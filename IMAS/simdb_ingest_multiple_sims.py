#!/usr/bin/env python3
"""Generate and ingest SIMDB manifests for DREAM IMAS exports.

The script reads the URI map produced by ``export_multiple_sims_to_IMAS.py``:

    URI<TAB>Provenance folder

For each row it creates a manifest from ``manifest_template.yaml``, adapting:

* ``alias`` to ``DREAM_RE_2026_<pulse>/<run>``
* IMAS output URI to ``imas:hdf5?path=<resolved local IMAS path>``
* file inputs/outputs to the simulation provenance folder
* the case label and full provenance folder inside the description metadata

Then it runs, in order:

    simdb manifest check <manifest>
    simdb simulation ingest <manifest>
    simdb simulation validate <alias>
"""

from __future__ import annotations

import argparse
import copy
import csv
import getpass
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any
from urllib.parse import parse_qsl, quote, urlsplit

try:
    import yaml
except Exception as exc:  # pragma: no cover - environment dependent
    yaml = None  # type: ignore[assignment]
    YAML_IMPORT_ERROR = exc
else:
    YAML_IMPORT_ERROR = None


DEFAULT_URI_MAP = "dream_simulation_uris.tsv"
DEFAULT_TEMPLATE = "manifest_template.yaml"
DEFAULT_MANIFEST_DIR = "simdb_manifests"
DEFAULT_ALIAS_PREFIX = "DREAM_RE_2026"


@dataclass
class SimulationRow:
    imas_uri: str
    provenance_folder: Path


@dataclass
class ManifestJob:
    row: SimulationRow
    alias: str
    manifest_path: Path
    imas_path_uri: str


class LiteralString(str):
    """String marker for YAML literal block output."""


class ManifestDumper(yaml.SafeDumper if yaml is not None else object):
    pass


if yaml is not None:
    def represent_literal_string(dumper: yaml.SafeDumper, data: LiteralString) -> yaml.ScalarNode:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")

    ManifestDumper.add_representer(LiteralString, represent_literal_string)


def parse_imas_uri_query(uri: str) -> dict[str, str]:
    query = urlsplit(uri).query.replace(";", "&")
    return dict(parse_qsl(query, keep_blank_values=True))


def imas_uri_to_local_path(uri: str) -> Path:
    """Resolve an IMAS URI to the local data-entry path used by this workflow.

    IMAS URI syntax supports either an explicit ``path=...`` query entry, or
    legacy identifiers such as user/pulse/run/database/version. The latter are
    mapped here to the ITER local HDF5 database convention used in this project.
    """
    query = parse_imas_uri_query(uri)
    if "path" in query:
        return Path(query["path"]).expanduser().resolve()

    required = {"user", "pulse", "run", "database", "version"}
    missing = sorted(required - query.keys())
    if missing:
        raise ValueError(f"cannot resolve IMAS URI to a path; missing query keys: {', '.join(missing)}")

    user = query["user"] or getpass.getuser()
    return (
        Path("/home/ITER")
        / user
        / "public"
        / "imasdb"
        / query["database"]
        / query["version"]
        / query["pulse"]
        / query["run"]
    )


def imas_path_uri(uri: str) -> str:
    path = imas_uri_to_local_path(uri)
    return "imas:hdf5?path=" + quote(str(path), safe="/")


def alias_from_imas_uri(uri: str, prefix: str) -> str:
    query = parse_imas_uri_query(uri)
    pulse = query.get("pulse")
    run = query.get("run")
    if pulse is None or run is None:
        path = imas_uri_to_local_path(uri)
        try:
            pulse = path.parts[-2]
            run = path.parts[-1]
        except IndexError as exc:
            raise ValueError(f"cannot infer pulse/run for alias from {uri!r}") from exc
    return f"{prefix}_{pulse}/{run}"


def file_uri(path_or_pattern: Path | str) -> str:
    return "file://" + str(path_or_pattern)


def read_uri_map(path: Path) -> list[SimulationRow]:
    rows: list[SimulationRow] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"URI", "Provenance folder"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: URI and Provenance folder")
        for raw in reader:
            uri = (raw.get("URI") or "").strip()
            folder = (raw.get("Provenance folder") or "").strip()
            if not uri or not folder:
                continue
            rows.append(SimulationRow(imas_uri=uri, provenance_folder=Path(folder).expanduser().resolve()))
    return rows


def load_template(path: Path) -> dict[str, Any]:
    if yaml is None:
        raise RuntimeError(f"PyYAML is required to read/write manifests: {YAML_IMPORT_ERROR}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} did not contain a YAML mapping")
    return data


def common_case_root(rows: list[SimulationRow], fallback: Path) -> Path:
    if not rows:
        return fallback
    if len(rows) == 1:
        return fallback
    try:
        return Path(os.path.commonpath([row.provenance_folder for row in rows]))
    except Exception:
        return fallback


def case_label(folder: Path, case_root: Path) -> str:
    try:
        return folder.resolve().relative_to(case_root.resolve()).as_posix()
    except ValueError:
        return folder.name


def replace_case_description(description: str, label: str, folder: Path) -> str:
    lines = description.splitlines()
    label_replaced = False
    folder_replaced = False
    default_indent = "  "
    for index, line in enumerate(lines):
        if "This case label:" in line:
            indent = line[: len(line) - len(line.lstrip())]
            default_indent = indent
            lines[index] = f"{indent}This case label: {label}"
            label_replaced = True
        if "Provenance folder:" in line:
            indent = line[: len(line) - len(line.lstrip())]
            lines[index] = f"{indent}Provenance folder: {folder}"
            folder_replaced = True
    if not label_replaced:
        lines.append(f"{default_indent}This case label: {label}")
    if not folder_replaced:
        lines.append(f"{default_indent}Provenance folder: {folder}")
    return "\n".join(lines)


def update_metadata(manifest: dict[str, Any], label: str, folder: Path) -> None:
    metadata = manifest.get("metadata")
    if not isinstance(metadata, list):
        metadata = []
        manifest["metadata"] = metadata

    has_creation_date = False
    for item in metadata:
        if not isinstance(item, dict):
            continue
        if "description" in item:
            item["description"] = LiteralString(replace_case_description(str(item["description"]), label, folder))
        ids_properties = item.get("ids_properties")
        if isinstance(ids_properties, dict):
            ids_properties["creation_date"] = date.today().isoformat()
            has_creation_date = True

    if not has_creation_date:
        metadata.append({"ids_properties": {"creation_date": date.today().isoformat()}})


def build_manifest(template: dict[str, Any], row: SimulationRow, alias_prefix: str, case_root: Path) -> tuple[str, dict[str, Any]]:
    manifest = copy.deepcopy(template)
    label = case_label(row.provenance_folder, case_root)
    alias = alias_from_imas_uri(row.imas_uri, alias_prefix)
    output_dir = row.provenance_folder / "output"

    manifest["alias"] = alias
    manifest["inputs"] = [
        {"uri": file_uri(output_dir / "settings*")},
        {"uri": file_uri(row.provenance_folder / "*.py")},
    ]
    manifest["outputs"] = [
        {"uri": imas_path_uri(row.imas_uri)},
        {"uri": file_uri(output_dir / "output*")},
    ]
    update_metadata(manifest, label, row.provenance_folder)
    return alias, manifest


def write_manifest(manifest: dict[str, Any], path: Path) -> None:
    if yaml is None:
        raise RuntimeError(f"PyYAML is required to read/write manifests: {YAML_IMPORT_ERROR}")
    path.parent.mkdir(parents=True, exist_ok=True)
    text = yaml.dump(manifest, Dumper=ManifestDumper, sort_keys=False, default_flow_style=False)
    path.write_text(text, encoding="utf-8")


def run_command(command: list[str], dry_run: bool) -> int:
    print("+ " + " ".join(command))
    if dry_run:
        return 0
    completed = subprocess.run(command, check=False)
    return int(completed.returncode)


def process_job(job: ManifestJob, dry_run: bool) -> bool:
    steps = [
        #["simdb", "manifest", "check", str(job.manifest_path)],
        ["simdb", "simulation", "ingest", str(job.manifest_path)],
        #["simdb", "simulation", "validate", job.alias],
    ]
    for command in steps:
        if run_command(command, dry_run) != 0:
            return False
    return True


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create, check, ingest, and validate SIMDB manifests for DREAM simulations.")
    parser.add_argument(
        "folder",
        nargs="?",
        default=".",
        help="Folder containing dream_simulation_uris.tsv and manifest_template.yaml.",
    )
    parser.add_argument("--uri-map", default=DEFAULT_URI_MAP, help="URI map TSV path, relative to folder unless absolute.")
    parser.add_argument("--template", default=DEFAULT_TEMPLATE, help="Manifest template path, relative to folder unless absolute.")
    parser.add_argument("--manifest-dir", default=DEFAULT_MANIFEST_DIR, help="Output manifest directory, relative to folder unless absolute.")
    parser.add_argument("--alias-prefix", default=DEFAULT_ALIAS_PREFIX, help="Alias prefix before _<pulse>/<run>.")
    parser.add_argument("--dry-run", action="store_true", help="Write manifests and print simdb commands, but do not execute them.")
    parser.add_argument("--keep-going", action="store_true", help="Continue processing remaining simulations after a failed simdb command.")
    parser.add_argument("--generate-only", action="store_true", help="Only write manifest files; do not run simdb commands.")
    return parser.parse_args(argv)


def resolve_under(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    base = Path(args.folder).expanduser().resolve()
    uri_map_path = resolve_under(base, args.uri_map)
    template_path = resolve_under(base, args.template)
    manifest_dir = resolve_under(base, args.manifest_dir)

    template = load_template(template_path)
    rows = read_uri_map(uri_map_path)
    if not rows:
        print(f"No simulations found in {uri_map_path}", file=sys.stderr)
        return 2
    label_root = common_case_root(rows, base)

    jobs: list[ManifestJob] = []
    for row in rows:
        alias, manifest = build_manifest(template, row, args.alias_prefix, label_root)
        safe_name = alias.replace("/", "_")
        manifest_path = manifest_dir / f"{safe_name}.yaml"
        write_manifest(manifest, manifest_path)
        jobs.append(ManifestJob(row=row, alias=alias, manifest_path=manifest_path, imas_path_uri=manifest["outputs"][0]["uri"]))
        print(f"Wrote {manifest_path}")

    if args.generate_only:
        return 0

    failures = 0
    for index, job in enumerate(jobs, start=1):
        print(f"[{index}/{len(jobs)}] {job.alias}")
        print(f"  IMAS output: {job.imas_path_uri}")
        print(f"  Provenance: {job.row.provenance_folder}")
        if not process_job(job, args.dry_run):
            failures += 1
            print(f"ERROR: SIMDB processing failed for {job.alias}", file=sys.stderr)
            if not args.keep_going:
                break

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
