import imas
from imas.util import get_data_dictionary_version

def read_uris(uri: str = None, uri_file: str = None) -> list[str]:
    if uri:
        return [uri]
    if uri_file:
        uris = []
        with open(uri_file, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                uris.append(s)
        if not uris:
            raise ValueError(f"No URIs found in file: {uri_file}")
        return uris
    raise ValueError("Provide either --uri or --uri_file")

def get_entries(uris: list[str], mode: str = "r", dd_version: str = "None"):
    entries = []
    
    for u in uris:
        if dd_version == "None":
          e = imas.DBEntry(uri=u, mode=mode)
        else:
          e = imas.DBEntry(uri=u, mode=mode, dd_version=dd_version)
        entries.append(e)
    return entries


def parse_imas_uri_kv(uri: str) -> dict:
    if "?" not in uri:
        return {}
    _, params = uri.split("?", 1)
    out = {}
    for part in params.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k.strip().lower()] = v.strip()
    return out

def make_entry_label(uri: str, idx: int) -> str:
    kv = parse_imas_uri_kv(uri)
    pulse = kv.get("pulse", "?")
    run = kv.get("run", "?")
    db = kv.get("database", "?")
    user = kv.get("user", "?")
    return f"{idx}: {pulse}/{run} {db} {user}"