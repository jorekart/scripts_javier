#!/usr/bin/env python3
"""Display JOREK grid and limiter data in preparation for patch editing.

This first implementation stage deliberately keeps the block geometry and file
parsers independent of Matplotlib.  A later interactive controller can build on
the data types and snapping candidates defined here without coupling persistent
state to plotting artists.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Any, Iterable, List, NamedTuple, Sequence

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_GRID_PATH = SCRIPT_DIR / "grid_polar.dat"
DEFAULT_LIMITER_PATH = SCRIPT_DIR / "simple_limiter.txt"
DEFAULT_OUTPUT_PATH = SCRIPT_DIR / "extension_blocks.in"

GRID_SOURCE = np.uint8(1)
LIMITER_SOURCE = np.uint8(2)


class Point(NamedTuple):
    """One point in the poloidal R-Z plane, in metres."""

    r: float
    z: float


@dataclass
class Block:
    """Geometry entered by the user for one JOREK extension block.

    ``left_points`` and ``right_points`` are ordered polylines.  They are not
    restricted to two points because patch jumps and shaped sides may require
    additional vertices.  The block's one-based JOREK index is intentionally
    not stored: it will be derived from its position in the document when the
    blocks are displayed or written.
    """

    n_ext_block: int
    left_points: List[Point] = field(default_factory=list)
    right_points: List[Point] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.n_ext_block <= 0:
            raise ValueError("n_ext_block must be a positive integer")


@dataclass(frozen=True)
class SnapCandidates:
    """Unique vertices available to a future point-snapping controller.

    ``points`` has shape ``(n, 2)``.  Each value in ``source_flags`` is a bit
    mask made from ``GRID_SOURCE`` and ``LIMITER_SOURCE``.  A coordinate shared
    by the grid and limiter is stored once with both bits set.
    """

    points: np.ndarray
    source_flags: np.ndarray

    def __post_init__(self) -> None:
        if self.points.ndim != 2 or self.points.shape[1] != 2:
            raise ValueError("snap candidate points must have shape (n, 2)")
        if self.source_flags.shape != (self.points.shape[0],):
            raise ValueError("source_flags must contain one value per point")


@dataclass
class BlockArtists:
    """Matplotlib artists belonging to one block, separate from its geometry."""

    left_line: Any
    right_line: Any
    block_label: Any
    point_labels: List[Any] = field(default_factory=list)


@dataclass(frozen=True)
class PointSelection:
    """Index-based reference to one editable point in the block model."""

    block_index: int
    side: str
    point_index: int


class SnapMatch(NamedTuple):
    """Nearest snap candidate found within the configured pixel radius."""

    candidate_index: int
    point: Point
    source_flags: int
    distance_pixels: float


class DataFormatError(ValueError):
    """Raised when an input file does not contain valid two-column data."""


class JorekFormatError(DataFormatError):
    """Raised when a saved JOREK block file is incomplete or inconsistent."""


_JOREK_SCALAR_NAMES = {
    "n_ext_block",
    "n_block_points_left",
    "n_block_points_right",
}
_JOREK_COORDINATE_NAMES = {
    "R_block_points_left",
    "Z_block_points_left",
    "R_block_points_right",
    "Z_block_points_right",
}
_JOREK_DEFINITION_RE = re.compile(
    r"^\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*"
    r"\(\s*(?P<block>[+-]?\d+)\s*"
    r"(?:(?P<comma>,)\s*(?P<point>[+-]?\d+)\s*)?\)"
    r"\s*=\s*(?P<value>[^,\s!#]+)\s*,?\s*(?:[!#].*)?$"
)


def format_jorek_blocks(blocks: Sequence[Block]) -> str:
    """Serialize blocks using exact one-based JOREK input variable names."""

    formatted_blocks: List[str] = []
    for block_index, block in enumerate(blocks, start=1):
        if block.n_ext_block <= 0:
            raise ValueError(f"block {block_index}: n_ext_block must be positive")

        lines = [
            f"n_ext_block         ({block_index:2d})    = {block.n_ext_block}",
            f"n_block_points_left ({block_index:2d})    = {len(block.left_points)}",
        ]
        for point_index, point in enumerate(block.left_points, start=1):
            if not (math.isfinite(point.r) and math.isfinite(point.z)):
                raise ValueError(
                    f"block {block_index} left point {point_index} is not finite"
                )
            lines.extend(
                (
                    f"R_block_points_left ({block_index:2d},{point_index})  = "
                    f"{point.r:+.8f}",
                    f"Z_block_points_left ({block_index:2d},{point_index})  = "
                    f"{point.z:+.8f}",
                )
            )

        lines.append(
            f"n_block_points_right({block_index:2d})    = {len(block.right_points)}"
        )
        for point_index, point in enumerate(block.right_points, start=1):
            if not (math.isfinite(point.r) and math.isfinite(point.z)):
                raise ValueError(
                    f"block {block_index} right point {point_index} is not finite"
                )
            lines.extend(
                (
                    f"R_block_points_right({block_index:2d},{point_index})  = "
                    f"{point.r:+.8f}",
                    f"Z_block_points_right({block_index:2d},{point_index})  = "
                    f"{point.z:+.8f}",
                )
            )
        formatted_blocks.append("\n".join(lines))

    if not formatted_blocks:
        return ""
    return "\n\n".join(formatted_blocks) + "\n"


def incomplete_block_warnings(blocks: Sequence[Block]) -> List[str]:
    """Describe blocks that do not yet have two points on each side."""

    if not blocks:
        return ["no blocks are defined"]
    warnings: List[str] = []
    for block_index, block in enumerate(blocks, start=1):
        if len(block.left_points) < 2:
            warnings.append(
                f"block {block_index} left side has {len(block.left_points)} point(s)"
            )
        if len(block.right_points) < 2:
            warnings.append(
                f"block {block_index} right side has {len(block.right_points)} point(s)"
            )
    return warnings


def write_jorek_blocks_atomic(path: Path | str, blocks: Sequence[Block]) -> None:
    """Atomically replace ``path`` with the formatted JOREK block document."""

    destination = Path(path)
    parent = destination.parent
    if not parent.is_dir():
        raise OSError(f"output directory does not exist: {parent}")

    content = format_jorek_blocks(blocks)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=str(parent)
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(
            file_descriptor, "w", encoding="utf-8", newline="\n"
        ) as temporary_file:
            temporary_file.write(content)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, destination)
    except BaseException:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass
        raise


def _parse_jorek_integer(
    path: Path, line_number: int, name: str, value: str, minimum: int
) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise JorekFormatError(
            f"{path}:{line_number}: {name} must be an integer, got {value!r}"
        ) from error
    if parsed < minimum:
        qualifier = "positive" if minimum == 1 else "non-negative"
        raise JorekFormatError(
            f"{path}:{line_number}: {name} must be {qualifier}, got {parsed}"
        )
    return parsed


def _parse_jorek_coordinate(
    path: Path, line_number: int, name: str, value: str
) -> float:
    try:
        parsed = float(value.replace("D", "E").replace("d", "e"))
    except ValueError as error:
        raise JorekFormatError(
            f"{path}:{line_number}: {name} must be a number, got {value!r}"
        ) from error
    if not math.isfinite(parsed):
        raise JorekFormatError(
            f"{path}:{line_number}: {name} must be finite, got {value!r}"
        )
    return parsed


def read_jorek_blocks(path: Path | str) -> List[Block]:
    """Load blocks from the exact JOREK assignment family written above.

    Whitespace and Fortran ``D`` exponents are tolerated.  Definitions must be
    unique, block and point indices must be contiguous and one-based, declared
    counts must match, and every point must contain both R and Z coordinates.
    """

    input_path = Path(path)
    definitions: dict[tuple[str, int, int | None], tuple[str, int]] = {}
    allowed_names = _JOREK_SCALAR_NAMES | _JOREK_COORDINATE_NAMES

    try:
        with input_path.open("r", encoding="utf-8") as input_file:
            for line_number, raw_line in enumerate(input_file, start=1):
                stripped = raw_line.strip()
                if not stripped or stripped.startswith(("#", "!")):
                    continue
                match = _JOREK_DEFINITION_RE.match(raw_line)
                if match is None:
                    raise JorekFormatError(
                        f"{input_path}:{line_number}: expected a JOREK block "
                        f"assignment, got {raw_line.rstrip()!r}"
                    )

                name = match.group("name")
                if name not in allowed_names:
                    raise JorekFormatError(
                        f"{input_path}:{line_number}: unsupported variable {name!r}"
                    )
                block_index = int(match.group("block"))
                if block_index <= 0:
                    raise JorekFormatError(
                        f"{input_path}:{line_number}: block index must be positive, "
                        f"got {block_index}"
                    )
                point_text = match.group("point")
                point_index = int(point_text) if point_text is not None else None

                if name in _JOREK_SCALAR_NAMES and point_index is not None:
                    raise JorekFormatError(
                        f"{input_path}:{line_number}: {name} accepts only a block index"
                    )
                if name in _JOREK_COORDINATE_NAMES:
                    if point_index is None:
                        raise JorekFormatError(
                            f"{input_path}:{line_number}: {name} requires block and "
                            "point indices"
                        )
                    if point_index <= 0:
                        raise JorekFormatError(
                            f"{input_path}:{line_number}: point index must be positive, "
                            f"got {point_index}"
                        )

                key = (name, block_index, point_index)
                if key in definitions:
                    previous_line = definitions[key][1]
                    raise JorekFormatError(
                        f"{input_path}:{line_number}: duplicate definition of {name} "
                        f"for block {block_index}"
                        + (
                            f", point {point_index}"
                            if point_index is not None
                            else ""
                        )
                        + f" (first defined on line {previous_line})"
                    )
                definitions[key] = (match.group("value"), line_number)
    except UnicodeError as error:
        raise JorekFormatError(f"{input_path}: file is not valid UTF-8") from error

    if not definitions:
        raise JorekFormatError(f"{input_path}: no JOREK block definitions found")

    block_indices = sorted({key[1] for key in definitions})
    expected_blocks = list(range(1, block_indices[-1] + 1))
    if block_indices != expected_blocks:
        raise JorekFormatError(
            f"{input_path}: block indices must be contiguous from 1; "
            f"found {block_indices}"
        )

    def required(
        name: str, block_index: int, point_index: int | None = None
    ) -> tuple[str, int]:
        key = (name, block_index, point_index)
        if key not in definitions:
            suffix = (
                f" for block {block_index}, point {point_index}"
                if point_index is not None
                else f" for block {block_index}"
            )
            raise JorekFormatError(f"{input_path}: missing {name}{suffix}")
        return definitions[key]

    blocks: List[Block] = []
    for block_index in expected_blocks:
        n_ext_value, n_ext_line = required("n_ext_block", block_index)
        n_ext_block = _parse_jorek_integer(
            input_path, n_ext_line, "n_ext_block", n_ext_value, minimum=1
        )
        count_left_value, count_left_line = required(
            "n_block_points_left", block_index
        )
        count_left = _parse_jorek_integer(
            input_path,
            count_left_line,
            "n_block_points_left",
            count_left_value,
            minimum=0,
        )
        count_right_value, count_right_line = required(
            "n_block_points_right", block_index
        )
        count_right = _parse_jorek_integer(
            input_path,
            count_right_line,
            "n_block_points_right",
            count_right_value,
            minimum=0,
        )

        side_points: dict[str, List[Point]] = {"left": [], "right": []}
        for side, count in (("left", count_left), ("right", count_right)):
            coordinate_names = {
                f"R_block_points_{side}",
                f"Z_block_points_{side}",
            }
            extra_indices = sorted(
                {
                    point_index
                    for name, defined_block, point_index in definitions
                    if defined_block == block_index
                    and name in coordinate_names
                    and point_index is not None
                    and point_index > count
                }
            )
            if extra_indices:
                extra_index = extra_indices[0]
                extra_name = next(
                    name
                    for name in coordinate_names
                    if (name, block_index, extra_index) in definitions
                )
                extra_line = definitions[(extra_name, block_index, extra_index)][1]
                raise JorekFormatError(
                    f"{input_path}:{extra_line}: {side} point index {extra_index} "
                    f"exceeds declared count {count} for block {block_index}"
                )

            for point_index in range(1, count + 1):
                r_name = f"R_block_points_{side}"
                z_name = f"Z_block_points_{side}"
                r_value, r_line = required(r_name, block_index, point_index)
                z_value, z_line = required(z_name, block_index, point_index)
                r_coordinate = _parse_jorek_coordinate(
                    input_path, r_line, r_name, r_value
                )
                z_coordinate = _parse_jorek_coordinate(
                    input_path, z_line, z_name, z_value
                )
                side_points[side].append(Point(r_coordinate, z_coordinate))

        blocks.append(
            Block(
                n_ext_block=n_ext_block,
                left_points=side_points["left"],
                right_points=side_points["right"],
            )
        )

    return blocks


def _coordinate_from_line(path: Path, line_number: int, raw_line: str) -> Point | None:
    """Parse one coordinate line, returning ``None`` for a comment line."""

    content = raw_line.split("#", 1)[0].strip()
    if not content:
        return None

    fields = content.split()
    if len(fields) != 2:
        raise DataFormatError(
            f"{path}:{line_number}: expected exactly 2 columns, got "
            f"{len(fields)}: {raw_line.rstrip()!r}"
        )

    try:
        r_value, z_value = (float(value) for value in fields)
    except ValueError as error:
        raise DataFormatError(
            f"{path}:{line_number}: invalid floating-point coordinate: "
            f"{raw_line.rstrip()!r}"
        ) from error

    if not (math.isfinite(r_value) and math.isfinite(z_value)):
        raise DataFormatError(
            f"{path}:{line_number}: coordinates must be finite: "
            f"{raw_line.rstrip()!r}"
        )

    return Point(r_value, z_value)


def read_segment_file(path: Path | str) -> List[np.ndarray]:
    """Read blank-line-separated, two-column polylines from ``path``.

    Any run of one or more blank lines terminates the current segment.  Full or
    trailing ``#`` comments are accepted and do not themselves split a segment.
    Every returned array has shape ``(n, 2)`` and dtype ``float64``.
    """

    input_path = Path(path)
    segments: List[np.ndarray] = []
    current_segment: List[Point] = []

    try:
        with input_path.open("r", encoding="utf-8") as input_file:
            for line_number, raw_line in enumerate(input_file, start=1):
                if not raw_line.strip():
                    if current_segment:
                        segments.append(np.asarray(current_segment, dtype=float))
                        current_segment = []
                    continue

                coordinate = _coordinate_from_line(input_path, line_number, raw_line)
                if coordinate is not None:
                    current_segment.append(coordinate)
    except UnicodeError as error:
        raise DataFormatError(f"{input_path}: file is not valid UTF-8") from error

    if current_segment:
        segments.append(np.asarray(current_segment, dtype=float))

    if not segments:
        raise DataFormatError(f"{input_path}: no coordinate data found")

    return segments


def read_limiter(path: Path | str) -> np.ndarray:
    """Read an ordered two-column limiter contour from ``path``.

    Blank and comment lines are ignored.  Repeated vertices and an explicit
    closing point are preserved because they are part of the plotted contour.
    """

    input_path = Path(path)
    points: List[Point] = []

    try:
        with input_path.open("r", encoding="utf-8") as input_file:
            for line_number, raw_line in enumerate(input_file, start=1):
                if not raw_line.strip():
                    continue
                coordinate = _coordinate_from_line(input_path, line_number, raw_line)
                if coordinate is not None:
                    points.append(coordinate)
    except UnicodeError as error:
        raise DataFormatError(f"{input_path}: file is not valid UTF-8") from error

    if not points:
        raise DataFormatError(f"{input_path}: no coordinate data found")

    return np.asarray(points, dtype=float)


def build_snap_candidates(
    grid_segments: Iterable[np.ndarray], limiter: np.ndarray
) -> SnapCandidates:
    """Return the unique union of grid and limiter vertices with source flags."""

    segment_arrays = [np.asarray(segment, dtype=float) for segment in grid_segments]
    for segment_number, segment in enumerate(segment_arrays, start=1):
        if segment.ndim != 2 or segment.shape[1] != 2:
            raise ValueError(f"grid segment {segment_number} must have shape (n, 2)")

    limiter_array = np.asarray(limiter, dtype=float)
    if limiter_array.ndim != 2 or limiter_array.shape[1] != 2:
        raise ValueError("limiter must have shape (n, 2)")

    nonempty_segments = [segment for segment in segment_arrays if len(segment)]
    grid_points = (
        np.concatenate(nonempty_segments, axis=0)
        if nonempty_segments
        else np.empty((0, 2), dtype=float)
    )

    point_sets = []
    flag_sets = []
    if len(grid_points):
        point_sets.append(grid_points)
        flag_sets.append(np.full(len(grid_points), GRID_SOURCE, dtype=np.uint8))
    if len(limiter_array):
        point_sets.append(limiter_array)
        flag_sets.append(np.full(len(limiter_array), LIMITER_SOURCE, dtype=np.uint8))

    if not point_sets:
        return SnapCandidates(
            points=np.empty((0, 2), dtype=float),
            source_flags=np.empty(0, dtype=np.uint8),
        )

    all_points = np.concatenate(point_sets, axis=0)
    all_flags = np.concatenate(flag_sets)
    unique_points, inverse = np.unique(all_points, axis=0, return_inverse=True)
    unique_flags = np.zeros(len(unique_points), dtype=np.uint8)
    np.bitwise_or.at(unique_flags, inverse, all_flags)

    return SnapCandidates(points=unique_points, source_flags=unique_flags)


class PatchBuilderApp:
    """Interactive Matplotlib controller for creating extension-block sides.

    Block geometry lives exclusively in :class:`Block` instances.  This class
    owns the transient selection state, widgets, callbacks, and Matplotlib
    artists.  It creates canvas callbacks once during initialization; artist
    refreshes never reconnect callbacks or recreate the figure.
    """

    LEFT_COLOR = "#1f77b4"
    RIGHT_COLOR = "#ff7f0e"
    ACTIVE_ALPHA = 1.0
    INACTIVE_ALPHA = 0.32

    def __init__(
        self,
        grid_segments: Sequence[np.ndarray],
        limiter: np.ndarray,
        snap_candidates: SnapCandidates,
        default_n_ext: int = 30,
        snap_pixels: float = 12.0,
        output_path: Path | str = DEFAULT_OUTPUT_PATH,
        initial_blocks: Sequence[Block] | None = None,
    ) -> None:
        try:
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.lines import Line2D
            from matplotlib.widgets import Button, RadioButtons, TextBox
        except ModuleNotFoundError as error:
            if error.name == "matplotlib" or (error.name or "").startswith(
                "matplotlib."
            ):
                raise RuntimeError(
                    "Matplotlib is required for the graphical display but is not "
                    "installed in this Python environment"
                ) from error
            raise

        if default_n_ext <= 0:
            raise ValueError("default_n_ext must be positive")
        if snap_pixels <= 0.0:
            raise ValueError("snap_pixels must be positive")

        self._plt = plt
        self._Button = Button
        self._RadioButtons = RadioButtons
        self._TextBox = TextBox
        plt.rcParams["keymap.save"] = [
            key
            for key in plt.rcParams.get("keymap.save", [])
            if key.lower() not in {"ctrl+s", "cmd+s", "super+s"}
        ]
        self.grid_segments = list(grid_segments)
        self.limiter = np.asarray(limiter, dtype=float)
        self.snap_candidates = snap_candidates
        self.default_n_ext = default_n_ext
        self.snap_pixels = snap_pixels
        self.output_path = Path(output_path)

        self.blocks = [
            Block(
                n_ext_block=block.n_ext_block,
                left_points=[Point(point.r, point.z) for point in block.left_points],
                right_points=[Point(point.r, point.z) for point in block.right_points],
            )
            for block in (initial_blocks or [])
        ]
        self.block_artists: List[BlockArtists] = []
        self.active_block_index: int | None = 0 if self.blocks else None
        self.active_side = "left"
        self.add_point_mode = False
        self.insert_point_mode = False
        self.selected_point: PointSelection | None = None
        self._dragging_point = False
        self._drag_has_moved = False
        self._drag_start_display: tuple[float, float] | None = None
        self._drag_last_match: SnapMatch | None = None
        self._updating_resolution_text = False
        self._updating_side_selector = False
        self._canvas_callback_ids: List[int] = []
        self._axes_callback_ids: List[int] = []
        self._snap_display_points = np.empty((0, 2), dtype=float)
        self._snap_transform_signature: tuple[float, ...] | None = None
        self._pending_incomplete_save_signature: str | None = None

        self.figure = plt.figure(figsize=(13.2, 9.0))
        self.axes = self.figure.add_axes([0.07, 0.08, 0.66, 0.85])

        grid_collection = LineCollection(
            self.grid_segments,
            colors="#53565a",
            linewidths=0.35,
            alpha=0.65,
            label="Initial grid",
            rasterized=True,
        )
        self.axes.add_collection(grid_collection)
        self.axes.plot(
            self.limiter[:, 0],
            self.limiter[:, 1],
            color="#d62728",
            linewidth=1.4,
            marker="o",
            markersize=3.5,
            markerfacecolor="white",
            markeredgewidth=0.8,
            label="Limiter / first wall",
            zorder=3,
        )
        self.axes.autoscale_view()
        self.axes.set_aspect("equal", adjustable="box")
        self.axes.set_xlabel("R [m]")
        self.axes.set_ylabel("Z [m]")
        self.axes.set_title(
            "JOREK extension-patch builder\n"
            f"{len(self.grid_segments):,} grid segments; "
            f"{len(self.snap_candidates.points):,} unique snap vertices"
        )
        self.axes.grid(True, color="#d9d9d9", linewidth=0.5, alpha=0.6)

        legend_handles = [
            Line2D([], [], color="#53565a", linewidth=0.8, label="Initial grid"),
            Line2D(
                [],
                [],
                color="#d62728",
                marker="o",
                markerfacecolor="white",
                markersize=4,
                label="Limiter / first wall",
            ),
            Line2D(
                [],
                [],
                color=self.LEFT_COLOR,
                marker="o",
                linewidth=2.0,
                label="Block left side",
            ),
            Line2D(
                [],
                [],
                color=self.RIGHT_COLOR,
                marker="s",
                linewidth=2.0,
                label="Block right side",
            ),
        ]
        self.axes.legend(handles=legend_handles, loc="best")

        self.snap_indicator = self.axes.plot(
            [],
            [],
            linestyle="none",
            marker="*",
            markersize=15,
            markeredgecolor="black",
            markeredgewidth=0.8,
            zorder=20,
            visible=False,
        )[0]
        self.snap_annotation = self.axes.annotate(
            "",
            xy=(0.0, 0.0),
            xytext=(12, 12),
            textcoords="offset points",
            fontsize=8.5,
            fontweight="bold",
            color="#202020",
            zorder=21,
            visible=False,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#555555",
                "alpha": 0.92,
            },
        )
        self.selected_indicator = self.axes.plot(
            [],
            [],
            linestyle="none",
            marker="o",
            markersize=14,
            markerfacecolor="none",
            markeredgecolor="#f4d03f",
            markeredgewidth=2.6,
            zorder=19,
            visible=False,
        )[0]

        self.block_artists = [self._new_block_artists() for _block in self.blocks]

        self._create_controls()
        self._connect_callbacks()
        self._sync_resolution_text()
        if self.blocks:
            self._update_display(
                f"Loaded {len(self.blocks)} block(s). Select a point or continue editing."
            )
        else:
            self._update_display("Ready. Add a block to begin.")

    @property
    def active_block(self) -> Block | None:
        """Return the active geometry block, if one exists."""

        if self.active_block_index is None:
            return None
        return self.blocks[self.active_block_index]

    def _create_controls(self) -> None:
        """Create the persistent control panel once."""

        panel_x = 0.77
        panel_width = 0.20
        button_height = 0.042

        self.figure.text(
            panel_x,
            0.965,
            "Block controls",
            fontsize=13,
            fontweight="bold",
            ha="left",
        )

        add_axes = self.figure.add_axes([panel_x, 0.905, panel_width, button_height])
        delete_axes = self.figure.add_axes(
            [panel_x, 0.855, panel_width, button_height]
        )
        save_axes = self.figure.add_axes([panel_x, 0.805, panel_width, button_height])
        previous_axes = self.figure.add_axes(
            [panel_x, 0.750, 0.095, button_height]
        )
        next_axes = self.figure.add_axes(
            [panel_x + 0.105, 0.750, 0.095, button_height]
        )
        resolution_axes = self.figure.add_axes(
            [panel_x, 0.690, panel_width, button_height]
        )
        side_axes = self.figure.add_axes([panel_x, 0.560, panel_width, 0.105])
        mode_axes = self.figure.add_axes([panel_x, 0.500, panel_width, 0.045])
        delete_point_axes = self.figure.add_axes(
            [panel_x, 0.450, panel_width, button_height]
        )
        insert_point_axes = self.figure.add_axes(
            [panel_x, 0.400, panel_width, button_height]
        )
        earlier_axes = self.figure.add_axes(
            [panel_x, 0.350, 0.095, button_height]
        )
        later_axes = self.figure.add_axes(
            [panel_x + 0.105, 0.350, 0.095, button_height]
        )

        self.add_block_button = self._Button(add_axes, "Add block")
        self.delete_block_button = self._Button(delete_axes, "Delete active block")
        self.save_button = self._Button(save_axes, "Save blocks (Ctrl+S)")
        self.previous_block_button = self._Button(previous_axes, "Previous")
        self.next_block_button = self._Button(next_axes, "Next")
        self.resolution_textbox = self._TextBox(
            resolution_axes,
            "n_ext_block  ",
            initial=str(self.default_n_ext),
        )
        side_axes.set_title("Active side", fontsize=10, loc="left", pad=3)
        self.side_selector = self._RadioButtons(side_axes, ("Left", "Right"), active=0)
        self.add_mode_button = self._Button(mode_axes, "Start adding points")
        self.delete_point_button = self._Button(
            delete_point_axes, "Delete selected point"
        )
        self.insert_point_button = self._Button(
            insert_point_axes, "Insert after selected"
        )
        self.move_earlier_button = self._Button(earlier_axes, "Move earlier")
        self.move_later_button = self._Button(later_axes, "Move later")

        self.state_text = self.figure.text(
            panel_x,
            0.295,
            "",
            fontsize=9.8,
            va="top",
            ha="left",
            linespacing=1.30,
        )
        self.message_text = self.figure.text(
            panel_x,
            0.180,
            "",
            fontsize=9.0,
            va="top",
            ha="left",
            color="#303030",
            wrap=True,
        )
        self.figure.text(
            panel_x,
            0.110,
            "Active-block points: select/drag.\n"
            "Centre labels change active block.\n"
            "Add/insert uses the active side.\n"
            "The snap star identifies its source.\n"
            "Pan/zoom is one-shot; Esc exits.\n"
            "Ctrl+S saves to the --output path.",
            fontsize=8.5,
            va="top",
            ha="left",
            color="#555555",
            linespacing=1.35,
        )

    def _connect_callbacks(self) -> None:
        """Connect each widget and canvas callback exactly once."""

        self.add_block_button.on_clicked(self._on_add_block)
        self.delete_block_button.on_clicked(self._on_delete_block)
        self.save_button.on_clicked(self._on_save_requested)
        self.previous_block_button.on_clicked(self._on_previous_block)
        self.next_block_button.on_clicked(self._on_next_block)
        self.resolution_textbox.on_submit(self._on_resolution_submitted)
        self.side_selector.on_clicked(self._on_side_selected)
        self.add_mode_button.on_clicked(self._on_toggle_add_mode)
        self.delete_point_button.on_clicked(self._on_delete_selected_point)
        self.insert_point_button.on_clicked(self._on_toggle_insert_mode)
        self.move_earlier_button.on_clicked(self._on_move_point_earlier)
        self.move_later_button.on_clicked(self._on_move_point_later)

        for event_name, callback in (
            ("button_press_event", self._on_canvas_click),
            ("motion_notify_event", self._on_canvas_motion),
            ("button_release_event", self._on_canvas_release),
            ("resize_event", self._mark_snap_cache_dirty),
            ("key_press_event", self._on_key_press),
        ):
            callback_id = self.figure.canvas.mpl_connect(event_name, callback)
            self._canvas_callback_ids.append(callback_id)

        for event_name in ("xlim_changed", "ylim_changed"):
            callback_id = self.axes.callbacks.connect(
                event_name, self._mark_snap_cache_dirty
            )
            self._axes_callback_ids.append(callback_id)

    def _new_block_artists(self) -> BlockArtists:
        left_line = self.axes.plot(
            [],
            [],
            color=self.LEFT_COLOR,
            marker="o",
            markerfacecolor="white",
            markeredgewidth=1.5,
            linewidth=2.2,
            zorder=7,
        )[0]
        right_line = self.axes.plot(
            [],
            [],
            color=self.RIGHT_COLOR,
            marker="s",
            markerfacecolor="white",
            markeredgewidth=1.5,
            linewidth=2.2,
            zorder=7,
        )[0]
        block_label = self.axes.text(
            0.0,
            0.0,
            "",
            fontsize=9.5,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=9,
            visible=False,
            bbox={
                "boxstyle": "round,pad=0.22",
                "facecolor": "white",
                "edgecolor": "#333333",
                "alpha": 0.85,
            },
        )
        return BlockArtists(left_line, right_line, block_label)

    def _on_add_block(self, _event: Any) -> None:
        self._stop_point_placement()
        self.blocks.append(Block(n_ext_block=self.default_n_ext))
        self.block_artists.append(self._new_block_artists())
        self.active_block_index = len(self.blocks) - 1
        self.selected_point = None
        self.insert_point_mode = False
        self._sync_resolution_text()
        self._update_display(f"Added block {self.active_block_index + 1}.")

    def _on_delete_block(self, _event: Any) -> None:
        if self.active_block_index is None:
            self._update_display("There is no block to delete.")
            return

        self._stop_point_placement()
        deleted_index = self.active_block_index
        artists = self.block_artists.pop(deleted_index)
        artists.left_line.remove()
        artists.right_line.remove()
        artists.block_label.remove()
        for label in artists.point_labels:
            label.remove()
        del self.blocks[deleted_index]

        if self.selected_point is not None:
            if self.selected_point.block_index == deleted_index:
                self.selected_point = None
            elif self.selected_point.block_index > deleted_index:
                self.selected_point = PointSelection(
                    block_index=self.selected_point.block_index - 1,
                    side=self.selected_point.side,
                    point_index=self.selected_point.point_index,
                )

        if self.blocks:
            self.active_block_index = min(deleted_index, len(self.blocks) - 1)
        else:
            self.active_block_index = None
            self.add_point_mode = False
            self.insert_point_mode = False
            self.selected_point = None
        self._dragging_point = False
        self._sync_resolution_text()
        self._update_display(f"Deleted block {deleted_index + 1}.")

    def _on_previous_block(self, _event: Any) -> None:
        if not self.blocks:
            self._update_display("There are no blocks to select.")
            return
        current = self.active_block_index if self.active_block_index is not None else 0
        self._stop_point_placement()
        self.active_block_index = (current - 1) % len(self.blocks)
        self.selected_point = None
        self.insert_point_mode = False
        self._sync_resolution_text()
        self._update_display(f"Selected block {self.active_block_index + 1}.")

    def _on_next_block(self, _event: Any) -> None:
        if not self.blocks:
            self._update_display("There are no blocks to select.")
            return
        current = self.active_block_index if self.active_block_index is not None else -1
        self._stop_point_placement()
        self.active_block_index = (current + 1) % len(self.blocks)
        self.selected_point = None
        self.insert_point_mode = False
        self._sync_resolution_text()
        self._update_display(f"Selected block {self.active_block_index + 1}.")

    def _on_resolution_submitted(self, value: str) -> None:
        if self._updating_resolution_text:
            return
        block = self.active_block
        if block is None:
            self._sync_resolution_text()
            self._update_display("Add a block before setting n_ext_block.")
            return
        try:
            resolution = int(value.strip())
            if resolution <= 0:
                raise ValueError
        except ValueError:
            self._sync_resolution_text()
            self._update_display("n_ext_block must be a positive integer.")
            return

        block.n_ext_block = resolution
        self._sync_resolution_text()
        self._update_display(
            f"Block {self.active_block_index + 1}: n_ext_block = {resolution}."
        )

    def _on_side_selected(self, label: str) -> None:
        selected_side = label.lower()
        if self._updating_side_selector:
            self.active_side = selected_side
            return
        side_changed = selected_side != self.active_side
        self.active_side = selected_side
        if side_changed:
            self._stop_point_placement()
        if (
            self.selected_point is not None
            and self.selected_point.side != self.active_side
        ):
            self.selected_point = None
        self.insert_point_mode = False
        self._update_display(
            f"Active side changed to {self.active_side}. "
            f"Press Start adding {self.active_side.upper()} points to continue."
        )

    def _stop_point_placement(self) -> None:
        """Return to selection/edit mode after changing block or side."""

        self.add_point_mode = False
        self.insert_point_mode = False
        self._hide_snap_indicator()

    def _on_save_requested(self, _event: Any = None) -> None:
        """Save complete data immediately, or require a second draft confirmation."""

        try:
            document_signature = format_jorek_blocks(self.blocks)
        except ValueError as error:
            self._pending_incomplete_save_signature = None
            self._update_display(f"Cannot save: {error}")
            return

        warnings = incomplete_block_warnings(self.blocks)
        if warnings and self._pending_incomplete_save_signature != document_signature:
            self._pending_incomplete_save_signature = document_signature
            shown_warnings = "; ".join(warnings[:3])
            if len(warnings) > 3:
                shown_warnings += f"; and {len(warnings) - 3} more"
            self._update_display(
                "Incomplete draft warning: "
                f"{shown_warnings}. Press Confirm incomplete save or Ctrl+S again "
                "to save it anyway."
            )
            return

        try:
            write_jorek_blocks_atomic(self.output_path, self.blocks)
        except (OSError, ValueError) as error:
            self._update_display(f"Save failed for {self.output_path}: {error}")
            return

        self._pending_incomplete_save_signature = None
        saved_kind = "incomplete draft" if warnings else "blocks"
        self._update_display(
            f"Saved {saved_kind} atomically to {self.output_path}."
        )

    def _on_toggle_add_mode(self, _event: Any) -> None:
        if not self.blocks:
            self.add_point_mode = False
            self._update_display("Add a block before entering add-point mode.")
            return
        self.add_point_mode = not self.add_point_mode
        if self.add_point_mode:
            self.insert_point_mode = False
        message = (
            f"Click empty space to append points to the {self.active_side} side."
            if self.add_point_mode
            else "Add-point mode stopped."
        )
        self._update_display(message)

    def _on_canvas_click(self, event: Any) -> None:
        if event.inaxes is not self.axes:
            return
        if event.button != 1 or event.xdata is None or event.ydata is None:
            return

        if self._toolbar_is_active():
            return

        display_x, display_y = self._event_display_coordinates(event)
        clicked_block_index = self._find_clicked_block_label(display_x, display_y)
        if clicked_block_index is not None:
            self._stop_point_placement()
            self.active_block_index = clicked_block_index
            self.selected_point = None
            self._sync_resolution_text()
            self._update_display(
                f"Selected block {clicked_block_index + 1} from its centre label."
            )
            return

        nearest_selection = self._find_nearest_block_point(display_x, display_y)
        if nearest_selection is not None:
            selection, _distance = nearest_selection
            if self.add_point_mode:
                self.add_point_mode = False
                self._hide_snap_indicator()
            self._select_point(selection)
            self._dragging_point = True
            self._drag_has_moved = False
            self._drag_start_display = (display_x, display_y)
            self._drag_last_match = None
            self._update_display(
                f"Selected block {selection.block_index + 1} "
                f"{selection.side} point {selection.point_index + 1}; drag to move."
            )
            return

        block = self.active_block
        if block is None and (self.add_point_mode or self.insert_point_mode):
            self.add_point_mode = False
            self.insert_point_mode = False
            self._update_display("Add a block before adding points.")
            return

        if self.insert_point_mode:
            selection = self.selected_point
            points = self._points_for_selection(selection)
            if selection is None or points is None:
                self.insert_point_mode = False
                self._update_display("Select a point before inserting after it.")
                return
            point, match = self._snap_or_free_point(event)
            insert_index = selection.point_index + 1
            points.insert(insert_index, point)
            self.selected_point = PointSelection(
                selection.block_index, selection.side, insert_index
            )
            self.insert_point_mode = False
            self._update_snap_indicator(match)
            self._update_display(
                f"Inserted {self._snap_description(match)} as block "
                f"{selection.block_index + 1} {selection.side} point "
                f"{insert_index + 1}."
            )
            return

        if self.add_point_mode and block is not None:
            point, match = self._snap_or_free_point(event)
            points = (
                block.left_points
                if self.active_side == "left"
                else block.right_points
            )
            points.append(point)
            self.selected_point = PointSelection(
                self.active_block_index, self.active_side, len(points) - 1
            )
            self._update_snap_indicator(match)
            self._update_display(
                f"Added {self._snap_description(match)} as block "
                f"{self.active_block_index + 1} {self.active_side} point "
                f"{len(points)}."
            )
            return

        self.selected_point = None
        self._update_display("Selection cleared.")

    def _toolbar_is_active(self) -> bool:
        toolbar = getattr(self.figure.canvas, "toolbar", None)
        if toolbar is None:
            return False
        mode = getattr(toolbar, "mode", "")
        return bool(mode) and str(mode).lower() not in {"none", "_mode.none"}

    def _deactivate_toolbar_mode(self) -> bool:
        """Turn a persistent Matplotlib pan/zoom tool back into edit mode."""

        toolbar = getattr(self.figure.canvas, "toolbar", None)
        if toolbar is None or not self._toolbar_is_active():
            return False

        mode_name = str(getattr(toolbar, "mode", "")).lower()
        if "pan" in mode_name:
            toolbar.pan()
        elif "zoom" in mode_name:
            toolbar.zoom()
        else:
            return False
        return True

    def _resume_editing_after_view_change(self) -> None:
        self._deactivate_toolbar_mode()
        self._mark_snap_cache_dirty()
        self._dragging_point = False
        self._drag_has_moved = False
        self._hide_snap_indicator()
        self._update_display(
            "View changed. Pan/zoom is off; point selection and snapping are ready."
        )

    def _on_key_press(self, event: Any) -> None:
        key = str(getattr(event, "key", "")).lower()
        if key in {"ctrl+s", "cmd+s", "super+s"}:
            self._on_save_requested()
            return
        if key == "escape" and self._toolbar_is_active():
            self._resume_editing_after_view_change()

    def _event_display_coordinates(self, event: Any) -> tuple[float, float]:
        display_x = getattr(event, "x", None)
        display_y = getattr(event, "y", None)
        if display_x is None or display_y is None:
            display_x, display_y = self.axes.transData.transform(
                (float(event.xdata), float(event.ydata))
            )
        return float(display_x), float(display_y)

    def _snap_transform_state(self) -> tuple[float, ...]:
        """Return the view properties affecting the data-to-pixel transform."""

        values = (
            *self.axes.get_xlim(),
            *self.axes.get_ylim(),
            *self.axes.bbox.bounds,
            self.figure.dpi,
        )
        return tuple(float(value) for value in values)

    def _mark_snap_cache_dirty(self, _event: Any = None) -> None:
        self._snap_transform_signature = None

    def _ensure_snap_display_points(self) -> None:
        """Transform snap candidates after any view, layout, or DPI change."""

        signature = self._snap_transform_state()
        if signature == self._snap_transform_signature:
            return
        if len(self.snap_candidates.points):
            self._snap_display_points = self.axes.transData.transform(
                self.snap_candidates.points
            )
        else:
            self._snap_display_points = np.empty((0, 2), dtype=float)
        self._snap_transform_signature = signature

    def _find_snap_candidate(
        self, display_x: float, display_y: float
    ) -> SnapMatch | None:
        """Find the nearest grid/limiter vertex within ``snap_pixels``."""

        self._ensure_snap_display_points()
        if not len(self._snap_display_points):
            return None

        delta = self._snap_display_points - np.asarray((display_x, display_y))
        squared_distances = np.einsum("ij,ij->i", delta, delta)
        candidate_index = int(np.argmin(squared_distances))
        distance_pixels = float(math.sqrt(squared_distances[candidate_index]))
        if distance_pixels > self.snap_pixels:
            return None

        r_value, z_value = self.snap_candidates.points[candidate_index]
        return SnapMatch(
            candidate_index=candidate_index,
            point=Point(float(r_value), float(z_value)),
            source_flags=int(self.snap_candidates.source_flags[candidate_index]),
            distance_pixels=distance_pixels,
        )

    @staticmethod
    def _snap_source_name(source_flags: int) -> str:
        has_grid = bool(source_flags & int(GRID_SOURCE))
        has_limiter = bool(source_flags & int(LIMITER_SOURCE))
        if has_grid and has_limiter:
            return "grid + limiter vertex"
        if has_limiter:
            return "limiter vertex"
        return "grid vertex"

    @staticmethod
    def _snap_source_color(source_flags: int) -> str:
        has_grid = bool(source_flags & int(GRID_SOURCE))
        has_limiter = bool(source_flags & int(LIMITER_SOURCE))
        if has_grid and has_limiter:
            return "#9467bd"
        if has_limiter:
            return "#d62728"
        return "#2ca02c"

    def _update_snap_indicator(self, match: SnapMatch | None) -> None:
        if match is None:
            self._hide_snap_indicator()
            return
        color = self._snap_source_color(match.source_flags)
        source_name = self._snap_source_name(match.source_flags)
        self.snap_indicator.set_data([match.point.r], [match.point.z])
        self.snap_indicator.set_markerfacecolor(color)
        self.snap_indicator.set_visible(True)
        self.snap_annotation.xy = match.point
        self.snap_annotation.set_text(
            f"Snap: {source_name}\n{match.distance_pixels:.1f} px"
        )
        self.snap_annotation.get_bbox_patch().set_edgecolor(color)
        self.snap_annotation.set_visible(True)

    def _hide_snap_indicator(self) -> None:
        self.snap_indicator.set_visible(False)
        self.snap_annotation.set_visible(False)

    def _snap_or_free_point(self, event: Any) -> tuple[Point, SnapMatch | None]:
        display_x, display_y = self._event_display_coordinates(event)
        match = self._find_snap_candidate(display_x, display_y)
        if match is not None:
            return match.point, match
        return Point(float(event.xdata), float(event.ydata)), None

    def _snap_description(self, match: SnapMatch | None) -> str:
        if match is None:
            return "free point"
        return f"point snapped to a {self._snap_source_name(match.source_flags)}"

    def _find_nearest_block_point(
        self, display_x: float, display_y: float
    ) -> tuple[PointSelection, float] | None:
        """Find a vertex of the active block using a pixel-space tolerance."""

        if self.active_block_index is None:
            return None
        nearest_selection: PointSelection | None = None
        nearest_squared_distance = math.inf
        cursor = np.asarray((display_x, display_y))
        block = self.blocks[self.active_block_index]
        other_side = "right" if self.active_side == "left" else "left"

        for side in (self.active_side, other_side):
            points = block.left_points if side == "left" else block.right_points
            if not points:
                continue
            display_points = self.axes.transData.transform(
                np.asarray(points, dtype=float)
            )
            delta = display_points - cursor
            squared_distances = np.einsum("ij,ij->i", delta, delta)
            point_index = int(np.argmin(squared_distances))
            squared_distance = float(squared_distances[point_index])
            if squared_distance < nearest_squared_distance:
                nearest_squared_distance = squared_distance
                nearest_selection = PointSelection(
                    self.active_block_index, side, point_index
                )

        selection_tolerance = max(7.0, self.snap_pixels)
        if (
            nearest_selection is None
            or nearest_squared_distance > selection_tolerance**2
        ):
            return None
        return nearest_selection, math.sqrt(nearest_squared_distance)

    def _find_clicked_block_label(
        self, display_x: float, display_y: float
    ) -> int | None:
        """Return the block whose visible centre label contains the click."""

        renderer = self.figure.canvas.get_renderer()
        containing_labels: List[tuple[float, int]] = []
        for block_index, artists in enumerate(self.block_artists):
            if not artists.block_label.get_visible():
                continue
            bounds = artists.block_label.get_window_extent(renderer=renderer)
            if bounds.contains(display_x, display_y):
                distance = math.hypot(
                    display_x - (bounds.x0 + bounds.x1) / 2.0,
                    display_y - (bounds.y0 + bounds.y1) / 2.0,
                )
                containing_labels.append((distance, block_index))

        if not containing_labels:
            return None
        return min(containing_labels)[1]

    def _points_for_selection(
        self, selection: PointSelection | None
    ) -> List[Point] | None:
        if selection is None or not 0 <= selection.block_index < len(self.blocks):
            return None
        block = self.blocks[selection.block_index]
        points = block.left_points if selection.side == "left" else block.right_points
        if not 0 <= selection.point_index < len(points):
            return None
        return points

    def _selected_model_point(self) -> Point | None:
        points = self._points_for_selection(self.selected_point)
        if points is None or self.selected_point is None:
            return None
        return points[self.selected_point.point_index]

    def _sync_side_selector(self) -> None:
        desired_index = 0 if self.active_side == "left" else 1
        self._updating_side_selector = True
        try:
            self.side_selector.set_active(desired_index)
        finally:
            self._updating_side_selector = False

    def _select_point(self, selection: PointSelection) -> None:
        self.selected_point = selection
        self.active_block_index = selection.block_index
        self.active_side = selection.side
        self._sync_side_selector()
        self._sync_resolution_text()

    def _on_delete_selected_point(self, _event: Any) -> None:
        selection = self.selected_point
        points = self._points_for_selection(selection)
        if selection is None or points is None:
            self._update_display("Select a point before deleting it.")
            return

        deleted_number = selection.point_index + 1
        del points[selection.point_index]
        if points:
            self.selected_point = PointSelection(
                selection.block_index,
                selection.side,
                min(selection.point_index, len(points) - 1),
            )
        else:
            self.selected_point = None
        self.insert_point_mode = False
        self._update_display(
            f"Deleted block {selection.block_index + 1} {selection.side} "
            f"point {deleted_number}."
        )

    def _on_toggle_insert_mode(self, _event: Any) -> None:
        if self._points_for_selection(self.selected_point) is None:
            self.insert_point_mode = False
            self._update_display("Select a point before inserting after it.")
            return
        self.insert_point_mode = not self.insert_point_mode
        if self.insert_point_mode:
            self.add_point_mode = False
            message = "Click empty space to insert one point after the selection."
        else:
            message = "Insert-point mode cancelled."
        self._update_display(message)

    def _on_move_point_earlier(self, _event: Any) -> None:
        selection = self.selected_point
        points = self._points_for_selection(selection)
        if selection is None or points is None:
            self._update_display("Select a point before reordering it.")
            return
        if selection.point_index == 0:
            self._update_display("The selected point is already first.")
            return
        new_index = selection.point_index - 1
        points[new_index], points[selection.point_index] = (
            points[selection.point_index],
            points[new_index],
        )
        self.selected_point = PointSelection(
            selection.block_index, selection.side, new_index
        )
        self._update_display(f"Moved the selected point to position {new_index + 1}.")

    def _on_move_point_later(self, _event: Any) -> None:
        selection = self.selected_point
        points = self._points_for_selection(selection)
        if selection is None or points is None:
            self._update_display("Select a point before reordering it.")
            return
        if selection.point_index == len(points) - 1:
            self._update_display("The selected point is already last.")
            return
        new_index = selection.point_index + 1
        points[new_index], points[selection.point_index] = (
            points[selection.point_index],
            points[new_index],
        )
        self.selected_point = PointSelection(
            selection.block_index, selection.side, new_index
        )
        self._update_display(f"Moved the selected point to position {new_index + 1}.")

    def _on_canvas_motion(self, event: Any) -> None:
        if event.inaxes is not self.axes or event.xdata is None or event.ydata is None:
            if not self._dragging_point:
                self._hide_snap_indicator()
                self.figure.canvas.draw_idle()
            return
        if self._toolbar_is_active():
            return

        display_x, display_y = self._event_display_coordinates(event)
        if self._dragging_point:
            selection = self.selected_point
            points = self._points_for_selection(selection)
            if selection is None or points is None:
                self._dragging_point = False
                return
            if self._drag_start_display is not None:
                motion = math.hypot(
                    display_x - self._drag_start_display[0],
                    display_y - self._drag_start_display[1],
                )
                if motion < 2.0 and not self._drag_has_moved:
                    return
            self._drag_has_moved = True
            point, match = self._snap_or_free_point(event)
            points[selection.point_index] = point
            self._drag_last_match = match
            self._update_snap_indicator(match)
            self._update_block_artists()
            self.figure.canvas.draw_idle()
            return

        if self.add_point_mode or self.insert_point_mode:
            match = self._find_snap_candidate(display_x, display_y)
            self._update_snap_indicator(match)
            self.figure.canvas.draw_idle()
        elif self.snap_indicator.get_visible():
            self._hide_snap_indicator()
            self.figure.canvas.draw_idle()

    def _on_canvas_release(self, event: Any) -> None:
        if self._toolbar_is_active():
            self._resume_editing_after_view_change()
            return
        if not self._dragging_point:
            return
        selection = self.selected_point
        self._dragging_point = False

        if (
            self._drag_has_moved
            and selection is not None
            and event.inaxes is self.axes
            and event.xdata is not None
            and event.ydata is not None
        ):
            points = self._points_for_selection(selection)
            if points is not None:
                point, match = self._snap_or_free_point(event)
                points[selection.point_index] = point
                self._drag_last_match = match
                self._update_snap_indicator(match)
                self._update_display(
                    f"Moved block {selection.block_index + 1} {selection.side} "
                    f"point {selection.point_index + 1} to a "
                    f"{self._snap_description(match)}."
                )
                return

        self._update_display(
            "Point selected." if not self._drag_has_moved else "Point move finished."
        )

    def _sync_resolution_text(self) -> None:
        value = (
            str(self.active_block.n_ext_block)
            if self.active_block is not None
            else str(self.default_n_ext)
        )
        if self.resolution_textbox.text == value:
            return
        self._updating_resolution_text = True
        try:
            self.resolution_textbox.set_val(value)
        finally:
            self._updating_resolution_text = False

    @staticmethod
    def _line_coordinates(points: Sequence[Point]) -> tuple[np.ndarray, np.ndarray]:
        if not points:
            return np.empty(0), np.empty(0)
        coordinates = np.asarray(points, dtype=float)
        return coordinates[:, 0], coordinates[:, 1]

    def _update_block_artists(self) -> None:
        """Update block lines and labels in place without changing axes limits."""

        for index, (block, artists) in enumerate(
            zip(self.blocks, self.block_artists)
        ):
            is_active = index == self.active_block_index
            alpha = self.ACTIVE_ALPHA if is_active else self.INACTIVE_ALPHA
            linewidth = 2.5 if is_active else 1.25
            markersize = 7.0 if is_active else 4.5
            zorder = 8 if is_active else 5

            left_r, left_z = self._line_coordinates(block.left_points)
            right_r, right_z = self._line_coordinates(block.right_points)
            artists.left_line.set_data(left_r, left_z)
            artists.right_line.set_data(right_r, right_z)
            for line in (artists.left_line, artists.right_line):
                line.set_alpha(alpha)
                line.set_linewidth(linewidth)
                line.set_markersize(markersize)
                line.set_zorder(zorder)

            for point_label in artists.point_labels:
                point_label.remove()
            artists.point_labels.clear()

            for side_name, points, color in (
                ("L", block.left_points, self.LEFT_COLOR),
                ("R", block.right_points, self.RIGHT_COLOR),
            ):
                for point_index, point in enumerate(points, start=1):
                    side = "left" if side_name == "L" else "right"
                    is_selected = self.selected_point == PointSelection(
                        index, side, point_index - 1
                    )
                    label = self.axes.annotate(
                        f"{side_name}{point_index}",
                        xy=point,
                        xytext=(5, 5),
                        textcoords="offset points",
                        fontsize=8.5 if is_active else 7.0,
                        fontweight="bold" if is_active else "normal",
                        color=color,
                        alpha=alpha,
                        zorder=zorder + 1,
                        bbox=(
                            {
                                "boxstyle": "round,pad=0.16",
                                "facecolor": "#fff3a6",
                                "edgecolor": "#9a7d0a",
                                "alpha": 0.95,
                            }
                            if is_selected
                            else None
                        ),
                    )
                    artists.point_labels.append(label)

            all_points = block.left_points + block.right_points
            if all_points:
                label_r = sum(point.r for point in all_points) / len(all_points)
                label_z = sum(point.z for point in all_points) / len(all_points)
                artists.block_label.set_position((label_r, label_z))
                artists.block_label.set_text(f"Block {index + 1}")
                artists.block_label.set_alpha(0.95 if is_active else 0.45)
                artists.block_label.set_zorder(zorder + 2)
                artists.block_label.set_visible(True)
            else:
                artists.block_label.set_visible(False)

        selected_model_point = self._selected_model_point()
        if selected_model_point is None:
            self.selected_indicator.set_visible(False)
        else:
            self.selected_indicator.set_data(
                [selected_model_point.r], [selected_model_point.z]
            )
            self.selected_indicator.set_visible(True)

    def _update_display(self, message: str) -> None:
        """Refresh controls and block artists while preserving the plot view."""

        if self._pending_incomplete_save_signature is not None:
            try:
                current_signature = format_jorek_blocks(self.blocks)
            except ValueError:
                current_signature = None
            if current_signature != self._pending_incomplete_save_signature:
                self._pending_incomplete_save_signature = None

        if self.active_block_index is None:
            active_description = "none"
            count_description = f"0 / {len(self.blocks)}"
        else:
            active_description = str(self.active_block_index + 1)
            count_description = f"{self.active_block_index + 1} / {len(self.blocks)}"

        if self.insert_point_mode:
            mode_description = "INSERT POINT"
        elif self.add_point_mode:
            mode_description = "ADD POINTS"
        else:
            mode_description = "SELECT / EDIT"
        if self.selected_point is None:
            selected_description = "none"
        else:
            selected_description = (
                f"B{self.selected_point.block_index + 1} "
                f"{self.selected_point.side[0].upper()}"
                f"{self.selected_point.point_index + 1}"
            )
        self.state_text.set_text(
            f"Active block: {count_description}\n"
            f"Active side: {self.active_side.upper()}\n"
            f"Mode: {mode_description}\n"
            f"Selected point: {selected_description}"
        )
        self.message_text.set_text(message)
        self.add_mode_button.label.set_text(
            (
                f"Stop adding {self.active_side.upper()} points"
                if self.add_point_mode
                else f"Start adding {self.active_side.upper()} points"
            )
        )
        self.add_mode_button.ax.set_facecolor(
            "#ffd59a" if self.add_point_mode else "#f0f0f0"
        )
        self.insert_point_button.label.set_text(
            "Cancel insertion"
            if self.insert_point_mode
            else "Insert after selected"
        )
        self.insert_point_button.ax.set_facecolor(
            "#c8e6c9" if self.insert_point_mode else "#f0f0f0"
        )
        if self._pending_incomplete_save_signature is None:
            self.save_button.label.set_text("Save blocks (Ctrl+S)")
            self.save_button.ax.set_facecolor("#f0f0f0")
        else:
            self.save_button.label.set_text("Confirm incomplete save")
            self.save_button.ax.set_facecolor("#ffcc80")
        manager = getattr(self.figure.canvas, "manager", None)
        if manager is not None:
            manager.set_window_title(
                f"JOREK patch builder — active block {active_description}"
            )
        self._update_block_artists()
        self.figure.canvas.draw_idle()

    def show(self) -> None:
        """Enter Matplotlib's GUI event loop."""

        self._plt.show()


def run_interactive_builder(
    grid_segments: Sequence[np.ndarray],
    limiter: np.ndarray,
    snap_candidates: SnapCandidates,
    default_n_ext: int,
    snap_pixels: float,
    output_path: Path | str,
    initial_blocks: Sequence[Block] | None = None,
) -> PatchBuilderApp:
    """Create, display, and return the interactive patch-builder application."""

    application = PatchBuilderApp(
        grid_segments=grid_segments,
        limiter=limiter,
        snap_candidates=snap_candidates,
        default_n_ext=default_n_ext,
        snap_pixels=snap_pixels,
        output_path=output_path,
        initial_blocks=initial_blocks,
    )
    application.show()
    return application


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise argparse.ArgumentTypeError("must be a positive finite number")
    return parsed


def build_argument_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser without importing Matplotlib."""

    parser = argparse.ArgumentParser(
        description=(
            "Display a JOREK grid and limiter in preparation for interactively "
            "defining extension blocks."
        )
    )
    parser.add_argument(
        "--grid",
        type=Path,
        default=DEFAULT_GRID_PATH,
        help=f"blank-line-separated grid file (default: {DEFAULT_GRID_PATH})",
    )
    parser.add_argument(
        "--limiter",
        type=Path,
        default=DEFAULT_LIMITER_PATH,
        help=f"ordered two-column limiter file (default: {DEFAULT_LIMITER_PATH})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"JOREK block-output file (default: {DEFAULT_OUTPUT_PATH})",
    )
    parser.add_argument(
        "--load",
        type=Path,
        default=None,
        help="existing JOREK block file to load and resume editing",
    )
    parser.add_argument(
        "--snap-pixels",
        type=_positive_float,
        default=12.0,
        help="screen-space snapping radius in pixels (default: 12)",
    )
    parser.add_argument(
        "--default-n-ext",
        type=_positive_int,
        default=30,
        help="default radial resolution for new blocks (default: 30)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Load the configured geometry and open the interactive Matplotlib view."""

    arguments = build_argument_parser().parse_args(argv)

    try:
        grid_segments = read_segment_file(arguments.grid)
        limiter = read_limiter(arguments.limiter)
        snap_candidates = build_snap_candidates(grid_segments, limiter)
        initial_blocks = (
            read_jorek_blocks(arguments.load) if arguments.load is not None else []
        )
        run_interactive_builder(
            grid_segments=grid_segments,
            limiter=limiter,
            snap_candidates=snap_candidates,
            default_n_ext=arguments.default_n_ext,
            snap_pixels=arguments.snap_pixels,
            output_path=arguments.output,
            initial_blocks=initial_blocks,
        )
    except (OSError, DataFormatError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
