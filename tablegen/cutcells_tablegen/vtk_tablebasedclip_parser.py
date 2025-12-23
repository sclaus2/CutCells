"""Parse VTK TableBasedClip case tables.

This module is intentionally lightweight: it only extracts the VTK_HEXAHEDRON
case stream from VTK's `Filters/General/vtkTableBasedClipCases.h`.

We do not depend on VTK headers/libraries; instead we parse the published
header text (either from a local file path or downloaded URL).

The extracted structure is converted into CutCells tablegen `ClipCase` objects.

License/provenance:
    The case tables parsed by this module are derived from VTK's
    `Filters/General/vtkTableBasedClipCases.h`.
    VTK is licensed under the BSD 3-Clause license; see
    `third_party/VTK-Copyright.txt` in this repository.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import os
import re
import urllib.request


_VTK_DEFAULT_URL = "https://raw.githubusercontent.com/Kitware/VTK/{ref}/Filters/General/vtkTableBasedClipCases.h"


@dataclass(frozen=True)
class VtkCaseRecord:
    case_id: int
    items: Sequence[Sequence[str]]
    # Each item is a token list: [ST_*, n, p0, p1, ...]


def _strip_comments(text: str) -> str:
    # Remove /* ... */ comments and // line comments
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)
    text = re.sub(r"//.*?$", " ", text, flags=re.MULTILINE)
    return text


def _read_header_text(path: str | None, ref: str) -> str:
    if path:
        with open(path, "r", encoding="utf-8") as f:
            return f.read()

    url = _VTK_DEFAULT_URL.format(ref=ref)
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


def _find_hex_block(text: str) -> str:
    # VTK's header contains multiple hexahedron-related tables *and* several
    # earlier occurrences of the string "VTK_HEXAHEDRON" (e.g. in topology
    # comments). We must select the actual hexahedron case stream inside the
    # big CellCases array.
    #
    # Robust approach:
    #  - search for the *case table marker* line: "// VTK_HEXAHEDRON"
    #  - require it to be followed by the 256-case stream starting at
    #    "/* case 0 */ 0,"
    #  - validate a couple of hexahedron-specific signatures to avoid matching
    #    other tables.
    marker_pat = re.compile(r"\n\s*//\s*VTK_HEXAHEDRON\s*\n")

    # Signature checks within the hexahedron volume clip case stream.
    case0_pat = re.compile(r"/\*\s*case\s*0\s*\*/\s*(\d+)\s*,")
    case1_sig = re.compile(
        r"/\*\s*case\s*1\s*\*/\s*1\s*,\s*ST_TET\s*,\s*4\s*,\s*P0\s*,\s*EA\s*,\s*ED\s*,\s*EI\s*,"
    )
    # A second signature that differs between various hex-related tables.
    # In the VTK_HEXAHEDRON clip table, case 153 is a single hex:
    #   ST_HEX, 8, EA, EE, EG, EC, P0, P4, P7, P3
    case153_sig = re.compile(
        r"/\*\s*case\s*153\s*\*/\s*1\s*,\s*ST_HEX\s*,\s*8\s*,\s*EA\s*,\s*EE\s*,\s*EG\s*,\s*EC\s*,\s*P0\s*,\s*P4\s*,\s*P7\s*,\s*P3\s*,"
    )

    for m in marker_pat.finditer(text):
        # Look ahead; the full 256-case table is large.
        window = text[m.end() : m.end() + 250000]

        m0 = case0_pat.search(window)
        if not m0 or m0.group(1) != "0":
            continue

        # Bind signatures to the specific candidate stream.
        tail = window[m0.start() : m0.start() + 60000]
        if not case1_sig.search(tail):
            continue
        if not case153_sig.search(tail):
            continue

        # Important: later we strip comments before tokenizing; the "// VTK_*"
        # marker would be removed, so we return a slice starting exactly at the
        # numeric case stream (the count following "/* case 0 */").
        start = m.end() + m0.start(1)
        return text[start : start + 250000]

    raise RuntimeError(
        "Could not locate the VTK_HEXAHEDRON clip case stream (expected '/* case 0 */ 0' and known signatures)"
    )


def _tokenize_case_block(block: str) -> List[str]:
    # Keep identifiers like ST_TET, P0, EA, N0 and integers.
    # Strip comments first to avoid commas in comments.
    block = _strip_comments(block)

    # Reduce to comma-separated tokens.
    tokens = [t.strip() for t in block.split(",")]
    return [t for t in tokens if t]


def _parse_cases_from_tokens(tokens: Sequence[str]) -> List[VtkCaseRecord]:
    # The provided token stream is expected to start at the numeric case stream
    # (i.e. the count for case 0). Parse n cases sequentially.
    return _parse_n_cases_from_tokens(tokens, 256)


def _parse_n_cases_from_tokens(
    tokens: Sequence[str], n_cases: int
) -> List[VtkCaseRecord]:
    i = 0
    if i >= len(tokens) or not re.fullmatch(r"\d+", tokens[i]):
        raise RuntimeError("Token stream does not start with a case count")

    cases: List[VtkCaseRecord] = []
    for case_id in range(n_cases):
        if i >= len(tokens):
            raise RuntimeError(f"Unexpected end of tokens while parsing case {case_id}")
        n_items = int(tokens[i])
        i += 1

        items: List[List[str]] = []
        for _ in range(n_items):
            st = tokens[i]
            i += 1
            n = int(tokens[i])
            i += 1
            pts = list(tokens[i : i + n])
            i += n
            items.append([st, str(n), *pts])

        cases.append(VtkCaseRecord(case_id=case_id, items=tuple(items)))

    return cases


def load_hexahedron_cases(
    *,
    vtk_header_path: str | None = None,
    vtk_ref: str | None = None,
) -> List[VtkCaseRecord]:
    """Load VTK_HEXAHEDRON case records.

    If `vtk_header_path` is not provided, the header is downloaded from GitHub.

    Environment variables:
      - CUT_CELLS_VTK_HEADER: local path override
      - CUT_CELLS_VTK_REF: Git ref (default: master)
    """

    if vtk_header_path is None:
        vtk_header_path = os.environ.get("CUT_CELLS_VTK_HEADER")
    if vtk_ref is None:
        vtk_ref = os.environ.get("CUT_CELLS_VTK_REF", "master")

    text = _read_header_text(vtk_header_path, vtk_ref)
    block = _find_hex_block(text)
    tokens = _tokenize_case_block(block)
    return _parse_n_cases_from_tokens(tokens, 256)


def _find_wedge_block(text: str) -> str:
    # Select the VTK_WEDGE case stream inside the CellCases array.
    marker_pat = re.compile(r"\n\s*//\s*VTK_WEDGE\s*\n")

    case0_pat = re.compile(r"/\*\s*case\s*0\s*\*/\s*(\d+)\s*,")
    # Signature anchored to the wedge case table: case 1 is a tet containing P0
    # and the three edge points around vertex 0.
    case1_sig = re.compile(
        r"/\*\s*case\s*1\s*\*/\s*1\s*,\s*ST_TET\s*,\s*4\s*,\s*EG\s*,\s*EA\s*,\s*EC\s*,\s*P0\s*,"
    )
    case63_sig = re.compile(r"/\*\s*case\s*63\s*\*/")

    for m in marker_pat.finditer(text):
        window = text[m.end() : m.end() + 80000]
        m0 = case0_pat.search(window)
        if not m0 or m0.group(1) != "0":
            continue

        # Bind checks to the candidate stream.
        tail = window[m0.start() : m0.start() + 20000]
        if not case1_sig.search(tail):
            continue
        if not case63_sig.search(tail):
            continue

        start = m.end() + m0.start(1)
        return text[start : start + 80000]

    raise RuntimeError(
        "Could not locate the VTK_WEDGE clip case stream (expected '/* case 0 */ 0' and known signatures)"
    )


def load_wedge_cases(
    *,
    vtk_header_path: str | None = None,
    vtk_ref: str | None = None,
) -> List[VtkCaseRecord]:
    """Load VTK_WEDGE (prism) case records (64 cases)."""

    if vtk_header_path is None:
        vtk_header_path = os.environ.get("CUT_CELLS_VTK_HEADER")
    if vtk_ref is None:
        vtk_ref = os.environ.get("CUT_CELLS_VTK_REF", "master")

    text = _read_header_text(vtk_header_path, vtk_ref)
    block = _find_wedge_block(text)
    tokens = _tokenize_case_block(block)
    return _parse_n_cases_from_tokens(tokens, 64)


def _find_pyramid_block(text: str) -> str:
    # Select the VTK_PYRAMID case stream inside the CellCases array.
    marker_pat = re.compile(r"\n\s*//\s*VTK_PYRAMID\s*\n")

    case0_pat = re.compile(r"/\*\s*case\s*0\s*\*/\s*(\d+)\s*,")
    # Pyramid signature: case 1 is a tet containing P0 and the three edge points around vertex 0.
    case1_sig = re.compile(
        r"/\*\s*case\s*1\s*\*/\s*1\s*,\s*ST_TET\s*,\s*4\s*,\s*P0\s*,\s*EA\s*,\s*ED\s*,\s*EE\s*,"
    )
    case31_sig = re.compile(r"/\*\s*case\s*31\s*\*/")

    for m in marker_pat.finditer(text):
        window = text[m.end() : m.end() + 120000]
        m0 = case0_pat.search(window)
        if not m0 or m0.group(1) != "0":
            continue

        tail = window[m0.start() : m0.start() + 30000]
        if not case1_sig.search(tail):
            continue
        if not case31_sig.search(tail):
            continue

        start = m.end() + m0.start(1)
        return text[start : start + 120000]

    raise RuntimeError(
        "Could not locate the VTK_PYRAMID clip case stream (expected '/* case 0 */ 0' and known signatures)"
    )


def load_pyramid_cases(
    *,
    vtk_header_path: str | None = None,
    vtk_ref: str | None = None,
) -> List[VtkCaseRecord]:
    """Load VTK_PYRAMID case records (32 cases)."""

    if vtk_header_path is None:
        vtk_header_path = os.environ.get("CUT_CELLS_VTK_HEADER")
    if vtk_ref is None:
        vtk_ref = os.environ.get("CUT_CELLS_VTK_REF", "master")

    text = _read_header_text(vtk_header_path, vtk_ref)
    block = _find_pyramid_block(text)
    tokens = _tokenize_case_block(block)
    return _parse_n_cases_from_tokens(tokens, 32)
