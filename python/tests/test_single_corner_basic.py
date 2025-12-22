import pytest

pytest.importorskip(
    "tools.cutcells_tables",
    reason="Legacy table-tooling package is not part of this repo layout.",
)

from tools.cutcells_tables.build.hex_variant import build_hex_variant


def test_single_inside_corner_triangle():
    mask = 1  # only corner 0 inside
    vp = build_hex_variant(mask, 0)
    assert len(vp.iso_loops) == 1, f"expected 1 iso loop, got {len(vp.iso_loops)}"
    L = vp.iso_loops[0]
    assert len(L.verts) == 3, f"expected triangle, got {L.verts}"
    assert all(v >= 8 for v in L.verts), "loop should use edge vertices only"
    assert len(set(L.verts)) == 3, "triangle vertices must be distinct"
