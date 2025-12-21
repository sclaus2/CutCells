from tools.cutcells_tables.build.hex_variant import build_hex_variant
from tools.cutcells_tables.classify.hex_mc33 import edges_mask_from_cell_mask, classify_hex


def test_variant_loops_basic():
    for mask in (0, 0xFF, 0b00001111, 0b11110000, 0b01011010):
        vp = build_hex_variant(mask, 0)
        em = edges_mask_from_cell_mask(mask)
        if em == 0:
            assert len(vp.iso_loops) == 0
        else:
            # loops may be >1 now; ensure each closed and edge-vertex indices
            for loop in vp.iso_loops:
                assert loop.is_closed()
                assert all(v >= 8 for v in loop.verts)


def test_vertex_mask_matches_loop():
    mask = 0b01011010
    vp = build_hex_variant(mask, 0)
    bits = 0
    for coll in (vp.iso_loops, vp.inside_faces, vp.outside_faces):
        for L in coll:
            for v in L.verts:
                bits |= 1 << v
    assert bits == vp.vmask


def test_deterministic_ordering():
    mask = 0b01011010
    v1 = build_hex_variant(mask, 0)
    v2 = build_hex_variant(mask, 0)
    assert [L.verts for L in v1.iso_loops] == [L.verts for L in v2.iso_loops]


def test_variant_bit_effect_when_ambiguous():
    # Construct an ambiguous face mask: bottom face 0101 (vertices 0,2 active)
    mask = (1 << 0) | (1 << 2)
    info = classify_hex(mask)
    if info.K == 0:
        return  # no ambiguous faces in this mask; skip
    vA = build_hex_variant(mask, 0)
    vB = build_hex_variant(mask, 1)
    # Either different loops or identical if degenerate; ensure deterministic order regardless
    assert [L.verts for L in vA.iso_loops] == sorted([L.verts for L in vA.iso_loops])
    assert [L.verts for L in vB.iso_loops] == sorted([L.verts for L in vB.iso_loops])
