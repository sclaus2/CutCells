from tools.cutcells_tables.export.hex_polyhedron import build_hex_polyhedra


def test_single_corner_tetra_inside_polyhedron():
    # Mask with only corner 0 inside
    m = 1
    data = build_hex_polyhedra(m, 0)
    assert data.inside is not None, "Inside tetra should be watertight"
    # Tetra has 4 faces (iso + 3 corner faces)
    assert data.inside.n_faces == 4
    assert sum(data.inside.face_sizes) == len(data.inside.faces_stream)


def test_full_or_empty_masks_no_polyhedra():
    # No iso surface => boundaries not mixed (pure hex/no cut) -> cannot form cut polyhedra
    for m in (0, 0xFF):
        data = build_hex_polyhedra(m, 0)
        # Either inside or outside degenerates; both can't be polyhedron with iso surface absent
        assert data.inside is None or data.outside is None
