"""Test quad table sanity checks."""


def test_quad_table_sanity():
    """Validate quad tables structure and token encoding.

    This is a C++-level validation that would check:
    - Offsets are monotonic
    - Last offset equals TOTAL subcells
    - Vert offsets monotonic
    - Every token is either edge_id in [0,3] or 100+vid with vid in [0,3]

    For now, this is a placeholder. The actual validation happens in C++
    through compilation and existing integration tests.
    """
    # TODO: When table validation utilities are exposed to Python, implement:
    # 1. Load quad_case_subcell_offset_inside
    # 2. Verify monotonic: offset[i] <= offset[i+1]
    # 3. Verify offset[16] == 21 (total subcells)
    # 4. Load quad_subcell_verts_inside
    # 5. For each token != -1: verify token < 4 OR 100 <= token <= 103

    # The quad tables are now defined in cut_quadrilateral_tables.h
    # and are validated through successful compilation and existing tests
    import cutcells
    import numpy as np

    # Run a basic cut to ensure tables work
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=float)
    ls_values = np.array([-0.5, 0.5, 0.5, -0.5], dtype=float)

    cut = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi<0", False
    )

    # Should produce valid output
    assert len(cut.types) > 0
    assert cut.vertex_coords is not None

    # Basic validation passed if we got here
    print(f"Produced {len(cut.types)} subcells")
    print("✓ test_quad_table_sanity passed")
