"""Test CutContext for edge-owned intersection vertices."""

import numpy as np


def test_cutcontext_edge_reuse():
    """Test that CutContext correctly deduplicates edge intersections.

    This is a placeholder test until CutContext is exposed to Python.
    The actual functionality will be tested in C++ or through integration tests.
    """
    # TODO: Once CutContext is exposed to Python bindings, implement:
    # 1. Create CutContext with num_original_vertices=4, gdim=2
    # 2. Insert edge (0,1) with phi values → get vid1
    # 3. Insert same edge again → get same vid (reuse)
    # 4. Insert edge with endpoint phi==0 → returns original vertex id
    # 5. Verify ip_coords has exactly one entry per unique intersection

    # For now, just verify the concept works through existing cut operations
    # which should internally use similar logic
    import cutcells

    # Create a quad and cut it
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=float)
    ls_values = np.array([-0.5, 0.5, 0.5, -0.5], dtype=float)

    cut = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi<0", False
    )

    # Basic sanity: cut should produce cells
    assert len(cut.types) > 0
    assert cut.vertex_coords is not None

    # TODO: When CutContext is integrated, verify:
    # - Same edge intersections produce same vertex ids
    # - No duplicate vertices for shared edges


if __name__ == "__main__":
    test_cutcontext_edge_reuse()
    print("✓ test_cutcontext_edge_reuse passed (placeholder)")
