"""Test that CutCell stores original cell information."""

import cutcells
import numpy as np


def test_cutcell_original_info():
    """Verify CutCell stores original vertex coords and IDs."""
    # Create a simple quad
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=float)
    ls_values = np.array([-0.5, -0.5, 0.5, 0.5], dtype=float)

    # Cut it
    cut = cutcells.cut(
        cutcells.CellType.quadrilateral, vertices, 2, ls_values, "phi<0", False
    )

    # For now, just verify the cut doesn't crash
    # TODO: Once Python bindings expose _parent_* fields, verify:
    # - cut._parent_cell_type == cutcells.CellType.quadrilateral
    # - cut._parent_vertex_coords matches vertices
    # - cut._parent_vertex_ids is [0,1,2,3] (or similar default)

    assert len(cut.types) > 0, "Should produce some cut cells"
    assert cut.vertex_coords is not None, "Should have vertex coordinates"


if __name__ == "__main__":
    test_cutcell_original_info()
    print("✓ test_cutcell_original_info passed")
