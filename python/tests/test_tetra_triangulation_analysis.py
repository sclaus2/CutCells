import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cutcells

pytestmark = pytest.mark.skipif(
    not hasattr(cutcells, "analyze_all_cases"),
    reason="triangulation_analysis helpers are not packaged",
)

CellType = cutcells.CellType
analyze_all_cases = getattr(cutcells, "analyze_all_cases", None)
summarize_analysis = getattr(cutcells, "summarize_analysis", None)


def test_tetra_parent_summary_excludes_zero_only_edges_from_interior_root_count():
    records, summary = analyze_all_cases(CellType.tetrahedron)

    assert summary["cases"] == 42
    assert summary["total_new_edges"] == 218
    assert summary["new_edges_with_interior_roots"] == 0
    assert summary["new_edges_with_zero_tag"] == 18
    assert summary["new_edges_with_zero_endpoint"] == 122
    assert summary["total_new_simplices"] == 116
    assert summary["negative_new_simplices"] == 0
    assert summary["degenerate_new_simplices"] == 0


def test_tetra_parent_interface_quads_have_no_interior_root_edges_or_negative_triangles():
    records, _ = analyze_all_cases(CellType.tetrahedron)
    interface_summary = summarize_analysis([record for record in records if record["expr"] == "phi=0"])

    assert interface_summary["cases"] == 14
    assert interface_summary["source_type_counts"] == {"quadrilateral": 12}
    assert interface_summary["new_edges_with_interior_roots"] == 0
    assert interface_summary["new_edges_with_zero_tag"] == 6
    assert interface_summary["new_edges_with_zero_endpoint"] == 6
    assert interface_summary["negative_new_simplices"] == 0
    assert interface_summary["degenerate_new_simplices"] == 0


def test_tetra_parent_prism_volume_cases_have_only_endpoint_roots_and_positive_tets():
    records, _ = analyze_all_cases(CellType.tetrahedron)

    negative_summary = summarize_analysis([record for record in records if record["expr"] == "phi<0"])
    positive_summary = summarize_analysis([record for record in records if record["expr"] == "phi>0"])

    for summary in (negative_summary, positive_summary):
        assert summary["cases"] == 14
        assert summary["source_type_counts"] == {"prism": 52}
        assert summary["new_edges_with_interior_roots"] == 0
        assert summary["new_edges_with_zero_tag"] == 6
        assert summary["new_edges_with_zero_endpoint"] == 58
        assert summary["negative_new_simplices"] == 0
        assert summary["degenerate_new_simplices"] == 0
