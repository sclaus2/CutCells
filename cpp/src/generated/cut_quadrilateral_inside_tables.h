#ifndef CUT_CELLS_QUADRILATERAL_INSIDE_TABLES_H
#define CUT_CELLS_QUADRILATERAL_INSIDE_TABLES_H

#include "cell_types.h"
#include <cstdint>

namespace cutcells::cell::quadrilateral {

// Topology
constexpr int edges[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};

// Intersected edges per case (1 = intersected, 0 = not intersected)
constexpr int intersected_edges[16][4] = {
    { 0, 0, 0, 0 },
    { 1, 0, 0, 1 },
    { 1, 1, 0, 0 },
    { 0, 1, 0, 1 },
    { 0, 1, 1, 0 },
    { 1, 1, 1, 1 },
    { 1, 0, 1, 0 },
    { 0, 0, 1, 1 },
    { 0, 0, 1, 1 },
    { 1, 0, 1, 0 },
    { 1, 1, 1, 1 },
    { 0, 1, 1, 0 },
    { 0, 1, 0, 1 },
    { 1, 1, 0, 0 },
    { 1, 0, 0, 1 },
    { 0, 0, 0, 0 }
};

// Ambiguity metadata (opposite-corner cases)
constexpr uint8_t case_is_ambiguous_tbl[16] = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
constexpr int8_t amb_case_id[16] = { -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1 };

// Number of subcells produced for each case (inside volume)
constexpr int num_subcells_inside[16] = { 0, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1, 2, 2, 1 };

// Offset into subcell array for each case
constexpr int case_subcell_offset_inside[17] = { 0, 0, 1, 2, 3, 4, 6, 7, 9, 10, 11, 13, 15, 16, 18, 20, 21 };

// Cell types for inside subcells
constexpr type subcell_type_inside[25] = { cell::type::triangle, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::triangle, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::triangle, cell::type::triangle, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral };

// Subcell vertices (max 4 vertices per subcell, -1 padding)
// Tokens: <100 = edge id, >=100 = 100+vertex_id
constexpr int subcell_verts_inside[25][4] = {
    { 100, 0, 3, -1 },
    { 0, 101, 1, -1 },
    { 100, 101, 1, 3 },
    { 1, 102, 2, -1 },
    { 100, 0, 3, -1 },
    { 102, 1, 2, -1 },
    { 0, 101, 102, 2 },
    { 3, 100, 101, -1 },
    { 3, 101, 102, 2 },
    { 2, 103, 3, -1 },
    { 100, 0, 2, 103 },
    { 101, 1, 0, -1 },
    { 103, 3, 2, -1 },
    { 2, 103, 100, -1 },
    { 2, 100, 101, 1 },
    { 1, 102, 103, 3 },
    { 1, 102, 103, -1 },
    { 1, 103, 100, 0 },
    { 0, 101, 102, -1 },
    { 0, 102, 103, 3 },
    { 100, 101, 102, 103 },
    { 100, 0, 1, 102 },
    { 100, 102, 2, 3 },
    { 101, 1, 2, 103 },
    { 101, 103, 3, 0 }
};

// Ambiguity override ranges: [begin0,end0,begin1,end1] per ambiguous case
constexpr int amb_range_inside[8] = { 21, 23, 4, 6, 11, 13, 23, 25 };

} // namespace cutcells::cell::quadrilateral

#endif // CUT_CELLS_QUADRILATERAL_INSIDE_TABLES_H