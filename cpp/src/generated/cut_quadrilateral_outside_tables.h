#ifndef CUT_CELLS_QUADRILATERAL_OUTSIDE_TABLES_H
#define CUT_CELLS_QUADRILATERAL_OUTSIDE_TABLES_H

#include "cell_types.h"

namespace cutcells::cell::quadrilateral {

// Number of subcells produced for each case (outside volume)
constexpr int num_subcells_outside[16] = { 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1, 1, 1, 0 };

// Offset into subcell array for each case
constexpr int case_subcell_offset_outside[17] = { 0, 1, 3, 5, 6, 8, 10, 11, 12, 14, 15, 17, 18, 19, 20, 21, 21 };

// Cell types for outside subcells
constexpr type subcell_type_outside[21] = { cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::triangle, cell::type::triangle, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::quadrilateral, cell::type::triangle, cell::type::quadrilateral, cell::type::triangle, cell::type::triangle };

// Subcell vertices (max 4 vertices per subcell, -1 padding)
// Tokens: <100 = edge id, >=100 = 100+vertex_id
constexpr int subcell_verts_outside[21][4] = {
    { 100, 101, 102, 103 },
    { 0, 101, 102, -1 },
    { 0, 102, 103, 3 },
    { 1, 102, 103, -1 },
    { 1, 103, 100, 0 },
    { 1, 102, 103, 3 },
    { 2, 103, 100, -1 },
    { 2, 100, 101, 1 },
    { 101, 1, 2, 103 },
    { 103, 3, 0, 101 },
    { 100, 0, 2, 103 },
    { 2, 103, 3, -1 },
    { 3, 100, 101, -1 },
    { 3, 101, 102, 2 },
    { 0, 101, 102, 2 },
    { 100, 0, 1, 102 },
    { 102, 2, 3, 100 },
    { 1, 102, 2, -1 },
    { 100, 101, 1, 3 },
    { 0, 101, 1, -1 },
    { 100, 0, 3, -1 }
};

} // namespace cutcells::cell::quadrilateral

#endif // CUT_CELLS_QUADRILATERAL_OUTSIDE_TABLES_H