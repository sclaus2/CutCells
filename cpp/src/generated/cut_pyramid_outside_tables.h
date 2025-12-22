#ifndef CUT_CELLS_PYRAMID_OUTSIDE_TABLES_H
#define CUT_CELLS_PYRAMID_OUTSIDE_TABLES_H

#include "cell_types.h"

namespace cutcells::cell::pyramid {

// Number of subcells produced for each case (outside volume)
constexpr int num_subcells_outside[32] = { 1, 2, 2, 6, 2, 2, 6, 2, 2, 6, 2, 2, 6, 2, 2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1, 1, 1, 0 };

// Offset into subcell array for each case
constexpr int case_subcell_offset_outside[33] = { 0, 1, 3, 5, 11, 13, 15, 21, 23, 25, 31, 33, 35, 41, 43, 45, 46, 47, 49, 51, 52, 54, 56, 57, 58, 60, 61, 63, 64, 65, 66, 67, 67 };

// Cell types for outside subcells
constexpr type subcell_type_outside[67] = { cell::type::pyramid, cell::type::prism, cell::type::tetrahedron, cell::type::prism, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::tetrahedron, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::prism, cell::type::tetrahedron, cell::type::prism, cell::type::prism, cell::type::pyramid, cell::type::pyramid, cell::type::tetrahedron, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::prism, cell::type::tetrahedron, cell::type::prism, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::tetrahedron, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::tetrahedron, cell::type::tetrahedron, cell::type::pyramid, cell::type::pyramid, cell::type::prism, cell::type::tetrahedron, cell::type::prism, cell::type::tetrahedron, cell::type::pyramid, cell::type::hexahedron, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::tetrahedron, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::prism, cell::type::tetrahedron, cell::type::prism, cell::type::tetrahedron, cell::type::tetrahedron };

// Subcell vertices (max 8 vertices per subcell, -1 padding)
// Tokens: <100 = edge id, 100..199 = 100+vertex_id, >=200 = 200+special_point_id
constexpr int subcell_verts_outside[67][8] = {
    { 100, 101, 102, 103, 104, -1, -1, -1 },
    { 101, 103, 104, 0, 3, 4, -1, -1 },
    { 101, 102, 103, 104, -1, -1, -1, -1 },
    { 102, 100, 104, 1, 0, 5, -1, -1 },
    { 102, 103, 100, 104, -1, -1, -1, -1 },
    { 4, 5, 1, 3, 200, -1, -1, -1 },
    { 1, 102, 103, 3, 200, -1, -1, -1 },
    { 102, 104, 103, 200, -1, -1, -1, -1 },
    { 5, 4, 104, 200, -1, -1, -1, -1 },
    { 5, 104, 102, 1, 200, -1, -1, -1 },
    { 104, 4, 3, 103, 200, -1, -1, -1 },
    { 103, 101, 104, 2, 1, 6, -1, -1 },
    { 103, 100, 101, 104, -1, -1, -1, -1 },
    { 3, 0, 4, 103, 101, 104, -1, -1 },
    { 1, 2, 6, 101, 103, 104, -1, -1 },
    { 5, 6, 2, 0, 200, -1, -1, -1 },
    { 2, 103, 100, 0, 200, -1, -1, -1 },
    { 103, 104, 100, 200, -1, -1, -1, -1 },
    { 6, 5, 104, 200, -1, -1, -1, -1 },
    { 6, 104, 103, 2, 200, -1, -1, -1 },
    { 104, 5, 0, 100, 200, -1, -1, -1 },
    { 104, 4, 6, 103, 3, 2, -1, -1 },
    { 4, 5, 6, 104, -1, -1, -1, -1 },
    { 100, 102, 104, 3, 2, 7, -1, -1 },
    { 100, 101, 102, 104, -1, -1, -1, -1 },
    { 7, 4, 0, 2, 200, -1, -1, -1 },
    { 0, 101, 102, 2, 200, -1, -1, -1 },
    { 101, 104, 102, 200, -1, -1, -1, -1 },
    { 4, 7, 104, 200, -1, -1, -1, -1 },
    { 4, 104, 101, 0, 200, -1, -1, -1 },
    { 104, 7, 2, 102, 200, -1, -1, -1 },
    { 0, 1, 5, 100, 102, 104, -1, -1 },
    { 2, 3, 7, 102, 100, 104, -1, -1 },
    { 104, 7, 5, 102, 2, 1, -1, -1 },
    { 7, 4, 5, 104, -1, -1, -1, -1 },
    { 6, 7, 3, 1, 200, -1, -1, -1 },
    { 3, 100, 101, 1, 200, -1, -1, -1 },
    { 100, 104, 101, 200, -1, -1, -1, -1 },
    { 7, 6, 104, 200, -1, -1, -1, -1 },
    { 7, 104, 100, 3, 200, -1, -1, -1 },
    { 104, 6, 1, 101, 200, -1, -1, -1 },
    { 104, 6, 4, 101, 1, 0, -1, -1 },
    { 6, 7, 4, 104, -1, -1, -1, -1 },
    { 104, 5, 7, 100, 0, 3, -1, -1 },
    { 5, 6, 7, 104, -1, -1, -1, -1 },
    { 4, 5, 6, 7, 104, -1, -1, -1 },
    { 100, 101, 102, 103, 4, 5, 6, 7 },
    { 5, 6, 7, 101, 102, 103, -1, -1 },
    { 103, 3, 7, 101, 0, 5, -1, -1 },
    { 6, 7, 4, 102, 103, 100, -1, -1 },
    { 100, 0, 4, 102, 1, 6, -1, -1 },
    { 3, 7, 103, 1, 6, 102, -1, -1 },
    { 7, 4, 5, 103, 100, 101, -1, -1 },
    { 101, 1, 5, 103, 2, 7, -1, -1 },
    { 103, 3, 7, 101, 0, 5, -1, -1 },
    { 2, 103, 7, 1, 101, 5, -1, -1 },
    { 0, 4, 100, 2, 7, 103, -1, -1 },
    { 103, 3, 2, 7, -1, -1, -1, -1 },
    { 4, 5, 6, 100, 101, 102, -1, -1 },
    { 102, 2, 6, 100, 3, 4, -1, -1 },
    { 2, 6, 102, 0, 5, 101, -1, -1 },
    { 100, 0, 4, 102, 1, 6, -1, -1 },
    { 3, 100, 4, 2, 102, 6, -1, -1 },
    { 102, 2, 1, 6, -1, -1, -1, -1 },
    { 1, 5, 101, 3, 4, 100, -1, -1 },
    { 101, 1, 0, 5, -1, -1, -1, -1 },
    { 100, 0, 3, 4, -1, -1, -1, -1 }
};

// Special points (e.g. centroid): per-case definition streams
constexpr int special_point_count_outside[32] = { 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
constexpr int special_point_offset_outside[33] = { 0, 0, 0, 0, 8, 8, 8, 16, 16, 16, 24, 24, 24, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 };
constexpr int special_point_data_outside[32] = { 7, 5, 4, 1, 3, 102, 103, 104, 7, 6, 5, 2, 0, 103, 100, 104, 7, 4, 7, 0, 2, 101, 102, 104, 7, 7, 6, 3, 1, 100, 101, 104 };

} // namespace cutcells::cell::pyramid

#endif // CUT_CELLS_PYRAMID_OUTSIDE_TABLES_H