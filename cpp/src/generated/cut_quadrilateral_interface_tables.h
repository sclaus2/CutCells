#ifndef CUT_CELLS_QUADRILATERAL_INTERFACE_TABLES_H
#define CUT_CELLS_QUADRILATERAL_INTERFACE_TABLES_H

#include "cell_types.h"

namespace cutcells::cell::quadrilateral {

// Number of subcells produced for each case (interface volume)
constexpr int num_subcells_interface[16] = { 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0 };

// Offset into subcell array for each case
constexpr int case_subcell_offset_interface[17] = { 0, 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 16 };

// Cell types for interface subcells
constexpr type subcell_type_interface[16] = { cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval, cell::type::interval };

// Subcell vertices (max 2 vertices per subcell, -1 padding)
// Tokens: <100 = edge id, >=100 = 100+vertex_id
constexpr int subcell_verts_interface[16][2] = {
    { 0, 3 },
    { 0, 1 },
    { 1, 3 },
    { 1, 2 },
    { 0, 3 },
    { 1, 2 },
    { 0, 2 },
    { 2, 3 },
    { 2, 3 },
    { 0, 2 },
    { 0, 1 },
    { 2, 3 },
    { 1, 2 },
    { 1, 3 },
    { 0, 1 },
    { 0, 3 }
};

} // namespace cutcells::cell::quadrilateral

#endif // CUT_CELLS_QUADRILATERAL_INTERFACE_TABLES_H