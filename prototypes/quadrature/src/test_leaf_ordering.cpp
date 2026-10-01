// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Writes one Lagrange hexahedron and one Lagrange quadrilateral per degree 1..4
// whose node coordinates equal their parametric coordinates (i/p, j/p, k/p).
// tools/check_leaves.py reads them with VTK and compares each node with VTK's own
// parametric coordinates for that local index, which validates the node order.

#include <cstdio>
#include <string>

#include "leaf_mesh.h"

using namespace cutcells::proto;

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : ".";
    for (int p = 1; p <= 4; ++p)
    {
        LeafMesh hex, quad;
        const int p1 = p + 1;
        std::vector<std::int32_t> conn(p1 * p1 * p1);
        for (int k = 0; k < p1; ++k)
            for (int j = 0; j < p1; ++j)
                for (int i = 0; i < p1; ++i)
                {
                    const int node = hex.n_points();
                    hex.points.insert(hex.points.end(), {double(i) / p, double(j) / p, double(k) / p});
                    conn[vtk_lagrange_hex_index(i, j, k, p)] = node;
                }
        hex.connectivity = conn;
        hex.offsets.push_back(static_cast<std::int32_t>(conn.size()));
        hex.types.push_back(vtk_lagrange_hexahedron);
        hex.parent.push_back(0);
        hex.degree.push_back(p);
        write_leaf_mesh(dir + "/lagrange_hex_p" + std::to_string(p) + ".vtu", hex);

        conn.assign(p1 * p1, 0);
        for (int j = 0; j < p1; ++j)
            for (int i = 0; i < p1; ++i)
            {
                const int node = quad.n_points();
                quad.points.insert(quad.points.end(), {double(i) / p, double(j) / p, 0.0});
                conn[vtk_lagrange_quad_index(i, j, p)] = node;
            }
        quad.connectivity = conn;
        quad.offsets.push_back(static_cast<std::int32_t>(conn.size()));
        quad.types.push_back(vtk_lagrange_quadrilateral);
        quad.parent.push_back(0);
        quad.degree.push_back(p);
        write_leaf_mesh(dir + "/lagrange_quad_p" + std::to_string(p) + ".vtu", quad);
    }
    std::printf("wrote Lagrange hexahedra and quadrilaterals, degrees 1 to 4, to %s\n", dir.c_str());
    return 0;
}
