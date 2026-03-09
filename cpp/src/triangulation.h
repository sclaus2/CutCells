// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cell_types.h"

#include <numeric>
#include <stdexcept>
#include <vector>

namespace cutcells::cell
{
    /// @brief Decompose a non-simplex cell into simplices (triangles or tetrahedra).
    ///
    /// Each entry in @p tris is a list of vertex indices (into the parent vertices[])
    /// array) that form one simplex.
    ///
    /// Supported cell types and their outputs:
    ///   quadrilateral  → 2 triangles  (3 vertices each)
    ///   hexahedron     → 6 tetrahedra (4 vertices each)
    ///   prism          → 3 tetrahedra (4 vertices each)
    ///   pyramid        → 2 tetrahedra (4 vertices each)
    ///
    /// VTK vertex ordering is assumed for all cell types.
    inline void triangulation(const type cell_type, int* vertices, std::vector<std::vector<int>>& tris)
    {
        switch(cell_type)
        {
            case type::quadrilateral:
                tris.resize(2, std::vector<int>(3));
                tris = {{vertices[0],vertices[1],vertices[2]}, {vertices[0],vertices[2],vertices[3]}};
                break;

            case type::hexahedron:
                // VTK hex: v0=(0,0,0) v1=(1,0,0) v2=(1,1,0) v3=(0,1,0)
                //          v4=(0,0,1) v5=(1,0,1) v6=(1,1,1) v7=(0,1,1)
                //
                // Kuhn triangulation (6 tets, fan through the v0→v6 diagonal).
                // Each tet has |det J| = 1 on the unit cube → total volume = 6 × 1/6 = 1 ✓
                tris.resize(6, std::vector<int>(4));
                tris = {{vertices[0],vertices[1],vertices[2],vertices[6]},
                        {vertices[0],vertices[1],vertices[5],vertices[6]},
                        {vertices[0],vertices[3],vertices[2],vertices[6]},
                        {vertices[0],vertices[3],vertices[7],vertices[6]},
                        {vertices[0],vertices[4],vertices[5],vertices[6]},
                        {vertices[0],vertices[4],vertices[7],vertices[6]}};
                break;

            case type::prism:
                // VTK wedge: v0,v1,v2 bottom △, v3,v4,v5 top △
                // Tetrahedron 0: { v0, v2, v1, v3 }
                // Tetrahedron 1: { v1, v3, v5, v4 }
                // Tetrahedron 2: { v1, v2, v5, v3 }
                tris.resize(3, std::vector<int>(4));
                tris = {{vertices[0],vertices[2],vertices[1],vertices[3]},
                        {vertices[1],vertices[3],vertices[5],vertices[4]},
                        {vertices[1],vertices[2],vertices[5],vertices[3]}};
                break;

            case type::pyramid:
                // VTK pyramid: v0=(0,0,0) v1=(1,0,0) v2=(1,1,0) v3=(0,1,0) v4=apex
                // Tetrahedron 0: { v0, v1, v3, v4 }
                // Tetrahedron 1: { v1, v2, v3, v4 }
                tris.resize(2, std::vector<int>(4));
                tris = {{vertices[0],vertices[1],vertices[3],vertices[4]},
                        {vertices[1],vertices[2],vertices[3],vertices[4]}};
                break;

            default:
                throw std::invalid_argument("triangulation not implemented for given cell type");
                break;
        }
    };

    /// @brief Canonical reference vertices for a cell type (VTK ordering, flat).
    ///
    /// Returns a flat vector of nv * tdim values representing the reference-element
    /// vertex coordinates for the given cell type.  Coordinates are in the
    /// topological-dimension space (i.e. tdim values per vertex, NOT padded to gdim).
    ///
    /// VTK canonical positions used (these make the affine Jacobian = I for the
    /// unit reference cell):
    ///   interval      v0=(0)   v1=(1)
    ///   triangle      v0=(0,0) v1=(1,0) v2=(0,1)
    ///   quadrilateral v0=(0,0) v1=(1,0) v2=(1,1) v3=(0,1)
    ///   tetrahedron   v0=(0,0,0) v1=(1,0,0) v2=(0,1,0) v3=(0,0,1)
    ///   hexahedron    v0=(0,0,0) v1=(1,0,0) v2=(1,1,0) v3=(0,1,0)
    ///                 v4=(0,0,1) v5=(1,0,1) v6=(1,1,1) v7=(0,1,1)
    ///   prism         v0=(0,0,0) v1=(1,0,0) v2=(0,1,0)
    ///                 v3=(0,0,1) v4=(1,0,1) v5=(0,1,1)
    ///   pyramid       v0=(0,0,0) v1=(1,0,0) v2=(1,1,0) v3=(0,1,0) v4=(0,0,1)
    template <std::floating_point T>
    inline std::vector<T> canonical_vertices(type cell_type)
    {
        switch (cell_type)
        {
            case type::interval:
                return {T(0), T(1)};
            case type::triangle:
                return {T(0),T(0),  T(1),T(0),  T(0),T(1)};
            case type::quadrilateral:
                return {T(0),T(0),  T(1),T(0),  T(1),T(1),  T(0),T(1)};
            case type::tetrahedron:
                return {T(0),T(0),T(0),  T(1),T(0),T(0),  T(0),T(1),T(0),  T(0),T(0),T(1)};
            case type::hexahedron:
                return {T(0),T(0),T(0),  T(1),T(0),T(0),  T(1),T(1),T(0),  T(0),T(1),T(0),
                        T(0),T(0),T(1),  T(1),T(0),T(1),  T(1),T(1),T(1),  T(0),T(1),T(1)};
            case type::prism:
                return {T(0),T(0),T(0),  T(1),T(0),T(0),  T(0),T(1),T(0),
                        T(0),T(0),T(1),  T(1),T(0),T(1),  T(0),T(1),T(1)};
            case type::pyramid:
                return {T(0),T(0),T(0),  T(1),T(0),T(0),  T(1),T(1),T(0),  T(0),T(1),T(0),  T(0),T(0),T(1)};
            default:
                throw std::invalid_argument("canonical_vertices: unsupported cell type");
        }
    }

} // namespace cutcells::cell
