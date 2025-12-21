// Copyright (c) 2025 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include <unordered_map>
#include <vector>
#include <cstddef>

namespace cutcells
{
    /// Edge key for identifying unique edges (always a < b)
    struct EdgeKey
    {
        int a, b;

        EdgeKey(int v0, int v1) : a(v0 < v1 ? v0 : v1), b(v0 < v1 ? v1 : v0) {}

        bool operator==(const EdgeKey& other) const noexcept
        {
            return a == other.a && b == other.b;
        }
    };

    /// Hash functor for EdgeKey
    struct EdgeKeyHash
    {
        std::size_t operator()(const EdgeKey& key) const noexcept
        {
            // Simple hash combination
            return std::hash<int>()(key.a) ^ (std::hash<int>()(key.b) << 1);
        }
    };

    /// Context for managing edge-owned intersection vertices
    /// Ensures no duplicate intersection vertices across cells
    struct CutContext
    {
        /// Number of original vertices in the context
        int num_original_vertices;

        /// Geometric dimension
        int gdim;

        /// Map from edge to global vertex id
        std::unordered_map<EdgeKey, int, EdgeKeyHash> edge_to_vid;

        /// Packed coordinates of new intersection vertices only (size = num_new_verts * gdim)
        std::vector<double> ip_coords;

        CutContext(int nv, int dim) : num_original_vertices(nv), gdim(dim) {}
    };

    /// Get or create intersection vertex on edge
    /// Returns global vertex id
    /// If endpoint has phi==0, returns that endpoint's id (no new vertex)
    /// Otherwise creates (or reuses) edge intersection vertex
    int get_or_create_intersection_vid(
        CutContext& ctx,
        int gv0, int gv1,
        const double* p0, const double* p1,
        double phi0, double phi1);

} // namespace cutcells
