// Copyright (c) 2025 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "cut_context.h"
#include <cmath>
#include <limits>

namespace cutcells
{
    int get_or_create_intersection_vid(
        CutContext& ctx,
        int gv0, int gv1,
        const double* p0, const double* p1,
        double phi0, double phi1)
    {
        constexpr double eps = 1e-12;

        // If endpoint has phi==0, return that endpoint (no new vertex)
        if (std::abs(phi0) < eps)
            return gv0;
        if (std::abs(phi1) < eps)
            return gv1;

        // Create edge key
        EdgeKey key(gv0, gv1);

        // Check if intersection already exists
        auto it = ctx.edge_to_vid.find(key);
        if (it != ctx.edge_to_vid.end())
        {
            return it->second;
        }

        // Create new intersection vertex
        int new_vid = ctx.num_original_vertices + static_cast<int>(ctx.edge_to_vid.size());

        // Compute intersection point using linear interpolation
        double t = phi0 / (phi0 - phi1);
        
        // Add coordinates
        for (int d = 0; d < ctx.gdim; ++d)
        {
            double coord = p0[d] + t * (p1[d] - p0[d]);
            ctx.ip_coords.push_back(coord);
        }

        // Store mapping
        ctx.edge_to_vid[key] = new_vid;

        return new_vid;
    }

} // namespace cutcells
