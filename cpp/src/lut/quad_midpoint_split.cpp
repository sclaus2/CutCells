// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "quad_midpoint_split.h"

#include <string>

namespace cutcells::cell::quad_midpoint
{

namespace
{

[[noreturn]] void throw_invalid_quad_tokens(std::span<const int> quad_tokens,
                                            const std::string& message)
{
    std::string full = "quad_midpoint_split: " + message + " tokens=[";
    for (std::size_t i = 0; i < quad_tokens.size(); ++i)
    {
        if (i > 0)
            full += ",";
        full += std::to_string(quad_tokens[i]);
    }
    full += "]";
    throw std::runtime_error(full);
}

} // namespace

bool is_root_token(int token)
{
    return token >= 0 && token < 100;
}

QuadRoleAnalysis analyze_triangle_derived_quadrilateral(
    std::span<const int> quad_tokens)
{
    if (quad_tokens.size() != 4)
    {
        throw std::invalid_argument(
            "analyze_triangle_derived_quadrilateral: expected 4 quadrilateral tokens");
    }

    QuadRoleAnalysis out;

    for (int local_id = 0; local_id < 4; ++local_id)
    {
        const int token = quad_tokens[static_cast<std::size_t>(local_id)];
        if (is_root_token(token))
        {
            out.root_local_ids[static_cast<std::size_t>(out.num_roots)] = local_id;
            ++out.num_roots;
        }
        else
        {
            out.non_root_local_ids[static_cast<std::size_t>(out.num_non_roots)] = local_id;
            ++out.num_non_roots;
        }
    }

    if (out.num_roots != 2 || out.num_non_roots != 2)
    {
        throw_invalid_quad_tokens(
            quad_tokens,
            "triangle-derived quadrilateral must contain exactly 2 root and 2 non-root vertices");
    }

    bool found_canonical_order = false;
    for (int start = 0; start < 4; ++start)
    {
        const int v0 = start;
        const int v1 = (start + 1) % 4;
        const int v2 = (start + 2) % 4;
        const int v3 = (start + 3) % 4;
        if (is_root_token(quad_tokens[static_cast<std::size_t>(v0)])
            && is_root_token(quad_tokens[static_cast<std::size_t>(v1)])
            && !is_root_token(quad_tokens[static_cast<std::size_t>(v2)])
            && !is_root_token(quad_tokens[static_cast<std::size_t>(v3)]))
        {
            out.canonical_to_actual = {v0, v1, v2, v3};
            out.midpoint_edge = {v2, v3};
            found_canonical_order = true;
            break;
        }
    }

    if (!found_canonical_order)
    {
        throw_invalid_quad_tokens(
            quad_tokens,
            "expected cyclic ordering with two consecutive root vertices and two consecutive non-root vertices");
    }

    return out;
}

} // namespace cutcells::cell::quad_midpoint
