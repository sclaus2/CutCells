// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT

#include "prism_midpoint_split.h"

#include <algorithm>
#include <string>

namespace cutcells::cell::prism_midpoint
{

namespace
{

constexpr std::array<LocalEdge, 9> prism_edge_list = {{
    {0, 1}, {1, 2}, {2, 0},
    {3, 4}, {4, 5}, {5, 3},
    {0, 3}, {1, 4}, {2, 5},
}};

[[noreturn]] void throw_invalid_prism_tokens(std::span<const int> prism_tokens,
                                             const std::string& message)
{
    std::string full = "prism_midpoint_split: " + message + " tokens=[";
    for (std::size_t i = 0; i < prism_tokens.size(); ++i)
    {
        if (i > 0)
            full += ",";
        full += std::to_string(prism_tokens[i]);
    }
    full += "]";
    throw std::runtime_error(full);
}

} // namespace

bool is_root_token(int token)
{
    return token >= 0 && token < 100;
}

std::span<const LocalEdge> prism_edges()
{
    return std::span<const LocalEdge>(prism_edge_list.data(), prism_edge_list.size());
}

family classify_family(std::span<const int> prism_tokens)
{
    if (prism_tokens.size() != 6)
        throw std::invalid_argument("classify_family: expected 6 prism tokens");

    int num_roots = 0;
    for (const int token : prism_tokens)
        num_roots += is_root_token(token) ? 1 : 0;

    if (num_roots == 3)
        return family::roots3;
    if (num_roots == 4)
        return family::roots4;

    throw_invalid_prism_tokens(
        prism_tokens,
        "tetra-derived prism must contain exactly 3 or 4 root vertices");
}

PrismRoleAnalysis analyze_tetra_derived_prism(std::span<const int> prism_tokens)
{
    if (prism_tokens.size() != 6)
        throw std::invalid_argument("analyze_tetra_derived_prism: expected 6 prism tokens");

    PrismRoleAnalysis out;
    out.split_family = classify_family(prism_tokens);

    for (int local_id = 0; local_id < 6; ++local_id)
    {
        const int token = prism_tokens[static_cast<std::size_t>(local_id)];
        if (is_root_token(token))
        {
            const std::size_t idx = static_cast<std::size_t>(out.num_roots);
            out.root_local_ids[idx] = local_id;
            ++out.num_roots;
        }
        else
        {
            const std::size_t idx = static_cast<std::size_t>(out.num_non_roots);
            out.non_root_local_ids[idx] = local_id;
            ++out.num_non_roots;
        }
    }

    const auto edges = prism_edges();
    for (const LocalEdge edge : edges)
    {
        const bool a_root = is_root_token(prism_tokens[static_cast<std::size_t>(edge.a)]);
        const bool b_root = is_root_token(prism_tokens[static_cast<std::size_t>(edge.b)]);
        if (!a_root && !b_root)
            out.midpoint_edges.push_back(edge);
    }

    const std::size_t expected_midpoints
        = (out.split_family == family::roots3) ? std::size_t(3) : std::size_t(1);
    if (out.midpoint_edges.size() != expected_midpoints)
    {
        throw_invalid_prism_tokens(
            prism_tokens,
            "unexpected number of uncut non-root prism edges for midpoint insertion");
    }

    return out;
}

} // namespace cutcells::cell::prism_midpoint
