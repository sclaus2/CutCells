// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "clipped_box.h"
#include "selection_expr.h"

namespace cutcells::proto
{

/// Quadrature rule for one cell and one selection: points in the parent cell's
/// reference coordinates (flat, 3 per point) and physical weights, as in
/// quadrature::QuadratureRules.
struct Rule
{
    std::vector<double> points;
    std::vector<double> weights;

    int n_points() const { return static_cast<int>(weights.size()); }
};

/// Level set given by its values in physical space and its polynomial degree.
/// Generators interpolate it on each box (the engine will start from the
/// cell's Bernstein coefficients instead).
struct LevelSet
{
    std::function<double(const Vec3&)> value;
    int degree = 2;
};

/// What a single-level-set selection term selects.
enum class PartKind
{
    negative,  ///< phi < 0
    positive,  ///< phi > 0
    interface, ///< phi = 0
    whole      ///< no condition
};

/// Classify a compiled single-level-set term; throws for terms on other level sets.
PartKind part_kind(const SelectionTerm& term);

enum class AxisChoice
{
    algoim, ///< algoim's score over the whole box
    alpha   ///< Cui et al.'s angle indicator on the clipped region
};

struct GeneratorOptions
{
    std::string name;
    bool gauss_legendre = true; ///< false: algoim's AutoMixed (tanh-sinh where it sees vertical tangents)
    bool cell_masks = false;    ///< restrict algoim's masks to sub-cells that meet the clipped region
    AxisChoice axes = AxisChoice::algoim;
    double split_alpha = 2.0;   ///< split a box while its best alpha >= split_alpha (> 1: never);
                                ///< alpha is computed for this even when axes == algoim
    int split_depth = 0;        ///< maximum number of splits along a branch
};

struct GeneratorStats
{
    int boxes = 0;  ///< boxes integrated
    int splits = 0; ///< boxes split
};

/// algoim's multi-polynomial engine on a clipped box: the level set plus one
/// linear polynomial per clip plane; points outside the clips are discarded.
void algoim_clipped_box(const ClippedBox& cell, const LevelSet& ls, const SelectionTerm& term,
                        int q, const GeneratorOptions& opt, Rule& rule, GeneratorStats& stats);

/// algoim's 2015 engine (interval arithmetic, bisection) on an unclipped box,
/// for the level set |x - centre|^2 - radius^2.
void algoim_quadgen_sphere(const ClippedBox& cell, const Vec3& centre, double radius,
                           const SelectionTerm& term, int q, Rule& rule);

/// Named generator presets used by the study.
GeneratorOptions generator_preset(const std::string& name);

} // namespace cutcells::proto
