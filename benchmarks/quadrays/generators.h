// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// algoim's quadrature engines on clipped boxes, for comparisons with quadrays.
// Built only with CUTCELLS_WITH_ALGOIM; the library never includes algoim.

#pragma once

#include <functional>
#include <string>

#include <cutcells/quadrature.h>
#include <cutcells/quadrays/clipped_box.h>
#include <cutcells/selection_expr.h>

namespace cutcells::quadrays::benchmarks
{

/// Level set given by its values in physical space and its polynomial degree;
/// the generators interpolate it on each box.
struct LevelSet
{
    std::function<double(const Vec3<double>&)> value;
    int degree = 2;
};

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
    double split_alpha = 2.0; ///< split a box while its best alpha >= split_alpha (> 1: never);
                              ///< alpha is computed for this even when axes == algoim
    int split_depth = 0;      ///< maximum number of splits along a branch
};

struct GeneratorStats
{
    int boxes = 0;  ///< boxes integrated
    int splits = 0; ///< boxes split
};

/// algoim's multi-polynomial engine on a clipped box: the level set plus one
/// linear polynomial per clip plane; points outside the clips are discarded.
/// Points and weights are appended to @p rule as one rule of parent 0.
void algoim_clipped_box(const ClippedBox<double>& cell, const LevelSet& ls, const SelectionTerm& term, int q,
                        const GeneratorOptions& opt, quadrature::QuadratureRules<double>& rule,
                        GeneratorStats& stats);

/// algoim's 2015 engine (interval arithmetic, bisection) on an unclipped box,
/// for the level set |x - centre|^2 - radius^2.
void algoim_quadgen_sphere(const ClippedBox<double>& cell, const Vec3<double>& centre, double radius,
                           const SelectionTerm& term, int q, quadrature::QuadratureRules<double>& rule);

/// Named generator presets: algoim-auto, algoim-gl, gl-cellmask, alpha, split, alpha-split.
GeneratorOptions generator_preset(const std::string& name);

} // namespace cutcells::quadrays::benchmarks
