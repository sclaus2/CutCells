// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "debug.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace cutcells::quadrays
{

namespace
{

// The most points a report holds, and the grid of its level sets by the cell's dimension.
constexpr int report_points = 300000;
constexpr int report_grid[4] = {0, 201, 201, 41};
// Every n-th point, n prime and large: the points come in tensor products of a
// few points per segment, of which a small n would keep only some heights.
constexpr int report_strides[] = {1, 97, 197, 397, 797, 1597, 3203, 6421};

void write_string(std::FILE* f, const std::string& s)
{
    std::fputc('"', f);
    for (const char c : s)
    {
        if (c == '"' || c == '\\')
            std::fputc('\\', f);
        std::fputc(c, f);
    }
    std::fputc('"', f);
}

void write_key(std::FILE* f, const char* key)
{
    std::fprintf(f, ",\n\"%s\": ", key);
}

/// JSON has no NaN: the plot reads null as one.
template <std::floating_point T>
void write_number(std::FILE* f, T v, const char* format = "%.10g")
{
    if (std::isfinite(v))
        std::fprintf(f, format, static_cast<double>(v));
    else
        std::fputs("null", f);
}

template <std::floating_point T>
void write_numbers(std::FILE* f, std::span<const T> values, const char* format = "%.10g", std::size_t stride = 1)
{
    std::fputc('[', f);
    for (std::size_t i = 0; i < values.size(); i += stride)
    {
        if (i > 0)
            std::fputc(',', f);
        write_number(f, values[i], format);
    }
    std::fputc(']', f);
}

template <std::integral I>
void write_ints(std::FILE* f, std::span<const I> values)
{
    std::fputc('[', f);
    for (std::size_t i = 0; i < values.size(); ++i)
        std::fprintf(f, i > 0 ? ",%d" : "%d", static_cast<int>(values[i]));
    std::fputc(']', f);
}

} // namespace

template <std::floating_point T>
void write_cell_report(const std::string& path, const CellLabel<T>& label, const ClippedBox<T>& box,
                       std::span<const Source<T>> phis, const Options& opt, const BoxTrace<T>& trace,
                       const CellPoints<T>& points, const Stats& stats)
{
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f == nullptr)
        throw std::runtime_error("quadrays: cannot write the cell report " + path);

    std::fprintf(f, "{\"cell\": %lld", static_cast<long long>(label.cell));
    write_key(f, "cell_type");
    write_string(f, label.cell_type);
    std::fprintf(f, ",\n\"gdim\": %d, \"tdim\": %d", label.gdim, box.tdim);
    write_key(f, "vertices");
    write_numbers(f, std::span<const T>(label.vertices));
    write_key(f, "origin");
    write_numbers(f, std::span<const T>(box.origin));
    write_key(f, "jacobian");
    std::fputc('[', f);
    for (int i = 0; i < 3; ++i)
    {
        std::fputs(i > 0 ? "," : "", f);
        write_numbers(f, std::span<const T>(box.jacobian[i]));
    }
    std::fputc(']', f);
    write_key(f, "clips");
    std::fputc('[', f);
    for (std::size_t i = 0; i < box.clips.size(); ++i)
    {
        std::fputs(i > 0 ? "," : "", f);
        const T row[4] = {box.clips[i].c[0], box.clips[i].c[1], box.clips[i].c[2], box.clips[i].d};
        write_numbers(f, std::span<const T>(row, 4));
    }
    std::fputc(']', f);
    write_key(f, "level_sets");
    std::fputc('[', f);
    for (std::size_t i = 0; i < label.level_sets.size(); ++i)
    {
        std::fputs(i > 0 ? "," : "", f);
        write_string(f, label.level_sets[i]);
    }
    std::fputc(']', f);
    write_key(f, "selection");
    write_string(f, label.selection);
    std::fprintf(f, ",\n\"dim\": %d, \"order\": %d, \"q\": %d, \"seconds\": %.6f", label.dim, label.order, label.q,
                 label.seconds);
    std::fprintf(f, ",\n\"options\": {\"margin\": %g, \"tolerance\": %g, \"max_depth\": %d, \"max_bisections\": %d, "
                    "\"max_points\": %d, \"rotation_depth\": %d, \"fold_splits\": %s, \"anchor_folds\": %d}",
                 opt.margin, opt.tolerance, opt.max_depth, opt.max_bisections, point_limit(opt, label.q, box.tdim),
                 opt.rotation_depth, opt.fold_splits ? "true" : "false", opt.anchor_folds);
    std::fprintf(f, ",\n\"stats\": {\"bisections\": %d, \"uncertified\": %d, \"rotations\": %d, \"fold_splits\": %d, "
                    "\"surfaces\": %d, \"points_limited\": %d, \"bisections_limited\": %d, \"curve_lost\": %d}",
                 stats.bisections, stats.uncertified, stats.rotations, stats.fold_splits, stats.surfaces,
                 stats.points_limited, stats.bisections_limited, stats.curve_lost);
    std::fprintf(f, ",\n\"n_points\": %d", points.n_points());

    // the level sets on the grid, u0 slowest
    const int n = report_grid[box.tdim];
    std::fprintf(f, ",\n\"grid\": %d", n);
    write_key(f, "values");
    std::fputc('[', f);
    std::vector<T> values;
    for (std::size_t l = 0; l < phis.size(); ++l)
    {
        values.clear();
        const int n1 = box.tdim >= 2 ? n : 1, n2 = box.tdim >= 3 ? n : 1;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n1; ++j)
                for (int k = 0; k < n2; ++k)
                {
                    const T u[3] = {static_cast<T>(i) / static_cast<T>(n - 1),
                                    box.tdim >= 2 ? static_cast<T>(j) / static_cast<T>(n - 1) : T(0),
                                    box.tdim >= 3 ? static_cast<T>(k) / static_cast<T>(n - 1) : T(0)};
                    values.push_back(evaluate(phis[l], std::span<const T>(u, 3)));
                }
        std::fputs(l > 0 ? ",\n" : "", f);
        write_numbers(f, std::span<const T>(values), "%.6g");
    }
    std::fputc(']', f);

    // the boxes
    std::fputs(",\n\"boxes\": {\"parent\": ", f);
    write_ints(f, std::span<const int>(trace.parent));
    write_key(f, "level");
    write_ints(f, std::span<const int>(trace.level));
    write_key(f, "depth");
    write_ints(f, std::span<const int>(trace.depth));
    write_key(f, "height");
    write_ints(f, std::span<const int>(trace.height));
    write_key(f, "fate");
    std::fputc('[', f);
    for (std::size_t i = 0; i < trace.fate.size(); ++i)
        std::fprintf(f, i > 0 ? ",%d" : "%d", static_cast<int>(trace.fate[i]));
    std::fputc(']', f);
    write_key(f, "rotated");
    write_ints(f, std::span<const std::uint8_t>(trace.rotated));
    write_key(f, "lo");
    write_numbers(f, std::span<const T>(trace.lo));
    write_key(f, "hi");
    write_numbers(f, std::span<const T>(trace.hi));
    write_key(f, "frame");
    write_numbers(f, std::span<const T>(trace.frame));
    write_key(f, "points");
    write_ints(f, std::span<const int>(trace.points));
    std::fputc('}', f);

    // the points, in box coordinates
    std::size_t stride = 1;
    for (const int n : report_strides)
    {
        stride = static_cast<std::size_t>(n);
        if (points.n_points() / n <= report_points)
            break;
    }
    std::fprintf(f, ",\n\"point_stride\": %d", static_cast<int>(stride));
    write_key(f, "points");
    std::fputc('[', f);
    for (std::size_t p = 0; p < points.weights.size(); p += stride)
    {
        std::fputs(p > 0 ? "," : "", f);
        write_number(f, points.points[3 * p], "%.6g");
        std::fputc(',', f);
        write_number(f, points.points[3 * p + 1], "%.6g");
        std::fputc(',', f);
        write_number(f, points.points[3 * p + 2], "%.6g");
    }
    std::fputs("]}\n", f);
    std::fclose(f);
}

template void write_cell_report<float>(const std::string&, const CellLabel<float>&, const ClippedBox<float>&,
                                       std::span<const Source<float>>, const Options&, const BoxTrace<float>&,
                                       const CellPoints<float>&, const Stats&);
template void write_cell_report<double>(const std::string&, const CellLabel<double>&, const ClippedBox<double>&,
                                        std::span<const Source<double>>, const Options&, const BoxTrace<double>&,
                                        const CellPoints<double>&, const Stats&);

} // namespace cutcells::quadrays
