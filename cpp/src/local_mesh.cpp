// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#include "local_mesh.h"
#include "green_split.h"
#include "cell_subdivision.h"
#include "cell_topology.h"
#include "edge_classification.h"
#include "mapping.h"
#include "cell_flags.h"
#include "cut_cell.h"
#include "triangulation_repair.h"
#include "interface_split.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace cutcells
{

namespace
{
inline bool debug_multiroot_enabled()
{
    const char* env = std::getenv("CUTCELLS_DEBUG_MULTIROOT");
    if (env == nullptr)
        return false;
    const std::string_view value(env);
    return value == "1" || value == "true" || value == "TRUE";
}

struct EdgeKey
{
    int32_t a = -1;
    int32_t b = -1;

    bool operator==(const EdgeKey& other) const noexcept
    {
        return a == other.a && b == other.b;
    }
};

struct EdgeKeyHash
{
    std::size_t operator()(const EdgeKey& k) const noexcept
    {
        const auto h0 = std::hash<int32_t>{}(k.a);
        const auto h1 = std::hash<int32_t>{}(k.b);
        return h0 ^ (h1 + 0x9e3779b97f4a7c15ULL + (h0 << 6) + (h0 >> 2));
    }
};

struct ZeroOwnerKey
{
    uint64_t zero_mask = 0;
    uint8_t dim = 0;
    int32_t parent_dim = -1;
    int32_t parent_id = -1;

    bool operator==(const ZeroOwnerKey& other) const noexcept
    {
        return zero_mask == other.zero_mask
            && dim == other.dim
            && parent_dim == other.parent_dim
            && parent_id == other.parent_id;
    }
};

struct ZeroOwnerKeyHash
{
    std::size_t operator()(const ZeroOwnerKey& k) const noexcept
    {
        std::size_t h = std::hash<uint64_t>{}(k.zero_mask);
        h ^= std::hash<uint8_t>{}(k.dim) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>{}(k.parent_dim) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int32_t>{}(k.parent_id) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

template <std::floating_point T>
int find_local_edge_by_vertices(
    const LocalMesh<T>& mesh,
    int32_t             va,
    int32_t             vb)
{
    for (int e = 0; e < mesh.n_edges(); ++e)
    {
        const int32_t ev0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int32_t ev1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if ((ev0 == va && ev1 == vb) || (ev0 == vb && ev1 == va))
            return e;
    }
    return -1;
}

template <std::floating_point T, typename EvalPhi, typename EvalGrad>
cell::edge_root::RootSolveInfo<T> solve_edge_root_on_segment(
    std::span<const T>      p0,
    std::span<const T>      p1,
    EvalPhi&&               eval_phi,
    EvalGrad&&              eval_grad,
    bool                    has_exact_grad,
    cell::edge_root::method root_method,
    T                       tol)
{
    using RootInfo = cell::edge_root::RootSolveInfo<T>;
    if (p0.size() != p1.size())
        throw std::invalid_argument("solve_edge_root_on_segment: inconsistent point dimensions");

    std::vector<T> x(static_cast<std::size_t>(p0.size()), T(0));
    std::vector<T> direction(static_cast<std::size_t>(p0.size()), T(0));
    for (std::size_t d = 0; d < p0.size(); ++d)
        direction[d] = p1[d] - p0[d];

    int eval_count = 0;
    auto eval = [&](T t) -> T
    {
        for (std::size_t d = 0; d < p0.size(); ++d)
            x[d] = p0[d] + t * direction[d];
        ++eval_count;
        return eval_phi(std::span<const T>(x.data(), x.size()));
    };

    RootInfo best;
    best.t = T(0.5);
    best.residual = std::numeric_limits<T>::max();

    const auto update_best = [&](T t, T value)
    {
        const T residual = std::abs(value);
        if (residual < best.residual)
        {
            best.t = std::clamp(t, T(0), T(1));
            best.residual = residual;
        }
    };

    const T f0 = eval(T(0));
    const T f1 = eval(T(1));
    best.evaluations = 2;
    update_best(T(0), f0);
    update_best(T(1), f1);

    if (std::abs(f0) <= tol)
    {
        best.t = T(0);
        best.residual = std::abs(f0);
        best.converged = true;
        return best;
    }
    if (std::abs(f1) <= tol)
    {
        best.t = T(1);
        best.residual = std::abs(f1);
        best.converged = true;
        return best;
    }

    const auto solve_on_bracket = [&](T a, T b, T fa, T fb) -> RootInfo
    {
        RootInfo info;
        const T t_linear = cell::edge_root::linear_root_parameter<T>(fa, fb, T(0));
        int iterations = 0;
        bool converged = false;
        T t = t_linear;

        if (root_method == cell::edge_root::method::linear)
        {
            t = t_linear;
        }
        else if (root_method == cell::edge_root::method::brent)
        {
            t = cell::edge_root::brent_parameter<T>(
                eval, a, b, fa, fb, 64, tol, tol, &iterations, &converged);
        }
        else if (root_method == cell::edge_root::method::newton)
        {
            t = cell::edge_root::newton_parameter<T>(
                eval, a, b, fa, fb, t_linear, 64, tol, tol, &iterations, &converged);
        }
        else
        {
            t = cell::edge_root::itp_parameter<T>(
                eval, a, b, fa, fb, t_linear, 64, tol, tol, &iterations, &converged);
        }

        const T f = eval(t);
        info.t = std::clamp(t, T(0), T(1));
        info.residual = std::abs(f);
        info.iterations = iterations;
        info.evaluations = 1;
        info.converged = converged || (info.residual <= tol);
        return info;
    };

    bool have_bracket = (f0 * f1 < T(0));
    T a = T(0);
    T b = T(1);
    T fa = f0;
    T fb = f1;

    if (!have_bracket)
    {
        const int n_scan = 32;
        T prev_t = T(0);
        T prev_f = f0;
        for (int i = 1; i <= n_scan; ++i)
        {
            const T t = static_cast<T>(i) / static_cast<T>(n_scan);
            const T f = eval(t);
            update_best(t, f);
            if (std::abs(f) <= tol)
            {
                best.t = t;
                best.residual = std::abs(f);
                best.evaluations = eval_count;
                best.converged = true;
                return best;
            }
            if (prev_f * f < T(0))
            {
                a = prev_t;
                b = t;
                fa = prev_f;
                fb = f;
                have_bracket = true;
                break;
            }
            prev_t = t;
            prev_f = f;
        }
    }

    if (have_bracket)
    {
        RootInfo info = solve_on_bracket(a, b, fa, fb);
        info.evaluations += eval_count;
        return info;
    }

    if (root_method == cell::edge_root::method::newton && has_exact_grad)
    {
        std::vector<T> grad(static_cast<std::size_t>(p0.size()), T(0));
        T t = best.t;
        for (int iter = 0; iter < 64; ++iter)
        {
            const T f = eval(t);
            update_best(t, f);
            if (std::abs(f) <= tol)
            {
                best.t = t;
                best.residual = std::abs(f);
                best.evaluations = eval_count;
                best.iterations = iter + 1;
                best.converged = true;
                return best;
            }

            for (std::size_t d = 0; d < p0.size(); ++d)
                x[d] = p0[d] + t * direction[d];
            eval_grad(
                std::span<const T>(x.data(), x.size()),
                std::span<T>(grad.data(), grad.size()));

            T dphidt = T(0);
            for (std::size_t d = 0; d < p0.size(); ++d)
                dphidt += grad[d] * direction[d];

            if (!std::isfinite(dphidt) || std::abs(dphidt) <= T(1e-14))
                break;

            const T t_newton = std::clamp(t - f / dphidt, T(0), T(1));
            if (std::abs(t_newton - t) <= tol)
                break;
            t = t_newton;
        }
    }

    RootInfo info = cell::edge_root::find_root_parameter_info<T>(
        p0,
        p1,
        [&](std::span<const T> xp) -> T { return eval_phi(xp); },
        root_method,
        T(0),
        64,
        tol,
        tol);
    info.evaluations += eval_count;
    if (info.residual < best.residual)
        return info;
    best.evaluations = eval_count;
    return best;
}

inline int background_edge_id_from_corners(cell::type ct, int pv0, int pv1)
{
    const auto parent_edges = cell::edges(ct);
    for (int pe = 0; pe < static_cast<int>(parent_edges.size()); ++pe)
    {
        const int e0 = parent_edges[static_cast<std::size_t>(pe)][0];
        const int e1 = parent_edges[static_cast<std::size_t>(pe)][1];
        if ((pv0 == e0 && pv1 == e1) || (pv0 == e1 && pv1 == e0))
            return pe;
    }
    return -1;
}

inline bool background_edge_on_face(cell::type ct, int edge_id, int face_id)
{
    const auto edges = cell::edges(ct);
    const auto faces = cell::faces(ct);
    const int n_faces = cell::num_faces(ct);
    if (edge_id < 0 || edge_id >= static_cast<int>(edges.size()))
        return false;
    if (face_id < 0 || face_id >= n_faces)
        return false;

    const auto edge = edges[static_cast<std::size_t>(edge_id)];
    const auto face = faces[static_cast<std::size_t>(face_id)];
    const auto on_face = [&](int v)
    {
        for (const int fv : face)
        {
            if (fv == v)
                return true;
        }
        return false;
    };
    return on_face(edge[0]) && on_face(edge[1]);
}

template <std::floating_point T>
bool vertex_lies_on_background_face(const LocalMesh<T>& mesh, int lv, int bg_face_id)
{
    const int n_bg_faces = cell::num_faces(mesh.parent_cell_type);
    if (bg_face_id < 0 || bg_face_id >= n_bg_faces)
        return false;

    const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
    const int pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
    if (pdim == 2 && pid == bg_face_id)
        return true;
    if (pdim == 1 && pid >= 0 && pid < cell::num_edges(mesh.parent_cell_type))
        return background_edge_on_face(mesh.parent_cell_type, pid, bg_face_id);
    if (pdim == 0)
    {
        const auto bg_faces = cell::faces(mesh.parent_cell_type);
        const auto face = bg_faces[static_cast<std::size_t>(bg_face_id)];
        for (const int fv : face)
        {
            if (pid == fv)
                return true;
        }
    }
    return false;
}

template <std::floating_point T>
std::pair<int32_t, int32_t> infer_background_face_parent(
    const LocalMesh<T>&             mesh,
    std::span<const int32_t>        face_vertices)
{
    const int n_bg_faces = cell::num_faces(mesh.parent_cell_type);
    for (int f = 0; f < n_bg_faces; ++f)
    {
        bool all_on_face = true;
        for (const int32_t lv : face_vertices)
        {
            if (!vertex_lies_on_background_face(mesh, lv, f))
            {
                all_on_face = false;
                break;
            }
        }
        if (all_on_face)
            return {2, f};
    }
    return {-1, -1};
}

template <std::floating_point T>
std::pair<int32_t, int32_t> infer_background_edge_or_face_parent(
    const LocalMesh<T>& mesh,
    int                 lv0,
    int                 lv1)
{
    const auto edge_parent = infer_background_edge_parent(mesh, lv0, lv1);
    if (edge_parent.first >= 0)
        return edge_parent;

    const std::array<int32_t, 2> verts = {
        static_cast<int32_t>(lv0),
        static_cast<int32_t>(lv1)};
    return infer_background_face_parent(
        mesh, std::span<const int32_t>(verts.data(), verts.size()));
}

template <std::floating_point T>
inline uint64_t bits_key(T v)
{
    if constexpr (sizeof(T) == sizeof(uint64_t))
        return std::bit_cast<uint64_t>(v);
    else
        return static_cast<uint64_t>(std::bit_cast<uint32_t>(v));
}

template <std::floating_point T>
struct VertexDedupKey
{
    int32_t parent_dim = -1;
    int32_t parent_id = -1;
    std::array<uint64_t, 3> x = {0, 0, 0};
    std::array<uint64_t, 3> xref = {0, 0, 0};
    uint8_t gdim = 0;
    uint8_t tdim = 0;
    uint8_t has_ref = 0;

    bool operator==(const VertexDedupKey& other) const noexcept
    {
        return parent_dim == other.parent_dim
               && parent_id == other.parent_id
               && x == other.x
               && xref == other.xref
               && gdim == other.gdim
               && tdim == other.tdim
               && has_ref == other.has_ref;
    }
};

template <std::floating_point T>
struct VertexDedupKeyHash
{
    std::size_t operator()(const VertexDedupKey<T>& k) const noexcept
    {
        std::size_t h = std::hash<int32_t>{}(k.parent_dim);
        auto mix = [&](std::size_t v) {
            h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        mix(std::hash<int32_t>{}(k.parent_id));
        mix(std::hash<uint64_t>{}(k.x[0]));
        mix(std::hash<uint64_t>{}(k.x[1]));
        mix(std::hash<uint64_t>{}(k.x[2]));
        mix(std::hash<uint64_t>{}(k.xref[0]));
        mix(std::hash<uint64_t>{}(k.xref[1]));
        mix(std::hash<uint64_t>{}(k.xref[2]));
        mix(std::hash<uint8_t>{}(k.gdim));
        mix(std::hash<uint8_t>{}(k.tdim));
        mix(std::hash<uint8_t>{}(k.has_ref));
        return h;
    }
};

} // namespace  // end of first anonymous namespace

template <std::floating_point T>
void build_local_edges(LocalMesh<T>& mesh)
{
    const int nc = mesh.n_cells();
    std::vector<std::pair<int32_t, int32_t>> unique_edges;
    std::unordered_map<EdgeKey, int32_t, EdgeKeyHash> edge_to_id;
    edge_to_id.reserve(static_cast<std::size_t>(nc * 8));
    mesh.cell_edge_offsets.resize(nc + 1);
    mesh.cell_edges_flat.clear();

    for (int i = 0; i < nc; ++i)
    {
        mesh.cell_edge_offsets[i] = static_cast<int32_t>(mesh.cell_edges_flat.size());
        const int v_start = mesh.cell_offsets[i];
        const auto cell_typ = mesh.cell_types[i];
        const auto edge_patterns = cell::edges(cell_typ);
        for (const auto& ep : edge_patterns)
        {
            int32_t v0 = mesh.cell_vertices[v_start + ep[0]];
            int32_t v1 = mesh.cell_vertices[v_start + ep[1]];
            if (v0 > v1)
                std::swap(v0, v1);

            const EdgeKey key{v0, v1};
            const auto it = edge_to_id.find(key);
            int32_t edge_idx = 0;
            if (it == edge_to_id.end())
            {
                edge_idx = static_cast<int32_t>(unique_edges.size());
                unique_edges.push_back({v0, v1});
                edge_to_id.emplace(key, edge_idx);
            }
            else
            {
                edge_idx = it->second;
            }
            mesh.cell_edges_flat.push_back(edge_idx);
        }
    }
    mesh.cell_edge_offsets[nc] = static_cast<int32_t>(mesh.cell_edges_flat.size());

    const int n_edges = static_cast<int>(unique_edges.size());
    mesh.edge_vertices.resize(n_edges * 2);
    for (int i = 0; i < n_edges; ++i)
    {
        mesh.edge_vertices[i * 2] = unique_edges[i].first;
        mesh.edge_vertices[i * 2 + 1] = unique_edges[i].second;
    }
    mesh.edge_parent_dim.assign(n_edges, -1);
    mesh.edge_parent_id.assign(n_edges, -1);
    for (int i = 0; i < n_edges; ++i)
    {
        const auto [pdim, pid] = infer_background_edge_or_face_parent(
            mesh,
            mesh.edge_vertices[static_cast<std::size_t>(2 * i)],
            mesh.edge_vertices[static_cast<std::size_t>(2 * i + 1)]);
        mesh.edge_parent_dim[static_cast<std::size_t>(i)] = pdim;
        mesh.edge_parent_id[static_cast<std::size_t>(i)] = pid;
    }
    const int n_ls = mesh.n_level_sets;
    const int ne_ls = n_edges * n_ls;
    mesh.edge_state.assign(static_cast<std::size_t>(ne_ls), static_cast<uint8_t>(EdgeState::uncertain));
    mesh.edge_cert.assign(
        static_cast<std::size_t>(ne_ls),
        static_cast<uint8_t>(EdgeCertification::unresolved));
    mesh.edge_root_vertex.assign(static_cast<std::size_t>(ne_ls), -1);
    mesh.edge_root_parameter.assign(static_cast<std::size_t>(ne_ls), T(-1));
    mesh.edge_root_iterations.assign(static_cast<std::size_t>(ne_ls), 0);
    mesh.edge_root_evaluations.assign(static_cast<std::size_t>(ne_ls), 0);
    mesh.edge_root_converged.assign(static_cast<std::size_t>(ne_ls), 0);
    mesh.edge_root_residual.assign(static_cast<std::size_t>(ne_ls), T(0));
}

inline std::span<const std::array<int, 2>> cut_edge_order(cell::type ct)
{
    // Edge numbering follows the clipping/cut tables (VTK-style), not Basix.
    static constexpr std::array<std::array<int, 2>, 1> interval = {{{0, 1}}};
    static constexpr std::array<std::array<int, 2>, 3> triangle = {{{0, 1}, {1, 2}, {2, 0}}};
    static constexpr std::array<std::array<int, 2>, 4> quadrilateral = {{{0, 1}, {1, 2}, {2, 3}, {3, 0}}};
    static constexpr std::array<std::array<int, 2>, 6> tetrahedron = {{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}}};
    static constexpr std::array<std::array<int, 2>, 12> hexahedron = {{
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {3, 7}, {2, 6}
    }};
    static constexpr std::array<std::array<int, 2>, 9> prism = {{
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3}, {1, 4}, {2, 5}
    }};
    static constexpr std::array<std::array<int, 2>, 8> pyramid = {{
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {0, 4}, {1, 4}, {2, 4}, {3, 4}
    }};

    switch (ct)
    {
    case cell::type::interval: return std::span(interval);
    case cell::type::triangle: return std::span(triangle);
    case cell::type::quadrilateral: return std::span(quadrilateral);
    case cell::type::tetrahedron: return std::span(tetrahedron);
    case cell::type::hexahedron: return std::span(hexahedron);
    case cell::type::prism: return std::span(prism);
    case cell::type::pyramid: return std::span(pyramid);
    default:
        throw std::invalid_argument("cut_edge_order: unsupported cell type");
    }
}

template <std::floating_point T>
int local_edge_for_cut_edge_token(const LocalMesh<T>& mesh,
                                  const std::span<const int32_t> parent_cell_vertices,
                                  const std::span<const int32_t> parent_cell_edges,
                                  const cell::type ct,
                                  const int token_edge)
{
    const auto cut_edges = cut_edge_order(ct);
    if (token_edge < 0 || token_edge >= static_cast<int>(cut_edges.size()))
        return -1;

    const int lv0 = cut_edges[static_cast<std::size_t>(token_edge)][0];
    const int lv1 = cut_edges[static_cast<std::size_t>(token_edge)][1];
    if (lv0 < 0 || lv1 < 0
        || lv0 >= static_cast<int>(parent_cell_vertices.size())
        || lv1 >= static_cast<int>(parent_cell_vertices.size()))
        return -1;

    const int32_t gv0 = parent_cell_vertices[static_cast<std::size_t>(lv0)];
    const int32_t gv1 = parent_cell_vertices[static_cast<std::size_t>(lv1)];

    for (const int32_t local_mesh_e : parent_cell_edges)
    {
        if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
            continue;
        const int32_t e0 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e)];
        const int32_t e1 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e + 1)];
        if ((e0 == gv0 && e1 == gv1) || (e0 == gv1 && e1 == gv0))
            return local_mesh_e;
    }

    return -1;
}

template <std::floating_point T>
void deduplicate_vertices_exact(LocalMesh<T>& mesh)
{
    const int old_nv = mesh.n_vertices();
    if (old_nv <= 1)
        return;

    const bool has_ref = !mesh.vertex_ref_x.empty();
    const int tdim = mesh.tdim;
    const int gdim = mesh.gdim;

    std::vector<int32_t> old_to_new(static_cast<std::size_t>(old_nv), -1);
    std::vector<T> new_x;
    std::vector<T> new_ref;
    std::vector<int32_t> new_pdim;
    std::vector<int32_t> new_pid;
    std::vector<int32_t> new_root_edge_id;
    new_x.reserve(mesh.vertex_x.size());
    if (has_ref)
        new_ref.reserve(mesh.vertex_ref_x.size());
    new_pdim.reserve(mesh.vertex_parent_dim.size());
    new_pid.reserve(mesh.vertex_parent_id.size());
    new_root_edge_id.reserve(mesh.vertex_root_edge_id.size());

    std::unordered_map<VertexDedupKey<T>, int32_t, VertexDedupKeyHash<T>> key_to_new;
    key_to_new.reserve(static_cast<std::size_t>(old_nv * 2));

    auto make_key = [&](int old_v) -> VertexDedupKey<T>
    {
        VertexDedupKey<T> key;
        key.parent_dim = mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)];
        key.parent_id = mesh.vertex_parent_id[static_cast<std::size_t>(old_v)];
        key.gdim = static_cast<uint8_t>(gdim);
        key.tdim = static_cast<uint8_t>(tdim);
        key.has_ref = static_cast<uint8_t>(has_ref ? 1 : 0);
        for (int d = 0; d < gdim; ++d)
            key.x[static_cast<std::size_t>(d)] = bits_key(
                mesh.vertex_x[static_cast<std::size_t>(old_v * gdim + d)]);
        if (has_ref)
        {
            for (int d = 0; d < tdim; ++d)
                key.xref[static_cast<std::size_t>(d)] = bits_key(
                    mesh.vertex_ref_x[static_cast<std::size_t>(old_v * tdim + d)]);
        }
        return key;
    };

    int new_nv = 0;
    for (int old_v = 0; old_v < old_nv; ++old_v)
    {
        const auto key = make_key(old_v);
        const auto it = key_to_new.find(key);
        int found = -1;
        if (it == key_to_new.end())
        {
            found = new_nv++;
            key_to_new.emplace(key, found);
            for (int d = 0; d < gdim; ++d)
                new_x.push_back(mesh.vertex_x[static_cast<std::size_t>(old_v * gdim + d)]);
            if (has_ref)
            {
                for (int d = 0; d < tdim; ++d)
                    new_ref.push_back(mesh.vertex_ref_x[static_cast<std::size_t>(old_v * tdim + d)]);
            }
            new_pdim.push_back(mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
            new_pid.push_back(mesh.vertex_parent_id[static_cast<std::size_t>(old_v)]);
            new_root_edge_id.push_back(mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)]);
        }
        else
            found = it->second;
        old_to_new[static_cast<std::size_t>(old_v)] = found;
    }

    if (new_nv == old_nv)
        return;

    for (auto& cv : mesh.cell_vertices)
        cv = old_to_new[static_cast<std::size_t>(cv)];

    mesh.vertex_x.swap(new_x);
    mesh.vertex_ref_x.swap(new_ref);
    mesh.vertex_parent_dim.swap(new_pdim);
    mesh.vertex_parent_id.swap(new_pid);
    mesh.vertex_root_edge_id.swap(new_root_edge_id);
}

template <std::floating_point T>
void rebuild_parent_entity_maps(LocalMesh<T>& mesh)
{
    const int n_parent_vertices = cell::get_num_vertices(mesh.parent_cell_type);
    const int n_parent_edges = cell::num_edges(mesh.parent_cell_type);

    mesh.parent_vertex_to_local_vertex.assign(static_cast<std::size_t>(n_parent_vertices), -1);
    mesh.parent_edge_to_local_edge.assign(static_cast<std::size_t>(n_parent_edges), -1);

    for (int lv = 0; lv < mesh.n_vertices(); ++lv)
    {
        if (mesh.vertex_parent_dim[static_cast<std::size_t>(lv)] != 0)
            continue;
        const int pv = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
        if (pv >= 0 && pv < n_parent_vertices && mesh.parent_vertex_to_local_vertex[static_cast<std::size_t>(pv)] < 0)
            mesh.parent_vertex_to_local_vertex[static_cast<std::size_t>(pv)] = lv;
    }

    const auto parent_edges = cell::edges(mesh.parent_cell_type);

    // First attempt: match local edges that connect two parent corner vertices.
    for (int le = 0; le < mesh.n_edges(); ++le)
    {
        const int lv0 = mesh.edge_vertices[static_cast<std::size_t>(2 * le)];
        const int lv1 = mesh.edge_vertices[static_cast<std::size_t>(2 * le + 1)];

        if (mesh.vertex_parent_dim[static_cast<std::size_t>(lv0)] != 0
            || mesh.vertex_parent_dim[static_cast<std::size_t>(lv1)] != 0)
            continue;

        const int pv0 = mesh.vertex_parent_id[static_cast<std::size_t>(lv0)];
        const int pv1 = mesh.vertex_parent_id[static_cast<std::size_t>(lv1)];
        for (int pe = 0; pe < static_cast<int>(parent_edges.size()); ++pe)
        {
            const int e0 = parent_edges[static_cast<std::size_t>(pe)][0];
            const int e1 = parent_edges[static_cast<std::size_t>(pe)][1];
            if ((pv0 == e0 && pv1 == e1) || (pv0 == e1 && pv1 == e0))
            {
                if (mesh.parent_edge_to_local_edge[static_cast<std::size_t>(pe)] < 0)
                    mesh.parent_edge_to_local_edge[static_cast<std::size_t>(pe)] = le;
                break;
            }
        }
    }

    // Rewrite only parent-cell interior ids so they remain stable across later
    // local refinement. Edge ancestry stays in background-edge numbering.
    for (int lv = 0; lv < mesh.n_vertices(); ++lv)
    {
        const int pdim = mesh.vertex_parent_dim[static_cast<std::size_t>(lv)];
        int& pid = mesh.vertex_parent_id[static_cast<std::size_t>(lv)];
        if (pdim == mesh.tdim)
        {
            pid = mesh.parent_cell_id;
        }
    }
}

int infer_lagrange_order(cell::type ct, int n_vertices)
{
    if (n_vertices == cell::get_num_vertices(ct))
        return 1;

    for (int p = 2; p <= 4; ++p)
    {
        if (iso_p1_template(ct, p).n_vertices == n_vertices)
            return p;
    }
    throw std::invalid_argument("init_local_mesh_from_cell: could not infer interpolation order from node count");
}

std::pair<int, int> infer_gdim_and_order(cell::type ct, std::size_t flat_coord_size)
{
    std::vector<std::pair<int, int>> matches;
    for (int gdim = 1; gdim <= 3; ++gdim)
    {
        if (flat_coord_size % static_cast<std::size_t>(gdim) != 0)
            continue;
        const int n_vertices = static_cast<int>(flat_coord_size / static_cast<std::size_t>(gdim));
        try
        {
            const int order = infer_lagrange_order(ct, n_vertices);
            matches.push_back({gdim, order});
        }
        catch (const std::invalid_argument&)
        {
            // Try next candidate.
        }
    }

    if (matches.empty())
        throw std::invalid_argument("init_local_mesh_from_cell: could not infer geometric dimension and interpolation order");
    if (matches.size() > 1)
        throw std::invalid_argument("init_local_mesh_from_cell: ambiguous geometric dimension/order from coordinate array");
    return matches[0];
}

template <std::floating_point T>
int32_t append_generated_vertex(
    std::vector<T>&        vertex_x,
    std::vector<T>&        vertex_ref_x,
    std::vector<int32_t>&  vertex_parent_dim,
    std::vector<int32_t>&  vertex_parent_id,
    std::vector<int32_t>&  vertex_root_edge_id,
    std::span<const T>     x_ref,
    std::span<const T>     x_phys,
    int32_t                parent_dim,
    int32_t                parent_id)
{
    const int32_t new_id = static_cast<int32_t>(vertex_parent_dim.size());
    vertex_x.insert(vertex_x.end(), x_phys.begin(), x_phys.end());
    vertex_ref_x.insert(vertex_ref_x.end(), x_ref.begin(), x_ref.end());
    vertex_parent_dim.push_back(parent_dim);
    vertex_parent_id.push_back(parent_id);
    vertex_root_edge_id.push_back(-1);
    return new_id;
}

template <std::floating_point T>
int32_t append_segment_midpoint(
    const LocalMesh<T>&    old,
    int32_t                va,
    int32_t                vb,
    std::vector<T>&        vertex_x,
    std::vector<T>&        vertex_ref_x,
    std::vector<int32_t>&  vertex_parent_dim,
    std::vector<int32_t>&  vertex_parent_id,
    std::vector<int32_t>&  vertex_root_edge_id)
{
    std::array<T, 3> x_ref = {T(0), T(0), T(0)};
    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    for (int d = 0; d < old.tdim; ++d)
    {
        x_ref[static_cast<std::size_t>(d)]
            = T(0.5) * (old.vertex_ref_x[static_cast<std::size_t>(va * old.tdim + d)]
                        + old.vertex_ref_x[static_cast<std::size_t>(vb * old.tdim + d)]);
    }
    for (int d = 0; d < old.gdim; ++d)
    {
        x_phys[static_cast<std::size_t>(d)]
            = T(0.5) * (old.vertex_x[static_cast<std::size_t>(va * old.gdim + d)]
                        + old.vertex_x[static_cast<std::size_t>(vb * old.gdim + d)]);
    }

    const auto [pdim, pid] = infer_background_edge_or_face_parent(old, va, vb);
    return append_generated_vertex<T>(
        vertex_x,
        vertex_ref_x,
        vertex_parent_dim,
        vertex_parent_id,
        vertex_root_edge_id,
        std::span<const T>(x_ref.data(), static_cast<std::size_t>(old.tdim)),
        std::span<const T>(x_phys.data(), static_cast<std::size_t>(old.gdim)),
        pdim,
        pid);
}

template <std::floating_point T>
int32_t append_average_vertex(
    const LocalMesh<T>&    old,
    std::span<const int32_t> verts,
    int32_t                parent_dim,
    int32_t                parent_id,
    std::vector<T>&        vertex_x,
    std::vector<T>&        vertex_ref_x,
    std::vector<int32_t>&  vertex_parent_dim,
    std::vector<int32_t>&  vertex_parent_id,
    std::vector<int32_t>&  vertex_root_edge_id)
{
    std::array<T, 3> x_ref = {T(0), T(0), T(0)};
    std::array<T, 3> x_phys = {T(0), T(0), T(0)};
    const T inv = T(1) / static_cast<T>(verts.size());
    for (const int32_t v : verts)
    {
        for (int d = 0; d < old.tdim; ++d)
            x_ref[static_cast<std::size_t>(d)]
                += inv * old.vertex_ref_x[static_cast<std::size_t>(v * old.tdim + d)];
        for (int d = 0; d < old.gdim; ++d)
            x_phys[static_cast<std::size_t>(d)]
                += inv * old.vertex_x[static_cast<std::size_t>(v * old.gdim + d)];
    }

    return append_generated_vertex<T>(
        vertex_x,
        vertex_ref_x,
        vertex_parent_dim,
        vertex_parent_id,
        vertex_root_edge_id,
        std::span<const T>(x_ref.data(), static_cast<std::size_t>(old.tdim)),
        std::span<const T>(x_phys.data(), static_cast<std::size_t>(old.gdim)),
        parent_dim,
        parent_id);
}

template <typename Table>
void append_subdivided_cells(
    const Table&               table,
    std::span<const int32_t>   local_to_global,
    cell::type                 child_cell_type,
    std::vector<int32_t>&      cell_vertices,
    std::vector<int32_t>&      cell_offsets,
    std::vector<cell::type>&   cell_types)
{
    for (const auto& child : table)
    {
        for (const int lid : child)
            cell_vertices.push_back(local_to_global[static_cast<std::size_t>(lid)]);
        cell_offsets.push_back(static_cast<int32_t>(cell_vertices.size()));
        cell_types.push_back(child_cell_type);
    }
}
struct FaceKey
{
    std::array<int32_t, 4> verts = {-1, -1, -1, -1};

    bool operator==(const FaceKey& other) const noexcept
    {
        return verts == other.verts;
    }
};

struct FaceKeyHash
{
    std::size_t operator()(const FaceKey& k) const noexcept
    {
        std::size_t h = 0;
        for (int32_t v : k.verts)
        {
            const auto hv = std::hash<int32_t>{}(v);
            h ^= hv + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

// ============================================================================
// LocalMesh Implementation
// ============================================================================

template <std::floating_point T>
void init_local_mesh_from_template(
    LocalMesh<T>&             mesh,
    const RefinementTemplate& tpl,
    std::span<const T>        parent_cell_coords,
    cell::type                parent_cell_type,
    int                       parent_cell_id,
    int                       n_level_sets)
{
    if (parent_cell_type != tpl.bg_cell_type)
        throw std::invalid_argument("init_local_mesh_from_template: parent_cell_type does not match template background cell type");

    const int nv = tpl.n_vertices;
    const int nc = tpl.n_cells;
    const int tdim = tpl.tdim;
    if (nv <= 0)
        throw std::invalid_argument("init_local_mesh_from_template: template has no vertices");
    if (parent_cell_coords.size() % static_cast<std::size_t>(nv) != 0)
        throw std::invalid_argument("init_local_mesh_from_template: parent_cell_coords size is not divisible by n_vertices");

    const int gdim = static_cast<int>(parent_cell_coords.size() / static_cast<std::size_t>(nv));
    if (gdim <= 0)
        throw std::invalid_argument("init_local_mesh_from_template: invalid geometric dimension");

    mesh.gdim = gdim;
    mesh.tdim = tdim;
    mesh.parent_cell_id = parent_cell_id;
    mesh.parent_cell_type = parent_cell_type;
    mesh.n_level_sets = n_level_sets;
    const int p1_nv = cell::get_num_vertices(parent_cell_type);
    mesh.parent_cell_coords_p1.resize(static_cast<std::size_t>(p1_nv * gdim));
    for (int i = 0; i < p1_nv; ++i)
    {
        for (int d = 0; d < gdim; ++d)
            mesh.parent_cell_coords_p1[static_cast<std::size_t>(i * gdim + d)] = parent_cell_coords[static_cast<std::size_t>(i * gdim + d)];
    }

    // Reset and resize vertices
    mesh.vertex_ref_x.resize(static_cast<std::size_t>(nv * tdim));
    mesh.vertex_x.resize(nv * gdim);
    mesh.vertex_parent_dim.assign(tpl.vertex_parent_dim.begin(), tpl.vertex_parent_dim.end());
    mesh.vertex_parent_id.assign(tpl.vertex_parent_id.begin(), tpl.vertex_parent_id.end());
    mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(nv), -1);

    const std::span<const double> ref_coords = tpl.ref_vertex_coords.empty()
        ? ((nv == cell::get_num_vertices(parent_cell_type))
               ? p1_ref_coords(parent_cell_type)
               : iso_p1_ref_coords(parent_cell_type, infer_lagrange_order(parent_cell_type, nv)))
        : std::span<const double>(tpl.ref_vertex_coords.data(), tpl.ref_vertex_coords.size());
    if (static_cast<int>(ref_coords.size()) != nv * tdim)
        throw std::invalid_argument("init_local_mesh_from_template: template reference coordinates have invalid size");

    for (int i = 0; i < nv; ++i)
    {
        for (int d = 0; d < tdim; ++d)
            mesh.vertex_ref_x[static_cast<std::size_t>(i * tdim + d)] = static_cast<T>(
                ref_coords[static_cast<std::size_t>(i * tdim + d)]);
        for (int d = 0; d < gdim; ++d)
            mesh.vertex_x[i * gdim + d] = parent_cell_coords[i * gdim + d];
    }

    // Build cells
    mesh.cell_offsets.resize(nc + 1);
    mesh.cell_vertices.assign(tpl.cell_connectivity.begin(), tpl.cell_connectivity.end());
    for (int i = 0; i <= nc; ++i) {
        mesh.cell_offsets[i] = i * tpl.vertices_per_cell;
    }
    mesh.cell_types.assign(nc, tpl.child_cell_type);
    mesh.cell_domain.assign(nc, static_cast<uint8_t>(cell::domain::unset));

    deduplicate_vertices_exact(mesh);

    // Classification masks
    const int nvd = mesh.n_vertices();
    mesh.vertex_zero_mask.assign(nvd, 0);
    mesh.vertex_inside_mask.assign(nvd, 0);
    mesh.vertex_phi.assign(static_cast<std::size_t>(nvd * n_level_sets), T(0));

    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void init_local_mesh_from_cell(
    LocalMesh<T>&      mesh,
    std::span<const T> parent_cell_coords,
    cell::type         parent_cell_type,
    int                parent_cell_id,
    int                n_level_sets)
{
    const auto [gdim, order] = infer_gdim_and_order(parent_cell_type, parent_cell_coords.size());

    const RefinementTemplate& meta_tpl = (order == 1) ? p1_template(parent_cell_type) : iso_p1_template(parent_cell_type, order);
    const RefinementTemplate& p1_tpl = p1_template(parent_cell_type);

    mesh.gdim = gdim;
    mesh.tdim = meta_tpl.tdim;
    mesh.parent_cell_id = parent_cell_id;
    mesh.parent_cell_type = parent_cell_type;
    mesh.n_level_sets = n_level_sets;

    mesh.vertex_x.assign(parent_cell_coords.begin(), parent_cell_coords.end());
    const int p1_nv = cell::get_num_vertices(parent_cell_type);
    mesh.parent_cell_coords_p1.resize(static_cast<std::size_t>(p1_nv * gdim));
    for (int i = 0; i < p1_nv; ++i)
    {
        for (int d = 0; d < gdim; ++d)
            mesh.parent_cell_coords_p1[static_cast<std::size_t>(i * gdim + d)] = parent_cell_coords[static_cast<std::size_t>(i * gdim + d)];
    }
    mesh.vertex_ref_x.resize(static_cast<std::size_t>(meta_tpl.n_vertices * meta_tpl.tdim));
    if (order == 1)
    {
        const auto xref = p1_ref_coords(parent_cell_type);
        for (std::size_t i = 0; i < xref.size(); ++i)
            mesh.vertex_ref_x[i] = static_cast<T>(xref[i]);
    }
    else
    {
        const auto xref = iso_p1_ref_coords(parent_cell_type, order);
        for (std::size_t i = 0; i < xref.size(); ++i)
            mesh.vertex_ref_x[i] = static_cast<T>(xref[i]);
    }

    mesh.vertex_parent_dim.assign(meta_tpl.vertex_parent_dim.begin(), meta_tpl.vertex_parent_dim.end());
    mesh.vertex_parent_id.assign(meta_tpl.vertex_parent_id.begin(), meta_tpl.vertex_parent_id.end());
    mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(meta_tpl.n_vertices), -1);
    mesh.cell_offsets = {0, p1_tpl.vertices_per_cell};
    mesh.cell_vertices.assign(p1_tpl.cell_connectivity.begin(), p1_tpl.cell_connectivity.end());
    mesh.cell_types = {parent_cell_type};
    mesh.cell_domain = {static_cast<uint8_t>(cell::domain::unset)};

    deduplicate_vertices_exact(mesh);

    const int nvd = mesh.n_vertices();
    mesh.vertex_zero_mask.assign(nvd, 0);
    mesh.vertex_inside_mask.assign(nvd, 0);
    mesh.vertex_phi.assign(static_cast<std::size_t>(nvd * n_level_sets), T(0));

    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void refine_local_mesh_from_template(
    LocalMesh<T>&             mesh,
    const RefinementTemplate& tpl)
{
    const cell::type parent_cell_type = mesh.parent_cell_type;
    const int parent_cell_id = mesh.parent_cell_id;
    const int n_level_sets = mesh.n_level_sets;
    const int p1_nv = cell::get_num_vertices(parent_cell_type);
    if (p1_nv <= 0)
        throw std::invalid_argument("refine_local_mesh_from_template: invalid parent cell type");

    // If current local-mesh vertices already match template-node count, reuse them
    // directly as parent-cell interpolation-node coordinates.
    if (mesh.gdim > 0
        && static_cast<int>(mesh.vertex_x.size()) == tpl.n_vertices * mesh.gdim)
    {
        const std::vector<T> current_coords = mesh.vertex_x;
        init_local_mesh_from_template(
            mesh,
            tpl,
            std::span<const T>(current_coords.data(), current_coords.size()),
            parent_cell_type,
            parent_cell_id,
            n_level_sets);
        return;
    }

    const std::span<const T> parent_cell_coords(mesh.parent_cell_coords_p1.data(), mesh.parent_cell_coords_p1.size());
    int gdim = -1;
    if (parent_cell_coords.size() % static_cast<std::size_t>(p1_nv) == 0)
        gdim = static_cast<int>(parent_cell_coords.size() / static_cast<std::size_t>(p1_nv));
    if (gdim <= 0)
        throw std::invalid_argument("refine_local_mesh_from_template: missing valid stored parent-cell P1 coordinates");

    std::vector<T> mapped_template_coords(static_cast<std::size_t>(tpl.n_vertices * gdim), T(0));
    std::vector<T> parent_p1(parent_cell_coords.begin(), parent_cell_coords.end());
    const std::span<const double> ref_coords = tpl.ref_vertex_coords.empty()
        ? ((tpl.n_vertices == cell::get_num_vertices(parent_cell_type))
               ? p1_ref_coords(parent_cell_type)
               : iso_p1_ref_coords(parent_cell_type, infer_lagrange_order(parent_cell_type, tpl.n_vertices)))
        : std::span<const double>(tpl.ref_vertex_coords.data(), tpl.ref_vertex_coords.size());
    if (static_cast<int>(ref_coords.size()) != tpl.n_vertices * tpl.tdim)
        throw std::invalid_argument("refine_local_mesh_from_template: template reference coordinates have invalid size");
    const std::vector<T> ref_coords_t(ref_coords.begin(), ref_coords.end());

    cell::push_forward_affine(
        parent_cell_type, parent_p1, gdim,
        std::span<const T>(ref_coords_t.data(), ref_coords_t.size()),
        std::span<T>(mapped_template_coords.data(), mapped_template_coords.size()));

    init_local_mesh_from_template(
        mesh,
        tpl,
        std::span<const T>(mapped_template_coords.data(), mapped_template_coords.size()),
        parent_cell_type,
        parent_cell_id,
        n_level_sets);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void red_refine_marked_cells(
    LocalMesh<T>&            mesh,
    std::span<const uint8_t> marked_cells)
{
    if (static_cast<int>(marked_cells.size()) != mesh.n_cells())
        throw std::invalid_argument("red_refine_marked_cells: marked_cells size mismatch");

    const LocalMesh<T> old = mesh;

    std::vector<int32_t> old_to_kept(static_cast<std::size_t>(old.n_vertices()), -1);
    std::vector<T> new_vertex_x;
    std::vector<T> new_vertex_ref_x;
    std::vector<int32_t> new_vertex_parent_dim;
    std::vector<int32_t> new_vertex_parent_id;
    std::vector<int32_t> new_vertex_root_edge_id;
    new_vertex_x.reserve(old.vertex_x.size());
    new_vertex_ref_x.reserve(old.vertex_ref_x.size());
    new_vertex_parent_dim.reserve(old.vertex_parent_dim.size());
    new_vertex_parent_id.reserve(old.vertex_parent_id.size());
    new_vertex_root_edge_id.reserve(old.vertex_root_edge_id.size());

    for (int32_t old_v = 0; old_v < old.n_vertices(); ++old_v)
    {
        if (old.vertex_root_edge_id[static_cast<std::size_t>(old_v)] >= 0)
            continue;
        const int32_t new_v = static_cast<int32_t>(new_vertex_parent_dim.size());
        old_to_kept[static_cast<std::size_t>(old_v)] = new_v;
        for (int d = 0; d < old.gdim; ++d)
            new_vertex_x.push_back(old.vertex_x[static_cast<std::size_t>(old_v * old.gdim + d)]);
        for (int d = 0; d < old.tdim; ++d)
            new_vertex_ref_x.push_back(old.vertex_ref_x[static_cast<std::size_t>(old_v * old.tdim + d)]);
        new_vertex_parent_dim.push_back(old.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
        new_vertex_parent_id.push_back(old.vertex_parent_id[static_cast<std::size_t>(old_v)]);
        // Drop all cached roots from rejected candidate cuts.
        new_vertex_root_edge_id.push_back(-1);
    }

    std::vector<int32_t> new_cell_vertices;
    std::vector<int32_t> new_cell_offsets = {0};
    std::vector<cell::type> new_cell_types;
    new_cell_vertices.reserve(old.cell_vertices.size() * 4);
    new_cell_types.reserve(old.n_cells() * 4);

    for (int c = 0; c < old.n_cells(); ++c)
    {
        const int c0 = old.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = old.cell_offsets[static_cast<std::size_t>(c + 1)];
        const std::span<const int32_t> old_cell_verts(
            old.cell_vertices.data() + static_cast<std::size_t>(c0),
            static_cast<std::size_t>(c1 - c0));
        const cell::type ct = old.cell_types[static_cast<std::size_t>(c)];
        std::vector<int32_t> cell_verts(static_cast<std::size_t>(old_cell_verts.size()), -1);
        for (std::size_t i = 0; i < old_cell_verts.size(); ++i)
        {
            const int32_t old_v = old_cell_verts[i];
            if (old_v < 0 || old_v >= old.n_vertices())
                throw std::invalid_argument("red_refine_marked_cells: invalid cell vertex id");
            cell_verts[i] = old_to_kept[static_cast<std::size_t>(old_v)];
            if (cell_verts[i] < 0)
                throw std::invalid_argument("red_refine_marked_cells: cell references dropped root vertex");
        }

        if (!marked_cells[static_cast<std::size_t>(c)])
        {
            new_cell_vertices.insert(new_cell_vertices.end(), cell_verts.begin(), cell_verts.end());
            new_cell_offsets.push_back(static_cast<int32_t>(new_cell_vertices.size()));
            new_cell_types.push_back(ct);
            continue;
        }

        switch (ct)
        {
        case cell::type::interval:
        {
            std::array<int32_t, 3> local = {
                cell_verts[0],
                cell_verts[1],
                append_segment_midpoint(
                    old, old_cell_verts[0], old_cell_verts[1],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id)};
            append_subdivided_cells(
                cell::interval_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::interval,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::triangle:
        {
            std::array<int32_t, 6> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], -1, -1, -1};
            const auto edges = cell::edges(cell::type::triangle);
            for (int e = 0; e < 3; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(3 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            append_subdivided_cells(
                cell::triangle_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::triangle,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::quadrilateral:
        {
            std::array<int32_t, 9> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3],
                -1, -1, -1, -1, -1};
            const auto edges = cell::edges(cell::type::quadrilateral);
            for (int e = 0; e < 4; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(4 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            local[8] = append_average_vertex(
                old,
                old_cell_verts,
                old.tdim,
                old.parent_cell_id,
                new_vertex_x, new_vertex_ref_x,
                new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            append_subdivided_cells(
                cell::quadrilateral_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::quadrilateral,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::tetrahedron:
        {
            std::array<int32_t, 10> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3],
                -1, -1, -1, -1, -1, -1};
            const auto edges = cell::edges(cell::type::tetrahedron);
            for (int e = 0; e < 6; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(4 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            append_subdivided_cells(
                cell::tetrahedron_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::tetrahedron,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        case cell::type::hexahedron:
        {
            std::array<int32_t, 27> local = {
                cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3],
                cell_verts[4], cell_verts[5], cell_verts[6], cell_verts[7],
                -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1, -1};
            const auto edges = cell::edges(cell::type::hexahedron);
            for (int e = 0; e < 12; ++e)
            {
                const int a = edges[static_cast<std::size_t>(e)][0];
                const int b = edges[static_cast<std::size_t>(e)][1];
                local[static_cast<std::size_t>(8 + e)] = append_segment_midpoint(
                    old, old_cell_verts[static_cast<std::size_t>(a)], old_cell_verts[static_cast<std::size_t>(b)],
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            for (int f = 0; f < 6; ++f)
            {
                const auto face = cell::hexahedron_faces[static_cast<std::size_t>(f)];
                std::array<int32_t, 4> face_verts = {
                    old_cell_verts[static_cast<std::size_t>(face[0])],
                    old_cell_verts[static_cast<std::size_t>(face[1])],
                    old_cell_verts[static_cast<std::size_t>(face[2])],
                    old_cell_verts[static_cast<std::size_t>(face[3])]};
                const auto [pdim, pid] = infer_background_face_parent(
                    old, std::span<const int32_t>(face_verts.data(), face_verts.size()));
                local[static_cast<std::size_t>(20 + f)] = append_average_vertex(
                    old,
                    std::span<const int32_t>(face_verts.data(), face_verts.size()),
                    pdim,
                    pid,
                    new_vertex_x, new_vertex_ref_x,
                    new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            }
            local[26] = append_average_vertex(
                old,
                old_cell_verts,
                old.tdim,
                old.parent_cell_id,
                new_vertex_x, new_vertex_ref_x,
                new_vertex_parent_dim, new_vertex_parent_id, new_vertex_root_edge_id);
            append_subdivided_cells(
                cell::hexahedron_subdivision_table,
                std::span<const int32_t>(local.data(), local.size()),
                cell::type::hexahedron,
                new_cell_vertices,
                new_cell_offsets,
                new_cell_types);
            break;
        }
        default:
            throw std::invalid_argument("red_refine_marked_cells: unsupported cell type");
        }
    }

    mesh.vertex_x = std::move(new_vertex_x);
    mesh.vertex_ref_x = std::move(new_vertex_ref_x);
    mesh.vertex_parent_dim = std::move(new_vertex_parent_dim);
    mesh.vertex_parent_id = std::move(new_vertex_parent_id);
    mesh.vertex_root_edge_id = std::move(new_vertex_root_edge_id);
    mesh.cell_vertices = std::move(new_cell_vertices);
    mesh.cell_offsets = std::move(new_cell_offsets);
    mesh.cell_types = std::move(new_cell_types);
    mesh.cell_domain.assign(static_cast<std::size_t>(mesh.cell_types.size()),
                            static_cast<uint8_t>(cell::domain::unset));

    deduplicate_vertices_exact(mesh);

    const int nv = mesh.n_vertices();
    mesh.vertex_zero_mask.assign(static_cast<std::size_t>(nv), 0);
    mesh.vertex_inside_mask.assign(static_cast<std::size_t>(nv), 0);
    mesh.vertex_phi.assign(static_cast<std::size_t>(nv * mesh.n_level_sets), T(0));

    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void build_local_faces(LocalMesh<T>& mesh)
{
    // No-op for non-volume meshes
    if (mesh.tdim < 3)
    {
        mesh.face_vertices.clear();
        mesh.face_offsets.clear();
        mesh.face_types.clear();
        mesh.face_parent_dim.clear();
        mesh.face_parent_id.clear();
        mesh.cell_face_offsets.clear();
        mesh.cell_faces_flat.clear();
        mesh.parent_face_to_local_face.clear();
        return;
    }

    const int nc = mesh.n_cells();
    std::unordered_map<FaceKey, int32_t, FaceKeyHash> face_to_id;
    face_to_id.reserve(static_cast<std::size_t>(nc * 8));

    std::vector<std::array<int32_t, 4>> unique_faces;
    std::vector<cell::type> unique_face_types;
    std::vector<int32_t>    unique_face_parent_dim;
    std::vector<int32_t>    unique_face_parent_id;

    mesh.cell_face_offsets.resize(static_cast<std::size_t>(nc + 1));
    mesh.cell_faces_flat.clear();

    // Determine number of background faces for parent_face_to_local_face
    const int n_bg_faces = cell::num_faces(mesh.parent_cell_type);
    mesh.parent_face_to_local_face.assign(static_cast<std::size_t>(std::max(0, n_bg_faces)), -1);

    for (int c = 0; c < nc; ++c)
    {
        mesh.cell_face_offsets[static_cast<std::size_t>(c)] =
            static_cast<int32_t>(mesh.cell_faces_flat.size());

        const int cv0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const cell::type ct = mesh.cell_types[static_cast<std::size_t>(c)];

        // Get face definitions for this cell type
        const int nf_cell = cell::num_faces(ct);
        if (nf_cell == 0)
            continue;

        const auto face_defs = cell::faces(ct);
        const auto face_sizes = cell::face_vertex_counts(ct);

        for (int f = 0; f < nf_cell; ++f)
        {
            const int fsz = face_sizes[static_cast<std::size_t>(f)];
            const auto& fdef = face_defs[static_cast<std::size_t>(f)];

            // Gather global vertex ids for this face
            std::array<int32_t, 4> gv = {-1, -1, -1, -1};
            for (int k = 0; k < fsz; ++k)
            {
                const int local_v = fdef[static_cast<std::size_t>(k)];
                gv[static_cast<std::size_t>(k)] =
                    mesh.cell_vertices[static_cast<std::size_t>(cv0 + local_v)];
            }

            // Canonical key: sorted vertex ids
            FaceKey key;
            key.verts = gv;
            std::sort(key.verts.begin(), key.verts.begin() + fsz);

            const auto it = face_to_id.find(key);
            int32_t face_idx = 0;
            if (it == face_to_id.end())
            {
                face_idx = static_cast<int32_t>(unique_faces.size());
                unique_faces.push_back(gv);

                // Face type
                const cell::type ftype = (fsz == 3) ? cell::type::triangle : cell::type::quadrilateral;
                unique_face_types.push_back(ftype);

                int32_t bg_face_dim = -1;
                int32_t bg_face_id = -1;
                const auto [pdim, pid] = infer_background_face_parent(
                    mesh,
                    std::span<const int32_t>(gv.data(), static_cast<std::size_t>(fsz)));
                bg_face_dim = pdim;
                bg_face_id = pid;

                unique_face_parent_dim.push_back(bg_face_dim);
                unique_face_parent_id.push_back(bg_face_id);

                if (bg_face_id >= 0 && bg_face_id < n_bg_faces)
                {
                    if (mesh.parent_face_to_local_face[static_cast<std::size_t>(bg_face_id)] < 0)
                        mesh.parent_face_to_local_face[static_cast<std::size_t>(bg_face_id)] = face_idx;
                }

                face_to_id.emplace(key, face_idx);
            }
            else
            {
                face_idx = it->second;
            }
            mesh.cell_faces_flat.push_back(face_idx);
        }
    }
    mesh.cell_face_offsets[static_cast<std::size_t>(nc)] =
        static_cast<int32_t>(mesh.cell_faces_flat.size());

    // Build CSR face arrays
    const int n_faces = static_cast<int>(unique_faces.size());
    mesh.face_offsets.resize(static_cast<std::size_t>(n_faces + 1));
    mesh.face_vertices.clear();
    mesh.face_offsets[0] = 0;
    for (int f = 0; f < n_faces; ++f)
    {
        const auto& gv = unique_faces[static_cast<std::size_t>(f)];
        const cell::type ftype = unique_face_types[static_cast<std::size_t>(f)];
        const int fsz = (ftype == cell::type::triangle) ? 3 : 4;
        for (int k = 0; k < fsz; ++k)
            mesh.face_vertices.push_back(gv[static_cast<std::size_t>(k)]);
        mesh.face_offsets[static_cast<std::size_t>(f + 1)] =
            static_cast<int32_t>(mesh.face_vertices.size());
    }
    mesh.face_types = std::move(unique_face_types);
    mesh.face_parent_dim = std::move(unique_face_parent_dim);
    mesh.face_parent_id = std::move(unique_face_parent_id);
}

// ============================================================================
// build_zero_entities / topology caches
// ============================================================================

namespace
{

template <std::floating_point T>
void clear_zero_entity_compat_views(LocalMesh<T>& mesh)
{
    mesh.iface_vertices.clear();
    mesh.iface_offsets.clear();
    mesh.iface_level_set_id.clear();
    mesh.iface_parent_cell.clear();
    mesh.curved_iface_ref_nodes.clear();
    mesh.curved_iface_offsets.clear();
    mesh.curved_iface_converged.clear();
    mesh.curved_iface_status.clear();
}

template <std::floating_point T>
void sync_iface_compatibility_view(LocalMesh<T>& mesh)
{
    clear_zero_entity_compat_views(mesh);
    mesh.iface_offsets.push_back(0);

    const int codim1_dim = std::max(0, mesh.tdim - 1);
    const int n_zero = mesh.n_zero_entities();
    for (int ze = 0; ze < n_zero; ++ze)
    {
        if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
            continue;
        const uint64_t zm = mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)];
        if (zm == 0 || (zm & (zm - 1)) != 0)
            continue;

        const int z0 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze)];
        const int z1 = mesh.zero_entity_offsets[static_cast<std::size_t>(ze + 1)];
        for (int i = z0; i < z1; ++i)
            mesh.iface_vertices.push_back(mesh.zero_entity_vertices[static_cast<std::size_t>(i)]);
        mesh.iface_offsets.push_back(static_cast<int32_t>(mesh.iface_vertices.size()));
        mesh.iface_level_set_id.push_back(static_cast<int32_t>(std::countr_zero(zm)));
        mesh.iface_parent_cell.push_back(mesh.zero_entity_parent_cell[static_cast<std::size_t>(ze)]);
    }
}

template <std::floating_point T>
void append_zero_entity(
    LocalMesh<T>&         mesh,
    uint8_t               dim,
    uint64_t              zero_mask,
    std::span<const int32_t> verts,
    int32_t               parent_cell,
    int32_t               parent_dim = -1,
    int32_t               parent_id = -1,
    std::span<const int32_t> face_edge_vertices = {})
{
    mesh.zero_entity_dim.push_back(dim);
    mesh.zero_entity_zero_mask.push_back(zero_mask);
    mesh.zero_entity_sign_mask.push_back(0);
    mesh.zero_entity_parent_cell.push_back(parent_cell);
    mesh.zero_entity_parent_dim.push_back(parent_dim);
    mesh.zero_entity_parent_id.push_back(parent_id);
    mesh.zero_entity_is_owned.push_back(uint8_t(1));

    if (mesh.zero_entity_offsets.empty())
        mesh.zero_entity_offsets.push_back(0);
    for (int32_t v : verts)
        mesh.zero_entity_vertices.push_back(v);
    mesh.zero_entity_offsets.push_back(static_cast<int32_t>(mesh.zero_entity_vertices.size()));

    if (dim == 1 && verts.size() >= 2)
    {
        mesh.zero_entity_endpoint_v0.push_back(verts[0]);
        mesh.zero_entity_endpoint_v1.push_back(verts[1]);
    }
    else
    {
        mesh.zero_entity_endpoint_v0.push_back(-1);
        mesh.zero_entity_endpoint_v1.push_back(-1);
    }

    if (mesh.zero_face_edge_offsets.empty())
        mesh.zero_face_edge_offsets.push_back(0);
    for (int32_t ev : face_edge_vertices)
        mesh.zero_face_edge_vertices.push_back(ev);
    mesh.zero_face_edge_offsets.push_back(static_cast<int32_t>(mesh.zero_face_edge_vertices.size()));
}

template <std::floating_point T>
void clear_zero_topology_caches(LocalMesh<T>& mesh)
{
    mesh.zero_entity_dim.clear();
    mesh.zero_entity_zero_mask.clear();
    mesh.zero_entity_sign_mask.clear();
    mesh.zero_entity_vertices.clear();
    mesh.zero_entity_offsets.clear();
    mesh.zero_entity_parent_cell.clear();
    mesh.zero_entity_parent_dim.clear();
    mesh.zero_entity_parent_id.clear();
    mesh.zero_entity_is_owned.clear();
    mesh.zero_entity_endpoint_v0.clear();
    mesh.zero_entity_endpoint_v1.clear();
    mesh.zero_face_edge_vertices.clear();
    mesh.zero_face_edge_offsets.clear();

    clear_zero_entity_compat_views(mesh);

    mesh.zero_chain_offsets.clear();
    mesh.zero_chain_entity_ids.clear();
    mesh.zero_chain_entity_reversed.clear();
    mesh.zero_chain_is_closed.clear();
    mesh.zero_chain_zero_mask.clear();

    mesh.zero_patch_offsets.clear();
    mesh.zero_patch_entity_ids.clear();
    mesh.zero_patch_entity_oriented.clear();
    mesh.zero_patch_is_closed.clear();
    mesh.zero_patch_zero_mask.clear();
}

template <std::floating_point T>
void sync_curved_zero_compatibility_view(LocalMesh<T>& mesh)
{
    mesh.curved_zero_offsets.assign(
        static_cast<std::size_t>(mesh.n_zero_entities() + 1), 0);
    mesh.curved_zero_ref_nodes.clear();
    mesh.curved_zero_converged.clear();
    mesh.curved_zero_status.clear();

    std::vector<int32_t> iface_to_zero;
    const int n_zero = mesh.n_zero_entities();
    const int codim1_dim = std::max(0, mesh.tdim - 1);
    for (int ze = 0; ze < n_zero; ++ze)
    {
        if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != codim1_dim)
            continue;
        const uint64_t zm = mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)];
        if (zm == 0 || (zm & (zm - 1)) != 0)
            continue;
        iface_to_zero.push_back(ze);
    }

    if (mesh.curved_iface_offsets.empty()
        || static_cast<int>(mesh.curved_iface_offsets.size()) != static_cast<int>(iface_to_zero.size()) + 1)
        return;

    for (std::size_t iface_id = 0; iface_id < iface_to_zero.size(); ++iface_id)
    {
        const int ze = iface_to_zero[iface_id];
        const int c0 = mesh.curved_iface_offsets[iface_id];
        const int c1 = mesh.curved_iface_offsets[iface_id + 1];
        const int n_nodes = c1 - c0;
        mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)] = n_nodes;
    }

    for (int ze = 0; ze < n_zero; ++ze)
        mesh.curved_zero_offsets[static_cast<std::size_t>(ze + 1)] += mesh.curved_zero_offsets[static_cast<std::size_t>(ze)];

    const int total_nodes = mesh.curved_zero_offsets.back();
    mesh.curved_zero_ref_nodes.resize(static_cast<std::size_t>(total_nodes * mesh.tdim), T(0));
    mesh.curved_zero_converged.resize(static_cast<std::size_t>(total_nodes), 0);
    mesh.curved_zero_status.resize(
        static_cast<std::size_t>(total_nodes),
        static_cast<uint8_t>(CurveStatus::failed));

    for (std::size_t iface_id = 0; iface_id < iface_to_zero.size(); ++iface_id)
    {
        const int ze = iface_to_zero[iface_id];
        const int src0 = mesh.curved_iface_offsets[iface_id];
        const int src1 = mesh.curved_iface_offsets[iface_id + 1];
        const int dst0 = mesh.curved_zero_offsets[static_cast<std::size_t>(ze)];
        const int n_nodes = src1 - src0;

        for (int i = 0; i < n_nodes; ++i)
        {
            mesh.curved_zero_converged[static_cast<std::size_t>(dst0 + i)]
                = mesh.curved_iface_converged[static_cast<std::size_t>(src0 + i)];
            if (!mesh.curved_iface_status.empty())
            {
                mesh.curved_zero_status[static_cast<std::size_t>(dst0 + i)]
                    = mesh.curved_iface_status[static_cast<std::size_t>(src0 + i)];
            }
            for (int d = 0; d < mesh.tdim; ++d)
            {
                mesh.curved_zero_ref_nodes[static_cast<std::size_t>((dst0 + i) * mesh.tdim + d)]
                    = mesh.curved_iface_ref_nodes[static_cast<std::size_t>((src0 + i) * mesh.tdim + d)];
            }
        }
    }
}

} // namespace

template <std::floating_point T>
void build_zero_entities(LocalMesh<T>& mesh, int level_set_id)
{
    const uint64_t zero_mask = uint64_t(1) << level_set_id;

    if (mesh.zero_entity_offsets.empty())
        mesh.zero_entity_offsets.push_back(0);
    if (mesh.zero_face_edge_offsets.empty())
        mesh.zero_face_edge_offsets.push_back(0);

    clear_zero_topology_caches(mesh);

    constexpr uint8_t D_INSIDE  = static_cast<uint8_t>(cell::domain::inside);
    constexpr uint8_t D_OUTSIDE = static_cast<uint8_t>(cell::domain::outside);

    if (mesh.tdim == 2)
    {
        const int ne = mesh.n_edges();
        const int nc = mesh.n_cells();
        std::vector<int32_t> edge_adj_cell0(static_cast<std::size_t>(ne), -1);
        std::vector<int32_t> edge_adj_cell1(static_cast<std::size_t>(ne), -1);
        std::vector<uint8_t> edge_added(static_cast<std::size_t>(ne), uint8_t(0));

        for (int c = 0; c < nc; ++c)
        {
            const int e_start = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
            const int e_end   = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
            for (int ei = e_start; ei < e_end; ++ei)
            {
                const int eid = mesh.cell_edges_flat[static_cast<std::size_t>(ei)];
                if (edge_adj_cell0[static_cast<std::size_t>(eid)] < 0)
                    edge_adj_cell0[static_cast<std::size_t>(eid)] = static_cast<int32_t>(c);
                else
                    edge_adj_cell1[static_cast<std::size_t>(eid)] = static_cast<int32_t>(c);
            }
        }

        // Touch-only zero edges: keep them as interface entities directly.
        for (int e = 0; e < ne; ++e)
        {
            const auto st = static_cast<EdgeState>(mesh.edge_state_for(e, level_set_id));
            if (st != EdgeState::on_interface)
                continue;

            const std::array<int32_t, 2> verts = {
                mesh.edge_vertices[static_cast<std::size_t>(2 * e)],
                mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)]};
            const int32_t owner_cell = (edge_adj_cell0[static_cast<std::size_t>(e)] >= 0)
                ? edge_adj_cell0[static_cast<std::size_t>(e)]
                : edge_adj_cell1[static_cast<std::size_t>(e)];
            append_zero_entity(mesh, 1, zero_mask,
                               std::span<const int32_t>(verts.data(), verts.size()),
                               owner_cell,
                               mesh.edge_parent_dim[static_cast<std::size_t>(e)],
                               mesh.edge_parent_id[static_cast<std::size_t>(e)]);
            edge_added[static_cast<std::size_t>(e)] = uint8_t(1);
        }

        for (int e = 0; e < ne; ++e)
        {
            if (edge_added[static_cast<std::size_t>(e)] != 0)
                continue;
            const int32_t c0 = edge_adj_cell0[static_cast<std::size_t>(e)];
            const int32_t c1 = edge_adj_cell1[static_cast<std::size_t>(e)];
            if (c0 < 0 || c1 < 0)
                continue;

            const uint8_t d0 = mesh.cell_domain[static_cast<std::size_t>(c0)];
            const uint8_t d1 = mesh.cell_domain[static_cast<std::size_t>(c1)];
            if (!((d0 == D_INSIDE && d1 == D_OUTSIDE) ||
                  (d0 == D_OUTSIDE && d1 == D_INSIDE)))
                continue;

            const std::array<int32_t, 2> verts = {
                mesh.edge_vertices[static_cast<std::size_t>(2 * e)],
                mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)]};
            const int32_t parent_c = (d0 == D_INSIDE) ? c0 : c1;
            const int32_t carrier_dim = mesh.edge_parent_dim[static_cast<std::size_t>(e)];
            const int32_t carrier_id = mesh.edge_parent_id[static_cast<std::size_t>(e)];
            append_zero_entity(mesh, 1, zero_mask,
                               std::span<const int32_t>(verts.data(), verts.size()),
                               parent_c,
                               carrier_dim,
                               carrier_id);
            edge_added[static_cast<std::size_t>(e)] = uint8_t(1);
        }
    }
    else if (mesh.tdim == 3)
    {
        const int nf = mesh.n_faces();
        const int nc = mesh.n_cells();
        std::vector<int32_t> face_adj_cell0(static_cast<std::size_t>(nf), -1);
        std::vector<int32_t> face_adj_cell1(static_cast<std::size_t>(nf), -1);
        std::vector<uint8_t> face_added(static_cast<std::size_t>(nf), uint8_t(0));
        std::unordered_map<EdgeKey, int32_t, EdgeKeyHash> zero_edge_entity_by_key;

        auto normalize_edge = [](int32_t a, int32_t b) -> EdgeKey
        {
            return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
        };

        auto append_zero_edge_if_missing = [&](int32_t a, int32_t b, int32_t owner_cell)
        {
            const EdgeKey key = normalize_edge(a, b);
            if (zero_edge_entity_by_key.find(key) != zero_edge_entity_by_key.end())
                return;

            int32_t parent_dim = -1;
            int32_t parent_id = -1;
            const int edge_id = find_local_edge_by_vertices(mesh, a, b);
            if (edge_id >= 0)
            {
                parent_dim = mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)];
                parent_id = mesh.edge_parent_id[static_cast<std::size_t>(edge_id)];
            }

            const std::array<int32_t, 2> verts = {a, b};
            append_zero_entity(
                mesh, 1, zero_mask,
                std::span<const int32_t>(verts.data(), verts.size()),
                owner_cell,
                parent_dim,
                parent_id);
            zero_edge_entity_by_key.emplace(
                key,
                static_cast<int32_t>(mesh.n_zero_entities() - 1));
        };

        for (int c = 0; c < nc; ++c)
        {
            const int f_start = mesh.cell_face_offsets[static_cast<std::size_t>(c)];
            const int f_end   = mesh.cell_face_offsets[static_cast<std::size_t>(c + 1)];
            for (int fi = f_start; fi < f_end; ++fi)
            {
                const int fid = mesh.cell_faces_flat[static_cast<std::size_t>(fi)];
                if (face_adj_cell0[static_cast<std::size_t>(fid)] < 0)
                    face_adj_cell0[static_cast<std::size_t>(fid)] = static_cast<int32_t>(c);
                else
                    face_adj_cell1[static_cast<std::size_t>(fid)] = static_cast<int32_t>(c);
            }
        }

        // Touch-only zero faces: if all face vertices are zero for this level set,
        // keep the face directly as a zero entity without volume cutting.
        for (int f = 0; f < nf; ++f)
        {
            const int f0 = mesh.face_offsets[static_cast<std::size_t>(f)];
            const int f1 = mesh.face_offsets[static_cast<std::size_t>(f + 1)];
            const int n_face_verts = f1 - f0;
            bool all_zero = (n_face_verts > 0);
            std::vector<int32_t> verts(static_cast<std::size_t>(n_face_verts));
            for (int i = 0; i < n_face_verts; ++i)
            {
                const int32_t lv = mesh.face_vertices[static_cast<std::size_t>(f0 + i)];
                verts[static_cast<std::size_t>(i)] = lv;
                if ((mesh.vertex_zero_mask[static_cast<std::size_t>(lv)] & zero_mask) == 0)
                    all_zero = false;
            }
            if (!all_zero)
                continue;

            const int32_t owner_cell = (face_adj_cell0[static_cast<std::size_t>(f)] >= 0)
                ? face_adj_cell0[static_cast<std::size_t>(f)]
                : face_adj_cell1[static_cast<std::size_t>(f)];
            std::vector<int32_t> edge_pairs;
            edge_pairs.reserve(static_cast<std::size_t>(2 * n_face_verts));
            for (int i = 0; i < n_face_verts; ++i)
            {
                const int32_t a = verts[static_cast<std::size_t>(i)];
                const int32_t b = verts[static_cast<std::size_t>((i + 1) % n_face_verts)];
                append_zero_edge_if_missing(a, b, owner_cell);
                edge_pairs.push_back(a);
                edge_pairs.push_back(b);
            }
            append_zero_entity(
                mesh, 2, zero_mask,
                std::span<const int32_t>(verts.data(), verts.size()),
                owner_cell,
                mesh.face_parent_dim[static_cast<std::size_t>(f)],
                mesh.face_parent_id[static_cast<std::size_t>(f)],
                std::span<const int32_t>(edge_pairs.data(), edge_pairs.size()));
            face_added[static_cast<std::size_t>(f)] = uint8_t(1);
        }

        for (int f = 0; f < nf; ++f)
        {
            if (face_added[static_cast<std::size_t>(f)] != 0)
                continue;
            const int32_t c0 = face_adj_cell0[static_cast<std::size_t>(f)];
            const int32_t c1 = face_adj_cell1[static_cast<std::size_t>(f)];
            if (c0 < 0 || c1 < 0)
                continue;

            const uint8_t d0 = mesh.cell_domain[static_cast<std::size_t>(c0)];
            const uint8_t d1 = mesh.cell_domain[static_cast<std::size_t>(c1)];
            if (!((d0 == D_INSIDE && d1 == D_OUTSIDE) ||
                  (d0 == D_OUTSIDE && d1 == D_INSIDE)))
                continue;

            const int f0 = mesh.face_offsets[static_cast<std::size_t>(f)];
            const int f1 = mesh.face_offsets[static_cast<std::size_t>(f + 1)];
            const int n_face_verts = f1 - f0;
            std::vector<int32_t> verts(static_cast<std::size_t>(n_face_verts));
            for (int i = 0; i < n_face_verts; ++i)
                verts[static_cast<std::size_t>(i)] =
                    mesh.face_vertices[static_cast<std::size_t>(f0 + i)];

            const int32_t parent_c = (d0 == D_INSIDE) ? c0 : c1;
            std::vector<int32_t> edge_pairs;
            edge_pairs.reserve(static_cast<std::size_t>(2 * n_face_verts));
            for (int i = 0; i < n_face_verts; ++i)
            {
                const int32_t a = verts[static_cast<std::size_t>(i)];
                const int32_t b = verts[static_cast<std::size_t>((i + 1) % n_face_verts)];
                append_zero_edge_if_missing(a, b, parent_c);
                edge_pairs.push_back(a);
                edge_pairs.push_back(b);
            }
            append_zero_entity(
                mesh, 2, zero_mask,
                std::span<const int32_t>(verts.data(), verts.size()),
                parent_c,
                mesh.face_parent_dim[static_cast<std::size_t>(f)],
                mesh.face_parent_id[static_cast<std::size_t>(f)],
                std::span<const int32_t>(edge_pairs.data(), edge_pairs.size()));
            face_added[static_cast<std::size_t>(f)] = uint8_t(1);
        }
    }

    sync_iface_compatibility_view(mesh);
    ++mesh.zero_entity_version;
    mesh.curved_cache_version = 0;
    mesh.curved_cache_zero_version = 0;
}

template <std::floating_point T>
void build_interface_entities(LocalMesh<T>& mesh, int level_set_id)
{
    build_zero_entities(mesh, level_set_id);
}

template <std::floating_point T>
void assign_zero_entity_ownership(
    std::span<LocalMesh<T>*> local_meshes,
    uint64_t                 zero_mask)
{
    struct Candidate
    {
        LocalMesh<T>* mesh = nullptr;
        int entity_id = -1;
        int parent_cell_id = -1;
    };

    std::unordered_map<ZeroOwnerKey, std::vector<Candidate>, ZeroOwnerKeyHash> groups;

    for (LocalMesh<T>* mesh : local_meshes)
    {
        if (mesh == nullptr)
            continue;
        const int n_zero = mesh->n_zero_entities();
        for (int ze = 0; ze < n_zero; ++ze)
        {
            if (mesh->zero_entity_zero_mask[static_cast<std::size_t>(ze)] != zero_mask)
                continue;

            const int32_t parent_dim = mesh->zero_entity_parent_dim[static_cast<std::size_t>(ze)];
            const int32_t parent_id = mesh->zero_entity_parent_id[static_cast<std::size_t>(ze)];
            mesh->zero_entity_is_owned[static_cast<std::size_t>(ze)] = uint8_t(1);

            if (parent_dim < 0 || parent_id < 0)
                continue;

            const ZeroOwnerKey key{
                zero_mask,
                mesh->zero_entity_dim[static_cast<std::size_t>(ze)],
                parent_dim,
                parent_id};
            groups[key].push_back(Candidate{mesh, ze, mesh->parent_cell_id});
        }
    }

    for (auto& [key, candidates] : groups)
    {
        (void)key;
        if (candidates.empty())
            continue;

        Candidate* owner = &candidates.front();
        for (auto& cand : candidates)
        {
            if (cand.parent_cell_id < owner->parent_cell_id)
                owner = &cand;
        }

        for (auto& cand : candidates)
        {
            cand.mesh->zero_entity_is_owned[static_cast<std::size_t>(cand.entity_id)]
                = (&cand == owner) ? uint8_t(1) : uint8_t(0);
        }
    }
}

template <std::floating_point T>
void build_zero_chains(LocalMesh<T>& mesh, uint64_t zero_mask)
{
    mesh.zero_chain_offsets.clear();
    mesh.zero_chain_entity_ids.clear();
    mesh.zero_chain_entity_reversed.clear();
    mesh.zero_chain_is_closed.clear();
    mesh.zero_chain_zero_mask.clear();
    mesh.zero_chain_offsets.push_back(0);

    const int n_zero = mesh.n_zero_entities();
    std::vector<int32_t> entities;
    entities.reserve(static_cast<std::size_t>(n_zero));
    std::unordered_map<int32_t, std::vector<int32_t>> endpoint_to_entities;

    for (int ze = 0; ze < n_zero; ++ze)
    {
        if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 1)
            continue;
        if (mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] != zero_mask)
            continue;
        const int32_t v0 = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
        const int32_t v1 = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
        if (v0 < 0 || v1 < 0)
            continue;
        entities.push_back(ze);
        endpoint_to_entities[v0].push_back(ze);
        endpoint_to_entities[v1].push_back(ze);
    }

    for (const auto& [vertex, incident] : endpoint_to_entities)
    {
        (void)vertex;
        if (incident.size() > 2)
            throw std::runtime_error("build_zero_chains: endpoint degree > 2");
    }

    std::unordered_map<int32_t, uint8_t> visited;
    auto other_endpoint = [&](int32_t ze, int32_t v) -> int32_t
    {
        const int32_t a = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
        const int32_t b = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(ze)];
        return (a == v) ? b : a;
    };

    auto append_chain = [&](int32_t start_entity, int32_t start_vertex, bool closed)
    {
        int32_t current_entity = start_entity;
        int32_t current_vertex = start_vertex;

        while (true)
        {
            if (visited[current_entity])
                break;
            visited[current_entity] = 1;

            const int32_t a = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(current_entity)];
            const int32_t b = mesh.zero_entity_endpoint_v1[static_cast<std::size_t>(current_entity)];
            const bool reversed = (current_vertex == b);
            mesh.zero_chain_entity_ids.push_back(current_entity);
            mesh.zero_chain_entity_reversed.push_back(reversed ? uint8_t(1) : uint8_t(0));

            const int32_t next_vertex = reversed ? a : b;
            const auto it = endpoint_to_entities.find(next_vertex);
            if (it == endpoint_to_entities.end())
                break;

            int32_t next_entity = -1;
            for (const int32_t cand : it->second)
            {
                if (!visited[cand])
                {
                    next_entity = cand;
                    break;
                }
            }

            if (next_entity < 0)
                break;

            current_entity = next_entity;
            current_vertex = next_vertex;
        }

        mesh.zero_chain_offsets.push_back(static_cast<int32_t>(mesh.zero_chain_entity_ids.size()));
        mesh.zero_chain_is_closed.push_back(closed ? uint8_t(1) : uint8_t(0));
        mesh.zero_chain_zero_mask.push_back(zero_mask);
    };

    for (const auto& [vertex, incident] : endpoint_to_entities)
    {
        if (incident.size() != 1)
            continue;
        const int32_t ze = incident.front();
        if (visited[ze])
            continue;
        append_chain(ze, vertex, false);
    }

    for (const int32_t ze : entities)
    {
        if (visited[ze])
            continue;
        const int32_t start_v = mesh.zero_entity_endpoint_v0[static_cast<std::size_t>(ze)];
        append_chain(ze, start_v, true);
    }
}

template <std::floating_point T>
void build_zero_patches(LocalMesh<T>& mesh, uint64_t zero_mask)
{
    mesh.zero_patch_offsets.clear();
    mesh.zero_patch_entity_ids.clear();
    mesh.zero_patch_entity_oriented.clear();
    mesh.zero_patch_is_closed.clear();
    mesh.zero_patch_zero_mask.clear();
    mesh.zero_patch_offsets.push_back(0);

    const int n_zero = mesh.n_zero_entities();
    std::vector<int32_t> faces;
    std::unordered_map<EdgeKey, std::vector<int32_t>, EdgeKeyHash> edge_to_faces;

    auto normalize_edge = [](int32_t a, int32_t b) -> EdgeKey
    {
        return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
    };

    for (int ze = 0; ze < n_zero; ++ze)
    {
        if (mesh.zero_entity_dim[static_cast<std::size_t>(ze)] != 2)
            continue;
        if (mesh.zero_entity_zero_mask[static_cast<std::size_t>(ze)] != zero_mask)
            continue;
        faces.push_back(ze);

        const int e0 = mesh.zero_face_edge_offsets[static_cast<std::size_t>(ze)];
        const int e1 = mesh.zero_face_edge_offsets[static_cast<std::size_t>(ze + 1)];
        for (int i = e0; i + 1 < e1; i += 2)
        {
            const EdgeKey key = normalize_edge(
                mesh.zero_face_edge_vertices[static_cast<std::size_t>(i)],
                mesh.zero_face_edge_vertices[static_cast<std::size_t>(i + 1)]);
            edge_to_faces[key].push_back(ze);
        }
    }

    std::unordered_map<int32_t, std::vector<int32_t>> face_adj;
    std::unordered_map<int32_t, uint8_t> face_boundary;
    for (const auto& [edge, incident] : edge_to_faces)
    {
        (void)edge;
        if (incident.size() == 1)
        {
            face_boundary[incident.front()] = 1;
            continue;
        }

        for (std::size_t i = 0; i < incident.size(); ++i)
        {
            const int32_t fi = incident[i];
            auto& adj = face_adj[fi];
            for (std::size_t j = 0; j < incident.size(); ++j)
            {
                if (i == j)
                    continue;
                const int32_t fj = incident[j];
                if (std::find(adj.begin(), adj.end(), fj) == adj.end())
                    adj.push_back(fj);
            }
        }
    }

    std::unordered_map<int32_t, uint8_t> visited;
    for (const int32_t ze : faces)
    {
        if (visited[ze])
            continue;

        std::vector<int32_t> stack = {ze};
        bool closed = true;
        while (!stack.empty())
        {
            const int32_t cur = stack.back();
            stack.pop_back();
            if (visited[cur])
                continue;
            visited[cur] = 1;

            mesh.zero_patch_entity_ids.push_back(cur);
            mesh.zero_patch_entity_oriented.push_back(0);
            if (face_boundary[cur])
                closed = false;

            for (const int32_t nb : face_adj[cur])
            {
                if (!visited[nb])
                    stack.push_back(nb);
            }
        }

        mesh.zero_patch_offsets.push_back(static_cast<int32_t>(mesh.zero_patch_entity_ids.size()));
        mesh.zero_patch_is_closed.push_back(closed ? uint8_t(1) : uint8_t(0));
        mesh.zero_patch_zero_mask.push_back(zero_mask);
    }
}

template <std::floating_point T>
void evaluate_levelsets_on_vertices(
    LocalMesh<T>&                           mesh,
    const std::vector<LevelSetFunction<T>>& phi,
    int                                     n_ls,
    T                                       tol)
{
    int nv = mesh.n_vertices();
    int gdim = mesh.gdim;
    mesh.vertex_phi.assign(nv * n_ls, 0.0);
    mesh.vertex_zero_mask.assign(nv, 0);
    mesh.vertex_inside_mask.assign(nv, 0);

    for (int i = 0; i < nv; ++i) {
        const T* x = &mesh.vertex_x[i * gdim];
        for (int ls = 0; ls < n_ls; ++ls) {
            T val = 0.0;
            val = phi[ls].value(x, mesh.parent_cell_id);
            
            mesh.vertex_phi[i * n_ls + ls] = val;
            if (std::abs(val) < tol) {
                mesh.vertex_zero_mask[i] |= (1ULL << ls);
            } else if (val < 0) {
                mesh.vertex_inside_mask[i] |= (1ULL << ls);
            }
        }
    }
}

template <std::floating_point T>
void classify_local_edges(LocalMesh<T>& mesh, int level_set_id)
{
    int ne = mesh.n_edges();
    uint64_t mask = (1ULL << level_set_id);

    for (int i = 0; i < ne; ++i) {
        int v0 = mesh.edge_vertices[i * 2];
        int v1 = mesh.edge_vertices[i * 2 + 1];

        bool zero0 = (mesh.vertex_zero_mask[v0] & mask);
        bool zero1 = (mesh.vertex_zero_mask[v1] & mask);
        bool inside0 = (mesh.vertex_inside_mask[v0] & mask);
        bool inside1 = (mesh.vertex_inside_mask[v1] & mask);

        if (zero0 && zero1) {
            mesh.edge_state_for(i, level_set_id) = static_cast<uint8_t>(EdgeState::zero_edge);
        } else if (inside0 == inside1) {
            mesh.edge_state_for(i, level_set_id) = static_cast<uint8_t>(EdgeState::uncut);
        } else {
            mesh.edge_state_for(i, level_set_id) = static_cast<uint8_t>(EdgeState::single_cross);
        }
        mesh.edge_cert_for(i, level_set_id) = static_cast<uint8_t>(EdgeCertification::certified);
    }
}

template <std::floating_point T>
void clear_edge_root_cache(LocalMesh<T>& mesh)
{
    std::fill(mesh.edge_root_vertex.begin(), mesh.edge_root_vertex.end(), -1);
    std::fill(mesh.edge_root_parameter.begin(), mesh.edge_root_parameter.end(), T(-1));
    std::fill(mesh.edge_root_iterations.begin(), mesh.edge_root_iterations.end(), 0);
    std::fill(mesh.edge_root_evaluations.begin(), mesh.edge_root_evaluations.end(), 0);
    std::fill(mesh.edge_root_converged.begin(), mesh.edge_root_converged.end(), uint8_t(0));
    std::fill(mesh.edge_root_residual.begin(), mesh.edge_root_residual.end(), T(0));
}

template <std::floating_point T>
void drop_unreferenced_root_vertices(LocalMesh<T>& mesh)
{
    const int old_nv = mesh.n_vertices();
    if (old_nv == 0 || static_cast<int>(mesh.vertex_root_edge_id.size()) != old_nv)
        return;

    std::vector<uint8_t> referenced(static_cast<std::size_t>(old_nv), uint8_t(0));
    for (const int32_t v : mesh.cell_vertices)
    {
        if (v >= 0 && v < old_nv)
            referenced[static_cast<std::size_t>(v)] = 1;
    }

    std::vector<int32_t> old_to_new(static_cast<std::size_t>(old_nv), -1);
    std::vector<T> new_x;
    std::vector<T> new_ref_x;
    std::vector<int32_t> new_parent_dim;
    std::vector<int32_t> new_parent_id;
    std::vector<int32_t> new_root_edge_id;
    std::vector<T> new_phi;
    std::vector<uint64_t> new_zero_mask;
    std::vector<uint64_t> new_inside_mask;

    new_x.reserve(mesh.vertex_x.size());
    new_ref_x.reserve(mesh.vertex_ref_x.size());
    new_parent_dim.reserve(mesh.vertex_parent_dim.size());
    new_parent_id.reserve(mesh.vertex_parent_id.size());
    new_root_edge_id.reserve(mesh.vertex_root_edge_id.size());
    new_phi.reserve(mesh.vertex_phi.size());
    new_zero_mask.reserve(mesh.vertex_zero_mask.size());
    new_inside_mask.reserve(mesh.vertex_inside_mask.size());

    int new_nv = 0;
    for (int old_v = 0; old_v < old_nv; ++old_v)
    {
        const bool is_root = mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)] >= 0;
        if (is_root && !referenced[static_cast<std::size_t>(old_v)])
            continue;

        old_to_new[static_cast<std::size_t>(old_v)] = new_nv++;
        for (int d = 0; d < mesh.gdim; ++d)
        {
            new_x.push_back(
                mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)]);
        }
        if (!mesh.vertex_ref_x.empty())
        {
            for (int d = 0; d < mesh.tdim; ++d)
            {
                new_ref_x.push_back(
                    mesh.vertex_ref_x[static_cast<std::size_t>(old_v * mesh.tdim + d)]);
            }
        }
        new_parent_dim.push_back(mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
        new_parent_id.push_back(mesh.vertex_parent_id[static_cast<std::size_t>(old_v)]);
        new_root_edge_id.push_back(mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)]);
        if (!mesh.vertex_phi.empty())
        {
            for (int ls = 0; ls < mesh.n_level_sets; ++ls)
            {
                new_phi.push_back(
                    mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + ls)]);
            }
        }
        if (!mesh.vertex_zero_mask.empty())
            new_zero_mask.push_back(mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)]);
        if (!mesh.vertex_inside_mask.empty())
            new_inside_mask.push_back(mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)]);
    }

    if (new_nv == old_nv)
        return;

    for (int32_t& v : mesh.cell_vertices)
        v = old_to_new[static_cast<std::size_t>(v)];

    for (int32_t& v : mesh.edge_vertices)
        v = old_to_new[static_cast<std::size_t>(v)];

    mesh.vertex_x.swap(new_x);
    mesh.vertex_ref_x.swap(new_ref_x);
    mesh.vertex_parent_dim.swap(new_parent_dim);
    mesh.vertex_parent_id.swap(new_parent_id);
    mesh.vertex_root_edge_id.swap(new_root_edge_id);
    mesh.vertex_phi.swap(new_phi);
    mesh.vertex_zero_mask.swap(new_zero_mask);
    mesh.vertex_inside_mask.swap(new_inside_mask);
}

template <std::floating_point T>
void deduplicate_zero_vertices(LocalMesh<T>& mesh,
                               int           level_set_id,
                               T             coord_tol,
                               T)
{
    const int old_nv = mesh.n_vertices();
    if (old_nv <= 1
        || static_cast<int>(mesh.vertex_phi.size()) != old_nv * mesh.n_level_sets
        || level_set_id < 0
        || level_set_id >= mesh.n_level_sets)
    {
        return;
    }

    std::vector<int32_t> old_to_new(static_cast<std::size_t>(old_nv), -1);
    std::vector<T> new_x;
    std::vector<T> new_ref_x;
    std::vector<int32_t> new_parent_dim;
    std::vector<int32_t> new_parent_id;
    std::vector<int32_t> new_root_edge_id;
    std::vector<T> new_phi;
    std::vector<uint64_t> new_zero_mask;
    std::vector<uint64_t> new_inside_mask;

    new_x.reserve(mesh.vertex_x.size());
    new_ref_x.reserve(mesh.vertex_ref_x.size());
    new_parent_dim.reserve(mesh.vertex_parent_dim.size());
    new_parent_id.reserve(mesh.vertex_parent_id.size());
    new_root_edge_id.reserve(mesh.vertex_root_edge_id.size());
    new_phi.reserve(mesh.vertex_phi.size());
    new_zero_mask.reserve(mesh.vertex_zero_mask.size());
    new_inside_mask.reserve(mesh.vertex_inside_mask.size());

    const auto append_vertex = [&](int old_v)
    {
        for (int d = 0; d < mesh.gdim; ++d)
        {
            new_x.push_back(mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)]);
        }
        if (!mesh.vertex_ref_x.empty())
        {
            for (int d = 0; d < mesh.tdim; ++d)
            {
                new_ref_x.push_back(
                    mesh.vertex_ref_x[static_cast<std::size_t>(old_v * mesh.tdim + d)]);
            }
        }
        new_parent_dim.push_back(mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)]);
        new_parent_id.push_back(mesh.vertex_parent_id[static_cast<std::size_t>(old_v)]);
        new_root_edge_id.push_back(mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)]);
        for (int ls = 0; ls < mesh.n_level_sets; ++ls)
        {
            new_phi.push_back(mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + ls)]);
        }
        new_zero_mask.push_back(mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)]);
        new_inside_mask.push_back(mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)]);
    };

    auto overwrite_vertex = [&](int new_v, int old_v)
    {
        for (int d = 0; d < mesh.gdim; ++d)
        {
            new_x[static_cast<std::size_t>(new_v * mesh.gdim + d)]
                = mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)];
        }
        if (!mesh.vertex_ref_x.empty())
        {
            for (int d = 0; d < mesh.tdim; ++d)
            {
                new_ref_x[static_cast<std::size_t>(new_v * mesh.tdim + d)]
                    = mesh.vertex_ref_x[static_cast<std::size_t>(old_v * mesh.tdim + d)];
            }
        }
        new_parent_dim[static_cast<std::size_t>(new_v)]
            = mesh.vertex_parent_dim[static_cast<std::size_t>(old_v)];
        new_parent_id[static_cast<std::size_t>(new_v)]
            = mesh.vertex_parent_id[static_cast<std::size_t>(old_v)];
        new_root_edge_id[static_cast<std::size_t>(new_v)]
            = mesh.vertex_root_edge_id[static_cast<std::size_t>(old_v)];
        for (int ls = 0; ls < mesh.n_level_sets; ++ls)
        {
            new_phi[static_cast<std::size_t>(new_v * mesh.n_level_sets + ls)]
                = mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + ls)];
        }
        new_zero_mask[static_cast<std::size_t>(new_v)]
            = mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)];
        new_inside_mask[static_cast<std::size_t>(new_v)]
            = mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)];
    };

    int new_nv = 0;
    for (int old_v = 0; old_v < old_nv; ++old_v)
    {
        int found = -1;
        for (int new_v = 0; new_v < new_nv; ++new_v)
        {
            bool same = true;
            for (int d = 0; d < mesh.gdim; ++d)
            {
                const T dv
                    = new_x[static_cast<std::size_t>(new_v * mesh.gdim + d)]
                      - mesh.vertex_x[static_cast<std::size_t>(old_v * mesh.gdim + d)];
                if (std::abs(dv) > coord_tol)
                {
                    same = false;
                    break;
                }
            }
            if (same)
            {
                found = new_v;
                break;
            }
        }

        if (found < 0)
        {
            old_to_new[static_cast<std::size_t>(old_v)] = new_nv++;
            append_vertex(old_v);
            continue;
        }

        old_to_new[static_cast<std::size_t>(old_v)] = found;
        const T old_abs_phi = std::abs(
            mesh.vertex_phi[static_cast<std::size_t>(old_v * mesh.n_level_sets + level_set_id)]);
        const T new_abs_phi = std::abs(
            new_phi[static_cast<std::size_t>(found * mesh.n_level_sets + level_set_id)]);
        if (old_abs_phi < new_abs_phi)
            overwrite_vertex(found, old_v);

        new_zero_mask[static_cast<std::size_t>(found)]
            |= mesh.vertex_zero_mask[static_cast<std::size_t>(old_v)];
        new_inside_mask[static_cast<std::size_t>(found)]
            |= mesh.vertex_inside_mask[static_cast<std::size_t>(old_v)];
    }

    if (new_nv == old_nv)
        return;

    for (int32_t& v : mesh.cell_vertices)
        v = old_to_new[static_cast<std::size_t>(v)];

    mesh.vertex_x.swap(new_x);
    mesh.vertex_ref_x.swap(new_ref_x);
    mesh.vertex_parent_dim.swap(new_parent_dim);
    mesh.vertex_parent_id.swap(new_parent_id);
    mesh.vertex_root_edge_id.swap(new_root_edge_id);
    mesh.vertex_phi.swap(new_phi);
    mesh.vertex_zero_mask.swap(new_zero_mask);
    mesh.vertex_inside_mask.swap(new_inside_mask);
}

template <std::floating_point T>
constexpr T root_vertex_merge_tol()
{
    if constexpr (std::is_same_v<T, float>)
        return static_cast<T>(1e-5);
    else
        return static_cast<T>(1e-11);
}

template <std::floating_point T>
int compute_edge_root_linear(LocalMesh<T>& mesh, int edge_id, int level_set_id)
{
    if (edge_id < 0 || edge_id >= mesh.n_edges())
        throw std::invalid_argument("compute_edge_root_linear: invalid edge id");
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("compute_edge_root_linear: invalid level_set_id");

    if (mesh.edge_root_vertex_for(edge_id, level_set_id) >= 0)
        return mesh.edge_root_vertex_for(edge_id, level_set_id);

    int v0_idx = mesh.edge_vertices[edge_id * 2];
    int v1_idx = mesh.edge_vertices[edge_id * 2 + 1];

    int n_ls = mesh.n_level_sets;
    T val0 = mesh.vertex_phi[v0_idx * n_ls + level_set_id];
    T val1 = mesh.vertex_phi[v1_idx * n_ls + level_set_id];

    // Linear interpolation factor s in [0, 1]
    const T denom = (val1 - val0);
    T s = T(0.5);
    if (std::abs(denom) > std::numeric_limits<T>::epsilon())
        s = -val0 / denom;
    s = std::clamp(s, T(0), T(1));

    // Create new vertex
    int new_v_idx = mesh.n_vertices();
    int gdim = mesh.gdim;
    int tdim = mesh.tdim;

    for (int d = 0; d < gdim; ++d) {
        mesh.vertex_x.push_back(mesh.vertex_x[v0_idx * gdim + d] + s * (mesh.vertex_x[v1_idx * gdim + d] - mesh.vertex_x[v0_idx * gdim + d]));
    }
    for (int d = 0; d < tdim; ++d) {
        mesh.vertex_ref_x.push_back(mesh.vertex_ref_x[v0_idx * tdim + d] + s * (mesh.vertex_ref_x[v1_idx * tdim + d] - mesh.vertex_ref_x[v0_idx * tdim + d]));
    }

    // Parent info for the new root
    mesh.vertex_parent_dim.push_back(mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_parent_id.push_back(mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_root_edge_id.push_back(edge_id);

    // Level set values at root are approx 0
    mesh.vertex_phi.resize(mesh.vertex_phi.size() + n_ls, 0.0);
    mesh.vertex_phi[(new_v_idx * n_ls) + level_set_id] = 0.0;
    mesh.vertex_zero_mask.push_back(1ULL << level_set_id);
    mesh.vertex_inside_mask.push_back(0);

    mesh.edge_root_vertex_for(edge_id, level_set_id) = new_v_idx;
    mesh.edge_root_parameter_for(edge_id, level_set_id) = s;
    mesh.edge_root_iterations_for(edge_id, level_set_id) = 0;
    mesh.edge_root_evaluations_for(edge_id, level_set_id) = 0;
    mesh.edge_root_converged_for(edge_id, level_set_id) = 1;
    mesh.edge_root_residual_for(edge_id, level_set_id) = std::abs(val0 + s * (val1 - val0));
    return new_v_idx;
}

template <std::floating_point T, std::integral I>
int compute_edge_root_bernstein(
    LocalMesh<T>&                      mesh,
    const LocalLevelSetFunction<T, I>& level_set,
    int                                edge_id,
    int                                level_set_id,
    cell::edge_root::method            root_method,
    T                                  tol)
{
    if (edge_id < 0 || edge_id >= mesh.n_edges())
        throw std::invalid_argument("compute_edge_root_bernstein: invalid edge id");
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("compute_edge_root_bernstein: invalid level_set_id");
    if (mesh.edge_root_vertex_for(edge_id, level_set_id) >= 0)
        return mesh.edge_root_vertex_for(edge_id, level_set_id);
    if (mesh.vertex_ref_x.size() != static_cast<std::size_t>(mesh.n_vertices() * mesh.tdim))
        throw std::invalid_argument("compute_edge_root_bernstein: missing reference coordinates");

    const int v0_idx = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id)];
    const int v1_idx = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id + 1)];

    const std::span<const T> x0_ref(
        mesh.vertex_ref_x.data() + static_cast<std::size_t>(v0_idx * mesh.tdim),
        static_cast<std::size_t>(mesh.tdim));
    const std::span<const T> x1_ref(
        mesh.vertex_ref_x.data() + static_cast<std::size_t>(v1_idx * mesh.tdim),
        static_cast<std::size_t>(mesh.tdim));

    auto phi = [&](std::span<const T> x) -> T
    {
        return level_set.value(x.data());
    };
    auto grad = [&](std::span<const T> x, std::span<T> g)
    {
        level_set.grad(x.data(), g.data());
    };
    const auto root_info = solve_edge_root_on_segment<T>(
        x0_ref, x1_ref, phi, grad, level_set.has_gradient(), root_method, tol);
    const T t = std::clamp(root_info.t, T(0), T(1));

    const int new_v_idx = mesh.n_vertices();
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;

    for (int d = 0; d < gdim; ++d)
    {
        mesh.vertex_x.push_back(mesh.vertex_x[static_cast<std::size_t>(v0_idx * gdim + d)]
                                + t * (mesh.vertex_x[static_cast<std::size_t>(v1_idx * gdim + d)]
                                       - mesh.vertex_x[static_cast<std::size_t>(v0_idx * gdim + d)]));
    }
    for (int d = 0; d < tdim; ++d)
    {
        mesh.vertex_ref_x.push_back(mesh.vertex_ref_x[static_cast<std::size_t>(v0_idx * tdim + d)]
                                    + t * (mesh.vertex_ref_x[static_cast<std::size_t>(v1_idx * tdim + d)]
                                           - mesh.vertex_ref_x[static_cast<std::size_t>(v0_idx * tdim + d)]));
    }

    mesh.vertex_parent_dim.push_back(mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_parent_id.push_back(mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_root_edge_id.push_back(edge_id);

    mesh.vertex_phi.resize(mesh.vertex_phi.size() + static_cast<std::size_t>(n_ls), T(0));
    mesh.vertex_zero_mask.push_back(0);
    mesh.vertex_inside_mask.push_back(0);

    const T* x_root_ref = mesh.vertex_ref_x.data() + static_cast<std::size_t>(new_v_idx * tdim);
    const T val = level_set.value(x_root_ref);
    mesh.vertex_phi[static_cast<std::size_t>(new_v_idx * n_ls + level_set_id)] = val;
    const uint64_t mask = (uint64_t(1) << level_set_id);
    if (std::abs(val) <= tol)
        mesh.vertex_zero_mask.back() |= mask;
    else if (val < T(0))
        mesh.vertex_inside_mask.back() |= mask;

    mesh.edge_root_vertex_for(edge_id, level_set_id) = new_v_idx;
    mesh.edge_root_parameter_for(edge_id, level_set_id) = t;
    mesh.edge_root_iterations_for(edge_id, level_set_id) = root_info.iterations;
    mesh.edge_root_evaluations_for(edge_id, level_set_id) = root_info.evaluations;
    mesh.edge_root_converged_for(edge_id, level_set_id) = root_info.converged ? 1 : 0;
    mesh.edge_root_residual_for(edge_id, level_set_id) = root_info.residual;
    return new_v_idx;
}

template <std::floating_point T, std::integral I>
int compute_edge_root(LocalMesh<T>& mesh,
                      const LevelSetFunction<T, I>& level_set,
                      int edge_id,
                      int level_set_id,
                      cell::edge_root::method root_method)
{
    if (root_method == cell::edge_root::method::linear || !level_set.has_value())
        return compute_edge_root_linear(mesh, edge_id, level_set_id);

    if (edge_id < 0 || edge_id >= mesh.n_edges())
        throw std::invalid_argument("compute_edge_root: invalid edge id");
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("compute_edge_root: invalid level_set_id");
    if (mesh.edge_root_vertex_for(edge_id, level_set_id) >= 0)
        return mesh.edge_root_vertex_for(edge_id, level_set_id);

    const int v0_idx = mesh.edge_vertices[static_cast<std::size_t>(edge_id * 2)];
    const int v1_idx = mesh.edge_vertices[static_cast<std::size_t>(edge_id * 2 + 1)];
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int n_ls = mesh.n_level_sets;

    const std::span<const T> p0(mesh.vertex_x.data() + static_cast<std::size_t>(v0_idx * gdim),
                                static_cast<std::size_t>(gdim));
    const std::span<const T> p1(mesh.vertex_x.data() + static_cast<std::size_t>(v1_idx * gdim),
                                static_cast<std::size_t>(gdim));
    auto phi = [&](std::span<const T> x) -> T
    {
        return level_set.value(x.data(), static_cast<I>(mesh.parent_cell_id));
    };
    auto grad = [&](std::span<const T> x, std::span<T> g)
    {
        level_set.grad(x.data(), static_cast<I>(mesh.parent_cell_id), g.data());
    };

    const auto root_info = solve_edge_root_on_segment<T>(
        p0, p1, phi, grad, level_set.has_gradient(), root_method, static_cast<T>(1e-12));
    const T t = root_info.t;

    const int new_v_idx = mesh.n_vertices();
    for (int d = 0; d < gdim; ++d)
    {
        mesh.vertex_x.push_back(mesh.vertex_x[static_cast<std::size_t>(v0_idx * gdim + d)]
                                + t * (mesh.vertex_x[static_cast<std::size_t>(v1_idx * gdim + d)]
                                       - mesh.vertex_x[static_cast<std::size_t>(v0_idx * gdim + d)]));
    }
    for (int d = 0; d < tdim; ++d)
    {
        mesh.vertex_ref_x.push_back(mesh.vertex_ref_x[static_cast<std::size_t>(v0_idx * tdim + d)]
                                    + t * (mesh.vertex_ref_x[static_cast<std::size_t>(v1_idx * tdim + d)]
                                           - mesh.vertex_ref_x[static_cast<std::size_t>(v0_idx * tdim + d)]));
    }

    mesh.vertex_parent_dim.push_back(mesh.edge_parent_dim[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_parent_id.push_back(mesh.edge_parent_id[static_cast<std::size_t>(edge_id)]);
    mesh.vertex_root_edge_id.push_back(edge_id);

    mesh.vertex_phi.resize(mesh.vertex_phi.size() + static_cast<std::size_t>(n_ls), T(0));
    mesh.vertex_zero_mask.push_back(0);
    mesh.vertex_inside_mask.push_back(0);
    const T val = phi(std::span<const T>(
        mesh.vertex_x.data() + static_cast<std::size_t>(new_v_idx * gdim),
        static_cast<std::size_t>(gdim)));
    mesh.vertex_phi[static_cast<std::size_t>(new_v_idx * n_ls + level_set_id)] = val;
    const uint64_t mask = (uint64_t(1) << level_set_id);
    const T tol = static_cast<T>(1e-12);
    if (std::abs(val) <= tol)
        mesh.vertex_zero_mask.back() |= mask;
    else if (val < T(0))
        mesh.vertex_inside_mask.back() |= mask;

    mesh.edge_root_vertex_for(edge_id, level_set_id) = new_v_idx;
    mesh.edge_root_parameter_for(edge_id, level_set_id) = t;
    mesh.edge_root_iterations_for(edge_id, level_set_id) = root_info.iterations;
    mesh.edge_root_evaluations_for(edge_id, level_set_id) = root_info.evaluations;
    mesh.edge_root_converged_for(edge_id, level_set_id) = root_info.converged ? 1 : 0;
    mesh.edge_root_residual_for(edge_id, level_set_id) = root_info.residual;
    return new_v_idx;
}

template <std::floating_point T>
void compute_all_roots_linear(LocalMesh<T>& mesh, int level_set_id)
{
    drop_unreferenced_root_vertices(mesh);
    clear_edge_root_cache(mesh);
    int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i) {
        if (mesh.edge_state_for(i, level_set_id) == static_cast<uint8_t>(EdgeState::one_root)) {
            compute_edge_root_linear(mesh, i, level_set_id);
        }
    }
}

template <std::floating_point T, std::integral I>
void compute_all_roots(LocalMesh<T>& mesh,
                       const LevelSetFunction<T, I>& level_set,
                       int level_set_id,
                       cell::edge_root::method root_method)
{
    if (root_method == cell::edge_root::method::linear)
    {
        compute_all_roots_linear(mesh, level_set_id);
        return;
    }

    drop_unreferenced_root_vertices(mesh);
    clear_edge_root_cache(mesh);

    const int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i)
    {
        if (mesh.edge_state_for(i, level_set_id) != static_cast<uint8_t>(EdgeState::one_root))
            continue;
        compute_edge_root<T, I>(mesh, level_set, i, level_set_id, root_method);
    }
}

template <std::floating_point T, std::integral I>
void compute_all_roots_with_backend(LocalMesh<T>& mesh,
                                    const LevelSetFunction<T, I>& level_set,
                                    LocalLevelSetBackend backend,
                                    int level_set_id,
                                    cell::edge_root::method root_method,
                                    T tol)
{
    if (backend != LocalLevelSetBackend::bernstein)
    {
        compute_all_roots<T, I>(mesh, level_set, level_set_id, root_method);
        return;
    }

    if (!level_set.has_nodal_values())
        throw std::invalid_argument("compute_all_roots_with_backend: Bernstein backend requires nodal values");
    const auto local_level_set
        = level_set.mesh != nullptr
              ? make_local_level_set_function_bernstein<T, I>(
                    *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
              : make_local_level_set_function_bernstein<T, I>(
                    mesh.parent_cell_type, mesh.gdim, level_set,
                    static_cast<I>(mesh.parent_cell_id));

    drop_unreferenced_root_vertices(mesh);
    clear_edge_root_cache(mesh);

    const int ne = mesh.n_edges();
    for (int i = 0; i < ne; ++i)
    {
        if (mesh.edge_state_for(i, level_set_id) != static_cast<uint8_t>(EdgeState::one_root))
            continue;
        compute_edge_root_bernstein<T, I>(
            mesh, local_level_set, i, level_set_id, root_method, tol);
    }
}

namespace
{
template <std::floating_point T>
int find_vertex_by_coords(
    const LocalMesh<T>& mesh,
    const std::span<const T> x,
    T tol = root_vertex_merge_tol<T>())
{
    const int nv = mesh.n_vertices();
    const int gdim = mesh.gdim;
    for (int v = 0; v < nv; ++v)
    {
        bool same = true;
        for (int d = 0; d < gdim; ++d)
        {
            const T dv = mesh.vertex_x[static_cast<std::size_t>(v * gdim + d)] - x[static_cast<std::size_t>(d)];
            if (std::abs(dv) > tol)
            {
                same = false;
                break;
            }
        }
        if (same)
            return v;
    }
    return -1;
}

template <std::floating_point T>
int map_cut_vertex_to_local_mesh(
    LocalMesh<T>&               mesh,
    const cell::CutCell<T>&     cut_cell,
    const int                   lv,
    const cell::type            parent_cell_type,
    const std::span<const int32_t> parent_cell_vertices,
    const std::span<const int32_t> parent_cell_edges,
    const int                   level_set_id)
{
    const int gdim = mesh.gdim;
    const int tdim = mesh.tdim;
    const int nls = mesh.n_level_sets;
    const int ncutv = static_cast<int>(cut_cell._vertex_coords.size()) / gdim;

    if (lv < 0 || lv >= ncutv)
        throw std::invalid_argument("map_cut_vertex_to_local_mesh: invalid local cut-vertex id");

    const bool has_tokens = static_cast<int>(cut_cell._vertex_parent_entity.size()) == ncutv;
    const int32_t token = has_tokens ? cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv)] : -1;

    if (token >= 100 && token < 200)
    {
        const int local_parent_v = static_cast<int>(token - 100);
        if (local_parent_v >= 0 && local_parent_v < static_cast<int>(parent_cell_vertices.size()))
            return parent_cell_vertices[static_cast<std::size_t>(local_parent_v)];
    }
    else if (token >= 0 && token < 100)
    {
        const int local_parent_e = static_cast<int>(token);
        const int local_mesh_e = local_edge_for_cut_edge_token(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type, local_parent_e);
        if (local_mesh_e >= 0 && local_mesh_e < mesh.n_edges())
        {
            const auto edge_state = static_cast<EdgeState>(
                mesh.edge_state_for(local_mesh_e, level_set_id));
            int rv = mesh.edge_root_vertex_for(local_mesh_e, level_set_id);
            if (rv < 0 && edge_state == EdgeState::one_root)
            {
                rv = compute_edge_root_linear(mesh, local_mesh_e, level_set_id);
            }
            if (rv >= 0)
                return rv;

            // The legacy cut tables still emit edge-root tokens in zero-touch
            // cases. Reuse exact zero endpoints instead of fabricating a
            // spurious interior root vertex.
            if (edge_state != EdgeState::single_cross)
            {
                const int ev0 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e)];
                const int ev1 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e + 1)];
                const uint64_t mask = (uint64_t(1) << level_set_id);
                const bool zero0 = (mesh.vertex_zero_mask[static_cast<std::size_t>(ev0)] & mask) != 0;
                const bool zero1 = (mesh.vertex_zero_mask[static_cast<std::size_t>(ev1)] & mask) != 0;

                if (zero0 && !zero1)
                    return ev0;
                if (zero1 && !zero0)
                    return ev1;
                if (zero0 && zero1)
                {
                    const std::span<const T> xv(
                        cut_cell._vertex_coords.data() + static_cast<std::size_t>(lv * gdim),
                        static_cast<std::size_t>(gdim));
                    T dist0 = T(0);
                    T dist1 = T(0);
                    for (int d = 0; d < gdim; ++d)
                    {
                        const T dx0 = xv[static_cast<std::size_t>(d)]
                            - mesh.vertex_x[static_cast<std::size_t>(ev0 * gdim + d)];
                        const T dx1 = xv[static_cast<std::size_t>(d)]
                            - mesh.vertex_x[static_cast<std::size_t>(ev1 * gdim + d)];
                        dist0 += dx0 * dx0;
                        dist1 += dx1 * dx1;
                    }
                    return (dist0 <= dist1) ? ev0 : ev1;
                }
            }
        }
    }

    // Fallback (rare): lookup by coordinates and add if not found.
    const std::span<const T> xv(
        cut_cell._vertex_coords.data() + static_cast<std::size_t>(lv * gdim),
        static_cast<std::size_t>(gdim));
    const int existing = find_vertex_by_coords(mesh, xv);
    if (existing >= 0)
        return existing;

    const int new_v = mesh.n_vertices();
    for (int d = 0; d < gdim; ++d)
        mesh.vertex_x.push_back(xv[static_cast<std::size_t>(d)]);

    // Reconstruct reference coordinates for vertices that originate from
    // a parent-cell vertex or cut-edge root token. This avoids corrupting
    // later curved-interface reconstruction with zero reference points.
    if (!mesh.vertex_ref_x.empty())
    {
        std::vector<T> x_ref(static_cast<std::size_t>(tdim), T(0));

        if (token >= 100 && token < 200)
        {
            const int local_parent_v = static_cast<int>(token - 100);
            if (local_parent_v >= 0 && local_parent_v < static_cast<int>(parent_cell_vertices.size()))
            {
                const int gv = parent_cell_vertices[static_cast<std::size_t>(local_parent_v)];
                for (int d = 0; d < tdim; ++d)
                {
                    x_ref[static_cast<std::size_t>(d)]
                        = mesh.vertex_ref_x[static_cast<std::size_t>(gv * tdim + d)];
                }
            }
        }
        else if (token >= 0 && token < 100)
        {
            const int local_parent_e = static_cast<int>(token);
            const int local_mesh_e = local_edge_for_cut_edge_token(
                mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type, local_parent_e);
            if (local_mesh_e >= 0 && local_mesh_e < mesh.n_edges())
            {
                const int ev0 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e)];
                const int ev1 = mesh.edge_vertices[static_cast<std::size_t>(2 * local_mesh_e + 1)];

                T s = T(0.5);
                T den = T(0);
                T num = T(0);
                for (int d = 0; d < gdim; ++d)
                {
                    const T dx = mesh.vertex_x[static_cast<std::size_t>(ev1 * gdim + d)]
                               - mesh.vertex_x[static_cast<std::size_t>(ev0 * gdim + d)];
                    den += dx * dx;
                    num += (xv[static_cast<std::size_t>(d)]
                            - mesh.vertex_x[static_cast<std::size_t>(ev0 * gdim + d)]) * dx;
                }
                if (den > std::numeric_limits<T>::epsilon())
                    s = std::clamp(num / den, T(0), T(1));

                for (int d = 0; d < tdim; ++d)
                {
                    const T r0 = mesh.vertex_ref_x[static_cast<std::size_t>(ev0 * tdim + d)];
                    const T r1 = mesh.vertex_ref_x[static_cast<std::size_t>(ev1 * tdim + d)];
                    x_ref[static_cast<std::size_t>(d)] = r0 + s * (r1 - r0);
                }
            }
        }

        mesh.vertex_ref_x.insert(mesh.vertex_ref_x.end(), x_ref.begin(), x_ref.end());
    }

    int32_t parent_dim = -1;
    int32_t parent_id = -1;
    int32_t root_edge_id = -1;
    if (token >= 100 && token < 200)
    {
        parent_dim = 0;
        const int local_parent_v = static_cast<int>(token - 100);
        if (local_parent_v >= 0 && local_parent_v < static_cast<int>(parent_cell_vertices.size()))
            parent_id = parent_cell_vertices[static_cast<std::size_t>(local_parent_v)];
    }
    else if (token >= 0 && token < 100)
    {
        parent_dim = 1;
        const int local_parent_e = static_cast<int>(token);
        parent_id = local_edge_for_cut_edge_token(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type, local_parent_e);
        root_edge_id = parent_id;
    }
    mesh.vertex_parent_dim.push_back(parent_dim);
    mesh.vertex_parent_id.push_back(parent_id);
    mesh.vertex_root_edge_id.push_back(root_edge_id);

    mesh.vertex_phi.resize(mesh.vertex_phi.size() + static_cast<std::size_t>(nls), T(0));
    mesh.vertex_zero_mask.push_back(0);
    mesh.vertex_inside_mask.push_back(0);
    return new_v;
}

template <std::floating_point T>
void append_cut_fragments(
    LocalMesh<T>&                     mesh,
    const cell::CutCell<T>&           cut_cell,
    const cell::type                  parent_cell_type,
    const std::span<const int32_t>    parent_cell_vertices,
    const std::span<const int32_t>    parent_cell_edges,
    const cell::domain                dom,
    const int                         level_set_id,
    std::vector<int32_t>&             out_cell_vertices,
    std::vector<int32_t>&             out_cell_offsets,
    std::vector<cell::type>&          out_cell_types,
    std::vector<uint8_t>&             out_cell_domain)
{
    const int ncutv = static_cast<int>(cut_cell._vertex_coords.size()) / mesh.gdim;
    if (ncutv == 0 || cut_cell._offset.empty())
        return;

    std::vector<int32_t> cut_to_local(static_cast<std::size_t>(ncutv), -1);
    for (int lv = 0; lv < ncutv; ++lv)
    {
        cut_to_local[static_cast<std::size_t>(lv)] = map_cut_vertex_to_local_mesh(
            mesh, cut_cell, lv, parent_cell_type, parent_cell_vertices, parent_cell_edges, level_set_id);
    }

    const int nsub = cell::num_cells(cut_cell);
    for (int i = 0; i < nsub; ++i)
    {
        const auto verts = cell::cell_vertices(cut_cell, i);
        out_cell_types.push_back(cut_cell._types[static_cast<std::size_t>(i)]);
        out_cell_domain.push_back(static_cast<uint8_t>(dom));
        for (const int lv : verts)
            out_cell_vertices.push_back(cut_to_local[static_cast<std::size_t>(lv)]);
        out_cell_offsets.push_back(static_cast<int32_t>(out_cell_vertices.size()));
    }
}

inline int sign_with_cut_tol(double v, double tol)
{
    if (v > tol)
        return 1;
    if (v < -tol)
        return -1;
    return 0;
}

template <std::floating_point T>
inline int sign_with_cut_tol(T v, T tol)
{
    if (v > tol)
        return 1;
    if (v < -tol)
        return -1;
    return 0;
}

inline bool is_cut_original_vertex_token(cell::type parent_cell_type, int32_t token)
{
    return token >= 100
           && token < 100 + cell::get_num_vertices(parent_cell_type);
}

inline bool is_cut_root_token(cell::type parent_cell_type, int32_t token)
{
    return token >= 0 && token < cell::num_edges(parent_cell_type);
}

enum class CutRejectReason : uint8_t
{
    none = 0,
    empty_fragment = 1,
    new_volume_edge_intersects = 2,
    root_split_edge_extra_root = 3,
    interface_edge_extra_root = 4,
    reference_reconstruction_failed = 5,
    max_depth_uncertified = 6
};

enum class CutSubcellEdgeRole : uint8_t
{
    original_edge = 0,
    interface_edge = 1,
    root_split_edge = 2,
    new_volume_edge = 3
};

inline const char* cut_reject_reason_name(CutRejectReason reason)
{
    switch (reason)
    {
        case CutRejectReason::none:
            return "none";
        case CutRejectReason::empty_fragment:
            return "empty_fragment";
        case CutRejectReason::new_volume_edge_intersects:
            return "new_volume_edge_intersects";
        case CutRejectReason::root_split_edge_extra_root:
            return "root_split_edge_extra_root";
        case CutRejectReason::interface_edge_extra_root:
            return "interface_edge_extra_root";
        case CutRejectReason::reference_reconstruction_failed:
            return "reference_reconstruction_failed";
        case CutRejectReason::max_depth_uncertified:
            return "max_depth_uncertified";
    }
    return "unknown";
}

inline int first_marked_local_cell(std::span<const uint8_t> marked_cells)
{
    for (int c = 0; c < static_cast<int>(marked_cells.size()); ++c)
    {
        if (marked_cells[static_cast<std::size_t>(c)] != 0)
            return c;
    }
    return -1;
}

template <std::floating_point T>
int owning_local_cell_for_edge(const LocalMesh<T>& mesh, int edge_id)
{
    for (int c = 0; c < mesh.n_cells(); ++c)
    {
        const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int ce = ce0; ce < ce1; ++ce)
        {
            if (mesh.cell_edges_flat[static_cast<std::size_t>(ce)] == edge_id)
                return c;
        }
    }
    return -1;
}

template <std::floating_point T>
[[noreturn]] void throw_max_depth_uncertified(
    const LocalMesh<T>& mesh,
    int                 local_cell_id,
    CutRejectReason     reason)
{
    std::ostringstream oss;
    oss << "decompose_local_mesh_with_backend: parent cell " << mesh.parent_cell_id
        << " local cell " << local_cell_id
        << " uncertified at max_refine_depth (" << cut_reject_reason_name(reason) << ")";
    throw std::runtime_error(oss.str());
}

inline void record_first_reject(
    int             cell_id,
    CutRejectReason reason,
    int*            first_invalid_cell,
    CutRejectReason* first_reject_reason)
{
    if (first_invalid_cell != nullptr && *first_invalid_cell < 0)
        *first_invalid_cell = cell_id;
    if (first_reject_reason != nullptr && *first_reject_reason == CutRejectReason::none)
        *first_reject_reason = reason;
}

template <std::floating_point T>
EdgeOrigin classify_cut_edge_origin(
    cell::type                parent_cell_type,
    const cell::CutCell<T>&   cut_cell,
    int                       lv0,
    int                       lv1)
{
    if (static_cast<int>(cut_cell._vertex_parent_entity.size())
        != static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim)
    {
        return EdgeOrigin::lut_new;
    }

    const int32_t token0 = cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv0)];
    const int32_t token1 = cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv1)];

    if (is_cut_original_vertex_token(parent_cell_type, token0)
        && is_cut_original_vertex_token(parent_cell_type, token1))
    {
        const int pv0 = token0 - 100;
        const int pv1 = token1 - 100;
        return background_edge_id_from_corners(parent_cell_type, pv0, pv1) >= 0
                   ? EdgeOrigin::original
                   : EdgeOrigin::lut_new;
    }

    if (is_cut_root_token(parent_cell_type, token0)
        && is_cut_original_vertex_token(parent_cell_type, token1))
    {
        const auto edge = cut_edge_order(parent_cell_type)[static_cast<std::size_t>(token0)];
        const int pv = token1 - 100;
        return (pv == edge[0] || pv == edge[1])
                   ? EdgeOrigin::root_split
                   : EdgeOrigin::lut_new;
    }

    if (is_cut_root_token(parent_cell_type, token1)
        && is_cut_original_vertex_token(parent_cell_type, token0))
    {
        const auto edge = cut_edge_order(parent_cell_type)[static_cast<std::size_t>(token1)];
        const int pv = token0 - 100;
        return (pv == edge[0] || pv == edge[1])
                   ? EdgeOrigin::root_split
                   : EdgeOrigin::lut_new;
    }

    return EdgeOrigin::lut_new;
}

template <std::floating_point T>
CutSubcellEdgeRole classify_cut_subcell_edge_role(
    cell::type              parent_cell_type,
    const cell::CutCell<T>& cut_cell,
    int                     lv0,
    int                     lv1)
{
    const auto origin = classify_cut_edge_origin(parent_cell_type, cut_cell, lv0, lv1);
    if (origin == EdgeOrigin::original)
        return CutSubcellEdgeRole::original_edge;
    if (origin == EdgeOrigin::root_split)
        return CutSubcellEdgeRole::root_split_edge;

    if (static_cast<int>(cut_cell._vertex_parent_entity.size())
        != static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim)
    {
        return CutSubcellEdgeRole::new_volume_edge;
    }

    const int32_t token0 = cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv0)];
    const int32_t token1 = cut_cell._vertex_parent_entity[static_cast<std::size_t>(lv1)];
    if (is_cut_root_token(parent_cell_type, token0)
        && is_cut_root_token(parent_cell_type, token1))
    {
        return CutSubcellEdgeRole::interface_edge;
    }

    return CutSubcellEdgeRole::new_volume_edge;
}

template <std::floating_point T>
void build_subcell_edges(
    const cell::CutCell<T>&   cut_cell,
    cell::type                parent_cell_type,
    int                       cell_id,
    std::vector<int32_t>&     edges,
    std::vector<uint8_t>&     edge_origin)
{
    const auto subcell_type = cut_cell._types[static_cast<std::size_t>(cell_id)];
    const auto subcell_vertices = cell::cell_vertices(cut_cell, cell_id);
    const auto edge_patterns = cell::edges(subcell_type);

    edges.clear();
    edge_origin.clear();
    edges.reserve(edge_patterns.size() * 2);
    edge_origin.reserve(edge_patterns.size());

    for (const auto& ep : edge_patterns)
    {
        const int lv0 = subcell_vertices[static_cast<std::size_t>(ep[0])];
        const int lv1 = subcell_vertices[static_cast<std::size_t>(ep[1])];
        edges.push_back(lv0);
        edges.push_back(lv1);
        edge_origin.push_back(static_cast<uint8_t>(
            classify_cut_edge_origin(parent_cell_type, cut_cell, lv0, lv1)));
    }
}

template <std::floating_point T, std::integral I>
bool cut_vertex_level_set_value(
    const cell::CutCell<T>&    cut_cell,
    std::span<const T>         parent_ls_values,
    const LevelSetFunction<T, I>& level_set,
    int                        parent_cell_id,
    int                        local_vertex_id,
    T&                         value)
{
    if (static_cast<int>(cut_cell._vertex_parent_entity.size())
        == static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim)
    {
        const int32_t token = cut_cell._vertex_parent_entity[static_cast<std::size_t>(local_vertex_id)];
        if (is_cut_root_token(cut_cell._parent_cell_type, token))
        {
            value = T(0);
            return true;
        }
        if (token >= 100
            && token < 100 + static_cast<int32_t>(parent_ls_values.size()))
        {
            value = parent_ls_values[static_cast<std::size_t>(token - 100)];
            return true;
        }
    }

    if (level_set.has_value())
    {
        const T* x = cut_cell._vertex_coords.data()
                     + static_cast<std::size_t>(local_vertex_id * cut_cell._gdim);
        value = level_set.value(x, static_cast<I>(parent_cell_id));
        return true;
    }

    return false;
}

template <std::floating_point T, std::integral I>
bool edge_contains_root(
    const LevelSetFunction<T, I>&  level_set,
    const cell::CutCell<T>&        cut_cell,
    std::span<const T>             parent_ls_values,
    int                            parent_cell_id,
    int                            lv0,
    int                            lv1,
    T                              tol)
{
    T phi0 = T(0);
    T phi1 = T(0);
    if (!cut_vertex_level_set_value(
            cut_cell, parent_ls_values, level_set, parent_cell_id, lv0, phi0)
        || !cut_vertex_level_set_value(
            cut_cell, parent_ls_values, level_set, parent_cell_id, lv1, phi1))
    {
        return true;
    }

    const int s0 = sign_with_cut_tol(phi0, tol);
    const int s1 = sign_with_cut_tol(phi1, tol);
    if (s0 * s1 < 0)
        return true;

    if (!level_set.has_value())
        return false;

    std::vector<T> xmid(static_cast<std::size_t>(cut_cell._gdim), T(0));
    const T* x0 = cut_cell._vertex_coords.data()
                  + static_cast<std::size_t>(lv0 * cut_cell._gdim);
    const T* x1 = cut_cell._vertex_coords.data()
                  + static_cast<std::size_t>(lv1 * cut_cell._gdim);
    for (int d = 0; d < cut_cell._gdim; ++d)
        xmid[static_cast<std::size_t>(d)] = T(0.5) * (x0[d] + x1[d]);

    const T phim = level_set.value(xmid.data(), static_cast<I>(parent_cell_id));
    const int sm = sign_with_cut_tol(phim, tol);
    if (sm == 0)
        return true;
    if (s0 != 0 && s0 * sm < 0)
        return true;
    if (s1 != 0 && sm * s1 < 0)
        return true;
    return false;
}

template <std::floating_point T>
bool cut_vertex_reference_coords(
    const LocalMesh<T>&              mesh,
    std::span<const int32_t>         parent_cell_vertices,
    std::span<const int32_t>         parent_cell_edges,
    cell::type                       parent_cell_type,
    const cell::CutCell<T>&          cut_cell,
    int                              local_vertex_id,
    int                              level_set_id,
    std::span<T>                     x_ref)
{
    if (static_cast<int>(cut_cell._vertex_parent_entity.size())
        != static_cast<int>(cut_cell._vertex_coords.size()) / cut_cell._gdim)
    {
        return false;
    }
    if (static_cast<int>(x_ref.size()) != mesh.tdim)
        return false;

    const int32_t token = cut_cell._vertex_parent_entity[static_cast<std::size_t>(local_vertex_id)];
    if (token >= 100 && token < 200)
    {
        const int local_parent_v = static_cast<int>(token - 100);
        if (local_parent_v < 0 || local_parent_v >= static_cast<int>(parent_cell_vertices.size()))
            return false;
        const int gv = parent_cell_vertices[static_cast<std::size_t>(local_parent_v)];
        for (int d = 0; d < mesh.tdim; ++d)
            x_ref[static_cast<std::size_t>(d)]
                = mesh.vertex_ref_x[static_cast<std::size_t>(gv * mesh.tdim + d)];
        return true;
    }

    if (token >= 0 && token < 100)
    {
        const int local_parent_e = static_cast<int>(token);
        const int local_mesh_e = local_edge_for_cut_edge_token(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type, local_parent_e);
        if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
            return false;
        const int rv = mesh.edge_root_vertex_for(local_mesh_e, level_set_id);
        if (rv < 0 || rv >= mesh.n_vertices())
            return false;
        for (int d = 0; d < mesh.tdim; ++d)
            x_ref[static_cast<std::size_t>(d)]
                = mesh.vertex_ref_x[static_cast<std::size_t>(rv * mesh.tdim + d)];
        return true;
    }

    return false;
}

template <std::floating_point T, std::integral I>
bool certify_cut_edge_bernstein(
    const LocalMesh<T>&              mesh,
    std::span<const int32_t>         parent_cell_vertices,
    std::span<const int32_t>         parent_cell_edges,
    cell::type                       parent_cell_type,
    const ParentPolynomialContext<T>& parent_poly,
    const cell::CutCell<T>&          cut_cell,
    int                              level_set_id,
    int                              lv0,
    int                              lv1,
    T                                tol,
    CutRejectReason&                 reject_reason)
{
    reject_reason = CutRejectReason::none;
    const auto role = classify_cut_subcell_edge_role(
        parent_cell_type, cut_cell, lv0, lv1);
    if (role == CutSubcellEdgeRole::original_edge)
        return true;

    std::vector<T> x0_ref(static_cast<std::size_t>(mesh.tdim), T(0));
    std::vector<T> x1_ref(static_cast<std::size_t>(mesh.tdim), T(0));
    if (!cut_vertex_reference_coords(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type,
            cut_cell, lv0, level_set_id, std::span<T>(x0_ref.data(), x0_ref.size()))
        || !cut_vertex_reference_coords(
            mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type,
            cut_cell, lv1, level_set_id, std::span<T>(x1_ref.data(), x1_ref.size())))
    {
        reject_reason = CutRejectReason::reference_reconstruction_failed;
        return false;
    }

    const T ev0 = evaluate_bernstein_cell<T>(
        parent_poly.bernstein, std::span<const T>(x0_ref.data(), x0_ref.size()));
    const T ev1 = evaluate_bernstein_cell<T>(
        parent_poly.bernstein, std::span<const T>(x1_ref.data(), x1_ref.size()));
    std::vector<T> edge_coeffs;
    restrict_bernstein_to_segment_1d<T>(
        parent_poly.bernstein,
        std::span<const T>(x0_ref.data(), x0_ref.size()),
        std::span<const T>(x1_ref.data(), x1_ref.size()),
        edge_coeffs);
    const auto raw = classify_edge_from_full_cell_bernstein<T>(
        std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), tol);
    const auto interior_samples = sample_interior_bernstein_edge<T>(
        std::span<const T>(edge_coeffs.data(), edge_coeffs.size()));
    const T root_endpoint_tol = std::max(T(32) * tol, T(1e-10));
    const auto topo = classify_edge_zero_topology<T>(
        ev0, ev1,
        std::span<const std::pair<T, T>>(interior_samples.data(), interior_samples.size()),
        raw.state, tol, root_endpoint_tol);
    const int nonzero_endpoint_sign = sign_with_tol(
        sign_with_tol(ev0, tol) == 0 ? ev1 : ev0, tol);
    bool root_split_interior_same_sign = true;
    for (const auto& sample : interior_samples)
    {
        const int si = sign_with_tol(sample.second, tol);
        if (si != 0 && nonzero_endpoint_sign != 0 && si != nonzero_endpoint_sign)
        {
            root_split_interior_same_sign = false;
            break;
        }
    }

    switch (role)
    {
        case CutSubcellEdgeRole::interface_edge:
        {
            const bool ok = topo.endpoint_zero_count == 2
                            && topo.interior_root_count == 0;
            if (!ok)
                reject_reason = CutRejectReason::interface_edge_extra_root;
            return ok;
        }
        case CutSubcellEdgeRole::root_split_edge:
        {
            const bool ok = topo.endpoint_zero_count == 1
                            && topo.interior_root_count == 0
                            && root_split_interior_same_sign;
            if (!ok)
                reject_reason = CutRejectReason::root_split_edge_extra_root;
            return ok;
        }
        case CutSubcellEdgeRole::new_volume_edge:
        {
            const bool ok = (topo.endpoint_zero_count == 1)
                                ? (topo.interior_root_count == 0
                                   && root_split_interior_same_sign)
                                : (topo.state == EdgeState::uncut);
            if (!ok)
                reject_reason = CutRejectReason::new_volume_edge_intersects;
            return ok;
        }
        case CutSubcellEdgeRole::original_edge:
            return true;
    }

    reject_reason = CutRejectReason::max_depth_uncertified;
    return false;
}

template <std::floating_point T, std::integral I>
bool certify_subcell(
    const LevelSetFunction<T, I>&  level_set,
    const cell::CutCell<T>&        cut_cell,
    cell::type                     parent_cell_type,
    std::span<const T>             parent_ls_values,
    int                            parent_cell_id,
    int                            cell_id,
    T                              tol,
    bool                           debug)
{
    std::vector<int32_t> edges;
    std::vector<uint8_t> edge_origin;
    build_subcell_edges(cut_cell, parent_cell_type, cell_id, edges, edge_origin);

    for (std::size_t e = 0; e < edge_origin.size(); ++e)
    {
        const auto origin = static_cast<EdgeOrigin>(edge_origin[e]);
        if (origin == EdgeOrigin::original)
            continue;

        const int lv0 = edges[2 * e];
        const int lv1 = edges[2 * e + 1];
        if (edge_contains_root(
                level_set, cut_cell, parent_ls_values, parent_cell_id, lv0, lv1, tol))
        {
            if (debug)
                std::cout << "invalid edge detected\n";
            return false;
        }
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool certify_subcell_bernstein(
    const LocalMesh<T>&              mesh,
    std::span<const int32_t>         parent_cell_vertices,
    std::span<const int32_t>         parent_cell_edges,
    const ParentPolynomialContext<T>& parent_poly,
    const cell::CutCell<T>&          cut_cell,
    cell::type                       parent_cell_type,
    int                              level_set_id,
    int                              cell_id,
    T                                tol,
    bool                             debug,
    CutRejectReason*                 reject_reason)
{
    std::vector<int32_t> edges;
    std::vector<uint8_t> edge_origin;
    build_subcell_edges(cut_cell, parent_cell_type, cell_id, edges, edge_origin);

    for (std::size_t e = 0; e < edge_origin.size(); ++e)
    {
        const auto origin = static_cast<EdgeOrigin>(edge_origin[e]);
        if (origin == EdgeOrigin::original)
            continue;

        const int lv0 = edges[2 * e];
        const int lv1 = edges[2 * e + 1];
        CutRejectReason edge_reject_reason = CutRejectReason::none;
        if (!certify_cut_edge_bernstein<T, I>(
                mesh, parent_cell_vertices, parent_cell_edges, parent_cell_type,
                parent_poly, cut_cell, level_set_id, lv0, lv1, tol, edge_reject_reason))
        {
            if (reject_reason != nullptr && *reject_reason == CutRejectReason::none)
                *reject_reason = edge_reject_reason;
            if (debug)
                std::cout << "invalid bernstein edge detected: "
                          << cut_reject_reason_name(edge_reject_reason) << "\n";
            return false;
        }
    }

    return true;
}

template <std::floating_point T, std::integral I>
bool certify_cut(
    const LevelSetFunction<T, I>&  level_set,
    const cell::CutCell<T>&        cut_cell,
    cell::type                     parent_cell_type,
    std::span<const T>             parent_ls_values,
    int                            parent_cell_id,
    T                              tol,
    bool                           debug)
{
    const int nsub = cell::num_cells(cut_cell);
    if (nsub == 0)
    {
        if (debug)
            std::cout << "empty cut detected\n";
        return false;
    }
    for (int i = 0; i < nsub; ++i)
    {
        if (!certify_subcell(
                level_set, cut_cell, parent_cell_type, parent_ls_values,
                parent_cell_id, i, tol, debug))
        {
            return false;
        }
    }
    return true;
}

template <std::floating_point T, std::integral I>
bool certify_cut_bernstein(
    const LocalMesh<T>&              mesh,
    std::span<const int32_t>         parent_cell_vertices,
    std::span<const int32_t>         parent_cell_edges,
    const ParentPolynomialContext<T>& parent_poly,
    const cell::CutCell<T>&          cut_cell,
    cell::type                       parent_cell_type,
    int                              level_set_id,
    T                                tol,
    bool                             debug,
    CutRejectReason*                 reject_reason)
{
    const int nsub = cell::num_cells(cut_cell);
    if (nsub == 0)
    {
        if (reject_reason != nullptr)
            *reject_reason = CutRejectReason::empty_fragment;
        if (debug)
            std::cout << "empty bernstein cut detected\n";
        return false;
    }
    for (int i = 0; i < nsub; ++i)
    {
        CutRejectReason subcell_reject_reason = CutRejectReason::none;
        if (!certify_subcell_bernstein<T, I>(
                mesh, parent_cell_vertices, parent_cell_edges, parent_poly,
                cut_cell, parent_cell_type, level_set_id, i, tol, debug,
                &subcell_reject_reason))
        {
            if (reject_reason != nullptr && *reject_reason == CutRejectReason::none)
                *reject_reason = subcell_reject_reason;
            return false;
        }
    }
    return true;
}

template <std::floating_point T>
cell::domain classify_uncut_touch_cell_domain(std::span<const T> ls_values, T tol);

template <std::floating_point T>
int count_cell_one_root_edges(
    const LocalMesh<T>&      mesh,
    std::span<const int32_t> parent_cell_edges,
    int                      level_set_id);

template <std::floating_point T, std::integral I>
bool certify_cut_on_local_mesh(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    int                           level_set_id,
    bool                          triangulate,
    T                             tol,
    std::vector<uint8_t>&         marked_cells,
    int&                          invalid_cells,
    bool                          debug,
    int*                          first_invalid_cell,
    CutRejectReason*              first_reject_reason)
{
    marked_cells.assign(static_cast<std::size_t>(mesh.n_cells()), uint8_t(0));
    invalid_cells = 0;
    if (first_invalid_cell != nullptr)
        *first_invalid_cell = -1;
    if (first_reject_reason != nullptr)
        *first_reject_reason = CutRejectReason::none;
    const bool use_bernstein_cert
        = level_set.has_nodal_values() && level_set.kind == LevelSetFunction<T, I>::Kind::fem_nodal;
    ParentPolynomialContext<T> parent_poly;
    if (use_bernstein_cert)
        parent_poly = make_parent_polynomial_context<T, I>(mesh.parent_cell_type, level_set);

    const int nc = mesh.n_cells();
    for (int c = 0; c < nc; ++c)
    {
        const int c0 = mesh.cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = mesh.cell_offsets[static_cast<std::size_t>(c + 1)];
        const int nv = c1 - c0;
        const cell::type ct = mesh.cell_types[static_cast<std::size_t>(c)];

        const std::span<const int32_t> parent_cell_vertices(
            mesh.cell_vertices.data() + static_cast<std::size_t>(c0),
            static_cast<std::size_t>(nv));

        const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        const std::span<const int32_t> parent_cell_edges(
            mesh.cell_edges_flat.data() + static_cast<std::size_t>(ce0),
            static_cast<std::size_t>(ce1 - ce0));

        std::vector<T> ls_values(static_cast<std::size_t>(nv), T(0));
        std::vector<T> cell_vertex_coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
        for (int j = 0; j < nv; ++j)
        {
            const int gv = parent_cell_vertices[static_cast<std::size_t>(j)];
            ls_values[static_cast<std::size_t>(j)]
                = mesh.vertex_phi[static_cast<std::size_t>(gv * mesh.n_level_sets + level_set_id)];
            for (int d = 0; d < mesh.gdim; ++d)
            {
                cell_vertex_coords[static_cast<std::size_t>(j * mesh.gdim + d)]
                    = mesh.vertex_x[static_cast<std::size_t>(gv * mesh.gdim + d)];
            }
        }

        const cell::domain dom = cell::classify_cell_domain<T>(
            std::span<const T>(ls_values.data(), ls_values.size()));
        if (dom != cell::domain::intersected)
            continue;
        if (count_cell_one_root_edges(mesh, parent_cell_edges, level_set_id) == 0)
            continue;

        const int n_cut_edges = cell::num_edges(ct);
        std::vector<T> edge_root_coords(static_cast<std::size_t>(n_cut_edges * mesh.gdim), T(0));
        std::vector<uint8_t> edge_has_root(static_cast<std::size_t>(n_cut_edges), uint8_t(0));
        bool missing_root = false;
        for (int te = 0; te < n_cut_edges; ++te)
        {
            const int local_mesh_e = local_edge_for_cut_edge_token(
                mesh, parent_cell_vertices, parent_cell_edges, ct, te);
            if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
                continue;

            if (mesh.edge_state_for(local_mesh_e, level_set_id)
                == static_cast<uint8_t>(EdgeState::one_root))
            {
                const int rv = mesh.edge_root_vertex_for(local_mesh_e, level_set_id);
                if (rv < 0 || rv >= mesh.n_vertices())
                {
                    missing_root = true;
                    continue;
                }
                edge_has_root[static_cast<std::size_t>(te)] = 1;
                for (int d = 0; d < mesh.gdim; ++d)
                {
                    edge_root_coords[static_cast<std::size_t>(te * mesh.gdim + d)]
                        = mesh.vertex_x[static_cast<std::size_t>(rv * mesh.gdim + d)];
                }
            }
        }

        if (missing_root)
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            ++invalid_cells;
            record_first_reject(
                c, CutRejectReason::reference_reconstruction_failed,
                first_invalid_cell, first_reject_reason);
            continue;
        }

        cell::CutCell<T> cut_inside;
        cell::CutCell<T> cut_outside;
        cell::cut_from_cached_roots<T>(
            ct,
            std::span<const T>(cell_vertex_coords.data(), cell_vertex_coords.size()),
            mesh.gdim,
            std::span<const T>(ls_values.data(), ls_values.size()),
            "phi<0",
            std::span<const T>(edge_root_coords.data(), edge_root_coords.size()),
            std::span<const uint8_t>(edge_has_root.data(), edge_has_root.size()),
            cut_inside,
            triangulate);
        cell::cut_from_cached_roots<T>(
            ct,
            std::span<const T>(cell_vertex_coords.data(), cell_vertex_coords.size()),
            mesh.gdim,
            std::span<const T>(ls_values.data(), ls_values.size()),
            "phi>0",
            std::span<const T>(edge_root_coords.data(), edge_root_coords.size()),
            std::span<const uint8_t>(edge_has_root.data(), edge_has_root.size()),
            cut_outside,
            triangulate);

        if (cell::num_cells(cut_inside) == 0 || cell::num_cells(cut_outside) == 0)
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            ++invalid_cells;
            record_first_reject(
                c, CutRejectReason::empty_fragment,
                first_invalid_cell, first_reject_reason);
            continue;
        }

        CutRejectReason inside_reason = CutRejectReason::none;
        const bool inside_ok = use_bernstein_cert
            ? certify_cut_bernstein<T, I>(
                mesh, parent_cell_vertices, parent_cell_edges, parent_poly,
                cut_inside, ct, level_set_id, tol, debug, &inside_reason)
            : certify_cut(
                level_set, cut_inside, ct,
                std::span<const T>(ls_values.data(), ls_values.size()),
                mesh.parent_cell_id, tol, debug);
        CutRejectReason outside_reason = CutRejectReason::none;
        const bool outside_ok = use_bernstein_cert
            ? certify_cut_bernstein<T, I>(
                mesh, parent_cell_vertices, parent_cell_edges, parent_poly,
                cut_outside, ct, level_set_id, tol, debug, &outside_reason)
            : certify_cut(
                level_set, cut_outside, ct,
                std::span<const T>(ls_values.data(), ls_values.size()),
                mesh.parent_cell_id, tol, debug);
        if (!inside_ok || !outside_ok)
        {
            marked_cells[static_cast<std::size_t>(c)] = 1;
            ++invalid_cells;
            record_first_reject(
                c,
                inside_ok ? outside_reason : inside_reason,
                first_invalid_cell,
                first_reject_reason);
        }
    }

    return invalid_cells > 0;
}

template <std::floating_point T, std::integral I>
bool certify_final_local_mesh_edges_bernstein(
    const LocalMesh<T>&              mesh,
    const LevelSetFunction<T, I>&    level_set,
    int                              level_set_id,
    T                                tol,
    std::vector<uint8_t>&            marked_cells,
    int*                             first_invalid_cell,
    CutRejectReason*                 first_reject_reason)
{
    marked_cells.assign(static_cast<std::size_t>(mesh.n_cells()), uint8_t(0));
    if (first_invalid_cell != nullptr)
        *first_invalid_cell = -1;
    if (first_reject_reason != nullptr)
        *first_reject_reason = CutRejectReason::none;

    if (!(level_set.has_nodal_values()
          && level_set.kind == LevelSetFunction<T, I>::Kind::fem_nodal))
    {
        return false;
    }

    const auto parent_poly = make_parent_polynomial_context<T, I>(
        mesh.parent_cell_type, level_set);
    const int ne = mesh.n_edges();
    const T root_endpoint_tol = std::max(T(32) * tol, T(1e-10));
    bool any_invalid = false;

    for (int e = 0; e < ne; ++e)
    {
        const int lv0 = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int lv1 = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if (lv0 < 0 || lv1 < 0 || lv0 >= mesh.n_vertices() || lv1 >= mesh.n_vertices())
            continue;

        std::vector<T> x0_ref(static_cast<std::size_t>(mesh.tdim), T(0));
        std::vector<T> x1_ref(static_cast<std::size_t>(mesh.tdim), T(0));
        for (int d = 0; d < mesh.tdim; ++d)
        {
            x0_ref[static_cast<std::size_t>(d)]
                = mesh.vertex_ref_x[static_cast<std::size_t>(lv0 * mesh.tdim + d)];
            x1_ref[static_cast<std::size_t>(d)]
                = mesh.vertex_ref_x[static_cast<std::size_t>(lv1 * mesh.tdim + d)];
        }

        const T ev0 = evaluate_bernstein_cell<T>(
            parent_poly.bernstein, std::span<const T>(x0_ref.data(), x0_ref.size()));
        const T ev1 = evaluate_bernstein_cell<T>(
            parent_poly.bernstein, std::span<const T>(x1_ref.data(), x1_ref.size()));
        std::vector<T> edge_coeffs;
        restrict_bernstein_to_segment_1d<T>(
            parent_poly.bernstein,
            std::span<const T>(x0_ref.data(), x0_ref.size()),
            std::span<const T>(x1_ref.data(), x1_ref.size()),
            edge_coeffs);
        const auto raw = classify_edge_from_full_cell_bernstein<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), tol);
        const auto interior_samples = sample_interior_bernstein_edge<T>(
            std::span<const T>(edge_coeffs.data(), edge_coeffs.size()));
        const auto topo = classify_edge_zero_topology<T>(
            ev0, ev1,
            std::span<const std::pair<T, T>>(interior_samples.data(), interior_samples.size()),
            raw.state, tol, root_endpoint_tol);

        const bool zero0 = sign_with_tol(ev0, tol) == 0;
        const bool zero1 = sign_with_tol(ev1, tol) == 0;
        const int nonzero_endpoint_sign = sign_with_tol(zero0 ? ev1 : ev0, tol);
        bool root_split_interior_same_sign = true;
        for (const auto& sample : interior_samples)
        {
            const int si = sign_with_tol(sample.second, tol);
            if (si != 0 && nonzero_endpoint_sign != 0 && si != nonzero_endpoint_sign)
            {
                root_split_interior_same_sign = false;
                break;
            }
        }
        bool ok = false;
        CutRejectReason reject_reason = CutRejectReason::none;
        if (zero0 && zero1)
        {
            ok = topo.endpoint_zero_count == 2 && topo.interior_root_count == 0;
            if (!ok)
                reject_reason = CutRejectReason::interface_edge_extra_root;
        }
        else if (zero0 || zero1)
        {
            ok = topo.endpoint_zero_count == 1
                 && topo.interior_root_count == 0
                 && root_split_interior_same_sign;
            if (!ok)
                reject_reason = CutRejectReason::root_split_edge_extra_root;
        }
        else
        {
            ok = topo.state == EdgeState::uncut;
            if (!ok)
                reject_reason = CutRejectReason::new_volume_edge_intersects;
        }

        if (ok)
            continue;

        any_invalid = true;
        const int owner_cell = owning_local_cell_for_edge(mesh, e);
        if (owner_cell >= 0)
        {
            marked_cells[static_cast<std::size_t>(owner_cell)] = 1;
            record_first_reject(
                owner_cell, reject_reason, first_invalid_cell, first_reject_reason);
        }
    }

    return any_invalid;
}

template <std::floating_point T>
cell::domain classify_uncut_touch_cell_domain(std::span<const T> ls_values, T tol)
{
    bool has_neg = false;
    bool has_pos = false;
    for (const T v : ls_values)
    {
        if (std::abs(v) <= tol)
            continue;
        if (v < T(0))
            has_neg = true;
        else
            has_pos = true;
    }

    if (has_neg && !has_pos)
        return cell::domain::inside;
    if (has_pos && !has_neg)
        return cell::domain::outside;
    if (!has_neg && !has_pos)
        return cell::domain::inside;
    return cell::domain::intersected;
}

template <std::floating_point T>
int count_cell_one_root_edges(
    const LocalMesh<T>&      mesh,
    std::span<const int32_t> parent_cell_edges,
    int                      level_set_id)
{
    int count = 0;
    for (const int32_t eid : parent_cell_edges)
    {
        if (eid < 0 || eid >= mesh.n_edges())
            continue;
        if (mesh.edge_state_for(eid, level_set_id)
            == static_cast<uint8_t>(EdgeState::one_root))
        {
            ++count;
        }
    }
    return count;
}

template <std::floating_point T>
bool mark_zero_topology_cells(
    const LocalMesh<T>&  mesh,
    int                  level_set_id,
    std::vector<uint8_t>& marked_cells)
{
    (void)mesh;
    (void)level_set_id;
    marked_cells.clear();
    return false;
}

template <std::floating_point T>
bool bernstein_parent_has_interior_sign_contradiction(
    cell::type         parent_cell_type,
    int                degree,
    std::span<const T> nodal_values,
    T                  tol)
{
    if (degree <= 1)
        return false;

    const RefinementTemplate& tpl = iso_p1_template(parent_cell_type, degree);
    if (static_cast<int>(nodal_values.size()) != tpl.n_vertices)
        return false;

    const int n_corner_vertices = cell::get_num_vertices(parent_cell_type);
    int corner_sign = 0;
    for (int i = 0; i < n_corner_vertices; ++i)
    {
        const int si = sign_with_tol(nodal_values[static_cast<std::size_t>(i)], tol);
        if (si == 0)
            continue;
        if (corner_sign == 0)
            corner_sign = si;
        else if (corner_sign != si)
            return false;
    }
    if (corner_sign == 0)
        return false;

    for (int i = n_corner_vertices; i < tpl.n_vertices; ++i)
    {
        if (tpl.vertex_parent_dim[static_cast<std::size_t>(i)] != tpl.tdim)
            continue;
        const int si = sign_with_tol(nodal_values[static_cast<std::size_t>(i)], tol);
        if (si != 0 && si != corner_sign)
            return true;
    }

    return false;
}

template <std::floating_point T, std::integral I>
void seed_local_mesh_from_template_if_needed(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    LocalLevelSetBackend          backend,
    int                           level_set_id,
    T                             tol)
{
    (void)level_set_id;
    if (backend != LocalLevelSetBackend::bernstein)
        return;
    if (!level_set.has_nodal_values())
        return;
    if (mesh.n_cells() != 1)
        return;
    if (mesh.cell_types.size() != 1
        || mesh.cell_types[0] != mesh.parent_cell_type)
    {
        return;
    }

    const int degree = bernstein_level_set_degree(mesh, level_set);
    const bool interior_contradiction
        = bernstein_parent_has_interior_sign_contradiction<T>(
            mesh.parent_cell_type, degree, level_set.nodal_values, tol);
    if (!interior_contradiction)
        return;

    // Keep template seeding only for true parent-cell interior contradictions.
    // Parent-edge multi-cross ambiguity is handled later by the normal shared
    // local-cell marking/refinement path so nearby uncut children are not
    // created up front by a whole-parent iso subdivision.
    refine_local_mesh_from_template(
        mesh, iso_p1_template(mesh.parent_cell_type, degree));
}

// Green split functions are in green_split.h / green_split.cpp.

} // namespace
// ============================================================================
// resolve_multi_root_edges — green-split loop for multi-cross resolution
// ============================================================================

template <std::floating_point T, std::integral I>
bool resolve_multi_root_edges(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    LocalLevelSetBackend          backend,
    int                           level_set_id,
    T                             tol,
    int                           max_iterations)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("resolve_multi_root_edges: invalid level_set_id");
    const bool debug_multiroot = debug_multiroot_enabled();

    // Build Bernstein local level-set function once (for polynomial backend)
    // The segment_restriction closure captures the parent polynomial and
    // can produce 1D Bernstein coefficients for any edge given its ref coords.
    const bool use_bernstein = (backend == LocalLevelSetBackend::bernstein);
    LocalLevelSetFunction<T, I> local_ls;
    if (use_bernstein)
    {
        local_ls = level_set.mesh != nullptr
            ? make_local_level_set_function_bernstein<T, I>(
                  *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
            : make_local_level_set_function_bernstein<T, I>(
                  mesh.parent_cell_type, mesh.gdim, level_set,
                  static_cast<I>(mesh.parent_cell_id));
    }

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // (Re-)classify all edges
        if (use_bernstein)
        {
            classify_edges_from_local_level_set<T, I>(
                mesh, local_ls, level_set_id, tol, false, nullptr);
        }
        else
        {
            classify_edges_with_backend<T, I>(
                mesh, level_set, backend, level_set_id, tol, false);
        }

        // Find the first multi_cross or uncertain edge
        int target_edge = -1;
        const int ne = mesh.n_edges();
        for (int e = 0; e < ne; ++e)
        {
            const auto st = static_cast<EdgeState>(
                mesh.edge_state_for(e, level_set_id));
            if (st == EdgeState::multi_cross || st == EdgeState::uncertain)
            {
                target_edge = e;
                break;
            }
        }

        if (target_edge < 0)
        {
            if (debug_multiroot)
            {
                std::cerr << "Debug: resolve_multi_root_edges parent cell "
                          << mesh.parent_cell_id
                          << " resolved after " << iter
                          << " iteration(s); n_cells=" << mesh.n_cells()
                          << ", n_edges=" << mesh.n_edges() << "\n";
            }
            return true; // all edges are at most single_cross
        }

        if (debug_multiroot)
        {
            const auto st = static_cast<EdgeState>(
                mesh.edge_state_for(target_edge, level_set_id));
            const int va = mesh.edge_vertices[static_cast<std::size_t>(2 * target_edge)];
            const int vb = mesh.edge_vertices[static_cast<std::size_t>(2 * target_edge + 1)];
            std::cerr << "Debug: resolve_multi_root_edges parent cell "
                      << mesh.parent_cell_id
                      << " iter " << iter
                      << " target_edge=" << target_edge
                      << " state=" << static_cast<int>(st)
                      << " verts=(" << va << "," << vb << ")"
                      << " n_cells=" << mesh.n_cells()
                      << " n_edges=" << mesh.n_edges() << "\n";
        }

        // Find separator parameter for this edge
        T t_sep = T(0.5);
        if (use_bernstein
            && mesh.vertex_ref_x.size() == static_cast<std::size_t>(mesh.n_vertices() * mesh.tdim))
        {
            const int v0 = mesh.edge_vertices[static_cast<std::size_t>(2 * target_edge)];
            const int v1 = mesh.edge_vertices[static_cast<std::size_t>(2 * target_edge + 1)];
            const std::span<const T> x0_ref(
                mesh.vertex_ref_x.data() + static_cast<std::size_t>(v0 * mesh.tdim),
                static_cast<std::size_t>(mesh.tdim));
            const std::span<const T> x1_ref(
                mesh.vertex_ref_x.data() + static_cast<std::size_t>(v1 * mesh.tdim),
                static_cast<std::size_t>(mesh.tdim));
            std::vector<T> edge_coeffs;
            local_ls.segment_restriction(x0_ref, x1_ref, edge_coeffs);
            t_sep = find_separator_on_edge_bernstein<T>(
                std::span<const T>(edge_coeffs.data(), edge_coeffs.size()), tol);
        }
        else
        {
            t_sep = find_separator_on_edge_midpoint<T>(
                mesh, target_edge, level_set_id, tol);
        }

        // Attempt green split
        const int n_cells_before = mesh.n_cells();
        const bool green_ok = green_split_one_edge(mesh, target_edge, t_sep);
        if (!green_ok)
        {
            // Fall back to red refinement of incident cells
            const int nc = mesh.n_cells();
            std::vector<uint8_t> marks(static_cast<std::size_t>(nc), 0);
            // Mark cells incident on this edge
            // Re-find incident cells (green_split_one_edge may have failed
            // before modifying the mesh, so the edge is still valid)
            const auto incident = cells_incident_on_edge(mesh, target_edge);
            for (const int c : incident)
                marks[static_cast<std::size_t>(c)] = 1;
            red_refine_marked_cells(
                mesh, std::span<const uint8_t>(marks.data(), marks.size()));
        }

        if (debug_multiroot)
        {
            std::cerr << "Debug: resolve_multi_root_edges parent cell "
                      << mesh.parent_cell_id
                      << " iter " << iter
                      << " separator_t=" << t_sep
                      << " green_ok=" << (green_ok ? 1 : 0)
                      << " n_cells_before=" << n_cells_before
                      << " n_cells_after=" << mesh.n_cells()
                      << "\n";
        }

        // Evaluate level-set values at new vertices (needed for next
        // classification round). The classify functions populate
        // vertex_phi from the parent polynomial, so this is handled
        // implicitly in the next iteration.
    }

    // max_iterations reached — check if any multi-cross edges remain
    const int ne = mesh.n_edges();
    for (int e = 0; e < ne; ++e)
    {
        const auto st = static_cast<EdgeState>(
            mesh.edge_state_for(e, level_set_id));
        if (st == EdgeState::multi_cross || st == EdgeState::uncertain)
        {
            if (debug_multiroot)
            {
                std::cerr << "Debug: resolve_multi_root_edges parent cell "
                          << mesh.parent_cell_id
                          << " hit max_iterations with unresolved edge "
                          << e << " state=" << static_cast<int>(st) << "\n";
            }
            return false;
        }
    }
    return true;
}

template <std::floating_point T, std::integral I>
bool certify_cut_subcells(
    const LevelSetFunction<T, I>& level_set,
    const cell::CutCell<T>&       cut_cell,
    cell::type                    parent_cell_type,
    std::span<const T>            parent_ls_values,
    int                           parent_cell_id,
    T                             tol,
    bool                          debug)
{
    return certify_cut(
        level_set, cut_cell, parent_cell_type,
        parent_ls_values, parent_cell_id, tol, debug);
}

template <std::floating_point T>
void decompose_local_mesh_from_cached_roots(
    LocalMesh<T>& mesh,
    int           level_set_id,
    bool          triangulate,
    bool          fill_missing_linear)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh_linear: invalid level_set_id");
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    if (fill_missing_linear)
        compute_all_roots_linear(mesh, level_set_id);

    const std::vector<int32_t> old_cell_vertices = mesh.cell_vertices;
    const std::vector<int32_t> old_cell_offsets = mesh.cell_offsets;
    const std::vector<cell::type> old_cell_types = mesh.cell_types;
    const std::vector<int32_t> old_cell_edge_offsets = mesh.cell_edge_offsets;
    const std::vector<int32_t> old_cell_edges_flat = mesh.cell_edges_flat;

    std::vector<int32_t> new_cell_vertices;
    std::vector<int32_t> new_cell_offsets(1, 0);
    std::vector<cell::type> new_cell_types;
    std::vector<uint8_t> new_cell_domain;

    const int old_nc = static_cast<int>(old_cell_types.size());
    for (int c = 0; c < old_nc; ++c)
    {
        const int c0 = old_cell_offsets[static_cast<std::size_t>(c)];
        const int c1 = old_cell_offsets[static_cast<std::size_t>(c + 1)];
        const int nv = c1 - c0;
        const cell::type ct = old_cell_types[static_cast<std::size_t>(c)];

        const std::span<const int32_t> parent_cell_vertices(
            old_cell_vertices.data() + static_cast<std::size_t>(c0),
            static_cast<std::size_t>(nv));

        const int ce0 = old_cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = old_cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        const std::span<const int32_t> parent_cell_edges(
            old_cell_edges_flat.data() + static_cast<std::size_t>(ce0),
            static_cast<std::size_t>(ce1 - ce0));

        std::vector<T> ls_values(static_cast<std::size_t>(nv), T(0));
        std::vector<T> cell_vertex_coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));
        for (int j = 0; j < nv; ++j)
        {
            const int gv = parent_cell_vertices[static_cast<std::size_t>(j)];
            ls_values[static_cast<std::size_t>(j)] = mesh.vertex_phi[static_cast<std::size_t>(gv * mesh.n_level_sets + level_set_id)];
            for (int d = 0; d < mesh.gdim; ++d)
            {
                cell_vertex_coords[static_cast<std::size_t>(j * mesh.gdim + d)]
                    = mesh.vertex_x[static_cast<std::size_t>(gv * mesh.gdim + d)];
            }
        }

        const cell::domain dom = cell::classify_cell_domain<T>(
            std::span<const T>(ls_values.data(), ls_values.size()));
        if (dom == cell::domain::inside || dom == cell::domain::outside)
        {
            new_cell_types.push_back(ct);
            new_cell_domain.push_back(static_cast<uint8_t>(dom));
            for (int j = 0; j < nv; ++j)
                new_cell_vertices.push_back(parent_cell_vertices[static_cast<std::size_t>(j)]);
            new_cell_offsets.push_back(static_cast<int32_t>(new_cell_vertices.size()));
            continue;
        }

        const int n_one_root = count_cell_one_root_edges(mesh, parent_cell_edges, level_set_id);
        if (n_one_root == 0)
        {
            // Bernstein edge analysis found no root on any edge of this cell.
            // Trust the certification: classify as inside/outside based on the
            // vertex with the largest absolute level-set value.
            T max_abs = T(0);
            T dominant_val = T(0);
            for (const T v : ls_values)
            {
                if (std::abs(v) > max_abs)
                {
                    max_abs = std::abs(v);
                    dominant_val = v;
                }
            }
            const cell::domain no_root_dom = (dominant_val <= T(0))
                                                 ? cell::domain::inside
                                                 : cell::domain::outside;
            new_cell_types.push_back(ct);
            new_cell_domain.push_back(static_cast<uint8_t>(no_root_dom));
            for (int j = 0; j < nv; ++j)
                new_cell_vertices.push_back(parent_cell_vertices[static_cast<std::size_t>(j)]);
            new_cell_offsets.push_back(static_cast<int32_t>(new_cell_vertices.size()));
            continue;
        }

        const int n_cut_edges = cell::num_edges(ct);
        std::vector<T> edge_root_coords(static_cast<std::size_t>(n_cut_edges * mesh.gdim), T(0));
        std::vector<uint8_t> edge_has_root(static_cast<std::size_t>(n_cut_edges), 0);
        for (int te = 0; te < n_cut_edges; ++te)
        {
            const int local_mesh_e = local_edge_for_cut_edge_token(
                mesh, parent_cell_vertices, parent_cell_edges, ct, te);
            if (local_mesh_e < 0 || local_mesh_e >= mesh.n_edges())
                continue;
            if (mesh.edge_state_for(local_mesh_e, level_set_id)
                != static_cast<uint8_t>(EdgeState::one_root))
            {
                continue;
            }
            int rv = mesh.edge_root_vertex_for(local_mesh_e, level_set_id);
            if (rv < 0 && fill_missing_linear)
                rv = compute_edge_root_linear(mesh, local_mesh_e, level_set_id);
            if (rv < 0 || rv >= mesh.n_vertices())
                continue;
            edge_has_root[static_cast<std::size_t>(te)] = 1;
            for (int d = 0; d < mesh.gdim; ++d)
            {
                edge_root_coords[static_cast<std::size_t>(te * mesh.gdim + d)]
                    = mesh.vertex_x[static_cast<std::size_t>(rv * mesh.gdim + d)];
            }
        }

        cell::CutCell<T> cut_inside;
        cell::CutCell<T> cut_outside;
        cell::cut_from_cached_roots<T>(
            ct,
            std::span<const T>(cell_vertex_coords.data(), cell_vertex_coords.size()),
            mesh.gdim,
            std::span<const T>(ls_values.data(), ls_values.size()),
            "phi<0",
            std::span<const T>(edge_root_coords.data(), edge_root_coords.size()),
            std::span<const uint8_t>(edge_has_root.data(), edge_has_root.size()),
            cut_inside,
            triangulate);
        cell::cut_from_cached_roots<T>(
            ct,
            std::span<const T>(cell_vertex_coords.data(), cell_vertex_coords.size()),
            mesh.gdim,
            std::span<const T>(ls_values.data(), ls_values.size()),
            "phi>0",
            std::span<const T>(edge_root_coords.data(), edge_root_coords.size()),
            std::span<const uint8_t>(edge_has_root.data(), edge_has_root.size()),
            cut_outside,
            triangulate);

        if (cell::num_cells(cut_inside) == 0 || cell::num_cells(cut_outside) == 0)
        {
            std::ostringstream oss;
            oss << "decompose_local_mesh_from_cached_roots: parent cell "
                << mesh.parent_cell_id
                << " local cell " << c
                << " produced an empty cut fragment";
            throw std::runtime_error(oss.str());
        }

        append_cut_fragments(
            mesh, cut_inside, ct, parent_cell_vertices, parent_cell_edges,
            cell::domain::inside, level_set_id,
            new_cell_vertices, new_cell_offsets, new_cell_types, new_cell_domain);
        append_cut_fragments(
            mesh, cut_outside, ct, parent_cell_vertices, parent_cell_edges,
            cell::domain::outside, level_set_id,
            new_cell_vertices, new_cell_offsets, new_cell_types, new_cell_domain);
    }

    mesh.cell_vertices.swap(new_cell_vertices);
    mesh.cell_offsets.swap(new_cell_offsets);
    mesh.cell_types.swap(new_cell_types);
    mesh.cell_domain.swap(new_cell_domain);

    deduplicate_zero_vertices(
        mesh, level_set_id, root_vertex_merge_tol<T>(), root_vertex_merge_tol<T>());
    build_local_edges(mesh);
    build_local_faces(mesh);
    rebuild_parent_entity_maps(mesh);
}

template <std::floating_point T>
void decompose_local_mesh_linear(LocalMesh<T>& mesh, int level_set_id, bool triangulate)
{
    decompose_local_mesh_from_cached_roots(mesh, level_set_id, triangulate, true);
}

template <std::floating_point T, std::integral I>
void decompose_local_mesh(LocalMesh<T>& mesh,
                          const LevelSetFunction<T, I>& level_set,
                          int level_set_id,
                          cell::edge_root::method root_method,
                          bool triangulate)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh: invalid level_set_id");
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    compute_all_roots<T, I>(mesh, level_set, level_set_id, root_method);
    decompose_local_mesh_linear(mesh, level_set_id, triangulate);
}

// ============================================================================
// repair_triangulation_diagonals
// ============================================================================

namespace
{

inline void append_interface_split_log_line(
    int         parent_cell_id,
    int         local_cell_id,
    const char* split_name,
    int         all_zero_child_tets)
{
    const char* path = std::getenv("CUTCELLS_INTERFACE_SPLIT_LOG");
    if (path == nullptr || path[0] == '\0')
        return;

    std::ofstream out(path, std::ios::app);
    if (!out)
        return;

    out << parent_cell_id << " " << local_cell_id << " "
        << split_name << " " << all_zero_child_tets << "\n";
}

/// Build the set of edges (as canonically-ordered vertex pairs) from a mesh.
inline std::set<std::pair<int32_t, int32_t>>
collect_edge_set(const std::vector<int32_t>& edge_vertices, int n_edges)
{
    std::set<std::pair<int32_t, int32_t>> result;
    for (int e = 0; e < n_edges; ++e)
    {
        int32_t a = edge_vertices[static_cast<std::size_t>(2 * e)];
        int32_t b = edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if (a > b) std::swap(a, b);
        result.insert({a, b});
    }
    return result;
}

/// Test whether a post-cut edge is a triangulation diagonal.
///
/// A diagonal is an edge that:
///   1. was NOT present in the pre-cut mesh (i.e. introduced by LUT cutting), AND
///   2. is NOT a root-to-root interface edge (both endpoints are root vertices).
template <std::floating_point T>
bool is_local_mesh_diagonal(
    const LocalMesh<T>& mesh, int edge_id,
    const std::set<std::pair<int32_t, int32_t>>& pre_cut_edges)
{
    int32_t va = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id)];
    int32_t vb = mesh.edge_vertices[static_cast<std::size_t>(2 * edge_id + 1)];
    if (va > vb) std::swap(va, vb);

    // Edge existed before the cut — not a diagonal
    if (pre_cut_edges.count({va, vb}))
        return false;

    // Both endpoints are root vertices — this is an interface edge, not a diagonal
    const int32_t root_a = mesh.vertex_root_edge_id[static_cast<std::size_t>(va)];
    const int32_t root_b = mesh.vertex_root_edge_id[static_cast<std::size_t>(vb)];
    if (root_a >= 0 && root_b >= 0)
        return false;

    return true;
}

/// Check whether a diagonal in the decomposed local mesh crosses the
/// interface using multi-point sampling.  Accepts vertex ids directly
/// so that post-swap checks can be done without rebuilding edge arrays.
template <std::floating_point T>
bool local_mesh_diagonal_crosses_by_verts(
    const LocalMesh<T>& mesh,
    int va, int vb,
    int expected_sign,
    std::function<T(const T*, int)> eval_phi,
    T tol)
{
    const int gdim = mesh.gdim;
    static constexpr std::array<double, 5> samples = {0.2, 0.35, 0.5, 0.65, 0.8};

    for (const double t : samples)
    {
        std::array<T, 3> pt = {T(0), T(0), T(0)};
        for (int d = 0; d < gdim; ++d)
        {
            const T xa = mesh.vertex_x[static_cast<std::size_t>(va * gdim + d)];
            const T xb = mesh.vertex_x[static_cast<std::size_t>(vb * gdim + d)];
            pt[static_cast<std::size_t>(d)] = xa + static_cast<T>(t) * (xb - xa);
        }

        const T phi_val = eval_phi(pt.data(), -1);

        if (expected_sign < 0 && phi_val > tol)
            return true;
        if (expected_sign > 0 && phi_val < -tol)
            return true;
    }
    return false;
}

/// Swap a diagonal between two triangles in the local mesh cell arrays.
/// Returns true if the swap was performed.
template <std::floating_point T>
bool swap_local_mesh_diagonal_2d(
    LocalMesh<T>& mesh,
    int cell_i, int cell_j,
    int shared_a, int shared_b)
{
    if (mesh.cell_types[static_cast<std::size_t>(cell_i)] != cell::type::triangle
        || mesh.cell_types[static_cast<std::size_t>(cell_j)] != cell::type::triangle)
        return false;

    const int off_i = mesh.cell_offsets[static_cast<std::size_t>(cell_i)];
    const int off_j = mesh.cell_offsets[static_cast<std::size_t>(cell_j)];

    // Find the unshared vertex in each triangle
    int unshared_i = -1, unshared_j = -1;
    for (int k = 0; k < 3; ++k)
    {
        int v = mesh.cell_vertices[static_cast<std::size_t>(off_i + k)];
        if (v != shared_a && v != shared_b) { unshared_i = v; break; }
    }
    for (int k = 0; k < 3; ++k)
    {
        int v = mesh.cell_vertices[static_cast<std::size_t>(off_j + k)];
        if (v != shared_a && v != shared_b) { unshared_j = v; break; }
    }
    if (unshared_i < 0 || unshared_j < 0)
        return false;

    // Rewrite: T_i -> {shared_a, unshared_i, unshared_j}
    //          T_j -> {shared_b, unshared_j, unshared_i}
    mesh.cell_vertices[static_cast<std::size_t>(off_i + 0)] = shared_a;
    mesh.cell_vertices[static_cast<std::size_t>(off_i + 1)] = unshared_i;
    mesh.cell_vertices[static_cast<std::size_t>(off_i + 2)] = unshared_j;

    mesh.cell_vertices[static_cast<std::size_t>(off_j + 0)] = shared_b;
    mesh.cell_vertices[static_cast<std::size_t>(off_j + 1)] = unshared_j;
    mesh.cell_vertices[static_cast<std::size_t>(off_j + 2)] = unshared_i;

    return true;
}

/// Find the non-crossing edge of a cell for interface-aware green refinement.
///
/// A non-crossing edge has both endpoints on the same side of the interface
/// (both phi < 0 or both phi > 0).  Splitting such an edge creates a new
/// edge from the vertex on the opposite side to the midpoint that passes
/// through the interface in a quasi-normal direction.
///
/// This is the correct edge to split because it guarantees the interface
/// is subdivided quasi-normally, regardless of cell type or cut case.
///
/// Among multiple non-crossing candidates, picks the longest.
///
/// @param mesh          the local mesh (must have vertex_phi populated)
/// @param cell_id       which cell to inspect
/// @param level_set_id  which level set
/// @param out_v0,out_v1 output: vertex ids of the selected edge
/// @param tol           zero tolerance for sign classification
/// @return local-mesh edge id, or -1 if not found
template <std::floating_point T>
int find_interface_split_edge(
    const LocalMesh<T>& mesh, int cell_id, int level_set_id,
    int& out_v0, int& out_v1, T tol)
{
    const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(cell_id)];
    const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(cell_id + 1)];

    int best_edge = -1;
    T max_len = T(-1);

    for (int idx = ce0; idx < ce1; ++idx)
    {
        const int eid = mesh.cell_edges_flat[static_cast<std::size_t>(idx)];
        const int va  = mesh.edge_vertices[static_cast<std::size_t>(2 * eid)];
        const int vb  = mesh.edge_vertices[static_cast<std::size_t>(2 * eid + 1)];

        const T phi_a = mesh.vertex_phi[static_cast<std::size_t>(
            va * mesh.n_level_sets + level_set_id)];
        const T phi_b = mesh.vertex_phi[static_cast<std::size_t>(
            vb * mesh.n_level_sets + level_set_id)];

        // Skip edges that cross the interface (endpoints on different sides)
        if ((phi_a < -tol && phi_b > tol) || (phi_a > tol && phi_b < -tol))
            continue;

        // Skip edges where an endpoint is on the interface (phi ≈ 0)
        if (std::abs(phi_a) <= tol || std::abs(phi_b) <= tol)
            continue;

        // Both endpoints strictly on the same side: candidate non-crossing edge
        T len_sq = T(0);
        for (int d = 0; d < mesh.gdim; ++d)
        {
            const T diff = mesh.vertex_x[static_cast<std::size_t>(va * mesh.gdim + d)]
                         - mesh.vertex_x[static_cast<std::size_t>(vb * mesh.gdim + d)];
            len_sq += diff * diff;
        }
        if (len_sq > max_len)
        {
            max_len = len_sq;
            best_edge = eid;
            out_v0 = va;
            out_v1 = vb;
        }
    }
    return best_edge;
}

/// Find local-mesh edge id given two vertex ids.  Returns -1 if not found.
template <std::floating_point T>
int find_edge_by_vertices(const LocalMesh<T>& mesh, int v0, int v1)
{
    for (int e = 0; e < mesh.n_edges(); ++e)
    {
        const int a = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int b = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if ((a == v0 && b == v1) || (a == v1 && b == v0))
            return e;
    }
    return -1;
}

/// Map a post-cut cell back to its parent pre-cut cell.
///
/// Uses non-root vertices (original mesh vertices) and root-edge incidence
/// to uniquely identify which pre-cut cell produced the given post-cut cell.
/// Returns -1 if no parent can be found.
template <std::floating_point T>
int find_parent_cell(const LocalMesh<T>& post_cut, int cell_id,
                     const LocalMesh<T>& pre_cut)
{
    const int off = post_cut.cell_offsets[static_cast<std::size_t>(cell_id)];
    const int nv  = post_cut.cell_offsets[static_cast<std::size_t>(cell_id + 1)] - off;

    // Gather non-root vertices and one root edge id for disambiguation
    std::vector<int32_t> orig_verts;
    int any_root_edge = -1;
    for (int j = 0; j < nv; ++j)
    {
        const int v = post_cut.cell_vertices[static_cast<std::size_t>(off + j)];
        if (post_cut.vertex_root_edge_id[static_cast<std::size_t>(v)] < 0)
            orig_verts.push_back(static_cast<int32_t>(v));
        else if (any_root_edge < 0)
            any_root_edge = post_cut.vertex_root_edge_id[static_cast<std::size_t>(v)];
    }
    if (orig_verts.empty())
        return -1;

    for (int pc = 0; pc < pre_cut.n_cells(); ++pc)
    {
        const int poff = pre_cut.cell_offsets[static_cast<std::size_t>(pc)];
        const int pnv  = pre_cut.cell_offsets[static_cast<std::size_t>(pc + 1)] - poff;

        // Check all non-root vertices are in this pre-cut cell
        bool all_found = true;
        for (const int32_t ov : orig_verts)
        {
            bool found = false;
            for (int j = 0; j < pnv; ++j)
            {
                if (pre_cut.cell_vertices[static_cast<std::size_t>(poff + j)] == ov)
                {
                    found = true;
                    break;
                }
            }
            if (!found) { all_found = false; break; }
        }
        if (!all_found)
            continue;

        // If we have a root edge, verify it belongs to this cell for disambiguation
        if (any_root_edge >= 0)
        {
            const int ce0 = pre_cut.cell_edge_offsets[static_cast<std::size_t>(pc)];
            const int ce1 = pre_cut.cell_edge_offsets[static_cast<std::size_t>(pc + 1)];
            bool has_edge = false;
            for (int idx = ce0; idx < ce1; ++idx)
            {
                if (pre_cut.cell_edges_flat[static_cast<std::size_t>(idx)] == any_root_edge)
                {
                    has_edge = true;
                    break;
                }
            }
            if (!has_edge)
                continue;
        }

        return pc;
    }
    return -1;
}

} // anonymous namespace

template <std::floating_point T>
bool repair_triangulation_diagonals(
    LocalMesh<T>&                          mesh,
    const std::set<std::pair<int32_t, int32_t>>& pre_cut_edges,
    std::function<T(const T*, int)>        eval_phi,
    int                                    level_set_id,
    T                                      tol)
{
    // Build edge-to-cell adjacency
    const int ne = mesh.n_edges();
    const int nc = mesh.n_cells();
    std::vector<std::vector<int>> edge_to_cells(static_cast<std::size_t>(ne));

    for (int c = 0; c < nc; ++c)
    {
        const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int idx = ce0; idx < ce1; ++idx)
        {
            const int eid = mesh.cell_edges_flat[static_cast<std::size_t>(idx)];
            edge_to_cells[static_cast<std::size_t>(eid)].push_back(c);
        }
    }

    bool all_resolved = true;

    for (int e = 0; e < ne; ++e)
    {
        const auto& inc_cells = edge_to_cells[static_cast<std::size_t>(e)];
        if (inc_cells.size() != 2)
            continue;

        if (!is_local_mesh_diagonal(mesh, e, pre_cut_edges))
            continue;

        // Determine expected sign from the domain of incident cells
        const int ci = inc_cells[0];
        const int cj = inc_cells[1];
        const auto dom_i = static_cast<cell::domain>(mesh.cell_domain[static_cast<std::size_t>(ci)]);
        const auto dom_j = static_cast<cell::domain>(mesh.cell_domain[static_cast<std::size_t>(cj)]);
        if (dom_i != dom_j)
            continue; // different domains: not a triangulation diagonal

        int expected_sign = 0;
        if (dom_i == cell::domain::inside)
            expected_sign = -1;
        else if (dom_i == cell::domain::outside)
            expected_sign = +1;
        else
            continue;

        const int va = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int vb = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];

        if (!local_mesh_diagonal_crosses_by_verts(mesh, va, vb, expected_sign, eval_phi, tol))
            continue; // diagonal is fine

        // Diagonal crosses interface — try 2D swap
        if (mesh.tdim == 2
            && mesh.cell_types[static_cast<std::size_t>(ci)] == cell::type::triangle
            && mesh.cell_types[static_cast<std::size_t>(cj)] == cell::type::triangle)
        {
            if (!swap_local_mesh_diagonal_2d(mesh, ci, cj, va, vb))
            {
                all_resolved = false;
                continue;
            }

            // After swap: T_i = {shared_a, unshared_i, unshared_j}
            // The new diagonal connects unshared_i and unshared_j
            const int off_i = mesh.cell_offsets[static_cast<std::size_t>(ci)];
            const int new_va = mesh.cell_vertices[static_cast<std::size_t>(off_i + 1)];
            const int new_vb = mesh.cell_vertices[static_cast<std::size_t>(off_i + 2)];

            if (local_mesh_diagonal_crosses_by_verts(
                    mesh, new_va, new_vb, expected_sign, eval_phi, tol))
            {
                // New diagonal also bad — revert the swap
                swap_local_mesh_diagonal_2d(mesh, ci, cj, new_va, new_vb);
                all_resolved = false;
            }
        }
        else
        {
            // 3D or non-triangle: mark as unresolved
            all_resolved = false;
        }
    }

    // Rebuild edge/face arrays to reflect any swaps
    build_local_edges(mesh);
    build_local_faces(mesh);

    return all_resolved;
}

// ============================================================================
// decompose_local_mesh_with_backend
// ============================================================================

template <std::floating_point T, std::integral I>
LUTCertificationResult<T, I> decompose_local_mesh_with_backend(
    LocalMesh<T>&                 mesh,
    const LevelSetFunction<T, I>& level_set,
    LocalLevelSetBackend          backend,
    int                           level_set_id,
    cell::edge_root::method       root_method,
    bool                          triangulate,
    int                           max_refine_depth,
    T                             tol,
    bool                          debug,
    bool                          repair_diagonals,
    int                           max_repair_depth)
{
    if (level_set_id < 0 || level_set_id >= mesh.n_level_sets)
        throw std::invalid_argument("decompose_local_mesh_with_backend: invalid level_set_id");
    if (max_refine_depth < 0)
        throw std::invalid_argument("decompose_local_mesh_with_backend: max_refine_depth must be >= 0");

    LUTCertificationResult<T, I> result;
    const auto count_intersected_subcells = [&](const LocalMesh<T>& local_mesh) -> int
    {
        int n_intersected = 0;
        for (const uint8_t dom : local_mesh.cell_domain)
        {
            if (static_cast<cell::domain>(dom) == cell::domain::intersected)
                ++n_intersected;
        }
        return n_intersected;
    };
    const auto count_all_zero_tets_in_range = [&](const LocalMesh<T>& local_mesh,
                                                  int cell_begin,
                                                  int cell_end) -> int
    {
        int count = 0;
        for (int ci = std::max(0, cell_begin);
             ci < std::min(cell_end, local_mesh.n_cells());
             ++ci)
        {
            if (local_mesh.cell_types[static_cast<std::size_t>(ci)]
                != cell::type::tetrahedron)
                continue;
            const int c0 = local_mesh.cell_offsets[static_cast<std::size_t>(ci)];
            const int c1 = local_mesh.cell_offsets[static_cast<std::size_t>(ci + 1)];
            if (c1 - c0 != 4)
                continue;

            bool all_zero = true;
            for (int j = c0; j < c1; ++j)
            {
                const int32_t vid = local_mesh.cell_vertices[static_cast<std::size_t>(j)];
                const T phi = local_mesh.vertex_phi[static_cast<std::size_t>(
                    vid * local_mesh.n_level_sets + level_set_id)];
                if (std::abs(phi) > tol)
                {
                    all_zero = false;
                    break;
                }
            }
            if (all_zero)
                ++count;
        }
        return count;
    };
    if (static_cast<int>(mesh.vertex_root_edge_id.size()) != mesh.n_vertices())
        mesh.vertex_root_edge_id.assign(static_cast<std::size_t>(mesh.n_vertices()), -1);

    // Step 1: Seed the local mesh from a refinement template if needed
    seed_local_mesh_from_template_if_needed<T, I>(
        mesh, level_set, backend, level_set_id, tol);

    // Step 2: Initial classification and red refinement for zero-topology
    // or interior-sign conflicts (pre-LUT refinement loop)
    for (int depth = 0; depth <= max_refine_depth; ++depth)
    {
        if (backend == LocalLevelSetBackend::bernstein)
        {
            certify_local_mesh_polynomial<T, I>(
                mesh, level_set, level_set_id, 4, tol, true);
        }
        else
        {
            classify_edges_with_backend<T, I>(
                mesh, level_set, backend, level_set_id, tol, false);
        }

        std::vector<uint8_t> zero_topology_marks;
        if (mark_zero_topology_cells(mesh, level_set_id, zero_topology_marks))
        {
            result.refine_iterations = depth + 1;
            if (depth >= max_refine_depth)
            {
                // Hit max depth during zero-topology refinement
                decompose_local_mesh_linear(mesh, level_set_id, triangulate);
                result.certified = false;
                result.hit_max_depth = true;
                return result;
            }

            red_refine_marked_cells(
                mesh, std::span<const uint8_t>(
                    zero_topology_marks.data(), zero_topology_marks.size()));
            continue;
        }

        // No zero-topology issues — proceed to multi-root resolution
        break;
    }

    // Step 3: Resolve multi-root edges via green splits
    const bool all_resolved = resolve_multi_root_edges<T, I>(
        mesh, level_set, backend, level_set_id, tol, max_refine_depth * 4);

    if (!all_resolved && debug)
    {
        std::cerr << "decompose_local_mesh_with_backend: parent cell "
                  << mesh.parent_cell_id
                  << " has unresolved multi-root edges after green splits\n";
    }

    // Step 4: Final edge classification (after green splits may have changed topology)
    if (backend == LocalLevelSetBackend::bernstein)
    {
        classify_edges_with_backend<T, I>(
            mesh, level_set, backend, level_set_id, tol, false);
    }

    // Step 5: Compute roots on all single_cross edges
    compute_all_roots_with_backend<T, I>(
        mesh, level_set, backend, level_set_id, root_method, tol);

    // Step 6: LUT cutting with optional diagonal repair
    if (!repair_diagonals || !triangulate)
    {
        decompose_local_mesh_from_cached_roots(mesh, level_set_id, triangulate, false);
        result.invalid_cells = count_intersected_subcells(mesh);
        result.certified = all_resolved && (result.invalid_cells == 0);
        result.refine_iterations = 0;
        return result;
    }

    // Build a physical-space evaluator for diagonal certification.
    // For FEM/bernstein-backed data this evaluates the local polynomial after
    // pulling the sample point back to the parent reference cell.
    std::function<T(const T*, int)> eval_phi;
    LocalLevelSetFunction<T, I> local_phi;
    if (backend == LocalLevelSetBackend::bernstein)
    {
        local_phi = make_local_level_set_function_bernstein<T, I>(
            mesh.parent_cell_type, mesh.gdim, level_set,
            static_cast<I>(mesh.parent_cell_id));
        eval_phi = [&mesh, &local_phi](const T* x, int /*cell_id*/) -> T
        {
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};
            std::array<T, 3> x_ref = {T(0), T(0), T(0)};
            for (int d = 0; d < mesh.gdim; ++d)
                x_phys[static_cast<std::size_t>(d)] = x[d];
            cell::pull_back_affine<T>(
                mesh.parent_cell_type,
                mesh.parent_cell_coords_p1,
                mesh.gdim,
                std::span<const T>(x_phys.data(), static_cast<std::size_t>(mesh.gdim)),
                std::span<T>(x_ref.data(), static_cast<std::size_t>(mesh.tdim)));
            return local_phi.value(x_ref.data());
        };
    }
    else
    {
        const I parent_cell_id = static_cast<I>(mesh.parent_cell_id);
        eval_phi = [&level_set, parent_cell_id](const T* x, int /*cell_id*/) -> T {
            return level_set.value(x, parent_cell_id);
        };
    }

    // Save pre-decomposition state for potential rollback
    const LocalMesh<T> pre_cut_mesh = mesh;

    // Collect the pre-cut edge set so we can identify new edges (diagonals)
    const auto pre_cut_edges = collect_edge_set(
        mesh.edge_vertices, mesh.n_edges());

    // Trial decomposition with triangulation
    decompose_local_mesh_from_cached_roots(mesh, level_set_id, true, false);

    // Build edge-to-cell adjacency for the post-cut mesh
    const int ne_post = mesh.n_edges();
    std::vector<std::vector<int>> edge_to_cells(
        static_cast<std::size_t>(ne_post));
    for (int c = 0; c < mesh.n_cells(); ++c)
    {
        const int ce0 = mesh.cell_edge_offsets[static_cast<std::size_t>(c)];
        const int ce1 = mesh.cell_edge_offsets[static_cast<std::size_t>(c + 1)];
        for (int idx = ce0; idx < ce1; ++idx)
        {
            const int eid = mesh.cell_edges_flat[static_cast<std::size_t>(idx)];
            edge_to_cells[static_cast<std::size_t>(eid)].push_back(c);
        }
    }

    // Reclassify the post-cut mesh and inspect every new edge. This catches
    // unresolved face-subtriangulation edges in addition to the earlier
    // diagonal-crossing heuristic, so the stronger face-enriched split can be
    // triggered when a child-face edge is still intersected after LUT cutting.
    classify_edges_with_backend<T, I>(
        mesh, level_set, backend, level_set_id, tol, false);

    // Check all new triangulation diagonals for interface crossing.
    // For each crossing diagonal, map back to the pre-cut parent cell.
    std::set<int> bad_parent_cells;
    for (int e = 0; e < ne_post; ++e)
    {
        if (!is_local_mesh_diagonal(mesh, e, pre_cut_edges))
            continue;

        const auto& inc = edge_to_cells[static_cast<std::size_t>(e)];
        const auto edge_state = static_cast<EdgeState>(
            mesh.edge_state_for(e, level_set_id));
        const bool unresolved_new_edge =
            edge_state == EdgeState::single_cross
            || edge_state == EdgeState::multi_cross
            || edge_state == EdgeState::uncertain;

        if (unresolved_new_edge)
        {
            for (const int child_cell : inc)
            {
                const int parent = find_parent_cell(mesh, child_cell, pre_cut_mesh);
                if (parent >= 0)
                    bad_parent_cells.insert(parent);
            }
        }

        if (inc.size() != 2)
            continue;

        const int ci = inc[0];
        const int cj = inc[1];
        const auto dom_i = static_cast<cell::domain>(
            mesh.cell_domain[static_cast<std::size_t>(ci)]);
        const auto dom_j = static_cast<cell::domain>(
            mesh.cell_domain[static_cast<std::size_t>(cj)]);
        if (dom_i != dom_j)
            continue;

        int expected_sign = 0;
        if (dom_i == cell::domain::inside)  expected_sign = -1;
        else if (dom_i == cell::domain::outside) expected_sign = +1;
        else continue;

        const int va = mesh.edge_vertices[static_cast<std::size_t>(2 * e)];
        const int vb = mesh.edge_vertices[static_cast<std::size_t>(2 * e + 1)];
        if (!local_mesh_diagonal_crosses_by_verts(
                mesh, va, vb, expected_sign, eval_phi, tol))
            continue;

        // Diagonal crosses the interface — find the parent pre-cut cell
        const int parent = find_parent_cell(mesh, ci, pre_cut_mesh);
        if (parent >= 0)
            bad_parent_cells.insert(parent);
    }

    if (bad_parent_cells.empty())
    {
        // No crossing diagonals — accept the trial cut
        result.certified = all_resolved;
        result.invalid_cells = 0;
        result.refine_iterations = 0;
        return result;
    }

    // Crossing diagonals detected — try the certified interface-split
    // fallbacks on the affected cells. Start with the cheaper 10/12-child
    // split and only escalate to the median-based 19/24-child face-enriched
    // split if the cheaper fallback cannot certify its children.
    mesh = pre_cut_mesh;

    // Apply interface_split on bad cells in reverse cell-id order so that
    // replace_cell_with_children does not invalidate lower indices.
    std::vector<int> bad_cells_vec(bad_parent_cells.rbegin(),
                                   bad_parent_cells.rend());
    int n_interface_splits = 0;
    const char* face_only_env = std::getenv("CUTCELLS_FORCE_FACE_INTERFACE_SPLIT");
    const bool force_face_split_only = face_only_env != nullptr
        && (std::string_view(face_only_env) == "1"
            || std::string_view(face_only_env) == "true"
            || std::string_view(face_only_env) == "TRUE");
    for (int pc : bad_cells_vec)
    {
        if (mesh.cell_types[static_cast<std::size_t>(pc)]
            != cell::type::tetrahedron)
            continue;

        // Try the cheaper certified split first, then the stronger
        // face-enriched certified split for the matching topology.
        bool ok = false;
        if (!force_face_split_only)
        {
            ok = interface_split_topology1_tet(
                mesh, pc, level_set, level_set_id, tol);
        }
        if (ok)
        {
            const int all_zero_child_tets =
                count_all_zero_tets_in_range(mesh, pc, pc + 10);
            std::cerr << "Debug: repair used interface_split_topology1_tet"
                      << " on parent cell " << mesh.parent_cell_id
                      << ", local cell " << pc
                      << ", all_zero_child_tets=" << all_zero_child_tets
                      << "\n";
            append_interface_split_log_line(
                mesh.parent_cell_id, pc, "interface_split_topology1_tet",
                all_zero_child_tets);
        }
        if (!ok && !force_face_split_only)
        {
            ok = interface_split_topology2_tet(
                mesh, pc, level_set, level_set_id, tol);
            if (ok)
            {
                const int all_zero_child_tets =
                    count_all_zero_tets_in_range(mesh, pc, pc + 12);
                std::cerr << "Debug: repair used interface_split_topology2_tet"
                          << " on parent cell " << mesh.parent_cell_id
                          << ", local cell " << pc
                          << ", all_zero_child_tets=" << all_zero_child_tets
                          << "\n";
                append_interface_split_log_line(
                    mesh.parent_cell_id, pc, "interface_split_topology2_tet",
                    all_zero_child_tets);
            }
        }
        if (!ok)
        {
            ok = interface_split_topology1_tet_faces(
                mesh, pc, level_set, level_set_id, tol);
            if (ok)
            {
                const int all_zero_child_tets =
                    count_all_zero_tets_in_range(mesh, pc, pc + 19);
                std::cerr << "Debug: repair used interface_split_topology1_tet_faces"
                          << " on parent cell " << mesh.parent_cell_id
                          << ", local cell " << pc
                          << ", all_zero_child_tets=" << all_zero_child_tets
                          << "\n";
                append_interface_split_log_line(
                    mesh.parent_cell_id, pc, "interface_split_topology1_tet_faces",
                    all_zero_child_tets);
            }
        }
        if (!ok)
        {
            ok = interface_split_topology2_tet_faces(
                mesh, pc, level_set, level_set_id, tol);
            if (ok)
            {
                const int all_zero_child_tets =
                    count_all_zero_tets_in_range(mesh, pc, pc + 24);
                std::cerr << "Debug: repair used interface_split_topology2_tet_faces"
                          << " on parent cell " << mesh.parent_cell_id
                          << ", local cell " << pc
                          << ", all_zero_child_tets=" << all_zero_child_tets
                          << "\n";
                append_interface_split_log_line(
                    mesh.parent_cell_id, pc, "interface_split_topology2_tet_faces",
                    all_zero_child_tets);
            }
        }
        if (ok)
            ++n_interface_splits;
    }

    if (n_interface_splits > 0)
    {
        // Edge arrays were rebuilt by interface_split. Re-classify edges,
        // rerun the multi-root green-split resolver on any new child edges,
        // then recompute roots for the remaining single-cross edges.
        classify_edges_with_backend<T, I>(
            mesh, level_set, backend, level_set_id, tol, false);

        const bool repair_resolved = resolve_multi_root_edges<T, I>(
            mesh, level_set, backend, level_set_id, tol, max_repair_depth * 4);
        if (!repair_resolved)
        {
            std::cerr << "Warning: decompose_local_mesh_with_backend: parent cell "
                      << mesh.parent_cell_id
                      << " still has unresolved multi-root edges after post-interface-split green refinement\n";
        }

        classify_edges_with_backend<T, I>(
            mesh, level_set, backend, level_set_id, tol, false);
        compute_all_roots_with_backend<T, I>(
            mesh, level_set, backend, level_set_id, root_method, tol);
    }

    // Final LUT cut: interface_split children are automatically classified
    // as inside/outside (all non-zero vertices share the same sign), so
    // the decompose step passes them through without cutting.  Remaining
    // intersected cells get normal LUT cutting.
    decompose_local_mesh_from_cached_roots(mesh, level_set_id, true, false);

    result.invalid_cells = count_intersected_subcells(mesh);
    result.certified = all_resolved && (result.invalid_cells == 0);
    result.refine_iterations = 0;
    if (result.invalid_cells > 0)
    {
        std::cerr << "Warning: decompose_local_mesh_with_backend: parent cell "
                  << mesh.parent_cell_id << " still has "
                  << result.invalid_cells
                  << " intersected subcell(s) after repair attempt\n";
    }
    return result;
}

// ============================================================================
// curve_interface_entities_with_backend
// ============================================================================

template <std::floating_point T, std::integral I>
void curve_interface_entities_with_backend(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id,
    int                                geom_order,
    T                                  tol)
{
    mesh.curved_fallback_count = 0;

    if (mesh.tdim == 3)
    {
        curve_zero_entities_with_backend(
            mesh, level_set, backend, level_set_id, geom_order, tol);
        return;
    }

    if (backend == LocalLevelSetBackend::bernstein)
    {
        if (!level_set.has_nodal_values())
            throw std::invalid_argument(
                "curve_interface_entities_with_backend: Bernstein backend requires nodal values");

        const auto local_phi =
            level_set.mesh != nullptr
                ? make_local_level_set_function_bernstein<T, I>(
                      *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
                : make_local_level_set_function_bernstein<T, I>(
                      mesh.parent_cell_type, mesh.gdim, level_set,
                      static_cast<I>(mesh.parent_cell_id));

        auto eval_phi = [&local_phi](const T* x_ref, int /*tdim*/) -> T
        {
            return local_phi.value_fn(x_ref);
        };

        auto eval_grad = [&local_phi](const T* x_ref, int /*tdim*/, T* grad_out)
        {
            local_phi.grad_fn(x_ref, grad_out);
        };

        curve_interface_entities(mesh, eval_phi, eval_grad,
                                 level_set_id, geom_order, tol);
    }
    else if ((backend == LocalLevelSetBackend::analytical_callbacks
              || backend == LocalLevelSetBackend::nodal_signs)
             && level_set.has_value() && level_set.has_gradient())
    {
        const int tdim = mesh.tdim;
        const int gdim = mesh.gdim;
        const int cell_id = mesh.parent_cell_id;

        auto eval_phi = [&](const T* x_ref, int /*td*/) -> T
        {
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};
            cell::ref_to_phys_affine<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref,
                x_phys.data());
            return level_set.value(x_phys.data(), cell_id);
        };

        auto eval_grad = [&](const T* x_ref, int /*td*/, T* grad_out)
        {
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};
            cell::ref_to_phys_affine<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref,
                x_phys.data());

            std::array<T, 3> grad_phys = {T(0), T(0), T(0)};
            level_set.grad(x_phys.data(), cell_id, grad_phys.data());

            std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));
            cell::compute_jacobian<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref,
                J.data());

            for (int k = 0; k < tdim; ++k)
            {
                T val = T(0);
                for (int d = 0; d < gdim; ++d)
                    val += J[static_cast<std::size_t>(d * tdim + k)] * grad_phys[static_cast<std::size_t>(d)];
                grad_out[k] = val;
            }
        };

        curve_interface_entities(mesh, eval_phi, eval_grad,
                                 level_set_id, geom_order, tol);
    }
    else
    {
        throw std::invalid_argument(
            "curve_interface_entities_with_backend: unsupported backend or "
            "level set missing value/gradient functions");
    }
}

template <std::floating_point T, std::integral I>
void curve_zero_entities_with_backend(
    LocalMesh<T>&                      mesh,
    const LevelSetFunction<T, I>&      level_set,
    LocalLevelSetBackend               backend,
    int                                level_set_id,
    int                                geom_order,
    T                                  tol,
    cell::edge_root::method            line_search_method,
    CurveSearchDirection               search_direction,
    std::span<const T>                 custom_direction_ref)
{
    mesh.curved_fallback_count = 0;

    if (backend == LocalLevelSetBackend::bernstein)
    {
        if (!level_set.has_nodal_values())
            throw std::invalid_argument(
                "curve_zero_entities_with_backend: Bernstein backend requires nodal values");

        const auto local_phi =
            level_set.mesh != nullptr
                ? make_local_level_set_function_bernstein<T, I>(
                      *level_set.mesh, level_set, static_cast<I>(mesh.parent_cell_id))
                : make_local_level_set_function_bernstein<T, I>(
                      mesh.parent_cell_type, mesh.gdim, level_set,
                      static_cast<I>(mesh.parent_cell_id));

        auto eval_phi = [&local_phi](const T* x_ref, int /*tdim*/) -> T
        {
            return local_phi.value_fn(x_ref);
        };

        auto eval_grad = [&local_phi](const T* x_ref, int /*tdim*/, T* grad_out)
        {
            local_phi.grad_fn(x_ref, grad_out);
        };

        curve_zero_entities_impl(mesh, eval_phi, eval_grad,
                                 level_set_id, geom_order, tol,
                                 line_search_method, search_direction,
                                 custom_direction_ref);
    }
    else if ((backend == LocalLevelSetBackend::analytical_callbacks
              || backend == LocalLevelSetBackend::nodal_signs)
             && level_set.has_value() && level_set.has_gradient())
    {
        const int tdim = mesh.tdim;
        const int gdim = mesh.gdim;
        const int cell_id = mesh.parent_cell_id;

        auto eval_phi = [&](const T* x_ref, int /*td*/) -> T
        {
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};
            cell::ref_to_phys_affine<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref,
                x_phys.data());
            return level_set.value(x_phys.data(), cell_id);
        };

        auto eval_grad = [&](const T* x_ref, int /*td*/, T* grad_out)
        {
            std::array<T, 3> x_phys = {T(0), T(0), T(0)};
            cell::ref_to_phys_affine<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref,
                x_phys.data());

            std::array<T, 3> grad_phys = {T(0), T(0), T(0)};
            level_set.grad(x_phys.data(), cell_id, grad_phys.data());

            std::vector<T> J(static_cast<std::size_t>(gdim * tdim), T(0));
            cell::compute_jacobian<T>(
                mesh.parent_cell_type,
                std::span<const T>(mesh.parent_cell_coords_p1.data(),
                                   mesh.parent_cell_coords_p1.size()),
                gdim,
                tdim,
                x_ref,
                J.data());

            for (int k = 0; k < tdim; ++k)
            {
                T val = T(0);
                for (int d = 0; d < gdim; ++d)
                    val += J[static_cast<std::size_t>(d * tdim + k)] * grad_phys[static_cast<std::size_t>(d)];
                grad_out[k] = val;
            }
        };

        curve_zero_entities_impl(mesh, eval_phi, eval_grad,
                                 level_set_id, geom_order, tol,
                                 line_search_method, search_direction,
                                 custom_direction_ref);
    }
    else
    {
        throw std::invalid_argument(
            "curve_zero_entities_with_backend: unsupported backend or "
            "level set missing value/gradient functions");
    }
}

// Explicit instantiations for common types
template void build_local_edges<double>(LocalMesh<double>&);
template void build_local_edges<float>(LocalMesh<float>&);
template void build_local_faces<double>(LocalMesh<double>&);
template void build_local_faces<float>(LocalMesh<float>&);
template void rebuild_parent_entity_maps<double>(LocalMesh<double>&);
template void rebuild_parent_entity_maps<float>(LocalMesh<float>&);
template void build_zero_entities<double>(LocalMesh<double>&, int);
template void build_zero_entities<float>(LocalMesh<float>&, int);
template void assign_zero_entity_ownership<double>(std::span<LocalMesh<double>*>, uint64_t);
template void assign_zero_entity_ownership<float>(std::span<LocalMesh<float>*>, uint64_t);
template void build_interface_entities<double>(LocalMesh<double>&, int);
template void build_interface_entities<float>(LocalMesh<float>&, int);
template void build_zero_chains<double>(LocalMesh<double>&, uint64_t);
template void build_zero_chains<float>(LocalMesh<float>&, uint64_t);
template void build_zero_patches<double>(LocalMesh<double>&, uint64_t);
template void build_zero_patches<float>(LocalMesh<float>&, uint64_t);
template void init_local_mesh_from_template<double>(LocalMesh<double>&, const RefinementTemplate&, std::span<const double>, cell::type, int, int);
template void init_local_mesh_from_template<float>(LocalMesh<float>&, const RefinementTemplate&, std::span<const float>, cell::type, int, int);
template void init_local_mesh_from_cell<double>(LocalMesh<double>&, std::span<const double>, cell::type, int, int);
template void init_local_mesh_from_cell<float>(LocalMesh<float>&, std::span<const float>, cell::type, int, int);
template void refine_local_mesh_from_template<double>(LocalMesh<double>&, const RefinementTemplate&);
template void refine_local_mesh_from_template<float>(LocalMesh<float>&, const RefinementTemplate&);
template void red_refine_marked_cells<double>(LocalMesh<double>&, std::span<const uint8_t>);
template void red_refine_marked_cells<float>(LocalMesh<float>&, std::span<const uint8_t>);
template void evaluate_levelsets_on_vertices<double>(LocalMesh<double>&, const std::vector<LevelSetFunction<double>>&, int, double);
template void evaluate_levelsets_on_vertices<float>(LocalMesh<float>&, const std::vector<LevelSetFunction<float>>&, int, float);
template void classify_local_edges<double>(LocalMesh<double>&, int);
template void classify_local_edges<float>(LocalMesh<float>&, int);
template int compute_edge_root_linear<double>(LocalMesh<double>&, int, int);
template int compute_edge_root_linear<float>(LocalMesh<float>&, int, int);
template int compute_edge_root_bernstein<double, int>(LocalMesh<double>&, const LocalLevelSetFunction<double, int>&, int, int, cell::edge_root::method, double);
template int compute_edge_root_bernstein<float, int>(LocalMesh<float>&, const LocalLevelSetFunction<float, int>&, int, int, cell::edge_root::method, float);
template int compute_edge_root<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, int, cell::edge_root::method);
template int compute_edge_root<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, int, cell::edge_root::method);
template void compute_all_roots_with_backend<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, cell::edge_root::method, double);
template void compute_all_roots_with_backend<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, cell::edge_root::method, float);
template bool certify_cut_subcells<double, int>(const LevelSetFunction<double, int>&, const cell::CutCell<double>&, cell::type, std::span<const double>, int, double, bool);
template bool certify_cut_subcells<float, int>(const LevelSetFunction<float, int>&, const cell::CutCell<float>&, cell::type, std::span<const float>, int, float, bool);
template void compute_all_roots_linear<double>(LocalMesh<double>&, int);
template void compute_all_roots_linear<float>(LocalMesh<float>&, int);
template void compute_all_roots<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, cell::edge_root::method);
template void compute_all_roots<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, cell::edge_root::method);
template void decompose_local_mesh_linear<double>(LocalMesh<double>&, int, bool);
template void decompose_local_mesh_linear<float>(LocalMesh<float>&, int, bool);
template void decompose_local_mesh<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, int, cell::edge_root::method, bool);
template void decompose_local_mesh<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, int, cell::edge_root::method, bool);
template LUTCertificationResult<double, int> decompose_local_mesh_with_backend<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, cell::edge_root::method, bool, int, double, bool, bool, int);
template LUTCertificationResult<float, int> decompose_local_mesh_with_backend<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, cell::edge_root::method, bool, int, float, bool, bool, int);
template bool repair_triangulation_diagonals<double>(LocalMesh<double>&, const std::set<std::pair<int32_t, int32_t>>&, std::function<double(const double*, int)>, int, double);
template bool repair_triangulation_diagonals<float>(LocalMesh<float>&, const std::set<std::pair<int32_t, int32_t>>&, std::function<float(const float*, int)>, int, float);
template bool resolve_multi_root_edges<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, double, int);
template bool resolve_multi_root_edges<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, float, int);
template void curve_interface_entities_with_backend<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, int, double);
template void curve_interface_entities_with_backend<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, int, float);
template void curve_zero_entities_with_backend<double, int>(LocalMesh<double>&, const LevelSetFunction<double, int>&, LocalLevelSetBackend, int, int, double, cell::edge_root::method, CurveSearchDirection, std::span<const double>);
template void curve_zero_entities_with_backend<float, int>(LocalMesh<float>&, const LevelSetFunction<float, int>&, LocalLevelSetBackend, int, int, float, cell::edge_root::method, CurveSearchDirection, std::span<const float>);

} // namespace cutcells
