// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT
#pragma once

#include "cell_topology.h"
#include "reference_cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace cutcells::geom
{

/// Local parent entity in the reference parent cell.
/// dim = 0, 1, 2, tdim denotes vertex, edge, face, cell interior.
/// id is local to that entity dimension, except for cell interiors where it is
/// conventionally -1.
struct ParentEntity
{
    int dim = -1;
    int id = -1;

    bool valid() const noexcept { return dim >= 0; }
};

enum class Degeneracy
{
    none = 0,
    zero_input,
    zero_frame,
    zero_projection,
    invalid_parent_entity
};

template <std::floating_point T>
struct VectorQuantity
{
    std::vector<T> value;
    T norm = T(0);
    Degeneracy degeneracy = Degeneracy::none;

    bool degenerate() const noexcept
    {
        return degeneracy != Degeneracy::none;
    }
};

template <std::floating_point T>
struct DisplacementComponents
{
    std::vector<T> normal;
    std::vector<T> tangential;
    T signed_normal_magnitude = T(0);
    bool degenerate_normal = false;
};

template <std::floating_point T>
struct Alignment
{
    T cosine = T(0);
    T angle = T(0);
    bool degenerate = false;
};

template <std::floating_point T>
struct LineInterval
{
    bool valid = false;
    T t0 = T(0);
    T t1 = T(0);
};

template <std::floating_point T>
inline T default_tolerance()
{
    return T(128) * std::numeric_limits<T>::epsilon();
}

template <std::floating_point T>
inline void require_same_size(std::span<const T> a, std::span<const T> b,
                              const char* name)
{
    if (a.size() != b.size())
        throw std::invalid_argument(std::string(name) + ": dimension mismatch");
}

template <std::floating_point T>
inline void require_size(std::span<const T> a, std::size_t n, const char* name)
{
    if (a.size() != n)
        throw std::invalid_argument(std::string(name) + ": unexpected dimension");
}

template <std::floating_point T>
inline T dot(std::span<const T> a, std::span<const T> b)
{
    require_same_size<T>(a, b, "dot");
    T out = T(0);
    for (std::size_t i = 0; i < a.size(); ++i)
        out += a[i] * b[i];
    return out;
}

template <std::floating_point T>
inline T squared_norm(std::span<const T> a)
{
    return dot<T>(a, a);
}

template <std::floating_point T>
inline T norm(std::span<const T> a)
{
    return std::sqrt(squared_norm<T>(a));
}

template <std::floating_point T>
inline std::vector<T> subtract(std::span<const T> a, std::span<const T> b)
{
    require_same_size<T>(a, b, "subtract");
    std::vector<T> out(a.size(), T(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        out[i] = a[i] - b[i];
    return out;
}

template <std::floating_point T>
inline std::array<T, 3> cross(std::span<const T> a, std::span<const T> b)
{
    require_size<T>(a, 3, "cross");
    require_size<T>(b, 3, "cross");
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

template <std::floating_point T>
inline VectorQuantity<T> make_vector_quantity(std::vector<T> value,
                                               T tol,
                                               Degeneracy zero_reason)
{
    const T n = norm<T>(std::span<const T>(value.data(), value.size()));
    return {std::move(value), n, n <= tol ? zero_reason : Degeneracy::none};
}

template <std::floating_point T>
inline std::vector<T> reference_vertex(cell::type parent_cell_type,
                                       int vertex_id)
{
    const int tdim = cell::get_tdim(parent_cell_type);
    const int nverts = cell::get_num_vertices(parent_cell_type);
    if (vertex_id < 0 || vertex_id >= nverts)
        throw std::invalid_argument("reference_vertex: vertex id out of bounds");
    if (parent_cell_type == cell::type::point)
        return {};

    const auto vertices = cell::reference_vertices<T>(parent_cell_type);
    std::vector<T> out(static_cast<std::size_t>(tdim), T(0));
    for (int d = 0; d < tdim; ++d)
        out[static_cast<std::size_t>(d)] =
            vertices[static_cast<std::size_t>(vertex_id * tdim + d)];
    return out;
}

template <std::floating_point T>
inline VectorQuantity<T> segment_tangent(std::span<const T> a,
                                         std::span<const T> b,
                                         bool unit = true,
                                         T tol = default_tolerance<T>())
{
    auto out = subtract<T>(b, a);
    const T n = norm<T>(std::span<const T>(out.data(), out.size()));
    if (n <= tol)
        return {std::move(out), n, Degeneracy::zero_frame};

    if (unit)
    {
        for (T& x : out)
            x /= n;
        return {std::move(out), T(1), Degeneracy::none};
    }

    return {std::move(out), n, Degeneracy::none};
}

template <std::floating_point T>
inline VectorQuantity<T> parent_edge_tangent(cell::type parent_cell_type,
                                             int parent_edge_id,
                                             bool unit = true,
                                             T tol = default_tolerance<T>())
{
    const auto edges = cell::edges(parent_cell_type);
    if (parent_edge_id < 0 || parent_edge_id >= static_cast<int>(edges.size()))
        return {{}, T(0), Degeneracy::invalid_parent_entity};

    const auto edge = edges[static_cast<std::size_t>(parent_edge_id)];
    const auto a = reference_vertex<T>(parent_cell_type, edge[0]);
    const auto b = reference_vertex<T>(parent_cell_type, edge[1]);
    return segment_tangent<T>(
        std::span<const T>(a.data(), a.size()),
        std::span<const T>(b.data(), b.size()),
        unit,
        tol);
}

template <std::floating_point T>
inline VectorQuantity<T> face_normal(std::span<const T> a,
                                     std::span<const T> b,
                                     std::span<const T> c,
                                     bool unit = true,
                                     T tol = default_tolerance<T>())
{
    require_size<T>(a, 3, "face_normal");
    require_size<T>(b, 3, "face_normal");
    require_size<T>(c, 3, "face_normal");

    const auto ab = subtract<T>(b, a);
    const auto ac = subtract<T>(c, a);
    const auto n_array = cross<T>(
        std::span<const T>(ab.data(), ab.size()),
        std::span<const T>(ac.data(), ac.size()));

    std::vector<T> out(n_array.begin(), n_array.end());
    const T n = norm<T>(std::span<const T>(out.data(), out.size()));
    if (n <= tol)
        return {std::move(out), n, Degeneracy::zero_frame};

    if (unit)
    {
        for (T& x : out)
            x /= n;
        return {std::move(out), T(1), Degeneracy::none};
    }

    return {std::move(out), n, Degeneracy::none};
}

template <std::floating_point T>
inline VectorQuantity<T> parent_face_normal(cell::type parent_cell_type,
                                            int parent_face_id,
                                            bool unit = true,
                                            T tol = default_tolerance<T>())
{
    if (cell::get_tdim(parent_cell_type) != 3)
        return {{}, T(0), Degeneracy::invalid_parent_entity};
    if (parent_face_id < 0 || parent_face_id >= cell::num_faces(parent_cell_type))
        return {{}, T(0), Degeneracy::invalid_parent_entity};

    const auto face_vertices = cell::face_vertices(parent_cell_type, parent_face_id);
    if (face_vertices.size() < 3)
        return {{}, T(0), Degeneracy::invalid_parent_entity};

    const auto a = reference_vertex<T>(parent_cell_type, face_vertices[0]);
    const auto b = reference_vertex<T>(parent_cell_type, face_vertices[1]);
    const auto c = reference_vertex<T>(parent_cell_type, face_vertices[2]);
    return face_normal<T>(
        std::span<const T>(a.data(), a.size()),
        std::span<const T>(b.data(), b.size()),
        std::span<const T>(c.data(), c.size()),
        unit,
        tol);
}

template <std::floating_point T>
inline VectorQuantity<T> segment_normal_2d(std::span<const T> a,
                                           std::span<const T> b,
                                           bool unit = true,
                                           T tol = default_tolerance<T>())
{
    require_size<T>(a, 2, "segment_normal_2d");
    require_size<T>(b, 2, "segment_normal_2d");

    std::vector<T> out = {-(b[1] - a[1]), b[0] - a[0]};
    const T n = norm<T>(std::span<const T>(out.data(), out.size()));
    if (n <= tol)
        return {std::move(out), n, Degeneracy::zero_frame};

    if (unit)
    {
        out[0] /= n;
        out[1] /= n;
        return {std::move(out), T(1), Degeneracy::none};
    }

    return {std::move(out), n, Degeneracy::none};
}

template <std::floating_point T>
inline VectorQuantity<T> segment_normal(std::span<const T> a,
                                        std::span<const T> b,
                                        bool unit = true,
                                        T tol = default_tolerance<T>())
{
    return segment_normal_2d<T>(a, b, unit, tol);
}

/// Unit normal to a segment, constrained to a 3D face tangent plane.
/// The orientation is face_normal x segment_tangent.
template <std::floating_point T>
inline VectorQuantity<T> in_face_segment_normal(std::span<const T> a,
                                                std::span<const T> b,
                                                std::span<const T> face_normal_vector,
                                                bool unit = true,
                                                T tol = default_tolerance<T>())
{
    require_size<T>(a, 3, "in_face_segment_normal");
    require_size<T>(b, 3, "in_face_segment_normal");
    require_size<T>(face_normal_vector, 3, "in_face_segment_normal");

    auto tangent = segment_tangent<T>(a, b, false, tol);
    if (tangent.degenerate())
        return {std::move(tangent.value), tangent.norm, tangent.degeneracy};

    const auto n_array = cross<T>(
        face_normal_vector,
        std::span<const T>(tangent.value.data(), tangent.value.size()));
    std::vector<T> out(n_array.begin(), n_array.end());
    const T n = norm<T>(std::span<const T>(out.data(), out.size()));
    if (n <= tol)
        return {std::move(out), n, Degeneracy::zero_frame};

    if (unit)
    {
        for (T& x : out)
            x /= n;
        return {std::move(out), T(1), Degeneracy::none};
    }

    return {std::move(out), n, Degeneracy::none};
}

template <std::floating_point T>
inline VectorQuantity<T> project_onto_line(std::span<const T> vector,
                                           std::span<const T> line_direction,
                                           T tol = default_tolerance<T>())
{
    require_same_size<T>(vector, line_direction, "project_onto_line");
    const T dd = squared_norm<T>(line_direction);
    std::vector<T> out(vector.size(), T(0));
    if (dd <= tol * tol)
        return {std::move(out), T(0), Degeneracy::zero_frame};

    const T scale = dot<T>(vector, line_direction) / dd;
    for (std::size_t i = 0; i < vector.size(); ++i)
        out[i] = scale * line_direction[i];
    return make_vector_quantity<T>(std::move(out), tol, Degeneracy::zero_projection);
}

template <std::floating_point T>
inline VectorQuantity<T> project_into_plane(std::span<const T> vector,
                                            std::span<const T> plane_normal,
                                            T tol = default_tolerance<T>())
{
    require_same_size<T>(vector, plane_normal, "project_into_plane");
    const T nn = squared_norm<T>(plane_normal);
    std::vector<T> out(vector.begin(), vector.end());
    if (nn <= tol * tol)
    {
        std::fill(out.begin(), out.end(), T(0));
        return {std::move(out), T(0), Degeneracy::zero_frame};
    }

    const T scale = dot<T>(vector, plane_normal) / nn;
    for (std::size_t i = 0; i < vector.size(); ++i)
        out[i] -= scale * plane_normal[i];
    return make_vector_quantity<T>(std::move(out), tol, Degeneracy::zero_projection);
}

template <std::floating_point T>
inline VectorQuantity<T> project_into_parent_face_tangent(cell::type parent_cell_type,
                                                          int parent_face_id,
                                                          std::span<const T> vector,
                                                          T tol = default_tolerance<T>())
{
    const auto normal = parent_face_normal<T>(parent_cell_type, parent_face_id, false, tol);
    if (normal.degenerate())
        return {{}, T(0), normal.degeneracy};
    return project_into_plane<T>(
        vector,
        std::span<const T>(normal.value.data(), normal.value.size()),
        tol);
}

template <std::floating_point T>
inline VectorQuantity<T> project_onto_parent_edge(cell::type parent_cell_type,
                                                  int parent_edge_id,
                                                  std::span<const T> vector,
                                                  T tol = default_tolerance<T>())
{
    const auto tangent = parent_edge_tangent<T>(parent_cell_type, parent_edge_id, false, tol);
    if (tangent.degenerate())
        return {{}, T(0), tangent.degeneracy};
    return project_onto_line<T>(
        vector,
        std::span<const T>(tangent.value.data(), tangent.value.size()),
        tol);
}

/// Project a raw direction into the tangent frame admitted by the smallest
/// parent entity that hosts the point. The caller supplies that host entity;
/// use smallest_parent_entity_containing_point when it must be inferred from
/// reference coordinates.
template <std::floating_point T>
inline VectorQuantity<T> admissible_direction_in_parent_frame(
    cell::type parent_cell_type,
    ParentEntity host,
    std::span<const T> raw_direction,
    T tol = default_tolerance<T>())
{
    const int tdim = cell::get_tdim(parent_cell_type);
    if (static_cast<int>(raw_direction.size()) != tdim)
        throw std::invalid_argument(
            "admissible_direction_in_parent_frame: direction dimension mismatch");

    if (!host.valid() || host.dim > tdim)
        return {{}, T(0), Degeneracy::invalid_parent_entity};

    if (host.dim == tdim)
    {
        std::vector<T> out(raw_direction.begin(), raw_direction.end());
        return make_vector_quantity<T>(std::move(out), tol, Degeneracy::zero_input);
    }

    if (host.dim == 2 && tdim == 3)
        return project_into_parent_face_tangent<T>(
            parent_cell_type, host.id, raw_direction, tol);

    if (host.dim == 1)
        return project_onto_parent_edge<T>(
            parent_cell_type, host.id, raw_direction, tol);

    if (host.dim == 0)
        return {std::vector<T>(raw_direction.size(), T(0)),
                T(0),
                Degeneracy::zero_projection};

    return {{}, T(0), Degeneracy::invalid_parent_entity};
}

template <std::floating_point T>
inline VectorQuantity<T> restricted_level_set_gradient_in_parent_frame(
    cell::type parent_cell_type,
    ParentEntity host,
    std::span<const T> raw_gradient,
    T tol = default_tolerance<T>())
{
    return admissible_direction_in_parent_frame<T>(
        parent_cell_type, host, raw_gradient, tol);
}

template <std::floating_point T>
inline DisplacementComponents<T> decompose_displacement(
    std::span<const T> displacement,
    std::span<const T> normal_vector,
    T tol = default_tolerance<T>())
{
    require_same_size<T>(displacement, normal_vector, "decompose_displacement");
    DisplacementComponents<T> out;
    out.normal.assign(displacement.size(), T(0));
    out.tangential.assign(displacement.begin(), displacement.end());

    const T nn = squared_norm<T>(normal_vector);
    if (nn <= tol * tol)
    {
        out.degenerate_normal = true;
        return out;
    }

    const T scale = dot<T>(displacement, normal_vector) / nn;
    const T n = std::sqrt(nn);
    out.signed_normal_magnitude = dot<T>(displacement, normal_vector) / n;
    for (std::size_t i = 0; i < displacement.size(); ++i)
    {
        out.normal[i] = scale * normal_vector[i];
        out.tangential[i] = displacement[i] - out.normal[i];
    }
    return out;
}

template <std::floating_point T>
inline Alignment<T> alignment(std::span<const T> a,
                              std::span<const T> b,
                              T tol = default_tolerance<T>())
{
    require_same_size<T>(a, b, "alignment");
    const T na = norm<T>(a);
    const T nb = norm<T>(b);
    if (na <= tol || nb <= tol)
        return {T(0), T(0), true};

    const T c = std::clamp(dot<T>(a, b) / (na * nb), T(-1), T(1));
    return {c, std::acos(c), false};
}

template <std::floating_point T>
inline T cosine_alignment(std::span<const T> a,
                          std::span<const T> b,
                          T tol = default_tolerance<T>())
{
    return alignment<T>(a, b, tol).cosine;
}

template <std::floating_point T>
inline T absolute_alignment(std::span<const T> a,
                            std::span<const T> b,
                            T tol = default_tolerance<T>())
{
    return std::fabs(cosine_alignment<T>(a, b, tol));
}

template <std::floating_point T>
inline T angle_between(std::span<const T> a,
                       std::span<const T> b,
                       T tol = default_tolerance<T>())
{
    return alignment<T>(a, b, tol).angle;
}

template <std::floating_point T>
inline bool point_in_parent_cell(cell::type parent_cell_type,
                                 std::span<const T> x,
                                 T tol = default_tolerance<T>())
{
    const int tdim = cell::get_tdim(parent_cell_type);
    if (static_cast<int>(x.size()) != tdim)
        throw std::invalid_argument("point_in_parent_cell: point dimension mismatch");

    switch (parent_cell_type)
    {
    case cell::type::point:
        return x.empty();
    case cell::type::interval:
        return x[0] >= -tol && x[0] <= T(1) + tol;
    case cell::type::triangle:
        return x[0] >= -tol && x[1] >= -tol && x[0] + x[1] <= T(1) + tol;
    case cell::type::quadrilateral:
        return x[0] >= -tol && x[0] <= T(1) + tol
            && x[1] >= -tol && x[1] <= T(1) + tol;
    case cell::type::tetrahedron:
        return x[0] >= -tol && x[1] >= -tol && x[2] >= -tol
            && x[0] + x[1] + x[2] <= T(1) + tol;
    case cell::type::hexahedron:
        return x[0] >= -tol && x[0] <= T(1) + tol
            && x[1] >= -tol && x[1] <= T(1) + tol
            && x[2] >= -tol && x[2] <= T(1) + tol;
    case cell::type::prism:
        return x[0] >= -tol && x[1] >= -tol && x[0] + x[1] <= T(1) + tol
            && x[2] >= -tol && x[2] <= T(1) + tol;
    case cell::type::pyramid:
        return x[0] >= -tol && x[1] >= -tol && x[2] >= -tol
            && x[0] + x[2] <= T(1) + tol
            && x[1] + x[2] <= T(1) + tol;
    default:
        return false;
    }
}

template <std::floating_point T>
inline bool point_on_segment(std::span<const T> x,
                             std::span<const T> a,
                             std::span<const T> b,
                             T tol = default_tolerance<T>(),
                             T* parameter = nullptr)
{
    require_same_size<T>(x, a, "point_on_segment");
    require_same_size<T>(a, b, "point_on_segment");

    const auto ab = subtract<T>(b, a);
    const T ab2 = squared_norm<T>(std::span<const T>(ab.data(), ab.size()));
    if (ab2 <= tol * tol)
        return false;

    T ax_ab = T(0);
    for (std::size_t i = 0; i < x.size(); ++i)
        ax_ab += (x[i] - a[i]) * ab[i];
    const T t = ax_ab / ab2;
    if (parameter != nullptr)
        *parameter = t;

    if (t < -tol || t > T(1) + tol)
        return false;

    T distance2 = T(0);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const T closest = a[i] + t * ab[i];
        const T r = x[i] - closest;
        distance2 += r * r;
    }
    return distance2 <= tol * tol;
}

template <std::floating_point T>
inline bool point_on_parent_edge(cell::type parent_cell_type,
                                 int parent_edge_id,
                                 std::span<const T> x,
                                 T tol = default_tolerance<T>(),
                                 T* edge_parameter = nullptr)
{
    const auto edges = cell::edges(parent_cell_type);
    if (parent_edge_id < 0 || parent_edge_id >= static_cast<int>(edges.size()))
        return false;
    const auto edge = edges[static_cast<std::size_t>(parent_edge_id)];
    const auto a = reference_vertex<T>(parent_cell_type, edge[0]);
    const auto b = reference_vertex<T>(parent_cell_type, edge[1]);
    return point_on_segment<T>(
        x,
        std::span<const T>(a.data(), a.size()),
        std::span<const T>(b.data(), b.size()),
        tol,
        edge_parameter);
}

template <std::floating_point T>
inline bool point_on_parent_face(cell::type parent_cell_type,
                                 int parent_face_id,
                                 std::span<const T> x,
                                 T tol = default_tolerance<T>())
{
    if (cell::get_tdim(parent_cell_type) != 3)
        return false;
    if (parent_face_id < 0 || parent_face_id >= cell::num_faces(parent_cell_type))
        return false;

    const auto fv = cell::face_vertices(parent_cell_type, parent_face_id);
    if (fv.size() < 3)
        return false;

    const auto a = reference_vertex<T>(parent_cell_type, fv[0]);
    const auto normal = parent_face_normal<T>(parent_cell_type, parent_face_id, false, tol);
    if (normal.degenerate())
        return false;

    T signed_distance_num = T(0);
    for (std::size_t i = 0; i < x.size(); ++i)
        signed_distance_num += (x[i] - a[i]) * normal.value[i];
    if (std::fabs(signed_distance_num) > tol * normal.norm)
        return false;

    return point_in_parent_cell<T>(parent_cell_type, x, tol);
}

template <std::floating_point T>
inline bool point_in_parent_entity(cell::type parent_cell_type,
                                   ParentEntity entity,
                                   std::span<const T> x,
                                   T tol = default_tolerance<T>())
{
    const int tdim = cell::get_tdim(parent_cell_type);
    if (!entity.valid() || entity.dim > tdim)
        return false;

    if (entity.dim == tdim)
        return point_in_parent_cell<T>(parent_cell_type, x, tol);

    if (entity.dim == 2 && tdim == 3)
        return point_on_parent_face<T>(parent_cell_type, entity.id, x, tol);

    if (entity.dim == 1)
        return point_on_parent_edge<T>(parent_cell_type, entity.id, x, tol);

    if (entity.dim == 0)
    {
        const auto v = reference_vertex<T>(parent_cell_type, entity.id);
        const auto delta = subtract<T>(x, std::span<const T>(v.data(), v.size()));
        return norm<T>(std::span<const T>(delta.data(), delta.size())) <= tol;
    }

    return false;
}

template <std::floating_point T>
inline ParentEntity smallest_parent_entity_containing_point(
    cell::type parent_cell_type,
    std::span<const T> x,
    T tol = default_tolerance<T>())
{
    const int nverts = cell::get_num_vertices(parent_cell_type);
    for (int v = 0; v < nverts; ++v)
    {
        const ParentEntity entity{0, v};
        if (point_in_parent_entity<T>(parent_cell_type, entity, x, tol))
            return entity;
    }

    if (cell::get_tdim(parent_cell_type) >= 1)
    {
        const auto edges = cell::edges(parent_cell_type);
        for (int e = 0; e < static_cast<int>(edges.size()); ++e)
        {
            if (point_on_parent_edge<T>(parent_cell_type, e, x, tol))
                return {1, e};
        }
    }

    if (cell::get_tdim(parent_cell_type) == 3)
    {
        for (int f = 0; f < cell::num_faces(parent_cell_type); ++f)
        {
            if (point_on_parent_face<T>(parent_cell_type, f, x, tol))
                return {2, f};
        }
    }

    if (point_in_parent_cell<T>(parent_cell_type, x, tol))
        return {cell::get_tdim(parent_cell_type), -1};

    return {-1, -1};
}

template <std::floating_point T>
inline bool clip_halfspace(std::span<const T> x0,
                           std::span<const T> direction,
                           std::span<const T> normal,
                           T rhs,
                           LineInterval<T>& interval,
                           T tol)
{
    const T value0 = dot<T>(normal, x0);
    const T rate = dot<T>(normal, direction);
    if (std::fabs(rate) <= tol)
        return value0 <= rhs + tol;

    const T hit = (rhs - value0) / rate;
    if (rate > T(0))
        interval.t1 = std::min(interval.t1, hit);
    else
        interval.t0 = std::max(interval.t0, hit);

    return interval.t0 <= interval.t1 + tol;
}

template <std::floating_point T>
inline bool clip_coordinate_lower(int axis,
                                  std::span<const T> x0,
                                  std::span<const T> direction,
                                  LineInterval<T>& interval,
                                  T tol)
{
    std::vector<T> n(x0.size(), T(0));
    n[static_cast<std::size_t>(axis)] = T(-1);
    return clip_halfspace<T>(
        x0, direction, std::span<const T>(n.data(), n.size()), T(0), interval, tol);
}

template <std::floating_point T>
inline bool clip_coordinate_upper(int axis,
                                  T upper,
                                  std::span<const T> x0,
                                  std::span<const T> direction,
                                  LineInterval<T>& interval,
                                  T tol)
{
    std::vector<T> n(x0.size(), T(0));
    n[static_cast<std::size_t>(axis)] = T(1);
    return clip_halfspace<T>(
        x0, direction, std::span<const T>(n.data(), n.size()), upper, interval, tol);
}

template <std::floating_point T>
inline bool clip_sum_upper(std::span<const int> axes,
                           T upper,
                           std::span<const T> x0,
                           std::span<const T> direction,
                           LineInterval<T>& interval,
                           T tol)
{
    std::vector<T> n(x0.size(), T(0));
    for (const int axis : axes)
        n[static_cast<std::size_t>(axis)] = T(1);
    return clip_halfspace<T>(
        x0, direction, std::span<const T>(n.data(), n.size()), upper, interval, tol);
}

template <std::floating_point T>
inline LineInterval<T> clip_line_interval_in_parent_cell(
    cell::type parent_cell_type,
    std::span<const T> x0,
    std::span<const T> direction,
    T t0 = -std::numeric_limits<T>::infinity(),
    T t1 = std::numeric_limits<T>::infinity(),
    T tol = default_tolerance<T>())
{
    const int tdim = cell::get_tdim(parent_cell_type);
    if (static_cast<int>(x0.size()) != tdim
        || static_cast<int>(direction.size()) != tdim)
    {
        throw std::invalid_argument(
            "clip_line_interval_in_parent_cell: dimension mismatch");
    }

    LineInterval<T> interval{true, t0, t1};
    auto lower = [&](int axis)
    {
        return clip_coordinate_lower<T>(axis, x0, direction, interval, tol);
    };
    auto upper = [&](int axis, T value)
    {
        return clip_coordinate_upper<T>(axis, value, x0, direction, interval, tol);
    };
    auto sum_upper = [&](std::span<const int> axes, T value)
    {
        return clip_sum_upper<T>(axes, value, x0, direction, interval, tol);
    };

    switch (parent_cell_type)
    {
    case cell::type::point:
        interval.valid = x0.empty() && direction.empty();
        return interval;
    case cell::type::interval:
        interval.valid = lower(0) && upper(0, T(1));
        return interval;
    case cell::type::triangle:
    {
        const std::array<int, 2> axes = {0, 1};
        interval.valid = lower(0) && lower(1)
                      && sum_upper(std::span<const int>(axes), T(1));
        return interval;
    }
    case cell::type::quadrilateral:
        interval.valid = lower(0) && upper(0, T(1))
                      && lower(1) && upper(1, T(1));
        return interval;
    case cell::type::tetrahedron:
    {
        const std::array<int, 3> axes = {0, 1, 2};
        interval.valid = lower(0) && lower(1) && lower(2)
                      && sum_upper(std::span<const int>(axes), T(1));
        return interval;
    }
    case cell::type::hexahedron:
        interval.valid = lower(0) && upper(0, T(1))
                      && lower(1) && upper(1, T(1))
                      && lower(2) && upper(2, T(1));
        return interval;
    case cell::type::prism:
    {
        const std::array<int, 2> axes = {0, 1};
        interval.valid = lower(0) && lower(1)
                      && sum_upper(std::span<const int>(axes), T(1))
                      && lower(2) && upper(2, T(1));
        return interval;
    }
    case cell::type::pyramid:
    {
        const std::array<int, 2> xz = {0, 2};
        const std::array<int, 2> yz = {1, 2};
        interval.valid = lower(0) && lower(1) && lower(2)
                      && sum_upper(std::span<const int>(xz), T(1))
                      && sum_upper(std::span<const int>(yz), T(1));
        return interval;
    }
    default:
        interval.valid = false;
        return interval;
    }
}

template <std::floating_point T>
inline LineInterval<T> clip_line_interval_in_parent_entity(
    cell::type parent_cell_type,
    ParentEntity entity,
    std::span<const T> x0,
    std::span<const T> direction,
    T t0 = -std::numeric_limits<T>::infinity(),
    T t1 = std::numeric_limits<T>::infinity(),
    T tol = default_tolerance<T>())
{
    const int tdim = cell::get_tdim(parent_cell_type);
    if (!entity.valid() || entity.dim > tdim)
        return {false, T(0), T(0)};

    if (entity.dim == tdim)
        return clip_line_interval_in_parent_cell<T>(
            parent_cell_type, x0, direction, t0, t1, tol);

    if (entity.dim == 2 && tdim == 3)
    {
        if (!point_on_parent_face<T>(parent_cell_type, entity.id, x0, tol))
            return {false, T(0), T(0)};
        const auto normal = parent_face_normal<T>(parent_cell_type, entity.id, false, tol);
        if (normal.degenerate())
            return {false, T(0), T(0)};
        const T normal_rate = dot<T>(
            direction,
            std::span<const T>(normal.value.data(), normal.value.size()));
        if (std::fabs(normal_rate) > tol * std::max<T>(T(1), normal.norm))
            return {false, T(0), T(0)};
        return clip_line_interval_in_parent_cell<T>(
            parent_cell_type, x0, direction, t0, t1, tol);
    }

    if (entity.dim == 1)
    {
        T edge_s = T(0);
        if (!point_on_parent_edge<T>(parent_cell_type, entity.id, x0, tol, &edge_s))
            return {false, T(0), T(0)};

        const auto tangent = parent_edge_tangent<T>(parent_cell_type, entity.id, false, tol);
        if (tangent.degenerate())
            return {false, T(0), T(0)};

        const T tt = squared_norm<T>(
            std::span<const T>(tangent.value.data(), tangent.value.size()));
        const T ds = dot<T>(
            direction,
            std::span<const T>(tangent.value.data(), tangent.value.size())) / tt;
        std::vector<T> parallel(direction.size(), T(0));
        for (std::size_t i = 0; i < direction.size(); ++i)
            parallel[i] = ds * tangent.value[i];
        const auto residual = subtract<T>(
            direction, std::span<const T>(parallel.data(), parallel.size()));
        if (norm<T>(std::span<const T>(residual.data(), residual.size())) > tol)
            return {false, T(0), T(0)};

        if (std::fabs(ds) <= tol)
            return {true, t0, t1};
        const T a = (T(0) - edge_s) / ds;
        const T b = (T(1) - edge_s) / ds;
        return {std::max(t0, std::min(a, b)) <= std::min(t1, std::max(a, b)) + tol,
                std::max(t0, std::min(a, b)),
                std::min(t1, std::max(a, b))};
    }

    if (entity.dim == 0)
    {
        if (!point_in_parent_entity<T>(parent_cell_type, entity, x0, tol))
            return {false, T(0), T(0)};
        if (norm<T>(direction) <= tol)
            return {true, t0, t1};
        return {true, T(0), T(0)};
    }

    return {false, T(0), T(0)};
}

} // namespace cutcells::geom
