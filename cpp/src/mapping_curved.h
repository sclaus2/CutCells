// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

#pragma once

#include "local_level_set.h"

#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <span>
#include <vector>

namespace cutcells::mapping
{

/// Backend selection for curved volume mappings.
enum class CurvedMappingBackend : uint8_t
{
    collapsed    = 0,   ///< Current collapsed-triangle blend: F = (1-t)*gamma(s) + t*v_opp
    gordon_hall  = 1,   ///< Gordon–Hall transfinite interpolation
};

// ============================================================================
// Pk Lagrange basis on [0,1] with Gauss-Lobatto nodes
// ============================================================================

/// Evaluate all p+1 Lagrange basis functions at parameter s in [0,1].
/// gl_interior: the p-1 interior GL points in (0,1) (ascending order).
/// out: output array of size p+1.
template <std::floating_point T>
inline void eval_lagrange_basis_gl(T s,
                                   const T* gl_interior,
                                   int p,
                                   T* out)
{
    const int n_nodes = p + 1;

    auto node_pos = [&](int j) -> T {
        if (j == 0)   return T(0);
        if (j == p)   return T(1);
        return gl_interior[j - 1];
    };

    for (int j = 0; j < n_nodes; ++j)
    {
        T tj = node_pos(j);
        T num = T(1);
        T den = T(1);
        for (int k = 0; k < n_nodes; ++k)
        {
            if (k == j) continue;
            T tk = node_pos(k);
            num *= (s - tk);
            den *= (tj - tk);
        }
        out[j] = (std::abs(den) < T(1e-30)) ? T(0) : num / den;
    }
}

/// Evaluate derivative of all p+1 Lagrange basis functions at parameter s.
template <std::floating_point T>
inline void eval_lagrange_deriv_gl(T s,
                                   const T* gl_interior,
                                   int p,
                                   T* dout)
{
    const int n_nodes = p + 1;

    auto node_pos = [&](int j) -> T {
        if (j == 0)   return T(0);
        if (j == p)   return T(1);
        return gl_interior[j - 1];
    };

    for (int j = 0; j < n_nodes; ++j)
    {
        T tj = node_pos(j);
        T den = T(1);
        for (int k = 0; k < n_nodes; ++k)
        {
            if (k == j) continue;
            den *= (tj - node_pos(k));
        }

        T sum = T(0);
        for (int m = 0; m < n_nodes; ++m)
        {
            if (m == j) continue;
            T prod = T(1);
            for (int k = 0; k < n_nodes; ++k)
            {
                if (k == j || k == m) continue;
                prod *= (s - node_pos(k));
            }
            sum += prod;
        }
        dout[j] = (std::abs(den) < T(1e-30)) ? T(0) : sum / den;
    }
}

// ============================================================================
// Pk curve evaluation
// ============================================================================

/// Evaluate the Pk-curved position x(s) on one interface entity.
/// s in [0,1] is the parameter along the segment.
/// Nodes: p+1 total — v0 (endpoint), interior GL nodes, v1 (endpoint).
///
/// @param s          parameter in [0,1]
/// @param x_v0      coordinates of endpoint v0 (dim components)
/// @param x_v1      coordinates of endpoint v1 (dim components)
/// @param x_int     flat coordinates of p-1 interior nodes (dim each)
/// @param gl_pts    interior GL points in (0,1)
/// @param p         polynomial order (>= 1)
/// @param dim       number of coordinate components
/// @param x_out     output: dim coordinates
template <std::floating_point T>
inline void eval_pk_map(T s,
                        const T* x_v0,
                        const T* x_v1,
                        const T* x_int,
                        const T* gl_pts,
                        int p,
                        int dim,
                        T* x_out)
{
    const int n_nodes = p + 1;
    T L[16] = {};
    assert(n_nodes <= 16);

    eval_lagrange_basis_gl(s, gl_pts, p, L);

    for (int d = 0; d < dim; ++d)
        x_out[d] = T(0);

    for (int d = 0; d < dim; ++d)
        x_out[d] += L[0] * x_v0[d];

    for (int j = 1; j < p; ++j)
        for (int d = 0; d < dim; ++d)
            x_out[d] += L[j] * x_int[(j - 1) * dim + d];

    for (int d = 0; d < dim; ++d)
        x_out[d] += L[p] * x_v1[d];
}

/// Evaluate dx/ds (derivative of Pk map w.r.t. parameter s).
template <std::floating_point T>
inline void eval_pk_deriv(T s,
                          const T* x_v0,
                          const T* x_v1,
                          const T* x_int,
                          const T* gl_pts,
                          int p,
                          int dim,
                          T* dx_out)
{
    const int n_nodes = p + 1;
    T dL[16] = {};
    assert(n_nodes <= 16);

    eval_lagrange_deriv_gl(s, gl_pts, p, dL);

    for (int d = 0; d < dim; ++d)
        dx_out[d] = T(0);

    for (int d = 0; d < dim; ++d)
        dx_out[d] += dL[0] * x_v0[d];

    for (int j = 1; j < p; ++j)
        for (int d = 0; d < dim; ++d)
            dx_out[d] += dL[j] * x_int[(j - 1) * dim + d];

    for (int d = 0; d < dim; ++d)
        dx_out[d] += dL[p] * x_v1[d];
}

// ============================================================================
// Collapsed triangle mapping (existing approach)
// ============================================================================

/// @brief Evaluate the collapsed triangle mapping in reference space.
///
/// F_ref(s, t) = (1 - t) * gamma_ref(s) + t * r_opp
///
/// where gamma_ref is the Pk curved interface map parametrised by s in [0,1],
/// and t in [0,1] blends from the interface (t=0) to the opposite vertex (t=1).
///
/// @param s             parameter along interface edge [0,1]
/// @param t             blending parameter [0,1]: 0=interface, 1=opposite vertex
/// @param r_va          reference coords of interface endpoint v_a (tdim)
/// @param r_vb          reference coords of interface endpoint v_b (tdim)
/// @param x_int_ref     flat interior curved node coords in reference space
/// @param gl_pts        interior GL points in (0,1)
/// @param geom_order    polynomial order of Pk curve
/// @param r_opp         reference coords of opposite vertex (tdim)
/// @param tdim          topological dimension
/// @param F_ref         output: mapped reference position (tdim)
/// @param col_s_ref     output: dF/ds in reference space (tdim)
/// @param col_t_ref     output: dF/dt in reference space (tdim)
template <std::floating_point T>
inline void eval_collapsed_triangle(T s, T t,
                                    const T* r_va,
                                    const T* r_vb,
                                    const T* x_int_ref,
                                    const T* gl_pts,
                                    int geom_order,
                                    const T* r_opp,
                                    int tdim,
                                    T* F_ref,
                                    T* col_s_ref,
                                    T* col_t_ref)
{
    T gamma_ref[3] = {};
    T dgamma_ref[3] = {};
    eval_pk_map(s, r_va, r_vb, x_int_ref, gl_pts,
                geom_order, tdim, gamma_ref);
    eval_pk_deriv(s, r_va, r_vb, x_int_ref, gl_pts,
                  geom_order, tdim, dgamma_ref);

    for (int k = 0; k < tdim; ++k)
    {
        F_ref[k]     = (T(1) - t) * gamma_ref[k] + t * r_opp[k];
        col_s_ref[k] = (T(1) - t) * dgamma_ref[k];
        col_t_ref[k] = r_opp[k] - gamma_ref[k];
    }
}

// ============================================================================
// Gordon-Hall transfinite interpolation for 2D cells
// ============================================================================

/// @brief Evaluate the Gordon-Hall mapping for a triangle with one curved edge.
///
/// The triangle has three edges. Edge 0 (from v0 to v1) is curved via the Pk
/// map gamma(s); edges 1 and 2 are straight. The Gordon-Hall transfinite
/// interpolation for this case is:
///
///   F(xi, eta) = (1 - eta) * gamma(xi) + xi * v1 + eta * v2
///              - (1 - eta) * [(1 - xi) * v0 + xi * v1]
///              - eta * [(1 - xi) * v2 + xi * v2]   ... (simplified)
///
/// Using standard triangle coordinates (xi, eta) in [0,1] x [0,1-xi]:
///   - Edge 0: eta = 0, xi in [0,1]   → v0 to v1 (curved)
///   - Edge 1: xi + eta = 1 (on tri)  → v1 to v2 (straight)
///   - Edge 2: xi = 0                  → v0 to v2 (straight)
///
/// The GH formula for a triangle with one curved edge (edge 0) reduces to:
///
///   F(xi, eta) = (1 - eta) * gamma(xi) + eta * v2 + (1 - eta) * [xi * (v1 - gamma(1)) + (1 - xi) * (v0 - gamma(0))]
///
/// Since gamma(0) = v0 and gamma(1) = v1 (interpolatory endpoints), this
/// simplifies to:
///
///   F(xi, eta) = (1 - eta) * gamma(xi) + eta * v2
///
/// Jacobian columns:
///   dF/dxi  = (1 - eta) * dgamma/dxi
///   dF/deta = v2 - gamma(xi)
///
/// This is identical to the collapsed mapping when we identify:
///   xi = s, eta = t, v2 = v_opp
///
/// @param xi            first barycentric-like coordinate [0,1]
/// @param eta           second coordinate [0,1-xi]
/// @param r_v0          reference coords of vertex 0 (start of curved edge)
/// @param r_v1          reference coords of vertex 1 (end of curved edge)
/// @param x_int_ref     flat interior curved node coords in ref space
/// @param gl_pts        interior GL points in (0,1)
/// @param geom_order    polynomial order of Pk curve
/// @param r_v2          reference coords of vertex 2 (opposite to curved edge)
/// @param tdim          topological dimension
/// @param F_ref         output: mapped reference position
/// @param dF_dxi        output: dF/dxi in reference space
/// @param dF_deta       output: dF/deta in reference space
template <std::floating_point T>
inline void eval_gh_triangle_one_curved_edge(
    T xi, T eta,
    const T* r_v0,
    const T* r_v1,
    const T* x_int_ref,
    const T* gl_pts,
    int geom_order,
    const T* r_v2,
    int tdim,
    T* F_ref,
    T* dF_dxi,
    T* dF_deta)
{
    // For a triangle with one curved edge (edge 0: v0→v1),
    // F(xi, eta) = (1 - eta) * gamma(xi) + eta * v2
    //
    // This is mathematically identical to the collapsed mapping,
    // with the correspondence: s = xi, t = eta, v_opp = v2.
    //
    // The GH formulation makes the generalisation to quads and
    // multiple curved edges cleaner.

    T gamma_ref[3] = {};
    T dgamma_ref[3] = {};
    eval_pk_map(xi, r_v0, r_v1, x_int_ref, gl_pts,
                geom_order, tdim, gamma_ref);
    eval_pk_deriv(xi, r_v0, r_v1, x_int_ref, gl_pts,
                  geom_order, tdim, dgamma_ref);

    for (int k = 0; k < tdim; ++k)
    {
        F_ref[k]    = (T(1) - eta) * gamma_ref[k] + eta * r_v2[k];
        dF_dxi[k]   = (T(1) - eta) * dgamma_ref[k];
        dF_deta[k]  = r_v2[k] - gamma_ref[k];
    }
}

/// @brief Evaluate the Gordon-Hall mapping for a quadrilateral with one
///        curved edge.
///
/// Quadrilateral with vertices v0, v1, v2, v3 in order:
///   - Edge 0: v0 → v1 (curved via Pk map gamma)
///   - Edge 1: v1 → v2 (straight)
///   - Edge 2: v3 → v2 (straight, opposite to edge 0)
///   - Edge 3: v0 → v3 (straight)
///
/// Using (u, v) in [0,1]^2:
///   u runs along edges 0 and 2 (from v0-side to v1-side)
///   v runs along edges 3 and 1 (from edge-0-side to edge-2-side)
///
/// Gordon-Hall formula with one curved edge (edge 0 at v=0):
///
///   F(u, v) = (1-v) * gamma(u)        + v * edge2(u)
///           + (1-u) * edge3(v)         + u * edge1(v)
///           - (1-u)*(1-v)*v0 - u*(1-v)*v1 - u*v*v2 - (1-u)*v*v3
///
/// where:
///   gamma(u) = Pk map from v0 to v1
///   edge2(u) = (1-u)*v3 + u*v2    (straight)
///   edge3(v) = (1-v)*v0 + v*v3    (straight)
///   edge1(v) = (1-v)*v1 + v*v2    (straight)
///
/// Substituting and simplifying (since gamma(0)=v0, gamma(1)=v1):
///
///   F(u, v) = (1-v) * gamma(u) + v * [(1-u)*v3 + u*v2]
///
/// dF/du = (1-v) * dgamma/du  + v * (v2 - v3)
/// dF/dv = -gamma(u)  + (1-u)*v3 + u*v2
///       = [(1-u)*v3 + u*v2] - gamma(u)
///
/// @param u             parameter along curved edge [0,1]
/// @param v             parameter across quad [0,1]: 0=curved edge, 1=opposite
/// @param r_v0          reference coords of v0 (start of curved edge)
/// @param r_v1          reference coords of v1 (end of curved edge)
/// @param x_int_ref     flat interior curved node coords in ref space
/// @param gl_pts        interior GL points in (0,1)
/// @param geom_order    polynomial order of Pk curve
/// @param r_v2          reference coords of v2 (diagonal from v0)
/// @param r_v3          reference coords of v3 (opposite v1 on curved-edge side)
/// @param tdim          topological dimension
/// @param F_ref         output: mapped reference position
/// @param dF_du         output: dF/du in reference space
/// @param dF_dv         output: dF/dv in reference space
template <std::floating_point T>
inline void eval_gh_quad_one_curved_edge(
    T u, T v,
    const T* r_v0,
    const T* r_v1,
    const T* x_int_ref,
    const T* gl_pts,
    int geom_order,
    const T* r_v2,
    const T* r_v3,
    int tdim,
    T* F_ref,
    T* dF_du,
    T* dF_dv)
{
    T gamma_ref[3] = {};
    T dgamma_ref[3] = {};
    eval_pk_map(u, r_v0, r_v1, x_int_ref, gl_pts,
                geom_order, tdim, gamma_ref);
    eval_pk_deriv(u, r_v0, r_v1, x_int_ref, gl_pts,
                  geom_order, tdim, dgamma_ref);

    // Opposite edge position: straight from v3 to v2
    // edge2(u) = (1-u)*v3 + u*v2
    for (int k = 0; k < tdim; ++k)
    {
        const T opp_k = (T(1) - u) * r_v3[k] + u * r_v2[k];

        F_ref[k]  = (T(1) - v) * gamma_ref[k] + v * opp_k;
        dF_du[k]  = (T(1) - v) * dgamma_ref[k] + v * (r_v2[k] - r_v3[k]);
        dF_dv[k]  = opp_k - gamma_ref[k];
    }
}

// ============================================================================
// Gordon-Hall transfinite interpolation for 3D cells
// ============================================================================

template <std::floating_point T>
inline void eval_gh_face_triangle(
    std::span<const BernsteinCell<T>> coord_polys,
    T                                 xi,
    T                                 eta,
    T*                                F_ref,
    T*                                dF_dxi,
    T*                                dF_deta);

/// @brief Evaluate the GH mapping for a tetrahedron with one curved face.
///
/// The tetrahedron has face 0 (v0, v1, v2) curved. The opposite vertex is v3.
/// Using barycentric-like coordinates (xi, eta, zeta) in the reference tet:
///   - Face 0: zeta = 0 → triangle (v0, v1, v2) in (xi, eta)
///   - Vertex 3 at zeta = 1
///
/// GH blending with one curved face:
///   F(xi, eta, zeta) = (1 - zeta) * F_face(xi', eta') + zeta * v3
///
/// where (xi', eta') = (xi / (1 - zeta), eta / (1 - zeta)) are the face
/// coordinates rescaled to [0,1] range, and F_face is the curved face map.
///
/// For a triangular face curved via a Pk edge (edge v0→v1) with the face
/// opposite vertex v2, we have:
///   F_face(xi', eta') = (1 - eta') * gamma(xi') + eta' * v2
///
/// So the full tet mapping is:
///   F(xi, eta, zeta) = (1 - zeta) * [(1 - eta') * gamma(xi') + eta' * v2] + zeta * v3
///
/// where xi' = xi / (1 - zeta), eta' = eta / (1 - zeta), and the denominator
/// is clamped to avoid division by zero at the apex.
///
/// @param xi, eta, zeta  reference coordinates in [0,1]^3 (tet domain)
/// @param r_v0, r_v1     endpoints of curved edge on face 0
/// @param x_int_ref      interior GL nodes of curved edge (p-1 nodes, tdim coords each)
/// @param gl_pts         interior GL points in (0,1)
/// @param geom_order     polynomial order of curved edge
/// @param r_v2           third vertex of curved face (face-opposite to curved edge)
/// @param r_v3           opposite vertex (apex)
/// @param tdim           topological dimension (3)
/// @param F_ref          output: mapped position (tdim)
/// @param dF_dxi         output: dF/dxi (tdim)
/// @param dF_deta        output: dF/deta (tdim)
/// @param dF_dzeta       output: dF/dzeta (tdim)
template <std::floating_point T>
inline void eval_gh_tet_one_curved_face(
    T xi, T eta, T zeta,
    const T* r_v0,
    const T* r_v1,
    const T* x_int_ref,
    const T* gl_pts,
    int geom_order,
    const T* r_v2,
    const T* r_v3,
    int tdim,
    T* F_ref,
    T* dF_dxi,
    T* dF_deta,
    T* dF_dzeta)
{
    const T eps = static_cast<T>(1e-14);
    const T one_m_zeta = std::max(T(1) - zeta, eps);

    // Rescaled face coordinates
    const T xip  = xi  / one_m_zeta;
    const T etap = eta / one_m_zeta;

    // Evaluate Pk curved edge map on face: gamma(xip) and dgamma/dxip
    T gamma[3] = {};
    T dgamma[3] = {};
    eval_pk_map(xip, r_v0, r_v1, x_int_ref, gl_pts, geom_order, tdim, gamma);
    eval_pk_deriv(xip, r_v0, r_v1, x_int_ref, gl_pts, geom_order, tdim, dgamma);

    // F_face(xip, etap) = (1 - etap) * gamma(xip) + etap * v2
    // F(xi, eta, zeta)  = (1 - zeta) * F_face(xip, etap) + zeta * v3

    for (int k = 0; k < tdim; ++k)
    {
        const T face_k = (T(1) - etap) * gamma[k] + etap * r_v2[k];
        F_ref[k] = (T(1) - zeta) * face_k + zeta * r_v3[k];

        // Partial derivatives via chain rule.
        // xip = xi / (1-zeta), etap = eta / (1-zeta)
        //
        // dF/dxi = (1-zeta) * dF_face/dxip * dxip/dxi
        //        = (1-zeta) * (1 - etap) * dgamma(xip) * (1/(1-zeta))
        //        = (1 - etap) * dgamma(xip)
        dF_dxi[k] = (T(1) - etap) * dgamma[k];

        // dF/deta = (1-zeta) * dF_face/detap * detap/deta
        //         = (1-zeta) * (v2 - gamma(xip)) * (1/(1-zeta))
        //         = v2 - gamma(xip)
        dF_deta[k] = r_v2[k] - gamma[k];

        // dF/dzeta = -F_face + v3
        //          + (1-zeta) * [ dF_face/dxip * dxip/dzeta + dF_face/detap * detap/dzeta ]
        //
        // dxip/dzeta  = xi / (1-zeta)^2
        // detap/dzeta = eta / (1-zeta)^2
        //
        // dF_face/dxip  = (1 - etap) * dgamma(xip)
        // dF_face/detap = v2 - gamma(xip)
        const T dface_dxip_k  = (T(1) - etap) * dgamma[k];
        const T dface_detap_k = r_v2[k] - gamma[k];
        const T dxip_dzeta  = xi  / (one_m_zeta * one_m_zeta);
        const T detap_dzeta = eta / (one_m_zeta * one_m_zeta);

        dF_dzeta[k] = r_v3[k] - face_k
                     + one_m_zeta * (dface_dxip_k * dxip_dzeta
                                   + dface_detap_k * detap_dzeta);
    }
}

template <std::floating_point T>
inline void eval_gh_tet_one_curved_tri_face(
    T                                 xi,
    T                                 eta,
    T                                 zeta,
    std::span<const BernsteinCell<T>> face_coord_polys,
    const T*                          r_apex,
    T*                                F_ref,
    T*                                dF_dxi,
    T*                                dF_deta,
    T*                                dF_dzeta)
{
    const T eps = static_cast<T>(1e-14);
    const T one_m_zeta = std::max(T(1) - zeta, eps);
    const T xip = xi / one_m_zeta;
    const T etap = eta / one_m_zeta;

    std::array<T, 3> face_ref = {T(0), T(0), T(0)};
    std::array<T, 3> dface_dxi_p = {T(0), T(0), T(0)};
    std::array<T, 3> dface_deta_p = {T(0), T(0), T(0)};
    eval_gh_face_triangle<T>(
        face_coord_polys, xip, etap,
        face_ref.data(), dface_dxi_p.data(), dface_deta_p.data());

    const T dxip_dzeta  = xi / (one_m_zeta * one_m_zeta);
    const T detap_dzeta = eta / (one_m_zeta * one_m_zeta);

    for (std::size_t k = 0; k < face_coord_polys.size(); ++k)
    {
        F_ref[k] = one_m_zeta * face_ref[k] + zeta * r_apex[k];
        dF_dxi[k] = dface_dxi_p[k];
        dF_deta[k] = dface_deta_p[k];
        dF_dzeta[k] = r_apex[k] - face_ref[k]
                    + one_m_zeta
                        * (dface_dxi_p[k] * dxip_dzeta
                           + dface_deta_p[k] * detap_dzeta);
    }
}

/// @brief Evaluate the GH mapping for a prism with one curved triangular face.
///
/// Prism with triangle base face (v0, v1, v2) at w=0 and top face (v3, v4, v5)
/// at w=1. The bottom triangular face is curved via the Pk edge v0→v1
/// (with face-opposite vertex v2). Top face is straight.
///
/// GH blending:
///   F(xi, eta, w) = (1-w) * F_face_bot(xi, eta) + w * F_face_top(xi, eta)
///
/// where F_face_bot = (1-eta) * gamma(xi) + eta * v2
///       F_face_top = (1-eta) * [(1-xi)*v3 + xi*v4] + eta * v5  (straight)
///
/// @param xi, eta, w     reference coordinates: (xi, eta) on triangle, w in [0,1]
/// @param r_v0..r_v5     reference coords of 6 prism vertices
/// @param x_int_ref      interior GL nodes of curved edge (bottom face)
/// @param gl_pts         interior GL points in (0,1)
/// @param geom_order     polynomial order
/// @param tdim           topological dimension (3)
/// @param F_ref          output: mapped position (tdim)
/// @param dF_dxi         output: dF/dxi (tdim)
/// @param dF_deta        output: dF/deta (tdim)
/// @param dF_dw          output: dF/dw (tdim)
template <std::floating_point T>
inline void eval_gh_prism_one_curved_tri_face(
    T xi, T eta, T w,
    const T* r_v0,
    const T* r_v1,
    const T* r_v2,
    const T* r_v3,
    const T* r_v4,
    const T* r_v5,
    const T* x_int_ref,
    const T* gl_pts,
    int geom_order,
    int tdim,
    T* F_ref,
    T* dF_dxi,
    T* dF_deta,
    T* dF_dw)
{
    // Bottom face (curved): F_bot = (1-eta)*gamma(xi) + eta*v2
    T gamma[3] = {};
    T dgamma[3] = {};
    eval_pk_map(xi, r_v0, r_v1, x_int_ref, gl_pts, geom_order, tdim, gamma);
    eval_pk_deriv(xi, r_v0, r_v1, x_int_ref, gl_pts, geom_order, tdim, dgamma);

    for (int k = 0; k < tdim; ++k)
    {
        const T bot_k = (T(1) - eta) * gamma[k] + eta * r_v2[k];
        const T top_k = (T(1) - eta) * ((T(1) - xi) * r_v3[k] + xi * r_v4[k])
                       + eta * r_v5[k];

        F_ref[k] = (T(1) - w) * bot_k + w * top_k;

        // dF/dxi = (1-w) * (1-eta) * dgamma/dxi + w * (1-eta) * (v4 - v3)
        dF_dxi[k] = (T(1) - w) * (T(1) - eta) * dgamma[k]
                   + w * (T(1) - eta) * (r_v4[k] - r_v3[k]);

        // dF/deta = (1-w) * (v2 - gamma(xi)) + w * (v5 - (1-xi)*v3 - xi*v4)
        dF_deta[k] = (T(1) - w) * (r_v2[k] - gamma[k])
                    + w * (r_v5[k] - (T(1) - xi) * r_v3[k] - xi * r_v4[k]);

        // dF/dw = top - bot
        dF_dw[k] = top_k - bot_k;
    }
}

// ============================================================================
// 3D determinant helper
// ============================================================================

/// @brief Compute 3x3 determinant from three column vectors (each 3 components).
template <std::floating_point T>
inline T det_3d(const T* c0, const T* c1, const T* c2)
{
    return c0[0] * (c1[1] * c2[2] - c1[2] * c2[1])
         - c0[1] * (c1[0] * c2[2] - c1[2] * c2[0])
         + c0[2] * (c1[0] * c2[1] - c1[1] * c2[0]);
}

// ============================================================================
// Determinant helpers
// ============================================================================

/// @brief Compute 2D determinant of two column vectors (each gdim components).
template <std::floating_point T>
inline T det_2d(const T* col0, const T* col1, int gdim)
{
    if (gdim == 2)
    {
        return col0[0] * col1[1] - col0[1] * col1[0];
    }
    else if (gdim == 3)
    {
        // Cross product magnitude for embedded 2D surface in 3D
        const T cx = col0[1] * col1[2] - col0[2] * col1[1];
        const T cy = col0[2] * col1[0] - col0[0] * col1[2];
        const T cz = col0[0] * col1[1] - col0[1] * col1[0];
        return std::sqrt(cx * cx + cy * cy + cz * cz);
    }
    return T(0);
}

// ============================================================================
// 3D curved face maps
// ============================================================================

template <std::floating_point T>
inline void eval_gh_face_triangle(
    std::span<const BernsteinCell<T>> coord_polys,
    T                                 xi,
    T                                 eta,
    T*                                F_ref,
    T*                                dF_dxi,
    T*                                dF_deta)
{
    assert(!coord_polys.empty());
    assert(coord_polys[0].cell_type == cell::type::triangle);

    const std::array<T, 2> x_face = {xi, eta};
    for (std::size_t d = 0; d < coord_polys.size(); ++d)
    {
        F_ref[d] = evaluate_bernstein_cell<T>(
            coord_polys[d], std::span<const T>(x_face.data(), x_face.size()));

        std::array<T, 2> grad = {T(0), T(0)};
        evaluate_bernstein_cell_gradient<T>(
            coord_polys[d],
            std::span<const T>(x_face.data(), x_face.size()),
            std::span<T>(grad.data(), grad.size()));
        dF_dxi[d] = grad[0];
        dF_deta[d] = grad[1];
    }
}

template <std::floating_point T>
inline void eval_gh_face_quad(
    std::span<const BernsteinCell<T>> coord_polys,
    T                                 u,
    T                                 v,
    T*                                F_ref,
    T*                                dF_du,
    T*                                dF_dv)
{
    assert(!coord_polys.empty());
    assert(coord_polys[0].cell_type == cell::type::quadrilateral);

    const std::array<T, 2> x_face = {u, v};
    for (std::size_t d = 0; d < coord_polys.size(); ++d)
    {
        F_ref[d] = evaluate_bernstein_cell<T>(
            coord_polys[d], std::span<const T>(x_face.data(), x_face.size()));

        std::array<T, 2> grad = {T(0), T(0)};
        evaluate_bernstein_cell_gradient<T>(
            coord_polys[d],
            std::span<const T>(x_face.data(), x_face.size()),
            std::span<T>(grad.data(), grad.size()));
        dF_du[d] = grad[0];
        dF_dv[d] = grad[1];
    }
}

template <std::floating_point T>
inline T face_surface_metric(
    std::span<const T> J_parent,
    const T*           dF_du_ref,
    const T*           dF_dv_ref,
    int                gdim,
    int                tdim)
{
    std::array<T, 3> col0 = {T(0), T(0), T(0)};
    std::array<T, 3> col1 = {T(0), T(0), T(0)};
    for (int d = 0; d < gdim; ++d)
    {
        for (int k = 0; k < tdim; ++k)
        {
            col0[static_cast<std::size_t>(d)]
                += J_parent[static_cast<std::size_t>(d * tdim + k)] * dF_du_ref[k];
            col1[static_cast<std::size_t>(d)]
                += J_parent[static_cast<std::size_t>(d * tdim + k)] * dF_dv_ref[k];
        }
    }
    return det_2d(col0.data(), col1.data(), gdim);
}

} // namespace cutcells::mapping
