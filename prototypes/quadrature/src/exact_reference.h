// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier: MIT

// Exact reference values (about 1e-14) for a sphere and a convex polytope:
//   area of the sphere inside the polytope, and volume of the ball inside it.
//   area:   spherical coordinates in a frame whose poles are far from the polytope; for each
//           meridian phi the admissible theta-set is a union of intervals in closed form, so
//           the inner integral of sin(theta) is exact; the phi integral is split at every
//           corner and tangency of the boundary arcs and done with tanh-sinh.
//   volume: 3 Vol(B n P) = r Area(S n P) + sum_i d_i Area(F_i n B) (divergence theorem with x),
//           with exact polygon-disk areas.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace cutcells::proto::exact
{
using V3 = std::array<double, 3>;
inline V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline V3 add(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline V3 mul(double s, const V3& a) { return {s * a[0], s * a[1], s * a[2]}; }
inline double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline V3 cross(const V3& a, const V3& b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(const V3& a) { return std::sqrt(dot(a, a)); }
inline V3 unit(const V3& a) { return mul(1.0 / norm(a), a); }

// Convex polytope: faces as vertex loops (any orientation); outward normals are computed from the
// polytope centroid. Coordinates are relative to the sphere centre.
struct Face
{
  std::vector<V3> loop;
  V3 n;     // outward unit normal
  double d; // n . x = d on the face
};

inline std::vector<Face> make_faces(const std::vector<std::vector<V3>>& loops)
{
  V3 c{0, 0, 0};
  int cnt = 0;
  for (const auto& l : loops)
    for (const auto& v : l)
    {
      c = add(c, v);
      ++cnt;
    }
  c = mul(1.0 / cnt, c);
  std::vector<Face> faces;
  for (const auto& l : loops)
  {
    V3 nn{0, 0, 0}; // Newell normal
    for (size_t i = 0; i < l.size(); ++i)
      nn = add(nn, cross(l[i], l[(i + 1) % l.size()]));
    V3 n = unit(nn);
    double d = dot(n, l[0]);
    if (dot(n, c) > d) // centroid must be inside: n.c < d
    {
      n = mul(-1, n);
      d = -d;
    }
    faces.push_back({l, n, d});
  }
  return faces;
}

// ---- tanh-sinh on [a, b] ----
template <typename F>
double tanh_sinh(F&& f, double a, double b)
{
  const double hs = 1.0 / 32.0, half = 0.5 * (b - a), mid = 0.5 * (a + b);
  double s = 0;
  for (int k = -120; k <= 120; ++k)
  {
    const double t = k * hs;
    const double u = 0.5 * M_PI * std::sinh(t);
    const double ch = std::cosh(u);
    const double x = std::tanh(u);
    const double w = 0.5 * M_PI * std::cosh(t) / (ch * ch);
    if (w < 1e-300)
      continue;
    // distance to the endpoints computed without cancellation
    const double one_minus = 1.0 / (std::exp(u) * ch); // 1 - tanh(u)
    const double one_plus = std::exp(u) / ch;          // 1 + tanh(u)
    double xx;
    if (x > 0)
      xx = b - half * one_minus;
    else
      xx = a + half * one_plus;
    if (xx <= a || xx >= b)
      continue;
    (void)mid;
    s += w * f(xx);
  }
  return s * hs * half;
}

// ---- intervals ----
using Iv = std::pair<double, double>;
inline std::vector<Iv> intersect(const std::vector<Iv>& A, const std::vector<Iv>& B)
{
  std::vector<Iv> r;
  for (const auto& a : A)
    for (const auto& b : B)
    {
      const double lo = std::max(a.first, b.first), hi = std::min(a.second, b.second);
      if (hi > lo)
        r.push_back({lo, hi});
    }
  return r;
}

// Area of the sphere |x| = r inside the polytope.
inline double sphere_area(const std::vector<Face>& faces, double r)
{
  // frame: e_x towards the polytope centroid, poles far away
  V3 c{0, 0, 0};
  int cnt = 0;
  for (const auto& f : faces)
    for (const auto& v : f.loop)
    {
      c = add(c, v);
      ++cnt;
    }
  c = mul(1.0 / cnt, c);
  const V3 ex = unit(c);
  V3 tmp = std::abs(ex[0]) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
  const V3 ez = unit(cross(ex, tmp));
  const V3 ey = cross(ez, ex);
  struct C
  {
    double nx, ny, nz, s;
  };
  std::vector<C> cs;
  for (const auto& f : faces)
    cs.push_back({dot(f.n, ex), dot(f.n, ey), dot(f.n, ez), f.d / r});
  auto phi_of = [&](const V3& x) { return std::atan2(dot(x, ey), dot(x, ex)); };
  auto inside = [&](const V3& x, double tol)
  {
    for (const auto& f : faces)
      if (dot(f.n, x) > f.d + tol)
        return false;
    return true;
  };

  // admissible theta-set on meridian phi
  auto theta_set = [&](double ph)
  {
    std::vector<Iv> S = {{0.0, M_PI}};
    const double cp = std::cos(ph), sp = std::sin(ph);
    for (const auto& k : cs)
    {
      const double a = k.nx * cp + k.ny * sp, b = k.nz; // a sin(t) + b cos(t) <= s
      const double R = std::hypot(a, b);
      std::vector<Iv> allow;
      if (k.s >= R)
        continue;
      if (k.s < -R)
        return std::vector<Iv>{};
      const double beta = std::atan2(a, b), gam = std::acos(k.s / R);
      for (int m = -2; m <= 2; ++m)
      {
        const double lo = beta + gam + 2 * M_PI * m, hi = beta + 2 * M_PI - gam + 2 * M_PI * m;
        const double l2 = std::max(lo, 0.0), h2 = std::min(hi, M_PI);
        if (h2 > l2)
          allow.push_back({l2, h2});
      }
      S = intersect(S, allow);
      if (S.empty())
        return S;
    }
    return S;
  };
  auto fphi = [&](double ph)
  {
    double s = 0;
    for (const auto& iv : theta_set(ph))
      s += std::cos(iv.first) - std::cos(iv.second);
    return s;
  };

  // breakpoints: corners (sphere on an edge line) and tangencies of each circle with meridians
  std::vector<double> bp;
  const double tol = 1e-12 * r;
  for (size_t i = 0; i < faces.size(); ++i)
    for (size_t j = i + 1; j < faces.size(); ++j)
    {
      const V3 u = cross(faces[i].n, faces[j].n);
      const double uu = dot(u, u);
      if (uu < 1e-24)
        continue;
      // point on both planes: p = (d_i (n_j x u) + d_j (u x n_i)) / |u|^2
      const V3 p = mul(1.0 / uu, add(mul(faces[i].d, cross(faces[j].n, u)), mul(faces[j].d, cross(u, faces[i].n))));
      // |p + t u| = r
      const double A = uu, B = 2 * dot(p, u), Cc = dot(p, p) - r * r;
      const double disc = B * B - 4 * A * Cc;
      if (disc < 0)
        continue;
      for (int sg = -1; sg <= 1; sg += 2)
      {
        const double t = (-B + sg * std::sqrt(disc)) / (2 * A);
        const V3 x = add(p, mul(t, u));
        if (inside(x, tol))
          bp.push_back(phi_of(x));
      }
    }
  for (const auto& f : faces)
  {
    if (std::abs(f.d) >= r)
      continue;
    const double rho = std::sqrt(r * r - f.d * f.d);
    V3 t0 = std::abs(f.n[0]) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
    const V3 u = unit(cross(f.n, t0)), v = cross(f.n, u);
    const V3 c0 = mul(f.d, f.n);
    // X(t) = a0 + a1 cos t + a2 sin t, Y(t) = b0 + b1 cos t + b2 sin t
    const double a0 = dot(c0, ex), a1 = rho * dot(u, ex), a2 = rho * dot(v, ex);
    const double b0 = dot(c0, ey), b1 = rho * dot(u, ey), b2 = rho * dot(v, ey);
    const double P = a0 * b2 - b0 * a2, Q = -a0 * b1 + b0 * a1, K = a1 * b2 - a2 * b1;
    const double RR = std::hypot(P, Q);
    if (RR < 1e-300 || std::abs(K) > RR)
      continue;
    const double tau = std::atan2(Q, P), g = std::acos(-K / RR);
    for (int sg = -1; sg <= 1; sg += 2)
    {
      const double t = tau + sg * g;
      const V3 x = add(c0, add(mul(rho * std::cos(t), u), mul(rho * std::sin(t), v)));
      if (inside(x, tol))
        bp.push_back(phi_of(x));
    }
  }
  if (bp.empty())
    return 0.0;
  std::sort(bp.begin(), bp.end());
  double area = 0;
  for (size_t k = 0; k + 1 < bp.size(); ++k)
    if (bp[k + 1] > bp[k])
      area += tanh_sinh(fphi, bp[k], bp[k + 1]);
  return r * r * area;
}

// Area of polygon (planar loop) intersected with the disk of radius rho centred at c (same plane).
inline double polygon_disk_area(const std::vector<V3>& loop, const V3& n, const V3& c, double rho)
{
  V3 t0 = std::abs(n[0]) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
  const V3 e1 = unit(cross(n, t0)), e2 = cross(n, e1);
  auto p2 = [&](const V3& x) { return std::array<double, 2>{dot(sub(x, c), e1), dot(sub(x, c), e2)}; };
  auto seg = [&](std::array<double, 2> p, std::array<double, 2> q)
  {
    // signed area of disk n triangle (0, p, q)
    const double dx = q[0] - p[0], dy = q[1] - p[1];
    const double A = dx * dx + dy * dy, B = 2 * (p[0] * dx + p[1] * dy), C = p[0] * p[0] + p[1] * p[1] - rho * rho;
    std::vector<double> ts = {0.0};
    const double disc = B * B - 4 * A * C;
    if (A > 0 && disc > 0)
    {
      const double s = std::sqrt(disc);
      for (double t : {(-B - s) / (2 * A), (-B + s) / (2 * A)})
        if (t > 0 && t < 1)
          ts.push_back(t);
    }
    ts.push_back(1.0);
    double a = 0;
    for (size_t k = 0; k + 1 < ts.size(); ++k)
    {
      const std::array<double, 2> u = {p[0] + ts[k] * dx, p[1] + ts[k] * dy};
      const std::array<double, 2> w = {p[0] + ts[k + 1] * dx, p[1] + ts[k + 1] * dy};
      const double tm = 0.5 * (ts[k] + ts[k + 1]);
      const double mx = p[0] + tm * dx, my = p[1] + tm * dy;
      const double cr = u[0] * w[1] - u[1] * w[0];
      if (mx * mx + my * my <= rho * rho)
        a += 0.5 * cr; // triangle
      else
        a += 0.5 * rho * rho * std::atan2(cr, u[0] * w[0] + u[1] * w[1]); // sector
    }
    return a;
  };
  double a = 0;
  for (size_t i = 0; i < loop.size(); ++i)
    a += seg(p2(loop[i]), p2(loop[(i + 1) % loop.size()]));
  return std::abs(a);
}

inline double ball_volume(const std::vector<Face>& faces, double r, double area)
{
  double s = r * area;
  for (const auto& f : faces)
  {
    if (std::abs(f.d) >= r)
      continue;
    s += f.d * polygon_disk_area(f.loop, f.n, mul(f.d, f.n), std::sqrt(r * r - f.d * f.d));
  }
  return s / 3.0;
}

// Tetrahedron helper
inline std::vector<Face> tet_faces(const std::array<V3, 4>& X)
{
  return make_faces({{X[1], X[2], X[3]}, {X[0], X[2], X[3]}, {X[0], X[1], X[3]}, {X[0], X[1], X[2]}});
}
inline std::vector<Face> box_faces(const V3& lo, double h)
{
  auto P = [&](int i, int j, int k) { return V3{lo[0] + i * h, lo[1] + j * h, lo[2] + k * h}; };
  return make_faces({{P(0, 0, 0), P(0, 1, 0), P(0, 1, 1), P(0, 0, 1)},
                     {P(1, 0, 0), P(1, 1, 0), P(1, 1, 1), P(1, 0, 1)},
                     {P(0, 0, 0), P(1, 0, 0), P(1, 0, 1), P(0, 0, 1)},
                     {P(0, 1, 0), P(1, 1, 0), P(1, 1, 1), P(0, 1, 1)},
                     {P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(0, 1, 0)},
                     {P(0, 0, 1), P(1, 0, 1), P(1, 1, 1), P(0, 1, 1)}});
}
} // namespace cutcells::proto::exact
