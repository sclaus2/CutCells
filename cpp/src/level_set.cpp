// Copyright (c) 2026 ONERA
// Authors: Susanne Claus
// This file is part of CutCells
// SPDX-License-Identifier:    MIT

#include "level_set.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "cell_topology.h"
#include "cell_types.h"
#include "reference_cell.h"

namespace cutcells
{
namespace
{

template <std::integral I>
struct EdgeKey
{
  I v0 = -1;
  I v1 = -1;

  bool operator==(const EdgeKey&) const = default;
};

template <std::integral I>
struct EdgeKeyHash
{
  std::size_t operator()(const EdgeKey<I>& key) const noexcept
  {
    std::size_t seed = 0;
    seed ^= std::hash<I>{}(key.v0) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= std::hash<I>{}(key.v1) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

template <std::integral I>
struct FaceKey
{
  int nverts = 0;
  std::array<I, 4> verts = {-1, -1, -1, -1};

  bool operator==(const FaceKey&) const = default;
};

template <std::integral I>
struct FaceKeyHash
{
  std::size_t operator()(const FaceKey<I>& key) const noexcept
  {
    std::size_t seed = std::hash<int>{}(key.nverts);
    for (int i = 0; i < key.nverts; ++i)
      seed ^= std::hash<I>{}(key.verts[static_cast<std::size_t>(i)])
           + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct Param2Key
{
  int a = 0;
  int b = 0;

  bool operator==(const Param2Key&) const = default;
};

struct Param2KeyHash
{
  std::size_t operator()(const Param2Key& key) const noexcept
  {
    std::size_t seed = 0;
    seed ^= std::hash<int>{}(key.a) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= std::hash<int>{}(key.b) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

template <std::integral I>
struct FaceData
{
  cell::type face_type = cell::type::point;
  int nverts = 0;
  std::array<I, 4> verts = {-1, -1, -1, -1};
  std::unordered_map<Param2Key, I, Param2KeyHash> dof_map;
};

template <std::floating_point T, std::integral I>
I append_dof(LevelSetMeshData<T, I>& out, std::span<const T> x,
             int8_t parent_dim, int32_t parent_id, std::span<const T> params)
{
  const I dof = out.num_dofs();
  out.dof_coordinates.insert(out.dof_coordinates.end(), x.begin(), x.end());
  out.dof_parent_dim.push_back(parent_dim);
  out.dof_parent_id.push_back(parent_id);
  out.dof_parent_param.insert(out.dof_parent_param.end(), params.begin(), params.end());
  out.dof_parent_param_offset.push_back(
      static_cast<int32_t>(out.dof_parent_param.size()));
  return dof;
}

template <std::floating_point T, std::integral I>
void initialize_vertex_dofs(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out)
{
  const I ndofs = mesh.num_nodes();
  out.dof_coordinates.resize(static_cast<std::size_t>(ndofs)
                             * static_cast<std::size_t>(mesh.gdim));
  for (I i = 0; i < ndofs; ++i)
  {
    const T* x = mesh.node(i);
    for (int d = 0; d < mesh.gdim; ++d)
    {
      out.dof_coordinates[static_cast<std::size_t>(i)
                              * static_cast<std::size_t>(mesh.gdim)
                          + static_cast<std::size_t>(d)] = x[d];
    }
  }
  out.dof_parent_dim.assign(static_cast<std::size_t>(ndofs), 0);
  out.dof_parent_id.resize(static_cast<std::size_t>(ndofs));
  for (I i = 0; i < ndofs; ++i)
    out.dof_parent_id[static_cast<std::size_t>(i)] = static_cast<int32_t>(i);
  out.dof_parent_param_offset.assign(static_cast<std::size_t>(ndofs) + 1, 0);
}

template <std::floating_point T, std::integral I>
cell::type infer_cell_type(const MeshView<T, I>& mesh, I cell_id)
{
  if (mesh.has_cell_types())
  {
    return mesh.cell_type(cell_id);
  }

  const I n = mesh.cell_num_nodes(cell_id);
  switch (mesh.tdim)
  {
    case 1:
      if (n == 2)
        return cell::type::interval;
      break;
    case 2:
      if (n == 3)
        return cell::type::triangle;
      if (n == 4)
        return cell::type::quadrilateral;
      break;
    case 3:
      if (n == 4)
        return cell::type::tetrahedron;
      if (n == 5)
        return cell::type::pyramid;
      if (n == 6)
        return cell::type::prism;
      if (n == 8)
        return cell::type::hexahedron;
      break;
  }

  throw std::runtime_error("create_level_set_mesh_data: unsupported cell type inference from MeshView");
}

template <std::floating_point T, std::integral I>
std::vector<I> cell_vertices_in_basix_order(const MeshView<T, I>& mesh, I cell_id,
                                            cell::type ctype)
{
  std::vector<I> node_scratch;
  const std::span<const I> nodes = mesh.cell_nodes(cell_id, node_scratch);
  const int nverts = cell::get_num_vertices(ctype);
  if (static_cast<int>(nodes.size()) != nverts)
  {
    throw std::runtime_error(
        "create_level_set_mesh_data: MeshView cell connectivity must contain only corner vertices");
  }

  std::vector<I> verts;
  verts.reserve(static_cast<std::size_t>(nverts));
  if (!mesh.vtk_vertex_order)
  {
    verts.insert(verts.end(), nodes.begin(), nodes.end());
  }
  else
  {
    const auto perm = cell::vtk_to_basix_vertex_permutation(ctype);
    for (int i = 0; i < nverts; ++i)
      verts.push_back(nodes[static_cast<std::size_t>(perm[static_cast<std::size_t>(i)])]);
  }
  return verts;
}

template <std::floating_point T, std::integral I>
std::vector<T> gather_vertex_coordinates(const MeshView<T, I>& mesh,
                                         const std::vector<I>& vertex_ids)
{
  std::vector<T> verts;
  verts.reserve(static_cast<std::size_t>(vertex_ids.size() * mesh.gdim));
  for (I v : vertex_ids)
  {
    const T* x = mesh.node(v);
    verts.insert(verts.end(), x, x + mesh.gdim);
  }
  return verts;
}

template <std::integral I>
struct FaceDef
{
  cell::type face_type = cell::type::point;
  int nverts = 0;
  std::array<int, 4> verts = {-1, -1, -1, -1};
};

inline std::span<const FaceDef<int>> basix_faces(cell::type ctype)
{
  static constexpr std::array<FaceDef<int>, 4> tetrahedron = {{
      {cell::type::triangle, 3, {1, 2, 3, -1}},
      {cell::type::triangle, 3, {0, 2, 3, -1}},
      {cell::type::triangle, 3, {0, 1, 3, -1}},
      {cell::type::triangle, 3, {0, 1, 2, -1}},
  }};
  static constexpr std::array<FaceDef<int>, 6> hexahedron = {{
      {cell::type::quadrilateral, 4, {0, 1, 2, 3}},
      {cell::type::quadrilateral, 4, {4, 5, 6, 7}},
      {cell::type::quadrilateral, 4, {0, 1, 4, 5}},
      {cell::type::quadrilateral, 4, {0, 2, 4, 6}},
      {cell::type::quadrilateral, 4, {1, 3, 5, 7}},
      {cell::type::quadrilateral, 4, {2, 3, 6, 7}},
  }};
  static constexpr std::array<FaceDef<int>, 5> prism = {{
      {cell::type::triangle, 3, {0, 1, 2, -1}},
      {cell::type::triangle, 3, {3, 4, 5, -1}},
      {cell::type::quadrilateral, 4, {0, 1, 3, 4}},
      {cell::type::quadrilateral, 4, {0, 2, 3, 5}},
      {cell::type::quadrilateral, 4, {1, 2, 4, 5}},
  }};
  static constexpr std::array<FaceDef<int>, 5> pyramid = {{
      {cell::type::quadrilateral, 4, {0, 1, 2, 3}},
      {cell::type::triangle, 3, {0, 1, 4, -1}},
      {cell::type::triangle, 3, {1, 2, 4, -1}},
      {cell::type::triangle, 3, {2, 3, 4, -1}},
      {cell::type::triangle, 3, {0, 3, 4, -1}},
  }};

  switch (ctype)
  {
    case cell::type::tetrahedron: return std::span(tetrahedron);
    case cell::type::hexahedron: return std::span(hexahedron);
    case cell::type::prism: return std::span(prism);
    case cell::type::pyramid: return std::span(pyramid);
    default: return {};
  }
}

template <typename I>
std::array<I, 3> canonical_triangle(const std::array<I, 3>& v)
{
  std::array<std::array<I, 3>, 6> perms = {{
      {{v[0], v[1], v[2]}},
      {{v[0], v[2], v[1]}},
      {{v[1], v[0], v[2]}},
      {{v[1], v[2], v[0]}},
      {{v[2], v[0], v[1]}},
      {{v[2], v[1], v[0]}}
  }};
  return *std::min_element(perms.begin(), perms.end());
}

template <typename I>
std::array<I, 4> canonical_quad(const std::array<I, 4>& v)
{
  std::array<std::array<I, 4>, 8> perms = {{
      {{v[0], v[1], v[2], v[3]}},
      {{v[1], v[3], v[0], v[2]}},
      {{v[3], v[2], v[1], v[0]}},
      {{v[2], v[0], v[3], v[1]}},
      {{v[0], v[2], v[1], v[3]}},
      {{v[2], v[3], v[0], v[1]}},
      {{v[3], v[1], v[2], v[0]}},
      {{v[1], v[0], v[3], v[2]}}
  }};
  return *std::min_element(perms.begin(), perms.end());
}

template <typename I>
FaceKey<I> canonical_face_key(cell::type face_type, const std::vector<I>& local_vertices)
{
  FaceKey<I> key;
  key.nverts = static_cast<int>(local_vertices.size());
  if (face_type == cell::type::triangle)
  {
    const std::array<I, 3> loc = {
        local_vertices[0], local_vertices[1], local_vertices[2]};
    const std::array<I, 3> can = canonical_triangle(loc);
    std::copy(can.begin(), can.end(), key.verts.begin());
  }
  else if (face_type == cell::type::quadrilateral)
  {
    const std::array<I, 4> loc = {
        local_vertices[0], local_vertices[1], local_vertices[2], local_vertices[3]};
    const std::array<I, 4> can = canonical_quad(loc);
    std::copy(can.begin(), can.end(), key.verts.begin());
  }
  else
  {
    throw std::runtime_error("create_level_set_mesh_data: unsupported face type");
  }
  return key;
}

template <typename I>
std::vector<int> local_to_canonical_positions(const std::vector<I>& local_vertices,
                                              const std::array<I, 4>& canonical_vertices,
                                              int nverts)
{
  std::vector<int> pos(static_cast<std::size_t>(nverts), -1);
  for (int p = 0; p < nverts; ++p)
  {
    for (int q = 0; q < nverts; ++q)
    {
      if (local_vertices[static_cast<std::size_t>(p)] == canonical_vertices[static_cast<std::size_t>(q)])
      {
        pos[static_cast<std::size_t>(p)] = q;
        break;
      }
    }
    if (pos[static_cast<std::size_t>(p)] < 0)
      throw std::runtime_error("create_level_set_mesh_data: invalid face permutation");
  }
  return pos;
}

template <std::floating_point T>
void affine_edge_point(std::span<const T> x0, std::span<const T> x1, T t, std::vector<T>& out)
{
  out.resize(x0.size());
  for (std::size_t d = 0; d < x0.size(); ++d)
    out[d] = (T(1) - t) * x0[d] + t * x1[d];
}

template <std::floating_point T>
void bary_triangle_point(std::span<const T> x0, std::span<const T> x1, std::span<const T> x2,
                         T u, T v, std::vector<T>& out)
{
  const T w0 = T(1) - u - v;
  out.resize(x0.size());
  for (std::size_t d = 0; d < x0.size(); ++d)
    out[d] = w0 * x0[d] + u * x1[d] + v * x2[d];
}

template <std::floating_point T>
void bilinear_quad_point(std::span<const T> x0, std::span<const T> x1,
                         std::span<const T> x2, std::span<const T> x3,
                         T u, T v, std::vector<T>& out)
{
  const T w0 = (T(1) - u) * (T(1) - v);
  const T w1 = u * (T(1) - v);
  const T w2 = (T(1) - u) * v;
  const T w3 = u * v;
  out.resize(x0.size());
  for (std::size_t d = 0; d < x0.size(); ++d)
    out[d] = w0 * x0[d] + w1 * x1[d] + w2 * x2[d] + w3 * x3[d];
}

template <std::floating_point T>
void bary_tetra_point(std::span<const T> x0, std::span<const T> x1,
                      std::span<const T> x2, std::span<const T> x3,
                      T u, T v, T w, std::vector<T>& out)
{
  const T w0 = T(1) - u - v - w;
  out.resize(x0.size());
  for (std::size_t d = 0; d < x0.size(); ++d)
    out[d] = w0 * x0[d] + u * x1[d] + v * x2[d] + w * x3[d];
}

template <std::floating_point T>
void multilinear_hex_point(const std::vector<T>& verts, int gdim,
                           T u, T v, T w, std::vector<T>& out)
{
  out.resize(static_cast<std::size_t>(gdim));
  for (int d = 0; d < gdim; ++d)
  {
    out[static_cast<std::size_t>(d)] =
        (T(1) - u) * (T(1) - v) * (T(1) - w) * verts[static_cast<std::size_t>(d)] +
        u * (T(1) - v) * (T(1) - w) * verts[static_cast<std::size_t>(gdim + d)] +
        (T(1) - u) * v * (T(1) - w) * verts[static_cast<std::size_t>(2 * gdim + d)] +
        u * v * (T(1) - w) * verts[static_cast<std::size_t>(3 * gdim + d)] +
        (T(1) - u) * (T(1) - v) * w * verts[static_cast<std::size_t>(4 * gdim + d)] +
        u * (T(1) - v) * w * verts[static_cast<std::size_t>(5 * gdim + d)] +
        (T(1) - u) * v * w * verts[static_cast<std::size_t>(6 * gdim + d)] +
        u * v * w * verts[static_cast<std::size_t>(7 * gdim + d)];
  }
}

template <std::floating_point T>
void multilinear_prism_point(const std::vector<T>& verts, int gdim,
                             T u, T v, T w, std::vector<T>& out)
{
  const T l0 = T(1) - u - v;
  const T l1 = u;
  const T l2 = v;
  out.resize(static_cast<std::size_t>(gdim));
  for (int d = 0; d < gdim; ++d)
  {
    out[static_cast<std::size_t>(d)] =
        l0 * (T(1) - w) * verts[static_cast<std::size_t>(d)] +
        l1 * (T(1) - w) * verts[static_cast<std::size_t>(gdim + d)] +
        l2 * (T(1) - w) * verts[static_cast<std::size_t>(2 * gdim + d)] +
        l0 * w * verts[static_cast<std::size_t>(3 * gdim + d)] +
        l1 * w * verts[static_cast<std::size_t>(4 * gdim + d)] +
        l2 * w * verts[static_cast<std::size_t>(5 * gdim + d)];
  }
}

template <std::integral I>
EdgeKey<I> make_edge_key(I a, I b)
{
  return (a <= b) ? EdgeKey<I>{a, b} : EdgeKey<I>{b, a};
}

template <std::integral I>
bool edge_matches_canonical(I a, I b, const EdgeKey<I>& key)
{
  return a == key.v0 && b == key.v1;
}

inline std::array<int, 4> transform_quad_indices(int i, int j, int p,
                                                 const std::vector<int>& local_to_can)
{
  static constexpr std::array<std::array<int, 4>, 8> perms = {{
      {{0, 1, 2, 3}},
      {{1, 3, 0, 2}},
      {{3, 2, 1, 0}},
      {{2, 0, 3, 1}},
      {{0, 2, 1, 3}},
      {{2, 3, 0, 1}},
      {{3, 1, 2, 0}},
      {{1, 0, 3, 2}},
  }};

  for (const auto& perm : perms)
  {
    bool match = true;
    for (int k = 0; k < 4; ++k)
    {
      if (local_to_can[static_cast<std::size_t>(k)] != perm[static_cast<std::size_t>(k)])
      {
        match = false;
        break;
      }
    }
    if (!match)
      continue;

    switch (&perm - &perms[0])
    {
      case 0: return {i, j, 0, 0};
      case 1: return {j, p - i, 0, 0};
      case 2: return {p - i, p - j, 0, 0};
      case 3: return {p - j, i, 0, 0};
      case 4: return {j, i, 0, 0};
      case 5: return {p - i, j, 0, 0};
      case 6: return {p - j, p - i, 0, 0};
      case 7: return {i, p - j, 0, 0};
    }
  }

  throw std::runtime_error("create_level_set_mesh_data: unsupported quad face permutation");
}

template <std::floating_point T, std::integral I>
I get_or_create_edge_block(const MeshView<T, I>& mesh,
                           LevelSetMeshData<T, I>& out,
                           std::unordered_map<EdgeKey<I>, I, EdgeKeyHash<I>>& edge_ids,
                           std::vector<I>& edge_first_dof,
                           I a, I b, int degree)
{
  const EdgeKey<I> key = make_edge_key(a, b);
  auto it = edge_ids.find(key);
  if (it != edge_ids.end())
    return it->second;

  const I edge_id = static_cast<I>(edge_ids.size());
  edge_ids.emplace(key, edge_id);
  edge_first_dof.push_back(static_cast<I>(-1));
  if (degree == 1)
    return edge_id;

  edge_first_dof[static_cast<std::size_t>(edge_id)] = out.num_dofs();
  std::vector<T> x;
  std::array<T, 1> param;
  for (int k = 1; k < degree; ++k)
  {
    const T t = static_cast<T>(k) / static_cast<T>(degree);
    affine_edge_point(
        std::span<const T>(mesh.node(key.v0), static_cast<std::size_t>(mesh.gdim)),
        std::span<const T>(mesh.node(key.v1), static_cast<std::size_t>(mesh.gdim)),
        t, x);
    param[0] = t;
    append_dof(out, std::span<const T>(x.data(), x.size()), 1,
               static_cast<int32_t>(edge_id),
               std::span<const T>(param.data(), param.size()));
  }
  return edge_id;
}

template <std::floating_point T, std::integral I>
I get_or_create_face_block(const MeshView<T, I>& mesh,
                           LevelSetMeshData<T, I>& out,
                           std::unordered_map<FaceKey<I>, I, FaceKeyHash<I>>& face_ids,
                           std::vector<FaceData<I>>& faces,
                           cell::type face_type,
                           const std::vector<I>& local_vertices,
                           int degree)
{
  const FaceKey<I> key = canonical_face_key(face_type, local_vertices);
  auto it = face_ids.find(key);
  if (it != face_ids.end())
    return it->second;

  const I face_id = static_cast<I>(faces.size());
  face_ids.emplace(key, face_id);

  FaceData<I> face;
  face.face_type = face_type;
  face.nverts = key.nverts;
  face.verts = key.verts;

  if (degree > 2 || face_type == cell::type::quadrilateral)
  {
    std::vector<T> x;
    if (face_type == cell::type::triangle)
    {
      std::array<T, 2> param;
      const T* x0 = mesh.node(face.verts[0]);
      const T* x1 = mesh.node(face.verts[1]);
      const T* x2 = mesh.node(face.verts[2]);
      for (int j = 1; j <= degree - 2; ++j)
      {
        for (int i = 1; i <= degree - j - 1; ++i)
        {
          const T u = static_cast<T>(i) / static_cast<T>(degree);
          const T v = static_cast<T>(j) / static_cast<T>(degree);
          bary_triangle_point(
              std::span<const T>(x0, static_cast<std::size_t>(mesh.gdim)),
              std::span<const T>(x1, static_cast<std::size_t>(mesh.gdim)),
              std::span<const T>(x2, static_cast<std::size_t>(mesh.gdim)),
              u, v, x);
          param[0] = u;
          param[1] = v;
          const I dof = append_dof(out, std::span<const T>(x.data(), x.size()), 2,
                                   static_cast<int32_t>(face_id),
                                   std::span<const T>(param.data(), param.size()));
          face.dof_map.emplace(Param2Key{i, j}, dof);
        }
      }
    }
    else
    {
      std::array<T, 2> param;
      const T* x0 = mesh.node(face.verts[0]);
      const T* x1 = mesh.node(face.verts[1]);
      const T* x2 = mesh.node(face.verts[2]);
      const T* x3 = mesh.node(face.verts[3]);
      for (int j = 1; j < degree; ++j)
      {
        for (int i = 1; i < degree; ++i)
        {
          const T u = static_cast<T>(i) / static_cast<T>(degree);
          const T v = static_cast<T>(j) / static_cast<T>(degree);
          bilinear_quad_point(
              std::span<const T>(x0, static_cast<std::size_t>(mesh.gdim)),
              std::span<const T>(x1, static_cast<std::size_t>(mesh.gdim)),
              std::span<const T>(x2, static_cast<std::size_t>(mesh.gdim)),
              std::span<const T>(x3, static_cast<std::size_t>(mesh.gdim)),
              u, v, x);
          param[0] = u;
          param[1] = v;
          const I dof = append_dof(out, std::span<const T>(x.data(), x.size()), 2,
                                   static_cast<int32_t>(face_id),
                                   std::span<const T>(param.data(), param.size()));
          face.dof_map.emplace(Param2Key{i, j}, dof);
        }
      }
    }
  }

  faces.push_back(std::move(face));
  return face_id;
}

template <std::floating_point T, std::integral I>
void append_interval_cell_dofs(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                               I cell_id, const std::vector<I>& verts, int degree)
{
  out.cell_dofs.push_back(verts[0]);
  out.cell_dofs.push_back(verts[1]);
  if (degree == 1)
    return;

  std::vector<T> x;
  std::array<T, 1> param;
  for (int k = 1; k < degree; ++k)
  {
    const T t = static_cast<T>(k) / static_cast<T>(degree);
    affine_edge_point(
        std::span<const T>(mesh.node(verts[0]), static_cast<std::size_t>(mesh.gdim)),
        std::span<const T>(mesh.node(verts[1]), static_cast<std::size_t>(mesh.gdim)),
        t, x);
    param[0] = t;
    out.cell_dofs.push_back(append_dof(
        out, std::span<const T>(x.data(), x.size()), 1,
        static_cast<int32_t>(cell_id),
        std::span<const T>(param.data(), param.size())));
  }
}

template <std::floating_point T, std::integral I>
void append_triangle_cell_interior(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                                   I cell_id, const std::vector<I>& verts, int degree)
{
  if (degree <= 2)
    return;

  std::vector<T> x;
  std::array<T, 2> param;
  for (int j = 1; j <= degree - 2; ++j)
  {
    for (int i = 1; i <= degree - j - 1; ++i)
    {
      const T u = static_cast<T>(i) / static_cast<T>(degree);
      const T v = static_cast<T>(j) / static_cast<T>(degree);
      bary_triangle_point(
          std::span<const T>(mesh.node(verts[0]), static_cast<std::size_t>(mesh.gdim)),
          std::span<const T>(mesh.node(verts[1]), static_cast<std::size_t>(mesh.gdim)),
          std::span<const T>(mesh.node(verts[2]), static_cast<std::size_t>(mesh.gdim)),
          u, v, x);
      param[0] = u;
      param[1] = v;
      out.cell_dofs.push_back(append_dof(
          out, std::span<const T>(x.data(), x.size()), 2,
          static_cast<int32_t>(cell_id),
          std::span<const T>(param.data(), param.size())));
    }
  }
}

template <std::floating_point T, std::integral I>
void append_quadrilateral_cell_interior(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                                        I cell_id, const std::vector<I>& verts, int degree)
{
  if (degree <= 1)
    return;

  std::vector<T> x;
  std::array<T, 2> param;
  for (int j = 1; j < degree; ++j)
  {
    for (int i = 1; i < degree; ++i)
    {
      const T u = static_cast<T>(i) / static_cast<T>(degree);
      const T v = static_cast<T>(j) / static_cast<T>(degree);
      bilinear_quad_point(
          std::span<const T>(mesh.node(verts[0]), static_cast<std::size_t>(mesh.gdim)),
          std::span<const T>(mesh.node(verts[1]), static_cast<std::size_t>(mesh.gdim)),
          std::span<const T>(mesh.node(verts[2]), static_cast<std::size_t>(mesh.gdim)),
          std::span<const T>(mesh.node(verts[3]), static_cast<std::size_t>(mesh.gdim)),
          u, v, x);
      param[0] = u;
      param[1] = v;
      out.cell_dofs.push_back(append_dof(
          out, std::span<const T>(x.data(), x.size()), 2,
          static_cast<int32_t>(cell_id),
          std::span<const T>(param.data(), param.size())));
    }
  }
}

template <std::floating_point T, std::integral I>
void append_tetrahedron_cell_interior(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                                      I cell_id, const std::vector<I>& verts, int degree)
{
  if (degree <= 3)
    return;

  std::vector<T> x;
  std::array<T, 3> param;
  for (int k = 1; k <= degree - 3; ++k)
  {
    for (int j = 1; j <= degree - k - 2; ++j)
    {
      for (int i = 1; i <= degree - j - k - 1; ++i)
      {
        const T u = static_cast<T>(i) / static_cast<T>(degree);
        const T v = static_cast<T>(j) / static_cast<T>(degree);
        const T w = static_cast<T>(k) / static_cast<T>(degree);
        bary_tetra_point(
            std::span<const T>(mesh.node(verts[0]), static_cast<std::size_t>(mesh.gdim)),
            std::span<const T>(mesh.node(verts[1]), static_cast<std::size_t>(mesh.gdim)),
            std::span<const T>(mesh.node(verts[2]), static_cast<std::size_t>(mesh.gdim)),
            std::span<const T>(mesh.node(verts[3]), static_cast<std::size_t>(mesh.gdim)),
            u, v, w, x);
        param[0] = u;
        param[1] = v;
        param[2] = w;
        out.cell_dofs.push_back(append_dof(
            out, std::span<const T>(x.data(), x.size()), 3,
            static_cast<int32_t>(cell_id),
            std::span<const T>(param.data(), param.size())));
      }
    }
  }
}

template <std::floating_point T, std::integral I>
void append_hexahedron_cell_interior(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                                     I cell_id, const std::vector<I>& verts, int degree)
{
  if (degree <= 1)
    return;

  const std::vector<T> xverts = gather_vertex_coordinates(mesh, verts);
  std::vector<T> x;
  std::array<T, 3> param;
  for (int k = 1; k < degree; ++k)
  {
    for (int j = 1; j < degree; ++j)
    {
      for (int i = 1; i < degree; ++i)
      {
        const T u = static_cast<T>(i) / static_cast<T>(degree);
        const T v = static_cast<T>(j) / static_cast<T>(degree);
        const T w = static_cast<T>(k) / static_cast<T>(degree);
        multilinear_hex_point(xverts, mesh.gdim, u, v, w, x);
        param[0] = u;
        param[1] = v;
        param[2] = w;
        out.cell_dofs.push_back(append_dof(
            out, std::span<const T>(x.data(), x.size()), 3,
            static_cast<int32_t>(cell_id),
            std::span<const T>(param.data(), param.size())));
      }
    }
  }
}

template <std::floating_point T, std::integral I>
void append_prism_cell_interior(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                                I cell_id, const std::vector<I>& verts, int degree)
{
  if (degree <= 2)
    return;

  const std::vector<T> xverts = gather_vertex_coordinates(mesh, verts);
  std::vector<T> x;
  std::array<T, 3> param;
  for (int k = 1; k < degree; ++k)
  {
    for (int j = 1; j <= degree - 2; ++j)
    {
      for (int i = 1; i <= degree - j - 1; ++i)
      {
        const T u = static_cast<T>(i) / static_cast<T>(degree);
        const T v = static_cast<T>(j) / static_cast<T>(degree);
        const T w = static_cast<T>(k) / static_cast<T>(degree);
        multilinear_prism_point(xverts, mesh.gdim, u, v, w, x);
        param[0] = u;
        param[1] = v;
        param[2] = w;
        out.cell_dofs.push_back(append_dof(
            out, std::span<const T>(x.data(), x.size()), 3,
            static_cast<int32_t>(cell_id),
            std::span<const T>(param.data(), param.size())));
      }
    }
  }
}

template <std::floating_point T, std::integral I>
void append_cell_dofs(const MeshView<T, I>& mesh, LevelSetMeshData<T, I>& out,
                      I cell_id, cell::type ctype, const std::vector<I>& verts,
                      int degree,
                      std::unordered_map<EdgeKey<I>, I, EdgeKeyHash<I>>& edge_ids,
                      std::vector<I>& edge_first_dof,
                      std::unordered_map<FaceKey<I>, I, FaceKeyHash<I>>& face_ids,
                      std::vector<FaceData<I>>& faces)
{
  for (I v : verts)
    out.cell_dofs.push_back(v);

  if (ctype == cell::type::interval)
  {
    append_interval_cell_dofs(mesh, out, cell_id, verts, degree);
    return;
  }

  const auto edges = cell::edges(ctype);
  for (std::size_t e = 0; e < edges.size(); ++e)
  {
    const I a = verts[static_cast<std::size_t>(edges[e][0])];
    const I b = verts[static_cast<std::size_t>(edges[e][1])];
    const EdgeKey<I> key = make_edge_key(a, b);
    const bool same_orientation = edge_matches_canonical(a, b, key);
    const I edge_id = get_or_create_edge_block(
        mesh, out, edge_ids, edge_first_dof, a, b, degree);
    if (degree == 1)
      continue;

    const I first = edge_first_dof[static_cast<std::size_t>(edge_id)];
    for (int k = 1; k < degree; ++k)
    {
      const int canonical_k = same_orientation ? k : (degree - k);
      out.cell_dofs.push_back(first + static_cast<I>(canonical_k - 1));
    }
  }

  if (mesh.tdim == 3)
  {
    const auto faces_local = basix_faces(ctype);
    for (std::size_t f = 0; f < faces_local.size(); ++f)
    {
      const auto& fdef = faces_local[f];
      std::vector<I> fverts;
      fverts.reserve(static_cast<std::size_t>(fdef.nverts));
      for (int j = 0; j < fdef.nverts; ++j)
        fverts.push_back(verts[static_cast<std::size_t>(fdef.verts[static_cast<std::size_t>(j)])]);

      const I face_id = get_or_create_face_block(
          mesh, out, face_ids, faces, fdef.face_type, fverts, degree);
      const FaceData<I>& face = faces[static_cast<std::size_t>(face_id)];
      const std::vector<int> local_to_can =
          local_to_canonical_positions(fverts, face.verts, fdef.nverts);

      if (fdef.face_type == cell::type::triangle)
      {
        for (int j = 1; j <= degree - 2; ++j)
        {
          for (int i = 1; i <= degree - j - 1; ++i)
          {
            int local_bary[3] = {degree - i - j, i, j};
            int can_bary[3] = {0, 0, 0};
            for (int p = 0; p < 3; ++p)
              can_bary[local_to_can[static_cast<std::size_t>(p)]] = local_bary[p];

            const auto it = face.dof_map.find(Param2Key{can_bary[1], can_bary[2]});
            if (it == face.dof_map.end())
              throw std::runtime_error("create_level_set_mesh_data: missing triangle face dof");
            out.cell_dofs.push_back(it->second);
          }
        }
      }
      else
      {
        for (int j = 1; j < degree; ++j)
        {
          for (int i = 1; i < degree; ++i)
          {
            const auto can = transform_quad_indices(i, j, degree, local_to_can);
            const auto it = face.dof_map.find(Param2Key{can[0], can[1]});
            if (it == face.dof_map.end())
              throw std::runtime_error("create_level_set_mesh_data: missing quadrilateral face dof");
            out.cell_dofs.push_back(it->second);
          }
        }
      }
    }
  }

  switch (ctype)
  {
    case cell::type::triangle:
      append_triangle_cell_interior(mesh, out, cell_id, verts, degree);
      break;
    case cell::type::quadrilateral:
      append_quadrilateral_cell_interior(mesh, out, cell_id, verts, degree);
      break;
    case cell::type::tetrahedron:
      append_tetrahedron_cell_interior(mesh, out, cell_id, verts, degree);
      break;
    case cell::type::hexahedron:
      append_hexahedron_cell_interior(mesh, out, cell_id, verts, degree);
      break;
    case cell::type::prism:
      append_prism_cell_interior(mesh, out, cell_id, verts, degree);
      break;
    case cell::type::pyramid:
      if (degree > 1)
      {
        throw std::runtime_error(
            "create_level_set_mesh_data: pyramid degree > 1 is not implemented yet");
      }
      break;
    default:
      break;
  }
}

template <std::floating_point T, std::integral I>
void validate_mesh_data_args(int gdim, int tdim, int degree,
                             std::span<const T> dof_coordinates,
                             std::span<const I> cell_dofs,
                             std::span<const I> cell_offsets,
                             std::span<const cell::type> cell_types)
{
  if (gdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data: gdim must be positive");
  if (tdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data: tdim must be positive");
  if (degree < 1)
    throw std::runtime_error("create_level_set_mesh_data: degree must be >= 1");
  if (dof_coordinates.size() % static_cast<std::size_t>(gdim) != 0)
    throw std::runtime_error("create_level_set_mesh_data: dof_coordinates size must be divisible by gdim");
  if (cell_offsets.empty())
    throw std::runtime_error("create_level_set_mesh_data: cell_offsets must not be empty");
  if (cell_offsets.front() != 0)
    throw std::runtime_error("create_level_set_mesh_data: cell_offsets must start at 0");
  if (static_cast<std::size_t>(cell_offsets.back()) != cell_dofs.size())
    throw std::runtime_error("create_level_set_mesh_data: cell_offsets.back() must equal cell_dofs.size()");
  if (!cell_types.empty() && cell_types.size() + 1 != cell_offsets.size())
    throw std::runtime_error("create_level_set_mesh_data: cell_types size must equal number of cells");
}

template <std::floating_point T, std::integral I>
void validate_mesh_data_view_args(int gdim, int tdim, int degree,
                                  std::span<const T> dof_coordinates,
                                  std::span<const I> cell_dofs,
                                  I num_cells,
                                  I dofs_per_cell,
                                  I dof_stride,
                                  std::span<const cell::type> cell_types)
{
  if (gdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data_view: gdim must be positive");
  if (tdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data_view: tdim must be positive");
  if (degree < 1)
    throw std::runtime_error("create_level_set_mesh_data_view: degree must be >= 1");
  if (dof_coordinates.size() % static_cast<std::size_t>(gdim) != 0)
    throw std::runtime_error("create_level_set_mesh_data_view: dof_coordinates size must be divisible by gdim");
  if (num_cells < 0)
    throw std::runtime_error("create_level_set_mesh_data_view: num_cells must be non-negative");
  if (dofs_per_cell <= 0)
    throw std::runtime_error("create_level_set_mesh_data_view: dofs_per_cell must be positive");
  if (dof_stride < dofs_per_cell)
    throw std::runtime_error("create_level_set_mesh_data_view: dof_stride must be at least dofs_per_cell");
  if (cell_dofs.size()
      < static_cast<std::size_t>(num_cells) * static_cast<std::size_t>(dof_stride))
  {
    throw std::runtime_error(
        "create_level_set_mesh_data_view: cell_dofs view is too small for the requested layout");
  }
  if (!cell_types.empty()
      && cell_types.size() != static_cast<std::size_t>(num_cells))
  {
    throw std::runtime_error(
        "create_level_set_mesh_data_view: cell_types size must equal number of cells");
  }
}

template <std::floating_point T, std::integral I>
void validate_mesh_data_block_view_args(
    int gdim, int tdim, int degree,
    std::span<const T> dof_coordinates,
    std::span<const StridedRowBlock<I>> cell_dof_blocks,
    I num_cells,
    std::span<const cell::type> cell_types)
{
  if (gdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data_view: gdim must be positive");
  if (tdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data_view: tdim must be positive");
  if (degree < 1)
    throw std::runtime_error("create_level_set_mesh_data_view: degree must be >= 1");
  if (dof_coordinates.size() % static_cast<std::size_t>(gdim) != 0)
    throw std::runtime_error("create_level_set_mesh_data_view: dof_coordinates size must be divisible by gdim");
  if (num_cells < 0)
    throw std::runtime_error("create_level_set_mesh_data_view: num_cells must be non-negative");
  if (cell_dof_blocks.empty() && num_cells > 0)
    throw std::runtime_error("create_level_set_mesh_data_view: cell_dof_blocks must not be empty");

  std::vector<std::uint8_t> covered(static_cast<std::size_t>(num_cells), 0);
  for (const auto& block : cell_dof_blocks)
  {
    if (block.first_row < 0)
      throw std::runtime_error("create_level_set_mesh_data_view: block first_row must be non-negative");
    if (block.row_count < 0)
      throw std::runtime_error("create_level_set_mesh_data_view: block row_count must be non-negative");
    if (block.row_width <= 0)
      throw std::runtime_error("create_level_set_mesh_data_view: block row_width must be positive");
    const I stride = block.row_stride > 0 ? block.row_stride : block.row_width;
    if (stride < block.row_width)
      throw std::runtime_error("create_level_set_mesh_data_view: block row_stride must be at least row_width");
    const std::size_t required = static_cast<std::size_t>(block.row_count)
                               * static_cast<std::size_t>(stride);
    if (block.values.size() < required)
      throw std::runtime_error("create_level_set_mesh_data_view: block values view is too small");
    if (block.first_row + block.row_count > num_cells)
      throw std::runtime_error("create_level_set_mesh_data_view: block exceeds num_cells");
    for (I i = 0; i < block.row_count; ++i)
    {
      const std::size_t row = static_cast<std::size_t>(block.first_row + i);
      if (covered[row] != 0)
        throw std::runtime_error("create_level_set_mesh_data_view: overlapping blocks");
      covered[row] = 1;
    }
  }
  for (std::uint8_t row_is_covered : covered)
  {
    if (row_is_covered == 0)
      throw std::runtime_error("create_level_set_mesh_data_view: blocks must cover all cells");
  }

  if (!cell_types.empty()
      && cell_types.size() != static_cast<std::size_t>(num_cells))
  {
    throw std::runtime_error(
        "create_level_set_mesh_data_view: cell_types size must equal number of cells");
  }
}

template <std::integral I>
void validate_cell_type_blocks(std::span<const CellTypeBlock<I>> cell_type_blocks,
                               I num_cells)
{
  if (cell_type_blocks.empty() && num_cells > 0)
    throw std::runtime_error("create_level_set_mesh_data_view: cell_type_blocks must not be empty");

  std::vector<std::uint8_t> covered(static_cast<std::size_t>(num_cells), 0);
  for (const auto& block : cell_type_blocks)
  {
    if (block.first_cell < 0)
      throw std::runtime_error("create_level_set_mesh_data_view: cell type block first_cell must be non-negative");
    if (block.cell_count < 0)
      throw std::runtime_error("create_level_set_mesh_data_view: cell type block cell_count must be non-negative");
    if (block.first_cell + block.cell_count > num_cells)
      throw std::runtime_error("create_level_set_mesh_data_view: cell type block exceeds num_cells");
    for (I i = 0; i < block.cell_count; ++i)
    {
      const std::size_t row = static_cast<std::size_t>(block.first_cell + i);
      if (covered[row] != 0)
        throw std::runtime_error("create_level_set_mesh_data_view: overlapping cell type blocks");
      covered[row] = 1;
    }
  }
  for (std::uint8_t row_is_covered : covered)
  {
    if (row_is_covered == 0)
      throw std::runtime_error("create_level_set_mesh_data_view: cell type blocks must cover all cells");
  }
}

} // namespace

template <std::floating_point T, std::integral I>
LevelSetMeshData<T, I> create_level_set_mesh_data(
    const MeshView<T, I>& mesh, int degree, T merge_tol)
{
  (void)merge_tol;
  if (degree < 1)
    throw std::runtime_error("create_level_set_mesh_data: degree must be >= 1");
  if (mesh.gdim <= 0 || mesh.tdim <= 0)
    throw std::runtime_error("create_level_set_mesh_data: MeshView dimensions must be positive");

  LevelSetMeshData<T, I> out;
  out.gdim = mesh.gdim;
  out.tdim = mesh.tdim;
  out.degree = degree;
  out.cell_offsets.reserve(static_cast<std::size_t>(mesh.num_cells()) + 1);
  out.cell_types.reserve(static_cast<std::size_t>(mesh.num_cells()));
  out.cell_offsets.push_back(0);

  initialize_vertex_dofs(mesh, out);

  std::unordered_map<EdgeKey<I>, I, EdgeKeyHash<I>> edge_ids;
  std::vector<I> edge_first_dof;
  std::unordered_map<FaceKey<I>, I, FaceKeyHash<I>> face_ids;
  std::vector<FaceData<I>> faces;

  for (I cell_id = 0; cell_id < mesh.num_cells(); ++cell_id)
  {
    const cell::type ctype = infer_cell_type(mesh, cell_id);
    const std::vector<I> verts = cell_vertices_in_basix_order(mesh, cell_id, ctype);
    append_cell_dofs(mesh, out, cell_id, ctype, verts, degree,
                     edge_ids, edge_first_dof, face_ids, faces);

    out.cell_types.push_back(ctype);
    out.cell_offsets.push_back(static_cast<I>(out.cell_dofs.size()));
  }

  return out;
}

template <std::floating_point T, std::integral I>
LevelSetMeshData<T, I> create_level_set_mesh_data(
    int gdim,
    int tdim,
    int degree,
    std::span<const T> dof_coordinates,
    std::span<const I> cell_dofs,
    std::span<const I> cell_offsets,
    std::span<const cell::type> cell_types)
{
  validate_mesh_data_args(gdim, tdim, degree,
                          dof_coordinates, cell_dofs, cell_offsets, cell_types);

  LevelSetMeshData<T, I> out;
  out.gdim = gdim;
  out.tdim = tdim;
  out.degree = degree;
  out.dof_coordinates.assign(dof_coordinates.begin(), dof_coordinates.end());
  out.cell_dofs.assign(cell_dofs.begin(), cell_dofs.end());
  out.cell_offsets.assign(cell_offsets.begin(), cell_offsets.end());
  out.cell_types.assign(cell_types.begin(), cell_types.end());
  return out;
}

template <std::floating_point T, std::integral I>
LevelSetMeshData<T, I> create_level_set_mesh_data_view(
    int gdim,
    int tdim,
    int degree,
    std::vector<T>&& dof_coordinates,
    std::span<const I> cell_dofs,
    I num_cells,
    I dofs_per_cell,
    I dof_stride,
    std::span<const cell::type> cell_types)
{
  validate_mesh_data_view_args(gdim, tdim, degree,
                               std::span<const T>(dof_coordinates.data(),
                                                  dof_coordinates.size()),
                               cell_dofs, num_cells, dofs_per_cell,
                               dof_stride, cell_types);

  LevelSetMeshData<T, I> out;
  out.gdim = gdim;
  out.tdim = tdim;
  out.degree = degree;
  out.dof_coordinates = std::move(dof_coordinates);
  out.cell_dofs_view = cell_dofs;
  out.cell_count = num_cells;
  out.cell_dofs_per_cell = dofs_per_cell;
  out.cell_dof_stride = dof_stride;
  out.cell_types.assign(cell_types.begin(), cell_types.end());
  return out;
}

template <std::floating_point T, std::integral I>
LevelSetMeshData<T, I> create_level_set_mesh_data_view(
    int gdim,
    int tdim,
    int degree,
    std::vector<T>&& dof_coordinates,
    std::span<const StridedRowBlock<I>> cell_dof_blocks,
    I num_cells,
    std::span<const cell::type> cell_types)
{
  validate_mesh_data_block_view_args(gdim, tdim, degree,
                                     std::span<const T>(dof_coordinates.data(),
                                                        dof_coordinates.size()),
                                     cell_dof_blocks, num_cells, cell_types);

  LevelSetMeshData<T, I> out;
  out.gdim = gdim;
  out.tdim = tdim;
  out.degree = degree;
  out.dof_coordinates = std::move(dof_coordinates);
  out.cell_dof_blocks = cell_dof_blocks;
  out.cell_count = num_cells;
  out.cell_types.assign(cell_types.begin(), cell_types.end());
  return out;
}

template <std::floating_point T, std::integral I>
LevelSetMeshData<T, I> create_level_set_mesh_data_view(
    int gdim,
    int tdim,
    int degree,
    std::vector<T>&& dof_coordinates,
    std::span<const StridedRowBlock<I>> cell_dof_blocks,
    I num_cells,
    std::span<const CellTypeBlock<I>> cell_type_blocks)
{
  validate_mesh_data_block_view_args(gdim, tdim, degree,
                                     std::span<const T>(dof_coordinates.data(),
                                                        dof_coordinates.size()),
                                     cell_dof_blocks, num_cells,
                                     std::span<const cell::type>());
  validate_cell_type_blocks(cell_type_blocks, num_cells);

  LevelSetMeshData<T, I> out;
  out.gdim = gdim;
  out.tdim = tdim;
  out.degree = degree;
  out.dof_coordinates = std::move(dof_coordinates);
  out.cell_dof_blocks = cell_dof_blocks;
  out.cell_count = num_cells;
  out.cell_type_blocks = cell_type_blocks;
  return out;
}

template <std::floating_point T, std::integral I>
LevelSetFunction<T, I> create_level_set_function(
    LevelSetMeshData<T, I> mesh_data,
    std::span<const T> dof_values,
    std::string name)
{
  const std::size_t expected = static_cast<std::size_t>(mesh_data.num_dofs());
  if (dof_values.size() != expected)
  {
    throw std::runtime_error(
        "create_level_set_function: dof_values size mismatch (got "
        + std::to_string(dof_values.size()) + ", expected "
        + std::to_string(expected) + ")");
  }

  auto owned_values = std::make_shared<std::vector<T>>(dof_values.begin(), dof_values.end());

  LevelSetFunction<T, I> ls;
  ls.name = std::move(name);
  ls.type = LevelSetType::Polynomial;
  ls.gdim = mesh_data.gdim;
  ls.mesh_data = std::move(mesh_data);
  ls.has_mesh_data_storage = true;
  ls.dof_values = std::span<const T>(owned_values->data(), owned_values->size());
  ls.owner = std::static_pointer_cast<const void>(owned_values);
  return ls;
}

template <std::floating_point T, std::integral I>
LevelSetFunction<T, I> create_level_set_function_view(
    LevelSetMeshData<T, I> mesh_data,
    std::span<const T> dof_values,
    std::string name,
    std::shared_ptr<const void> owner)
{
  const std::size_t expected = static_cast<std::size_t>(mesh_data.num_dofs());
  if (dof_values.size() != expected)
  {
    throw std::runtime_error(
        "create_level_set_function_view: dof_values size mismatch (got "
        + std::to_string(dof_values.size()) + ", expected "
        + std::to_string(expected) + ")");
  }

  LevelSetFunction<T, I> ls;
  ls.name = std::move(name);
  ls.type = LevelSetType::Polynomial;
  ls.gdim = mesh_data.gdim;
  ls.mesh_data = std::move(mesh_data);
  ls.has_mesh_data_storage = true;
  ls.dof_values = dof_values;
  ls.owner = std::move(owner);
  return ls;
}

template LevelSetMeshData<float, int> create_level_set_mesh_data(
    const MeshView<float, int>& mesh, int degree, float merge_tol);
template LevelSetMeshData<double, int> create_level_set_mesh_data(
    const MeshView<double, int>& mesh, int degree, double merge_tol);

template LevelSetMeshData<float, int> create_level_set_mesh_data(
    int gdim, int tdim, int degree,
    std::span<const float> dof_coordinates,
    std::span<const int> cell_dofs,
    std::span<const int> cell_offsets,
    std::span<const cell::type> cell_types);
template LevelSetMeshData<double, int> create_level_set_mesh_data(
    int gdim, int tdim, int degree,
    std::span<const double> dof_coordinates,
    std::span<const int> cell_dofs,
    std::span<const int> cell_offsets,
    std::span<const cell::type> cell_types);

template LevelSetMeshData<float, int> create_level_set_mesh_data_view(
    int gdim, int tdim, int degree,
    std::vector<float>&& dof_coordinates,
    std::span<const int> cell_dofs,
    int num_cells,
    int dofs_per_cell,
    int dof_stride,
    std::span<const cell::type> cell_types);
template LevelSetMeshData<double, int> create_level_set_mesh_data_view(
    int gdim, int tdim, int degree,
    std::vector<double>&& dof_coordinates,
    std::span<const int> cell_dofs,
    int num_cells,
    int dofs_per_cell,
    int dof_stride,
    std::span<const cell::type> cell_types);

template LevelSetMeshData<float, int> create_level_set_mesh_data_view(
    int gdim, int tdim, int degree,
    std::vector<float>&& dof_coordinates,
    std::span<const StridedRowBlock<int>> cell_dof_blocks,
    int num_cells,
    std::span<const cell::type> cell_types);
template LevelSetMeshData<double, int> create_level_set_mesh_data_view(
    int gdim, int tdim, int degree,
    std::vector<double>&& dof_coordinates,
    std::span<const StridedRowBlock<int>> cell_dof_blocks,
    int num_cells,
    std::span<const cell::type> cell_types);

template LevelSetMeshData<float, int> create_level_set_mesh_data_view(
    int gdim, int tdim, int degree,
    std::vector<float>&& dof_coordinates,
    std::span<const StridedRowBlock<int>> cell_dof_blocks,
    int num_cells,
    std::span<const CellTypeBlock<int>> cell_type_blocks);
template LevelSetMeshData<double, int> create_level_set_mesh_data_view(
    int gdim, int tdim, int degree,
    std::vector<double>&& dof_coordinates,
    std::span<const StridedRowBlock<int>> cell_dof_blocks,
    int num_cells,
    std::span<const CellTypeBlock<int>> cell_type_blocks);

template LevelSetFunction<float, int> create_level_set_function(
    LevelSetMeshData<float, int> mesh_data,
    std::span<const float> dof_values,
    std::string name);
template LevelSetFunction<double, int> create_level_set_function(
    LevelSetMeshData<double, int> mesh_data,
    std::span<const double> dof_values,
    std::string name);

template LevelSetFunction<float, int> create_level_set_function_view(
    LevelSetMeshData<float, int> mesh_data,
    std::span<const float> dof_values,
    std::string name,
    std::shared_ptr<const void> owner);
template LevelSetFunction<double, int> create_level_set_function_view(
    LevelSetMeshData<double, int> mesh_data,
    std::span<const double> dof_values,
    std::string name,
    std::shared_ptr<const void> owner);

} // namespace cutcells
