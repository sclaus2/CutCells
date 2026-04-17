// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT


#include <iostream>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/shared_ptr.h>
#include <span>
#include <type_traits>
#include <stdexcept>
#include <memory>
#include <optional>
#include <numeric>
#include <unordered_set>

#include "../../cpp/src/cell_types.h"
#include "../../cpp/src/cut_cell.h"
#include "../../cpp/src/cut_mesh.h"
#include "../../cpp/src/write_vtk.h"
#include "../../cpp/src/bernstein.h"
#include "../../cpp/src/mapping.h"
#include "../../cpp/src/quadrature.h"
#include "../../cpp/src/quadrature_tables.h"
#include "../../cpp/src/mesh_view.h"
#include "../../cpp/src/level_set.h"
#include "../../cpp/src/adapt_cell.h"
#include "../../cpp/src/cell_topology.h"
#include "../../cpp/src/level_set_cell.h"
#include "../../cpp/src/reference_cell.h"
#include "../../cpp/src/triangulation.h"
#include "../../cpp/src/edge_certification.h"
#include "../../cpp/src/cell_certification.h"
#include "../../cpp/src/refine_cell.h"
#include "../../cpp/src/ho_cut_mesh.h"
#include "../../cpp/src/ho_mesh_part_output.h"

namespace nb = nanobind;

using namespace cutcells;

namespace
{
const std::string& cell_domain_to_str(cell::domain domain_id)
{
  static const std::map<cell::domain, std::string> type_to_name
      = {{cell::domain::inside, "inside"},
         {cell::domain::intersected, "intersected"},
         {cell::domain::outside, "outside"}};

  auto it = type_to_name.find(domain_id);
  if (it == type_to_name.end())
    throw std::runtime_error("Can't find type");

  return it->second;
}

template <typename V>
auto as_nbarray(V&& x, std::size_t ndim, const std::size_t* shape)
{
  using _V = std::decay_t<V>;
  _V* ptr = new _V(std::move(x));
  return nb::ndarray<typename _V::value_type, nb::numpy>(
      ptr->data(), ndim, shape,
      nb::capsule(ptr, [](void* p) noexcept { delete (_V*)p; }));
}

template <typename V>
auto as_nbarray(V&& x, const std::initializer_list<std::size_t> shape)
{
  return as_nbarray(x, shape.size(), shape.begin());
}

template <typename V>
auto as_nbarray(V&& x)
{
  return as_nbarray(std::move(x), {x.size()});
}

template <typename V, std::size_t U>
auto as_nbarrayp(std::pair<V, std::array<std::size_t, U>>&& x)
{
  return as_nbarray(std::move(x.first), x.second.size(), x.second.data());
}

template <typename T>
using ndarray1 = nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>;

template <typename T>
using ndarray2 = nb::ndarray<const T, nb::numpy, nb::shape<-1, -1>, nb::c_contig>;

std::shared_ptr<void> make_owner_from_objects(nb::object a = nb::object(),
                                              nb::object b = nb::object(),
                                              nb::object c = nb::object(),
                                              nb::object d = nb::object())
{
  struct Owner
  {
    nb::object a, b, c, d;
  };
  auto owner = std::make_shared<Owner>();
  owner->a = std::move(a);
  owner->b = std::move(b);
  owner->c = std::move(c);
  owner->d = std::move(d);
  return owner;
}

// Convert CSR connectivity+offset arrays into VTK packed cells layout
// [n0, v0_0, v0_1, ..., n1, v1_0, ...].
static std::vector<int> csr_to_vtk_cells_impl(std::span<const int> connectivity,
                                             std::span<const int> offsets)
{
  std::vector<int> out;
  if (offsets.size() == 0)
    return out;
  const std::size_t ncells = (offsets.size() > 0) ? (offsets.size() - 1) : 0;
  out.reserve(connectivity.size() + ncells);
  for (std::size_t i = 0; i < ncells; ++i)
  {
    int start = offsets[i];
    int end = offsets[i + 1];
    int n = end - start;
    out.push_back(n);
    for (int j = start; j < end; ++j)
      out.push_back(connectivity[static_cast<std::size_t>(j)]);
  }
  return out;
}

template <typename T>
cutcells::MeshView<T, int> make_mesh_view_from_numpy(
    const ndarray2<T>& coordinates,
    const ndarray1<int>& connectivity,
    const ndarray1<int>& offsets,
    const std::optional<ndarray1<int>>& cell_types,
    int tdim)
{
  cutcells::MeshView<T, int> mesh;
  mesh.gdim = static_cast<int>(coordinates.shape(1));
  mesh.tdim = tdim;

  mesh.coordinates = std::span<const T>(coordinates.data(),
                                        static_cast<std::size_t>(coordinates.size()));
  mesh.connectivity = std::span<const int>(connectivity.data(),
                                           static_cast<std::size_t>(connectivity.size()));
  mesh.offsets = std::span<const int>(offsets.data(),
                                      static_cast<std::size_t>(offsets.size()));

  // Keep-alive bundle: holds numpy arrays and the converted cell-type vector.
  struct MeshViewOwner
  {
    nb::object coords, conn, offs, types_numpy;
    std::vector<cutcells::cell::type> cell_types_converted;
  };
  auto owner_data = std::make_shared<MeshViewOwner>();
  owner_data->coords = nb::cast(coordinates);
  owner_data->conn   = nb::cast(connectivity);
  owner_data->offs   = nb::cast(offsets);

  if (cell_types.has_value())
  {
    // Convert VTK integer codes to cutcells cell::type enum at the Python boundary.
    const int*        raw    = cell_types->data();
    const std::size_t ncells = static_cast<std::size_t>(cell_types->size());
    owner_data->types_numpy  = nb::cast(*cell_types);
    owner_data->cell_types_converted.reserve(ncells);
    for (std::size_t i = 0; i < ncells; ++i)
      owner_data->cell_types_converted.push_back(
          cutcells::cell::map_vtk_type_to_cell_type(
              static_cast<cutcells::cell::vtk_types>(raw[i])));
    mesh.cell_types = std::span<const cutcells::cell::type>(
        owner_data->cell_types_converted.data(), ncells);
    // Connectivity from Python/VTK uses VTK vertex ordering.
    mesh.vtk_vertex_order = true;
  }

  mesh.owner = owner_data;
  return mesh;
}

template <typename T>
std::vector<T> parent_cell_vertex_coords_vtk(
    const cutcells::MeshView<T, int>& mesh, int cell_id)
{
  const auto ctype = mesh.cell_type(cell_id);
  const int nv = cutcells::cell::get_num_vertices(ctype);
  std::vector<T> coords(static_cast<std::size_t>(nv * mesh.gdim), T(0));

  for (int vtk_v = 0; vtk_v < nv; ++vtk_v)
  {
    const int local_v = mesh.vtk_vertex_order
                            ? vtk_v
                            : cutcells::cell::vtk_to_basix_vertex(ctype, vtk_v);
    const int node_id = mesh.cell_node(cell_id, local_v);
    const T* x = mesh.node(node_id);
    for (int d = 0; d < mesh.gdim; ++d)
      coords[static_cast<std::size_t>(vtk_v * mesh.gdim + d)] = x[d];
  }
  return coords;
}

template <typename T>
bool vertex_is_zero_for_level_set(const cutcells::AdaptCell<T>& adapt_cell,
                                  int vertex_id,
                                  int level_set_id)
{
  const std::uint64_t bit = std::uint64_t(1) << level_set_id;
  return (adapt_cell.zero_mask_per_vertex[static_cast<std::size_t>(vertex_id)] & bit) != 0;
}

struct SelectedEntity
{
  cutcells::cell::type type = cutcells::cell::type::point;
  std::vector<int> vertices;
};

template <typename T>
std::vector<SelectedEntity> part_selected_entities(
    const cutcells::HOMeshPart<T, int>& part,
    const cutcells::AdaptCell<T>& adapt_cell)
{
  if (part.expr.clauses.size() != 1 || part.expr.clauses.front().level_set_index != 0)
  {
    throw std::runtime_error(
        "HOMeshPart direct straight output currently supports only "
        "one-clause single-level-set selections");
  }

  const auto relation = part.expr.clauses.front().relation;
  std::vector<SelectedEntity> entities;

  if (part.dim == adapt_cell.tdim)
  {
    if (relation == cutcells::Relation::EqualTo)
      throw std::runtime_error("HOMeshPart: phi = 0 is not a volume selection");

    if (adapt_cell.cell_cert_tag_num_level_sets <= 0)
      throw std::runtime_error("HOMeshPart: missing cell certification tags");

    const auto target =
        (relation == cutcells::Relation::LessThan)
            ? cutcells::CellCertTag::negative
            : cutcells::CellCertTag::positive;

    const int n_cells = adapt_cell.n_entities(adapt_cell.tdim);
    entities.reserve(static_cast<std::size_t>(n_cells));
    for (int c = 0; c < n_cells; ++c)
    {
      if (adapt_cell.get_cell_cert_tag(/*level_set_id=*/0, c) == target)
      {
        auto verts = adapt_cell.entity_to_vertex[adapt_cell.tdim][static_cast<std::int32_t>(c)];
        SelectedEntity entity;
        entity.type = adapt_cell.entity_types[adapt_cell.tdim][static_cast<std::size_t>(c)];
        entity.vertices.assign(verts.begin(), verts.end());
        entities.push_back(std::move(entity));
      }
    }
    return entities;
  }

  if (relation != cutcells::Relation::EqualTo)
  {
    throw std::runtime_error(
        "HOMeshPart: lower-dimensional direct export currently supports only phi = 0");
  }

  if (part.dim < 1 || part.dim >= adapt_cell.tdim)
  {
    throw std::runtime_error("HOMeshPart: unsupported selection dimension");
  }

  if (part.dim != adapt_cell.tdim - 1)
  {
    throw std::runtime_error(
        "HOMeshPart direct straight output currently supports only codim-1 phi = 0 selections");
  }

  std::map<std::vector<int>, SelectedEntity> unique_entities;
  const int n_cells = adapt_cell.n_entities(adapt_cell.tdim);
  for (int c = 0; c < n_cells; ++c)
  {
    const auto cell_type = adapt_cell.entity_types[adapt_cell.tdim][static_cast<std::size_t>(c)];
    auto cell_verts = adapt_cell.entity_to_vertex[adapt_cell.tdim][static_cast<std::int32_t>(c)];

    if (adapt_cell.tdim == 2)
    {
      for (const auto& edge : cutcells::cell::edges(cell_type))
      {
        SelectedEntity entity;
        entity.type = cutcells::cell::type::interval;
        entity.vertices = {
            static_cast<int>(cell_verts[static_cast<std::size_t>(edge[0])]),
            static_cast<int>(cell_verts[static_cast<std::size_t>(edge[1])])};

        bool all_zero = true;
        for (int gv : entity.vertices)
        {
          if (!vertex_is_zero_for_level_set(adapt_cell, gv, /*level_set_id=*/0))
          {
            all_zero = false;
            break;
          }
        }
        if (!all_zero)
          continue;

        auto key = entity.vertices;
        std::sort(key.begin(), key.end());
        unique_entities.try_emplace(std::move(key), std::move(entity));
      }
      continue;
    }

    const int n_faces = cutcells::cell::num_faces(cell_type);
    for (int fi = 0; fi < n_faces; ++fi)
    {
      auto local_face = cutcells::cell::face_vertices(cell_type, fi);
      SelectedEntity entity;
      entity.type = cutcells::cell::face_type(cell_type, fi);
      entity.vertices.reserve(local_face.size());
      for (auto lv : local_face)
      {
        entity.vertices.push_back(
            static_cast<int>(cell_verts[static_cast<std::size_t>(lv)]));
      }

      bool all_zero = true;
      for (int gv : entity.vertices)
      {
        if (!vertex_is_zero_for_level_set(adapt_cell, gv, /*level_set_id=*/0))
        {
          all_zero = false;
          break;
        }
      }
      if (!all_zero)
        continue;

      auto key = entity.vertices;
      std::sort(key.begin(), key.end());
      unique_entities.try_emplace(std::move(key), std::move(entity));
    }
  }

  entities.reserve(unique_entities.size());
  for (auto& [key, entity] : unique_entities)
    entities.push_back(std::move(entity));
  return entities;
}

inline bool cell_type_is_simplex(cutcells::cell::type cell_type)
{
  using cutcells::cell::type;
  return cell_type == type::interval
      || cell_type == type::triangle
      || cell_type == type::tetrahedron;
}

inline cutcells::cell::type simplex_type_for_dim(int dim)
{
  using cutcells::cell::type;
  switch (dim)
  {
    case 1:
      return type::interval;
    case 2:
      return type::triangle;
    case 3:
      return type::tetrahedron;
    default:
      throw std::runtime_error("Unsupported simplex dimension");
  }
}

template <typename T>
std::vector<T> entity_reference_coords(const cutcells::AdaptCell<T>& adapt_cell,
                                       std::span<const int> entity_vertices)
{
  std::vector<T> coords(
      static_cast<std::size_t>(entity_vertices.size() * adapt_cell.tdim), T(0));

  for (std::size_t j = 0; j < entity_vertices.size(); ++j)
  {
    const int gv = entity_vertices[j];
    for (int d = 0; d < adapt_cell.tdim; ++d)
    {
      coords[static_cast<std::size_t>(j * adapt_cell.tdim + d)] =
          adapt_cell.vertex_coords[static_cast<std::size_t>(gv * adapt_cell.tdim + d)];
    }
  }

  return coords;
}

template <typename T>
std::vector<T> push_forward_parent_reference_coords(
    cutcells::cell::type parent_cell_type,
    const std::vector<T>& parent_vertex_coords_vtk,
    int parent_tdim,
    std::span<const T> reference_coords)
{
  std::vector<T> physical_coords(reference_coords.size(), T(0));
  cutcells::cell::push_forward_affine(
      parent_cell_type,
      parent_vertex_coords_vtk,
      parent_tdim,
      reference_coords,
      std::span<T>(physical_coords.data(), physical_coords.size()));
  return physical_coords;
}

template <typename T>
void gather_subcell_vertices(std::span<const T> coords,
                             int coord_dim,
                             std::span<const int> vertex_ids,
                             std::vector<T>& out)
{
  out.resize(static_cast<std::size_t>(vertex_ids.size() * coord_dim));
  for (std::size_t j = 0; j < vertex_ids.size(); ++j)
  {
    const int local_v = vertex_ids[j];
    for (int d = 0; d < coord_dim; ++d)
    {
      out[static_cast<std::size_t>(j * coord_dim + d)] =
          coords[static_cast<std::size_t>(local_v * coord_dim + d)];
    }
  }
}

template <typename T>
void map_canonical_to_subcell_points(const T* canonical_points,
                                     int num_points,
                                     int simplex_dim,
                                     const T* subcell_vertices,
                                     int parent_tdim,
                                     T* out_points)
{
  const T* v0 = subcell_vertices;

  for (int q = 0; q < num_points; ++q)
  {
    const T* X = canonical_points + q * simplex_dim;
    T* x = out_points + q * parent_tdim;

    for (int d = 0; d < parent_tdim; ++d)
      x[d] = v0[d];

    for (int i = 1; i <= simplex_dim; ++i)
    {
      const T* vi = subcell_vertices + i * parent_tdim;
      for (int d = 0; d < parent_tdim; ++d)
        x[d] += X[i - 1] * (vi[d] - v0[d]);
    }
  }
}

template <typename T>
T simplex_physical_measure(const T* vertices,
                           int simplex_dim,
                           int gdim)
{
  T J[9] = {};
  const T* v0 = vertices;

  for (int col = 0; col < simplex_dim; ++col)
  {
    const T* vi = vertices + (col + 1) * gdim;
    for (int row = 0; row < gdim; ++row)
      J[col * gdim + row] = vi[row] - v0[row];
  }

  if (simplex_dim == gdim)
  {
    if (simplex_dim == 1)
      return std::abs(J[0]);
    if (simplex_dim == 2)
      return std::abs(J[0] * J[3] - J[2] * J[1]);

    const T det =
        J[0] * (J[4] * J[8] - J[7] * J[5])
      - J[3] * (J[1] * J[8] - J[7] * J[2])
      + J[6] * (J[1] * J[5] - J[4] * J[2]);
    return std::abs(det);
  }

  T G[9] = {};
  for (int i = 0; i < simplex_dim; ++i)
  {
    for (int j = 0; j < simplex_dim; ++j)
    {
      T sum = 0;
      for (int k = 0; k < gdim; ++k)
        sum += J[i * gdim + k] * J[j * gdim + k];
      G[i * simplex_dim + j] = sum;
    }
  }

  if (simplex_dim == 1)
    return std::sqrt(G[0]);
  if (simplex_dim == 2)
    return std::sqrt(G[0] * G[3] - G[1] * G[2]);

  const T det =
      G[0] * (G[4] * G[8] - G[7] * G[5])
    - G[3] * (G[1] * G[8] - G[7] * G[2])
    + G[6] * (G[1] * G[5] - G[4] * G[2]);
  return std::sqrt(det);
}

template <typename T>
void append_mesh_entity(cutcells::mesh::CutMesh<T>& out,
                        std::span<const T> physical_coords,
                        int gdim,
                        cutcells::cell::type cell_type,
                        int parent_cell_id,
                        bool triangulate)
{
  if (out._gdim == 0)
    out._gdim = gdim;
  if (out._tdim == 0)
    out._tdim = cutcells::cell::get_tdim(cell_type);

  const int nv = static_cast<int>(physical_coords.size()) / gdim;
  const int vertex_base = out._num_vertices;
  out._vertex_coords.insert(
      out._vertex_coords.end(), physical_coords.begin(), physical_coords.end());
  out._num_vertices += nv;

  if (triangulate && !cell_type_is_simplex(cell_type))
  {
    std::vector<int> local_ids(static_cast<std::size_t>(nv));
    std::iota(local_ids.begin(), local_ids.end(), 0);

    std::vector<std::vector<int>> simplices;
    cutcells::cell::triangulation(cell_type, local_ids.data(), simplices);
    const auto simplex_type =
        simplex_type_for_dim(cutcells::cell::get_tdim(cell_type));

    for (const auto& simplex : simplices)
    {
      for (int lv : simplex)
        out._connectivity.push_back(vertex_base + lv);
      out._offset.push_back(static_cast<int>(out._connectivity.size()));
      out._types.push_back(simplex_type);
      out._parent_map.push_back(parent_cell_id);
      out._num_cells += 1;
    }
    return;
  }

  for (int lv = 0; lv < nv; ++lv)
    out._connectivity.push_back(vertex_base + lv);
  out._offset.push_back(static_cast<int>(out._connectivity.size()));
  out._types.push_back(cell_type);
  out._parent_map.push_back(parent_cell_id);
  out._num_cells += 1;
}

template <typename T>
void append_simplex_quadrature(cutcells::quadrature::QuadratureRules<T>& rules,
                               cutcells::cell::type simplex_type,
                               std::span<const T> ref_vertices,
                               std::span<const T> physical_vertices,
                               int parent_tdim,
                               int gdim,
                               int order)
{
  const auto ref_rule = cutcells::quadrature::get_reference_rule<T>(simplex_type, order);
  const int num_points = ref_rule._num_points;
  const int simplex_dim = ref_rule._tdim;

  std::vector<T> mapped_ref_points(static_cast<std::size_t>(num_points * parent_tdim), T(0));
  map_canonical_to_subcell_points(
      ref_rule._points.data(),
      num_points,
      simplex_dim,
      ref_vertices.data(),
      parent_tdim,
      mapped_ref_points.data());

  rules._points.insert(
      rules._points.end(), mapped_ref_points.begin(), mapped_ref_points.end());

  const T measure = simplex_physical_measure(
      physical_vertices.data(), simplex_dim, gdim);
  for (int q = 0; q < num_points; ++q)
    rules._weights.push_back(ref_rule._weights[q] * measure);
}

template <typename T>
void append_entity_quadrature(cutcells::quadrature::QuadratureRules<T>& rules,
                              cutcells::cell::type cell_type,
                              std::span<const T> ref_vertices,
                              std::span<const T> physical_vertices,
                              int parent_tdim,
                              int gdim,
                              int parent_cell_id,
                              int order,
                              bool triangulate)
{
  if (rules._tdim == 0)
    rules._tdim = parent_tdim;
  if (rules._offset.empty())
    rules._offset.push_back(0);

  const int cell_dim = cutcells::cell::get_tdim(cell_type);
  if (triangulate && !cell_type_is_simplex(cell_type))
  {
    const int nv = static_cast<int>(ref_vertices.size()) / parent_tdim;
    std::vector<int> local_ids(static_cast<std::size_t>(nv));
    std::iota(local_ids.begin(), local_ids.end(), 0);

    std::vector<std::vector<int>> simplices;
    cutcells::cell::triangulation(cell_type, local_ids.data(), simplices);
    const auto simplex_type = simplex_type_for_dim(cell_dim);

    std::vector<T> ref_simplex;
    std::vector<T> phys_simplex;
    for (const auto& simplex : simplices)
    {
      gather_subcell_vertices(
          ref_vertices,
          parent_tdim,
          std::span<const int>(simplex.data(), simplex.size()),
          ref_simplex);
      gather_subcell_vertices(
          physical_vertices,
          gdim,
          std::span<const int>(simplex.data(), simplex.size()),
          phys_simplex);
      append_simplex_quadrature(
          rules,
          simplex_type,
          std::span<const T>(ref_simplex.data(), ref_simplex.size()),
          std::span<const T>(phys_simplex.data(), phys_simplex.size()),
          parent_tdim,
          gdim,
          order);
    }
  }
  else if (cell_type_is_simplex(cell_type))
  {
    append_simplex_quadrature(
        rules,
        cell_type,
        ref_vertices,
        physical_vertices,
        parent_tdim,
        gdim,
        order);
  }
  else
  {
    const auto ref_rule = cutcells::quadrature::get_reference_rule<T>(cell_type, order);
    const T measure = cutcells::cell::affine_volume_factor<T>(
        cell_type, physical_vertices.data(), gdim);
    rules._points.insert(
        rules._points.end(), ref_rule._points.begin(), ref_rule._points.end());
    for (int q = 0; q < ref_rule._num_points; ++q)
      rules._weights.push_back(ref_rule._weights[q] * measure);
  }

  rules._parent_map.push_back(parent_cell_id);
  rules._offset.push_back(static_cast<int32_t>(rules._weights.size()));
}

inline bool part_mode_is_cut_only(std::string_view mode)
{
  if (mode == "cut_only")
    return true;
  if (mode == "full")
    return false;
  throw std::runtime_error("HOMeshPart mode must be 'cut_only' or 'full'");
}

template <typename T>
cutcells::mesh::CutMesh<T> part_visualization_mesh(
    const cutcells::HOMeshPart<T, int>& part,
    std::string_view mode,
    bool triangulate)
{
  const bool cut_only = part_mode_is_cut_only(mode);
  return cutcells::output::visualization_mesh(
      part, /*include_uncut_cells=*/!cut_only, triangulate);
}

template <typename T>
cutcells::quadrature::QuadratureRules<T> part_quadrature(
    const cutcells::HOMeshPart<T, int>& part,
    int order,
    std::string_view mode,
    bool triangulate)
{
  const bool cut_only = part_mode_is_cut_only(mode);
  return cutcells::output::quadrature_rules(
      part, order, /*include_uncut_cells=*/!cut_only, triangulate);
}

template <typename T>
void declare_meshview_and_levelset(nb::module_& m, const std::string& suffix)
{
  using MeshViewT = cutcells::MeshView<T, int>;
  using LevelSetMeshDataT = cutcells::LevelSetMeshData<T, int>;
  using LevelSetT = cutcells::LevelSetFunction<T, int>;

  const std::string mesh_name = "MeshView_" + suffix;
  nb::class_<MeshViewT>(m, mesh_name.c_str(), "Lightweight mesh view")
      .def(
          "__init__",
          [](MeshViewT* self,
             const ndarray2<T>& coordinates,
             const ndarray1<int>& connectivity,
             const ndarray1<int>& offsets,
             nb::object cell_types_obj,
             int tdim)
          {
            std::optional<ndarray1<int>> cell_types;
            if (!cell_types_obj.is_none())
              cell_types = nb::cast<ndarray1<int>>(cell_types_obj);

            new (self) MeshViewT(
                make_mesh_view_from_numpy<T>(coordinates, connectivity, offsets, cell_types, tdim));
          },
          nb::arg("coordinates"),
          nb::arg("connectivity"),
          nb::arg("offsets"),
          nb::arg("cell_types") = nb::none(),
          nb::arg("tdim"))
      .def_prop_ro("gdim", [](const MeshViewT& self) { return self.gdim; })
      .def_prop_ro("tdim", [](const MeshViewT& self) { return self.tdim; })
      .def("num_nodes", &MeshViewT::num_nodes)
      .def("num_cells", &MeshViewT::num_cells)
      .def("has_cell_types", &MeshViewT::has_cell_types)
      .def("cell_num_nodes", &MeshViewT::cell_num_nodes)
      .def("cell_node", &MeshViewT::cell_node)
      .def(
          "node",
          [](const MeshViewT& self, int node_id)
          {
            const T* x = self.node(node_id);
            return nb::ndarray<const T, nb::numpy>(
                x, {static_cast<std::size_t>(self.gdim)}, nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "coordinates",
          [](const MeshViewT& self)
          {
            return nb::ndarray<const T, nb::numpy>(
                self.coordinates.data(),
                {static_cast<std::size_t>(self.num_nodes()),
                 static_cast<std::size_t>(self.gdim)},
                nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "connectivity",
          [](const MeshViewT& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.connectivity.data(),
                {self.connectivity.size()},
                nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "offsets",
          [](const MeshViewT& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.offsets.data(),
                {self.offsets.size()},
                nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_types",
          [](const MeshViewT& self)
          {
            if (self.cell_types.empty())
              return nb::ndarray<const int, nb::numpy>(nullptr, {0}, nb::handle());
            // cell::type has explicit underlying type int; safe to expose as int array.
            const int* data = reinterpret_cast<const int*>(self.cell_types.data());
            return nb::ndarray<const int, nb::numpy>(
                data,
                {self.cell_types.size()},
                nb::handle());
          },
          nb::rv_policy::reference_internal);

  const std::string ls_mesh_name = "LevelSetMeshData_" + suffix;
  nb::class_<LevelSetMeshDataT>(m, ls_mesh_name.c_str(), "Discrete level-set mesh data")
      .def(nb::init<>())
      .def_prop_ro("gdim", [](const LevelSetMeshDataT& self) { return self.gdim; })
      .def_prop_ro("tdim", [](const LevelSetMeshDataT& self) { return self.tdim; })
      .def_prop_ro("degree", [](const LevelSetMeshDataT& self) { return self.degree; })
      .def("num_dofs", &LevelSetMeshDataT::num_dofs)
      .def("num_cells", &LevelSetMeshDataT::num_cells)
      .def("cell_num_dofs", &LevelSetMeshDataT::cell_num_dofs, nb::arg("cell_id"))
      .def_prop_ro(
          "dof_coordinates",
          [](const LevelSetMeshDataT& self)
          {
            const std::size_t gdim = static_cast<std::size_t>(self.gdim);
            if (gdim == 0)
              return nb::ndarray<const T, nb::numpy>(nullptr, {0, 0}, nb::handle());
            const std::size_t n = self.dof_coordinates.size() / gdim;
            return nb::ndarray<const T, nb::numpy>(
                self.dof_coordinates.data(),
                {n, gdim},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_dofs",
          [](const LevelSetMeshDataT& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.cell_dofs.data(),
                {self.cell_dofs.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_offsets",
          [](const LevelSetMeshDataT& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.cell_offsets.data(),
                {self.cell_offsets.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_types",
          [](const LevelSetMeshDataT& self)
          {
            // Convert cell::type enum values to int for Python
            auto* owner = new std::vector<int>();
            owner->reserve(self.cell_types.size());
            for (auto ct : self.cell_types)
              owner->push_back(static_cast<int>(ct));
            return nb::ndarray<int, nb::numpy>(
                owner->data(),
                {owner->size()},
                nb::capsule(owner, [](void* p) noexcept {
                  delete static_cast<std::vector<int>*>(p);
                }));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "dof_parent_dim",
          [](const LevelSetMeshDataT& self)
          {
            return nb::ndarray<const int8_t, nb::numpy>(
                self.dof_parent_dim.data(),
                {self.dof_parent_dim.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "dof_parent_id",
          [](const LevelSetMeshDataT& self)
          {
            return nb::ndarray<const int32_t, nb::numpy>(
                self.dof_parent_id.data(),
                {self.dof_parent_id.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "dof_parent_param",
          [](const LevelSetMeshDataT& self)
          {
            return nb::ndarray<const T, nb::numpy>(
                self.dof_parent_param.data(),
                {self.dof_parent_param.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "dof_parent_param_offset",
          [](const LevelSetMeshDataT& self)
          {
            return nb::ndarray<const int32_t, nb::numpy>(
                self.dof_parent_param_offset.data(),
                {self.dof_parent_param_offset.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal);

  m.def(
      "create_level_set_mesh_data",
      [](const MeshViewT& mesh, int degree, T merge_tol)
      {
        nb::gil_scoped_release release;
        return cutcells::create_level_set_mesh_data<T, int>(mesh, degree, merge_tol);
      },
      nb::arg("mesh"),
      nb::arg("degree"),
      nb::arg("merge_tol") = T(-1));

  m.def(
      "create_level_set_mesh_data",
      [](const ndarray2<T>& dof_coordinates,
         const ndarray1<int>& cell_dofs,
         const ndarray1<int>& cell_offsets,
         int degree,
         int tdim,
         nb::object cell_types_obj)
      {
        // Convert int cell types to cell::type enum
        // Python passes VTK type integers; convert to cell::type.
        std::vector<cutcells::cell::type> cell_types_vec;
        if (!cell_types_obj.is_none())
        {
          auto arr = nb::cast<ndarray1<int>>(cell_types_obj);
          cell_types_vec.reserve(static_cast<std::size_t>(arr.size()));
          for (std::size_t i = 0; i < static_cast<std::size_t>(arr.size()); ++i)
            cell_types_vec.push_back(
                cutcells::cell::map_vtk_type_to_cell_type(
                    static_cast<cutcells::cell::vtk_types>(arr.data()[i])));
        }

        std::span<const cutcells::cell::type> cell_types_span;
        if (!cell_types_vec.empty())
          cell_types_span = std::span<const cutcells::cell::type>(
              cell_types_vec.data(), cell_types_vec.size());

        return cutcells::create_level_set_mesh_data<T, int>(
            static_cast<int>(dof_coordinates.shape(1)),
            tdim,
            degree,
            std::span<const T>(dof_coordinates.data(),
                               static_cast<std::size_t>(dof_coordinates.size())),
            std::span<const int>(cell_dofs.data(),
                                 static_cast<std::size_t>(cell_dofs.size())),
            std::span<const int>(cell_offsets.data(),
                                 static_cast<std::size_t>(cell_offsets.size())),
            cell_types_span);
      },
      nb::arg("dof_coordinates"),
      nb::arg("cell_dofs"),
      nb::arg("cell_offsets"),
      nb::arg("degree"),
      nb::arg("tdim"),
      nb::arg("cell_types") = nb::none());

  const std::string ls_name = "LevelSetFunction_" + suffix;
  nb::class_<LevelSetT>(m, ls_name.c_str(), "Level-set function")
      .def(
          "__init__",
          [](LevelSetT* self,
             nb::object value_obj,
             nb::object grad_obj,
             nb::object nodal_values_obj,
             int gdim)
          {
            using ValueFn = std::function<T(const T*, int)>;
            using GradFn = std::function<void(const T*, int, T*)>;

            ValueFn value_fn;
            GradFn grad_fn;
            std::span<const T> nodal_values;
            nb::object nodal_owner;

            if (!value_obj.is_none())
            {
              nb::callable value_callable = nb::cast<nb::callable>(value_obj);

              value_fn = [value_callable](const T* x, int cell_id) -> T
              {
                nb::gil_scoped_acquire gil;
                nb::ndarray<const T, nb::numpy> x_arr(x, {static_cast<std::size_t>(3)}, nb::handle());
                try
                {
                  return nb::cast<T>(value_callable(x_arr, cell_id));
                }
                catch (const nb::python_error&)
                {
                  return nb::cast<T>(value_callable(x_arr));
                }
              };
            }

            if (!grad_obj.is_none())
            {
              nb::callable grad_callable = nb::cast<nb::callable>(grad_obj);

              grad_fn = [grad_callable](const T* x, int cell_id, T* g)
              {
                nb::gil_scoped_acquire gil;
                nb::ndarray<const T, nb::numpy> x_arr(x, {static_cast<std::size_t>(3)}, nb::handle());

                nb::object result;
                try
                {
                  result = grad_callable(x_arr, cell_id);
                }
                catch (const nb::python_error&)
                {
                  result = grad_callable(x_arr);
                }

                auto grad_arr = nb::cast<ndarray1<T>>(result);
                for (std::size_t i = 0; i < grad_arr.size(); ++i)
                  g[i] = grad_arr(i);
              };
            }

            if (!nodal_values_obj.is_none())
            {
              auto nodal_array = nb::cast<ndarray1<T>>(nodal_values_obj);
              nodal_values = std::span<const T>(nodal_array.data(),
                                                static_cast<std::size_t>(nodal_array.size()));
              nodal_owner = nodal_values_obj;
            }

            if (!value_fn && nodal_values.empty())
              throw std::runtime_error("LevelSetFunction requires at least one of value or nodal_values.");

            if (grad_fn && !value_fn)
              throw std::runtime_error("LevelSetFunction: grad requires value.");

            new (self) LevelSetT{};
            self->value_fn = std::move(value_fn);
            self->grad_fn = std::move(grad_fn);
            self->nodal_values = nodal_values;
            self->gdim = gdim;
            self->owner = make_owner_from_objects(value_obj, grad_obj, nodal_owner);
          },
          nb::arg("value") = nb::none(),
          nb::arg("grad") = nb::none(),
          nb::arg("nodal_values") = nb::none(),
          nb::arg("gdim") = 0)
      .def_prop_ro("gdim", [](const LevelSetT& self) { return self.gdim; })
      .def("has_value", &LevelSetT::has_value)
      .def("has_gradient", &LevelSetT::has_gradient)
      .def("has_nodal_values", &LevelSetT::has_nodal_values)
      .def("has_mesh_data", &LevelSetT::has_mesh_data)
      .def("has_dof_values", &LevelSetT::has_dof_values)
      .def(
          "value",
          [](const LevelSetT& self, const ndarray1<T>& x, int cell_id)
          {
            return self.value(x.data(), cell_id);
          },
          nb::arg("x"),
          nb::arg("cell_id") = -1)
      .def(
          "grad",
          [](const LevelSetT& self, const ndarray1<T>& x, int cell_id)
          {
            std::vector<T> g(static_cast<std::size_t>(self.gdim), T(0));
            self.grad(x.data(), cell_id, g.data());
            return nb::ndarray<const T, nb::numpy>(
                g.data(),
                {g.size()},
                nb::capsule(new std::vector<T>(std::move(g)),
                            [](void* p) noexcept { delete static_cast<std::vector<T>*>(p); }));
          },
          nb::arg("x"),
          nb::arg("cell_id") = -1)
      .def("value_at_node", &LevelSetT::value_at_node)
      .def_prop_ro(
          "nodal_values",
          [](const LevelSetT& self)
          {
            if (self.nodal_values.empty())
              return nb::ndarray<const T, nb::numpy>(nullptr, {0}, nb::handle());
            return nb::ndarray<const T, nb::numpy>(
                self.nodal_values.data(),
                {self.nodal_values.size()},
                nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "dof_values",
          [](const LevelSetT& self)
          {
            if (self.dof_values.empty())
              return nb::ndarray<const T, nb::numpy>(nullptr, {0}, nb::handle());
            return nb::ndarray<const T, nb::numpy>(
                self.dof_values.data(),
                {self.dof_values.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "mesh_data",
          [](const LevelSetT& self) -> nb::object
          {
            if (!self.mesh_data)
              return nb::none();
            return nb::cast(self.mesh_data);
          },
          nb::rv_policy::reference_internal);

  m.def(
      "create_level_set_function",
      [](const LevelSetMeshDataT& mesh_data, const ndarray1<T>& dof_values,
         const std::string& name)
      {
        auto mesh_data_ptr = std::make_shared<LevelSetMeshDataT>(mesh_data);
        return cutcells::create_level_set_function<T, int>(
            std::move(mesh_data_ptr),
            std::span<const T>(
                dof_values.data(),
                static_cast<std::size_t>(dof_values.size())),
            name);
      },
      nb::arg("mesh_data"),
      nb::arg("dof_values"),
      nb::arg("name") = "phi",
      "Create a polynomial LevelSetFunction from mesh_data and global dof values.");

  m.def(
      "create_level_set",
      [](const MeshViewT& mesh, nb::callable phi, int degree,
         const std::string& name)
      {
        LevelSetMeshDataT mesh_data;
        {
          nb::gil_scoped_release release;
          mesh_data = cutcells::create_level_set_mesh_data<T, int>(mesh, degree, T(-1));
        }

        const std::size_t num_dofs = static_cast<std::size_t>(mesh_data.num_dofs());
        const std::size_t gdim = static_cast<std::size_t>(mesh_data.gdim);
        // Expose as (ndim, npoints) so x[0] gives all x-coords, etc.
        // Underlying storage is row-major (npoints, gdim), so use transposed strides.
        const std::size_t shape[2] = {gdim, num_dofs};
        // Strides are in element counts (nanobind/DLPack convention).
        // data layout is (num_dofs, gdim) row-major, so to view as (gdim, num_dofs):
        //   stride[0] = 1  (adjacent coords of the same point)
        //   stride[1] = gdim  (first coord of the next point)
        const int64_t strides[2] = {1LL, static_cast<int64_t>(gdim)};
        nb::ndarray<const T, nb::numpy> x(
            mesh_data.dof_coordinates.data(),
            2, shape, nb::handle(), strides);

        // Intentional single batched callback invocation.
        nb::object values_obj = phi(x);
        auto values = nb::cast<ndarray1<T>>(values_obj);
        if (static_cast<std::size_t>(values.size()) != num_dofs)
        {
          throw std::runtime_error(
              "create_level_set: callback must return a 1D array with length num_dofs");
        }

        auto mesh_data_ptr = std::make_shared<LevelSetMeshDataT>(std::move(mesh_data));
        return cutcells::create_level_set_function<T, int>(
            std::move(mesh_data_ptr),
            std::span<const T>(
                values.data(),
                static_cast<std::size_t>(values.size())),
            name);
      },
      nb::arg("mesh"),
      nb::arg("phi"),
      nb::arg("degree"),
      nb::arg("name") = "phi",
      "Interpolate a batched callable phi(X) at higher-order level-set dof coordinates.");

  // ---- cut_mesh_view ----
  // Cuts a MeshView using a LevelSetFunction, returning a CutMesh.
  // Nodal level-set values are taken from ls.nodal_values if available,
  // otherwise ls.value() is evaluated at every mesh node.
  m.def(
      "cut_mesh_view",
      [](const MeshViewT& mesh, const LevelSetT& ls,
         const std::string& cut_type_str, bool triangulate)
      {
        if (!mesh.has_cell_types())
          throw std::runtime_error(
              "cut_mesh_view: MeshView must have cell_types (VTK type IDs)");

        // Build nodal level-set values while GIL is held
        // (ls.value() may call back into Python)
        const std::size_t n = static_cast<std::size_t>(mesh.num_nodes());
        std::vector<T> ls_vals(n);
        if (ls.has_nodal_values())
        {
          if (ls.nodal_values.size() != n)
            throw std::runtime_error(
                "cut_mesh_view: nodal_values size does not match MeshView num_nodes");
          std::copy(ls.nodal_values.begin(), ls.nodal_values.end(), ls_vals.begin());
        }
        else if (ls.has_value())
        {
          for (std::size_t i = 0; i < n; ++i)
            ls_vals[i] = ls.value(mesh.node(static_cast<int>(i)), -1);
        }
        else
        {
          throw std::runtime_error(
              "cut_mesh_view: LevelSetFunction has neither value nor nodal_values");
        }

        // Cut mesh with GIL released
        // Convert cell::type back to VTK int codes for the legacy cut_vtk_mesh API.
        std::vector<int> vtk_types_vec;
        vtk_types_vec.reserve(mesh.cell_types.size());
        for (auto ct : mesh.cell_types)
            vtk_types_vec.push_back(
                static_cast<int>(cutcells::cell::map_cell_type_to_vtk(ct)));

        nb::gil_scoped_release release;
        return mesh::cut_vtk_mesh<T>(
            std::span<const T>(ls_vals.data(), ls_vals.size()),
            mesh.coordinates,
            mesh.connectivity,
            mesh.offsets,
            std::span<const int>(vtk_types_vec.data(), vtk_types_vec.size()),
            cut_type_str,
            triangulate);
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("cut_type"),
      nb::arg("triangulate") = true,
      "Cut a MeshView with a LevelSetFunction.\n"
      "Returns a CutMesh containing cells classified by cut_type (\"phi<0\", \"phi=0\", \"phi>0\").\n"
      "Level-set values are taken from nodal_values if set, otherwise evaluated via value().");

  // Simple aliases for Python
  if constexpr (std::is_same_v<T, double>)
  {
    m.def(
        "write_level_set_vtu",
        [](const std::string& filename, const LevelSetT& ls, const std::string& field_name)
        {
          nb::gil_scoped_release release;
          io::write_level_set_vtu(filename, ls, field_name);
        },
        nb::arg("filename"),
        nb::arg("level_set"),
        nb::arg("field_name") = "phi");

    m.attr("MeshView") = m.attr(mesh_name.c_str());
    m.attr("LevelSetMeshData") = m.attr(ls_mesh_name.c_str());
    m.attr("LevelSetFunction") = m.attr(ls_name.c_str());
  }
}

template <typename T>
void declare_float(nb::module_& m, std::string type)
{
    m.def("classify_cell_domain", [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_values){
          cell::domain domain_id = cell::classify_cell_domain<T>(std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())});
          auto domain_str = cell_domain_to_str(domain_id);
          return domain_str;
        }
        , "classify a cell domain");

    std::string name = "CutCell_" + type;
    nb::class_<cell::CutCell<T>>(m, name.c_str(), "Cut Cell")
        .def(nb::init<>())
        //shape and classes for cutcell vertex coords, connectivity and types are for visualization with pyvista
        .def_prop_ro(
          "vertex_coords",
          [](const cell::CutCell<T>& self) {
            const std::size_t gdim = static_cast<std::size_t>(self._gdim);
            if (gdim == 0)
              return nb::ndarray<const T, nb::numpy>(nullptr, {0, 0}, nb::handle());
            const std::size_t n = self._vertex_coords.size() / gdim;
            return nb::ndarray<const T, nb::numpy>(
              self._vertex_coords.data(),
              {n, gdim},
              nb::handle());
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of cut-cell vertex coordinates as shape (num_vertices, gdim).")
        .def_prop_ro(
          "parent_vertex_coords",
          [](const cell::CutCell<T>& self) {
            const std::size_t gdim = static_cast<std::size_t>(self._gdim);
            if (gdim == 0)
              return nb::ndarray<const T, nb::numpy>(nullptr, {0, 0}, nb::handle());
            const std::size_t n = self._parent_vertex_coords.size() / gdim;
            return nb::ndarray<const T, nb::numpy>(
              self._parent_vertex_coords.data(),
              {n, gdim},
              nb::handle());
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of parent-cell vertex coordinates as shape (num_parent_vertices, gdim).")
        .def_prop_ro(
          "parent_vertex_ids",
          [](const cell::CutCell<T>& self) {
            return nb::ndarray<const int, nb::numpy>(
              self._parent_vertex_ids.data(),
              {self._parent_vertex_ids.size()},
              nb::handle());
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of parent vertex ids (context-global indices when available).")
        .def_prop_ro(
          "connectivity",
          [](const cell::CutCell<T>& self) {
            return nb::ndarray<const int, nb::numpy, nb::shape<-1>, nb::c_contig>(
              self._connectivity.data(),
              {self._connectivity.size()},
              nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of flat CSR connectivity array.")
        .def_prop_ro(
          "offsets",
          [](const cell::CutCell<T>& self) {
            return nb::ndarray<const int, nb::numpy, nb::shape<-1>, nb::c_contig>(
              self._offset.data(),
              {self._offset.size()},
              nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of CSR offsets array.")
        .def_prop_ro(
          "cells",
          [](const cell::CutCell<T>& self) {
            return as_nbarray(csr_to_vtk_cells_impl(
              std::span<const int>(self._connectivity.data(), self._connectivity.size()),
              std::span<const int>(self._offset.data(), self._offset.size())));
          },
          nb::rv_policy::move,
          "Packed VTK cells array [n0, v0..., n1, v1..., ...] built from connectivity+offsets.")
        .def_prop_ro(
          "types",
          [](const cell::CutCell<T>& self) {
            using type_id_t = std::underlying_type_t<cell::type>;
            static_assert(std::is_integral_v<type_id_t>, "cell::type must have integral underlying type");
            return nb::ndarray<const type_id_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              reinterpret_cast<const type_id_t*>(self._types.data()),
              {self._types.size()},
              nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of cut-cell type ids (cell::type enum underlying values).")
        .def_prop_ro(
          "vtk_types",
          [](const cell::CutCell<T>& self) {
            std::vector<uint8_t> vtk;
            vtk.reserve(self._types.size());
            for (const auto t : self._types)
              vtk.push_back(static_cast<uint8_t>(cell::map_cell_type_to_vtk(t)));
            return as_nbarray(std::move(vtk));
          },
          nb::rv_policy::move,
          "VTK type IDs for each sub-cell (uint8), suitable for pv.UnstructuredGrid.")
        .def_prop_ro(
          "vertex_parent_entity",
          [](const cell::CutCell<T>& self) {
            return nb::ndarray<const int32_t, nb::numpy>(
              self._vertex_parent_entity.data(),
              {self._vertex_parent_entity.size()},
              nb::handle());
          },
          nb::rv_policy::reference_internal,
          "Return parent entity token for each cut-cell vertex.\n"
          "Tokens encode the origin: edge intersections use edge id, original vertices use 100+vid, special points use 200+sid.")
        .def_prop_ro(
          "vertex_coords_phys",
          [](const cell::CutCell<T>& self) {
            const std::size_t gdim = static_cast<std::size_t>(self._gdim);
            if (gdim == 0 || self._vertex_coords_phys.empty())
              return nb::ndarray<const T, nb::numpy>(nullptr, {0, 0}, nb::handle());
            const std::size_t n = self._vertex_coords_phys.size() / gdim;
            return nb::ndarray<const T, nb::numpy>(
              self._vertex_coords_phys.data(),
              {n, gdim},
              nb::handle());
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of cut-cell physical vertex coordinates as shape (num_vertices, gdim). "
          "Empty until compute_physical_cut_vertices() or complete_from_physical() is called.")
        .def("str", [](const cell::CutCell<T>& self) {cell::str(self); return ;})
        .def("volume", [](const cell::CutCell<T>& self) {return cell::volume(self);})
        .def("write_vtk", [](cell::CutCell<double>& self, std::string fname) {io::write_vtk(fname,self); return ;});

    name = "CutCells_" + type;
    nb::class_<mesh::CutCells<T>>(m, name.c_str(), "Cut Cells")
        .def(nb::init<>())
        .def_prop_ro(
          "cut_cells",
          [](const mesh::CutCells<T>& self)
          {
            return self._cut_cells;
          },
          nb::rv_policy::reference_internal,
          "Return vector of cut cells.")
        .def_prop_ro(
          "types",
          [](const mesh::CutCells<T>& self) {
            using type_id_t = std::underlying_type_t<cell::type>;
            static_assert(std::is_integral_v<type_id_t>, "cell::type must have integral underlying type");
            return nb::ndarray<const type_id_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              reinterpret_cast<const type_id_t*>(self._types.data()),
              {self._types.size()},
              nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of cut-cell type ids (cell::type enum underlying values).")
        .def_prop_ro(
          "parent_map",
          [](const mesh::CutCells<T>& self) {
            return nb::ndarray<const int32_t, nb::numpy>(self._parent_map.data(),{self._parent_map.size()}, nb::handle());
          },
          nb::rv_policy::reference_internal,
          " Return parent map of cut cells.");

  name = "CutMesh_" + type;
  nb::class_<mesh::CutMesh<T>>(m, name.c_str(), "Cut Mesh")
        .def(nb::init<>())
        //shape and classes for cutcell vertex coords, connectivity and types are for visualization with pyvista
        .def_prop_ro(
          "vertex_coords",
          [](const mesh::CutMesh<T>& self) {
            const std::size_t gdim = static_cast<std::size_t>(self._gdim);
            if (gdim == 0)
              return nb::ndarray<T, nb::numpy>(nullptr, {0, 0}, nb::handle());
            const std::size_t n = self._vertex_coords.size() / gdim;
            // Return an *owned* copy so pyvista (which retains the array) does not
            // keep the CutMesh alive via a numpy base reference at interpreter shutdown.
            std::vector<T> copy = self._vertex_coords;
            return as_nbarray(std::move(copy), {n, gdim});
          },
          nb::rv_policy::move,
          "Copy of mesh vertex coordinates as shape (num_vertices, gdim).")
        .def_prop_ro(
          "connectivity",
          [](const mesh::CutMesh<T>& self) {
            return nb::ndarray<const int, nb::numpy>(self._connectivity.data(),{self._connectivity.size()}, nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          " Return connectivity vector.")
        .def_prop_ro(
          "offset",
          [](const mesh::CutMesh<T>& self) {
            return nb::ndarray<const int, nb::numpy>(self._offset.data(),{self._offset.size()}, nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          " Return offset vector.")
        .def_prop_ro(
          "types",
          [](const mesh::CutMesh<T>& self) {
            using type_id_t = std::underlying_type_t<cell::type>;
            static_assert(std::is_integral_v<type_id_t>, "cell::type must have integral underlying type");
            return nb::ndarray<const type_id_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              reinterpret_cast<const type_id_t*>(self._types.data()),
              {self._types.size()},
              nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          "Zero-copy view of cut-mesh type ids (cell::type enum underlying values).")
        .def_prop_ro(
          "vtk_types",
          [](const mesh::CutMesh<T>& self) {
            std::vector<uint8_t> vtk;
            vtk.reserve(self._types.size());
            for (const auto t : self._types)
              vtk.push_back(static_cast<uint8_t>(cell::map_cell_type_to_vtk(t)));
            return as_nbarray(std::move(vtk));
          },
          nb::rv_policy::move,
          "VTK type IDs for each sub-cell (uint8), suitable for pv.UnstructuredGrid.")
        .def_prop_ro(
          "cells",
          [](const mesh::CutMesh<T>& self) {
            return as_nbarray(csr_to_vtk_cells_impl(
              std::span<const int>(self._connectivity.data(), self._connectivity.size()),
              std::span<const int>(self._offset.data(), self._offset.size())));
          },
          nb::rv_policy::move,
          "Convenience packed cells view [n, v0, ...] for VTK-style consumers (allocates)."
          " For high-performance workflows, prefer zero-copy 'connectivity' + 'offset'.")
        .def_prop_ro(
          "parent_map",
          [](const mesh::CutMesh<T>& self) {
            return nb::ndarray<const int32_t, nb::numpy>(self._parent_map.data(),{self._parent_map.size()}, nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          " Return parent map of cut mesh.");

  // ---- frame-completion helpers ----
  m.def("compute_physical_cut_vertices",
    [](cell::CutCell<T>& cut_cell) {
      nb::gil_scoped_release release;
      cell::compute_physical_cut_vertices(cut_cell);
    },
    nb::arg("cut_cell"),
    "Fill _vertex_coords_phys via affine push-forward of _vertex_coords (reference-frame cut path).");

  m.def("complete_from_physical",
    [](cell::CutCell<T>& cut_cell) {
      nb::gil_scoped_release release;
      cell::complete_from_physical(cut_cell);
    },
    nb::arg("cut_cell"),
    "Copy vertex_coords to _vertex_coords_phys, then pull back to fill _vertex_coords "
    "with parent reference coordinates (physical-frame cut path).");

  // ---- QuadratureRules class ----
  {
    std::string qr_name = "QuadratureRules_" + type;
    nb::class_<quadrature::QuadratureRules<T>>(m, qr_name.c_str(),
        "Flat batch quadrature rules for a collection of cut cells.\n"
        "Points are in parent reference space; weights incorporate the physical Jacobian determinant.")
      .def(nb::init<>())
      .def_prop_ro("tdim",
        [](const quadrature::QuadratureRules<T>& self) { return self._tdim; },
        "Topological dimension of the point coordinates.")
      .def_prop_ro("points",
        [](const quadrature::QuadratureRules<T>& self) {
          // Always return a flat 1-D array; caller reshapes with [:, tdim] if needed.
          // Use nb::cast(self, ...) as the ndarray owner so NumPy holds a proper
          // strong reference to the parent — avoids keep_alive cycles at shutdown.
          return nb::ndarray<const T, nb::numpy>(
            self._points.data(), {self._points.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal,
        "Quadrature points in parent reference space, flat array of length (total_points * tdim). "
        "Reshape to (-1, tdim) to get shape (total_points, tdim).")
      .def_prop_ro("weights",
        [](const quadrature::QuadratureRules<T>& self) {
          return nb::ndarray<const T, nb::numpy>(
            self._weights.data(), {self._weights.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal,
        "Physical integration weights, shape (total_points,).")
      .def_prop_ro("offset",
        [](const quadrature::QuadratureRules<T>& self) {
          return nb::ndarray<const int32_t, nb::numpy>(
            self._offset.data(), {self._offset.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal,
        "Offsets into points/weights per cut-cell rule, shape (num_rules+1,).")
      .def_prop_ro("parent_map",
        [](const quadrature::QuadratureRules<T>& self) {
          return nb::ndarray<const int32_t, nb::numpy>(
            self._parent_map.data(), {self._parent_map.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal,
        "Index of the originating cut-cell in the input list, shape (num_rules,).");
  }

  // ---- make_quadrature ----
  m.def("make_quadrature",
    [](const std::vector<cell::CutCell<T>>& cut_cells, int order) {
      nb::gil_scoped_release release;
      return quadrature::make_quadrature(cut_cells, order);
    },
    nb::arg("cut_cells"), nb::arg("order"),
    "Generate flat quadrature rules for a list of enriched CutCells.\n"
    "Both vertex_coords (_vertex_coords, reference) and vertex_coords_phys must be populated on each cell.\n"
    "Returns a QuadratureRules_<T> object.");

  m.def("create_cut_mesh", [](mesh::CutCells<T>& cut_cells){
              nb::gil_scoped_release release;
              return mesh::create_cut_mesh(cut_cells);
             }
             , "Creating a cut mesh");
  m.def("cut", [](cell::type cell_type,
                   const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& vertex_coordinates,
                   const int gdim,
                   const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_values,
                   const std::string& cut_type_str,
                   bool triangulate){
              cell::CutCell<T> cut_cell;
              nb::gil_scoped_release release;
              cell::cut<T>(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, cut_cell, triangulate);
              return cut_cell;
             }
             , "cut a cell");

  m.def("higher_order_cut", [](cell::type cell_type,
             const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& vertex_coordinates,
             const int gdim,
             const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_values,
             const std::string& cut_type_str,
             bool triangulate){
              nb::gil_scoped_release release;
              cell::CutCell<T> cut_cell = cell::higher_order_cut<T>(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, triangulate);
              return cut_cell;
             }
             , "cut a second order cell");

    m.def("locate_cells", [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_vals,
                             const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& points,
                             const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
                             const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offset,
                             const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& vtk_type,
                             const std::string& cut_type_str){
              std::vector<int> located_cells;
              {
                nb::gil_scoped_release release;
                located_cells = mesh::locate_cells<T>(std::span(ls_vals.data(),ls_vals.size()),
                              std::span(points.data(),points.size()),
                              std::span(connectivity.data(),connectivity.size()),
                              std::span(offset.data(),offset.size()),
                              std::span(vtk_type.data(),vtk_type.size()),
                              cell::string_to_cut_type(cut_type_str));
              }
              return as_nbarray(std::move(located_cells));
             }
             , "locate cells in vtk mesh");

    m.def("cut_vtk_mesh", [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_vals,
                             const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& points,
                             const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
                             const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offset,
                             const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& vtk_type,
                             const std::string& cut_type_str,
                             bool triangulate){
              nb::gil_scoped_release release;
              return  mesh::cut_vtk_mesh<T>(std::span(ls_vals.data(),ls_vals.size()),
                            std::span(points.data(),points.size()),
                            std::span(connectivity.data(),connectivity.size()),
                            std::span(offset.data(),offset.size()),
                            std::span(vtk_type.data(),vtk_type.size()),
                            cut_type_str,
                            triangulate);
             }
             , nb::arg("ls_vals"), nb::arg("points"), nb::arg("connectivity"), nb::arg("offset"), nb::arg("vtk_type"),
               nb::arg("cut_type_str"), nb::arg("triangulate") = true
             , "cut vtk mesh");

    m.def("runtime_quadrature",
      [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_vals,
         const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& points,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offset,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& vtk_type,
         const std::string& cut_type_str,
         bool triangulate,
         int order) {
          {
            nb::gil_scoped_release release;
            return quadrature::runtime_quadrature<T>(
                std::span(ls_vals.data(),     ls_vals.size()),
                std::span(points.data(),      points.size()),
                std::span(connectivity.data(),connectivity.size()),
                std::span(offset.data(),      offset.size()),
                std::span(vtk_type.data(),    vtk_type.size()),
                cut_type_str,
                triangulate,
                order);
          }
      },
      nb::arg("ls_vals"), nb::arg("points"), nb::arg("connectivity"),
      nb::arg("offset"), nb::arg("vtk_type"), nb::arg("cut_type_str"),
      nb::arg("triangulate") = true, nb::arg("order") = 3,
      "Generate flat quadrature rules for all mesh cells in the requested "
      "level-set domain.\n"
      "Full cells (entirely inside/outside) use a direct reference rule scaled "
      "by |det J|.\n"
      "Cut cells (intersected) are cut and integrated via append_quadrature.\n"
      "Returns a QuadratureRules object with reference-space points and "
      "physical weights.");

    m.def("physical_points",
      [](const quadrature::QuadratureRules<T>& rules,
         const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& points,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offset,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& vtk_type) {
          std::vector<T> pts;
          {
            nb::gil_scoped_release release;
            pts = quadrature::physical_points<T>(
                rules,
                std::span(points.data(),      points.size()),
                std::span(connectivity.data(),connectivity.size()),
                std::span(offset.data(),      offset.size()),
                std::span(vtk_type.data(),    vtk_type.size()));
          }
          return as_nbarray(std::move(pts));
      },
      nb::arg("rules"), nb::arg("points"), nb::arg("connectivity"),
      nb::arg("offset"), nb::arg("vtk_type"),
      "Map reference-space quadrature points to physical space.\n"
      "Returns a flat numpy array of shape (total_num_points * 3,).");
}

// ---- HO cut types (BackgroundMeshData, HOCutCells, HOMeshPart) ----

template <typename T>
void declare_ho_cut(nb::module_& m, const std::string& type)
{
    using MeshViewT = cutcells::MeshView<T, int>;
    using LevelSetT = cutcells::LevelSetFunction<T, int>;
    using BGDataT = cutcells::BackgroundMeshData<T, int>;
    using HOCutT = cutcells::HOCutCells<T, int>;
    using PartT = cutcells::HOMeshPart<T, int>;

    // --- Wrapper that owns both HOCutCells + BackgroundMeshData ---
    struct HOCutResult
    {
        HOCutT cut_cells;
        BGDataT bg;
        std::shared_ptr<void> level_set_owner;
    };

    std::string result_name = "HOCutResult_" + type;
    auto py_result = nb::class_<HOCutResult>(m, result_name.c_str(),
        "Result of ho cut(): holds HOCutCells and BackgroundMeshData.");

    py_result
        .def_prop_ro("num_cut_cells",
            [](const HOCutResult& self) { return self.cut_cells.num_cut_cells(); })
        .def_prop_ro("num_cells",
            [](const HOCutResult& self) { return self.bg.num_cells; })
        .def_prop_ro("num_level_sets",
            [](const HOCutResult& self) { return self.bg.num_level_sets; })
        .def_prop_ro("level_set_names",
            [](const HOCutResult& self) { return self.bg.level_set_names; })
        .def_prop_ro("parent_cell_ids",
            [](const HOCutResult& self) {
                return nb::ndarray<const int, nb::numpy>(
                    self.cut_cells.parent_cell_ids.data(),
                    {self.cut_cells.parent_cell_ids.size()},
                    nb::handle());
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro("cell_domains",
            [](const HOCutResult& self) {
                const int* data = reinterpret_cast<const int*>(
                    self.bg.cell_domains.data());
                return nb::ndarray<const int, nb::numpy>(
                    data,
                    {static_cast<std::size_t>(self.bg.num_level_sets),
                     static_cast<std::size_t>(self.bg.num_cells)},
                    nb::handle());
            },
            nb::rv_policy::reference_internal,
            "Per-level-set domain classification, shape (num_level_sets, num_cells).")
        .def("__getitem__",
            [](const HOCutResult& self, const std::string& expr_str) {
                return cutcells::select_part(
                    self.cut_cells,
                    self.bg,
                    std::string_view(expr_str));
            },
            nb::arg("expr"),
            nb::keep_alive<0, 1>(),
            "Select a mesh part via expression, e.g. result[\"phi < 0\"].")
        .def(
            "adapt_cell",
            [](const HOCutResult& self, int cut_cell_id) -> const cutcells::AdaptCell<T>&
            {
                if (cut_cell_id < 0
                    || cut_cell_id >= static_cast<int>(self.cut_cells.adapt_cells.size()))
                {
                    throw std::out_of_range("adapt_cell: cut_cell_id out of range");
                }
                return self.cut_cells.adapt_cells[static_cast<std::size_t>(cut_cell_id)];
            },
            nb::arg("cut_cell_id"),
            nb::rv_policy::reference_internal,
            "Return the AdaptCell for a cut-cell index.");

    // --- HOMeshPart ---
    std::string part_name = "HOMeshPart_" + type;
    nb::class_<PartT>(m, part_name.c_str(), "Mesh part selected by expression")
        .def_prop_ro("dim",
            [](const PartT& self) { return self.dim; })
        .def_prop_ro("cut_only",
            [](const PartT& self) { return self.cut_only; })
        .def_prop_ro("num_cut_cells",
            [](const PartT& self) {
                return static_cast<int>(self.cut_cell_ids.size());
            })
        .def_prop_ro("num_uncut_cells",
            [](const PartT& self) {
                return static_cast<int>(self.uncut_cell_ids.size());
            })
        .def_prop_ro("cut_cell_ids",
            [](const PartT& self) {
                return nb::ndarray<const std::int32_t, nb::numpy>(
                    self.cut_cell_ids.data(),
                    {self.cut_cell_ids.size()},
                    nb::handle());
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro("uncut_cell_ids",
            [](const PartT& self) {
                return nb::ndarray<const int, nb::numpy>(
                    self.uncut_cell_ids.data(),
                    {self.uncut_cell_ids.size()},
                    nb::handle());
            },
            nb::rv_policy::reference_internal)
        .def(
            "visualization_mesh",
            [](const PartT& self, const std::string& mode, bool triangulate) {
                nb::gil_scoped_release release;
                return part_visualization_mesh(self, mode, triangulate);
            },
            nb::arg("mode") = "full",
            nb::arg("triangulate") = false,
            "Return a straight visualization mesh for a one-level-set HOMeshPart.\n"
            "This bridge currently supports only single-level-set selections.")
        .def(
            "quadrature",
            [](const PartT& self, int order, const std::string& mode, bool triangulate) {
                nb::gil_scoped_release release;
                return part_quadrature(self, order, mode, triangulate);
            },
            nb::arg("order") = 3,
            nb::arg("mode") = "full",
            nb::arg("triangulate") = false,
            "Return straight quadrature rules for a one-level-set HOMeshPart.\n"
            "This bridge currently supports only single-level-set selections.")
        .def(
            "write_vtu",
            [](const PartT& self,
               const std::string& filename,
               const std::string& mode,
               bool triangulate) {
                nb::gil_scoped_release release;
                auto vis = part_visualization_mesh(self, mode, triangulate);
                std::vector<double> coords(
                    vis._vertex_coords.begin(), vis._vertex_coords.end());
                io::write_vtk(
                    filename,
                    std::span<const double>(coords.data(), coords.size()),
                    std::span<const int>(vis._connectivity.data(), vis._connectivity.size()),
                    std::span<const int>(vis._offset.data(), vis._offset.size()),
                    std::span<cell::type>(vis._types.data(), vis._types.size()),
                    vis._gdim);
            },
            nb::arg("filename"),
            nb::arg("mode") = "full",
            nb::arg("triangulate") = false,
            "Write a straight visualization VTU file for a one-level-set HOMeshPart.");

    // --- ho_cut() factory ---
    m.def("ho_cut",
        [](const MeshViewT& mesh, const LevelSetT& ls) {
            nb::gil_scoped_release release;
            auto owned_ls = std::make_shared<LevelSetT>(ls);
            auto [hc, bg] = cutcells::cut(mesh, *owned_ls);
            return HOCutResult{std::move(hc), std::move(bg), owned_ls};
        },
        nb::arg("mesh"), nb::arg("level_set"),
        "Cut a MeshView with a single LevelSetFunction (HO pipeline).\n"
        "Returns an HOCutResult; use result[\"phi < 0\"] to select parts.");

    m.def("ho_cut",
        [](const MeshViewT& mesh, const std::vector<LevelSetT>& level_sets) {
            nb::gil_scoped_release release;
            auto owned_ls = std::make_shared<std::vector<LevelSetT>>(level_sets);
            auto [hc, bg] = cutcells::cut(mesh, *owned_ls);
            return HOCutResult{std::move(hc), std::move(bg), owned_ls};
        },
        nb::arg("mesh"), nb::arg("level_sets"),
        "Cut a MeshView with multiple LevelSetFunctions (HO pipeline).\n"
        "Returns an HOCutResult; use result[\"phi1 < 0 and phi2 = 0\"] to select parts.");

    m.def("cut",
        [](const MeshViewT& mesh, const LevelSetT& ls) {
            nb::gil_scoped_release release;
            auto owned_ls = std::make_shared<LevelSetT>(ls);
            auto [hc, bg] = cutcells::cut(mesh, *owned_ls);
            return HOCutResult{std::move(hc), std::move(bg), owned_ls};
        },
        nb::arg("mesh"), nb::arg("level_set"),
        "Cut a MeshView with a single LevelSetFunction (HO pipeline).\n"
        "Returns an HOCutResult; use result[\"phi < 0\"] to select parts.");

    m.def("cut",
        [](const MeshViewT& mesh, const std::vector<LevelSetT>& level_sets) {
            nb::gil_scoped_release release;
            auto owned_ls = std::make_shared<std::vector<LevelSetT>>(level_sets);
            auto [hc, bg] = cutcells::cut(mesh, *owned_ls);
            return HOCutResult{std::move(hc), std::move(bg), owned_ls};
        },
        nb::arg("mesh"), nb::arg("level_sets"),
        "Cut a MeshView with multiple LevelSetFunctions (HO pipeline).\n"
        "Returns an HOCutResult; use result[\"phi1 < 0 and phi2 = 0\"] to select parts.");

    // Simple aliases for Python
    if constexpr (std::is_same_v<T, double>)
    {
        m.attr("HOCutResult") = m.attr(result_name.c_str());
        m.attr("HOMeshPart") = m.attr(part_name.c_str());
    }
}

template <typename T>
void declare_certification(nb::module_& m, const std::string& suffix)
{
    using MeshViewT = cutcells::MeshView<T, int>;
    using LevelSetT = cutcells::LevelSetFunction<T, int>;
    using AdaptCellT = cutcells::AdaptCell<T>;
    using LevelSetCellT = cutcells::LevelSetCell<T, int>;

    const std::string adapt_name = "AdaptCell_" + suffix;
    nb::class_<AdaptCellT>(m, adapt_name.c_str(), "Adaptive local cell topology")
        .def(nb::init<>())
        .def_prop_ro("tdim", [](const AdaptCellT& self) { return self.tdim; })
        .def("num_vertices", &AdaptCellT::n_vertices)
        .def("num_edges", [](const AdaptCellT& self) { return self.n_entities(1); })
        .def("num_cells", [](const AdaptCellT& self) { return self.n_entities(self.tdim); })
        .def_prop_ro(
            "vertex_coords",
            [](const AdaptCellT& self)
            {
                return nb::ndarray<const T, nb::numpy>(
                    self.vertex_coords.data(),
                    {static_cast<std::size_t>(self.n_vertices()),
                     static_cast<std::size_t>(self.tdim)},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "vertex_source_edge_id",
            [](const AdaptCellT& self)
            {
                return nb::ndarray<const std::int32_t, nb::numpy>(
                    self.vertex_source_edge_id.data(),
                    {self.vertex_source_edge_id.size()},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "edge_connectivity",
            [](const AdaptCellT& self)
            {
                return nb::ndarray<const std::int32_t, nb::numpy>(
                    self.entity_to_vertex[1].indices.data(),
                    {self.entity_to_vertex[1].indices.size()},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "edge_offsets",
            [](const AdaptCellT& self)
            {
                return nb::ndarray<const std::int32_t, nb::numpy>(
                    self.entity_to_vertex[1].offsets.data(),
                    {self.entity_to_vertex[1].offsets.size()},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "cell_types",
            [](const AdaptCellT& self)
            {
                std::vector<int> types;
                const int n_cells = self.n_entities(self.tdim);
                types.reserve(static_cast<std::size_t>(n_cells));
                for (int c = 0; c < n_cells; ++c)
                    types.push_back(static_cast<int>(self.entity_types[self.tdim][static_cast<std::size_t>(c)]));
                return as_nbarray(std::move(types));
            },
            nb::rv_policy::move)
        .def_prop_ro(
            "cell_connectivity",
            [](const AdaptCellT& self)
            {
                return nb::ndarray<const std::int32_t, nb::numpy>(
                    self.entity_to_vertex[self.tdim].indices.data(),
                    {self.entity_to_vertex[self.tdim].indices.size()},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "cell_offsets",
            [](const AdaptCellT& self)
            {
                return nb::ndarray<const std::int32_t, nb::numpy>(
                    self.entity_to_vertex[self.tdim].offsets.data(),
                    {self.entity_to_vertex[self.tdim].offsets.size()},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal)
        .def(
            "edge_root_tags",
            [](const AdaptCellT& self, int level_set_id)
            {
                std::vector<int> tags;
                const int n_edges = self.n_entities(1);
                tags.reserve(static_cast<std::size_t>(n_edges));
                for (int e = 0; e < n_edges; ++e)
                    tags.push_back(static_cast<int>(self.get_edge_root_tag(level_set_id, e)));
                return as_nbarray(std::move(tags));
            },
            nb::arg("level_set_id"))
        .def(
            "cell_cert_tags",
            [](const AdaptCellT& self, int level_set_id)
            {
                std::vector<int> tags;
                const int n_cells = self.n_entities(self.tdim);
                tags.reserve(static_cast<std::size_t>(n_cells));
                for (int c = 0; c < n_cells; ++c)
                    tags.push_back(static_cast<int>(self.get_cell_cert_tag(level_set_id, c)));
                return as_nbarray(std::move(tags));
            },
            nb::arg("level_set_id"))
        .def(
            "edge_green_split_params",
            [](const AdaptCellT& self, int level_set_id)
            {
                const int n_edges = self.n_entities(1);
                std::vector<T> params(static_cast<std::size_t>(n_edges), T(0));
                for (int e = 0; e < n_edges; ++e)
                {
                    const auto idx = static_cast<std::size_t>(level_set_id * n_edges + e);
                    if (idx < self.edge_green_split_param.size())
                        params[static_cast<std::size_t>(e)] = self.edge_green_split_param[idx];
                }
                return as_nbarray(std::move(params));
            },
            nb::arg("level_set_id"))
        .def(
            "edge_green_split_mask",
            [](const AdaptCellT& self, int level_set_id)
            {
                const int n_edges = self.n_entities(1);
                std::vector<std::uint8_t> mask(static_cast<std::size_t>(n_edges), 0);
                for (int e = 0; e < n_edges; ++e)
                {
                    const auto idx = static_cast<std::size_t>(level_set_id * n_edges + e);
                    if (idx < self.edge_green_split_has_value.size())
                        mask[static_cast<std::size_t>(e)] = self.edge_green_split_has_value[idx];
                }
                return as_nbarray(std::move(mask));
            },
            nb::arg("level_set_id"));

    const std::string lsc_name = "LevelSetCell_" + suffix;
    nb::class_<LevelSetCellT>(m, lsc_name.c_str(), "Cell-local Bernstein level set")
        .def(nb::init<>())
        .def_prop_ro("cell_id", [](const LevelSetCellT& self) { return self.cell_id; })
        .def_prop_ro("bernstein_order", [](const LevelSetCellT& self) { return self.bernstein_order; })
        .def_prop_ro("cell_type", [](const LevelSetCellT& self) { return self.cell_type; })
        .def_prop_ro(
            "bernstein_coeffs",
            [](const LevelSetCellT& self)
            {
                return nb::ndarray<const T, nb::numpy>(
                    self.bernstein_coeffs.data(),
                    {self.bernstein_coeffs.size()},
                    nb::cast(self, nb::rv_policy::reference));
            },
            nb::rv_policy::reference_internal);

    m.def(
        "make_adapt_cell",
        [](const MeshViewT& mesh, int cell_id)
        {
            nb::gil_scoped_release release;
            return cutcells::make_adapt_cell(mesh, cell_id);
        },
        nb::arg("mesh"),
        nb::arg("cell_id"));

    m.def(
        "build_edges",
        [](AdaptCellT& adapt_cell)
        {
            nb::gil_scoped_release release;
            cutcells::build_edges(adapt_cell);
        },
        nb::arg("adapt_cell"));

    m.def(
        "build_faces",
        [](AdaptCellT& adapt_cell)
        {
            nb::gil_scoped_release release;
            cutcells::build_faces(adapt_cell);
        },
        nb::arg("adapt_cell"));

    m.def(
        "make_cell_level_set",
        [](const LevelSetT& ls, int cell_id)
        {
            nb::gil_scoped_release release;
            return cutcells::make_cell_level_set(ls, cell_id);
        },
        nb::arg("level_set"),
        nb::arg("cell_id"));

    m.def(
        "evaluate_bernstein",
        [](cell::type cell_type, int degree, const ndarray1<T>& coeffs, const ndarray1<T>& xi)
        {
            return cutcells::bernstein::evaluate<T>(
                cell_type,
                degree,
                std::span<const T>(coeffs.data(), static_cast<std::size_t>(coeffs.size())),
                std::span<const T>(xi.data(), static_cast<std::size_t>(xi.size())));
        },
        nb::arg("cell_type"),
        nb::arg("degree"),
        nb::arg("coeffs"),
        nb::arg("xi"));

    m.def(
        "extract_parent_edge_bernstein",
        [](cell::type parent_cell_type, int degree,
           const ndarray1<T>& parent_coeffs, int parent_local_edge_id)
        {
            std::vector<T> out;
            nb::gil_scoped_release release;
            cutcells::extract_parent_edge_bernstein<T>(
                parent_cell_type,
                degree,
                std::span<const T>(parent_coeffs.data(), static_cast<std::size_t>(parent_coeffs.size())),
                parent_local_edge_id,
                out);
            return out;
        },
        nb::arg("parent_cell_type"),
        nb::arg("degree"),
        nb::arg("parent_coeffs"),
        nb::arg("parent_local_edge_id"));

    m.def(
        "restrict_edge_bernstein_exact",
        [](cell::type parent_cell_type, int degree,
           const ndarray1<T>& parent_coeffs,
           const ndarray1<T>& xi_a,
           const ndarray1<T>& xi_b)
        {
            std::vector<T> out;
            nb::gil_scoped_release release;
            cutcells::restrict_edge_bernstein_exact<T>(
                parent_cell_type,
                degree,
                std::span<const T>(parent_coeffs.data(), static_cast<std::size_t>(parent_coeffs.size())),
                std::span<const T>(xi_a.data(), static_cast<std::size_t>(xi_a.size())),
                std::span<const T>(xi_b.data(), static_cast<std::size_t>(xi_b.size())),
                out);
            return out;
        },
        nb::arg("parent_cell_type"),
        nb::arg("degree"),
        nb::arg("parent_coeffs"),
        nb::arg("xi_a"),
        nb::arg("xi_b"));

    m.def(
        "classify_edge_roots",
        [](const ndarray1<T>& edge_coeffs,
           T zero_tol,
           T sign_tol,
           int max_depth)
        {
            T split_t = T(0);
            bool has_split = false;
            EdgeRootTag tag = cutcells::classify_edge_roots<T>(
                std::span<const T>(edge_coeffs.data(), static_cast<std::size_t>(edge_coeffs.size())),
                zero_tol,
                sign_tol,
                max_depth,
                split_t,
                has_split);
            return nb::make_tuple(tag, has_split ? nb::cast(split_t) : nb::none());
        },
        nb::arg("edge_coeffs"),
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12),
        nb::arg("max_depth") = 20);

    m.def(
        "classify_new_edges",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           T zero_tol, T sign_tol, int max_depth)
        {
            nb::gil_scoped_release release;
            cutcells::classify_new_edges(adapt_cell, ls_cell, level_set_id,
                                         zero_tol, sign_tol, max_depth);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12),
        nb::arg("max_depth") = 20);

    m.def(
        "fill_all_vertex_signs_from_level_set",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           T zero_tol)
        {
            nb::gil_scoped_release release;
            cutcells::fill_all_vertex_signs_from_level_set(
                adapt_cell, ls_cell, level_set_id, zero_tol);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("zero_tol") = T(1e-12));

    m.def(
        "classify_leaf_cells",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           T zero_tol, T sign_tol)
        {
            nb::gil_scoped_release release;
            cutcells::classify_leaf_cells(adapt_cell, ls_cell, level_set_id,
                                          zero_tol, sign_tol);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12));

    m.def(
        "classify_leaf_faces",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           T zero_tol, T sign_tol)
        {
            nb::gil_scoped_release release;
            cutcells::classify_leaf_faces(adapt_cell, ls_cell, level_set_id,
                                          zero_tol, sign_tol);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12));

    m.def(
        "process_ready_to_cut_cells",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           T zero_tol, T sign_tol, int edge_max_depth)
        {
            nb::gil_scoped_release release;
            cutcells::process_ready_to_cut_cells(
                adapt_cell, ls_cell, level_set_id,
                zero_tol, sign_tol, edge_max_depth);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12),
        nb::arg("edge_max_depth") = 20);

    m.def(
        "refine_green_on_multiple_root_edges",
        [](AdaptCellT& adapt_cell, int level_set_id)
        {
            nb::gil_scoped_release release;
            return cutcells::refine_green_on_multiple_root_edges(adapt_cell, level_set_id);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_id"));

    m.def(
        "refine_red_on_ambiguous_cells",
        [](AdaptCellT& adapt_cell, int level_set_id)
        {
            nb::gil_scoped_release release;
            return cutcells::refine_red_on_ambiguous_cells(adapt_cell, level_set_id);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_id"));

    m.def(
        "certify_refine_and_process_ready_cells",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           int max_iterations, T zero_tol, T sign_tol, int edge_max_depth)
        {
            nb::gil_scoped_release release;
            cutcells::certify_refine_and_process_ready_cells(
                adapt_cell, ls_cell, level_set_id,
                max_iterations, zero_tol, sign_tol, edge_max_depth);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("max_iterations") = 8,
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12),
        nb::arg("edge_max_depth") = 20);

    m.def(
        "certify_and_refine",
        [](AdaptCellT& adapt_cell, const LevelSetCellT& ls_cell, int level_set_id,
           int max_iterations, T zero_tol, T sign_tol, int edge_max_depth)
        {
            nb::gil_scoped_release release;
            cutcells::certify_and_refine(adapt_cell, ls_cell, level_set_id,
                                         max_iterations, zero_tol, sign_tol,
                                         edge_max_depth);
        },
        nb::arg("adapt_cell"),
        nb::arg("level_set_cell"),
        nb::arg("level_set_id"),
        nb::arg("max_iterations") = 8,
        nb::arg("zero_tol") = T(1e-12),
        nb::arg("sign_tol") = T(1e-12),
        nb::arg("edge_max_depth") = 20);

    if constexpr (std::is_same_v<T, double>)
    {
        m.attr("AdaptCell") = m.attr(adapt_name.c_str());
        m.attr("LevelSetCell") = m.attr(lsc_name.c_str());
    }
}

} // namespace

NB_MODULE(_cutcellscpp, m)
{
  // Create module for C++ wrappers
  m.doc() = "CutCells Python interface";

  nb::enum_<cell::type>(m, "CellType")
    .value("point", cell::type::point)
    .value("interval", cell::type::interval)
    .value("triangle", cell::type::triangle)
    .value("tetrahedron", cell::type::tetrahedron)
    .value("quadrilateral", cell::type::quadrilateral)
    .value("hexahedron", cell::type::hexahedron)
    .value("prism", cell::type::prism)
    .value("pyramid", cell::type::pyramid);

  nb::enum_<cutcells::EdgeRootTag>(m, "EdgeRootTag")
    .value("not_classified", cutcells::EdgeRootTag::not_classified)
    .value("no_root", cutcells::EdgeRootTag::no_root)
    .value("one_root", cutcells::EdgeRootTag::one_root)
    .value("multiple_roots", cutcells::EdgeRootTag::multiple_roots)
    .value("zero", cutcells::EdgeRootTag::zero);

  nb::enum_<cutcells::CellCertTag>(m, "CellCertTag")
    .value("not_classified", cutcells::CellCertTag::not_classified)
    .value("positive", cutcells::CellCertTag::positive)
    .value("negative", cutcells::CellCertTag::negative)
    .value("cut", cutcells::CellCertTag::cut)
    .value("zero", cutcells::CellCertTag::zero)
    .value("ambiguous", cutcells::CellCertTag::ambiguous)
    .value("ready_to_cut", cutcells::CellCertTag::ready_to_cut);

  nb::enum_<cutcells::FaceCertTag>(m, "FaceCertTag")
    .value("not_classified", cutcells::FaceCertTag::not_classified)
    .value("positive", cutcells::FaceCertTag::positive)
    .value("negative", cutcells::FaceCertTag::negative)
    .value("cut", cutcells::FaceCertTag::cut)
    .value("zero", cutcells::FaceCertTag::zero)
    .value("ambiguous", cutcells::FaceCertTag::ambiguous);

  declare_float<float>(m, "float32");
  declare_float<double>(m, "float64");

  declare_meshview_and_levelset<float>(m, "float32");
  declare_meshview_and_levelset<double>(m, "float64");
  declare_certification<float>(m, "float32");
  declare_certification<double>(m, "float64");

  declare_ho_cut<float>(m, "float32");
  declare_ho_cut<double>(m, "float64");

  m.def("csr_to_vtk_cells",
        [](const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
           const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offsets)
        {
          return as_nbarray(csr_to_vtk_cells_impl(
            std::span<const int>(connectivity.data(), connectivity.size()),
            std::span<const int>(offsets.data(), offsets.size())));
        },
        nb::arg("connectivity"),
        nb::arg("offsets"),
        "Pack CSR connectivity/offsets to VTK cells layout [n0, v0..., n1, v1..., ...].");

  m.def("write_vtk",
        [](std::string filename,
           const nb::ndarray<const double, nb::shape<-1>, nb::c_contig>& vertex_coordinates,
           const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
           const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offsets,
           const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& element_types,
           int gdim)
        {
          std::vector<cell::type> types;
          types.reserve(element_types.size());
          for (std::size_t i = 0; i < element_types.size(); ++i)
            types.push_back(static_cast<cell::type>(element_types.data()[i]));

          nb::gil_scoped_release release;
          io::write_vtk(
            std::move(filename),
            std::span<const double>(vertex_coordinates.data(), vertex_coordinates.size()),
            std::span<const int>(connectivity.data(), connectivity.size()),
            std::span<const int>(offsets.data(), offsets.size()),
            std::span<cell::type>(types.data(), types.size()),
            gdim);
        },
        nb::arg("filename"),
        nb::arg("vertex_coordinates"),
        nb::arg("connectivity"),
        nb::arg("offsets"),
        nb::arg("element_types"),
        nb::arg("gdim"),
        "Write an unstructured VTK XML file using the existing C++ writer.\n"
        "vertex_coordinates is a flat array of length num_points * gdim.");
}
