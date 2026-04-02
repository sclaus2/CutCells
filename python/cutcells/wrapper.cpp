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

#include <cutcells/cell_types.h>
#include <cutcells/cut_cell.h>
#include <cutcells/cut_mesh.h>
#include <cutcells/cut_mesh_curved.h>
#include <cutcells/write_vtk.h>
#include <cutcells/mapping.h>
#include <cutcells/quadrature.h>
#include <cutcells/quadrature_curved.h>
#include <cutcells/mesh_view.h>
#include <cutcells/level_set.h>
#include <cutcells/local_level_set.h>
#include <cutcells/iso_refine.h>
#include <cutcells/local_mesh.h>
#include <cutcells/curved_mesh.h>
#include <cutcells/mapping_curved.h>
#include <cutcells/edge_classification.h>
#include <cutcells/triangulation_repair.h>
#include <cutcells/green_split.h>
#include <cutcells/ray_split.h>
#include <cutcells/macro_split.h>
#include <cutcells/interface_split.h>

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

  nb::object coords_obj = nb::cast(coordinates);
  nb::object conn_obj = nb::cast(connectivity);
  nb::object offs_obj = nb::cast(offsets);
  nb::object types_obj;

  if (cell_types.has_value())
  {
    mesh.cell_types = std::span<const int>(cell_types->data(),
                                           static_cast<std::size_t>(cell_types->size()));
    types_obj = nb::cast(*cell_types);
  }

  mesh.owner = make_owner_from_objects(coords_obj, conn_obj, offs_obj, types_obj);
  return mesh;
}

template <typename T>
void declare_meshview_and_levelset(nb::module_& m, const std::string& suffix)
{
  using MeshViewT = cutcells::MeshView<T, int>;
  using LevelSetT = cutcells::LevelSetFunction<T, int>;
  using LocalLevelSetT = cutcells::LocalLevelSetFunction<T, int>;
  using LocalEdgeRestrictionT = cutcells::LocalLevelSetEdgeRestriction<T>;

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
            return nb::ndarray<const int, nb::numpy>(
                self.cell_types.data(),
                {self.cell_types.size()},
                nb::handle());
          },
          nb::rv_policy::reference_internal);

  const std::string edge_name = "LocalLevelSetEdgeRestriction_" + suffix;
  nb::class_<LocalEdgeRestrictionT>(m, edge_name.c_str(), "Local edge restriction")
      .def_prop_ro("degree", [](const LocalEdgeRestrictionT& self) { return self.degree; })
      .def_prop_ro(
          "coeffs",
          [](const LocalEdgeRestrictionT& self)
          {
            if (self.coeffs.empty())
              return nb::ndarray<const T, nb::numpy>(nullptr, {0}, nb::handle());
            return nb::ndarray<const T, nb::numpy>(
                self.coeffs.data(),
                {self.coeffs.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal);

  const std::string local_ls_name = "LocalLevelSetFunction_" + suffix;
  nb::class_<LocalLevelSetT>(m, local_ls_name.c_str(), "Cell-local level-set function")
      .def_prop_ro("parent_cell_id", [](const LocalLevelSetT& self) { return self.parent_cell_id; })
      .def_prop_ro("gdim", [](const LocalLevelSetT& self) { return self.gdim; })
      .def_prop_ro("tdim", [](const LocalLevelSetT& self) { return self.tdim; })
      .def_prop_ro("degree", [](const LocalLevelSetT& self) { return self.degree; })
      .def_prop_ro("backend", [](const LocalLevelSetT& self) { return self.backend; })
      .def("has_value", &LocalLevelSetT::has_value)
      .def("has_gradient", &LocalLevelSetT::has_gradient)
      .def("has_edge_restriction", &LocalLevelSetT::has_edge_restriction)
      .def(
          "value",
          [](const LocalLevelSetT& self, const ndarray1<T>& x_ref)
          {
            return self.value(x_ref.data());
          },
          nb::arg("x_ref"))
      .def(
          "grad",
          [](const LocalLevelSetT& self, const ndarray1<T>& x_ref)
          {
            auto* data = new std::vector<T>(static_cast<std::size_t>(self.tdim), T(0));
            self.grad(x_ref.data(), data->data());
            return nb::ndarray<const T, nb::numpy>(
                data->data(),
                {data->size()},
                nb::capsule(data, [](void* p) noexcept { delete static_cast<std::vector<T>*>(p); }));
          },
          nb::arg("x_ref"))
      .def("edge_restriction", &LocalLevelSetT::edge_restriction, nb::arg("edge_id"));

  const std::string ls_name = "LevelSetFunction_" + suffix;
  nb::class_<LevelSetT>(m, ls_name.c_str(), "Level-set function")
      .def(
          "__init__",
          [](LevelSetT* self,
             nb::object value_obj,
             nb::object grad_obj,
             nb::object nodal_values_obj,
             int gdim,
             int degree)
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
            self->degree = degree;
            self->kind = LevelSetT::Kind::callable;
            self->owner = make_owner_from_objects(value_obj, grad_obj, nodal_owner);
          },
          nb::arg("value") = nb::none(),
          nb::arg("grad") = nb::none(),
          nb::arg("nodal_values") = nb::none(),
          nb::arg("gdim") = 0,
          nb::arg("degree") = -1)
      .def_prop_ro("gdim", [](const LevelSetT& self) { return self.gdim; })
      .def_prop_ro("degree", [](const LevelSetT& self) { return self.degree; })
      .def_prop_ro("kind", [](const LevelSetT& self) { return self.kind; })
      .def("has_value", &LevelSetT::has_value)
      .def("has_gradient", &LevelSetT::has_gradient)
      .def("has_nodal_values", &LevelSetT::has_nodal_values)
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
      .def(
          "local",
          [](const LevelSetT& self, int cell_id)
          {
            return self.local(cell_id);
          },
          nb::arg("cell_id"))
      .def("value_at_node", &LevelSetT::value_at_node)
      .def_static(
          "from_fem",
          [](const MeshViewT& mesh,
             const ndarray1<T>& nodal_values,
             int degree)
          {
            return cutcells::make_level_set_function_from_fem<T, int>(
                mesh,
                std::span<const T>(nodal_values.data(), static_cast<std::size_t>(nodal_values.size())),
                degree,
                make_owner_from_objects(nb::cast(mesh), nb::cast(nodal_values)));
          },
          nb::arg("mesh"),
          nb::arg("nodal_values"),
          nb::arg("degree"))
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
          nb::rv_policy::reference_internal);

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
        nb::gil_scoped_release release;
        return mesh::cut_vtk_mesh<T>(
            std::span<const T>(ls_vals.data(), ls_vals.size()),
            mesh.coordinates,
            mesh.connectivity,
            mesh.offsets,
            mesh.cell_types,
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
    m.attr("MeshView") = m.attr(mesh_name.c_str());
    m.attr("LevelSetFunction") = m.attr(ls_name.c_str());
    m.attr("LocalLevelSetFunction") = m.attr(local_ls_name.c_str());
    m.attr("LocalLevelSetEdgeRestriction") = m.attr(edge_name.c_str());
  }
}

template <typename T>
void declare_float(nb::module_& m, std::string type)
{
    using LevelSetT = cutcells::LevelSetFunction<T, int>;
    m.def("classify_cell_domain", [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_values){
          cell::domain domain_id = cell::classify_cell_domain<T>(std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())});
          auto domain_str = cell_domain_to_str(domain_id);
          return domain_str;
        }
        , "classify a cell domain");

    using LocalMeshT = cutcells::LocalMesh<T>;
    const std::string lm_name = "LocalMesh_" + type;
    nb::class_<LocalMeshT>(m, lm_name.c_str(), "Cell-local mesh")
      .def(nb::init<>())
      .def_prop_ro("gdim", [](const LocalMeshT& self) { return self.gdim; })
      .def_prop_ro("tdim", [](const LocalMeshT& self) { return self.tdim; })
      .def_prop_ro("parent_cell_id", [](const LocalMeshT& self) { return self.parent_cell_id; })
      .def_prop_ro("n_level_sets", [](const LocalMeshT& self) { return self.n_level_sets; })
      .def("n_vertices", &LocalMeshT::n_vertices)
      .def("n_edges", &LocalMeshT::n_edges)
      .def("n_cells", &LocalMeshT::n_cells)
      .def("n_faces", &LocalMeshT::n_faces)
      .def_prop_ro("vertex_x",
        [](const LocalMeshT& self) {
          return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.vertex_x.data(), {self.vertex_x.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("cell_vertices",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.cell_vertices.data(), {self.cell_vertices.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("cell_offsets",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.cell_offsets.data(), {self.cell_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("cell_types",
        [](const LocalMeshT& self) {
          using type_id_t = std::underlying_type_t<cell::type>;
          return nb::ndarray<const type_id_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            reinterpret_cast<const type_id_t*>(self.cell_types.data()),
            {self.cell_types.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_state",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const uint8_t, nb::numpy>(
            self.edge_state.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_cert",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const uint8_t, nb::numpy>(
            self.edge_cert.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_root_vertex",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const int32_t, nb::numpy>(
            self.edge_root_vertex.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_root_parameter",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const T, nb::numpy>(
            self.edge_root_parameter.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_root_iterations",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const int32_t, nb::numpy>(
            self.edge_root_iterations.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_root_evaluations",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const int32_t, nb::numpy>(
            self.edge_root_evaluations.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_root_converged",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const uint8_t, nb::numpy>(
            self.edge_root_converged.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("edge_root_residual",
        [](const LocalMeshT& self) {
          const std::size_t shape[2] = {
            static_cast<std::size_t>(self.n_edges()),
            static_cast<std::size_t>(self.n_level_sets)};
          return nb::ndarray<const T, nb::numpy>(
            self.edge_root_residual.data(), 2, shape,
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("cell_domain",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.cell_domain.data(), {self.cell_domain.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      // ---- zero entities / topology caches ----
      .def_prop_ro("zero_entity_dim",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_dim.data(), {self.zero_entity_dim.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_zero_mask",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint64_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_zero_mask.data(), {self.zero_entity_zero_mask.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_vertices",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_vertices.data(), {self.zero_entity_vertices.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_offsets",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_offsets.data(), {self.zero_entity_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_parent_cell",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_parent_cell.data(), {self.zero_entity_parent_cell.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_parent_dim",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_parent_dim.data(), {self.zero_entity_parent_dim.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_parent_id",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_parent_id.data(), {self.zero_entity_parent_id.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_is_owned",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_is_owned.data(), {self.zero_entity_is_owned.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_endpoint_v0",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_endpoint_v0.data(), {self.zero_entity_endpoint_v0.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_entity_endpoint_v1",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_entity_endpoint_v1.data(), {self.zero_entity_endpoint_v1.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_face_edge_vertices",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_face_edge_vertices.data(), {self.zero_face_edge_vertices.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_face_edge_offsets",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_face_edge_offsets.data(), {self.zero_face_edge_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("curved_zero_ref_nodes",
        [](const LocalMeshT& self) {
          return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_zero_ref_nodes.data(), {self.curved_zero_ref_nodes.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("curved_zero_offsets",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_zero_offsets.data(), {self.curved_zero_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("curved_zero_converged",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_zero_converged.data(), {self.curved_zero_converged.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("curved_zero_status",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_zero_status.data(), {self.curved_zero_status.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_chain_offsets",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_chain_offsets.data(), {self.zero_chain_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_chain_entity_ids",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_chain_entity_ids.data(), {self.zero_chain_entity_ids.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_chain_entity_reversed",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_chain_entity_reversed.data(), {self.zero_chain_entity_reversed.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_patch_offsets",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_patch_offsets.data(), {self.zero_patch_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_patch_entity_ids",
        [](const LocalMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_patch_entity_ids.data(), {self.zero_patch_entity_ids.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def_prop_ro("zero_patch_entity_oriented",
        [](const LocalMeshT& self) {
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.zero_patch_entity_oriented.data(), {self.zero_patch_entity_oriented.size()},
            nb::cast(self, nb::rv_policy::reference));
        }, nb::rv_policy::reference_internal)
      .def("n_zero_entities", &LocalMeshT::n_zero_entities)
      .def("n_zero_chains", &LocalMeshT::n_zero_chains)
      .def("n_zero_patches", &LocalMeshT::n_zero_patches)
      // ---- straight interface entities ----
      .def_prop_ro("iface_vertices",
        [](const LocalMeshT& self) {
          if (self.iface_vertices.empty())
            return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.iface_vertices.data(), {self.iface_vertices.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("iface_offsets",
        [](const LocalMeshT& self) {
          if (self.iface_offsets.empty())
            return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.iface_offsets.data(), {self.iface_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("iface_level_set_id",
        [](const LocalMeshT& self) {
          if (self.iface_level_set_id.empty())
            return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.iface_level_set_id.data(), {self.iface_level_set_id.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("iface_parent_cell",
        [](const LocalMeshT& self) {
          if (self.iface_parent_cell.empty())
            return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.iface_parent_cell.data(), {self.iface_parent_cell.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def("n_iface_entities", &LocalMeshT::n_iface_entities)
      // ---- curved interface entity cache ----
      .def_prop_ro("curved_iface_ref_nodes",
        [](const LocalMeshT& self) {
          if (self.curved_iface_ref_nodes.empty())
            return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_iface_ref_nodes.data(),
            {self.curved_iface_ref_nodes.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("curved_iface_offsets",
        [](const LocalMeshT& self) {
          if (self.curved_iface_offsets.empty())
            return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_iface_offsets.data(),
            {self.curved_iface_offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("curved_iface_converged",
        [](const LocalMeshT& self) {
          if (self.curved_iface_converged.empty())
            return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_iface_converged.data(),
            {self.curved_iface_converged.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro("curved_iface_status",
        [](const LocalMeshT& self) {
          if (self.curved_iface_status.empty())
            return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
              nullptr, {std::size_t(0)}, nb::handle());
          return nb::ndarray<const uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.curved_iface_status.data(),
            {self.curved_iface_status.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      // ---- curved geometry cache ----
      .def_prop_ro("curved_geometry_order",
        [](const LocalMeshT& self) { return self.curved_geometry_order; })
      .def_prop_ro("curved_fallback_count",
        [](const LocalMeshT& self) { return self.curved_fallback_count; })
      .def_prop_ro("curved_ambiguous_count",
        [](const LocalMeshT& self) { return self.curved_ambiguous_count; })
      .def_prop_ro("curved_failed_count",
        [](const LocalMeshT& self) { return self.curved_failed_count; });

    // ---- BernsteinCell and lagrange_to_bernstein_cell ----
    using BernsteinCellT = cutcells::BernsteinCell<T>;
    const std::string bc_name = "BernsteinCell_" + type;
    nb::class_<BernsteinCellT>(m, bc_name.c_str(), "Bernstein polynomial on a reference cell")
      .def(nb::init<>())
      .def_prop_ro("cell_type", [](const BernsteinCellT& self) { return self.cell_type; })
      .def_prop_ro("degree", [](const BernsteinCellT& self) { return self.degree; })
      .def_prop_ro("tdim", [](const BernsteinCellT& self) { return self.tdim; })
      .def_prop_ro("coeffs",
        [](const BernsteinCellT& self) {
          return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.coeffs.data(), {self.coeffs.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal);

    const std::string ltb_name = "lagrange_to_bernstein_cell_" + type;
    m.def(
      ltb_name.c_str(),
      [](cell::type cell_type, int order,
         const ndarray1<T>& lagrange_values) {
        return cutcells::make_bernstein_cell<T>(
          cell_type, order,
          std::span<const T>(lagrange_values.data(), lagrange_values.size()));
      },
      nb::arg("cell_type"),
      nb::arg("order"),
      nb::arg("lagrange_values"),
      "Convert Lagrange nodal values to a BernsteinCell.");

    const std::string init_name = "init_local_mesh_from_template_" + type;
    m.def(
      init_name.c_str(),
      [](const RefinementTemplate& tpl,
         const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& parent_cell_coords,
         cell::type parent_cell_type,
         int parent_cell_id,
         int n_level_sets) {
        LocalMeshT mesh;
        cutcells::init_local_mesh_from_template<T>(
          mesh, tpl,
          std::span<const T>(parent_cell_coords.data(), parent_cell_coords.size()),
          parent_cell_type,
          parent_cell_id,
          n_level_sets);
        return mesh;
      },
      nb::arg("template"),
      nb::arg("parent_cell_coords"),
      nb::arg("parent_cell_type"),
      nb::arg("parent_cell_id"),
      nb::arg("n_level_sets") = 1,
      "Initialize and return a LocalMesh from a refinement template and one parent cell.");

    const std::string init_cell_name = "init_local_mesh_from_cell_" + type;
    m.def(
      init_cell_name.c_str(),
      [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& parent_cell_coords,
         cell::type parent_cell_type,
         int parent_cell_id,
         int n_level_sets) {
        LocalMeshT mesh;
        cutcells::init_local_mesh_from_cell<T>(
          mesh,
          std::span<const T>(parent_cell_coords.data(), parent_cell_coords.size()),
          parent_cell_type,
          parent_cell_id,
          n_level_sets);
        return mesh;
      },
      nb::arg("parent_cell_coords"),
      nb::arg("parent_cell_type"),
      nb::arg("parent_cell_id"),
      nb::arg("n_level_sets") = 1,
      "Initialize and return a LocalMesh from one parent mesh cell.");

    const std::string classify_name = "classify_edges_on_local_mesh_" + type;
    m.def(
      classify_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& level_set,
         int level_set_id,
         T tol,
         bool refine_on_uncertain,
         LocalLevelSetBackend backend) {
        return cutcells::classify_edges_with_backend<T, int>(
          mesh, level_set, backend, level_set_id, tol, refine_on_uncertain);
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("level_set_id") = 0,
      nb::arg("tol") = static_cast<T>(1e-14),
      nb::arg("refine_on_uncertain") = false,
      nb::arg("backend") = LocalLevelSetBackend::nodal_signs,
      "Classify all local-mesh edges for one level set.\n"
      "Returns True if refinement is recommended.");

    const std::string roots_name = "compute_roots_on_local_mesh_" + type;
    m.def(
      roots_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& level_set,
         int level_set_id,
         cell::edge_root::method root_method,
         T tol,
         LocalLevelSetBackend backend) {
        const LocalLevelSetBackend resolved_backend =
          (backend == LocalLevelSetBackend::nodal_signs && level_set.has_value())
            ? LocalLevelSetBackend::analytical_callbacks
            : backend;

        if (resolved_backend == LocalLevelSetBackend::bernstein)
        {
          cutcells::certify_local_mesh_polynomial<T, int>(
            mesh, level_set, level_set_id, 8, tol, true);
        }
        else
        {
          cutcells::classify_edges_with_backend<T, int>(
            mesh, level_set, resolved_backend, level_set_id, tol, false);
        }
        cutcells::compute_all_roots_with_backend<T, int>(
          mesh, level_set, resolved_backend, level_set_id, root_method, tol);
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("level_set_id") = 0,
      nb::arg("root_method") = cell::edge_root::method::linear,
      nb::arg("tol") = static_cast<T>(1e-14),
      nb::arg("backend") = LocalLevelSetBackend::nodal_signs,
      "Compute and cache roots on all one-root edges and store per-edge solver diagnostics.");

    const std::string decompose_name = "decompose_local_mesh_" + type;
    m.def(
      decompose_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& level_set,
         int level_set_id,
         cell::edge_root::method root_method,
         bool triangulate,
         T tol,
         LocalLevelSetBackend backend,
         bool repair_diagonals,
         int max_repair_depth) {
        const bool use_backend_path =
          backend == LocalLevelSetBackend::bernstein
          || level_set.has_value();

        if (use_backend_path)
        {
          const LocalLevelSetBackend resolved_backend =
            (backend == LocalLevelSetBackend::nodal_signs && level_set.has_value())
              ? LocalLevelSetBackend::analytical_callbacks
              : backend;
          const auto result = cutcells::decompose_local_mesh_with_backend<T, int>(
            mesh, level_set, resolved_backend, level_set_id, root_method, triangulate, 6, tol, false, repair_diagonals, max_repair_depth);
          if (result.invalid_cells > 0)
          {
            nb::module_ warnings = nb::module_::import_("warnings");
            nb::object runtime_warning = nb::module_::import_("builtins").attr("RuntimeWarning");
            warnings.attr("warn")(
              "LocalMesh decomposition still contains intersected subcells after the repair attempt.",
              runtime_warning);
          }
        }
        else
        {
          cutcells::decompose_local_mesh<T, int>(
            mesh, level_set, level_set_id, root_method, triangulate);
        }
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("level_set_id") = 0,
      nb::arg("root_method") = cell::edge_root::method::linear,
      nb::arg("triangulate") = true,
      nb::arg("tol") = static_cast<T>(1e-14),
      nb::arg("backend") = LocalLevelSetBackend::nodal_signs,
      nb::arg("repair_diagonals") = false,
      nb::arg("max_repair_depth") = 3,
      "Decompose a local mesh into phi<0 and phi>0 fragments using cached roots.\n"
      "When repair_diagonals=True, triangulation diagonals that cross the curved\n"
      "interface are detected and the parent cell is green-refined and the children\n"
      "are re-cut.  This repair is applied recursively up to max_repair_depth times.");

    const std::string build_iface_name = "build_interface_entities_" + type;
    const std::string build_zero_name = "build_zero_entities_" + type;
    m.def(
      build_zero_name.c_str(),
      [](LocalMeshT& mesh, int level_set_id) {
        nb::gil_scoped_release release;
        cutcells::build_zero_entities<T>(mesh, level_set_id);
      },
      nb::arg("mesh"),
      nb::arg("level_set_id") = 0,
      "Build straight zero entities after decomposition.");

    m.def(
      build_iface_name.c_str(),
      [](LocalMeshT& mesh, int level_set_id) {
        nb::gil_scoped_release release;
        cutcells::build_interface_entities<T>(mesh, level_set_id);
      },
      nb::arg("mesh"),
      nb::arg("level_set_id") = 0,
      "Build straight interface entities after decomposition.");

    const std::string green_split_name = "green_split_one_edge_" + type;
    m.def(
      green_split_name.c_str(),
      [](LocalMeshT& mesh, int edge_id, T t_sep) {
        nb::gil_scoped_release release;
        return cutcells::green_split_one_edge<T>(mesh, edge_id, t_sep);
      },
      nb::arg("mesh"),
      nb::arg("edge_id"),
      nb::arg("t_sep") = static_cast<T>(0.5),
      "Green-split all cells incident on the given edge at parameter t_sep.\n"
      "Returns True on success.");

    const std::string ray_split_name = "ray_split_one_tet_" + type;
    m.def(
      ray_split_name.c_str(),
      [](LocalMeshT& mesh, int cell_id,
         const LevelSetT& ls, int level_set_id, T tol) {
        nb::gil_scoped_release release;
        return cutcells::ray_split_one_tet<T>(
            mesh, cell_id, ls, level_set_id, tol);
      },
      nb::arg("mesh"),
      nb::arg("cell_id"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Split one intersected tet using a ray from the odd vertex to the\n"
      "opposite face centroid.  Returns True on success.");

    const std::string ray_refine_name = "ray_refine_local_mesh_" + type;
    m.def(
      ray_refine_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& ls, int level_set_id,
         int max_depth, T tol) {
        nb::gil_scoped_release release;
        return cutcells::ray_refine_local_mesh<T>(
            mesh, ls, level_set_id, max_depth, tol);
      },
      nb::arg("mesh"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("max_depth") = 5,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Iteratively ray-split intersected tets up to max_depth levels.\n"
      "Returns the total number of splits performed.");

    const std::string macro_split_name = "macro_split_topology1_tet_" + type;
    m.def(
      macro_split_name.c_str(),
      [](LocalMeshT& mesh, int cell_id,
         const LevelSetT& ls, int level_set_id, T tol) {
        nb::gil_scoped_release release;
        return cutcells::macro_split_topology1_tet<T>(
            mesh, cell_id, ls, level_set_id, tol);
      },
      nb::arg("mesh"),
      nb::arg("cell_id"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Split one topology-1 (1-vs-3) intersected tet into 6 subtets\n"
      "using the face midpoint + centroid macro pattern.");

    const std::string macro_refine_name = "macro_refine_topology1_" + type;
    m.def(
      macro_refine_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& ls, int level_set_id,
         int max_depth, T tol) {
        nb::gil_scoped_release release;
        return cutcells::macro_refine_topology1<T>(
            mesh, ls, level_set_id, max_depth, tol);
      },
      nb::arg("mesh"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("max_depth") = 5,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Iteratively apply topology-1 macro split to all 1-vs-3 tets.\n"
      "Returns the total number of splits performed.");

    const std::string iface_split_name = "interface_split_topology1_tet_" + type;
    m.def(
      iface_split_name.c_str(),
      [](LocalMeshT& mesh, int cell_id,
         const LevelSetT& ls, int level_set_id, T tol) {
        nb::gil_scoped_release release;
        return cutcells::interface_split_topology1_tet<T>(
            mesh, cell_id, ls, level_set_id, tol);
      },
      nb::arg("mesh"),
      nb::arg("cell_id"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Split one topology-1 (1-vs-3) intersected tet into 10 subtets\n"
      "using interface-adapted root vertices on crossing edges and an\n"
      "interior root on the apex-to-face-centroid ray.");

    const std::string iface_refine_name = "interface_refine_topology1_" + type;
    m.def(
      iface_refine_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& ls, int level_set_id,
         int max_depth, T tol) {
        nb::gil_scoped_release release;
        return cutcells::interface_refine_topology1<T>(
            mesh, ls, level_set_id, max_depth, tol);
      },
      nb::arg("mesh"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("max_depth") = 5,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Iteratively apply interface-adapted split to all 1-vs-3 tets.\n"
      "Returns the total number of splits performed.");

    const std::string iface_split2_name = "interface_split_topology2_tet_" + type;
    m.def(
      iface_split2_name.c_str(),
      [](LocalMeshT& mesh, int cell_id,
         const LevelSetT& ls, int level_set_id, T tol) {
        nb::gil_scoped_release release;
        return cutcells::interface_split_topology2_tet<T>(
            mesh, cell_id, ls, level_set_id, tol);
      },
      nb::arg("mesh"),
      nb::arg("cell_id"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Split one topology-2 (2-vs-2) intersected tet into 12 subtets\n"
      "using interface-adapted root vertices on all 4 crossing edges\n"
      "and an interior root vertex.");

    const std::string iface_refine_all_name = "interface_refine_" + type;
    m.def(
      iface_refine_all_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& ls, int level_set_id,
         int max_depth, T tol) {
        nb::gil_scoped_release release;
        return cutcells::interface_refine<T>(
            mesh, ls, level_set_id, max_depth, tol);
      },
      nb::arg("mesh"),
      nb::arg("ls"),
      nb::arg("level_set_id") = 0,
      nb::arg("max_depth") = 5,
      nb::arg("tol") = static_cast<T>(1e-13),
      "Iteratively apply interface-adapted split to all intersected tets.\n"
      "Handles both 1-vs-3 (10 children) and 2-vs-2 (12 children).\n"
      "Returns the total number of splits performed.");

    const std::string curve_iface_name = "curve_interface_entities_" + type;
    const std::string curve_zero_name = "curve_zero_entities_" + type;
    m.def(
      curve_zero_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& level_set,
         int level_set_id,
         int geom_order,
         T tol,
         LocalLevelSetBackend backend) {
        const LocalLevelSetBackend resolved_backend =
          (backend == LocalLevelSetBackend::nodal_signs && level_set.has_value())
            ? LocalLevelSetBackend::analytical_callbacks
            : backend;
        nb::gil_scoped_release release;
        cutcells::curve_zero_entities_with_backend<T, int>(
          mesh, level_set, resolved_backend, level_set_id, geom_order, tol);
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("level_set_id") = 0,
      nb::arg("geom_order") = 2,
      nb::arg("tol") = static_cast<T>(1e-12),
      nb::arg("backend") = LocalLevelSetBackend::nodal_signs,
      "Curve zero entities by projecting high-order nodes onto the isosurface.");

    m.def(
      curve_iface_name.c_str(),
      [](LocalMeshT& mesh,
         const LevelSetT& level_set,
         int level_set_id,
         int geom_order,
         T tol,
         LocalLevelSetBackend backend) {
        const LocalLevelSetBackend resolved_backend =
          (backend == LocalLevelSetBackend::nodal_signs && level_set.has_value())
            ? LocalLevelSetBackend::analytical_callbacks
            : backend;
        nb::gil_scoped_release release;
        cutcells::curve_interface_entities_with_backend<T, int>(
          mesh, level_set, resolved_backend, level_set_id, geom_order, tol);
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("level_set_id") = 0,
      nb::arg("geom_order") = 2,
      nb::arg("tol") = static_cast<T>(1e-12),
      nb::arg("backend") = LocalLevelSetBackend::nodal_signs,
      "Curve interface entities by projecting GL nodes onto the isosurface.");

    const std::string build_chain_name = "build_zero_chains_" + type;
    m.def(
      build_chain_name.c_str(),
      [](LocalMeshT& mesh, uint64_t zero_mask) {
        nb::gil_scoped_release release;
        cutcells::build_zero_chains<T>(mesh, zero_mask);
      },
      nb::arg("mesh"),
      nb::arg("zero_mask"),
      "Build connected chains of dim-1 zero entities for the given zero mask.");

    const std::string build_patch_name = "build_zero_patches_" + type;
    m.def(
      build_patch_name.c_str(),
      [](LocalMeshT& mesh, uint64_t zero_mask) {
        nb::gil_scoped_release release;
        cutcells::build_zero_patches<T>(mesh, zero_mask);
      },
      nb::arg("mesh"),
      nb::arg("zero_mask"),
      "Build connected patches of dim-2 zero entities for the given zero mask.");

    const std::string assign_zero_owner_name = "assign_zero_entity_ownership_" + type;
    m.def(
      assign_zero_owner_name.c_str(),
      [](nb::list meshes_list, uint64_t zero_mask) {
        std::vector<LocalMesh<T>*> ptrs;
        ptrs.reserve(meshes_list.size());
        for (auto item : meshes_list)
          ptrs.push_back(&nb::cast<LocalMesh<T>&>(item));
        nb::gil_scoped_release release;
        cutcells::assign_zero_entity_ownership<T>(
          std::span<LocalMesh<T>*>(ptrs.data(), ptrs.size()),
          zero_mask);
      },
      nb::arg("local_meshes"),
      nb::arg("zero_mask"),
      "Assign deterministic quadrature ownership for zero entities sharing a background carrier.");

    const std::string root_segment_name = "find_root_on_segment_" + type;
    m.def(
      root_segment_name.c_str(),
      [](const ndarray1<T>& p0_arr,
         const ndarray1<T>& p1_arr,
         const LevelSetT& level_set,
         cell::edge_root::method root_method,
         int cell_id,
         T level,
         int max_iter,
         T xtol,
         T ftol)
      {
        if (!level_set.has_value())
          throw std::runtime_error("find_root_on_segment: level_set.value(...) is required.");
        if (p0_arr.size() != p1_arr.size())
          throw std::invalid_argument("find_root_on_segment: p0 and p1 must have same dimension.");
        if (p0_arr.size() == 0 || p0_arr.size() > 3)
          throw std::invalid_argument("find_root_on_segment: dimension must be 1..3.");

        const std::span<const T> p0(p0_arr.data(), p0_arr.size());
        const std::span<const T> p1(p1_arr.data(), p1_arr.size());
        auto phi = [&](std::span<const T> x) -> T
        {
          return level_set.value(x.data(), static_cast<int>(cell_id));
        };

        const auto info = cell::edge_root::find_root_parameter_info<T>(
          p0, p1, phi, root_method, level, max_iter, xtol, ftol);

        std::vector<T> point(p0.size(), T(0));
        cell::edge_root::interpolate_point<T>(
          p0, p1, info.t, std::span<T>(point.data(), point.size()));

        nb::dict out;
        out["point"] = as_nbarray(std::move(point));
        out["t"] = nb::float_(info.t);
        out["residual"] = nb::float_(info.residual);
        out["iterations"] = nb::int_(info.iterations);
        out["evaluations"] = nb::int_(info.evaluations);
        out["converged"] = nb::bool_(info.converged);
        return out;
      },
      nb::arg("p0"),
      nb::arg("p1"),
      nb::arg("level_set"),
      nb::arg("root_method") = cell::edge_root::method::newton,
      nb::arg("cell_id") = -1,
      nb::arg("level") = static_cast<T>(0),
      nb::arg("max_iter") = 64,
      nb::arg("xtol") = static_cast<T>(1e-12),
      nb::arg("ftol") = static_cast<T>(1e-12),
      "Find a root on a segment p(t)=p0+t*(p1-p0), t in [0,1], and return "
      "{point, t, residual, iterations, evaluations, converged}.");

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

  m.def("cut_from_cached_roots",
        [](cell::type cell_type,
           const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& vertex_coordinates,
           const int gdim,
           const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_values,
           const std::string& cut_type_str,
           const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& edge_root_coords,
           const nb::ndarray<const uint8_t, nb::shape<-1>, nb::c_contig>& edge_has_root,
           bool triangulate) {
          cell::CutCell<T> cut_cell;
          nb::gil_scoped_release release;
          cell::cut_from_cached_roots<T>(
            cell_type,
            std::span{vertex_coordinates.data(), static_cast<unsigned long>(vertex_coordinates.size())},
            gdim,
            std::span{ls_values.data(), static_cast<unsigned long>(ls_values.size())},
            cut_type_str,
            std::span{edge_root_coords.data(), static_cast<unsigned long>(edge_root_coords.size())},
            std::span{edge_has_root.data(), static_cast<unsigned long>(edge_has_root.size())},
            cut_cell,
            triangulate);
          return cut_cell;
        },
        nb::arg("cell_type"),
        nb::arg("vertex_coordinates"),
        nb::arg("gdim"),
        nb::arg("ls_values"),
        nb::arg("cut_type_str"),
        nb::arg("edge_root_coords"),
        nb::arg("edge_has_root"),
        nb::arg("triangulate") = true,
        "Cut a cell while overriding edge-intersection coordinates from a root cache.");

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

    // ---- append_interface_quadrature_curved ----
    const std::string ai_curved_name = "append_interface_quadrature_curved_" + type;
    m.def(ai_curved_name.c_str(),
      [](const LocalMeshT& mesh, int zero_entity_id, int order,
         quadrature::QuadratureRules<T>& rules) {
          nb::gil_scoped_release release;
          quadrature::append_interface_quadrature_curved<T>(
              mesh, zero_entity_id, order, rules);
      },
      nb::arg("mesh"), nb::arg("zero_entity_id"), nb::arg("order"), nb::arg("rules"),
      "Append curved Pk quadrature for one zero entity to rules.");

    // ---- append_volume_quadrature_curved ----
    const std::string av_curved_name = "append_volume_quadrature_curved_" + type;
    m.def(av_curved_name.c_str(),
      [](const LocalMeshT& mesh, int cell_id, int level_set_id, int order,
         quadrature::QuadratureRules<T>& rules) {
          nb::gil_scoped_release release;
          int fallback_count = 0;
          quadrature::append_volume_quadrature_curved<T>(
              mesh, cell_id, level_set_id, order, rules, fallback_count);
      },
      nb::arg("mesh"), nb::arg("cell_id"), nb::arg("level_set_id"),
      nb::arg("order"), nb::arg("rules"),
      "Append curved volume quadrature for one inside sub-cell to rules.");

    // ---- make_quadrature_curved ----
    const std::string mqc_name = "make_quadrature_curved_" + type;
    m.def(mqc_name.c_str(),
      [](const LocalMeshT& mesh, int level_set_id, int order) {
          quadrature::QuadratureRules<T> vol_rules, iface_rules;
          {
              nb::gil_scoped_release release;
              quadrature::make_quadrature_curved<T>(
                  mesh, level_set_id, order, vol_rules, iface_rules);
          }
          return nb::make_tuple(std::move(vol_rules), std::move(iface_rules));
      },
      nb::arg("mesh"), nb::arg("level_set_id"), nb::arg("order"),
      "Build curved quadrature for all inside sub-cells and interface entities.\n"
      "Returns (volume_rules, interface_rules) as a tuple.");

    const std::string mqc_phys_name = "make_quadrature_curved_physical_points_" + type;
    m.def(mqc_phys_name.c_str(),
      [](const LocalMeshT& mesh, int level_set_id, int order) {
          quadrature::QuadratureRules<T> vol_rules, iface_rules;
          std::vector<T> vol_phys_points, iface_phys_points;
          {
              nb::gil_scoped_release release;
              quadrature::make_quadrature_curved<T>(
                  mesh, level_set_id, order,
                  vol_rules, vol_phys_points,
                  iface_rules, iface_phys_points);
          }
          return nb::make_tuple(
              std::move(vol_rules),
              as_nbarray(std::move(vol_phys_points)),
              std::move(iface_rules),
              as_nbarray(std::move(iface_phys_points)));
      },
      nb::arg("mesh"), nb::arg("level_set_id"), nb::arg("order"),
      "Build curved quadrature and return the matching physical quadrature points.\n"
      "Returns (volume_rules, volume_physical_points, interface_rules, interface_physical_points).");

    const std::string mqc_phys_backend_name =
      "make_quadrature_curved_physical_points_with_backend_" + type;
    m.def(mqc_phys_backend_name.c_str(),
      [](const LocalMeshT& mesh, int level_set_id, int order,
         mapping::CurvedMappingBackend volume_backend) {
          quadrature::QuadratureRules<T> vol_rules, iface_rules;
          std::vector<T> vol_phys_points, iface_phys_points;
          {
              nb::gil_scoped_release release;
              quadrature::make_quadrature_curved<T>(
                  mesh, level_set_id, order,
                  vol_rules, vol_phys_points,
                  iface_rules, iface_phys_points,
                  volume_backend);
          }
          return nb::make_tuple(
              std::move(vol_rules),
              as_nbarray(std::move(vol_phys_points)),
              std::move(iface_rules),
              as_nbarray(std::move(iface_phys_points)));
      },
      nb::arg("mesh"), nb::arg("level_set_id"), nb::arg("order"),
      nb::arg("volume_backend") = mapping::CurvedMappingBackend::collapsed,
      "Build curved quadrature with explicit volume mapping backend and return\n"
      "the matching physical quadrature points.\n"
      "Returns (volume_rules, volume_physical_points, interface_rules, interface_physical_points).");

    // ---- make_quadrature_curved with backend selection ----
    const std::string mqc_backend_name = "make_quadrature_curved_with_backend_" + type;
    m.def(mqc_backend_name.c_str(),
      [](const LocalMeshT& mesh, int level_set_id, int order,
         mapping::CurvedMappingBackend volume_backend) {
          quadrature::QuadratureRules<T> vol_rules, iface_rules;
          {
              nb::gil_scoped_release release;
              quadrature::make_quadrature_curved<T>(
                  mesh, level_set_id, order, vol_rules, iface_rules,
                  volume_backend);
          }
          return nb::make_tuple(std::move(vol_rules), std::move(iface_rules));
      },
      nb::arg("mesh"), nb::arg("level_set_id"), nb::arg("order"),
      nb::arg("volume_backend") = mapping::CurvedMappingBackend::collapsed,
      "Build curved quadrature with explicit volume mapping backend.\n"
      "Returns (volume_rules, interface_rules) as a tuple.");

    // ---- CurvedGlobalMesh binding ----
    using CurvedMeshT = mesh::CurvedGlobalMesh<T>;
    const std::string cgm_name = "CurvedGlobalMesh_" + type;

    nb::class_<CurvedMeshT>(m, cgm_name.c_str())
      .def(nb::init<>())
      .def_prop_ro("gdim", [](const CurvedMeshT& self) { return self.gdim; })
      .def_prop_ro("geom_order", [](const CurvedMeshT& self) { return self.geom_order; })
      .def_prop_ro("level_set_id", [](const CurvedMeshT& self) { return self.level_set_id; })
      .def_prop_ro("n_fallback_nodes", [](const CurvedMeshT& self) { return self.n_fallback_nodes; })
      .def("n_vertices", &CurvedMeshT::n_vertices)
      .def("n_cells", &CurvedMeshT::n_cells)
      .def(
        "write_vtu",
        [](const CurvedMeshT& self, const std::string& fname)
        {
          io::write_vtu(fname, self);
        },
        nb::arg("filename"))
      .def_prop_ro(
        "vertex_coords",
        [](const CurvedMeshT& self) {
          return nb::ndarray<const T, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.vertex_coords.data(), {self.vertex_coords.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "connectivity",
        [](const CurvedMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.connectivity.data(), {self.connectivity.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "offsets",
        [](const CurvedMeshT& self) {
          return nb::ndarray<const int32_t, nb::numpy, nb::shape<-1>, nb::c_contig>(
            self.offsets.data(), {self.offsets.size()},
            nb::cast(self, nb::rv_policy::reference));
        },
        nb::rv_policy::reference_internal);

    // ---- CurvedCutMeshResult binding ----
    using CurvedCutResultT = mesh::CurvedCutMeshResult<T>;
    const std::string ccmr_name = "CurvedCutMeshResult_" + type;
    nb::class_<CurvedCutResultT>(m, ccmr_name.c_str())
      .def(nb::init<>())
      .def_prop_ro(
        "inside",
        [](const CurvedCutResultT& self) -> const CurvedMeshT& { return self.inside; },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "interface",
        [](const CurvedCutResultT& self) -> const CurvedMeshT& { return self.interface; },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "outside",
        [](const CurvedCutResultT& self) -> const CurvedMeshT& { return self.outside; },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "inside_vis",
        [](const CurvedCutResultT& self) -> const mesh::CutMesh<T>& { return self.inside_vis; },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "interface_vis",
        [](const CurvedCutResultT& self) -> const mesh::CutMesh<T>& { return self.interface_vis; },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "interface_curved_vis",
        [](const CurvedCutResultT& self) -> const mesh::CutMesh<T>& { return self.interface_curved_vis; },
        nb::rv_policy::reference_internal)
      .def_prop_ro(
        "outside_vis",
        [](const CurvedCutResultT& self) -> const mesh::CutMesh<T>& { return self.outside_vis; },
        nb::rv_policy::reference_internal)
      .def("n_local_meshes", [](const CurvedCutResultT& self) {
          return static_cast<int>(self.local_meshes.size());
        })
      .def(
        "local_mesh",
        [](CurvedCutResultT& self, int i) -> LocalMeshT& {
          if (i < 0 || i >= static_cast<int>(self.local_meshes.size()))
            throw std::out_of_range("CurvedCutMeshResult.local_mesh: index out of range");
          return self.local_meshes[static_cast<std::size_t>(i)];
        },
        nb::arg("i"),
        nb::rv_policy::reference_internal)
      .def_prop_ro("n_parent_inside", [](const CurvedCutResultT& self) { return self.n_parent_inside; })
      .def_prop_ro("n_parent_outside", [](const CurvedCutResultT& self) { return self.n_parent_outside; })
      .def_prop_ro("n_parent_intersected", [](const CurvedCutResultT& self) { return self.n_parent_intersected; });

    const std::string ccm_cut_name = "cut_vtk_mesh_curved_" + type;
    m.def(
      ccm_cut_name.c_str(),
      [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& points,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offset,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& vtk_type,
         nb::callable value_callable,
         nb::object grad_callable_obj,
         const std::string& backend,
         int geom_order,
         int vis_subdivision,
         T tol,
         bool repair_diagonals,
         int max_repair_depth)
      {
        LevelSetT level_set;
        level_set.gdim = 3;
        level_set.degree = -1;
        level_set.kind = LevelSetT::Kind::callable;
        level_set.owner = make_owner_from_objects(value_callable, grad_callable_obj);
        level_set.value_fn = [value_callable](const T* x, int cell_id) -> T
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
        if (!grad_callable_obj.is_none())
        {
          nb::callable grad_callable = nb::cast<nb::callable>(grad_callable_obj);
          level_set.grad_fn = [grad_callable](const T* x, int cell_id, T* g) -> void
          {
            nb::gil_scoped_acquire gil;
            nb::ndarray<const T, nb::numpy> x_arr(
              x, {static_cast<std::size_t>(3)}, nb::handle());
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
            const T* grad_ptr = grad_arr.data();
            const std::size_t n = std::min<std::size_t>(3, grad_arr.size());
            for (std::size_t i = 0; i < n; ++i)
              g[i] = grad_ptr[i];
            for (std::size_t i = n; i < 3; ++i)
              g[i] = T(0);
          };
        }

        nb::gil_scoped_release release;
        return mesh::cut_vtk_mesh_curved<T>(
          std::span<const T>(points.data(), points.size()),
          std::span<const int>(connectivity.data(), connectivity.size()),
          std::span<const int>(offset.data(), offset.size()),
          std::span<const int>(vtk_type.data(), vtk_type.size()),
          level_set,
          geom_order,
          backend,
          vis_subdivision,
          tol,
          repair_diagonals,
          max_repair_depth);
      },
      nb::arg("points"),
      nb::arg("connectivity"),
      nb::arg("offset"),
      nb::arg("vtk_type"),
      nb::arg("f"),
      nb::arg("grad") = nb::none(),
      nb::arg("backend") = "bernstein",
      nb::arg("geom_order") = 4,
      nb::arg("vis_subdivision") = 3,
      nb::arg("tol") = static_cast<T>(1e-12),
      nb::arg("repair_diagonals") = false,
      nb::arg("max_repair_depth") = 3,
      "Whole-mesh curved cut pipeline using callable f(x) and optional grad(x).\n"
      "When repair_diagonals=True, the local tetrahedral decomposition enables\n"
      "the interface-split repair path for unresolved triangulation diagonals.");

    m.def(
      ccm_cut_name.c_str(),
      [](const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& points,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& connectivity,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& offset,
         const nb::ndarray<const int, nb::shape<-1>, nb::c_contig>& vtk_type,
         const LevelSetT& level_set,
         const std::string& backend,
         int geom_order,
         int vis_subdivision,
         T tol,
         bool repair_diagonals,
         int max_repair_depth)
      {
        nb::gil_scoped_release release;
        return mesh::cut_vtk_mesh_curved<T>(
          std::span<const T>(points.data(), points.size()),
          std::span<const int>(connectivity.data(), connectivity.size()),
          std::span<const int>(offset.data(), offset.size()),
          std::span<const int>(vtk_type.data(), vtk_type.size()),
          level_set,
          geom_order,
          backend,
          vis_subdivision,
          tol,
          repair_diagonals,
          max_repair_depth);
      },
      nb::arg("points"),
      nb::arg("connectivity"),
      nb::arg("offset"),
      nb::arg("vtk_type"),
      nb::arg("level_set"),
      nb::arg("backend") = "bernstein",
      nb::arg("geom_order") = 4,
      nb::arg("vis_subdivision") = 3,
      nb::arg("tol") = static_cast<T>(1e-12),
      nb::arg("repair_diagonals") = false,
      nb::arg("max_repair_depth") = 3,
      "Whole-mesh curved cut pipeline using a LevelSetFunction carrying either\n"
      "callable data or FEM nodal values.");

    // ---- assemble_curved_interface_mesh (convenience: builds cache if needed) ----
    const std::string aci_name = "assemble_curved_interface_mesh_" + type;
    m.def(
      aci_name.c_str(),
      [](nb::list meshes_list,
         const LevelSetT& level_set,
         LocalLevelSetBackend backend,
         int level_set_id,
         int geom_order,
         T tol)
      {
        std::vector<LocalMesh<T>*> ptrs;
        ptrs.reserve(meshes_list.size());
        for (auto item : meshes_list)
          ptrs.push_back(&nb::cast<LocalMesh<T>&>(item));
        const LocalLevelSetBackend resolved_backend =
          (backend == LocalLevelSetBackend::nodal_signs && level_set.has_value())
            ? LocalLevelSetBackend::analytical_callbacks
            : backend;
        nb::gil_scoped_release release;
        return mesh::assemble_curved_interface_mesh<T, int>(
          std::span<LocalMesh<T>*>(ptrs.data(), ptrs.size()),
          level_set, resolved_backend, level_set_id, geom_order, tol);
      },
      nb::arg("local_meshes"),
      nb::arg("level_set"),
      nb::arg("backend") = LocalLevelSetBackend::nodal_signs,
      nb::arg("level_set_id") = 0,
      nb::arg("geom_order") = 2,
      nb::arg("tol") = static_cast<T>(1e-12),
      "Assemble curved interface mesh. Builds interface entity cache if not done.");

    // ---- assemble_curved_interface_mesh_from_cache (no level set needed) ----
    const std::string aci_cache_name = "assemble_curved_interface_mesh_from_cache_" + type;
    m.def(
      aci_cache_name.c_str(),
      [](nb::list meshes_list, int level_set_id, int geom_order)
      {
        std::vector<const LocalMesh<T>*> ptrs;
        ptrs.reserve(meshes_list.size());
        for (auto item : meshes_list)
          ptrs.push_back(&nb::cast<const LocalMesh<T>&>(item));
        nb::gil_scoped_release release;
        return mesh::assemble_curved_interface_mesh<T>(
          std::span<const LocalMesh<T>*>(ptrs.data(), ptrs.size()),
          level_set_id, geom_order);
      },
      nb::arg("local_meshes"),
      nb::arg("level_set_id") = 0,
      nb::arg("geom_order") = 2,
      "Assemble curved interface mesh directly from the iface entity cache.");


    const std::string acv_name = "assemble_curved_volume_mesh_" + type;
    m.def(
      acv_name.c_str(),
      [](nb::list meshes_list, int level_set_id, int geom_order)
      {
        std::vector<const LocalMesh<T>*> ptrs;
        ptrs.reserve(meshes_list.size());
        for (auto item : meshes_list)
          ptrs.push_back(&nb::cast<const LocalMesh<T>&>(item));
        nb::gil_scoped_release release;
        return mesh::assemble_curved_volume_mesh<T>(
          std::span<const LocalMesh<T>*>(ptrs.data(), ptrs.size()),
          level_set_id, geom_order);
      },
      nb::arg("local_meshes"),
      nb::arg("level_set_id") = 0,
      nb::arg("geom_order") = 2,
      "Assemble curved volume mesh (phi<0 sub-cells) from decomposed LocalMeshes.");

    const std::string acv_backend_name = "assemble_curved_volume_mesh_with_backend_" + type;
    m.def(
      acv_backend_name.c_str(),
      [](nb::list meshes_list, int level_set_id, int geom_order,
         mapping::CurvedMappingBackend volume_backend, int vis_subdivision)
      {
        std::vector<const LocalMesh<T>*> ptrs;
        ptrs.reserve(meshes_list.size());
        for (auto item : meshes_list)
          ptrs.push_back(&nb::cast<const LocalMesh<T>&>(item));
        nb::gil_scoped_release release;
        return mesh::assemble_curved_volume_mesh<T>(
          std::span<const LocalMesh<T>*>(ptrs.data(), ptrs.size()),
          level_set_id, geom_order, volume_backend, vis_subdivision);
      },
      nb::arg("local_meshes"),
      nb::arg("level_set_id") = 0,
      nb::arg("geom_order") = 2,
      nb::arg("volume_backend") = mapping::CurvedMappingBackend::collapsed,
      nb::arg("vis_subdivision") = 1,
      "Assemble curved volume mesh with explicit mapping backend.\n"
      "When gordon_hall, interface-touching cells get GH-mapped geometry nodes.");

    // ---- triangulation repair bindings ----
    const std::string repair_name = "repair_cut_cell_diagonals_" + type;
    m.def(
      repair_name.c_str(),
      [](cell::CutCell<T>& cut_cell,
         cell::type parent_cell_type,
         nb::object eval_phi_py,
         int expected_sign,
         T tol,
         bool debug)
      {
        std::function<T(const T*)> eval_phi
            = [eval_phi_py](const T* x) -> T
        {
            nb::object result = eval_phi_py(
                nb::ndarray<const T, nb::numpy, nb::shape<-1>>(
                    x, {3}));
            return nb::cast<T>(result);
        };
        // GIL kept held: eval_phi calls back into Python
        return cell::repair_cut_cell_diagonals<T>(
            cut_cell, parent_cell_type, eval_phi,
            expected_sign, tol, debug);
      },
      nb::arg("cut_cell"),
      nb::arg("parent_cell_type"),
      nb::arg("eval_phi"),
      nb::arg("expected_sign"),
      nb::arg("tol") = static_cast<T>(1e-14),
      nb::arg("debug") = false,
      "Check and repair triangulation diagonals in a CutCell.\n"
      "Returns a TriangulationRepairInfo struct.");

    const std::string diag_check_name = "diagonal_crosses_interface_" + type;
    m.def(
      diag_check_name.c_str(),
      [](const cell::CutCell<T>& cut_cell,
         int va, int vb,
         nb::object eval_phi_py,
         int expected_sign,
         T tol)
      {
        std::function<T(const T*)> eval_phi
            = [eval_phi_py](const T* x) -> T
        {
            nb::object result = eval_phi_py(
                nb::ndarray<const T, nb::numpy, nb::shape<-1>>(
                    x, {3}));
            return nb::cast<T>(result);
        };
        // GIL kept held: eval_phi calls back into Python
        return cell::diagonal_crosses_interface<T>(
            cut_cell, va, vb, eval_phi, expected_sign, tol);
      },
      nb::arg("cut_cell"),
      nb::arg("va"),
      nb::arg("vb"),
      nb::arg("eval_phi"),
      nb::arg("expected_sign"),
      nb::arg("tol") = static_cast<T>(1e-14),
      "Check if a CutCell edge midpoint crosses the interface.");

    const std::string swap2d_name = "swap_diagonal_2d_" + type;
    m.def(
      swap2d_name.c_str(),
      [](cell::CutCell<T>& cut_cell, int cell_i, int cell_j) {
        return cell::swap_diagonal_2d<T>(cut_cell, cell_i, cell_j);
      },
      nb::arg("cut_cell"),
      nb::arg("cell_i"),
      nb::arg("cell_j"),
      "Swap the shared diagonal between two triangles in a CutCell.");

    if constexpr (std::is_same_v<T, double>)
    {
      m.attr("BernsteinCell") = m.attr(bc_name.c_str());
      m.attr("lagrange_to_bernstein_cell") = m.attr(ltb_name.c_str());
      m.attr("LocalMesh") = m.attr(lm_name.c_str());
      m.attr("init_local_mesh_from_template") = m.attr(init_name.c_str());
      m.attr("init_local_mesh_from_cell") = m.attr(init_cell_name.c_str());
      m.attr("classify_edges_on_local_mesh") = m.attr(classify_name.c_str());
      m.attr("compute_roots_on_local_mesh") = m.attr(roots_name.c_str());
      m.attr("decompose_local_mesh") = m.attr(decompose_name.c_str());
      m.attr("find_root_on_segment") = m.attr(root_segment_name.c_str());
      m.attr("CurvedGlobalMesh") = m.attr(cgm_name.c_str());
      m.attr("CurvedCutMeshResult") = m.attr(ccmr_name.c_str());
      m.attr("cut_vtk_mesh_curved") = m.attr(ccm_cut_name.c_str());
      m.attr("assemble_curved_interface_mesh") = m.attr(aci_name.c_str());
      m.attr("assemble_curved_interface_mesh_from_cache") = m.attr(aci_cache_name.c_str());
      m.attr("assemble_curved_volume_mesh") = m.attr(acv_name.c_str());
      m.attr("build_zero_entities") = m.attr(build_zero_name.c_str());
      m.attr("build_interface_entities") = m.attr(build_iface_name.c_str());
      m.attr("green_split_one_edge") = m.attr(green_split_name.c_str());
      m.attr("ray_split_one_tet") = m.attr(ray_split_name.c_str());
      m.attr("ray_refine_local_mesh") = m.attr(ray_refine_name.c_str());
      m.attr("macro_split_topology1_tet") = m.attr(macro_split_name.c_str());
      m.attr("macro_refine_topology1") = m.attr(macro_refine_name.c_str());
      m.attr("interface_split_topology1_tet") = m.attr(iface_split_name.c_str());
      m.attr("interface_refine_topology1") = m.attr(iface_refine_name.c_str());
      m.attr("interface_split_topology2_tet") = m.attr(iface_split2_name.c_str());
      m.attr("interface_refine") = m.attr(iface_refine_all_name.c_str());
      m.attr("curve_zero_entities") = m.attr(curve_zero_name.c_str());
      m.attr("curve_interface_entities") = m.attr(curve_iface_name.c_str());
      m.attr("build_zero_chains") = m.attr(build_chain_name.c_str());
      m.attr("build_zero_patches") = m.attr(build_patch_name.c_str());
      m.attr("assign_zero_entity_ownership") = m.attr(assign_zero_owner_name.c_str());
      m.attr("append_interface_quadrature_curved") = m.attr(ai_curved_name.c_str());
      m.attr("append_volume_quadrature_curved")    = m.attr(av_curved_name.c_str());
      m.attr("make_quadrature_curved")             = m.attr(mqc_name.c_str());
      m.attr("make_quadrature_curved_with_backend") = m.attr(mqc_backend_name.c_str());
      m.attr("make_quadrature_curved_physical_points") = m.attr(mqc_phys_name.c_str());
      m.attr("make_quadrature_curved_physical_points_with_backend")
        = m.attr(mqc_phys_backend_name.c_str());
      m.attr("assemble_curved_volume_mesh_with_backend") = m.attr(acv_backend_name.c_str());
      m.attr("repair_cut_cell_diagonals") = m.attr(repair_name.c_str());
      m.attr("diagonal_crosses_interface") = m.attr(diag_check_name.c_str());
      m.attr("swap_diagonal_2d") = m.attr(swap2d_name.c_str());
    }
}

} // namespace

NB_MODULE(_cutcellscpp, m)
{
  // Create module for C++ wrappers
  m.doc() = "CutCells Python interface";

  nb::class_<cell::TriangulationRepairInfo>(m, "TriangulationRepairInfo")
    .def(nb::init<>())
    .def_ro("checked_diagonals", &cell::TriangulationRepairInfo::checked_diagonals)
    .def_ro("invalid_diagonals", &cell::TriangulationRepairInfo::invalid_diagonals)
    .def_ro("swapped_diagonals", &cell::TriangulationRepairInfo::swapped_diagonals)
    .def_ro("unresolved", &cell::TriangulationRepairInfo::unresolved);

  nb::class_<RefinementTemplate>(m, "RefinementTemplate")
    .def_prop_ro("n_vertices", [](const RefinementTemplate& self) { return self.n_vertices; })
    .def_prop_ro("n_cells", [](const RefinementTemplate& self) { return self.n_cells; })
    .def_prop_ro("tdim", [](const RefinementTemplate& self) { return self.tdim; })
    .def_prop_ro("vertices_per_cell", [](const RefinementTemplate& self) { return self.vertices_per_cell; })
    .def_prop_ro("bg_cell_type", [](const RefinementTemplate& self) { return self.bg_cell_type; })
    .def_prop_ro("child_cell_type", [](const RefinementTemplate& self) { return self.child_cell_type; })
    .def_prop_ro(
      "vertex_parent_dim",
      [](const RefinementTemplate& self) {
        return nb::ndarray<const int, nb::numpy, nb::shape<-1>, nb::c_contig>(
          self.vertex_parent_dim.data(),
          {self.vertex_parent_dim.size()},
          nb::cast(self, nb::rv_policy::reference));
      },
      nb::rv_policy::reference_internal)
    .def_prop_ro(
      "vertex_parent_id",
      [](const RefinementTemplate& self) {
        return nb::ndarray<const int, nb::numpy, nb::shape<-1>, nb::c_contig>(
          self.vertex_parent_id.data(),
          {self.vertex_parent_id.size()},
          nb::cast(self, nb::rv_policy::reference));
      },
      nb::rv_policy::reference_internal)
    .def_prop_ro(
      "cell_connectivity",
      [](const RefinementTemplate& self) {
        return nb::ndarray<const int, nb::numpy, nb::shape<-1>, nb::c_contig>(
          self.cell_connectivity.data(),
          {self.cell_connectivity.size()},
          nb::cast(self, nb::rv_policy::reference));
      },
      nb::rv_policy::reference_internal);

  m.def(
    "iso_p1_template",
    [](cell::type cell_type, int order) -> const RefinementTemplate& {
      return cutcells::iso_p1_template(cell_type, order);
    },
    nb::arg("cell_type"), nb::arg("order"),
    nb::rv_policy::reference);

  m.def(
    "iso_p1_ref_coords",
    [](cell::type cell_type, int order) {
      auto x = cutcells::iso_p1_ref_coords(cell_type, order);
      return nb::ndarray<const double, nb::numpy>(
        x.data(), {x.size()}, nb::handle());
    },
    nb::arg("cell_type"), nb::arg("order"),
    "Return flat reference-node coordinates for iso-P1 (size = n_vertices * tdim).");

  nb::enum_<cell::type>(m, "CellType")
    .value("point", cell::type::point)
    .value("interval", cell::type::interval)
    .value("triangle", cell::type::triangle)
    .value("tetrahedron", cell::type::tetrahedron)
    .value("quadrilateral", cell::type::quadrilateral)
    .value("hexahedron", cell::type::hexahedron)
    .value("prism", cell::type::prism)
    .value("pyramid", cell::type::pyramid);

  nb::enum_<cell::edge_root::method>(m, "EdgeRootMethod")
    .value("linear", cell::edge_root::method::linear)
    .value("brent", cell::edge_root::method::brent)
    .value("itp", cell::edge_root::method::itp)
    .value("newton", cell::edge_root::method::newton);

  nb::enum_<LocalLevelSetBackend>(m, "LocalLevelSetBackend")
    .value("nodal_signs", LocalLevelSetBackend::nodal_signs)
    .value("bernstein", LocalLevelSetBackend::bernstein)
    .value("analytical_callbacks", LocalLevelSetBackend::analytical_callbacks);

  nb::enum_<mapping::CurvedMappingBackend>(m, "CurvedMappingBackend")
    .value("collapsed", mapping::CurvedMappingBackend::collapsed)
    .value("gordon_hall", mapping::CurvedMappingBackend::gordon_hall);

  nb::enum_<cutcells::LevelSetFunction<double, int>::Kind>(m, "LevelSetKind")
    .value("callable", cutcells::LevelSetFunction<double, int>::Kind::callable)
    .value("fem_nodal", cutcells::LevelSetFunction<double, int>::Kind::fem_nodal);

  nb::enum_<LocalLevelSetKind>(m, "LocalLevelSetKind")
    .value("callable", LocalLevelSetKind::callable)
    .value("bernstein", LocalLevelSetKind::bernstein)
    .value("taylor", LocalLevelSetKind::taylor);

  declare_float<float>(m, "float32");
  declare_float<double>(m, "float64");

  declare_meshview_and_levelset<float>(m, "float32");
  declare_meshview_and_levelset<double>(m, "float64");

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
}
