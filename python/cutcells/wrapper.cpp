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

#include "../../cpp/src/cell_types.h"
#include "../../cpp/src/cut_cell.h"
#include "../../cpp/src/cut_mesh.h"
#include "../../cpp/src/write_vtk.h"
#include "../../cpp/src/mapping.h"
#include "../../cpp/src/quadrature.h"
#include "../../cpp/src/mesh_view.h"
#include "../../cpp/src/level_set.h"

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
            return nb::ndarray<const int, nb::numpy>(
                self.cell_types.data(),
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
            return nb::ndarray<const int, nb::numpy>(
                self.cell_types.data(),
                {self.cell_types.size()},
                nb::cast(self, nb::rv_policy::reference));
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
        std::optional<ndarray1<int>> cell_types;
        if (!cell_types_obj.is_none())
          cell_types = nb::cast<ndarray1<int>>(cell_types_obj);

        std::span<const int> cell_types_span;
        if (cell_types.has_value())
        {
          cell_types_span = std::span<const int>(
              cell_types->data(), static_cast<std::size_t>(cell_types->size()));
        }

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
      [](const LevelSetMeshDataT& mesh_data, const ndarray1<T>& dof_values)
      {
        auto mesh_data_ptr = std::make_shared<LevelSetMeshDataT>(mesh_data);
        return cutcells::create_level_set_function<T, int>(
            std::move(mesh_data_ptr),
            std::span<const T>(
                dof_values.data(),
                static_cast<std::size_t>(dof_values.size())));
      },
      nb::arg("mesh_data"),
      nb::arg("dof_values"),
      "Create a polynomial LevelSetFunction from mesh_data and global dof values.");

  m.def(
      "interpolate_level_set",
      [](const MeshViewT& mesh, nb::callable phi, int degree)
      {
        LevelSetMeshDataT mesh_data;
        {
          nb::gil_scoped_release release;
          mesh_data = cutcells::create_level_set_mesh_data<T, int>(mesh, degree, T(-1));
        }

        const std::size_t num_dofs = static_cast<std::size_t>(mesh_data.num_dofs());
        const std::size_t gdim = static_cast<std::size_t>(mesh_data.gdim);
        nb::ndarray<const T, nb::numpy> x(
            mesh_data.dof_coordinates.data(),
            {num_dofs, gdim},
            nb::handle());

        // Intentional single batched callback invocation.
        nb::object values_obj = phi(x);
        auto values = nb::cast<ndarray1<T>>(values_obj);
        if (static_cast<std::size_t>(values.size()) != num_dofs)
        {
          throw std::runtime_error(
              "interpolate_level_set: callback must return a 1D array with length num_dofs");
        }

        auto mesh_data_ptr = std::make_shared<LevelSetMeshDataT>(std::move(mesh_data));
        return cutcells::create_level_set_function<T, int>(
            std::move(mesh_data_ptr),
            std::span<const T>(
                values.data(),
                static_cast<std::size_t>(values.size())));
      },
      nb::arg("mesh"),
      nb::arg("phi"),
      nb::arg("degree"),
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
