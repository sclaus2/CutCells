// Copyright (c) 2022 ONERA 
// Authors: Susanne Claus 
// This file is part of CutCells
//
// SPDX-License-Identifier:    MIT


#include <iostream>
#include <cmath>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/pair.h>
#include <algorithm>
#include <array>
#include <span>
#include <type_traits>
#include <stdexcept>
#include <memory>
#include <optional>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <cstdint>
#include <string_view>

#include <cutcells/bernstein.h>
#include <cutcells/cell_flags.h>
#include <cutcells/cell_topology.h>
#include <cutcells/cell_types.h>
#include <cutcells/compression/compress.h>
#include <cutcells/level_set.h>
#include <cutcells/level_set_cell.h>
#include <cutcells/lut/cut_cell.h>
#include <cutcells/lut/cut_mesh.h>
#include <cutcells/lut/iso_refine.h>
#include <cutcells/lut/triangulation.h>
#include <cutcells/mapping.h>
#include <cutcells/mesh_view.h>
#include <cutcells/quadrature.h>
#include <cutcells/quadrature_tables.h>
#include <cutcells/quadrays/adapters/shapeforest_tape.h>
#include <cutcells/quadrays/analytic.h>
#include <cutcells/quadrays/leaves.h>
#include <cutcells/quadrays/rules.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/reference_cell.h>
#include <cutcells/write_vtk.h>

namespace nb = nanobind;

using namespace cutcells;

namespace
{
template <typename T>
T default_sqrt_epsilon_tol()
{
  return std::sqrt(std::numeric_limits<T>::epsilon());
}

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

cell::TriangulationStrategy make_triangulation_strategy(
    bool triangulate, const std::string& triangulation)
{
  return triangulate ? cell::triangulation_strategy_from_string(triangulation)
                     : cell::TriangulationStrategy::none;
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
  return as_nbarray(std::forward<V>(x), shape.size(), shape.begin());
}

template <typename V>
auto as_nbarray(V&& x)
{
  const std::size_t size = x.size();
  return as_nbarray(std::forward<V>(x), {size});
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

// ============================================================================
// Analytic level sets (quadrays/analytic.h)
// ============================================================================

/// An analytic level set as Python holds it. The interface struct is shared
/// with the level-set functions and cut results made from it; whoever made it
/// keeps its context alive through the shared pointer.
struct PyAnalyticLevelSet
{
  std::shared_ptr<const quadrays::AnalyticLevelSet> phi;
};

/// The interface struct and the data its context points to, in one allocation.
template <typename Data>
struct AnalyticHolder
{
  quadrays::AnalyticLevelSet phi;
  Data data;
};

template <typename Data>
PyAnalyticLevelSet share_analytic(std::shared_ptr<AnalyticHolder<Data>> holder)
{
  const quadrays::AnalyticLevelSet* phi = &holder->phi;
  return PyAnalyticLevelSet{std::shared_ptr<const quadrays::AnalyticLevelSet>(std::move(holder), phi)};
}

/// Python objects, released with the GIL held wherever the last owner goes.
struct PythonObjects
{
  std::vector<nb::object> objects;

  PythonObjects() = default;
  PythonObjects(const PythonObjects&) = delete;
  PythonObjects& operator=(const PythonObjects&) = delete;
  ~PythonObjects()
  {
    if (!Py_IsInitialized())
    {
      for (nb::object& o : objects)
        o.release(); // the interpreter is gone
      return;
    }
    nb::gil_scoped_acquire gil;
    objects.clear();
  }
};

/// Python callables behind an AnalyticLevelSet: objects[0..3] are value,
/// gradient, box_bounds and taylor_bounds (None if absent), objects[4] is
/// numpy.ascontiguousarray.
struct PythonCallbacks
{
  PythonObjects keep;
};

nb::object numpy_copy(const double* data, std::size_t rows, std::size_t cols = 0)
{
  std::vector<double> v(data, data + rows * std::max<std::size_t>(cols, 1));
  if (cols == 0)
    return nb::cast(as_nbarray(std::move(v)));
  return nb::cast(as_nbarray(std::move(v), {rows, cols}));
}

/// The result of a Python callback as contiguous doubles: @p n of them, or
/// @p value_only of them if that is not 0. Returns 1, or 2 for value_only.
int read_doubles(const PythonCallbacks& cb, const nb::object& result, double* out, std::size_t n,
                 const char* what, std::size_t value_only = 0)
{
  const nb::object array = cb.keep.objects[4](result, "float64");
  const auto a = nb::cast<nb::ndarray<const double, nb::numpy, nb::c_contig>>(array);
  if (a.size() != n && (value_only == 0 || a.size() != value_only))
  {
    throw std::runtime_error(std::string("AnalyticLevelSet: ") + what + " must return "
                             + std::to_string(n) + " numbers"
                             + (value_only != 0 ? " (or " + std::to_string(value_only) + " for the value alone)" : ""));
  }
  std::copy(a.data(), a.data() + a.size(), out);
  return a.size() == n ? 1 : 2;
}

double python_value(const double* x, void* context)
{
  nb::gil_scoped_acquire gil;
  const auto& cb = *static_cast<const PythonCallbacks*>(context);
  try
  {
    return nb::cast<double>(cb.keep.objects[0](numpy_copy(x, 3)));
  }
  catch (const nb::python_error& e)
  {
    throw std::runtime_error(std::string("AnalyticLevelSet value: ") + e.what());
  }
}

double python_gradient(const double* x, double* grad, void* context)
{
  nb::gil_scoped_acquire gil;
  const auto& cb = *static_cast<const PythonCallbacks*>(context);
  try
  {
    read_doubles(cb, cb.keep.objects[1](numpy_copy(x, 3)), grad, 3, "gradient");
    return nb::cast<double>(cb.keep.objects[0](numpy_copy(x, 3)));
  }
  catch (const nb::python_error& e)
  {
    throw std::runtime_error(std::string("AnalyticLevelSet gradient: ") + e.what());
  }
}

int python_box_bounds(const double* lo, const double* hi, double* b, void* context)
{
  nb::gil_scoped_acquire gil;
  const auto& cb = *static_cast<const PythonCallbacks*>(context);
  try
  {
    const nb::object r = cb.keep.objects[2](numpy_copy(lo, 3), numpy_copy(hi, 3));
    if (r.is_none())
      return 0;
    return read_doubles(cb, r, b, 8, "box_bounds", 2);
  }
  catch (const nb::python_error& e)
  {
    throw std::runtime_error(std::string("AnalyticLevelSet box_bounds: ") + e.what());
  }
}

int python_taylor_bounds(const double* centre, const double* axes, int m, double* models, void* context)
{
  nb::gil_scoped_acquire gil;
  const auto& cb = *static_cast<const PythonCallbacks*>(context);
  try
  {
    const nb::object r = cb.keep.objects[3](numpy_copy(centre, 3), numpy_copy(axes, 3, m));
    if (r.is_none())
      return 0;
    return read_doubles(cb, r, models, static_cast<std::size_t>((m + 1) * (m + 2)), "taylor_bounds",
                        static_cast<std::size_t>(m + 2));
  }
  catch (const nb::python_error& e)
  {
    throw std::runtime_error(std::string("AnalyticLevelSet taylor_bounds: ") + e.what());
  }
}

/// A capsule's struct, copied, and the capsule that keeps its context alive.
struct CapsuleSource
{
  PythonObjects keep;
};

/// A ShapeForest tape and the functor that runs it.
struct TapeSource
{
  quadrays::shapeforest::Tape tape;
  quadrays::shapeforest::TapeLevelSet functor;
};

/// The sphere |x - c| - r (distance) or |x - c|^2 - r^2 as a C++ functor.
struct SphereFunctor
{
  std::array<double, 3> c = {0, 0, 0};
  double r = 0;
  bool distance = true;

  template <typename V>
  V operator()(const std::array<V, 3>& x) const
  {
    using std::sqrt;
    const V dx = x[0] - c[0], dy = x[1] - c[1], dz = x[2] - c[2];
    const V q = dx * dx + dy * dy + dz * dz;
    return distance ? sqrt(q) - r : q - r * r;
  }
};

constexpr const char* analytic_capsule_name = "cutcells.AnalyticLevelSet";

/// Point coordinates padded to 3.
std::array<double, 3> point3(const std::vector<double>& x)
{
  if (x.size() < 1 || x.size() > 3)
    throw std::invalid_argument("AnalyticLevelSet: a point has 1 to 3 coordinates");
  std::array<double, 3> p = {0, 0, 0};
  std::copy(x.data(), x.data() + x.size(), p.begin());
  return p;
}

void declare_analytic(nb::module_& m)
{
  nb::class_<PyAnalyticLevelSet>(m, "AnalyticLevelSet",
      "An analytic level set in physical coordinates for the quadrays backend: "
      "values, gradients, and bounds of both over boxes, as in algoim. Made from "
      "Python callables (slow, for experiments), from a capsule holding a "
      "cutcells.AnalyticLevelSet struct of function pointers (fast, for compiled "
      "geometry libraries), from a ShapeForest tape (cutcells.shapeforest) or by "
      "analytic_sphere. cut() accepts it.")
      .def(
          "__init__",
          [](PyAnalyticLevelSet* self, nb::callable value, nb::callable gradient, nb::callable box_bounds,
             nb::object taylor_bounds)
          {
            auto holder = std::make_shared<AnalyticHolder<PythonCallbacks>>();
            holder->data.keep.objects = {value, gradient, box_bounds, taylor_bounds,
                                         nb::module_::import_("numpy").attr("ascontiguousarray")};
            quadrays::AnalyticLevelSet& phi = holder->phi;
            phi.context = &holder->data;
            phi.value = python_value;
            phi.gradient = python_gradient;
            phi.box_bounds = python_box_bounds;
            phi.taylor_bounds = taylor_bounds.is_none() ? nullptr : python_taylor_bounds;
            new (self) PyAnalyticLevelSet(share_analytic(std::move(holder)));
          },
          nb::arg("value"), nb::arg("gradient"), nb::arg("box_bounds"), nb::arg("taylor_bounds") = nb::none(),
          "value(x) -> float and gradient(x) -> 3 numbers at a point x (3 coordinates); "
          "box_bounds(lo, hi) -> 8 numbers (phi in [b0, b1], d phi/dx_i in [b(2+2i), b(3+2i)]), "
          "2 if only the value has bounds, or None; optional taylor_bounds(centre, axes) -> "
          "(m+1, m+2) array of first-order Taylor models (alpha, beta_0..beta_{m-1}, eps) of "
          "phi and of d phi/dt_j over {centre + axes t : t in [-1, 1]^m}, axes of shape (3, m), "
          "its first row if only the value has a model, or None.")
      .def_static(
          "from_capsule",
          [](nb::capsule capsule)
          {
            const char* name = PyCapsule_GetName(capsule.ptr());
            if (name == nullptr || std::string_view(name) != analytic_capsule_name)
            {
              throw std::invalid_argument(std::string("AnalyticLevelSet.from_capsule: the capsule must be named ")
                                          + analytic_capsule_name);
            }
            const auto* src = static_cast<const quadrays::AnalyticLevelSet*>(PyCapsule_GetPointer(capsule.ptr(), name));
            if (src == nullptr || src->value == nullptr || src->gradient == nullptr || src->box_bounds == nullptr)
            {
              throw std::invalid_argument(
                  "AnalyticLevelSet.from_capsule: the struct needs value, gradient and box_bounds");
            }
            auto holder = std::make_shared<AnalyticHolder<CapsuleSource>>();
            holder->phi = *src;
            holder->data.keep.objects = {nb::borrow(capsule)};
            return share_analytic(std::move(holder));
          },
          nb::arg("capsule"),
          "An AnalyticLevelSet from a capsule named 'cutcells.AnalyticLevelSet' that points "
          "to the C struct of quadrays/analytic.h: void* context, then the function pointers "
          "value, gradient, box_bounds, taylor_bounds and hessian_bounds (the last two may be NULL). "
          "The struct is copied; "
          "the capsule is kept alive, so it may own the context.")
      .def_prop_ro(
          "capsule",
          [](const PyAnalyticLevelSet& self)
          {
            auto* keep = new std::shared_ptr<const quadrays::AnalyticLevelSet>(self.phi);
            PyObject* capsule = PyCapsule_New(
                const_cast<quadrays::AnalyticLevelSet*>(keep->get()), analytic_capsule_name,
                [](PyObject* c)
                {
                  delete static_cast<std::shared_ptr<const quadrays::AnalyticLevelSet>*>(PyCapsule_GetContext(c));
                });
            if (capsule == nullptr)
            {
              delete keep;
              throw nb::python_error();
            }
            PyCapsule_SetContext(capsule, keep);
            return nb::steal<nb::capsule>(capsule);
          },
          "A capsule named 'cutcells.AnalyticLevelSet' pointing to the C struct; it keeps "
          "this level set alive.")
      .def_prop_ro("has_taylor_bounds",
                   [](const PyAnalyticLevelSet& self) { return self.phi->taylor_bounds != nullptr; })
      .def(
          "value",
          [](const PyAnalyticLevelSet& self, const std::vector<double>& x)
          {
            const std::array<double, 3> p = point3(x);
            return self.phi->value(p.data(), self.phi->context);
          },
          nb::arg("x"), "phi at a point.")
      .def(
          "gradient",
          [](const PyAnalyticLevelSet& self, const std::vector<double>& x)
          {
            const std::array<double, 3> p = point3(x);
            std::vector<double> g(3);
            self.phi->gradient(p.data(), g.data(), self.phi->context);
            return as_nbarray(std::move(g));
          },
          nb::arg("x"), "Gradient of phi at a point.")
      .def(
          "box_bounds",
          [](const PyAnalyticLevelSet& self, const std::vector<double>& lo, const std::vector<double>& hi)
              -> nb::object
          {
            const std::array<double, 3> l = point3(lo), h = point3(hi);
            std::vector<double> b(8);
            const int status = self.phi->box_bounds(l.data(), h.data(), b.data(), self.phi->context);
            if (status == 0)
              return nb::none();
            if (status == 2)
              b.resize(2);
            return nb::cast(as_nbarray(std::move(b)));
          },
          nb::arg("lo"), nb::arg("hi"),
          "Bounds over the box [lo, hi]: phi in [b0, b1], d phi/dx_i in [b(2+2i), b(3+2i)]; "
          "only [b0, b1] if the gradient has no bound there, None if the value has none.")
      .def(
          "taylor_bounds",
          [](const PyAnalyticLevelSet& self, const nb::ndarray<const double, nb::numpy, nb::c_contig>& centre,
             const nb::ndarray<const double, nb::numpy, nb::c_contig>& axes) -> nb::object
          {
            if (centre.size() != 3 || axes.ndim() != 2 || axes.shape(0) != 3 || axes.shape(1) < 1
                || axes.shape(1) > 3)
              throw std::invalid_argument("taylor_bounds: centre has 3 entries, axes the shape (3, m), m <= 3");
            const int mm = static_cast<int>(axes.shape(1));
            std::vector<double> models(static_cast<std::size_t>((mm + 1) * (mm + 2)));
            const int status = quadrays::parallelepiped_bounds(*self.phi, centre.data(), axes.data(), mm,
                                                               models.data());
            if (status == 0)
              return nb::none();
            const std::size_t rows = status == 2 ? 1 : static_cast<std::size_t>(mm + 1);
            models.resize(rows * static_cast<std::size_t>(mm + 2));
            return nb::cast(as_nbarray(std::move(models), {rows, static_cast<std::size_t>(mm + 2)}));
          },
          nb::arg("centre"), nb::arg("axes"),
          "First-order Taylor models of phi and of d phi/dt_j over {centre + axes t : t in "
          "[-1, 1]^m}: rows (alpha, beta_0..beta_{m-1}, eps); from box_bounds if the level "
          "set has no taylor_bounds. Only the first row if the derivatives have no model "
          "there, None if the value has none.");

  m.def(
      "analytic_sphere",
      [](const std::vector<double>& centre, double radius, bool signed_distance)
      {
        auto holder = std::make_shared<AnalyticHolder<SphereFunctor>>();
        holder->data.c = point3(centre);
        holder->data.r = radius;
        holder->data.distance = signed_distance;
        holder->phi = quadrays::analytic_level_set(holder->data);
        return share_analytic(std::move(holder));
      },
      nb::arg("centre"), nb::arg("radius"), nb::arg("signed_distance") = true,
      "The sphere as a compiled analytic level set: |x - centre| - radius, or "
      "|x - centre|^2 - radius^2 with signed_distance=False.");

  m.def(
      "analytic_level_set_from_tape",
      [](const nb::ndarray<const std::uint8_t, nb::numpy, nb::shape<-1>, nb::c_contig>& op,
         const ndarray1<std::int32_t>& a, const ndarray1<std::int32_t>& b, const ndarray1<std::int32_t>& c,
         const ndarray1<std::int32_t>& out, const ndarray1<double>& imm, int n_registers, int output,
         const ndarray1<std::int32_t>& inputs, const ndarray1<std::int32_t>& extra_registers,
         const ndarray1<double>& extra_values)
      {
        auto holder = std::make_shared<AnalyticHolder<TapeSource>>();
        quadrays::shapeforest::Tape& t = holder->data.tape;
        auto copy = [](const auto& array, auto& vec) { vec.assign(array.data(), array.data() + array.size()); };
        copy(op, t.op);
        copy(a, t.a);
        copy(b, t.b);
        copy(c, t.c);
        copy(out, t.out);
        copy(imm, t.imm);
        copy(extra_registers, t.extra_registers);
        copy(extra_values, t.extra_values);
        if (inputs.size() != 3)
          throw std::invalid_argument("analytic_level_set_from_tape: inputs has the x, y and z registers");
        for (int i = 0; i < 3; ++i)
          t.inputs[i] = inputs(i);
        t.n_registers = n_registers;
        t.output = output;
        quadrays::shapeforest::prepare_tape(t);
        holder->data.functor.tape = &t;
        holder->phi = quadrays::analytic_level_set(holder->data.functor);
        return share_analytic(std::move(holder));
      },
      nb::arg("op"), nb::arg("a"), nb::arg("b"), nb::arg("c"), nb::arg("out"), nb::arg("imm"),
      nb::arg("n_registers"), nb::arg("output"), nb::arg("inputs"), nb::arg("extra_registers"),
      nb::arg("extra_values"),
      "An analytic level set from the arrays of a ShapeForest tape (one entry per "
      "instruction; inputs: the registers of x, y and z; extra_registers/extra_values: "
      "registers preset to constants). cutcells.shapeforest.analytic_level_set builds "
      "the arrays from a shape.");
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

/// Variant of csr_to_vtk_cells_impl for CutMesh data stored in basix ordering.
/// Non-simplex cells are permuted from basix to VTK vertex ordering.
static std::vector<int> csr_to_vtk_cells_basix_impl(
    std::span<const int> connectivity,
    std::span<const int> offsets,
    std::span<const cutcells::cell::type> types)
{
  std::vector<int> out;
  if (offsets.size() == 0)
    return out;
  const std::size_t ncells = offsets.size() - 1;
  out.reserve(connectivity.size() + ncells);
  for (std::size_t i = 0; i < ncells; ++i)
  {
    const int start = offsets[i];
    const int end   = offsets[i + 1];
    const int n     = end - start;
    out.push_back(n);
    const cutcells::cell::type ctype = types[i];
    if (ctype == cutcells::cell::type::point
        || ctype == cutcells::cell::type::interval
        || ctype == cutcells::cell::type::triangle
        || ctype == cutcells::cell::type::tetrahedron)
    {
      for (int j = start; j < end; ++j)
        out.push_back(connectivity[static_cast<std::size_t>(j)]);
    }
    else
    {
      const auto perm = cutcells::cell::basix_to_vtk_vertex_permutation(ctype);
      for (int j = 0; j < n; ++j)
        out.push_back(connectivity[static_cast<std::size_t>(
            start + perm[static_cast<std::size_t>(j)])]);
    }
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

inline bool part_mode_is_cut_only(std::string_view mode)
{
  if (mode == "cut_only")
    return true;
  if (mode == "full")
    return false;
  throw std::invalid_argument("mode must be 'cut_only' or 'full'");
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
            const std::span<const int> cell_dofs = self.cell_dofs_storage_span();
            return nb::ndarray<const int, nb::numpy>(
                cell_dofs.data(),
                {cell_dofs.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_offsets",
          [](const LevelSetMeshDataT& self)
          {
            const std::span<const int> cell_offsets = self.cell_offsets_span();
            return nb::ndarray<const int, nb::numpy>(
                cell_offsets.data(),
                {cell_offsets.size()},
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
          [](const LevelSetT& self) -> const LevelSetMeshDataT*
          {
            if (!self.has_mesh_data())
              return nullptr;
            return &self.mesh_data;
          },
          nb::rv_policy::reference_internal);

  m.def(
      "create_level_set_function",
      [](const LevelSetMeshDataT& mesh_data, const ndarray1<T>& dof_values,
         const std::string& name)
      {
        return cutcells::create_level_set_function<T, int>(
            mesh_data,
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

        return cutcells::create_level_set_function<T, int>(
            std::move(mesh_data),
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

  m.def(
      "create_level_set",
      [](const MeshViewT& mesh, const PyAnalyticLevelSet& phi, int degree, const std::string& name)
      {
        nb::gil_scoped_release release;
        return cutcells::create_level_set_function<T, int>(mesh, phi.phi, degree, name);
      },
      nb::arg("mesh"),
      nb::arg("phi"),
      nb::arg("degree"),
      nb::arg("name") = "phi",
      "A level set from an AnalyticLevelSet. cut() classifies tetrahedra and "
      "hexahedra by its own bounds and backend='quadrays' integrates it; its "
      "interpolant of the given degree feeds the straight and algoim backends.");

  m.def(
      "interpolate_level_set",
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

        // Expose as (num_dofs, gdim), i.e. one point per row (matches tests + docs).
        nb::ndarray<const T, nb::numpy> X(
            mesh_data.dof_coordinates.data(),
            {num_dofs, gdim},
            nb::handle());

        // Intentional single batched callback invocation.
        nb::object values_obj = phi(X);
        auto values = nb::cast<ndarray1<T>>(values_obj);
        if (static_cast<std::size_t>(values.size()) != num_dofs)
        {
          throw std::runtime_error(
              "interpolate_level_set: callback must return a 1D array with length num_dofs");
        }

        return cutcells::create_level_set_function<T, int>(
            std::move(mesh_data),
            std::span<const T>(
                values.data(),
                static_cast<std::size_t>(values.size())),
            name);
      },
      nb::arg("mesh"),
      nb::arg("phi"),
      nb::arg("degree"),
      nb::arg("name") = "phi",
      "Interpolate a batched callable phi(X) on higher-order level-set DOF coordinates.\n"
      "X is passed as an ndarray of shape (num_dofs, gdim).");

  // ---- cut_mesh_view ----
  // Cuts a MeshView using a LevelSetFunction, returning a CutMesh.
  // Nodal level-set values are taken from ls.nodal_values if available,
  // otherwise ls.value() is evaluated at every mesh node.
  m.def(
      "cut_mesh_view",
      [](const MeshViewT& mesh, const LevelSetT& ls,
         const std::string& cut_type_str, bool triangulate,
         const std::string& triangulation)
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
        vtk_types_vec.reserve(static_cast<std::size_t>(mesh.num_cells()));
        for (int c = 0; c < mesh.num_cells(); ++c)
        {
          const auto ct = mesh.cell_type(c);
          vtk_types_vec.push_back(
              static_cast<int>(cutcells::cell::map_cell_type_to_vtk(ct)));
        }

        const auto strategy = make_triangulation_strategy(triangulate, triangulation);
        nb::gil_scoped_release release;
        return mesh::cut_vtk_mesh<T>(
            std::span<const T>(ls_vals.data(), ls_vals.size()),
            mesh.coordinates,
            mesh.connectivity,
            mesh.offsets,
            std::span<const int>(vtk_types_vec.data(), vtk_types_vec.size()),
            cut_type_str,
            strategy);
      },
      nb::arg("mesh"),
      nb::arg("level_set"),
      nb::arg("cut_type"),
      nb::arg("triangulate") = true,
      nb::arg("triangulation") = "classical",
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
            return as_nbarray(csr_to_vtk_cells_basix_impl(
              std::span<const int>(self._connectivity.data(), self._connectivity.size()),
              std::span<const int>(self._offset.data(), self._offset.size()),
              std::span<const cell::type>(self._types.data(), self._types.size())));
          },
          nb::rv_policy::move,
          "Packed cells view [n, v0, ...] with basix-to-VTK vertex permutation applied."
          " Suitable for pv.UnstructuredGrid. Allocates on each call.")
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
                   bool triangulate,
                   const std::string& triangulation){
              cell::CutCell<T> cut_cell;
              const auto strategy = make_triangulation_strategy(triangulate, triangulation);
              nb::gil_scoped_release release;
              cell::cut<T>(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, cut_cell, strategy);
              return cut_cell;
             }
             , nb::arg("cell_type"), nb::arg("vertex_coordinates"), nb::arg("gdim"),
             nb::arg("ls_values"), nb::arg("cut_type_str"),
             nb::arg("triangulate") = false,
             nb::arg("triangulation") = "classical",
             "cut a cell");

  m.def("higher_order_cut", [](cell::type cell_type,
             const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& vertex_coordinates,
             const int gdim,
             const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& ls_values,
             const std::string& cut_type_str,
             bool triangulate,
             const std::string& triangulation){
              const auto strategy = make_triangulation_strategy(triangulate, triangulation);
              nb::gil_scoped_release release;
              cell::CutCell<T> cut_cell = cell::higher_order_cut<T>(cell_type, std::span{vertex_coordinates.data(),static_cast<unsigned long>(vertex_coordinates.size())}, gdim, std::span{ls_values.data(),static_cast<unsigned long>(ls_values.size())}, cut_type_str, strategy);
              return cut_cell;
             }
             , nb::arg("cell_type"), nb::arg("vertex_coordinates"), nb::arg("gdim"),
             nb::arg("ls_values"), nb::arg("cut_type_str"),
             nb::arg("triangulate") = false,
             nb::arg("triangulation") = "classical",
             "cut a second order cell");

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
                             bool triangulate,
                             const std::string& triangulation){
              const auto strategy = make_triangulation_strategy(triangulate, triangulation);
              nb::gil_scoped_release release;
              return  mesh::cut_vtk_mesh<T>(std::span(ls_vals.data(),ls_vals.size()),
                            std::span(points.data(),points.size()),
                            std::span(connectivity.data(),connectivity.size()),
                            std::span(offset.data(),offset.size()),
                            std::span(vtk_type.data(),vtk_type.size()),
                            cut_type_str,
                            strategy);
             }
             , nb::arg("ls_vals"), nb::arg("points"), nb::arg("connectivity"), nb::arg("offset"), nb::arg("vtk_type"),
               nb::arg("cut_type_str"), nb::arg("triangulate") = true,
               nb::arg("triangulation") = "classical"
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

template <typename T>
void declare_level_set_cell(nb::module_& m, const std::string& suffix)
{
    using LevelSetT = cutcells::LevelSetFunction<T, int>;
    using LevelSetCellT = cutcells::LevelSetCell<T, int>;

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

    if constexpr (std::is_same_v<T, double>)
        m.attr("LevelSetCell") = m.attr(lsc_name.c_str());
}

template <typename T>
void declare_write_vtk(nb::module_& m)
{
  m.def("write_vtk",
        [](std::string filename,
           const nb::ndarray<const T, nb::shape<-1>, nb::c_contig>& vertex_coordinates,
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
            std::span<const T>(vertex_coordinates.data(), vertex_coordinates.size()),
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

} // namespace

template <typename T>
void declare_quadrays(nb::module_& m, const std::string& type)
{
  namespace qr = cutcells::quadrays;
  using LeafMeshT = qr::LeafMesh<T>;

  std::string leaf_name = "QuadraysLeafMesh_" + type;
  nb::class_<LeafMeshT>(m, leaf_name.c_str(),
      "Cells for visualisation: the leaves of the quadrays decomposition as VTK "
      "Lagrange cells, and whole cells as linear VTK cells (CSR layout).")
      .def(nb::init<>())
      .def_prop_ro("points",
          [](const LeafMeshT& self) {
            return nb::ndarray<const T, nb::numpy>(
                self.points.data(), {self.points.size() / 3, 3},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal, "Physical node coordinates, shape (n_points, 3).")
      .def_prop_ro("connectivity",
          [](const LeafMeshT& self) {
            return nb::ndarray<const std::int32_t, nb::numpy>(
                self.connectivity.data(), {self.connectivity.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal, "Node indices of all cells, in VTK's node order.")
      .def_prop_ro("offsets",
          [](const LeafMeshT& self) {
            return nb::ndarray<const std::int32_t, nb::numpy>(
                self.offsets.data(), {self.offsets.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal, "CSR offsets, shape (n_cells + 1,).")
      .def_prop_ro("vtk_types",
          [](const LeafMeshT& self) {
            return nb::ndarray<const std::uint8_t, nb::numpy>(
                self.vtk_types.data(), {self.vtk_types.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal, "VTK cell type per cell.")
      .def_prop_ro("parent",
          [](const LeafMeshT& self) {
            return nb::ndarray<const std::int32_t, nb::numpy>(
                self.parent.data(), {self.parent.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal, "Background cell per cell.")
      .def_prop_ro("degree",
          [](const LeafMeshT& self) {
            return nb::ndarray<const std::int32_t, nb::numpy>(
                self.degree.data(), {self.degree.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal, "Polynomial degree per cell (1 for linear cells).")
      .def("n_points", &LeafMeshT::n_points)
      .def("n_cells", &LeafMeshT::n_cells);

  // one cell given by its type, vertices and level-set Bernstein coefficients
  auto cell_inputs = [](cell::type cell_type, const nb::ndarray<const T, nb::numpy, nb::c_contig>& vertex_coords,
                        int degree, const ndarray1<T>& coeffs, qr::ClippedBox<T>& box,
                        qr::BoxBernstein<T>& phi)
  {
    qr::make_clipped_box<T>(cell_type, std::span<const T>(vertex_coords.data(), vertex_coords.size()),
                            cell::get_tdim(cell_type), box);
    qr::cell_bernstein_on_box<T>(cell_type, degree, std::span<const T>(coeffs.data(), coeffs.size()), phi);
  };
  auto cell_part = [](const std::string& selection, const std::string& name)
  {
    cutcells::SelectionExpr expr = cutcells::parse_selection_expr(selection);
    cutcells::compile_selection_expr(expr, {name});
    if (expr.terms.size() != 1)
      throw std::runtime_error("quadrays: one selection term per call");
    return qr::part_of(expr.terms.front(), 0);
  };

  m.def(("quadrays_cell_rules_" + type).c_str(),
        [cell_inputs, cell_part](cell::type cell_type,
                                 const nb::ndarray<const T, nb::numpy, nb::c_contig>& vertex_coords, int degree,
                                 const ndarray1<T>& bernstein_coeffs, const std::string& selection, int q,
                                 const qr::Options& options, const std::string& level_set_name)
        {
          qr::ClippedBox<T> box;
          qr::BoxBernstein<T> phi;
          cell_inputs(cell_type, vertex_coords, degree, bernstein_coeffs, box, phi);
          const qr::Part part = cell_part(selection, level_set_name);
          quadrature::QuadratureRules<T> rules;
          qr::Stats stats;
          {
            nb::gil_scoped_release release;
            qr::append_rules(box, phi, part, q, options, 0, rules, stats);
          }
          return std::make_pair(std::move(rules), std::move(stats));
        },
        nb::arg("cell_type"), nb::arg("vertex_coords"), nb::arg("degree"), nb::arg("bernstein_coeffs"),
        nb::arg("selection"), nb::arg("q") = 3, nb::arg("options") = qr::Options{},
        nb::arg("level_set_name") = "phi",
        "Quadrature rule of one part of one cell (triangle, quadrilateral, tetrahedron, "
        "hexahedron, prism or pyramid) from its vertices (Basix order, tdim coordinates each) and the Bernstein "
        "coefficients of its level set (CutCells' order). Returns (QuadratureRules, QuadraysStats).");

  m.def(("quadrays_cell_leaves_" + type).c_str(),
        [cell_inputs, cell_part](cell::type cell_type,
                                 const nb::ndarray<const T, nb::numpy, nb::c_contig>& vertex_coords, int degree,
                                 const ndarray1<T>& bernstein_coeffs, const std::string& selection,
                                 int leaf_degree, const qr::Options& options, const std::string& level_set_name)
        {
          qr::ClippedBox<T> box;
          qr::BoxBernstein<T> phi;
          cell_inputs(cell_type, vertex_coords, degree, bernstein_coeffs, box, phi);
          const qr::Part part = cell_part(selection, level_set_name);
          LeafMeshT leaves;
          qr::Stats stats;
          {
            nb::gil_scoped_release release;
            qr::append_leaves(box, phi, part, leaf_degree, options, 0, leaves, stats);
          }
          return leaves;
        },
        nb::arg("cell_type"), nb::arg("vertex_coords"), nb::arg("degree"), nb::arg("bernstein_coeffs"),
        nb::arg("selection"), nb::arg("leaf_degree") = 3, nb::arg("options") = qr::Options{},
        nb::arg("level_set_name") = "phi",
        "Leaf cells of one part of one cell, as quadrays_cell_rules takes it.");

  m.def(("quadrays_cell_rules_" + type).c_str(),
        [cell_part](cell::type cell_type, const nb::ndarray<const T, nb::numpy, nb::c_contig>& vertex_coords,
                    const PyAnalyticLevelSet& level_set, const std::string& selection, int q,
                    const qr::Options& options, const std::string& level_set_name)
        {
          qr::ClippedBox<T> box;
          qr::make_clipped_box<T>(cell_type, std::span<const T>(vertex_coords.data(), vertex_coords.size()),
                                  cell::get_tdim(cell_type), box);
          const qr::Part part = cell_part(selection, level_set_name);
          quadrature::QuadratureRules<T> rules;
          qr::Stats stats;
          {
            nb::gil_scoped_release release;
            qr::append_rules(box, qr::analytic_source(*level_set.phi, box), part, q, options, 0, rules, stats);
          }
          return std::make_pair(std::move(rules), std::move(stats));
        },
        nb::arg("cell_type"), nb::arg("vertex_coords"), nb::arg("level_set"), nb::arg("selection"),
        nb::arg("q") = 3, nb::arg("options") = qr::Options{}, nb::arg("level_set_name") = "phi",
        "Quadrature rule of one part of one cell for an AnalyticLevelSet; any cell of quadrays "
        "(vertices with tdim coordinates each). Returns (QuadratureRules, QuadraysStats).");

  m.def(("quadrays_cell_leaves_" + type).c_str(),
        [cell_part](cell::type cell_type, const nb::ndarray<const T, nb::numpy, nb::c_contig>& vertex_coords,
                    const PyAnalyticLevelSet& level_set, const std::string& selection, int leaf_degree,
                    const qr::Options& options, const std::string& level_set_name)
        {
          qr::ClippedBox<T> box;
          qr::make_clipped_box<T>(cell_type, std::span<const T>(vertex_coords.data(), vertex_coords.size()),
                                  cell::get_tdim(cell_type), box);
          const qr::Part part = cell_part(selection, level_set_name);
          LeafMeshT leaves;
          qr::Stats stats;
          {
            nb::gil_scoped_release release;
            qr::append_leaves(box, qr::analytic_source(*level_set.phi, box), part, leaf_degree, options, 0,
                              leaves, stats);
          }
          return leaves;
        },
        nb::arg("cell_type"), nb::arg("vertex_coords"), nb::arg("level_set"), nb::arg("selection"),
        nb::arg("leaf_degree") = 3, nb::arg("options") = qr::Options{}, nb::arg("level_set_name") = "phi",
        "Leaf cells of one part of one cell for an AnalyticLevelSet.");

  m.def(("write_quadrays_leaves_" + type).c_str(),
        [](const std::string& filename, const LeafMeshT& leaves)
        {
          nb::gil_scoped_release release;
          qr::write_leaves(filename, leaves);
        },
        nb::arg("filename"), nb::arg("leaves"),
        "Write leaf cells to a .vtu file (VTK 9.1 node order, with HigherOrderDegrees "
        "and the cell data parent_id).");

  if constexpr (std::is_same_v<T, double>)
  {
    m.attr("QuadraysLeafMesh") = m.attr(leaf_name.c_str());
    m.attr("quadrays_cell_rules") = m.attr("quadrays_cell_rules_float64");
    m.attr("quadrays_cell_leaves") = m.attr("quadrays_cell_leaves_float64");
    m.attr("write_quadrays_leaves") = m.attr("write_quadrays_leaves_float64");
  }
}

// ============================================================================
// The front end: cutcells.cut, and the submodule cutcells.part
// ============================================================================

/// What Python holds for a cut: the mesh and level sets the result points to,
/// at stable addresses, and the backend its parts use unless a call names one,
/// with that backend's options (None: its defaults).
template <typename T>
struct PartCutResult
{
  std::shared_ptr<const cutcells::MeshView<T, int>> mesh;
  std::shared_ptr<const std::vector<cutcells::LevelSetFunction<T, int>>> level_sets;
  cutcells::part::CutResult<T, int> result;
  std::string backend = "quadrays";
  nb::object options = nb::none();
};

/// A part as Python holds it, with the backend and options of its result.
template <typename T>
struct PartSelection
{
  cutcells::part::MeshPart<T, int> part;
  std::string backend;
  nb::object options;
};

/// A backend's name; "straight" is the lookup tables' old name.
inline std::string part_backend(const std::string& backend)
{
  if (backend == "lut" || backend == "straight")
    return "lut";
  if (backend == "quadrays")
    return backend;
  if (backend == "algoim" || backend == "algoim_general")
    throw std::invalid_argument("part: the algoim backends moved to benchmarks/ (CUTCELLS_WITH_ALGOIM there)");
  throw std::invalid_argument("part: unknown backend '" + backend + "'; expected 'quadrays' or 'lut'");
}

/// The backend of a call: the one named, else the part's.
inline std::string call_backend(nb::handle backend, const std::string& part_default)
{
  return part_backend(backend.is_none() ? part_default : nb::cast<std::string>(backend));
}

/// The options of a call: the ones given, else the part's if they belong to
/// the backend, else the backend's defaults.
template <typename Options>
Options call_options(nb::handle given, nb::handle part_options, const std::string& backend, const char* type_name)
{
  if (!given.is_none())
  {
    if (!nb::isinstance<Options>(given))
      throw nb::type_error(("part: the backend '" + backend + "' takes " + type_name).c_str());
    return nb::cast<Options>(given);
  }
  if (part_options.is_valid() && nb::isinstance<Options>(part_options))
    return nb::cast<Options>(part_options);
  return Options{};
}

/// Checks options against a backend; None passes.
inline nb::object checked_options(nb::handle options, const std::string& backend)
{
  if (options.is_none())
    return nb::none();
  const bool lut = backend == "lut";
  if (lut ? !nb::isinstance<cutcells::lut::Options>(options) : !nb::isinstance<cutcells::quadrays::Options>(options))
    throw nb::type_error(("part: the backend '" + backend + "' takes " + (lut ? "LutOptions" : "QuadraysOptions")).c_str());
  return nb::borrow(options);
}

/// The level sets of a cut, a LevelSetFunction or AnalyticLevelSet or a list
/// of them, classified on the mesh.
template <typename T>
PartCutResult<T> cut_mesh_part(const cutcells::MeshView<T, int>& mesh, nb::handle level_sets, nb::handle names,
                               int max_depth)
{
  using MeshViewT = cutcells::MeshView<T, int>;
  using LevelSetT = cutcells::LevelSetFunction<T, int>;
  std::vector<nb::handle> items;
  if (nb::isinstance<nb::list>(level_sets) || nb::isinstance<nb::tuple>(level_sets))
  {
    for (nb::handle h : level_sets)
      items.push_back(h);
  }
  else
    items.push_back(level_sets);
  std::vector<std::string> given;
  if (!names.is_none())
    given = nb::cast<std::vector<std::string>>(names);
  if (!given.empty() && given.size() != items.size())
    throw std::invalid_argument("cut: give one name per level set");

  auto owned = std::make_shared<std::vector<LevelSetT>>();
  for (std::size_t i = 0; i < items.size(); ++i)
  {
    const std::string fallback = items.size() == 1 ? "phi" : "phi" + std::to_string(i + 1);
    if (nb::isinstance<LevelSetT>(items[i]))
    {
      owned->push_back(nb::cast<const LevelSetT&>(items[i]));
      if (!given.empty())
        owned->back().name = given[i];
    }
    else if (nb::isinstance<PyAnalyticLevelSet>(items[i]))
    {
      const PyAnalyticLevelSet& phi = nb::cast<const PyAnalyticLevelSet&>(items[i]);
      owned->push_back(
          cutcells::create_level_set_function<T, int>(phi.phi, mesh.gdim, given.empty() ? fallback : given[i]));
    }
    else
      throw nb::type_error("cut: level sets are LevelSetFunctions or AnalyticLevelSets");
  }
  for (std::size_t i = 0; i < owned->size(); ++i)
    for (std::size_t j = 0; j < i; ++j)
      if ((*owned)[i].name == (*owned)[j].name)
        throw std::invalid_argument("cut: two level sets are named '" + (*owned)[i].name + "'");

  PartCutResult<T> r;
  r.mesh = std::make_shared<const MeshViewT>(mesh);
  r.level_sets = owned;
  cutcells::part::ClassifyOptions options;
  options.max_depth = max_depth;
  {
    nb::gil_scoped_release release;
    r.result = cutcells::part::cut<T, int>(*r.mesh, std::span<const LevelSetT>(*r.level_sets), options);
  }
  return r;
}

/// The lookup tables' options from the keywords of the old cut().
inline cutcells::lut::Options legacy_lut_options(bool triangulate, const std::string& triangulation,
                                                 const std::string& cut_approximation, int cut_approximation_order,
                                                 nb::handle degree)
{
  cutcells::lut::Options o;
  o.triangulate = triangulate;
  o.triangulation = cell::triangulation_strategy_from_string(triangulation);
  if (cut_approximation == "auto")
    o.template_order = degree.is_none() ? 0 : nb::cast<int>(degree);
  else if (cut_approximation == "linear")
  {
    if (cut_approximation_order != 1)
      throw std::invalid_argument("cut: cut_approximation='linear' requires cut_approximation_order=1");
    o.template_order = 1;
  }
  else if (cut_approximation == "iso_p1")
    o.template_order = cut_approximation_order;
  else
    throw std::invalid_argument("cut: cut_approximation must be 'auto', 'linear', or 'iso_p1'");
  if (o.template_order < 0 || o.template_order > 4)
    throw std::invalid_argument("cut: the template order (cut_approximation_order, degree) goes from 1 to 4");
  return o;
}

template <typename T>
void declare_part(nb::module_& m, nb::module_& part_module, const std::string& type)
{
  namespace qr = cutcells::quadrays;
  using MeshViewT = cutcells::MeshView<T, int>;
  using ResultT = PartCutResult<T>;
  using PartT = PartSelection<T>;

  const std::string result_name = "CutResult_" + type;
  nb::class_<ResultT>(part_module, result_name.c_str(),
      "Every cell classified by every level set as inside, outside or cut, by "
      "the level sets' own bounds, and the faces lying in a zero set with the "
      "cell that owns each. result[\"phi1 < 0 and phi2 = 0\"] selects a MeshPart. "
      "Its parts integrate with result.backend ('quadrays' or 'lut') and "
      "result.options unless a call names others.")
      .def_prop_ro("level_set_names", [](const ResultT& self) { return self.result.level_set_names; })
      .def_prop_ro("num_cells", [](const ResultT& self) { return self.result.num_cells; })
      .def_prop_ro("num_level_sets", [](const ResultT& self) { return self.result.n_level_sets(); })
      .def_prop_ro("num_cut_cells",
                   [](const ResultT& self) { return static_cast<int>(self.result.cut_cells.size()); })
      .def_prop_ro(
          "cut_cells",
          [](const ResultT& self)
          {
            return nb::ndarray<const int, nb::numpy>(self.result.cut_cells.data(), {self.result.cut_cells.size()},
                                                     nb::handle());
          },
          nb::rv_policy::reference_internal, "Cells that some level set cuts, ascending.")
      .def_prop_ro(
          "parent_cell_ids",
          [](const ResultT& self)
          {
            return nb::ndarray<const int, nb::numpy>(self.result.cut_cells.data(), {self.result.cut_cells.size()},
                                                     nb::handle());
          },
          nb::rv_policy::reference_internal, "The cut cells (the name of HOCutResult).")
      .def_prop_ro(
          "domains",
          [](const ResultT& self)
          {
            std::vector<std::int8_t> d(self.result.domains.size());
            for (std::size_t i = 0; i < d.size(); ++i)
              d[i] = static_cast<std::int8_t>(self.result.domains[i]);
            return as_nbarray(std::move(d), {static_cast<std::size_t>(self.result.n_level_sets()),
                                             static_cast<std::size_t>(self.result.num_cells)});
          },
          nb::rv_policy::move, "Per level set and cell: 0 inside, 1 cut, 2 outside; shape (num_level_sets, num_cells).")
      .def_prop_ro(
          "cell_domains",
          [](const ResultT& self)
          {
            std::vector<int> d(self.result.domains.size());
            for (std::size_t i = 0; i < d.size(); ++i)
              d[i] = static_cast<int>(self.result.domains[i]);
            return as_nbarray(std::move(d), {static_cast<std::size_t>(self.result.n_level_sets()),
                                             static_cast<std::size_t>(self.result.num_cells)});
          },
          nb::rv_policy::move, "domains as int (the name of HOCutResult).")
      .def_prop_ro(
          "zero_faces",
          [](const ResultT& self)
          {
            const std::size_t n = static_cast<std::size_t>(self.result.n_zero_faces());
            std::vector<int> z(3 * n);
            for (std::size_t i = 0; i < n; ++i)
            {
              z[3 * i] = self.result.zero_face_level_sets[i];
              z[3 * i + 1] = self.result.zero_face_cells[i];
              z[3 * i + 2] = self.result.zero_face_local[i];
            }
            return as_nbarray(std::move(z), {n, std::size_t(3)});
          },
          nb::rv_policy::move,
          "Faces in a zero set, each once: (level set, owning cell, face of the cell in "
          "Basix numbering), shape (n, 3). The owner is the cell on the negative side, "
          "else the lower cell index.")
      .def_prop_rw(
          "backend", [](const ResultT& self) { return self.backend; },
          [](ResultT& self, const std::string& backend)
          {
            const std::string b = part_backend(backend);
            self.options = b == self.backend ? self.options : nb::none();
            self.backend = b;
          },
          "The backend of parts selected from now on: 'quadrays' or 'lut' ('straight').")
      .def_prop_rw(
          "options", [](const ResultT& self) { return self.options; },
          [](ResultT& self, nb::handle options) { self.options = checked_options(options, self.backend); },
          "Options of the backend for parts selected from now on (None: its defaults).")
      .def(
          "__getitem__",
          [](const ResultT& self, const std::string& expr)
          { return PartT{cutcells::part::select(self.result, std::string_view(expr)), self.backend, self.options}; },
          nb::arg("expr"), nb::keep_alive<0, 1>(),
          "The MeshPart that a selection expression such as \"phi < 0\" selects.");

  const std::string part_name = "MeshPart_" + type;
  nb::class_<PartT>(part_module, part_name.c_str(),
      "A part of the mesh: the cells wholly in it, the cut cells holding a piece "
      "of it, and the zero faces in it. Its quadrature and visualisation come from "
      "a backend: its result's unless a call names one.")
      .def_prop_ro("dim", [](const PartT& self) { return self.part.dim; })
      .def_prop_ro("num_cut_cells", [](const PartT& self) { return self.part.n_cut_cells(); })
      .def_prop_ro("num_uncut_cells", [](const PartT& self) { return self.part.n_uncut_cells(); })
      .def_prop_ro("backend", [](const PartT& self) { return self.backend; })
      .def_prop_ro("options", [](const PartT& self) { return self.options; })
      .def_prop_ro(
          "cut_cells",
          [](const PartT& self)
          {
            return nb::ndarray<const int, nb::numpy>(self.part.cut_cells.data(), {self.part.cut_cells.size()},
                                                     nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "uncut_cells",
          [](const PartT& self)
          {
            return nb::ndarray<const int, nb::numpy>(self.part.uncut_cells.data(), {self.part.uncut_cells.size()},
                                                     nb::handle());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "uncut_cell_ids",
          [](const PartT& self)
          {
            return nb::ndarray<const int, nb::numpy>(self.part.uncut_cells.data(), {self.part.uncut_cells.size()},
                                                     nb::handle());
          },
          nb::rv_policy::reference_internal, "uncut_cells (the name of HOMeshPart).")
      .def_prop_ro(
          "cut_cell_ids",
          [](const PartT& self)
          {
            // positions in the result's cut cells, as HOMeshPart numbered them
            const std::vector<int>& all = self.part.result->cut_cells;
            std::vector<int> ids;
            ids.reserve(self.part.cut_cells.size());
            for (const int c : self.part.cut_cells)
              ids.push_back(static_cast<int>(std::lower_bound(all.begin(), all.end(), c) - all.begin()));
            return as_nbarray(std::move(ids));
          },
          nb::rv_policy::move, "Positions of the part's cut cells in result.parent_cell_ids.")
      .def_prop_ro(
          "zero_faces",
          [](const PartT& self)
          {
            return nb::ndarray<const int, nb::numpy>(self.part.zero_faces.data(), {self.part.zero_faces.size()},
                                                     nb::handle());
          },
          nb::rv_policy::reference_internal, "Indices into the result's zero_faces.")
      .def(
          "quadrature",
          [](const PartT& self, int order, const std::string& mode, nb::handle backend, nb::handle options)
          {
            const std::string b = call_backend(backend, self.backend);
            const bool cut_only = part_mode_is_cut_only(mode);
            if (b == "lut")
            {
              const auto o = call_options<cutcells::lut::Options>(options, self.options, b, "LutOptions");
              nb::gil_scoped_release release;
              return cutcells::part::quadrature_rules(self.part, order, !cut_only, o);
            }
            const auto o = call_options<qr::Options>(options, self.options, b, "QuadraysOptions");
            nb::gil_scoped_release release;
            return cutcells::part::quadrature_rules(self.part, order, !cut_only, o);
          },
          nb::arg("order") = 3, nb::arg("mode") = "full", nb::arg("backend") = nb::none(),
          nb::arg("options") = nb::none(),
          "Quadrature rules, one per cell: the backend's on cut cells, rules on owned "
          "zero faces, and with mode 'full' those of the whole cells. backend: 'quadrays' "
          "(options: QuadraysOptions) or 'lut' ('straight'), the lookup tables on Pk-iso-P1 "
          "templates (options: LutOptions); None: the part's. quadrays takes order "
          "Gauss-Legendre points per segment; the lookup tables' straight pieces, whole "
          "cells and faces get rules exact for degree 2 order - 1 (at most 10).")
      .def(
          "visualization_mesh",
          [](const PartT& self, const std::string& mode, nb::handle backend, int degree,
             nb::handle options) -> nb::object
          {
            const std::string b = call_backend(backend, self.backend);
            const bool cut_only = part_mode_is_cut_only(mode);
            if (b == "lut")
            {
              const auto o = call_options<cutcells::lut::Options>(options, self.options, b, "LutOptions");
              cutcells::mesh::CutMesh<T> out;
              {
                nb::gil_scoped_release release;
                out = cutcells::part::visualization_mesh(self.part, !cut_only, o);
              }
              return nb::cast(std::move(out));
            }
            const auto o = call_options<qr::Options>(options, self.options, b, "QuadraysOptions");
            qr::LeafMesh<T> out;
            {
              nb::gil_scoped_release release;
              out = cutcells::part::visualization_mesh(self.part, degree, !cut_only, o);
            }
            return nb::cast(std::move(out));
          },
          nb::arg("mode") = "full", nb::arg("backend") = nb::none(), nb::arg("degree") = 3,
          nb::arg("options") = nb::none(),
          "Cells for visualisation: the backend's pieces of cut cells, zero faces, and with "
          "mode 'full' the whole cells. quadrays gives a QuadraysLeafMesh of Lagrange cells "
          "of the given degree; the lookup tables a CutMesh of straight cells.")
      .def(
          "write_vtu",
          [](const PartT& self, const std::string& filename, const std::string& mode, nb::handle backend,
             int degree, nb::handle options)
          {
            const std::string b = call_backend(backend, self.backend);
            const bool cut_only = part_mode_is_cut_only(mode);
            if (b == "lut")
            {
              const auto o = call_options<cutcells::lut::Options>(options, self.options, b, "LutOptions");
              nb::gil_scoped_release release;
              cutcells::part::write_vtu(filename, self.part, !cut_only, o);
              return;
            }
            const auto o = call_options<qr::Options>(options, self.options, b, "QuadraysOptions");
            nb::gil_scoped_release release;
            cutcells::part::write_vtu(filename, self.part, degree, !cut_only, o);
          },
          nb::arg("filename"), nb::arg("mode") = "full", nb::arg("backend") = nb::none(), nb::arg("degree") = 3,
          nb::arg("options") = nb::none(), "Write visualization_mesh to a .vtu file.");

  part_module.def(
      ("cut_" + type).c_str(),
      [](const MeshViewT& mesh, nb::handle level_sets, nb::handle names, int max_depth, const std::string& backend,
         nb::handle options)
      {
        ResultT r = cut_mesh_part<T>(mesh, level_sets, names, max_depth);
        r.backend = part_backend(backend);
        r.options = checked_options(options, r.backend);
        return r;
      },
      nb::arg("mesh"), nb::arg("level_sets"), nb::arg("names") = nb::none(), nb::arg("max_depth") = 12,
      nb::arg("backend") = "quadrays", nb::arg("options") = nb::none(),
      "Classify every cell of the mesh by every level set (LevelSetFunctions with dof "
      "values, or AnalyticLevelSets, alone or in a list) by their own bounds. Analytic "
      "level sets are named 'phi', or 'phi1', 'phi2', ... in a list, unless names are "
      "given. max_depth: bisections of a cell before an unproven sign counts as cut. "
      "backend and options: the default of the result's parts.");

  // cutcells.cut and cutcells.ho_cut: the same, with the lookup tables by
  // default and the keywords of the former cut()
  for (const char* name : {"cut", "ho_cut"})
  {
    m.def(
        name,
        [](const MeshViewT& mesh, nb::handle level_sets, nb::handle names, int max_depth, const std::string& backend,
           nb::handle options, bool triangulate, const std::string& triangulation,
           const std::string& cut_approximation, int cut_approximation_order, nb::handle degree, nb::handle name)
        {
          nb::object given = nb::borrow(names);
          if (!name.is_none())
          {
            if (!names.is_none())
              throw std::invalid_argument("cut: give name or names, not both");
            given = nb::make_tuple(name);
          }
          ResultT r = cut_mesh_part<T>(mesh, level_sets, given, max_depth);
          r.backend = part_backend(backend);
          if (!options.is_none())
            r.options = checked_options(options, r.backend);
          else if (r.backend == "lut")
            r.options = nb::cast(legacy_lut_options(triangulate, triangulation, cut_approximation,
                                                    cut_approximation_order, degree));
          return r;
        },
        nb::arg("mesh"), nb::arg("level_sets"), nb::arg("names") = nb::none(), nb::arg("max_depth") = 12,
        nb::arg("backend") = "lut", nb::arg("options") = nb::none(), nb::arg("triangulate") = false,
        nb::arg("triangulation") = "classical", nb::arg("cut_approximation") = "auto",
        nb::arg("cut_approximation_order") = 1, nb::arg("degree") = nb::none(), nb::arg("name") = nb::none(),
        "Cut a MeshView by level sets (LevelSetFunctions or AnalyticLevelSets, alone or in "
        "a list): cutcells.part.cut with the lookup tables ('lut', 'straight') as the "
        "default backend of the result's parts. Without options, the lookup tables take "
        "triangulate and triangulation ('classical', 'midpoint') and the template order "
        "from cut_approximation: 'auto' (the level sets' degree, or degree for analytic "
        "ones), 'linear' (1) or 'iso_p1' (cut_approximation_order). name: of a single "
        "level set. Returns an HOCutResult (cutcells.part.CutResult); "
        "result[\"phi1 < 0 and phi2 = 0\"] selects a part.");
  }

  m.def(("quadrays_quadrature_" + type).c_str(),
        [](const PartT& part, int order, const std::string& mode, nb::handle options)
        {
          const bool cut_only = part_mode_is_cut_only(mode);
          const auto o = call_options<qr::Options>(options, nb::handle(), "quadrays", "QuadraysOptions");
          nb::gil_scoped_release release;
          return cutcells::part::quadrature_rules(part.part, order, !cut_only, o);
        },
        nb::arg("part"), nb::arg("order") = 3, nb::arg("mode") = "full", nb::arg("options") = nb::none(),
        "part.quadrature(order, mode, backend='quadrays', options). "
        "order: Gauss-Legendre points per segment of each height line.");

  m.def(("quadrays_leaves_" + type).c_str(),
        [](const PartT& part, int degree, const std::string& mode, nb::handle options)
        {
          const bool cut_only = part_mode_is_cut_only(mode);
          const auto o = call_options<qr::Options>(options, nb::handle(), "quadrays", "QuadraysOptions");
          nb::gil_scoped_release release;
          return cutcells::part::visualization_mesh(part.part, degree, !cut_only, o);
        },
        nb::arg("part"), nb::arg("degree") = 3, nb::arg("mode") = "full", nb::arg("options") = nb::none(),
        "part.visualization_mesh(mode, backend='quadrays', degree, options): the pieces "
        "quadrays integrates in cut cells as Lagrange cells of the given degree; with mode "
        "'full', the whole cells of volume parts as linear cells.");

  m.attr(("HOCutResult_" + type).c_str()) = part_module.attr(result_name.c_str());
  m.attr(("HOMeshPart_" + type).c_str()) = part_module.attr(part_name.c_str());
  if constexpr (std::is_same_v<T, double>)
  {
    part_module.attr("CutResult") = part_module.attr(result_name.c_str());
    part_module.attr("MeshPart") = part_module.attr(part_name.c_str());
    part_module.attr("cut") = part_module.attr("cut_float64");
    m.attr("HOCutResult") = part_module.attr(result_name.c_str());
    m.attr("HOMeshPart") = part_module.attr(part_name.c_str());
    m.attr("quadrays_quadrature") = m.attr("quadrays_quadrature_float64");
    m.attr("quadrays_leaves") = m.attr("quadrays_leaves_float64");
  }
}

// ============================================================================
// Compression of quadrature rules (compression/compress.h)
// ============================================================================

template <typename T>
void declare_compression(nb::module_& m, const std::string& type)
{
  namespace cc = cutcells::compression;
  m.def(("compress_rules_" + type).c_str(),
        [](const quadrature::QuadratureRules<T>& rules, int degree, const std::string& space)
        {
          const cc::MomentSpace s = cc::string_to_moment_space(space);
          quadrature::QuadratureRules<T> out;
          cc::CompressionStats stats;
          {
            nb::gil_scoped_release release;
            cc::compress_rules(rules, degree, s, out, stats);
          }
          return std::make_pair(std::move(out), stats);
        },
        nb::arg("rules"), nb::arg("degree"), nb::arg("space") = "tensor",
        "Every rule replaced by at most as many of its points as the polynomial space has "
        "moments, with positive weights that integrate the space exactly as the original "
        "rule does: space 'tensor' (degree <= degree in each reference coordinate, (degree + 1)^tdim "
        "moments; hexahedra) or 'total' (total degree <= degree; tetrahedra). Q_k elements "
        "on affine hexahedra need degree 2k for stiffness and mass, P_k on affine tetrahedra "
        "space 'total' and degree 2k. Rules with a negative weight, or no more points than "
        "moments, are copied. Returns (QuadratureRules, CompressionStats).");

  if constexpr (std::is_same_v<T, double>)
    m.attr("compress_rules") = m.attr("compress_rules_float64");
}

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

  nb::class_<cutcells::IsoRefineTemplate>(
      m, "IsoRefineTemplate", "Topology-only Pk-iso-P1 refinement template.")
      .def_prop_ro("n_vertices",
                   [](const cutcells::IsoRefineTemplate& self) { return self.n_vertices; })
      .def_prop_ro("n_cells",
                   [](const cutcells::IsoRefineTemplate& self) { return self.n_cells; })
      .def_prop_ro("tdim",
                   [](const cutcells::IsoRefineTemplate& self) { return self.tdim; })
      .def_prop_ro("vertices_per_cell",
                   [](const cutcells::IsoRefineTemplate& self) { return self.vertices_per_cell; })
      .def_prop_ro("parent_cell_type",
                   [](const cutcells::IsoRefineTemplate& self) { return self.parent_cell_type; })
      .def_prop_ro("child_cell_type",
                   [](const cutcells::IsoRefineTemplate& self) { return self.child_cell_type; })
      .def_prop_ro(
          "ref_vertex_coords",
          [](const cutcells::IsoRefineTemplate& self)
          {
            return nb::ndarray<const double, nb::numpy>(
                self.ref_vertex_coords.data(),
                {self.ref_vertex_coords.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "vertex_parent_dim",
          [](const cutcells::IsoRefineTemplate& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.vertex_parent_dim.data(),
                {self.vertex_parent_dim.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "vertex_parent_id",
          [](const cutcells::IsoRefineTemplate& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.vertex_parent_id.data(),
                {self.vertex_parent_id.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_connectivity",
          [](const cutcells::IsoRefineTemplate& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.cell_connectivity.data(),
                {self.cell_connectivity.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "cell_offsets",
          [](const cutcells::IsoRefineTemplate& self)
          {
            return nb::ndarray<const int, nb::numpy>(
                self.cell_offsets.data(),
                {self.cell_offsets.size()},
                nb::cast(self, nb::rv_policy::reference));
          },
          nb::rv_policy::reference_internal,
          "Children in CSR layout: child c has cell_connectivity[cell_offsets[c]:cell_offsets[c + 1]].")
      .def_prop_ro(
          "cell_types",
          [](const cutcells::IsoRefineTemplate& self) { return self.cell_types; },
          "The type of each child (a pyramid's template mixes pyramids and tetrahedra).");

  m.attr("RefinementTemplate") = m.attr("IsoRefineTemplate");

  m.def("iso_p1_template",
        [](cell::type cell_type, int order) -> const cutcells::IsoRefineTemplate&
        {
          return cutcells::iso_p1_template(cell_type, order);
        },
        nb::arg("cell_type"),
        nb::arg("order"),
        nb::rv_policy::reference,
        "Return a topology-only Pk-iso-P1 refinement template.");

  m.def("iso_p1_ref_coords",
        [](cell::type cell_type, int order)
        {
          auto x = cutcells::iso_p1_ref_coords(cell_type, order);
          return nb::ndarray<const double, nb::numpy>(
              x.data(), {x.size()}, nb::handle());
        },
        nb::arg("cell_type"),
        nb::arg("order"),
        "Return flat reference coordinates for a Pk-iso-P1 template.");

  declare_float<float>(m, "float32");
  declare_float<double>(m, "float64");

  declare_meshview_and_levelset<float>(m, "float32");
  declare_meshview_and_levelset<double>(m, "float64");
  declare_level_set_cell<float>(m, "float32");
  declare_level_set_cell<double>(m, "float64");

  declare_analytic(m);

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

  declare_write_vtk<float>(m);
  declare_write_vtk<double>(m);

  // ---- quadrays: height-function quadrature ----
  nb::class_<cutcells::quadrays::Options>(m, "QuadraysOptions",
      "Options of the quadrays engine.")
      .def(nb::init<>())
      .def_rw("margin", &cutcells::quadrays::Options::margin,
              "Required |d_k psi| / |grad psi| on a box for a height direction; "
              "also the accuracy control.")
      .def_rw("max_depth", &cutcells::quadrays::Options::max_depth,
              "Bisections allowed per level along a branch.")
      .def_rw("max_bisections", &cutcells::quadrays::Options::max_bisections,
              "Bisections allowed per cell; beyond it boxes are integrated uncertified.")
      .def_rw("prune_bounds", &cutcells::quadrays::Options::prune_bounds,
              "Drop bounds of the height lines that are never active.")
      .def_rw("split_bounds", &cutcells::quadrays::Options::split_bounds,
              "Split the base where the active bounds of the height lines change.")
      .def_rw("diagonal_frames", &cutcells::quadrays::Options::diagonal_frames,
              "Level 2: try the diagonal frame before bisecting.")
      .def_rw("rotation_depth", &cutcells::quadrays::Options::rotation_depth,
              "Where no axis suits the zero sets of two level sets, try a frame rotated between "
              "their normals: at level 2 always, in 3D from this depth of bisection on.")
      .def_rw("mask_subdivisions", &cutcells::quadrays::Options::mask_subdivisions,
              "M > 1: margins from M^D sub-cells; 1: bounds on the whole box.")
      .def_rw("taylor_subdivisions", &cutcells::quadrays::Options::taylor_subdivisions,
              "Analytic level sets: M > 1: Taylor models over M^D sub-boxes on boxes a clip plane "
              "cuts; 1: over the whole box.")
      .def_rw("two_roots_depth", &cutcells::quadrays::Options::two_roots_depth,
              "From this depth of bisection on, accept a direction with two roots per height line "
              "that never merge in the box (two sheets of one level set).")
      .def_rw("diagnose", &cutcells::quadrays::Options::diagnose,
              "Record why each bisection happened in QuadraysStats.causes.");
  nb::class_<cutcells::lut::Options>(m, "LutOptions", "Options of the lookup-table backend.")
      .def(
          "__init__",
          [](cutcells::lut::Options* self, int template_order, bool triangulate, const std::string& triangulation)
          {
            new (self) cutcells::lut::Options{template_order, triangulate,
                                              cell::triangulation_strategy_from_string(triangulation)};
          },
          nb::arg("template_order") = 0, nb::arg("triangulate") = false, nb::arg("triangulation") = "classical")
      .def_rw("template_order", &cutcells::lut::Options::template_order,
              "Order k of the Pk-iso-P1 template that subdivides a cut cell, 1 to 4; 0: the "
              "highest degree of the level sets that cut the cell, 2 for analytic ones.")
      .def_rw("triangulate", &cutcells::lut::Options::triangulate, "Split the cut pieces into simplices.")
      .def_prop_rw(
          "triangulation",
          [](const cutcells::lut::Options& self)
          { return std::string(cell::triangulation_strategy_to_string(self.triangulation)); },
          [](cutcells::lut::Options& self, const std::string& triangulation)
          { self.triangulation = cell::triangulation_strategy_from_string(triangulation); },
          "How triangulate splits the pieces: 'classical' or 'midpoint'.");
  nb::class_<cutcells::quadrays::Stats>(m, "QuadraysStats",
      "Counters of quadrays engine runs.")
      .def(nb::init<>())
      .def_ro("bisections", &cutcells::quadrays::Stats::bisections)
      .def_ro("uncertified", &cutcells::quadrays::Stats::uncertified)
      .def_ro("rotations", &cutcells::quadrays::Stats::rotations)
      .def_ro("incomplete_leaves", &cutcells::quadrays::Stats::incomplete_leaves)
      .def_ro("two_roots", &cutcells::quadrays::Stats::two_roots)
      .def_ro("surfaces", &cutcells::quadrays::Stats::surfaces)
      .def_ro("causes", &cutcells::quadrays::Stats::causes);
  declare_quadrays<float>(m, "float32");
  declare_quadrays<double>(m, "float64");

  nb::class_<cutcells::compression::CompressionStats>(m, "CompressionStats",
      "Counters of compress_rules.")
      .def(nb::init<>())
      .def_ro("n_rules", &cutcells::compression::CompressionStats::n_rules, "Rules read.")
      .def_ro("n_compressed", &cutcells::compression::CompressionStats::n_compressed, "Rules reduced.")
      .def_ro("n_skipped", &cutcells::compression::CompressionStats::n_skipped,
              "Rules with a negative weight, copied unchanged.")
      .def_ro("points_before", &cutcells::compression::CompressionStats::points_before, "Points read.")
      .def_ro("points_after", &cutcells::compression::CompressionStats::points_after, "Points written.")
      .def_ro("max_residual", &cutcells::compression::CompressionStats::max_residual,
              "Largest moment error relative to a rule's total |weight|.");
  declare_compression<float>(m, "float32");
  declare_compression<double>(m, "float64");

  nb::module_ part_module = m.def_submodule(
      "part", "The front end: cut(mesh, level_sets) classifies cells by the level sets' "
              "own bounds; result[expr] selects a MeshPart, whose quadrature and "
              "visualisation come from a backend (quadrays, or the lookup tables: 'lut').");
  declare_part<float>(m, part_module, "float32");
  declare_part<double>(m, part_module, "float64");
}
