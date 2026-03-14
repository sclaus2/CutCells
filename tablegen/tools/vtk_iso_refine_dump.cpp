// Dump VTK tessellations for higher-order reference cells.
//
// Output format is JSON to stdout:
// {
//   "cell_type": "...",
//   "order": p,
//   "tdim": d,
//   "npoints": N,
//   "simplex_size": d+1,
//   "points": [[x,y,z], ...],
//   "simplices_vtk": [[i0,i1,...], ...]
// }
//
// Build:
//   cmake -S tablegen/tools -B tablegen/tools/build
//   cmake --build tablegen/tools/build --target vtk_iso_refine_dump

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkHigherOrderHexahedron.h>
#include <vtkHigherOrderQuadrilateral.h>
#include <vtkHigherOrderTetra.h>
#include <vtkHigherOrderTriangle.h>
#include <vtkHigherOrderWedge.h>
#include <vtkIdList.h>
#include <vtkLagrangeHexahedron.h>
#include <vtkLagrangeQuadrilateral.h>
#include <vtkLagrangeTetra.h>
#include <vtkLagrangeTriangle.h>
#include <vtkLagrangeWedge.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTriQuadraticPyramid.h>

#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

int tdim_for_cell(const std::string& ct)
{
  if (ct == "interval")
    return 1;
  if (ct == "triangle" || ct == "quadrilateral")
    return 2;
  if (ct == "tetrahedron" || ct == "hexahedron" || ct == "prism" || ct == "pyramid")
    return 3;
  throw std::invalid_argument("Unknown cell type: " + ct);
}

int lagrange_point_count(const std::string& ct, int p)
{
  if (p < 1)
    throw std::invalid_argument("Order must be >= 1");
  if (ct == "interval")
    return p + 1;
  if (ct == "triangle")
    return (p + 1) * (p + 2) / 2;
  if (ct == "quadrilateral")
    return (p + 1) * (p + 1);
  if (ct == "tetrahedron")
    return (p + 1) * (p + 2) * (p + 3) / 6;
  if (ct == "hexahedron")
    return (p + 1) * (p + 1) * (p + 1);
  if (ct == "prism")
    return (p + 1) * (p + 1) * (p + 2) / 2;
  if (ct == "pyramid")
    return (p + 1) * (p + 2) * (2 * p + 3) / 6;
  throw std::invalid_argument("Unknown cell type: " + ct);
}

vtkSmartPointer<vtkCell> make_cell(const std::string& ct, int p)
{
  if (ct == "triangle")
  {
    return vtkSmartPointer<vtkLagrangeTriangle>::New();
  }
  if (ct == "quadrilateral")
  {
    auto c = vtkSmartPointer<vtkLagrangeQuadrilateral>::New();
    c->SetUniformOrderFromNumPoints(lagrange_point_count(ct, p));
    return c;
  }
  if (ct == "tetrahedron")
  {
    return vtkSmartPointer<vtkLagrangeTetra>::New();
  }
  if (ct == "hexahedron")
  {
    auto c = vtkSmartPointer<vtkLagrangeHexahedron>::New();
    c->SetUniformOrderFromNumPoints(lagrange_point_count(ct, p));
    return c;
  }
  if (ct == "prism")
  {
    auto c = vtkSmartPointer<vtkLagrangeWedge>::New();
    c->SetUniformOrderFromNumPoints(lagrange_point_count(ct, p));
    return c;
  }
  if (ct == "pyramid")
  {
    if (p != 2)
      throw std::invalid_argument("VTK in this build supports pyramid tessellation only for order 2");
    auto c = vtkSmartPointer<vtkTriQuadraticPyramid>::New();
    return c;
  }
  throw std::invalid_argument("Unsupported cell type for VTK dump: " + ct);
}

int configure_cell(vtkCell* cell, const std::string& ct, int p)
{
  const int npts_expected = lagrange_point_count(ct, p);

  // In VTK 9.6, triangle/tetra infer order from point count.
  if (ct == "triangle")
  {
    auto* tri = vtkHigherOrderTriangle::SafeDownCast(cell);
    if (tri == nullptr)
      throw std::runtime_error("Internal error: triangle cell cast failed");
    tri->GetPointIds()->SetNumberOfIds(npts_expected);
    for (int i = 0; i < npts_expected; ++i)
      tri->GetPointIds()->SetId(i, i);
    tri->ComputeOrder();
    return npts_expected;
  }
  if (ct == "tetrahedron")
  {
    auto* tet = vtkHigherOrderTetra::SafeDownCast(cell);
    if (tet == nullptr)
      throw std::runtime_error("Internal error: tetrahedron cell cast failed");
    tet->GetPointIds()->SetNumberOfIds(npts_expected);
    for (int i = 0; i < npts_expected; ++i)
      tet->GetPointIds()->SetId(i, i);
    tet->ComputeOrder();
    return npts_expected;
  }

  // Quadrilateral/hexahedron/wedge are configured from expected Lagrange point count.
  if (ct == "quadrilateral" || ct == "hexahedron" || ct == "prism")
    return npts_expected;

  // Pyramid currently uses vtkTriQuadraticPyramid in this VTK build.
  const int npts_pyr = static_cast<int>(cell->GetNumberOfPoints());
  if (npts_pyr <= 0)
    throw std::runtime_error("VTK pyramid has non-positive point count");
  return npts_pyr;
}

std::vector<std::array<double, 3>> get_parametric_points(vtkCell* cell, int npts)
{
  const double* pc = cell->GetParametricCoords();
  if (pc == nullptr)
    throw std::runtime_error("VTK cell returned null parametric coordinates");

  std::vector<std::array<double, 3>> out(npts);
  for (int i = 0; i < npts; ++i)
  {
    out[i] = {pc[3 * i + 0], pc[3 * i + 1], pc[3 * i + 2]};
  }
  return out;
}

std::vector<std::vector<int>> triangulate_ids(vtkCell* cell, int tdim)
{
  vtkNew<vtkIdList> out_ids;
  vtkNew<vtkPoints> out_pts;
  const int ok = cell->Triangulate(0, out_ids, out_pts);
  if (!ok)
    throw std::runtime_error("VTK Triangulate failed");

  const int simplex_size = tdim + 1;
  if (out_ids->GetNumberOfIds() % simplex_size != 0)
    throw std::runtime_error("Triangulate returned unexpected id count");

  std::vector<std::vector<int>> simplices;
  simplices.reserve(out_ids->GetNumberOfIds() / simplex_size);
  for (vtkIdType i = 0; i < out_ids->GetNumberOfIds(); i += simplex_size)
  {
    std::vector<int> s(simplex_size);
    for (int j = 0; j < simplex_size; ++j)
      s[j] = static_cast<int>(out_ids->GetId(i + j));
    simplices.push_back(std::move(s));
  }
  return simplices;
}

void print_json(const std::string& ct,
                int order,
                int tdim,
                const std::vector<std::array<double, 3>>& points,
                const std::vector<std::vector<int>>& simplices)
{
  std::cout << "{\n";
  std::cout << "  \"cell_type\": \"" << ct << "\",\n";
  std::cout << "  \"order\": " << order << ",\n";
  std::cout << "  \"tdim\": " << tdim << ",\n";
  std::cout << "  \"npoints\": " << points.size() << ",\n";
  std::cout << "  \"simplex_size\": " << (tdim + 1) << ",\n";
  std::cout << "  \"points\": [\n";
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    const auto& p = points[i];
    std::cout << "    [" << p[0] << ", " << p[1] << ", " << p[2] << "]";
    if (i + 1 < points.size())
      std::cout << ",";
    std::cout << "\n";
  }
  std::cout << "  ],\n";
  std::cout << "  \"simplices_vtk\": [\n";
  for (std::size_t i = 0; i < simplices.size(); ++i)
  {
    std::cout << "    [";
    for (std::size_t j = 0; j < simplices[i].size(); ++j)
    {
      std::cout << simplices[i][j];
      if (j + 1 < simplices[i].size())
        std::cout << ", ";
    }
    std::cout << "]";
    if (i + 1 < simplices.size())
      std::cout << ",";
    std::cout << "\n";
  }
  std::cout << "  ]\n";
  std::cout << "}\n";
}

void usage(const char* exe)
{
  std::cerr << "Usage: " << exe << " --cell <cell_type> --order <p>\n";
  std::cerr << "  cell_type: triangle|quadrilateral|tetrahedron|hexahedron|prism|pyramid\n";
}

} // namespace

int main(int argc, char** argv)
{
  try
  {
    std::cout << std::setprecision(17);

    std::string cell_type;
    int order = -1;
    for (int i = 1; i < argc; ++i)
    {
      const std::string a = argv[i];
      if (a == "--cell" && i + 1 < argc)
      {
        cell_type = argv[++i];
      }
      else if (a == "--order" && i + 1 < argc)
      {
        order = std::stoi(argv[++i]);
      }
      else
      {
        usage(argv[0]);
        return 2;
      }
    }

    if (cell_type.empty() || order < 1)
    {
      usage(argv[0]);
      return 2;
    }

    const int tdim = tdim_for_cell(cell_type);
    auto cell = make_cell(cell_type, order);
    const int npts = configure_cell(cell, cell_type, order);

    vtkNew<vtkPoints> pts;
    pts->SetNumberOfPoints(npts);
    std::vector<vtkIdType> ids(npts);
    for (int i = 0; i < npts; ++i)
    {
      ids[i] = i;
      pts->SetPoint(i, 0.0, 0.0, 0.0);
    }
    cell->Initialize(npts, ids.data(), pts);
    cell->Initialize();

    // Get parametric coordinates now that the higher-order cell has been initialized.
    const auto param = get_parametric_points(cell, npts);
    for (int i = 0; i < npts; ++i)
      pts->SetPoint(i, param[i][0], param[i][1], param[i][2]);
    cell->Initialize(npts, ids.data(), pts);
    cell->Initialize();

    const auto simplices = triangulate_ids(cell, tdim);
    print_json(cell_type, order, tdim, param, simplices);
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
