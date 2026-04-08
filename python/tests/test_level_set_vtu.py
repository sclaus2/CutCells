from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import cutcells


TRIANGLE_BASIX_TO_VTK = {
    2: [0, 1, 2, 5, 3, 4],
    3: [0, 1, 2, 7, 8, 3, 4, 6, 5, 9],
    4: [0, 1, 2, 9, 10, 11, 3, 4, 5, 8, 7, 6, 12, 13, 14],
}

TETRA_BASIX_TO_VTK = {
    2: [0, 1, 2, 3, 9, 6, 8, 7, 5, 4],
    3: [0, 1, 2, 3, 14, 15, 8, 9, 13, 12, 10, 11, 6, 7, 4, 5, 18, 16, 17, 19],
    4: [
        0,
        1,
        2,
        3,
        19,
        20,
        21,
        10,
        11,
        12,
        18,
        17,
        16,
        13,
        14,
        15,
        7,
        8,
        9,
        4,
        5,
        6,
        28,
        29,
        30,
        23,
        24,
        22,
        25,
        27,
        26,
        31,
        33,
        32,
        34,
    ],
}


def _read_dataarray(path: Path, name: str) -> str:
    root = ET.parse(path).getroot()
    for data_array in root.iter("DataArray"):
        if data_array.attrib.get("Name") == name:
            return (data_array.text or "").strip()
    raise AssertionError(f"DataArray '{name}' not found in {path}")


def _find_dataarray(path: Path, name: str) -> ET.Element:
    root = ET.parse(path).getroot()
    for data_array in root.iter("DataArray"):
        if data_array.attrib.get("Name") == name:
            return data_array
    raise AssertionError(f"DataArray '{name}' not found in {path}")


def _make_single_cell_level_set(cell_family: str, degree: int):
    if cell_family == "triangle":
        local_dofs = (degree + 1) * (degree + 2) // 2
        vtk_type = 5
        tdim = 2
    elif cell_family == "tetra":
        local_dofs = (degree + 1) * (degree + 2) * (degree + 3) // 6
        vtk_type = 10
        tdim = 3
    else:
        raise ValueError(f"Unsupported cell_family={cell_family}")

    dof_coordinates = np.zeros((local_dofs, 3), dtype=np.float64)
    dof_coordinates[:, 0] = np.arange(local_dofs, dtype=np.float64)
    cell_dofs = np.arange(local_dofs, dtype=np.int32)
    cell_offsets = np.array([0, local_dofs], dtype=np.int32)
    cell_types = np.array([vtk_type], dtype=np.int32)
    dof_values = np.arange(local_dofs, dtype=np.float64) + 0.125

    mesh_data = cutcells.create_level_set_mesh_data(
        dof_coordinates=dof_coordinates,
        cell_dofs=cell_dofs,
        cell_offsets=cell_offsets,
        degree=degree,
        tdim=tdim,
        cell_types=cell_types,
    )
    ls = cutcells.create_level_set_function(mesh_data, dof_values)
    return ls, dof_values


@pytest.mark.parametrize("degree", [2, 3, 4])
def test_write_level_set_vtu_triangle_connectivity_permutation(tmp_path, degree):
    ls, dof_values = _make_single_cell_level_set("triangle", degree)
    output = tmp_path / f"triangle_p{degree}.vtu"
    cutcells.write_level_set_vtu(str(output), ls, field_name="phi")

    connectivity = [int(x) for x in _read_dataarray(output, "connectivity").split()]
    types = [int(x) for x in _read_dataarray(output, "types").split()]
    phi = np.fromstring(_read_dataarray(output, "phi"), sep=" ")

    assert connectivity == TRIANGLE_BASIX_TO_VTK[degree]
    assert types == [69]  # VTK_LAGRANGE_TRIANGLE
    np.testing.assert_allclose(phi, dof_values)


@pytest.mark.parametrize("degree", [2, 3, 4])
def test_write_level_set_vtu_tetra_connectivity_permutation(tmp_path, degree):
    ls, dof_values = _make_single_cell_level_set("tetra", degree)
    output = tmp_path / f"tetra_p{degree}.vtu"
    cutcells.write_level_set_vtu(str(output), ls, field_name="phi")

    connectivity = [int(x) for x in _read_dataarray(output, "connectivity").split()]
    types = [int(x) for x in _read_dataarray(output, "types").split()]
    phi = np.fromstring(_read_dataarray(output, "phi"), sep=" ")

    assert connectivity == TETRA_BASIX_TO_VTK[degree]
    assert types == [71]  # VTK_LAGRANGE_TETRAHEDRON
    np.testing.assert_allclose(phi, dof_values)


def test_write_level_set_vtu_xml_matches_vtk_expectations(tmp_path):
    ls, _ = _make_single_cell_level_set("triangle", 3)
    output = tmp_path / "triangle_p3.vtu"
    cutcells.write_level_set_vtu(str(output), ls, field_name="phi")

    root = ET.parse(output).getroot()
    assert root.attrib["type"] == "UnstructuredGrid"
    assert root.attrib["version"] == "0.1"
    assert root.attrib["byte_order"] in {"LittleEndian", "BigEndian"}

    types = _find_dataarray(output, "types")
    assert types.attrib["type"] == "UInt8"


@pytest.mark.parametrize(
    ("cell_family", "degree", "expected_type", "expected_local_dofs"),
    [
        ("triangle", 2, 69, 6),
        ("triangle", 3, 69, 10),
        ("triangle", 4, 69, 15),
        ("tetra", 2, 71, 10),
        ("tetra", 3, 71, 20),
        ("tetra", 4, 71, 35),
    ],
)
def test_write_level_set_vtu_cell_widths_match_degree(
    tmp_path, cell_family, degree, expected_type, expected_local_dofs
):
    ls, _ = _make_single_cell_level_set(cell_family, degree)
    output = tmp_path / f"{cell_family}_p{degree}.vtu"
    cutcells.write_level_set_vtu(str(output), ls, field_name="phi")

    connectivity = [int(x) for x in _read_dataarray(output, "connectivity").split()]
    offsets = [int(x) for x in _read_dataarray(output, "offsets").split()]
    types = [int(x) for x in _read_dataarray(output, "types").split()]

    prev = 0
    widths = []
    for off in offsets:
        widths.append(off - prev)
        prev = off

    assert widths == [expected_local_dofs]
    assert len(connectivity) == expected_local_dofs
    assert types == [expected_type]


@pytest.mark.parametrize("cell_family", ["triangle", "tetra"])
@pytest.mark.parametrize("degree", [2, 3, 4])
def test_write_level_set_vtu_roundtrip_with_pyvista(tmp_path, cell_family, degree):
    pv = pytest.importorskip("pyvista")

    if cell_family == "triangle":
        grid = cutcells.rectangle_triangle_mesh(-1.0, -1.0, 1.0, 1.0, 5, 5)
        mesh = cutcells.mesh_from_pyvista(grid, tdim=2)
        expected_type = 69
    else:
        grid = cutcells.box_tetrahedron_mesh(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 4, 4, 4)
        mesh = cutcells.mesh_from_pyvista(grid, tdim=3)
        expected_type = 71

    def phi(X):
        return np.sqrt((X[:, 0] - 0.15) ** 2 + (X[:, 1] + 0.1) ** 2 + (X[:, 2] - 0.05) ** 2) - 0.55

    ls = cutcells.interpolate_level_set(mesh, phi, degree=degree)
    output = tmp_path / f"{cell_family}_p{degree}.vtu"
    cutcells.write_level_set_vtu(str(output), ls, field_name="phi")

    out_grid = pv.read(output)
    assert out_grid.n_points == ls.mesh_data.num_dofs()
    assert "phi" in out_grid.point_data
    assert len(np.asarray(out_grid.point_data["phi"])) == ls.mesh_data.num_dofs()
    assert np.all(np.asarray(out_grid.celltypes) == expected_type)
    assert all(out_grid.GetCell(i).GetOrder() == degree for i in range(out_grid.n_cells))
    assert all(not out_grid.GetCell(i).IsLinear() for i in range(out_grid.n_cells))
