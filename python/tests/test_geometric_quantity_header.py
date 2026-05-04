import os
import shlex
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


def test_geometric_quantity_header_primitives(tmp_path):
    cxx = os.environ.get("CXX")
    if cxx:
        compiler = shlex.split(cxx)
    else:
        path = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if path is None:
            pytest.skip("C++ compiler not available")
        compiler = [path]

    repo = Path(__file__).resolve().parents[2]
    source = tmp_path / "geometric_quantity_check.cpp"
    exe = tmp_path / "geometric_quantity_check"
    source.write_text(
        textwrap.dedent(
            r"""
            #include "geometric_quantity.h"

            #include <array>
            #include <cassert>
            #include <cmath>
            #include <span>

            namespace geom = cutcells::geom;
            namespace cell = cutcells::cell;

            bool close(double a, double b)
            {
                return std::fabs(a - b) < 1.0e-12;
            }

            int main()
            {
                const auto tet = cell::type::tetrahedron;

                // 1. Projecting a 3D gradient into a face plane.
                // Tetrahedron face 3 is the z = 0 face with normal (0, 0, 1).
                std::array<double, 3> gradient = {1.0, 2.0, 3.0};
                const auto face_gradient =
                    geom::project_into_parent_face_tangent<double>(
                        tet, 3, std::span<const double>(gradient));
                assert(!face_gradient.degenerate());
                assert(close(face_gradient.value[0], 1.0));
                assert(close(face_gradient.value[1], 2.0));
                assert(close(face_gradient.value[2], 0.0));

                // 2. Computing an in-face normal to a segment.
                std::array<double, 3> a = {0.0, 0.0, 0.0};
                std::array<double, 3> b = {1.0, 0.0, 0.0};
                std::array<double, 3> face_normal = {0.0, 0.0, 1.0};
                const auto in_face_normal =
                    geom::in_face_segment_normal<double>(
                        std::span<const double>(a),
                        std::span<const double>(b),
                        std::span<const double>(face_normal));
                assert(!in_face_normal.degenerate());
                assert(close(in_face_normal.value[0], 0.0));
                assert(close(in_face_normal.value[1], 1.0));
                assert(close(in_face_normal.value[2], 0.0));

                // 3. Projecting a direction onto a parent edge.
                // Tetrahedron edge 5 is the x-axis edge from vertex 0 to 1.
                std::array<double, 3> direction = {3.0, 4.0, 5.0};
                const auto edge_direction =
                    geom::project_onto_parent_edge<double>(
                        tet, 5, std::span<const double>(direction));
                assert(!edge_direction.degenerate());
                assert(close(edge_direction.value[0], 3.0));
                assert(close(edge_direction.value[1], 0.0));
                assert(close(edge_direction.value[2], 0.0));

                // 4. Detecting a degenerate projected direction.
                std::array<double, 3> vertical = {0.0, 0.0, 7.0};
                const auto admissible =
                    geom::admissible_direction_in_parent_frame<double>(
                        tet,
                        geom::ParentEntity{2, 3},
                        std::span<const double>(vertical));
                assert(admissible.degenerate());
                assert(admissible.degeneracy == geom::Degeneracy::zero_projection);

                return 0;
            }
            """
        )
    )

    compile_result = subprocess.run(
        [
            *compiler,
            "-std=c++20",
            "-I",
            str(repo / "cpp/src"),
            str(source),
            "-o",
            str(exe),
        ],
        capture_output=True,
        text=True,
    )
    assert compile_result.returncode == 0, compile_result.stderr

    run_result = subprocess.run([str(exe)], capture_output=True, text=True)
    assert run_result.returncode == 0, run_result.stderr
