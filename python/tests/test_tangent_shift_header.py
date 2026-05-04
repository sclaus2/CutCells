import os
import shlex
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


def test_tangent_shift_header_reports(tmp_path):
    cxx = os.environ.get("CXX")
    if cxx:
        compiler = shlex.split(cxx)
    else:
        path = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if path is None:
            pytest.skip("C++ compiler not available")
        compiler = [path]

    repo = Path(__file__).resolve().parents[2]
    source = tmp_path / "tangent_shift_check.cpp"
    exe = tmp_path / "tangent_shift_check"
    source.write_text(
        textwrap.dedent(
            r"""
            #include "tangent_shift.h"

            #include <array>
            #include <cassert>
            #include <cmath>
            #include <span>
            #include <vector>

            namespace cell = cutcells::cell;
            namespace geom = cutcells::geom;
            namespace ts = cutcells::tangent_shift;

            bool close(double a, double b, double tol = 1.0e-8)
            {
                return std::fabs(a - b) <= tol;
            }

            std::span<const double> sp(const std::vector<double>& v)
            {
                return std::span<const double>(v.data(), v.size());
            }

            std::array<double, 3> circle_point(double fraction)
            {
                const double x = 0.2 + 0.6 * fraction;
                const double dx = x - 0.5;
                const double y = -0.5 + std::sqrt(0.36 - dx * dx);
                return {x, y, 0.0};
            }

            std::vector<double> make_points(
                std::initializer_list<std::array<double, 3>> points)
            {
                std::vector<double> out;
                for (const auto& point : points)
                    out.insert(out.end(), point.begin(), point.end());
                return out;
            }

            int main()
            {
                const auto tet = cell::type::tetrahedron;
                const geom::ParentEntity z_face{2, 3};
                const auto a = circle_point(0.0);
                const auto b = circle_point(1.0);
                const auto midpoint = circle_point(0.5);
                const auto drifted = circle_point(0.70);

                const std::vector<double> provisional = make_points({a, b});
                const std::vector<double> desired_mid = {0.0, 0.5, 1.0};

                ts::Options<double> opts;
                opts.max_relative_arclength_drift = 0.02;
                opts.max_absolute_arclength_drift = 1.0e-12;
                opts.max_relative_final_arclength_error = 0.03;
                opts.max_absolute_final_arclength_error = 1.0e-10;
                opts.min_relative_spacing = 1.0e-8;
                opts.phi_tolerance = 1.0e-10;
                opts.max_correction_iterations = 24;

                auto circle_phi = [](std::span<const double> x) -> double
                {
                    const double dx = x[0] - 0.5;
                    const double dy = x[1] + 0.5;
                    return dx * dx + dy * dy - 0.36;
                };
                auto circle_grad = [](std::span<const double> x,
                                      std::span<double> g)
                {
                    g[0] = 2.0 * (x[0] - 0.5);
                    g[1] = 2.0 * (x[1] + 0.5);
                    g[2] = 0.0;
                };

                // 1. Projected midpoint already has acceptable arclength drift.
                {
                    const std::vector<double> edge = make_points({a, midpoint, b});
                    const auto report =
                        ts::correct_projected_edge_node<double>(
                            tet, z_face, sp(provisional), sp(edge), 3,
                            sp(desired_mid), 1, circle_phi, circle_grad, opts);
                    assert(report.accepted);
                    assert(!report.corrected);
                    assert(report.failure_reason == ts::FailureReason::none);
                    assert(report.metrics.initial_arclength_error < 1.0e-12);
                }

                // 2. Excessive midpoint drift is shifted and re-corrected to phi = 0.
                {
                    const std::vector<double> edge = make_points({a, drifted, b});
                    const auto report =
                        ts::correct_projected_edge_node<double>(
                            tet, z_face, sp(provisional), sp(edge), 3,
                            sp(desired_mid), 1, circle_phi, circle_grad, opts);
                    assert(report.accepted);
                    assert(report.corrected);
                    assert(close(report.corrected_node[0], midpoint[0]));
                    assert(close(report.corrected_node[1], midpoint[1]));
                    assert(report.metrics.phi_residual <= opts.phi_tolerance);
                }

                // 3. If the correction would leave the admissible parent face, reject.
                {
                    auto outside_phi = [](std::span<const double> x) -> double
                    {
                        return x[0] + x[1] - 1.2;
                    };
                    auto outside_grad = [](std::span<const double>,
                                           std::span<double> g)
                    {
                        g[0] = 1.0;
                        g[1] = 1.0;
                        g[2] = 0.0;
                    };
                    const std::vector<double> edge = make_points({a, drifted, b});
                    auto fail_opts = opts;
                    fail_opts.max_correction_iterations = 8;
                    const auto report =
                        ts::correct_projected_edge_node<double>(
                            tet, z_face, sp(provisional), sp(edge), 3,
                            sp(desired_mid), 1, outside_phi, outside_grad, fail_opts);
                    assert(!report.accepted);
                    assert(report.request_refinement);
                    assert(report.failure_reason
                        == ts::FailureReason::correction_left_parent_entity);
                }

                // 4. A correction that crosses the next edge node destroys ordering.
                {
                    const auto before = circle_point(0.35);
                    const auto next = circle_point(0.45);
                    const std::vector<double> edge = make_points({a, before, next, b});
                    const std::vector<double> desired = {0.0, 0.5, 0.75, 1.0};
                    const auto report =
                        ts::correct_projected_edge_node<double>(
                            tet, z_face, sp(provisional), sp(edge), 3,
                            sp(desired), 1, circle_phi, circle_grad, opts);
                    assert(!report.accepted);
                    assert(report.request_refinement);
                    assert(report.failure_reason
                        == ts::FailureReason::ordering_destroyed);
                }

                // 5. If re-correction cannot satisfy phi = 0, reject.
                {
                    const std::vector<double> edge = make_points({a, drifted, b});
                    auto no_iter = opts;
                    no_iter.max_correction_iterations = 0;
                    const auto report =
                        ts::correct_projected_edge_node<double>(
                            tet, z_face, sp(provisional), sp(edge), 3,
                            sp(desired_mid), 1, circle_phi, circle_grad, no_iter);
                    assert(!report.accepted);
                    assert(report.request_refinement);
                    assert(report.failure_reason
                        == ts::FailureReason::phi_residual_too_large);
                }

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
