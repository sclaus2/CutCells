import os
import shlex
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


def test_graph_criteria_header_reports(tmp_path):
    cxx = os.environ.get("CXX")
    if cxx:
        compiler = shlex.split(cxx)
    else:
        path = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if path is None:
            pytest.skip("C++ compiler not available")
        compiler = [path]

    repo = Path(__file__).resolve().parents[2]
    source = tmp_path / "graph_criteria_check.cpp"
    exe = tmp_path / "graph_criteria_check"
    source.write_text(
        textwrap.dedent(
            r"""
            #include "graph_criteria.h"

            #include <array>
            #include <cassert>
            #include <span>
            #include <vector>

            namespace cell = cutcells::cell;
            namespace geom = cutcells::geom;
            namespace graph = cutcells::graph_criteria;

            template <std::size_t N>
            std::span<const double> sp(const std::array<double, N>& a)
            {
                return std::span<const double>(a.data(), a.size());
            }

            int main()
            {
                const auto tet = cell::type::tetrahedron;
                const geom::ParentEntity z_face{2, 3};

                graph::Options<double> opts;
                opts.min_restricted_gradient_strength = 1.0e-12;
                opts.min_transversality = 1.0e-6;
                opts.max_relative_correction_distance = 1.0;
                opts.max_relative_tangential_shift = 0.25;

                graph::HostFrame<double> edge_host;
                edge_host.dimension = graph::HostDimension::edge;
                edge_host.normal = {1.0, 0.0, 0.0};
                edge_host.tangent = {0.0, 1.0, 0.0};
                edge_host.h = 1.0;

                // 1. Candidate direction tangent to the true zero set.
                {
                    const std::array<double, 3> xh = {0.2, 0.2, 0.0};
                    const std::array<double, 3> grad = {0.0, 1.0, 0.0};
                    const std::array<double, 3> d = {1.0, 0.0, 0.0};
                    const std::array<double, 3> xc = {0.25, 0.2, 0.0};
                    const auto report = graph::evaluate_direction<double>(
                        tet, z_face, edge_host, sp(xh), sp(grad), sp(d), sp(xc),
                        graph::DirectionKind::projected_level_set_gradient, opts);
                    assert(!report.accepted);
                    assert(report.failure_reason == graph::FailureReason::tangent_to_zero_set);
                }

                // 2. Candidate direction too tangential to the straight host normal.
                {
                    graph::Options<double> drift_opts = opts;
                    drift_opts.min_host_normal_alignment = 0.0;
                    drift_opts.max_drift_amplification = 4.0;
                    const std::array<double, 3> xh = {0.2, 0.2, 0.0};
                    const std::array<double, 3> grad = {1.0, 0.01, 0.0};
                    const std::array<double, 3> d = {1.0, 0.01, 0.0};
                    const std::array<double, 3> xc = {0.25, 0.2005, 0.0};
                    graph::HostFrame<double> host = edge_host;
                    host.normal = {0.0, 1.0, 0.0};
                    host.tangent = {1.0, 0.0, 0.0};
                    const auto report = graph::evaluate_direction<double>(
                        tet, z_face, host, sp(xh), sp(grad), sp(d), sp(xc),
                        graph::DirectionKind::projected_level_set_gradient, drift_opts);
                    assert(!report.accepted);
                    assert(report.failure_reason
                        == graph::FailureReason::excessive_drift_amplification);
                }

                // 3. Root search segment leaves the admissible parent face.
                {
                    const std::array<double, 3> xh = {0.9, 0.05, 0.0};
                    const std::array<double, 3> grad = {1.0, 0.0, 0.0};
                    const std::array<double, 3> d = {1.0, 0.0, 0.0};
                    const std::array<double, 3> xc = {1.1, 0.05, 0.0};
                    const auto report = graph::evaluate_direction<double>(
                        tet, z_face, edge_host, sp(xh), sp(grad), sp(d), sp(xc),
                        graph::DirectionKind::projected_straight_host_normal, opts);
                    assert(!report.accepted);
                    assert(report.failure_reason
                        == graph::FailureReason::root_segment_leaves_parent_entity);
                }

                // 4. Parent-face edge using the in-face segment normal is valid.
                {
                    const std::array<double, 3> a = {0.25, 0.25, 0.0};
                    const std::array<double, 3> b = {0.75, 0.25, 0.0};
                    const auto face_normal =
                        geom::parent_face_normal<double>(tet, 3, true);
                    const auto in_face_normal =
                        geom::in_face_segment_normal<double>(
                            sp(a), sp(b),
                            std::span<const double>(
                                face_normal.value.data(), face_normal.value.size()));
                    const auto tangent =
                        geom::segment_tangent<double>(sp(a), sp(b), true);
                    graph::HostFrame<double> host;
                    host.dimension = graph::HostDimension::edge;
                    host.normal = in_face_normal.value;
                    host.tangent = tangent.value;
                    host.h = 0.5;
                    const std::array<double, 3> xh = {0.4, 0.25, 0.0};
                    const std::array<double, 3> grad = {0.0, 1.0, 0.0};
                    const std::array<double, 3> d = {0.0, 1.0, 0.0};
                    const std::array<double, 3> xc = {0.4, 0.30, 0.0};
                    const auto report = graph::evaluate_direction<double>(
                        tet, z_face, host, sp(xh), sp(grad), sp(d), sp(xc),
                        graph::DirectionKind::projected_straight_host_normal, opts);
                    assert(report.accepted);
                    assert(report.failure_reason == graph::FailureReason::none);
                }

                // 5. Excessive tangential shift is rejected even when the root lies on the ray.
                {
                    graph::Options<double> shift_opts = opts;
                    shift_opts.max_drift_amplification = 10.0;
                    shift_opts.max_relative_tangential_shift = 0.1;
                    graph::HostFrame<double> host = edge_host;
                    host.normal = {0.0, 1.0, 0.0};
                    host.tangent = {1.0, 0.0, 0.0};
                    const std::array<double, 3> xh = {0.1, 0.1, 0.0};
                    const std::array<double, 3> grad = {0.8, 1.0, 0.0};
                    const std::array<double, 3> d = {0.8, 1.0, 0.0};
                    const std::array<double, 3> xc = {0.26, 0.3, 0.0};
                    const auto report = graph::evaluate_direction<double>(
                        tet, z_face, host, sp(xh), sp(grad), sp(d), sp(xc),
                        graph::DirectionKind::projected_level_set_gradient, shift_opts);
                    assert(!report.accepted);
                    assert(report.failure_reason
                        == graph::FailureReason::excessive_tangential_shift);
                }

                // Preferred direction rule: accept the straight-host normal first.
                {
                    graph::HostFrame<double> host = edge_host;
                    host.normal = {0.0, 1.0, 0.0};
                    host.tangent = {1.0, 0.0, 0.0};
                    const std::array<double, 3> xh = {0.2, 0.2, 0.0};
                    const std::array<double, 3> grad = {0.0, 1.0, 0.0};
                    const std::array<double, 3> grad_d = {0.2, 1.0, 0.0};
                    const std::array<double, 3> grad_root = {0.21, 0.25, 0.0};
                    const std::array<double, 3> normal_d = {0.0, 1.0, 0.0};
                    const std::array<double, 3> normal_root = {0.2, 0.25, 0.0};
                    const auto selected = graph::select_preferred_direction<double>(
                        tet, z_face, host, sp(xh), sp(grad), sp(grad_d), sp(grad_root),
                        sp(normal_d), sp(normal_root), opts);
                    assert(selected.accepted);
                    assert(selected.selected_kind
                        == graph::DirectionKind::projected_straight_host_normal);
                }

                // Ordering and face quality helpers remain pure accept/reject checks.
                {
                    const std::array<double, 6> host_pts = {0.0, 0.0, 0.5, 0.0, 1.0, 0.0};
                    const std::array<double, 6> ok_pts = {0.0, 0.0, 0.5, 0.0, 1.0, 0.0};
                    const std::array<double, 6> folded_pts = {0.0, 0.0, 0.6, 0.0, 0.55, 0.0};
                    const std::array<double, 2> tangent = {1.0, 0.0};
                    const auto ok = graph::evaluate_projected_edge_ordering<double>(
                        sp(host_pts), sp(ok_pts), 2, sp(tangent), 1.0, opts);
                    const auto folded = graph::evaluate_projected_edge_ordering<double>(
                        sp(host_pts), sp(folded_pts), 2, sp(tangent), 1.0, opts);
                    assert(ok.accepted);
                    assert(!folded.accepted);
                    assert(folded.failure_reason == graph::FailureReason::edge_ordering_fold);

                    assert(graph::failure_reason_name(
                        graph::FailureReason::surface_jacobian_not_positive)
                        == "surface_jacobian_not_positive");
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
