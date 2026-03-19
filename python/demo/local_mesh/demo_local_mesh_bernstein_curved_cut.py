import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import cutcells


def _domain_color(domain_id: int) -> str:
    if domain_id == 0:
        return "tab:blue"
    if domain_id == 2:
        return "tab:orange"
    return "tab:gray"


def _collect_triangle_interface_edges(cells, offsets, cell_domain):
    edge_to_domains = {}
    for i in range(len(offsets) - 1):
        c0, c1 = offsets[i], offsets[i + 1]
        tri = cells[c0:c1]
        if tri.size != 3:
            continue
        dom = int(cell_domain[i])
        for a, b in ((0, 1), (1, 2), (2, 0)):
            v0 = int(tri[a])
            v1 = int(tri[b])
            if v0 > v1:
                v0, v1 = v1, v0
            edge_to_domains.setdefault((v0, v1), []).append(dom)

    iface = []
    for e, doms in edge_to_domains.items():
        uniq = set(doms)
        if 0 in uniq and 2 in uniq:
            iface.append(e)
    return iface


def _triangle_phi(x):
    # Use a radius large enough to be resolved by the sampled P4 nodal field.
    return (x[0] - 0.40) ** 2 + (x[1] - 0.35) ** 2 - 0.18**2


def _tetrahedron_phi(x):
    return (
        (x[0] - 0.35) ** 2
        + (x[1] - 0.30) ** 2
        + (x[2] - 0.28) ** 2
        - 0.23**2
    )


def _sample_nodal_values(parent_cell_coords, gdim, phi):
    x = np.asarray(parent_cell_coords, dtype=np.float64).reshape(-1, gdim)
    return np.asarray([phi(xi) for xi in x], dtype=np.float64)


def _print_root_stats(
    shape_name,
    order,
    edge_state,
    root_param,
    root_iter,
    root_eval,
    root_conv,
    root_res,
    max_preview=12,
):
    one_root_edges = np.flatnonzero(edge_state == 1)
    if not one_root_edges.size:
        print(f"[{shape_name}] order={order}, no one-root edges")
        return one_root_edges

    n_conv = int(np.count_nonzero(root_conv[one_root_edges]))
    print(
        f"[{shape_name}] order={order}, backend=bernstein, "
        f"one_root_edges={one_root_edges.size}, "
        f"converged={n_conv}/{one_root_edges.size}, "
        f"avg_it={root_iter[one_root_edges].mean():.2f}, "
        f"avg_eval={root_eval[one_root_edges].mean():.2f}, "
        f"max_res={root_res[one_root_edges].max():.3e}"
    )
    preview = one_root_edges[:max_preview]
    for e in preview:
        print(
            f"  edge {int(e):3d}: t={root_param[e]:.6f}, it={root_iter[e]:2d}, "
            f"eval={root_eval[e]:2d}, converged={bool(root_conv[e])}, residual={root_res[e]:.3e}"
        )
    if one_root_edges.size > preview.size:
        print(f"  ... {one_root_edges.size - preview.size} more one-root edges omitted")
    return one_root_edges


def demo_triangle_bernstein(order=4):
    ct = cutcells.CellType.triangle
    tpl = cutcells.iso_p1_template(ct, order)
    parent_cell_coords = np.asarray(cutcells.iso_p1_ref_coords(ct, order), dtype=np.float64)

    lm = cutcells.init_local_mesh_from_template(
        tpl, parent_cell_coords, ct, parent_cell_id=0, n_level_sets=1
    )

    nodal_values = _sample_nodal_values(parent_cell_coords, 2, _triangle_phi)
    ls = cutcells.LevelSetFunction(nodal_values=nodal_values, gdim=2, degree=order)

    cutcells.compute_roots_on_local_mesh(
        lm,
        ls,
        level_set_id=0,
        root_method=cutcells.EdgeRootMethod.linear,
        tol=1e-13,
        backend=cutcells.LocalLevelSetBackend.bernstein,
    )

    edge_state = np.asarray(lm.edge_state, dtype=np.uint8).copy()
    root_param = np.asarray(lm.edge_root_parameter, dtype=np.float64).copy()
    root_iter = np.asarray(lm.edge_root_iterations, dtype=np.int32).copy()
    root_eval = np.asarray(lm.edge_root_evaluations, dtype=np.int32).copy()
    root_conv = np.asarray(lm.edge_root_converged, dtype=np.uint8).copy()
    root_res = np.asarray(lm.edge_root_residual, dtype=np.float64).copy()
    one_root_edges = _print_root_stats(
        "triangle", order, edge_state, root_param, root_iter, root_eval, root_conv, root_res
    )

    cutcells.decompose_local_mesh(
        lm,
        ls,
        level_set_id=0,
        root_method=cutcells.EdgeRootMethod.linear,
        triangulate=True,
        tol=1e-13,
        backend=cutcells.LocalLevelSetBackend.bernstein,
    )

    x = np.asarray(lm.vertex_x, dtype=np.float64).reshape(-1, lm.gdim)
    offsets = np.asarray(lm.cell_offsets, dtype=np.int32)
    cells = np.asarray(lm.cell_vertices, dtype=np.int32)
    cell_domain = np.asarray(lm.cell_domain, dtype=np.uint8)
    edge_root_vertex = np.asarray(lm.edge_root_vertex, dtype=np.int32)

    fig, ax = plt.subplots(figsize=(6.5, 6))

    for i in range(lm.n_cells()):
        c0, c1 = offsets[i], offsets[i + 1]
        tri = cells[c0:c1]
        if tri.size != 3:
            continue
        poly = x[tri, :2]
        col = _domain_color(int(cell_domain[i]))
        ax.fill(poly[:, 0], poly[:, 1], alpha=0.35, color=col)
        seg = np.vstack([poly, poly[0]])
        ax.plot(seg[:, 0], seg[:, 1], color="black", lw=0.8)

    ax.scatter(x[:, 0], x[:, 1], s=10, c="black", zorder=3)

    root_vids = np.unique(edge_root_vertex[edge_root_vertex >= 0])
    if root_vids.size:
        xr = x[root_vids]
        ax.scatter(
            xr[:, 0],
            xr[:, 1],
            s=30,
            c="crimson",
            marker="o",
            zorder=4,
            label="Bernstein roots",
        )
        phi_at_roots = np.asarray([_triangle_phi(xi) for xi in xr], dtype=np.float64)
        print(
            f"[triangle] Bernstein roots={root_vids.size}, "
            f"max|phi_exact(root)|={np.max(np.abs(phi_at_roots)):.3e}"
        )

    iface_edges = _collect_triangle_interface_edges(cells, offsets, cell_domain)
    if iface_edges:
        segs = np.array([[x[e0, :2], x[e1, :2]] for e0, e1 in iface_edges], dtype=np.float64)
        ax.add_collection(
            LineCollection(
                segs,
                colors="tab:green",
                linewidths=2.0,
                zorder=5,
                label="Bernstein interface",
            )
        )

    n = 250
    gx = np.linspace(0.0, 1.0, n)
    gy = np.linspace(0.0, 1.0, n)
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    Phi = (X - 0.40) ** 2 + (Y - 0.35) ** 2 - 0.18**2
    Phi[X + Y > 1.0] = np.nan
    ax.contour(X, Y, Phi, levels=[0.0], colors="crimson", linewidths=1.6)

    ax.set_title(f"Triangle Curved Cut via Bernstein Backend (P{order})")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper right")
    fig.tight_layout()
    plt.show()

    return {
        "order": order,
        "one_root_edges": int(one_root_edges.size),
        "converged": int(np.count_nonzero(root_conv[one_root_edges])) if one_root_edges.size else 0,
        "avg_it": float(root_iter[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "avg_eval": float(root_eval[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "max_res": float(root_res[one_root_edges].max()) if one_root_edges.size else 0.0,
    }


def demo_tetrahedron_bernstein(order=4):
    ct = cutcells.CellType.tetrahedron
    tpl = cutcells.iso_p1_template(ct, order)
    parent_cell_coords = np.asarray(cutcells.iso_p1_ref_coords(ct, order), dtype=np.float64)

    lm = cutcells.init_local_mesh_from_template(
        tpl, parent_cell_coords, ct, parent_cell_id=1, n_level_sets=1
    )

    nodal_values = _sample_nodal_values(parent_cell_coords, 3, _tetrahedron_phi)
    ls = cutcells.LevelSetFunction(nodal_values=nodal_values, gdim=3, degree=order)

    cutcells.compute_roots_on_local_mesh(
        lm,
        ls,
        level_set_id=0,
        root_method=cutcells.EdgeRootMethod.linear,
        tol=1e-13,
        backend=cutcells.LocalLevelSetBackend.bernstein,
    )

    edge_state = np.asarray(lm.edge_state, dtype=np.uint8).copy()
    root_iter = np.asarray(lm.edge_root_iterations, dtype=np.int32).copy()
    root_eval = np.asarray(lm.edge_root_evaluations, dtype=np.int32).copy()
    root_res = np.asarray(lm.edge_root_residual, dtype=np.float64).copy()
    root_conv = np.asarray(lm.edge_root_converged, dtype=np.uint8).copy()
    root_param = np.asarray(lm.edge_root_parameter, dtype=np.float64).copy()
    one_root_edges = _print_root_stats(
        "tetrahedron", order, edge_state, root_param, root_iter, root_eval, root_conv, root_res
    )

    cutcells.decompose_local_mesh(
        lm,
        ls,
        level_set_id=0,
        root_method=cutcells.EdgeRootMethod.linear,
        triangulate=True,
        tol=1e-13,
        backend=cutcells.LocalLevelSetBackend.bernstein,
    )

    x = np.asarray(lm.vertex_x, dtype=np.float64).reshape(-1, lm.gdim)
    offsets = np.asarray(lm.cell_offsets, dtype=np.int32)
    cells = np.asarray(lm.cell_vertices, dtype=np.int32)
    cell_domain = np.asarray(lm.cell_domain, dtype=np.uint8)
    edge_root_vertex = np.asarray(lm.edge_root_vertex, dtype=np.int32)

    inside_segments = []
    outside_segments = []
    for i in range(lm.n_cells()):
        c0, c1 = offsets[i], offsets[i + 1]
        tet = cells[c0:c1]
        if tet.size != 4:
            continue
        verts = x[tet]
        segs = [
            [verts[0], verts[1]],
            [verts[0], verts[2]],
            [verts[0], verts[3]],
            [verts[1], verts[2]],
            [verts[1], verts[3]],
            [verts[2], verts[3]],
        ]
        if int(cell_domain[i]) == 0:
            inside_segments.extend(segs)
        elif int(cell_domain[i]) == 2:
            outside_segments.extend(segs)

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    if inside_segments:
        ax.add_collection3d(
            Line3DCollection(inside_segments, colors="tab:blue", linewidths=0.7, alpha=0.8)
        )
    if outside_segments:
        ax.add_collection3d(
            Line3DCollection(outside_segments, colors="tab:orange", linewidths=0.7, alpha=0.6)
        )

    ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=8, c="black", depthshade=False)

    root_vids = np.unique(edge_root_vertex[edge_root_vertex >= 0])
    if root_vids.size:
        xr = x[root_vids]
        ax.scatter(
            xr[:, 0], xr[:, 1], xr[:, 2],
            s=18, c="crimson", depthshade=False, label="Bernstein roots"
        )
        phi_at_roots = np.asarray([_tetrahedron_phi(xi) for xi in xr], dtype=np.float64)
        print(
            f"[tetrahedron] Bernstein roots={root_vids.size}, "
            f"max|phi_exact(root)|={np.max(np.abs(phi_at_roots)):.3e}"
        )

    ax.set_title(f"Tetrahedron Curved Cut via Bernstein Backend (P{order})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect((1, 1, 1))
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper right")
    fig.tight_layout()
    plt.show()

    return {
        "order": order,
        "one_root_edges": int(one_root_edges.size),
        "converged": int(np.count_nonzero(root_conv[one_root_edges])) if one_root_edges.size else 0,
        "avg_it": float(root_iter[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "avg_eval": float(root_eval[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "max_res": float(root_res[one_root_edges].max()) if one_root_edges.size else 0.0,
    }


def _print_summary(shape_name, stats):
    n = max(1, stats["one_root_edges"])
    conv_pct = 100.0 * stats["converged"] / n
    print(
        f"[{shape_name}] order={stats['order']}, one_root_edges={stats['one_root_edges']}, "
        f"converged={stats['converged']}/{stats['one_root_edges']} ({conv_pct:.1f}%), "
        f"avg_it={stats['avg_it']:.2f}, avg_eval={stats['avg_eval']:.2f}, "
        f"max_res={stats['max_res']:.3e}"
    )


def main():
    tri_stats = demo_triangle_bernstein(order=4)
    tet_stats = demo_tetrahedron_bernstein(order=4)
    _print_summary("triangle", tri_stats)
    _print_summary("tetrahedron", tet_stats)


if __name__ == "__main__":
    main()
