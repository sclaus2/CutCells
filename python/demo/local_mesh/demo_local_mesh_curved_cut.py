import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

import cutcells

_COMPARE_METHODS = (
    cutcells.EdgeRootMethod.brent,
    cutcells.EdgeRootMethod.itp,
    cutcells.EdgeRootMethod.newton,
)


def _domain_color(domain_id: int) -> str:
    # cell::domain: inside=0, intersected=1, outside=2
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


def demo_triangle_curved(root_method=cutcells.EdgeRootMethod.brent):
    order = 3
    ct = cutcells.CellType.triangle
    tpl = cutcells.iso_p1_template(ct, order)
    parent_cell_coords = np.asarray(cutcells.iso_p1_ref_coords(ct, order), dtype=np.float64)

    lm = cutcells.init_local_mesh_from_template(
        tpl, parent_cell_coords, ct, parent_cell_id=0, n_level_sets=1
    )

    # Curved interface: circle
    ls = cutcells.LevelSetFunction(
        value=lambda x, cell_id=-1: (x[0] - 0.40) ** 2 + (x[1] - 0.35) ** 2 - 0.13**2,
        gdim=2,
    )

    cutcells.compute_roots_on_local_mesh(
        lm, ls, level_set_id=0, root_method=root_method, tol=1e-13
    )
    edge_state = np.asarray(lm.edge_state, dtype=np.uint8).flatten()
    root_param = np.asarray(lm.edge_root_parameter, dtype=np.float64).flatten()
    root_iter = np.asarray(lm.edge_root_iterations, dtype=np.int32).flatten()
    root_eval = np.asarray(lm.edge_root_evaluations, dtype=np.int32).flatten()
    root_conv = np.asarray(lm.edge_root_converged, dtype=np.uint8).flatten()
    root_res = np.asarray(lm.edge_root_residual, dtype=np.float64).flatten()
    one_root_edges = np.flatnonzero(edge_state == 1)
    if one_root_edges.size:
        print(f"[triangle] root_method={root_method.name}, one_root_edges={one_root_edges.size}")
        for e in one_root_edges:
            print(
                f"  edge {e:3d}: t={root_param[e]:.6f}, it={root_iter[e]:2d}, "
                f"eval={root_eval[e]:2d}, converged={bool(root_conv[e])}, residual={root_res[e]:.3e}"
            )

    cutcells.decompose_local_mesh(
        lm,
        ls,
        level_set_id=0,
        root_method=root_method,
        triangulate=True,
        tol=1e-13,
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

    # Mark actual computed root vertices (these should lie on phi=0).
    root_vids = np.unique(edge_root_vertex[edge_root_vertex >= 0])
    if root_vids.size:
        xr = x[root_vids]
        ax.scatter(
            xr[:, 0], xr[:, 1], s=30, c="crimson", marker="o",
            zorder=4, label="computed roots"
        )
        phi_at_roots = (xr[:, 0] - 0.40) ** 2 + (xr[:, 1] - 0.35) ** 2 - 0.13**2
        print(
            f"[triangle] roots={root_vids.size}, "
            f"max|phi(root)|={np.max(np.abs(phi_at_roots)):.3e}"
        )
    else:
        phi_at_roots = np.zeros((0,), dtype=np.float64)

    # Draw only the actual inside/outside interface from decomposed cells.
    iface_edges = _collect_triangle_interface_edges(cells, offsets, cell_domain)
    if iface_edges:
        segs = np.array([[x[e0, :2], x[e1, :2]] for e0, e1 in iface_edges], dtype=np.float64)
        ax.add_collection(
            LineCollection(segs, colors="tab:green", linewidths=2.0, zorder=5, label="piecewise interface")
        )

    # Exact phi=0 contour in the reference triangle, for visual comparison.
    n = 250
    gx = np.linspace(0.0, 1.0, n)
    gy = np.linspace(0.0, 1.0, n)
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    Phi = (X - 0.40) ** 2 + (Y - 0.35) ** 2 - 0.13**2
    Phi[X + Y > 1.0] = np.nan
    ax.contour(X, Y, Phi, levels=[0.0], colors="crimson", linewidths=1.6)

    ax.set_title(f"Curved Triangle Cut (root={root_method.name})")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper right")
    fig.tight_layout()
    plt.show()
    n_conv = int(np.count_nonzero(root_conv[one_root_edges])) if one_root_edges.size else 0
    return {
        "method": root_method.name,
        "one_root_edges": int(one_root_edges.size),
        "converged": n_conv,
        "avg_it": float(root_iter[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "avg_eval": float(root_eval[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "max_res": float(root_res[one_root_edges].max()) if one_root_edges.size else 0.0,
        "max_abs_phi_root": float(np.max(np.abs(phi_at_roots))) if phi_at_roots.size else 0.0,
    }


def demo_tetrahedron_curved(root_method=cutcells.EdgeRootMethod.itp):
    order = 3
    ct = cutcells.CellType.tetrahedron
    tpl = cutcells.iso_p1_template(ct, order)
    parent_cell_coords = np.asarray(cutcells.iso_p1_ref_coords(ct, order), dtype=np.float64)

    lm = cutcells.init_local_mesh_from_template(
        tpl, parent_cell_coords, ct, parent_cell_id=1, n_level_sets=1
    )

    # Curved interface: sphere
    ls = cutcells.LevelSetFunction(
        value=lambda x, cell_id=-1: (x[0] - 0.35) ** 2
        + (x[1] - 0.30) ** 2
        + (x[2] - 0.28) ** 2
        - 0.23**2,
        gdim=3,
    )

    cutcells.compute_roots_on_local_mesh(
        lm, ls, level_set_id=0, root_method=root_method, tol=1e-13
    )
    edge_state = np.asarray(lm.edge_state, dtype=np.uint8).flatten()
    root_iter = np.asarray(lm.edge_root_iterations, dtype=np.int32).flatten()
    root_eval = np.asarray(lm.edge_root_evaluations, dtype=np.int32).flatten()
    root_res = np.asarray(lm.edge_root_residual, dtype=np.float64).flatten()
    root_conv = np.asarray(lm.edge_root_converged, dtype=np.uint8).flatten()
    root_param = np.asarray(lm.edge_root_parameter, dtype=np.float64).flatten()
    one_root_edges = np.flatnonzero(edge_state == 1)
    if one_root_edges.size:
        n_conv = int(np.count_nonzero(root_conv[one_root_edges]))
        n_total = int(one_root_edges.size)
        conv_rate = 100.0 * n_conv / max(1, n_total)
        print(
            f"[tetrahedron] root_method={root_method.name}, one_root_edges={n_total}, "
            f"converged={n_conv}/{n_total} ({conv_rate:.1f}%), "
            f"avg_it={root_iter[one_root_edges].mean():.2f}, "
            f"avg_eval={root_eval[one_root_edges].mean():.2f}, "
            f"max_res={root_res[one_root_edges].max():.3e}"
        )
        for e in one_root_edges:
            print(
                f"  edge {int(e):3d}: t={float(root_param[e]):.6f}, it={int(root_iter[e]):2d}, "
                f"eval={int(root_eval[e]):2d}, converged={bool(root_conv[e])}, residual={float(root_res[e]):.3e}"
            )

    cutcells.decompose_local_mesh(
        lm,
        ls,
        level_set_id=0,
        root_method=root_method,
        triangulate=True,
        tol=1e-13,
    )

    x = np.asarray(lm.vertex_x, dtype=np.float64).reshape(-1, lm.gdim)
    offsets = np.asarray(lm.cell_offsets, dtype=np.int32)
    cells = np.asarray(lm.cell_vertices, dtype=np.int32)
    cell_domain = np.asarray(lm.cell_domain, dtype=np.uint8)

    inside_segments = []
    outside_segments = []
    interface_tris = []   # triangular facets with domain==1

    for i in range(lm.n_cells()):
        c0, c1 = offsets[i], offsets[i + 1]
        verts_idx = cells[c0:c1]
        dom = int(cell_domain[i])
        if verts_idx.size == 4:
            # Tetrahedral cell — draw wireframe edges
            verts = x[verts_idx]
            segs = [
                [verts[0], verts[1]],
                [verts[0], verts[2]],
                [verts[0], verts[3]],
                [verts[1], verts[2]],
                [verts[1], verts[3]],
                [verts[2], verts[3]],
            ]
            if dom == 0:
                inside_segments.extend(segs)
            elif dom == 2:
                outside_segments.extend(segs)
        elif verts_idx.size == 3:
            # Triangular cell — collect interface facets (domain==1)
            if dom == 1:
                interface_tris.append(x[verts_idx])

    # ---------------------------------------------------------------
    # Figure
    # ---------------------------------------------------------------
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Inside / outside wireframe
    if inside_segments:
        ax.add_collection3d(
            Line3DCollection(inside_segments, colors="tab:blue", linewidths=0.6, alpha=0.55,
                             label="inside tets")
        )
    if outside_segments:
        ax.add_collection3d(
            Line3DCollection(outside_segments, colors="tab:orange", linewidths=0.6, alpha=0.35,
                             label="outside tets")
        )

    # Piecewise-linear interface patches
    if interface_tris:
        pc = Poly3DCollection(
            interface_tris,
            facecolor="tab:green",
            edgecolor="darkgreen",
            linewidths=0.8,
            alpha=0.70,
            label="PL interface",
            zorder=4,
        )
        ax.add_collection3d(pc)

    # Computed root vertices
    ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=8, c="black", depthshade=False, zorder=3)
    root_vids_tet = np.unique(np.asarray(lm.edge_root_vertex, dtype=np.int32).flatten())
    root_vids_tet = root_vids_tet[root_vids_tet >= 0]
    if root_vids_tet.size:
        xr = x[root_vids_tet]
        ax.scatter(xr[:, 0], xr[:, 1], xr[:, 2], s=25, c="crimson", marker="o",
                   depthshade=False, zorder=5, label="computed roots")

    # Exact sphere surface for comparison
    cx, cy, cz, radius = 0.35, 0.30, 0.28, 0.23
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    sx = cx + radius * np.outer(np.cos(u), np.sin(v))
    sy = cy + radius * np.outer(np.sin(u), np.sin(v))
    sz = cz + radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(sx, sy, sz, color="crimson", alpha=0.15, linewidth=0,
                    antialiased=True, label="exact sphere")
    # Add a proxy artist for the legend (plot_surface legend support is limited)
    from matplotlib.patches import Patch
    legend_proxy = [
        Patch(facecolor="tab:blue", alpha=0.55, label="inside tets"),
        Patch(facecolor="tab:orange", alpha=0.55, label="outside tets"),
        Patch(facecolor="tab:green", alpha=0.70, label="PL interface"),
        Patch(facecolor="crimson", alpha=0.35, label="exact sphere"),
    ]

    ax.set_title(f"Curved Tetrahedron Cut (root={root_method.name})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect((1, 1, 1))
    ax.legend(handles=legend_proxy, loc="upper left", fontsize=8)
    fig.tight_layout()
    plt.show()
    n_conv = int(np.count_nonzero(root_conv[one_root_edges])) if one_root_edges.size else 0
    return {
        "method": root_method.name,
        "one_root_edges": int(one_root_edges.size),
        "converged": n_conv,
        "avg_it": float(root_iter[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "avg_eval": float(root_eval[one_root_edges].mean()) if one_root_edges.size else 0.0,
        "max_res": float(root_res[one_root_edges].max()) if one_root_edges.size else 0.0,
    }


def _print_compare_table(shape_name, stats):
    print(f"\n[{shape_name}] root-method comparison")
    for s in stats:
        n = max(1, s["one_root_edges"])
        conv_pct = 100.0 * s["converged"] / n
        line = (
            f"  {s['method']:>7s}: one_root_edges={s['one_root_edges']:3d}, "
            f"converged={s['converged']:3d}/{s['one_root_edges']:3d} ({conv_pct:5.1f}%), "
            f"avg_it={s['avg_it']:.2f}, avg_eval={s['avg_eval']:.2f}, max_res={s['max_res']:.3e}"
        )
        if "max_abs_phi_root" in s:
            line += f", max|phi(root)|={s['max_abs_phi_root']:.3e}"
        print(line)


def main():
    tri_stats = [demo_triangle_curved(m) for m in _COMPARE_METHODS]
    tet_stats = [demo_tetrahedron_curved(m) for m in _COMPARE_METHODS]
    _print_compare_table("triangle", tri_stats)
    _print_compare_table("tetrahedron", tet_stats)


if __name__ == "__main__":
    main()
