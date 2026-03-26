"""
Demo 3 — Straight vs curved interface on one cell (Step 5)

Compare the straight (P1) and curved (Pk) interface on a single
background cell. Uses a circle level set on an iso-P1 refined triangle.

Shows:
- Straight sub-mesh edges (thin grey)
- Straight interface (red dashed line connecting root vertices)
- Curved interface Pk nodes (green dots at Gauss-Lobatto positions)
- Curved interface as smooth interpolated curve through the green dots
- Exact zero-contour of phi (black line)
- Subplots for P2, P3, P4 geometry orders

Verifies: Curved node placement, Gauss-Lobatto spacing, convergence with
geometry order.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import cutcells as cc


# ============================================================================
# Level-set definition: circle
# ============================================================================

CX, CY, R = 0.40, 0.35, 0.18


def phi_value(x, cell_id=-1):
    return float((x[0] - CX) ** 2 + (x[1] - CY) ** 2 - R**2)


def phi_grad(x, cell_id=-1):
    return np.array([2.0 * (x[0] - CX), 2.0 * (x[1] - CY), 0.0])


# ============================================================================
# Helpers
# ============================================================================


def collect_triangle_interface_edges(cells, offsets, cell_domain, vertex_x):
    """Find edges shared between inside (domain=0) and outside (domain=2) cells."""
    edge_to_domains = {}
    for i in range(len(offsets) - 1):
        c0, c1 = offsets[i], offsets[i + 1]
        tri = cells[c0:c1]
        if tri.size != 3:
            continue
        dom = int(cell_domain[i])
        for a, b in ((0, 1), (1, 2), (2, 0)):
            v0, v1 = int(tri[a]), int(tri[b])
            if v0 > v1:
                v0, v1 = v1, v0
            edge_to_domains.setdefault((v0, v1), []).append(dom)

    iface = []
    for e, doms in edge_to_domains.items():
        uniq = set(doms)
        if 0 in uniq and 2 in uniq:
            iface.append(e)
    return iface


def lagrange_interp_1d(nodes, values, t_fine):
    """Lagrange interpolation through (nodes[i], values[i]) at points t_fine."""
    n = len(nodes)
    result = np.zeros(len(t_fine))
    for i in range(n):
        l_i = np.ones(len(t_fine))
        for j in range(n):
            if j != i:
                l_i *= (t_fine - nodes[j]) / (nodes[i] - nodes[j])
        result += values[i] * l_i
    return result


def gauss_lobatto_nodes_1d(order):
    """Return Gauss-Lobatto nodes on [0, 1]."""
    if order <= 1:
        return np.array([0.0, 1.0], dtype=np.float64)
    P = np.polynomial.legendre.Legendre.basis(order)
    interior = np.sort(P.deriv().roots())
    return np.concatenate(([0.0], 0.5 * (interior + 1.0), [1.0]))


def plot_curved_interface_from_cache(ax, mesh, vertex_x, geom_order):
    """Plot curved interface entities directly from the per-entity cache."""
    iface_vertices = np.array(mesh.iface_vertices, dtype=np.int32)
    iface_offsets = np.array(mesh.iface_offsets, dtype=np.int32)
    curved_ref = np.array(mesh.curved_iface_ref_nodes, dtype=np.float64)
    curved_offsets = np.array(mesh.curved_iface_offsets, dtype=np.int32)
    curved_conv = np.array(mesh.curved_iface_converged, dtype=np.uint8)
    tdim = mesh.tdim

    if iface_offsets.size == 0:
        return 0.0, 0

    max_phi = 0.0
    n_nodes = 0

    for ie in range(mesh.n_iface_entities()):
        ev0 = iface_offsets[ie]
        ev1 = iface_offsets[ie + 1]
        ent_verts = iface_vertices[ev0:ev1]
        if ent_verts.size < 2:
            continue

        node_start = curved_offsets[ie] if curved_offsets.size else 0
        node_end = curved_offsets[ie + 1] if curved_offsets.size else 0

        pts = [vertex_x[int(ent_verts[0]), :2]]
        for ni in range(node_start, node_end):
            xr = curved_ref[ni * tdim : (ni + 1) * tdim]
            pts.append(xr[:2])
            n_nodes += 1
            if curved_conv.size > ni and curved_conv[ni]:
                max_phi = max(max_phi, abs(phi_value(np.array([xr[0], xr[1], 0.0]))))
        pts.append(vertex_x[int(ent_verts[-1]), :2])
        pts = np.asarray(pts, dtype=np.float64)

        if len(pts) >= 2:
            t_nodes = gauss_lobatto_nodes_1d(geom_order)
            if len(t_nodes) != len(pts):
                t_nodes = np.linspace(0.0, 1.0, len(pts))
            t_fine = np.linspace(0.0, 1.0, 200)
            x_fine = lagrange_interp_1d(t_nodes, pts[:, 0], t_fine)
            y_fine = lagrange_interp_1d(t_nodes, pts[:, 1], t_fine)
            ax.plot(x_fine, y_fine, "-", color="tab:green", lw=2.0, zorder=5)

        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            s=25,
            c="tab:green",
            marker="o",
            zorder=7,
            edgecolors="black",
            linewidths=0.5,
        )

    return max_phi, n_nodes


# ============================================================================
# Run pipeline for one geometry order
# ============================================================================


def run_single_order(geom_order, ax):
    """Run the full pipeline on one cell and plot on the given axes."""
    iso_order = max(geom_order, 3)
    ct = cc.CellType.triangle
    tpl = cc.iso_p1_template(ct, iso_order)
    parent_coords_flat = np.asarray(
        cc.iso_p1_ref_coords(ct, iso_order), dtype=np.float64
    )
    parent_coords = parent_coords_flat.reshape(-1, 2)
    nodal_values = np.asarray(
        [phi_value(np.array([x, y, 0.0])) for x, y in parent_coords[:, :2]],
        dtype=np.float64,
    )

    mesh = cc.init_local_mesh_from_template(
        tpl,
        parent_coords_flat,
        ct,
        parent_cell_id=0,
        n_level_sets=1,
    )

    ls = cc.LevelSetFunction(nodal_values=nodal_values, gdim=2, degree=iso_order)

    # Decompose using the polynomial/Bernstein path only.
    cc.decompose_local_mesh(
        mesh,
        ls,
        level_set_id=0,
        root_method=cc.EdgeRootMethod.itp,
        triangulate=True,
        tol=1e-14,
        backend=cc.LocalLevelSetBackend.bernstein,
    )

    # Build curved geometry (new API)
    cc.build_interface_entities(mesh, level_set_id=0)
    cc.curve_interface_entities(
        mesh,
        ls,
        level_set_id=0,
        geom_order=geom_order,
        backend=cc.LocalLevelSetBackend.bernstein,
        tol=1e-14,
    )

    vertex_x = np.array(mesh.vertex_x, dtype=np.float64).reshape(-1, mesh.gdim)
    cells = np.array(mesh.cell_vertices, dtype=np.int32)
    offsets = np.array(mesh.cell_offsets, dtype=np.int32)
    cell_domain = np.array(mesh.cell_domain, dtype=np.uint8)

    curved_ref = np.array(mesh.curved_iface_ref_nodes, dtype=np.float64)
    curved_conv = np.array(mesh.curved_iface_converged, dtype=np.uint8)
    tdim = mesh.tdim

    # ---- Draw reference triangle ----
    tri_border = np.array([[0, 0], [1, 0], [0, 1], [0, 0]], dtype=float)
    ax.plot(tri_border[:, 0], tri_border[:, 1], "k-", lw=1.5)

    # ---- Draw sub-mesh edges (thin grey) ----
    for i in range(mesh.n_cells()):
        c0, c1 = offsets[i], offsets[i + 1]
        tri = cells[c0:c1]
        if tri.size != 3:
            continue
        poly = vertex_x[tri, :2]
        seg = np.vstack([poly, poly[0]])
        ax.plot(seg[:, 0], seg[:, 1], color="grey", lw=0.4, alpha=0.5)

    # ---- Straight interface (red dashed) ----
    iface = collect_triangle_interface_edges(cells, offsets, cell_domain, vertex_x)
    if iface:
        segs = np.array([[vertex_x[e0, :2], vertex_x[e1, :2]] for e0, e1 in iface])
        ax.add_collection(
            LineCollection(
                segs,
                colors="red",
                linewidths=1.5,
                linestyles="dashed",
                zorder=4,
                label="straight interface",
            )
        )

    # ---- Exact zero contour (black) ----
    n_fine = 300
    gx = np.linspace(0.0, 1.0, n_fine)
    gy = np.linspace(0.0, 1.0, n_fine)
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    Phi = (X - CX) ** 2 + (Y - CY) ** 2 - R**2
    Phi[X + Y > 1.0] = np.nan
    ax.contour(X, Y, Phi, levels=[0.0], colors="black", linewidths=1.8)

    # ---- Curved Pk nodes from interface-entity cache (if any) ----
    curved_ref = np.array(mesh.curved_iface_ref_nodes, dtype=np.float64)
    curved_conv = np.array(mesh.curved_iface_converged, dtype=np.uint8)
    tdim = mesh.tdim

    for i in range(len(curved_ref) // tdim):
        xr = curved_ref[i * tdim : (i + 1) * tdim]
        converged = bool(curved_conv[i])
        if converged:
            ax.plot(
                xr[0], xr[1], "o", color="tab:green", markersize=3, zorder=6, alpha=0.4
            )
        else:
            ax.plot(
                xr[0], xr[1], "x", color="tab:orange", markersize=3, zorder=6, alpha=0.4
            )

    # ---- Max |phi| at curved_iface_ref_nodes (machine-precision check) ----
    n_iface_nodes = len(curved_ref) // tdim
    max_phi_iface = 0.0
    for i in range(n_iface_nodes):
        xr = curved_ref[i * tdim : (i + 1) * tdim]
        if curved_conv[i]:
            pv = abs(phi_value(np.array([xr[0], xr[1], 0.0])))
            max_phi_iface = max(max_phi_iface, pv)

    # ---- Curved interface directly from per-entity cache ----
    max_dist_from_zero, n_curved_nodes = plot_curved_interface_from_cache(
        ax, mesh, vertex_x, geom_order
    )

    ax.set_title(f"P{geom_order}  (iface max|ϕ|={max_phi_iface:.1e})")
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(alpha=0.2)

    return {
        "geom_order": geom_order,
        "n_curved_nodes": n_curved_nodes,
        "max_phi_at_converged": max_phi_iface,
        "fallback_count": mesh.curved_fallback_count,
        "n_interface_cells": mesh.n_iface_entities(),
    }


# ============================================================================
# Main
# ============================================================================


def main():
    orders = [2, 3, 4]
    fig, axes = plt.subplots(1, len(orders), figsize=(5 * len(orders), 5))
    if len(orders) == 1:
        axes = [axes]

    results = []
    for order, ax in zip(orders, axes):
        info = run_single_order(order, ax)
        results.append(info)

    # Add legend to first subplot
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(loc="upper right", fontsize=7)

    fig.suptitle("Demo 3: Straight vs Curved Interface on One Cell", fontsize=12)
    fig.tight_layout()
    plt.savefig("demo_single_cell_curved.png", dpi=150)
    print("Saved demo_single_cell_curved.png")
    plt.show()

    # Print convergence table
    print("\nConvergence with geometry order:")
    print(
        f"{'Order':>6s} {'Curved Nodes':>13s} {'Max|phi|':>12s} {'Fallback':>9s} {'Iface Cells':>12s}"
    )
    for r in results:
        print(
            f"{r['geom_order']:6d} {r['n_curved_nodes']:13d} "
            f"{r['max_phi_at_converged']:12.2e} {r['fallback_count']:9d} "
            f"{r['n_interface_cells']:12d}"
        )


if __name__ == "__main__":
    main()
