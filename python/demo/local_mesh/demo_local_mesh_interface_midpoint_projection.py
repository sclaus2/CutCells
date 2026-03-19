import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import cutcells


METHODS = (
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


def _build_projection_segment(mid, normal, phi, base_half_len, max_expand=8):
    half_len = float(base_half_len)
    for _ in range(max_expand):
        p0 = mid - half_len * normal
        p1 = mid + half_len * normal
        f0 = phi(p0)
        f1 = phi(p1)
        if f0 * f1 <= 0.0:
            return p0, p1
        half_len *= 1.8
    return None


def _project_midpoints(root_method):
    order = 3
    ct = cutcells.CellType.triangle
    tpl = cutcells.iso_p1_template(ct, order)
    parent_cell_coords = np.asarray(
        cutcells.iso_p1_ref_coords(ct, order), dtype=np.float64
    )

    lm = cutcells.init_local_mesh_from_template(
        tpl, parent_cell_coords, ct, parent_cell_id=0, n_level_sets=1
    )

    def phi(x):
        return (x[0] - 0.40) ** 2 + (x[1] - 0.35) ** 2 - 0.23**2

    def grad_phi(x):
        return np.array([2.0 * (x[0] - 0.40), 2.0 * (x[1] - 0.35)], dtype=np.float64)

    ls = cutcells.LevelSetFunction(
        value=lambda x, cell_id=-1: phi(x),
        grad=lambda x, cell_id=-1: grad_phi(x),
        gdim=2,
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

    iface_edges = _collect_triangle_interface_edges(cells, offsets, cell_domain)
    if not iface_edges:
        return {
            "method": root_method.name,
            "n_segments": 0,
            "n_projected": 0,
            "avg_shift": 0.0,
            "max_abs_phi_proj": 0.0,
            "avg_it": 0.0,
            "avg_eval": 0.0,
            "segments": np.zeros((0, 2, 2), dtype=np.float64),
            "mid_linear": np.zeros((0, 2), dtype=np.float64),
            "mid_projected": np.zeros((0, 2), dtype=np.float64),
            "normal": np.zeros((0, 2), dtype=np.float64),
            "x": x,
            "cells": cells,
            "offsets": offsets,
            "cell_domain": cell_domain,
        }

    segs = np.array(
        [[x[e0, :2], x[e1, :2]] for e0, e1 in iface_edges], dtype=np.float64
    )
    mid_linear = 0.5 * (segs[:, 0, :] + segs[:, 1, :])

    projected = []
    normals = []
    shifts = []
    phi_proj = []
    iters = []
    evals = []

    for i, m in enumerate(mid_linear):
        g = grad_phi(m)
        gnorm = np.linalg.norm(g)
        if gnorm <= 1e-14:
            continue
        n = g / gnorm

        seg_len = np.linalg.norm(segs[i, 1, :] - segs[i, 0, :])
        bracket = _build_projection_segment(
            m, n, phi, base_half_len=max(1e-3, 0.75 * seg_len)
        )
        if bracket is None:
            continue
        p0, p1 = bracket

        info = cutcells.find_root_on_segment(
            p0,
            p1,
            ls,
            root_method=root_method,
            cell_id=lm.parent_cell_id,
            level=0.0,
            max_iter=64,
            xtol=1e-12,
            ftol=1e-12,
        )

        if not bool(info["converged"]):
            continue

        mp = np.asarray(info["point"], dtype=np.float64)
        projected.append(mp)
        normals.append(n)
        shifts.append(float(np.linalg.norm(mp - m)))
        phi_proj.append(abs(phi(mp)))
        iters.append(int(info["iterations"]))
        evals.append(int(info["evaluations"]))

    projected = (
        np.asarray(projected, dtype=np.float64)
        if projected
        else np.zeros((0, 2), dtype=np.float64)
    )
    normals = (
        np.asarray(normals, dtype=np.float64)
        if normals
        else np.zeros((0, 2), dtype=np.float64)
    )

    return {
        "method": root_method.name,
        "n_segments": int(segs.shape[0]),
        "n_projected": int(projected.shape[0]),
        "avg_shift": float(np.mean(shifts)) if shifts else 0.0,
        "max_abs_phi_proj": float(np.max(phi_proj)) if phi_proj else 0.0,
        "avg_it": float(np.mean(iters)) if iters else 0.0,
        "avg_eval": float(np.mean(evals)) if evals else 0.0,
        "segments": segs,
        "mid_linear": mid_linear,
        "mid_projected": projected,
        "normal": normals,
        "x": x,
        "cells": cells,
        "offsets": offsets,
        "cell_domain": cell_domain,
    }


def _plot_projection_result(data):
    fig, ax = plt.subplots(figsize=(6.8, 6.2))

    x = data["x"]
    cells = data["cells"]
    offsets = data["offsets"]
    cell_domain = data["cell_domain"]
    for i in range(len(offsets) - 1):
        c0, c1 = offsets[i], offsets[i + 1]
        tri = cells[c0:c1]
        if tri.size != 3:
            continue
        poly = x[tri, :2]
        col = _domain_color(int(cell_domain[i]))
        ax.fill(poly[:, 0], poly[:, 1], alpha=0.18, color=col)
        seg = np.vstack([poly, poly[0]])
        ax.plot(seg[:, 0], seg[:, 1], color="black", lw=0.6, alpha=0.7)

    segs = data["segments"]
    if segs.size:
        ax.add_collection(
            LineCollection(
                segs,
                colors="tab:green",
                linewidths=2.0,
                label="straight interface segments",
            )
        )

    mlin = data["mid_linear"]
    if mlin.size:
        ax.scatter(
            mlin[:, 0], mlin[:, 1], s=26, c="tab:blue", label="segment midpoints"
        )

    mproj = data["mid_projected"]
    if mproj.size:
        ax.scatter(
            mproj[:, 0], mproj[:, 1], s=28, c="crimson", label="projected midpoints"
        )

    k = min(len(mproj), len(mlin))
    if k:
        for i in range(k):
            ax.plot(
                [mlin[i, 0], mproj[i, 0]],
                [mlin[i, 1], mproj[i, 1]],
                color="gray",
                lw=0.9,
                alpha=0.8,
            )

    gx = np.linspace(0.0, 1.0, 250)
    gy = np.linspace(0.0, 1.0, 250)
    X, Y = np.meshgrid(gx, gy, indexing="xy")
    Phi = (X - 0.40) ** 2 + (Y - 0.35) ** 2 - 0.13**2
    Phi[X + Y > 1.0] = np.nan
    ax.contour(X, Y, Phi, levels=[0.0], colors="black", linewidths=1.2)

    ax.set_title(f"Interface-midpoint projection (root={data['method']})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    plt.show()


def main():
    all_stats = []
    for method in METHODS:
        data = _project_midpoints(method)
        all_stats.append(data)
        print(
            f"[{data['method']}] intervals={data['n_segments']}, projected={data['n_projected']}, "
            f"avg_shift={data['avg_shift']:.3e}, max|phi(projected)|={data['max_abs_phi_proj']:.3e}, "
            f"avg_it={data['avg_it']:.2f}, avg_eval={data['avg_eval']:.2f}"
        )

    print("\n[midpoint projection comparison]")
    for s in all_stats:
        print(
            f"  {s['method']:>7s}: projected={s['n_projected']:2d}/{s['n_segments']:2d}, "
            f"avg_shift={s['avg_shift']:.3e}, max|phi|={s['max_abs_phi_proj']:.3e}, "
            f"avg_it={s['avg_it']:.2f}, avg_eval={s['avg_eval']:.2f}"
        )

    # Plot one representative method in detail.
    for s in all_stats:
        if s["method"] == "newton":
            _plot_projection_result(s)
            break


if __name__ == "__main__":
    main()
