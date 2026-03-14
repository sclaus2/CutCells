import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import cutcells


def demo_triangle(order: int):
    tpl = cutcells.iso_p1_template(cutcells.CellType.triangle, order)
    parent_cell_coords = np.asarray(
        cutcells.iso_p1_ref_coords(cutcells.CellType.triangle, order), dtype=np.float64
    )

    lm = cutcells.init_local_mesh_from_template(
        tpl,
        parent_cell_coords,
        cutcells.CellType.triangle,
        parent_cell_id=0,
        n_level_sets=1,
    )

    print("LocalMesh")
    print("  gdim:", lm.gdim)
    print("  tdim:", lm.tdim)
    print("  n_vertices:", lm.n_vertices())
    print("  n_edges:", lm.n_edges())
    print("  n_cells:", lm.n_cells())
    print("  cell_offsets:", np.asarray(lm.cell_offsets))
    print("  cell_vertices:", np.asarray(lm.cell_vertices))

    x = np.asarray(lm.vertex_x, dtype=np.float64).reshape(-1, lm.gdim)
    cells = np.asarray(lm.cell_vertices, dtype=np.int32).reshape(-1, 3)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.triplot(x[:, 0], x[:, 1], cells, color="black", lw=1.2)
    ax.scatter(x[:, 0], x[:, 1], s=40, color="tab:blue", zorder=3)

    for i, (xi, yi) in enumerate(x):
        ax.text(xi + 0.015, yi + 0.015, f"v{i}", fontsize=9)

    ax.set_title(f"Local Mesh from Iso-Refine Template (Triangle P{order})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    plt.show()

def demo_tetrahedron(order: int):
    tpl = cutcells.iso_p1_template(cutcells.CellType.tetrahedron, order)
    parent_cell_coords = np.asarray(
        cutcells.iso_p1_ref_coords(cutcells.CellType.tetrahedron, order), dtype=np.float64
    )

    lm = cutcells.init_local_mesh_from_template(
        tpl,
        parent_cell_coords,
        cutcells.CellType.tetrahedron,
        parent_cell_id=1,
        n_level_sets=1,
    )

    print("\nLocalMesh (tetrahedron)")
    print("  gdim:", lm.gdim)
    print("  tdim:", lm.tdim)
    print("  n_vertices:", lm.n_vertices())
    print("  n_edges:", lm.n_edges())
    print("  n_cells:", lm.n_cells())

    x = np.asarray(lm.vertex_x, dtype=np.float64).reshape(-1, lm.gdim)
    tet_cells = np.asarray(lm.cell_vertices, dtype=np.int32).reshape(-1, 4)

    # Draw unique tetra edges as 3D line segments.
    edge_pairs = set()
    for a, b, c, d in tet_cells:
        verts = [a, b, c, d]
        for i in range(4):
            for j in range(i + 1, 4):
                e0, e1 = sorted((verts[i], verts[j]))
                edge_pairs.add((e0, e1))
    segments = [[x[i], x[j]] for i, j in edge_pairs]

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Line3DCollection(segments, colors="black", linewidths=0.9))
    ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=35, color="tab:red", depthshade=False)

    for i, (xi, yi, zi) in enumerate(x):
        ax.text(xi + 0.02, yi + 0.02, zi + 0.02, f"v{i}", fontsize=8)

    ax.set_title(f"Local Mesh from Iso-Refine Template (Tetrahedron P{order})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect((1, 1, 1))
    fig.tight_layout()
    plt.show()

def main():
    for order in (2, 3):
        demo_triangle(order)
    for order in (2, 3):
        demo_tetrahedron(order)


if __name__ == "__main__":
    main()
