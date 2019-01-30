import matplotlib
matplotlib.use('Agg')

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from semel import Semel, PointCloud, Graph
from semel import CSM
from semel import hodge

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)

    semel = Semel()

    K = 1.8
    csm = CSM(K)
    csm.set_start_random()
    csm.iterate(num_iter=1000, track=True)
    P = csm.get_point_cloud()

    Graph().point_cloud(P, S=None, labels=False, xlim=[0, 2 * np.pi],
                        base_folder='output', name='chirikov_pc')

    pci = 0
    manifold = {
        'type': semel.MANIFOLD_SURFACE_TORUS,
        'distance': semel.DISTANCE_SURFACE_TORUS,
        'border': [2 * np.pi, 2 * np.pi]
    }
    semel.add_point_cloud(pci, manifold, P)
    semel.add_complex(pci)

    # Periodic Delaunay triangulation via shifted copies
    from pyhull.delaunay import DelaunayTri

    border = manifold['border']

    def shift(P, i, j, border):
        Q = np.copy(P)
        Q[:, 0] = P[:, 0] + i * border[0]
        Q[:, 1] = P[:, 1] + j * border[1]
        return Q

    P_ = P
    for i in [-1, 0, +1]:
        for j in [-1, 0, +1]:
            if (i != 0) or (j != 0):
                P_ = np.concatenate([P_, shift(P, i, j, border)])

    nP = P.shape[0]
    S = [[x % nP, y % nP, z % nP] for [x, y, z] in DelaunayTri(P_).vertices]
    S = [list(x) for x in set(tuple(x) for x in S)]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.scatter(P_[:, 0], P_[:, 1], color='blue')
    ax.scatter(P[:, 0], P[:, 1], color='red')
    plt.savefig('output/chirikov_delaunay_copies.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    typeS = 0
    semel.add_simplices(pci, typeS, S)

    dim_max = 2
    r_max = 0.35
    semel.build_complex(pci, dim_max, r_max)

    q_complex = semel.export_complex(pci)
    q_cohomology = semel.cohomology(pci, measures_only=False)
    q_cocycles = semel.export_cocyles(pci)

    # 3D torus visualisation with cocycles
    def points_on_torus(R, r, x_, y_):
        px = (R + r * np.cos(x_)) * np.cos(y_)
        py = (R + r * np.cos(x_)) * np.sin(y_)
        pz = r * np.sin(x_)
        return [px, py, pz]

    def line(ax, P0, P1, color='b', linewidth=0.6, alpha=0.4):
        [x0, y0] = P0; [x1, y1] = P1
        steps = 100
        t = np.linspace(0.0, 1.0, steps)
        tmp_x_1 = np.absolute(x0 - x1); tmp_x_2 = border[0] - tmp_x_1
        tmp_y_1 = np.absolute(y0 - y1); tmp_y_2 = border[1] - tmp_y_1
        if tmp_x_1 < tmp_x_2:
            x_ = ((1.0 - t) * x0) + (t * x1)
        else:
            x_ = (t * x0) + ((1.0 - t) * x1)
        if tmp_y_1 < tmp_y_2:
            y_ = ((1.0 - t) * y0) + (t * y1)
        else:
            y_ = (t * y0) + ((1.0 - t) * y1)
        [px, py, pz] = points_on_torus(R, r, x_, y_)
        ax.plot(px, py, pz, color=color, linewidth=linewidth, alpha=alpha)

    R = 1.0; r = 0.40

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-1.5, 1.5)
    ax.set_ylim3d(-1.5, 1.5)
    ax.set_zlim3d(-1.5, 1.5)

    [px, py, pz] = points_on_torus(R, r, P[:, 0], P[:, 1])
    ax.scatter(px, py, pz, s=10, c='b')

    for s in q_complex['complex']:
        s_ = s['simplex']
        if s['dimension'] == 1:
            line(ax, P[s_[0], :], P[s_[1], :])

    colors = ['blue', 'red', 'green', 'yellow', 'cyan', 'black']
    for i, z in enumerate(q_cocycles):
        if z['dimension'] == 1:
            q = z['expression']
            for tmp in q:
                s_ = q_complex['complex'][tmp[1]]['simplex']
                line(ax, P[s_[0], :], P[s_[1], :], color=colors[i % len(colors)], linewidth=2.0, alpha=1.0)

    ax.view_init(30, 60)
    plt.savefig('output/chirikov_torus_cocycles.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    Graph().bundle_0(dim_max, r_max, q_cohomology, base_folder='output', name='chirikov_bundle')

    print(q_cohomology['1d'])

    for dim in [0, 1]:
        q = semel.export_coboundary(pci, dim)
        total = q['m'] * q['n']
        nnz = len(q['matrix_sparse'])
        pct = 100 * nnz / total if total > 0 else 0.0
        print("d^%d: C^%d (R^%d) -> C^%d (R^%d), nnz=%d / %d = %.2f%%" % (
            dim, dim, q['n'], dim + 1, q['m'], nnz, total, pct))

    # Hodge decomposition
    q_complex = semel.export_complex(pci)
    q_cocycles = semel.export_cocyles(pci)
    q_coboundary = semel.export_coboundary(pci, 0)

    [q_cc, y_, b_] = hodge(q_complex, q_cocycles, q_coboundary)

    Graph().point_cloud_cc(P, q_cc[0], labels=False, xlim=[0, 2 * np.pi],
                           base_folder='output', name='chirikov_cc')

    # 3D torus visualisation with Hodge
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-1.2, 1.2)
    ax.set_ylim3d(-1.2, 1.2)
    ax.set_zlim3d(-1.2, 1.2)

    [px, py, pz] = points_on_torus(R, r, P[:, 0], P[:, 1])
    ax.scatter(px, py, pz, s=10, c='gray')

    for s in q_complex['complex']:
        s_ = s['simplex']
        if s['dimension'] == 1:
            line(ax, P[s_[0], :], P[s_[1], :], color='gray', linewidth=0.6, alpha=0.4)

    colors = ['cyan', 'black', 'green', 'yellow', 'cyan', 'black']
    for i, z in enumerate(y_):
        idx = np.argsort(np.absolute(z), axis=0)[-40:]
        for tmp in idx:
            for qtmp in q_complex['complex']:
                if qtmp['dimension'] == 1:
                    if qtmp['index'] == tmp[0]:
                        s_ = qtmp['simplex']
                        line(ax, P[s_[0], :], P[s_[1], :], color=colors[i % len(colors)], linewidth=2.0, alpha=1.0)

    ax.view_init(30, 0)
    plt.savefig('output/chirikov_torus_hodge.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # Bifurcation diagram from pre-computed data
    if os.path.exists('data/csm_00.dat'):
        with open('data/csm_00.dat', 'rb') as f:
            [K_, mass, entropy] = pickle.load(f)

        fig = plt.figure(figsize=(20, 20))

        ax0 = fig.add_subplot(221)
        ax0.scatter(K_, mass[:, 0], s=2)
        ax0.set_title('H^0 mass')
        ax0.grid()

        ax1 = fig.add_subplot(222)
        ax1.scatter(K_, entropy[:, 0], s=2)
        ax1.set_ylim([-1, 200])
        ax1.set_title('H^0 entropy')
        ax1.grid()

        ax2 = fig.add_subplot(223)
        ax2.scatter(K_, mass[:, 1], s=2)
        ax2.set_title('H^1 mass')
        ax2.grid()

        ax3 = fig.add_subplot(224)
        ax3.scatter(K_, entropy[:, 1], s=2)
        ax3.set_title('H^1 entropy')
        ax3.grid()

        plt.savefig('output/chirikov_bifurcation.jpg', dpi=150, bbox_inches='tight')
        plt.close()
    else:
        print("data/csm_00.dat not found — skipping bifurcation plot")
