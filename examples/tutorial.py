import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import matplotlib.pyplot as plt

from semel import Semel, PointCloud, Graph, TimeSeries

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)

    semel = Semel()
    pc = PointCloud()
    graph = Graph()

    # 1. Circle (annulus) — expected: one 1-cycle
    np.random.seed(42)

    nP = 80
    P_circle = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P_circle = (P_circle - P_circle.min()) / (P_circle.max() - P_circle.min())

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.scatter(P_circle[:, 0], P_circle[:, 1], s=12)
    ax.set_aspect('equal')
    ax.set_title('Circle point cloud (%d points)' % nP)
    ax.grid(True)
    plt.savefig('output/tutorial_circle.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    pci = 0
    manifold = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }
    semel.add_point_cloud(pci, manifold, P_circle)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)

    dim_max = 2
    r_max = 0.6
    semel.build_complex(pci, dim_max, r_max)

    q_circle = semel.cohomology(pci, measures_only=False)
    print('H^0 mass: %.4f, entropy: %.4f' % (q_circle['1d']['mass'][0], q_circle['1d']['entropy'][0]))
    print('H^1 mass: %.4f, entropy: %.4f' % (q_circle['1d']['mass'][1], q_circle['1d']['entropy'][1]))

    graph.bundle_0(dim_max, r_max, q_circle, base_folder='output', name='tutorial_circle')

    q_complex = semel.export_complex(pci)
    graph.complex(P_circle, q_complex, labels=False, xlim=[-0.1, 1.1], base_folder='output', name='tutorial_complex')

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # 2. Random point cloud — expected: no persistent features
    nP = 200
    d = 2
    P_random = pc.random(nP, d)
    P_random = (P_random - P_random.min()) / (P_random.max() - P_random.min())

    pci = 0
    manifold = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }
    semel.add_point_cloud(pci, manifold, P_random)
    semel.add_complex(pci)

    semel.generate_delaunay(pci)

    dim_max = 2
    r_max = 0.3
    semel.build_complex(pci, dim_max, r_max)

    q_random = semel.cohomology(pci, measures_only=False)
    print('H^0 mass: %.4f, entropy: %.4f' % (q_random['1d']['mass'][0], q_random['1d']['entropy'][0]))
    if 1 in q_random['1d']['mass']:
        print('H^1 mass: %.4f, entropy: %.4f' % (q_random['1d']['mass'][1], q_random['1d']['entropy'][1]))
    else:
        print('H^1 mass: 0 (no 1-cycles)')

    graph.bundle_0(dim_max, r_max, q_random, base_folder='output', name='tutorial_random')

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # 3. Time series embedding
    t = np.linspace(0, 4 * np.pi, 300)
    signal = np.sin(t) + 0.05 * np.random.randn(len(t))

    ts = TimeSeries(points=signal.tolist())
    P_ts = ts.embedding(dim=2, tau=15, method=0)
    P_ts = (P_ts - P_ts.min()) / (P_ts.max() - P_ts.min())

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].plot(t, signal)
    axes[0].set_title('Signal')
    axes[0].set_xlabel('t')
    axes[0].grid(True)

    axes[1].scatter(P_ts[:, 0], P_ts[:, 1], s=8)
    axes[1].set_title('Time-delay embedding (d=2, tau=15)')
    axes[1].set_aspect('equal')
    axes[1].grid(True)
    plt.tight_layout()
    plt.savefig('output/tutorial_timeseries_embedding.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    pci = 0
    manifold = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }
    semel.add_point_cloud(pci, manifold, P_ts)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)

    dim_max = 2
    r_max = 0.5
    semel.build_complex(pci, dim_max, r_max)

    q_ts = semel.cohomology(pci, measures_only=False)
    print('H^0 mass: %.4f' % q_ts['1d']['mass'][0])
    if 1 in q_ts['1d']['mass']:
        print('H^1 mass: %.4f  (1-cycle detected)' % q_ts['1d']['mass'][1])

    graph.bundle_0(dim_max, r_max, q_ts, base_folder='output', name='tutorial_timeseries')

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # 4. Coboundary operator sparsity
    pci = 0
    manifold = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }
    semel.add_point_cloud(pci, manifold, P_circle)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.6)

    for d in range(2):
        q_cob = semel.export_coboundary(pci, d)
        nnz = len(q_cob['matrix_sparse'])
        total = q_cob['m'] * q_cob['n']
        print('d^%d: C^%d (R^%d) -> C^%d (R^%d), nnz=%d / %d (%.2f%%)' % (
            d, d, q_cob['n'], d+1, q_cob['m'], nnz, total, 100 * nnz / total))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)
