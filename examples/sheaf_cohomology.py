import os
import numpy as np

from semel import Semel, PointCloud

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)

    semel = Semel()
    pc = PointCloud()

    np.random.seed(42)

    MANIFOLD = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }

    nP = 150
    P = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P = (P - P.min()) / (P.max() - P.min())

    pci = 0
    dim_max = 2
    r_max = 0.5

    # ── 1. Rank-1 sheaf matches scalar local_system ──────────────────

    print('=' * 60)
    print('1. Rank-1 sheaf matches scalar local_system')
    print('=' * 60)

    n = P.shape[0]
    edges = np.array([[i, j] for i in range(n) for j in range(i+1, n)],
                     dtype=np.uint64)

    # scalar local_system with weight 2 on first edge
    weights_scalar = np.ones(len(edges), dtype=np.int64)
    weights_scalar[0] = 2

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.set_local_system(pci, edges, weights_scalar)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_scalar = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # rank-1 sheaf with same weight
    matrices_r1 = weights_scalar.reshape(-1, 1, 1)

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.set_sheaf(pci, 1, edges, matrices_r1)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_sheaf_r1 = semel.sheaf_cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nScalar local_system betti:', q_scalar.get('betti', {}))
    print('Rank-1 sheaf betti:       ', q_sheaf_r1.get('betti', {}))

    for d in sorted(set(q_scalar['1d']['mass'].keys()) | set(q_sheaf_r1['1d']['mass'].keys())):
        m_s = float(q_scalar['1d']['mass'].get(d, 0))
        m_r1 = float(q_sheaf_r1['1d']['mass'].get(d, 0))
        print('  H%d: scalar=%.6f  sheaf_r1=%.6f  diff=%.2e' % (d, m_s, m_r1, abs(m_s - m_r1)))

    assert str(q_scalar.get('betti', {})) == str(q_sheaf_r1.get('betti', {})), \
        'rank-1 sheaf must match scalar local_system'

    # ── 2. Identity rank-2 sheaf doubles Betti numbers ───────────────

    print('\n' + '=' * 60)
    print('2. Identity rank-2 sheaf doubles Betti numbers')
    print('=' * 60)

    # baseline (no local system)
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_base = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # identity rank-2 sheaf
    matrices_id = np.array([np.eye(2, dtype=np.int64)] * len(edges), dtype=np.int64)

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.set_sheaf(pci, 2, edges, matrices_id)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_sheaf_id = semel.sheaf_cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nBaseline betti:          ', q_base.get('betti', {}))
    print('Identity rank-2 betti:   ', q_sheaf_id.get('betti', {}))

    betti_base = q_base.get('betti', {})
    betti_sheaf = q_sheaf_id.get('betti', {})
    for d in sorted(betti_base.keys()):
        b = int(betti_base[d])
        s = int(betti_sheaf.get(d, 0))
        print('  H%d: baseline=%d  sheaf=%d  ratio=%.1f' % (d, b, s, s / b if b > 0 else 0))
        assert s == 2 * b, 'identity rank-2 sheaf should double Betti number in H%d' % d

    # ── 3. Non-trivial rank-2 transport ──────────────────────────────

    print('\n' + '=' * 60)
    print('3. Non-trivial rank-2 transport (rotation on one edge)')
    print('=' * 60)

    matrices_twisted = np.array([np.eye(2, dtype=np.int64)] * len(edges), dtype=np.int64)
    # 90-degree integer rotation on first edge
    matrices_twisted[0] = np.array([[0, -1], [1, 0]], dtype=np.int64)

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.set_sheaf(pci, 2, edges, matrices_twisted)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_sheaf_tw = semel.sheaf_cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nTwisted rank-2 betti:    ', q_sheaf_tw.get('betti', {}))
    print('(compare: identity r2:   ', betti_sheaf, ')')

    for d in sorted(q_sheaf_tw['1d']['mass'].keys()):
        print('  H%d mass: %.6f' % (d, float(q_sheaf_tw['1d']['mass'][d])))

    # ── 4. Sheaf Laplacian eigenvalues ───────────────────────────────

    print('\n' + '=' * 60)
    print('4. Sheaf Laplacian eigenvalues')
    print('=' * 60)

    # use F_p for Laplacian (identity rank-2)
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.set_sheaf(pci, 2, edges, matrices_id)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_lap = semel.sheaf_Laplacian(pci, d=dim_max, nev=10)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    for i, eig in enumerate(q_lap['eigenvalues']):
        if len(eig) > 0:
            n_zero = np.sum(np.abs(eig) < 1e-8)
            print('  L_%d: %d eigenvalues, %d near-zero' % (i+1, len(eig), n_zero))
            print('    smallest:', eig[:min(5, len(eig))])

    print('\nall tests passed')
