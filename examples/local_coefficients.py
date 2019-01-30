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

    # ── 1. Circle — baseline vs trivial local system ──────────────────

    print('=' * 60)
    print('1. Circle — baseline vs trivial local system (weights = 1)')
    print('=' * 60)

    nP = 150
    P = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P = (P - P.min()) / (P.max() - P.min())

    pci = 0
    dim_max = 2
    r_max = 0.5

    # baseline (no local system)
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_base = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # trivial local system (all weights = 1)
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    n = P.shape[0]
    edges = np.array([[i, j] for i in range(n) for j in range(i+1, n)],
                     dtype=np.uint64)
    weights_trivial = np.ones(len(edges), dtype=np.int64)
    semel.set_local_system(pci, edges, weights_trivial)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_trivial = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nBaseline (no local system):')
    print('  betti:', q_base.get('betti', {}))
    for d in sorted(q_base['intervals'].keys()):
        print('  H%d: %d intervals' % (d, len(q_base['intervals'][d])))

    print('\nTrivial local system (all weights = 1):')
    print('  betti:', q_trivial.get('betti', {}))
    for d in sorted(q_trivial['intervals'].keys()):
        print('  H%d: %d intervals' % (d, len(q_trivial['intervals'][d])))

    print('\nCross-check (mass):')
    for d in sorted(set(q_base['1d']['mass'].keys()) | set(q_trivial['1d']['mass'].keys())):
        m_base = float(q_base['1d']['mass'].get(d, 0))
        m_trivial = float(q_trivial['1d']['mass'].get(d, 0))
        print('  H%d: baseline=%.6f  trivial=%.6f  diff=%.2e' % (d, m_base, m_trivial, abs(m_base - m_trivial)))

    assert str(q_base['betti']) == str(q_trivial['betti']), \
        'trivial local system must match baseline'

    # ── 2. Twisted local system — breaking the 1-cycle ────────────────

    print('\n' + '=' * 60)
    print('2. Twisted local system — weight 2 on one edge')
    print('=' * 60)
    print('\nA single non-unit weight on an edge through the circle')
    print('breaks the flat cocycle condition, changing the cohomology.')

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    weights_twisted = np.ones(len(edges), dtype=np.int64)
    weights_twisted[0] = 2
    semel.set_local_system(pci, edges, weights_twisted)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_twisted = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nTwisted local system:')
    print('  betti:', q_twisted.get('betti', {}))
    for d in sorted(q_twisted['intervals'].keys()):
        intervals = q_twisted['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d: %d intervals (%d free, %d torsion)' % (d, len(intervals), n_free, n_torsion))

    print('\nComparison:')
    print('  baseline betti: %s' % q_base.get('betti', {}))
    print('  twisted  betti: %s' % q_twisted.get('betti', {}))

    # ── 3. F_p path with local coefficients ───────────────────────────

    print('\n' + '=' * 60)
    print('3. F_p coefficients with local system')
    print('=' * 60)

    # baseline F_p
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_fp_base = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # twisted F_p
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.set_local_system(pci, edges, weights_twisted)
    semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
    q_fp_twisted = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nF_p baseline mass:')
    for d in sorted(q_fp_base['1d']['mass'].keys()):
        print('  H%d: %.6f' % (d, float(q_fp_base['1d']['mass'][d])))

    print('\nF_p twisted mass:')
    for d in sorted(q_fp_twisted['1d']['mass'].keys()):
        print('  H%d: %.6f' % (d, float(q_fp_twisted['1d']['mass'][d])))

    # ── 4. Varying weights ────────────────────────────────────────────

    print('\n' + '=' * 60)
    print('4. Varying the transport weight on a single edge')
    print('=' * 60)

    for w in [1, -1, 2, 3, 5]:
        semel.add_point_cloud(pci, MANIFOLD, P)
        semel.add_complex_ladic(pci, l=0)
        semel.generate_delaunay(pci)
        ww = np.ones(len(edges), dtype=np.int64)
        ww[0] = w
        semel.set_local_system(pci, edges, ww)
        semel.build_complex(pci, dim_max, r_max, filtration=semel.FILTRATION_RIPS)
        q = semel.cohomology(pci, measures_only=False)
        semel.unload_complex(pci)
        semel.unload_point_cloud(pci)
        print('  w=%2d  betti=%s  torsion=%s' % (w, q.get('betti', {}), q.get('torsion', {})))

    print('\nall tests passed')
