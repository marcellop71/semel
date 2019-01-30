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

    # ── 1. Circle — no torsion expected ──────────────────────────────

    print('=' * 60)
    print('1. Circle (annulus) — expected: beta_0=1, beta_1=1, no torsion')
    print('=' * 60)

    nP = 100
    P = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P = (P - P.min()) / (P.max() - P.min())

    pci = 0
    dim_max = 2
    r_max = 0.5

    # F_p path (baseline)
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max)
    q_fp = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # integral path
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max)
    q_z = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nF_p intervals:')
    for d in sorted(q_fp['intervals'].keys()):
        print('  H%d: %d intervals' % (d, len(q_fp['intervals'][d])))

    print('\nIntegral (Z) results:')
    print('  betti:', q_z.get('betti', {}))
    print('  torsion:', q_z.get('torsion', {}))
    for d in sorted(q_z['intervals'].keys()):
        intervals = q_z['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d: %d intervals (%d free, %d torsion)' % (d, len(intervals), n_free, n_torsion))

    # cross-check: persistence masses should agree
    print('\nCross-check (mass):')
    for d in sorted(set(q_fp['1d']['mass'].keys()) | set(q_z['1d']['mass'].keys())):
        m_fp = float(q_fp['1d']['mass'].get(d, 0))
        m_z = float(q_z['1d']['mass'].get(d, 0))
        print('  H%d: F_p=%.6f  Z=%.6f  diff=%.2e' % (d, m_fp, m_z, abs(m_fp - m_z)))

    # ── 2. Random point cloud — mostly noise ─────────────────────────

    print('\n' + '=' * 60)
    print('2. Random point cloud — expected: beta_0>0, few persistent features')
    print('=' * 60)

    nP = 150
    P = np.random.rand(nP, 2)

    pci = 0
    dim_max = 2
    r_max = 0.3

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max)
    q_z = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nIntegral (Z) results:')
    print('  betti:', q_z.get('betti', {}))
    print('  torsion:', q_z.get('torsion', {}))
    for d in sorted(q_z['intervals'].keys()):
        intervals = q_z['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d: %d intervals (%d free, %d torsion)' % (d, len(intervals), n_free, n_torsion))

    # ── 3. l-adic filtering ──────────────────────────────────────────

    print('\n' + '=' * 60)
    print('3. l-adic filtering (l=2) on same random point cloud')
    print('=' * 60)

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=2)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max)
    q_l2 = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\n2-adic results:')
    print('  l:', q_l2.get('l', 'N/A'))
    print('  betti:', q_l2.get('betti', {}))
    print('  torsion:', q_l2.get('torsion', {}))

    # ── 4. Two concentric circles ────────────────────────────────────

    print('\n' + '=' * 60)
    print('4. Two concentric circles — expected: beta_0=1, beta_1=2, no torsion')
    print('=' * 60)

    n_inner = 40
    n_outer = 60
    P_inner = pc.circle(n_inner, cx=0.0, cy=0.0, rmin=0.35, rmax=0.45, method=1)
    P_outer = pc.circle(n_outer, cx=0.0, cy=0.0, rmin=0.85, rmax=0.95, method=1)
    P = np.vstack([P_inner, P_outer])
    P = (P - P.min()) / (P.max() - P.min())

    pci = 0
    dim_max = 2
    r_max = 0.5

    # F_p baseline
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max)
    q_fp = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # integral
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, dim_max, r_max)
    q_z = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\nF_p intervals:')
    for d in sorted(q_fp['intervals'].keys()):
        print('  H%d: %d intervals' % (d, len(q_fp['intervals'][d])))

    print('\nIntegral (Z) results:')
    print('  betti:', q_z.get('betti', {}))
    print('  torsion:', q_z.get('torsion', {}))
    for d in sorted(q_z['intervals'].keys()):
        intervals = q_z['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d: %d intervals (%d free, %d torsion)' % (d, len(intervals), n_free, n_torsion))

    print('\nCross-check (mass):')
    for d in sorted(set(q_fp['1d']['mass'].keys()) | set(q_z['1d']['mass'].keys())):
        m_fp = float(q_fp['1d']['mass'].get(d, 0))
        m_z = float(q_z['1d']['mass'].get(d, 0))
        print('  H%d: F_p=%.6f  Z=%.6f  diff=%.2e' % (d, m_fp, m_z, abs(m_fp - m_z)))
