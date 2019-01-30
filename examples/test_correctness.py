"""Correctness regression test for semel.

Exercises all code paths (F_p, Z, local coefficients, sheaf, Laplacian,
coboundary export) on multiple topologies and prints deterministic numerical
outputs that can be diffed across versions.
"""
import os
import json
import numpy as np

from semel import Semel, PointCloud

def fmt(x, decimals=6):
    if isinstance(x, float) or isinstance(x, np.floating):
        return ('%.*f' % (decimals, x))
    return str(x)

if __name__ == '__main__':
    semel = Semel()
    pc = PointCloud()

    np.random.seed(42)

    MANIFOLD = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }

    # ── Test 1: Circle, F_p cohomology (default) ────────────────────────

    print('=' * 60)
    print('TEST 1: Circle F_p cohomology')
    print('=' * 60)

    nP = 80
    P = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P = (P - P.min()) / (P.max() - P.min())

    pci = 0
    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.6)
    q = semel.cohomology(pci, measures_only=False)

    for d in sorted(q['1d']['mass'].keys()):
        print('  H%d mass=%s entropy=%s' % (d, fmt(q['1d']['mass'][d]), fmt(q['1d']['entropy'][d])))
    for d in sorted(q['intervals'].keys()):
        intervals = q['intervals'][d]
        print('  H%d intervals=%d' % (d, len(intervals)))
        # print first 5 intervals sorted by persistence
        if len(intervals) > 0:
            arr = np.array(intervals)
            pers = arr[:, 1] - arr[:, 0]
            idx = np.argsort(-pers)[:5]
            for i in idx:
                print('    [%s, %s] pers=%s' % (fmt(arr[i, 0]), fmt(arr[i, 1]), fmt(pers[i])))
    for d in sorted(q['cocycles'].keys()):
        print('  H%d cocycles=%d' % (d, len(q['cocycles'][d])))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # ── Test 2: Circle, coboundary operator export ──────────────────────

    print('\n' + '=' * 60)
    print('TEST 2: Coboundary operator export')
    print('=' * 60)

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.6)

    for d in range(2):
        cob = semel.export_coboundary(pci, d)
        nnz = len(cob['matrix_sparse'])
        print('  d^%d: %dx%d, nnz=%d' % (d, cob['m'], cob['n'], nnz))
        # checksum of sparse entries
        if nnz > 0:
            arr = np.array(cob['matrix_sparse'])
            checksum = float(np.sum(arr[:, 0]) + np.sum(arr[:, 1]) + np.sum(arr[:, 2]))
            print('    checksum=%s' % fmt(checksum))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # ── Test 3: Circle, Hodge Laplacian ─────────────────────────────────

    print('\n' + '=' * 60)
    print('TEST 3: Hodge Laplacian eigenvalues')
    print('=' * 60)

    semel.add_point_cloud(pci, MANIFOLD, P)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.6)

    eig = semel.HodgeLaplace(pci, 2, nev=0)
    for d, vals in enumerate(eig['eigenvalues']):
        if vals is not None and len(vals) > 0:
            n_zero = int(np.sum(np.abs(vals) < 1e-8))
            print('  L_%d: %d eigenvalues, %d near-zero' % (d, len(vals), n_zero))
            print('    smallest 5: [%s]' % ', '.join(fmt(v) for v in vals[:5]))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # ── Test 4: Random cloud, F_p ───────────────────────────────────────

    print('\n' + '=' * 60)
    print('TEST 4: Random point cloud F_p')
    print('=' * 60)

    nP = 150
    P_rand = pc.random(nP, 2)
    P_rand = (P_rand - P_rand.min()) / (P_rand.max() - P_rand.min())

    semel.add_point_cloud(pci, MANIFOLD, P_rand)
    semel.add_complex(pci)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.3)
    q = semel.cohomology(pci, measures_only=False)

    for d in sorted(q['1d']['mass'].keys()):
        print('  H%d mass=%s' % (d, fmt(q['1d']['mass'][d])))
    for d in sorted(q['intervals'].keys()):
        print('  H%d intervals=%d' % (d, len(q['intervals'][d])))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # ── Test 5: Circle, Z (integral) cohomology ─────────────────────────

    print('\n' + '=' * 60)
    print('TEST 5: Circle integral (Z) cohomology')
    print('=' * 60)

    nP = 100
    P_circ = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P_circ = (P_circ - P_circ.min()) / (P_circ.max() - P_circ.min())

    semel.add_point_cloud(pci, MANIFOLD, P_circ)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q = semel.cohomology(pci, measures_only=False)

    print('  betti: %s' % q.get('betti', {}))
    print('  torsion: %s' % q.get('torsion', {}))
    for d in sorted(q['1d']['mass'].keys()):
        print('  H%d mass=%s' % (d, fmt(q['1d']['mass'][d])))
    for d in sorted(q['intervals'].keys()):
        intervals = q['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d intervals=%d (free=%d, torsion=%d)' % (d, len(intervals), n_free, n_torsion))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # ── Test 6: Local coefficients (trivial = baseline) ─────────────────

    print('\n' + '=' * 60)
    print('TEST 6: Local coefficients (trivial must match baseline)')
    print('=' * 60)

    semel.add_point_cloud(pci, MANIFOLD, P_circ)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q_base = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    semel.add_point_cloud(pci, MANIFOLD, P_circ)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    n = P_circ.shape[0]
    edges = np.array([[i, j] for i in range(n) for j in range(i+1, n)], dtype=np.uint64)
    weights_trivial = np.ones(len(edges), dtype=np.int64)
    semel.set_local_system(pci, edges, weights_trivial)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q_trivial = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    match = True
    for d in sorted(set(q_base['1d']['mass'].keys()) | set(q_trivial['1d']['mass'].keys())):
        m_base = float(q_base['1d']['mass'].get(d, 0))
        m_trivial = float(q_trivial['1d']['mass'].get(d, 0))
        diff = abs(m_base - m_trivial)
        ok = 'OK' if diff < 1e-10 else 'MISMATCH'
        if diff >= 1e-10:
            match = False
        print('  H%d: base=%s trivial=%s diff=%s %s' % (d, fmt(m_base), fmt(m_trivial), fmt(diff, 2), ok))
    print('  betti match: %s' % (str(q_base.get('betti', {})) == str(q_trivial.get('betti', {}))))
    print('  PASS' if match else '  FAIL')

    # ── Test 7: Twisted local coefficients ──────────────────────────────

    print('\n' + '=' * 60)
    print('TEST 7: Twisted local coefficients')
    print('=' * 60)

    semel.add_point_cloud(pci, MANIFOLD, P_circ)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    weights_twisted = np.ones(len(edges), dtype=np.int64)
    weights_twisted[0] = 2
    semel.set_local_system(pci, edges, weights_twisted)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q_tw = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('  betti: %s' % q_tw.get('betti', {}))
    for d in sorted(q_tw['1d']['mass'].keys()):
        print('  H%d mass=%s' % (d, fmt(q_tw['1d']['mass'][d])))
    for d in sorted(q_tw['intervals'].keys()):
        intervals = q_tw['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d intervals=%d (free=%d, torsion=%d)' % (d, len(intervals), n_free, n_torsion))

    # ── Test 8: Sheaf cohomology (rank-1 = local system) ────────────────

    print('\n' + '=' * 60)
    print('TEST 8: Sheaf cohomology (rank-1 matches local_system)')
    print('=' * 60)

    nP = 80
    P_sh = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P_sh = (P_sh - P_sh.min()) / (P_sh.max() - P_sh.min())

    # scalar local_system with weight 2 on edge 0
    semel.add_point_cloud(pci, MANIFOLD, P_sh)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    n = P_sh.shape[0]
    edges_sh = np.array([[i, j] for i in range(n) for j in range(i+1, n)], dtype=np.uint64)
    w_scalar = np.ones(len(edges_sh), dtype=np.int64)
    w_scalar[0] = 2
    semel.set_local_system(pci, edges_sh, w_scalar)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q_scalar = semel.cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # rank-1 sheaf with same weight
    semel.add_point_cloud(pci, MANIFOLD, P_sh)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    matrices_r1 = np.ones((len(edges_sh), 1, 1), dtype=np.int64)
    matrices_r1[0, 0, 0] = 2
    semel.set_sheaf(pci, 1, edges_sh, matrices_r1)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q_sheaf = semel.sheaf_cohomology(pci, measures_only=False)
    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('  scalar betti: %s' % q_scalar.get('betti', {}))
    print('  sheaf  betti: %s' % q_sheaf.get('betti', {}))
    match = True
    for d in sorted(set(q_scalar['1d']['mass'].keys()) | set(q_sheaf['1d']['mass'].keys())):
        m_s = float(q_scalar['1d']['mass'].get(d, 0))
        m_sh = float(q_sheaf['1d']['mass'].get(d, 0))
        diff = abs(m_s - m_sh)
        ok = 'OK' if diff < 1e-10 else 'MISMATCH'
        if diff >= 1e-10:
            match = False
        print('  H%d: scalar=%s sheaf=%s diff=%s %s' % (d, fmt(m_s), fmt(m_sh), fmt(diff, 2), ok))
    print('  PASS' if match else '  FAIL')

    # ── Test 9: Bottom-up Rips (no Delaunay) ────────────────────────────

    print('\n' + '=' * 60)
    print('TEST 9: Bottom-up Rips (no Delaunay)')
    print('=' * 60)

    nP = 60
    P_rips = pc.circle(nP, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
    P_rips = (P_rips - P_rips.min()) / (P_rips.max() - P_rips.min())

    semel.add_point_cloud(pci, MANIFOLD, P_rips)
    semel.add_complex(pci)
    # no generate_delaunay: bottom-up Rips construction
    semel.build_complex(pci, 2, 0.25, filtration=semel.FILTRATION_RIPS)
    q_rips = semel.cohomology(pci, measures_only=False)

    for d in sorted(q_rips['1d']['mass'].keys()):
        print('  H%d mass=%s' % (d, fmt(q_rips['1d']['mass'][d])))
    for d in sorted(q_rips['intervals'].keys()):
        print('  H%d intervals=%d' % (d, len(q_rips['intervals'][d])))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    # ── Test 10: Two concentric circles, Z ──────────────────────────────

    print('\n' + '=' * 60)
    print('TEST 10: Two concentric circles, integral Z')
    print('=' * 60)

    n_inner = 40
    n_outer = 60
    P_inner = pc.circle(n_inner, cx=0.0, cy=0.0, rmin=0.35, rmax=0.45, method=1)
    P_outer = pc.circle(n_outer, cx=0.0, cy=0.0, rmin=0.85, rmax=0.95, method=1)
    P_2c = np.vstack([P_inner, P_outer])
    P_2c = (P_2c - P_2c.min()) / (P_2c.max() - P_2c.min())

    semel.add_point_cloud(pci, MANIFOLD, P_2c)
    semel.add_complex_ladic(pci, l=0)
    semel.generate_delaunay(pci)
    semel.build_complex(pci, 2, 0.5, filtration=semel.FILTRATION_RIPS)
    q_2c = semel.cohomology(pci, measures_only=False)

    print('  betti: %s' % q_2c.get('betti', {}))
    for d in sorted(q_2c['1d']['mass'].keys()):
        print('  H%d mass=%s' % (d, fmt(q_2c['1d']['mass'][d])))
    for d in sorted(q_2c['intervals'].keys()):
        intervals = q_2c['intervals'][d]
        n_free = sum(1 for iv in intervals if iv[2] == 0)
        n_torsion = sum(1 for iv in intervals if iv[2] > 0)
        print('  H%d intervals=%d (free=%d, torsion=%d)' % (d, len(intervals), n_free, n_torsion))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)

    print('\n' + '=' * 60)
    print('ALL TESTS COMPLETE')
    print('=' * 60)
