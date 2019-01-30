import matplotlib
matplotlib.use('Agg')

import os
import time
import numpy as np
import matplotlib.pyplot as plt

from semel import Semel, PointCloud
from ripser import ripser
from persim import bottleneck

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)

    semel_ctx = Semel()
    pc = PointCloud()

    # ── helpers ──────────────────────────────────────────────────────────

    MANIFOLD = {
        'type': semel_ctx.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel_ctx.DISTANCE_L2,
        'border': []
    }

    def run_semel(P, dim_max=2, r_max=0.5):
        """Delaunay + Rips pipeline.
        r_max is in diameter units (same as ripser's thresh).
        Semel Rips stores 0.5*max_pairwise_dist (radius), so we pass r_max/2.
        Delaunay provides the sparse skeleton; Rips filtration values for comparison.
        """
        pci = 0
        t0 = time.perf_counter()
        semel_ctx.add_point_cloud(pci, MANIFOLD, P)
        semel_ctx.add_complex(pci)
        semel_ctx.generate_delaunay(pci)
        semel_ctx.build_complex(pci, dim_max, r_max / 2,
                                filtration=semel_ctx.FILTRATION_RIPS)
        q = semel_ctx.cohomology(pci, measures_only=False)
        elapsed = time.perf_counter() - t0
        semel_ctx.unload_complex(pci)
        semel_ctx.unload_point_cloud(pci)
        return q, elapsed

    def run_ripser(P, maxdim=1, thresh=0.5):
        t0 = time.perf_counter()
        result = ripser(P, maxdim=maxdim, thresh=thresh)
        elapsed = time.perf_counter() - t0
        return result, elapsed

    def semel_to_diagrams(q, r_max):
        """Convert semel finite intervals to {dim: ndarray([birth, death])}.
        Semel Rips stores radius; ripser stores diameter → scale ×2.
        Excludes surviving cocycles (ripser also filters them out).
        """
        dgms = {}
        for dim, intervals in q['intervals'].items():
            if len(intervals) > 0:
                dgms[dim] = np.array(intervals) * 2  # radius → diameter
            else:
                dgms[dim] = np.empty((0, 2))
        return dgms

    def ripser_to_diagrams(result):
        """Convert ripser output to {dim: ndarray([birth, death])}, finite only."""
        dgms = {}
        for dim, dgm in enumerate(result['dgms']):
            finite = dgm[np.isfinite(dgm[:, 1])]
            dgms[dim] = finite
        return dgms

    def compare_diagrams(dgm_semel, dgm_ripser):
        """Compute bottleneck distance per dimension."""
        dims = sorted(set(dgm_semel.keys()) | set(dgm_ripser.keys()))
        result = {}
        for d in dims:
            a = dgm_semel.get(d, np.empty((0, 2)))
            b = dgm_ripser.get(d, np.empty((0, 2)))
            if len(a) == 0 and len(b) == 0:
                result[d] = 0.0
            else:
                result[d] = bottleneck(a, b)
        return result

    # ── test case generators ─────────────────────────────────────────────

    def make_circle(n):
        P = pc.circle(n, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
        P = (P - P.min()) / (P.max() - P.min())
        return P

    def make_two_circles(n):
        n_inner = n // 3
        n_outer = n - n_inner
        P_inner = pc.circle(n_inner, cx=0.0, cy=0.0, rmin=0.35, rmax=0.45, method=1)
        P_outer = pc.circle(n_outer, cx=0.0, cy=0.0, rmin=0.85, rmax=0.95, method=1)
        P = np.vstack([P_inner, P_outer])
        P = (P - P.min()) / (P.max() - P.min())
        return P

    def make_random(n):
        P = pc.random(n, 2)
        P = (P - P.min()) / (P.max() - P.min())
        return P

    TOPOLOGIES = [
        ('circle',      make_circle,      0.5, 'H1=1 cycle'),
        ('two_circles', make_two_circles, 0.5, 'H1=2 cycles'),
        ('random',      make_random,      0.3, 'H1=noise'),
    ]

    dim_max = 2

    # ── CORRECTNESS SECTION ──────────────────────────────────────────────

    print('=' * 60)
    print('CORRECTNESS COMPARISON  (n = 200)')
    print('=' * 60)

    np.random.seed(42)
    N_CORRECT = 200

    correctness_data = []

    for name, gen_fn, r_max, expected in TOPOLOGIES:
        P = gen_fn(N_CORRECT)

        q_semel, t_semel = run_semel(P, dim_max=dim_max, r_max=r_max)
        q_ripser, t_ripser = run_ripser(P, maxdim=dim_max - 1, thresh=r_max)

        dgm_s = semel_to_diagrams(q_semel, r_max)
        dgm_r = ripser_to_diagrams(q_ripser)
        bn = compare_diagrams(dgm_s, dgm_r)

        correctness_data.append({
            'name': name, 'expected': expected,
            'P': P, 'dgm_s': dgm_s, 'dgm_r': dgm_r,
            'bn': bn, 'r_max': r_max,
            't_semel': t_semel, 't_ripser': t_ripser,
        })

        print('\n--- %s (expected: %s) ---' % (name, expected))
        print('  semel  %.4fs   ripser %.4fs' % (t_semel, t_ripser))
        for d in sorted(set(list(dgm_s.keys()) + list(dgm_r.keys()))):
            ns = len(dgm_s.get(d, []))
            nr = len(dgm_r.get(d, []))
            mass_s = float(np.sum(dgm_s[d][:, 1] - dgm_s[d][:, 0])) if d in dgm_s and len(dgm_s[d]) > 0 else 0.0
            mass_r = float(np.sum(dgm_r[d][:, 1] - dgm_r[d][:, 0])) if d in dgm_r and len(dgm_r[d]) > 0 else 0.0
            print('  H%d: semel %3d feat (mass %.4f)  ripser %3d feat (mass %.4f)  bottleneck %.6f'
                  % (d, ns, mass_s, nr, mass_r, bn.get(d, 0.0)))

    # ── correctness figure: panel grid (rows=topology, cols=semel/ripser/overlay) ──

    n_topo = len(TOPOLOGIES)
    fig, axes = plt.subplots(n_topo, 3, figsize=(15, 5 * n_topo))
    colors = {0: 'C0', 1: 'C1'}

    for row, cd in enumerate(correctness_data):
        r_max = cd['r_max']
        diag = [0, r_max]

        # semel diagram
        ax = axes[row, 0]
        for d, dgm in cd['dgm_s'].items():
            if len(dgm) > 0:
                ax.scatter(dgm[:, 0], dgm[:, 1], s=20, c=colors.get(d, 'C2'),
                           label='H%d' % d, alpha=0.7)
        ax.plot(diag, diag, 'k--', lw=0.5)
        ax.set_xlim([0, r_max * 1.1])
        ax.set_ylim([0, r_max * 1.1])
        ax.set_title('%s — semel' % cd['name'])
        ax.set_xlabel('birth'); ax.set_ylabel('death')
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        # ripser diagram
        ax = axes[row, 1]
        for d, dgm in cd['dgm_r'].items():
            if len(dgm) > 0:
                ax.scatter(dgm[:, 0], dgm[:, 1], s=20, c=colors.get(d, 'C2'),
                           label='H%d' % d, alpha=0.7)
        ax.plot(diag, diag, 'k--', lw=0.5)
        ax.set_xlim([0, r_max * 1.1])
        ax.set_ylim([0, r_max * 1.1])
        ax.set_title('%s — ripser' % cd['name'])
        ax.set_xlabel('birth'); ax.set_ylabel('death')
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        # overlay
        ax = axes[row, 2]
        for d, dgm in cd['dgm_s'].items():
            if len(dgm) > 0:
                ax.scatter(dgm[:, 0], dgm[:, 1], s=30, marker='o', facecolors='none',
                           edgecolors=colors.get(d, 'C2'), label='semel H%d' % d, alpha=0.7)
        for d, dgm in cd['dgm_r'].items():
            if len(dgm) > 0:
                ax.scatter(dgm[:, 0], dgm[:, 1], s=20, marker='x',
                           c=colors.get(d, 'C2'), label='ripser H%d' % d, alpha=0.7)
        ax.plot(diag, diag, 'k--', lw=0.5)
        ax.set_xlim([0, r_max * 1.1])
        ax.set_ylim([0, r_max * 1.1])
        bn_str = ', '.join('H%d=%.4f' % (d, v) for d, v in sorted(cd['bn'].items()))
        ax.set_title('%s — overlay\nbottleneck: %s' % (cd['name'], bn_str), fontsize=9)
        ax.set_xlabel('birth'); ax.set_ylabel('death')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('output/comparison_ripser_correctness.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    # ── feature counts bar chart ─────────────────────────────────────────

    fig, axes = plt.subplots(1, n_topo, figsize=(5 * n_topo, 4))
    for i, cd in enumerate(correctness_data):
        ax = axes[i]
        dims = sorted(set(list(cd['dgm_s'].keys()) + list(cd['dgm_r'].keys())))
        x = np.arange(len(dims))
        counts_s = [len(cd['dgm_s'].get(d, [])) for d in dims]
        counts_r = [len(cd['dgm_r'].get(d, [])) for d in dims]
        w = 0.35
        ax.bar(x - w / 2, counts_s, w, label='semel', color='C0')
        ax.bar(x + w / 2, counts_r, w, label='ripser', color='C1')
        ax.set_xticks(x)
        ax.set_xticklabels(['H%d' % d for d in dims])
        ax.set_title('%s (%s)' % (cd['name'], cd['expected']))
        ax.set_ylabel('feature count')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('output/comparison_ripser_features.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    # ── SCALING SECTION ──────────────────────────────────────────────────

    print('\n' + '=' * 60)
    print('SCALING COMPARISON')
    print('=' * 60)

    SIZES = [50, 100, 200, 400, 800, 1200, 1600, 2000]

    timing = {name: {'semel': [], 'ripser': []}
              for name, _, _, _ in TOPOLOGIES}

    for name, gen_fn, r_max, _ in TOPOLOGIES:
        print('\n--- %s ---' % name)
        print('%6s  %12s  %12s' % ('n', 'semel', 'ripser'))
        for n in SIZES:
            np.random.seed(42)
            P = gen_fn(n)
            _, ts = run_semel(P, dim_max=dim_max, r_max=r_max)
            _, tr = run_ripser(P, maxdim=dim_max - 1, thresh=r_max)
            timing[name]['semel'].append(ts)
            timing[name]['ripser'].append(tr)
            print('%6d  %11.4fs  %11.4fs' % (n, ts, tr))

    # ── scaling figure: log-log timing plot ──────────────────────────────

    fig, axes = plt.subplots(1, n_topo, figsize=(6 * n_topo, 5))
    for i, (name, _, _, _) in enumerate(TOPOLOGIES):
        ax = axes[i]
        ax.loglog(SIZES, timing[name]['semel'], 'o-', label='semel', color='C0')
        ax.loglog(SIZES, timing[name]['ripser'], 's-', label='ripser', color='C1')
        ax.set_xlabel('Number of points')
        ax.set_ylabel('Wall-clock time (s)')
        ax.set_title('%s — scaling' % name)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('output/comparison_ripser_scaling.jpg', dpi=150, bbox_inches='tight')
    plt.close()

    print('\nOutputs saved to output/comparison_ripser_*.jpg')
