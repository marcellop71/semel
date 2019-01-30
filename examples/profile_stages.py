"""Profile each stage of the semel pipeline separately for circle point clouds."""

import os
import sys
import time
import numpy as np

# Ensure the semel wrapper is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'wrap', 'python'))

from semel import Semel, PointCloud

if __name__ == '__main__':
    np.random.seed(42)

    semel_ctx = Semel()
    pc = PointCloud()

    MANIFOLD = {
        'type': 0,       # MANIFOLD_FLAT_UNBOUNDED
        'distance': 0,   # DISTANCE_L2
        'border': []
    }

    dim_max = 2
    r_max = 0.5

    SIZES = [500, 1000, 2000]

    STAGES = [
        'add_point_cloud_raw',
        'add_complex',
        'generate_delaunay',
        'build_complex',
        'cohomology',
        'cleanup',
    ]

    for n in SIZES:
        # Re-seed so each size gets the same base randomness pattern
        np.random.seed(42)

        # Generate circle point cloud
        P = pc.circle(n, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
        P = (P - P.min()) / (P.max() - P.min())
        P = np.ascontiguousarray(P, dtype=np.float64)

        pci = 0
        timings = {}

        # Stage 1: add_point_cloud_raw
        t0 = time.perf_counter()
        semel_ctx.add_point_cloud_raw(pci, MANIFOLD, P)
        t1 = time.perf_counter()
        timings['add_point_cloud_raw'] = t1 - t0

        # Stage 2: add_complex
        t0 = time.perf_counter()
        semel_ctx.add_complex(pci)
        t1 = time.perf_counter()
        timings['add_complex'] = t1 - t0

        # Stage 3: generate_delaunay (C-side qhull)
        t0 = time.perf_counter()
        semel_ctx.generate_delaunay(pci)
        t1 = time.perf_counter()
        timings['generate_delaunay'] = t1 - t0

        # Stage 4: build_complex (distance matrix + filtration + face/coface)
        t0 = time.perf_counter()
        semel_ctx.build_complex(pci, dim_max, r_max)
        t1 = time.perf_counter()
        timings['build_complex'] = t1 - t0

        # Stage 5: cohomology (persistent cohomology computation)
        t0 = time.perf_counter()
        q = semel_ctx.cohomology(pci, measures_only=False)
        t1 = time.perf_counter()
        timings['cohomology'] = t1 - t0

        # Stage 6: cleanup
        t0 = time.perf_counter()
        semel_ctx.unload_complex(pci)
        semel_ctx.unload_point_cloud(pci)
        t1 = time.perf_counter()
        timings['cleanup'] = t1 - t0

        total = sum(timings.values())

        # Print table
        print('=' * 64)
        print('  n = %d points (circle)' % n)
        print('=' * 64)
        print('  %-28s %10s %8s' % ('Stage', 'Time (s)', '% Total'))
        print('  ' + '-' * 50)
        for stage in STAGES:
            t = timings[stage]
            pct = 100.0 * t / total if total > 0 else 0.0
            print('  %-28s %10.6f %7.1f%%' % (stage, t, pct))
        print('  ' + '-' * 50)
        print('  %-28s %10.6f %7.1f%%' % ('TOTAL', total, 100.0))
        print()
