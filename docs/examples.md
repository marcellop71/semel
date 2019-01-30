# Examples

All examples are standalone Python scripts in `examples/`. Figures are saved to `examples/output/` as JPEG files.

| Script | Description |
| --- | --- |
| `tutorial.py` | Self-contained demo: circle, random cloud, time series, coboundary operators |
| `comparison_ripser.py` | Correctness + scaling benchmark vs ripser (bottleneck distance, log-log timing) |
| `profile_stages.py` | Pipeline profiling: per-stage timing at multiple point cloud sizes |
| `chirikov.py` | Chirikov standard map on the torus, periodic Delaunay, 3D visualisation |
| `timeseries_rolling.py` | Rolling-window cohomology on synthetic signals (sin, harmonics) |
| `ladic_cohomology.py` | Integral and l-adic cohomology: F_p vs Z comparison, torsion detection |
| `local_coefficients.py` | Rank-1 local coefficients: trivial vs twisted transport weights on edges |
| `grid_coboundary.py` | Manual simplices on a grid, coboundary sparsity, torch cocycle optimisation |
| `sheaf_cohomology.py` | Rank-r sheaf cohomology: identity sheaf, twisted transport, sheaf Laplacian |

Run examples with:

```bash
python examples/tutorial.py
python examples/ladic_cohomology.py
```
