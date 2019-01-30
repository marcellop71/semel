# semel

*semel* is a C library (with Python and Rust bindings) for persistent cohomology on
simplicial complexes.
A _simplex_ is something "folded-once", _semel_-_plico_ from Latin,
something without com _plica_ tions.

## Documentation

| | |
|---|---|
| [Approach](docs/approach.md) | Delaunay pipeline, filtration, clearing algorithm, Delaunay vs Rips, manifold support, Hodge Laplacian, limitations, references |
| [Coefficients](docs/coefficients.md) | F_p (finite field), Z (integral / l-adic), and rank-1 local coefficient systems |
| [Sheaf](docs/sheaf.md) | Rank-r cellular sheaf cohomology: transport matrices, unrolling, sheaf Laplacian |
| [Build](docs/build.md) | Nix build, development shell, Python setup, dependencies, environment variables |
| [Usage](docs/usage.md) | Python API, workflow, persistence diagrams, coboundary operators, Hodge Laplacian, time series |
| [Examples](docs/examples.md) | Runnable scripts: tutorial, ripser comparison, l-adic cohomology, local coefficients, Chirikov map, time series |
| [Architecture](docs/architecture.md) | Source tree layout (C and Python) |

## Quick start

```bash
nix build
nix develop
source .venv/bin/activate  # or: python -m venv .venv --system-site-packages && source .venv/bin/activate
export SEMEL_LIB_PATH=./result/lib LD_LIBRARY_PATH=./result/lib:$LD_LIBRARY_PATH SEMEL_ZLOG_CONF=./result/etc/semel/zlog.conf
pip install ./wrap/python
```

```python
from semel import Semel, PointCloud

s = Semel()
P = PointCloud().circle(100, cx=0, cy=0, rmin=0.9, rmax=1.1, method=1)
P = (P - P.min()) / (P.max() - P.min())

m = {'type': s.MANIFOLD_FLAT_UNBOUNDED, 'distance': s.DISTANCE_L2, 'border': []}
s.add_point_cloud(0, m, P)
s.add_complex(0)                    # F_p coefficients (default F_5)
# s.add_complex_ladic(0, l=0)       # or: Z coefficients (integral cohomology)
s.generate_delaunay(0)
s.build_complex(0, dim_max=2, r_max=0.5)

q = s.cohomology(0, measures_only=False)
print(q['1d']['mass'])              # persistence mass per dimension
print(q['intervals'])               # birth-death pairs

s.unload_complex(0)
s.unload_point_cloud(0)
```
## Rust bindings

The `wrap/rust/` directory contains two crates:

- **semel-sys** -- raw FFI declarations mapping 1:1 to the C API in `semel.h`.
- **semel** -- safe, idiomatic wrapper with RAII resource management.

### Building

```bash
export SEMEL_LIB_PATH=./result/lib LD_LIBRARY_PATH=./result/lib:$LD_LIBRARY_PATH
cd wrap/rust/semel && cargo build
```

### Usage

```rust
use semel::{SemelSlot, Filtration};

let mut slot = SemelSlot::new();
let mut points = vec![1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0];
slot.add_point_cloud(4, 2, &mut points).unwrap();
slot.add_complex(5, 1).unwrap();
slot.generate_delaunay().unwrap();
slot.build_complex(2, 2.0, Filtration::Alpha).unwrap();
let result = slot.cohomology(false).unwrap();
println!("{}", result);
// slot is automatically cleaned up when it goes out of scope
```

### API overview

**`SemelSlot`** -- RAII handle for a point-cloud / complex slot. Unloads the complex and point cloud on drop.

| Method | Description |
|---|---|
| `add_point_cloud`, `add_point_cloud_with` | Load points (flat/manifold, L2/cosine/surface metrics) |
| `add_complex`, `add_complex_ladic` | Attach cochain complex (F_p or integral/l-adic coefficients) |
| `generate_delaunay` | Delaunay triangulation via qhull |
| `build_complex` | Build filtered simplicial complex (Alpha or Rips) |
| `add_simplices` | Add simplices directly (bypass triangulation) |
| `add_graph` | Load a weighted edge list |
| `cohomology` | Persistent cohomology (returns JSON) |
| `export_cocycles`, `export_complex`, `export_coboundary` | JSON export |
| `laplace_eigen` | Hodge Laplacian eigendecomposition |
| `get_scales` | Filtration scales |
| `set_local_system` | Rank-1 local system (transport weights on edges) |
| `set_sheaf` | Rank-r sheaf transport matrices |
| `sheaf_cohomology` | Persistent sheaf cohomology |
| `sheaf_laplace_eigen` | Sheaf Laplacian eigendecomposition |

**`TimeSeriesSlot`** -- RAII handle for time-series data. Unloads on drop.

| Method | Description |
|---|---|
| `add_time_series` | Load time-series data |
| `time_delay_embeddings` | Compute delay embeddings, returns a `SemelSlot` for the resulting point cloud |

**Free functions** -- `entropy_dist`, `dist_words_seq_bin`, `nrps_seq_bin`, `nsrps`, `ds_logistic`, `unload_point_cloud_all`, `unload_complex_all`, `unload_time_series_all`.

**Enums** -- `Manifold` (`FlatUnbounded`, `SurfaceCylinder`, `SurfaceTorus`), `Distance` (`L2`, `Cosine`, ...), `Filtration` (`Alpha`, `Rips`).

## AI Assistance Disclosure

Parts of this repository were created with assistance from AI-powered coding tools, specifically Claude by Anthropic. Not all generated code may have been reviewed. Generated code may have been adapted by the author. Design choices, architectural decisions, and final validation were performed independently by the author.

## License

MIT -- see [LICENSE](LICENSE).
