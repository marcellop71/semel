# Usage

## Over F_p (default)

```python
from semel import Semel, PointCloud

semel = Semel()
pc = PointCloud()

# 1. generate a point cloud
P = pc.circle(100, cx=0.0, cy=0.0, rmin=0.9, rmax=1.1, method=1)
P = (P - P.min()) / (P.max() - P.min())

# 2. load it with a manifold specification
manifold = {
    'type': semel.MANIFOLD_FLAT_UNBOUNDED,
    'distance': semel.DISTANCE_L2,
    'border': []
}
semel.add_point_cloud(0, manifold, P)

# 3. create a cochain complex over F_5 and compute Delaunay
semel.add_complex(0)             # default: F_5
semel.generate_delaunay(0)

# 4. build the filtered complex and compute persistent cohomology
semel.build_complex(0, dim_max=2, r_max=0.5)
q = semel.cohomology(0, measures_only=False)

print(q['1d']['mass'])      # persistent mass per dimension
print(q['1d']['entropy'])   # persistent entropy per dimension
print(q['intervals'])       # birth-death pairs: {dim: [[b, d], ...]}
print(q['cocycles'])        # cocycle representatives

# 5. clean up
semel.unload_complex(0)
semel.unload_point_cloud(0)
```

## Over Z (integral / l-adic)

```python
semel.add_point_cloud(0, manifold, P)

# use add_complex_ladic instead of add_complex
semel.add_complex_ladic(0, l=0)  # l=0: full integral cohomology
                                  # l=2: 2-adic (keep only 2-power torsion)
semel.generate_delaunay(0)
semel.build_complex(0, dim_max=2, r_max=0.5)
q = semel.cohomology(0, measures_only=False)

print(q['1d']['mass'])      # same persistence mass
print(q['betti'])           # {0: 1, 1: 1}  free rank per dimension
print(q['torsion'])         # {1: [2, 4]}    invariant factors > 1 per dim
print(q['intervals'])       # {dim: [[b, d, torsion_order], ...]}
print(q['cocycles'])        # {dim: [[b, torsion_order], ...]}

semel.unload_complex(0)
semel.unload_point_cloud(0)
```

Each interval is `[birth, death, torsion_order]` where `torsion_order = 0`
means free (the cocycle dies by an exact relation) and `torsion_order > 1`
means the cocycle carries torsion of that order (it takes that many copies to
kill it). When `l > 0`, only torsion whose order is a power of l is reported.

See [coefficients](coefficients.md) for details on the two coefficient rings.

## With local coefficients (rank-1 local system)

```python
import numpy as np

semel.add_point_cloud(0, manifold, P)
semel.add_complex_ladic(0, l=0)       # or add_complex(0) for F_p
semel.generate_delaunay(0)

# define transport weights on edges
edges = np.array([[0, 1], [1, 2], [0, 2]], dtype=np.uint64)
weights = np.array([2, 1, 1], dtype=np.int64)
semel.set_local_system(0, edges, weights)

semel.build_complex(0, dim_max=2, r_max=0.5)
q = semel.cohomology(0, measures_only=False)
```

`set_local_system` must be called after `add_complex`/`add_complex_ladic` and
before `build_complex`. Edges not listed default to weight 1. Setting all
weights to 1 recovers standard (untwisted) cohomology.

See [coefficients](coefficients.md) for details on the coefficient rings and
local systems.

## With sheaf cohomology (rank-r local system)

```python
import numpy as np

semel.add_point_cloud(0, manifold, P)
semel.add_complex_ladic(0, l=0)       # or add_complex(0) for F_p
semel.generate_delaunay(0)

# define rank-2 transport matrices on edges
rank = 2
edges = np.array([[0, 1], [1, 2], [0, 2]], dtype=np.uint64)
matrices = np.array([
    [[1, 0], [0, 1]],   # identity
    [[0, -1], [1, 0]],  # rotation
    [[1, 0], [0, 1]],   # identity
], dtype=np.int64)
semel.set_sheaf(0, rank, edges, matrices)

semel.build_complex(0, dim_max=2, r_max=0.5)
q = semel.sheaf_cohomology(0, measures_only=False)

# sheaf Laplacian
q_lap = semel.sheaf_Laplacian(0, d=2)
```

`set_sheaf` must be called after `add_complex`/`add_complex_ladic` and before
`build_complex`. Edges not listed default to identity transport. See
[sheaf](sheaf.md) for full documentation.

## Workflow

```
Point cloud -> add_point_cloud()
                     |
          add_complex() or add_complex_ladic()
                     |
            generate_delaunay()
                     |
          (optional) set_local_system() or set_sheaf()
                     |
               build_complex(dim_max, r_max)
                     |
      +--------------+------------------+-----------------+
      |              |                  |                 |
 cohomology()  export_complex()  export_coboundary()  HodgeLaplace()  sheaf_cohomology()  sheaf_Laplacian()
      |              |                  |                 |
 mass/entropy   simplex data      sparse matrix    eigenvalues
 + intervals    + cocycles                         + eigenvectors
      |
 unload_complex() -> unload_point_cloud()
```

## Time series

```python
from semel import TimeSeries

signal = [...]  # 1D time series
ts = TimeSeries(points=signal)
P = ts.embedding(dim=2, tau=15, method=0)  # sliding-window embedding
# then pass P through the standard pipeline above
```

## Persistence diagrams

`cohomology()` returns summary statistics by default. With `measures_only=False`,
it also returns the full persistence barcode:

```python
q = semel.cohomology(0, measures_only=False)

q['1d']['mass']       # {0: 1.23, 1: 0.45}  total lifespan per dimension
q['1d']['entropy']    # {0: 0.98, 1: 0.67}  Shannon entropy of normalised lifespans
q['intervals']        # {0: [[b,d], ...], 1: [[b,d], ...]}  finite bars (birth, death)
q['cocycles']         # {0: [[b], ...], 1: [[b], ...]}      infinite bars (birth only)
```

With `add_complex_ladic`, the output additionally contains:

```python
q['betti']            # {0: 1, 1: 1}        free rank per dimension
q['torsion']          # {1: [2]}             invariant factors > 1 per dim
q['intervals']        # {0: [[b,d,t], ...]}  t = torsion order (0 = free)
q['cocycles']         # {0: [[b,t], ...]}    t = torsion order (0 = free)
q['l']                # 2  (only present when l > 0)
```

**Intervals** are features that are born and later die (finite bars).
**Cocycles** are features that are born and survive to `r_max` (infinite bars).
To merge them into standard persistence diagrams:

```python
from semel import coalesce

dgm = coalesce(q['cocycles'], q['intervals'], r_max=0.5)
# dgm[dim] = numpy array of [birth, death] pairs, with death=r_max for infinite bars
```

Visualization is provided by the `Graph` class:

```python
from semel import Graph

fig, axes = plt.subplots(1, 2)
Graph().persistent_barcode(axes[0], dim=1, r_max=0.5, q=q)   # barcode plot
Graph().persistent_diagram(axes[1], dim=1, r_max=0.5, q=q)   # birth-death diagram
```

The diagram supports horizontal mode (`horizontal=True`, plotting birth vs lifespan)
and weighted mode (`weighted=True`, circle area proportional to lifespan).

Additional tools in `measures.py`:
- `persistent_1d_distribution(d)` -- normalised lifespan distribution and total mass
- `persistent_landscape(t, pi)` -- persistence landscape (top-4 envelope functions)
- `transport(A, B)` -- optimal transport distance between two diagrams (requires `pot`)

## Coboundary operators

The coboundary operator delta^k : C^k -> C^{k+1} is available as a sparse matrix:

```python
q = semel.export_coboundary(0, 0)   # delta^0: vertices -> edges

q['m']               # number of (dim+1)-simplices (rows)
q['n']               # number of dim-simplices (columns)
q['matrix_sparse']   # list of [row, col, value] triplets
```

Values are integers in the symmetric range of F_p (e.g. {-2,-1,0,1,2} for p=5).
The sign convention is (-1)^i where i is the position of the omitted vertex in the
sorted vertex list of the coface -- the standard simplicial coboundary sign.

## Hodge Laplacian

The Hodge Laplacian Delta_k on k-cochains decomposes into an "up" and "down" part:

```
Delta_k = delta_{k-1} delta_{k-1}^T + delta_k^T delta_k
```

where delta_k is the coboundary operator from C^k to C^{k+1}. Semel computes this
via GSL BLAS (dense matrix multiply) and returns the eigendecomposition:

```python
q = semel.HodgeLaplace(0, d=2)          # all eigenvalues (dense, O(n_k^3))
q = semel.HodgeLaplace(0, d=2, nev=10)  # 10 smallest eigenvalues (sparse ARPACK)

q['eigenvalues']     # list of arrays, one per Laplacian
q['eigenvectors']    # list of matrices (None when nev > 0)
```

`HodgeLaplace(k, d)` returns `d-1` Laplacians. Index convention: `eigenvalues[0]`
and `eigenvectors[0]` correspond to the 1-Laplacian (on edges), `eigenvalues[1]` to
the 2-Laplacian (on triangles), and so on. The number of zero eigenvalues equals the
Betti number of that dimension (by the Hodge theorem).

When `nev > 0`, the sparse Lanczos path computes only the `nev` smallest eigenvalues
using Lanczos iteration with full reorthogonalization. This avoids assembling the
dense Laplacian matrix. Eigenvectors are not returned via the sparse path.
