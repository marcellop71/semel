# Sheaf cohomology

Semel supports **rank-r cellular sheaves** on simplicial complexes. Each edge
carries an r x r integer transport matrix instead of a scalar weight. The
sheaf coboundary is built via **unrolling**: each simplex becomes r virtual
simplices, and the existing persistent cohomology engine runs unmodified on
the unrolled complex.

## Mathematical background

A **cellular sheaf** F on a simplicial complex K assigns:
- A vector space F(sigma) = Z^r (the **stalk**) to each simplex sigma
- A linear map F(sigma <= tau) : F(sigma) -> F(tau) (the **restriction map**)
  to each face relation sigma <= tau

The **sheaf coboundary** delta_F : C^p(K; F) -> C^{p+1}(K; F) generalises
the simplicial coboundary by replacing scalar signs with r x r blocks.

### Steenrod convention

For a coface tau = [v_0, v_1, ..., v_{p+1}] with sorted vertices:
- **Position-0 face** (removing v_0): the r x r block is `sign * T(v_0, v_1)`
  where T is the transport matrix for edge (v_0, v_1)
- **All other faces**: the r x r block is `sign * I_r` (identity)

This generalises the rank-1 local system: setting r=1 with 1x1 scalar matrices
recovers `set_local_system` exactly.

### Unrolling scheme

Each simplex sigma with ID `sid` becomes r virtual simplices with IDs
`sid*r`, `sid*r+1`, ..., `sid*r+(r-1)`. All share the same dimension,
filtration age, and radius. The coboundary of virtual simplex `sid*r+i`
encodes column i of the block structure.

## API

### Python

```python
import numpy as np
from semel import Semel

s = Semel()

# ... add point cloud, complex, build complex as usual ...

# Define rank-2 transport matrices on edges
rank = 2
edges = np.array([[0, 1], [1, 2], [0, 2]], dtype=np.uint64)
matrices = np.array([
    [[1, 0], [0, 1]],   # identity on edge (0,1)
    [[0, -1], [1, 0]],  # rotation on edge (1,2)
    [[1, 0], [0, 1]],   # identity on edge (0,2)
], dtype=np.int64)

s.set_sheaf(0, rank, edges, matrices)
s.build_complex(0, dim_max=2, r_max=0.5)

# Persistent sheaf cohomology
q = s.sheaf_cohomology(0, measures_only=False)
print(q['1d']['mass'])
print(q['intervals'])

# Sheaf Laplacian
q = s.sheaf_Laplacian(0, d=2)
print(q['eigenvalues'])
```

### C

```c
int32_t _set_sheaf(uint64_t k, uint32_t rank,
    uint32_t n_edges, uint64_t * endpoints, int64_t * matrices);

char * _sheaf_cohomology(uint64_t k, uint8_t measures_only);

int32_t _sheaf_Laplace_eigen(uint64_t k, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v);
```

## Three computation modes

| Mode | Function | Description |
|------|----------|-------------|
| Persistent cohomology | `sheaf_cohomology()` | Birth-death intervals, mass, entropy |
| Hodge Laplacian (dense) | `sheaf_Laplacian(k, d)` | Full eigendecomposition |
| Hodge Laplacian (sparse) | `sheaf_Laplacian(k, d, nev=10)` | Smallest nev eigenvalues via Lanczos |

## Relationship to local_system

| | `set_local_system` | `set_sheaf` |
|---|---|---|
| **Rank** | 1 (scalar) | r (matrix) |
| **Storage** | `int64_t` per edge | `int64_t[r*r]` per edge |
| **Coboundary** | Modified in-place during `build_complex` | Built in temporary unrolled complex |
| **When r=1** | -- | Equivalent to `set_local_system` |

## Coefficient rings

Sheaf cohomology works with both coefficient paths:

- **F_p** (`add_complex`): field coboundary, detects Betti numbers
- **Z** (`add_complex_ladic`): integral coboundary, detects Betti numbers + torsion

## Edge cases

- Edges not in the `edges` array default to identity transport (I_r)
- `sheaf_cohomology` returns NULL if no sheaf is set (rank=0)
- Identity sheaf (all matrices = I_r) gives Betti numbers = r x standard Betti numbers

## Correctness invariants

| Test | Expected |
|------|----------|
| Rank-1 sheaf with scalar [w] | Matches `set_local_system` with weight w |
| Identity rank-r sheaf | Betti numbers = r x standard Betti numbers |
| delta^2 = 0 (flat system) | Trivial transport gives consistent coboundary |
| Virtual ID ordering | `sid*r+c` preserves filtration order |
