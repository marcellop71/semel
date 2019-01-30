# Approach

Semel computes persistent cohomology through a Delaunay-based pipeline that
differs fundamentally from the Vietoris-Rips approach used by libraries like
Ripser.

## The pipeline

```
Point cloud (n points in R^d)
        |
  distance matrix  O(n^2)
        |
  Delaunay triangulation  O(n log n) in 2D
        |
  filtered simplicial complex
  (top-down face decomposition with radius filtration)
        |
  persistent cohomology
  (clearing algorithm over finite fields or integers, via FLINT)
        |
  cocycles + intervals + mass + entropy
```

**1. Delaunay seed complex.**
Instead of considering all possible simplices (as Rips does), semel starts from
the Delaunay triangulation of the point cloud. In 2D this produces O(n)
triangles — a sparse skeleton that already captures the topology. The Delaunay
triangulation is computed in C via qhull.

**2. Top-down filtration.**
Given the Delaunay simplices, semel decomposes each into all its faces
recursively and assigns each simplex a filtration value. Two modes are
available (selected via `filtration` parameter, default: alpha):

- **Alpha** (`FILTRATION_ALPHA = 0`): `radius(sigma) = circumradius(sigma)`.
  This is the true alpha complex filtration (Čech on Delaunay): a simplex
  enters when its circumball radius reaches r.
- **Rips** (`FILTRATION_RIPS = 1`): `radius(sigma) = 0.5 * max pairwise distance`.
  A simplex enters when its vertices are pairwise within diameter 2r.

For edges both modes coincide (circumradius = 0.5 * edge length).
Simplices are then ordered by (radius, lexicographic serialization) to build
the filtered complex.

**3. Persistent cohomology with clearing.**
The cohomology algorithm processes simplices in filtration order. For each
simplex sigma_k:

- If no existing cocycle has sigma_k in its coboundary: **birth** — a new
  cocycle is created.
- If one or more cocycles have sigma_k in their coboundary: **death** — the
  latest-born cocycle is killed (creating an interval), and earlier cocycles
  are reduced via row operations.

Two coefficient rings are supported (see [coefficients](coefficients.md)).

**4. Output.**
Beyond persistence diagrams, semel returns:

- **Cocycle representatives**: explicit cochains (F_p or Z coefficients)
- **Persistence mass**: total lifespan per dimension (sum of death - birth)
- **Persistence entropy**: Shannon entropy of the normalised lifespans
- **Coboundary operators**: sparse matrices (for Hodge Laplacians, etc.)
- **Betti numbers** (integral path): free rank per dimension
- **Torsion invariant factors** (integral path): list of orders > 1 per dimension

## Why Delaunay instead of Rips in low dimension

| | Semel (Delaunay) | Ripser (Vietoris-Rips) |
|---|---|---|
| **Simplices** | O(n) in 2D from Delaunay | O(n^{k+1}) for k-simplices, all pairs considered |
| **Filtration** | circumradius (alpha, default) or 0.5 * max pairwise dist (Rips) | max pairwise dist |
| **Scaling** | Dominated by Delaunay O(n log n) + cohomology on sparse complex | O(n^3) in dense regime |
| **Crossover** | Slower for n < ~400 (fixed overhead) | Faster for small n, cubic blowup for large n |
| **Algebra** | Cohomology (yields cocycles) | Homology (yields cycles) |
| **Manifolds** | Euclidean, torus, cylinder with geodesic distances | Euclidean only |
| **Arithmetic** | Exact (finite field F_p, or integers Z) | Exact (Z/2Z) |

The key tradeoff: Rips is fast for small point clouds because it avoids the
Delaunay overhead, but its cubic scaling makes it impractical beyond ~1000
points. Semel pays a fixed cost for Delaunay triangulation and complex setup,
but the resulting sparse complex keeps the cohomology computation efficient at
scale.

## Manifold support

Semel supports point clouds on non-Euclidean manifolds:

| Manifold | Distance | Use case |
|---|---|---|
| `MANIFOLD_FLAT_UNBOUNDED` | Euclidean L2 | Standard R^d |
| `MANIFOLD_FLAT_UNBOUNDED` | Cosine | Angular data, embeddings |
| `MANIFOLD_SURFACE_CYLINDER` | Geodesic on S^1 x R | Periodic in one direction |
| `MANIFOLD_SURFACE_TORUS` | Geodesic on S^1 x S^1 | Periodic in both directions |

The distance matrix and filtration adapt to the manifold geometry, so
persistent cohomology correctly detects topology on periodic domains (e.g.
Chirikov map on the torus).

## Hodge Laplacian

The Hodge Laplacian Delta_k on k-cochains decomposes as:

```
Delta_k = delta_{k-1} delta_{k-1}^T + delta_k^T delta_k
```

where delta_k is the coboundary operator C^k -> C^{k+1}. The discrete Hodge
theorem guarantees that dim(ker Delta_k) = beta_k (the k-th Betti number), so
the multiplicity of the zero eigenvalue recovers the homological rank. The
nonzero eigenvalues encode further spectral geometry of the complex.

Semel provides two computational paths:

**Dense (full spectrum).** Assembles the coboundary matrices as dense GSL
matrices, forms Delta_k = D D^T + D^T D via BLAS `dgemm`, and computes the
full eigendecomposition via `gsl_eigen_symmv`. This is O(n_k^3) in the number
of k-simplices and returns all eigenvalues and eigenvectors.

**Sparse (partial spectrum).** For the `nev` smallest eigenvalues, semel uses
the Lanczos algorithm with full reorthogonalization. The coboundary operators
are stored in CSR (compressed sparse row) format and the Laplacian is never
assembled explicitly — only the matrix-vector product y = Delta_k x is
computed, via four sparse matrix-vector multiplications:

```
w1 = D_k^T x       w2 = D_k w1       (down part)
w3 = D_{k+1} x     w4 = D_{k+1}^T w3 (up part)
y  = w2 + w4
```

The Lanczos iteration builds an orthonormal Krylov basis Q and a tridiagonal
matrix T such that Delta_k ~ Q T Q^T. Full reorthogonalization (Gram-Schmidt
against all previous Lanczos vectors at each step) prevents spurious
eigenvalue copies from loss of orthogonality. The tridiagonal T is then
eigendecomposed via GSL to obtain Ritz values. The subspace dimension is set
to max(n_k/2, 20 * nev) to ensure convergence, falling back to the dense path
when nev >= n_k.

## Limitations

### 1. Ambient dimension ceiling

Delaunay triangulation has complexity O(n^{ceil(d/2)}) in d dimensions. In 2D
it's O(n log n), in 3D it's O(n^2), but by d=10 it's O(n^5) and qhull will
struggle or run out of memory. In practice, the library is effective for
ambient dimension d <= ~6. High-dimensional point clouds need to be projected
first (PCA, UMAP, etc.), which can distort topology.

Rips doesn't have this problem -- it works directly on a distance matrix
regardless of ambient dimension. That's the real tradeoff: semel is faster for
low-d, but Rips generalises to arbitrary metric spaces.

### 2. Requires coordinates, not just distances

Delaunay triangulation needs an embedding in R^d. You can't feed in an abstract
distance matrix (e.g. from a graph shortest-path metric, or a kernel matrix).
The manifold support (torus, cylinder) partially addresses this but only for
specific geometries. A Rips-based library only needs the pairwise distance
matrix.

### 3. Coefficient growth over Z

Each GCD elimination step does `alpha_i <- (c_j/g) * alpha_i - (c_i/g) * alpha_j`.
The multiplier `c_j/g` can be > 1, so cocycle coefficients grow with each
reduction. Content reduction (dividing all coefficients of a cocycle by their
GCD) is performed after each reduction step to keep coefficient bit-lengths
bounded. This does not change the cohomology class. For large complexes,
performance still degrades compared to the F_p path where coefficients are
bounded by p, but content reduction keeps the growth manageable.

### 4. Dense Hodge Laplacian (full spectrum)

The full eigendecomposition (`nev=0`) uses GSL dense matrices, so it's O(n_k^3)
in the number of k-simplices. For the `nev` smallest eigenvalues, a sparse
Lanczos path is available (`nev > 0`) that uses only sparse matrix-vector
products, avoiding explicit dense matrix assembly. Eigenvectors are not returned
via the sparse path.

### 5. Not a full persistence module decomposition over Z

Over a field, the persistence module decomposes into interval modules (the
barcode) -- this is a theorem. Over Z, the structure is richer: you get
interval modules for the free part plus torsion modules whose structure depends
on the ring. The current implementation captures the torsion order at each
death event (from the pivot), which is correct for identifying torsion in
individual persistence pairs, but it doesn't compute the full Smith normal form
of the boundary matrix or decompose the persistence module into its
indecomposable summands over Z. For most practical purposes the per-interval
torsion order is what you want, but it's worth noting the theoretical gap.

### 6. Single-threaded

The entire pipeline is single-threaded. The cohomology loop is sequential by
algorithmic necessity (each death event depends on the current state of all
active cocycles). The Delaunay step is a single qhull call. The
complex-building step inserts into shared mutable data structures (Judy arrays)
that would require synchronization to access concurrently.

## References

- Eckmann, B. (1945). "Harmonische Funktionen und Randwertaufgaben in einem
  Komplex." *Commentarii Mathematici Helvetici*, 17(1), 240-255.
  Discrete Hodge theory: harmonic cochains represent cohomology classes.

- Dodziuk, J. (1976). "Finite-difference approach to the Hodge theory of
  harmonic forms." *American Journal of Mathematics*, 98(1), 79-104.
  Discrete Hodge theorem: ker(Delta_k) = H^k for simplicial complexes.

- Lanczos, C. (1950). "An iteration method for the solution of the eigenvalue
  problem of linear differential and integral operators." *Journal of Research
  of the National Bureau of Standards*, 45(4), 255-282.
  The Lanczos algorithm for symmetric eigenvalue problems.

- Parlett, B. N. and Scott, D. S. (1979). "The Lanczos algorithm with
  selective orthogonalization." *Mathematics of Computation*, 33(145), 217-238.
  Reorthogonalization strategies for numerical stability in Lanczos.

- Golub, G. H. and Van Loan, C. F. (1996). *Matrix Computations* (3rd ed.).
  Johns Hopkins University Press.
  Standard reference for Lanczos iteration (Section 9.3) and sparse eigensolvers.

- de Silva, V., Morozov, D., and Vejdemo-Johansson, M. (2011). "Dualities in
  persistent (co)homology." *Inverse Problems*, 27(12), 124003.
  Duality between persistent homology and cohomology; algorithmic equivalence.

- Chen, C. and Kerber, M. (2011). "Persistent homology computation with a
  twist." *Proceedings of the 27th European Workshop on Computational Geometry*,
  197-200.
  The clearing algorithm for persistent cohomology.

- de Silva, V. and Carlsson, G. (2004). "Topological estimation using witness
  complexes." *Symposium on Point-Based Graphics*, 157-166.
  Sparse complexes (witness, Delaunay) as alternatives to Vietoris-Rips.
