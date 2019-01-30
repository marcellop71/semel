# Cohomology coefficients

Semel supports two coefficient rings for persistent cohomology, selected at
complex creation time.

## F_p (finite field — default)

```python
semel.add_complex(k)             # default: F_5
semel.add_complex(k, fp=7)       # F_7
semel.add_complex(k, fp=2, fd=1) # F_2
```

All arithmetic is exact, performed in a finite field F_p using FLINT. Row
reduction uses field division: `scale = -c_i / c_j`.

This detects Betti numbers and persistence intervals but is blind to torsion —
H^k(K; Z) may contain torsion subgroups (e.g. Z/2Z in RP^2) invisible over
any single F_p.

Output per interval: `[birth, death]`.

## Z (integral / l-adic)

```python
semel.add_complex_ladic(k, l=0)  # full integral cohomology
semel.add_complex_ladic(k, l=2)  # 2-adic (see below)
semel.add_complex_ladic(k, l=3)  # 3-adic
```

Row reduction uses GCD-based elimination instead of field division:

```
g = gcd(c_i, c_j)
alpha_i <- (c_j/g) * alpha_i - (c_i/g) * alpha_j
```

This eliminates the pivot simplex from alpha_i's coboundary while minimising
coefficient growth. When a cocycle dies, the absolute value of the pivot
coefficient is the **torsion order** of that interval:

- |pivot| = 1: the death is free (no torsion), `torsion_order = 0`
- |pivot| > 1: the interval carries torsion of that order

Output per interval: `[birth, death, torsion_order]`.

### l-adic filtering

With `l = 0` (the default), all torsion orders are reported. With `l > 0`,
the output is filtered to the **l-primary component**: only torsion orders
that are a power of l (i.e. l, l^2, l^3, ...) are kept. This corresponds to
localising the cohomology at the prime l — extracting the Sylow l-subgroup of
each torsion part.

The full integral reduction is always performed regardless of l; the filtering
happens at export time. For example, if the computed torsion orders are
{2, 3, 4, 6, 8}, then:

- `l=0` reports all: {2, 3, 4, 6, 8}
- `l=2` reports {2, 4, 8} (powers of 2; 3 and 6 are discarded)
- `l=3` reports {3} (powers of 3 only)

This is useful when you know which prime carries the torsion you care about
(e.g. Z/2Z in RP^2) and want to suppress noise from other primes.

## Local coefficients (rank-1 local systems)

Instead of constant coefficients C^k(K; R), semel can compute C^k(K; L)
where L is a rank-1 local system -- a representation pi_1(K) -> R*
specified by scalar transport weights on edges.

```python
edges = np.array([[0, 1], [1, 2], [0, 2]], dtype=np.uint64)
weights = np.array([2, 1, 1], dtype=np.int64)

semel.add_complex_ladic(k, l=0)   # or add_complex(k) for F_p
semel.generate_delaunay(k)
semel.set_local_system(k, edges, weights)  # after add_complex, before build_complex
semel.build_complex(k, dim_max=2, r_max=0.5)
```

The coboundary operator is twisted using the **Steenrod convention**: for a
coface tau = [v_0, v_1, ..., v_{p+1}] with sorted vertices, only the
position-0 face (removing the smallest vertex v_0) is multiplied by the
transport factor w(v_0, v_1). All other faces keep their standard sign
(-1)^i. This gives delta^2 = 0 when the local system is flat (cocycle
condition on triangles).

- Edges not in the `edges` array default to weight 1 (identity transport)
- Weight 1 on all edges recovers standard (untwisted) cohomology exactly
- Works with both F_p and Z coefficients
- All downstream code (persistent cohomology, row reduction, Hodge Laplacian)
  operates on stored coboundary coefficients and needs no changes

## Sheaf cohomology (rank-r local systems)

For higher-rank local systems, semel supports **rank-r cellular sheaves** where
each edge carries an r x r integer transport matrix. See [sheaf](sheaf.md) for
the full API and mathematical background.

```python
rank = 2
edges = np.array([[0, 1], [1, 2]], dtype=np.uint64)
matrices = np.array([
    [[1, 0], [0, 1]],
    [[0, -1], [1, 0]],
], dtype=np.int64)

semel.add_complex_ladic(k, l=0)   # or add_complex(k) for F_p
semel.generate_delaunay(k)
semel.set_sheaf(k, rank, edges, matrices)
semel.build_complex(k, dim_max=2, r_max=0.5)
q = semel.sheaf_cohomology(k, measures_only=False)
```

## Comparison

| | Finite field F_p | Integers Z |
|---|---|---|
| **Setup** | `add_complex(k)` | `add_complex_ladic(k, l=0)` |
| **Arithmetic** | Field division: `scale = -c_i / c_j` | GCD elimination: `alpha_i <- (c_j/g)*alpha_i - (c_i/g)*alpha_j` |
| **Default** | p=5, degree 1 | -- |
| **Detects** | Betti numbers (free rank) | Betti numbers + torsion subgroups |
| **Interval** | `[birth, death]` | `[birth, death, torsion_order]` |
| **Output** | mass, entropy, cocycles, intervals | same + `betti`, `torsion` invariant factors |
| **Performance** | Coefficients bounded by p | Coefficients can grow (arbitrary precision) |
