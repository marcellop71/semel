# Architecture

```
src/
  apis.c                 Public C API (raw pointer entry points)
  complex_simplicial.c   Complex construction (top-down, bottom-up, filtration)
  cohomology.c           Persistent cohomology with clearing (over F_p)
  cohomology_ladic.c     Persistent cohomology over Z (GCD-based elimination)
  complex_cochain.c      Coboundary operators, cochain arithmetic (F_p and Z)
  geometry.c             Distance metrics, circumradius, circumsphere
  simplex.c              Simplex data structure (JudySL-based)
  export.c               Result export (cohomology, complex, coboundary)
  sheaf.c               Sheaf unrolling (rank-r transport → virtual complex)
  hodge.c                Hodge Laplacian eigendecomposition
  pointclouds.c          Point cloud init, distance matrix
  serialization.c        Simplex serialization for hashing
  tools.c                Judy array utilities, bit operations
  ctx.c                  Global context management
  time.c                 Monotonic clock
  structs.c              Data structure init/free

include/
  structs.h              Core data structures (point_cloud_t, complex_cochain_t)
  semel.h                Function declarations
  enums.h                Manifold, distance, and cohomology algorithm enumerations
  defs.h                 Constants (PRECISION=100000, MAX_DIM=1024)
  common.h               System and library includes
  fwddecls.h             Forward type declarations
  ctx.h                  Global context structure

wrap/python/semel/
  apis.py                Python bindings (ctypes)
  measures.py            coalesce(), entropy(), persistence landscape
  pointclouds.py         Point cloud generators (circle, random, grid, sphere)
  graphs.py              Visualization (persistence diagrams, barcodes)
  timeseries.py          Time series embedding

wrap/rust/
  semel-sys/             Raw FFI crate
    build.rs             Linker configuration (SEMEL_LIB_PATH)
    src/lib.rs           extern "C" declarations + enum constants
  semel/                 Safe wrapper crate
    src/lib.rs           SemelSlot, TimeSeriesSlot, enums, free functions
```

## Rust bindings

The Rust wrapping lives in `wrap/rust/` and is split into two crates:

### semel-sys (FFI layer)

`semel-sys` contains raw `extern "C"` declarations that map 1:1 to the
underscore-prefixed C API in `semel.h` (e.g. `_add_point_cloud_raw`,
`_cohomology`, `_Laplace_eigen`). It also re-exports the enum constants
(`MANIFOLD_*`, `DISTANCE_*`, `FILTRATION_*`, `COHOMOLOGY_*`).

`build.rs` tells Cargo where to find `libsemel.so` — it reads
`SEMEL_LIB_PATH` or falls back to a default path.

### semel (safe wrapper)

`semel` provides an idiomatic Rust API on top of `semel-sys`:

- **`SemelSlot`** — RAII handle for a point-cloud / complex slot. Holds
  the slot index `k` and boolean flags (`has_point_cloud`, `has_complex`).
  On `Drop`, it calls `_unload_complex` then `_unload_point_cloud` in that
  order. All slot-scoped operations (add points, build complex, cohomology,
  export, Hodge Laplacian, sheaf, local systems) are methods on this type.

- **`TimeSeriesSlot`** — RAII handle for a time-series slot. On `Drop`, it
  calls `_unload_time_series`. The `time_delay_embeddings` method returns a
  `SemelSlot` wrapping the point cloud produced by the embedding.

- **Enums** — `Manifold`, `Distance`, `Filtration` are `#[repr(u32/u8)]`
  enums that replace bare integer constants.

- **Result types** — `EigenDecomposition` (eigenvalues + optional
  eigenvectors per dimension), `EntropyMeasures` (Shannon, collision, min,
  Renyi), `FiltrationScales` (scale values + precision).

- **Free functions** — `entropy_dist`, `dist_words_seq_bin`,
  `nrps_seq_bin`, `nsrps`, `ds_logistic` for operations not tied to a slot,
  plus `unload_*_all` for global cleanup.

- **JSON handling** — functions returning C `char *` (malloc'd JSON) are
  wrapped by `parse_and_free`, which copies into a `serde_json::Value` and
  immediately calls `libc::free`.

- **Global init** — `s_init()` is called exactly once via `std::sync::Once`
  (triggered by `SemelSlot::new()`, `TimeSeriesSlot::new()`, or any free
  function).

### Ownership boundary

The C library owns all heavy data (point clouds, simplicial complexes,
cochain complexes, time series) in global Judy-array maps keyed by slot
index. The Rust wrapper never holds pointers into C-managed memory beyond
the scope of a single method call — output arrays from `_Laplace_eigen`,
`_get_scales`, etc. are immediately copied into owned `Vec`s.

## Judy arrays

Nearly every mutable data structure in semel is a Judy array (libJudy). A Judy
array is a compressed digital trie that maps machine-word keys to machine-word
values with O(log_{256} n) access time and cache-line-aligned internal nodes.
Semel uses two variants:

- **JudyL** — uint64 keys. Used for simplex IDs, filtration indices, dimension
  counters, coboundary maps, cocycle storage, and the global context.
- **JudySL** — null-terminated byte-string keys. Used for simplex serialization
  lookups (vertex list -> ID) and filtration ordering (serialization provides
  lexicographic tiebreaking within each radius bucket).

### Why Judy arrays

The core data structures in persistent cohomology are sparse and dynamic:
cocycles gain and lose nonzero coefficients during row reduction, the inverted
index must track which cocycles touch each simplex, and the face/coface maps
are irregular. Hash tables would work but Judy arrays give two properties that
hash tables do not:

1. **Sorted iteration.** `JLF`/`JLN` (first/next) traverse keys in ascending
   order with no sorting step. This is used throughout: filtration-order
   iteration over simplices, dimension-order traversal of the complex, and
   ordered enumeration of cocycle coefficients during linear combination.

2. **Predictable memory.** Judy arrays allocate in cache-line-sized chunks and
   compress empty branches. For the sparse coefficient vectors in cohomology
   (which may have a few dozen nonzero entries out of thousands of simplices),
   this is substantially more compact than a pre-allocated dense array.

### Layout of key structures

```
semel_ctx                             (ctx.h)
  .P   JL: pointcloud_index -> point_cloud_t *
  .Cc  JL: complex_index -> complex_cochain_t *
  .S   JL: timeseries_index -> timeseries_t *

complex_simplicial                    (structs.h)
  .filtration    JL: radius -> { JSL: serialization -> simplex }
  .simplices     JL: id -> simplex
  .sx2id         JSL: serialization -> id
  .faces         JL: id -> { JL: face_id -> index }
  .cofaces       JL: id -> { JL: coface_id -> index }
  .scd[dim]      JL: counter -> id        (enumerate simplices per dimension)
  .sci[dim]      JL: id -> counter        (reverse index)

complex_cochain                       (structs.h)
  .elementary_coboundaries      JL: id -> { JL: coface_id -> fq_t * }
  .elementary_coboundaries_z    JL: id -> { JL: coface_id -> fmpz_t * }
  .cocycles                     JL: cocycle_id -> { JL: simplex_id -> fq_t * }
  .cocycles_z                   JL: cocycle_id -> { JL: simplex_id -> fmpz_t * }
  .intervals                    JL: birth_id -> death_id
  .intervals_torsion            JL: interval_id -> torsion_order
```

Nesting is common: the filtration is a JudyL of JudySL arrays (one per radius
bucket), each coboundary is a JudyL of JudyL arrays (one per simplex), and
each cocycle is a JudyL mapping simplex IDs to coefficient pointers.

### Deallocation

Because Judy arrays store `uint64_t` values that are often cast pointers to
FLINT field elements (`fq_t *`, `fmpz_t *`) or to inner Judy arrays,
deallocation must walk the structure and free contents before freeing the
array itself. `tools.c` provides typed recursive freers:

| Function | Structure |
|---|---|
| `Judy_free_1_a` | flat JL of pointers |
| `Judy_free_2_a` | nested JL of JL |
| `Judy_free_1_b` / `Judy_free_2_b` | flat / nested JL of `fq_t *` |
| `Judy_free_1_z` / `Judy_free_2_z` | flat / nested JL of `fmpz_t *` |
