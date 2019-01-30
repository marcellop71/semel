//! Raw FFI declarations for libsemel (persistent cohomology library).
//!
//! These map 1:1 to the C API in `semel.h`. Prefer the safe wrapper
//! in the `semel` crate over calling these directly.

use std::os::raw::c_char;

// ── Manifold types (enums.h) ────────────────────────────────────────
pub const MANIFOLD_FLAT_UNBOUNDED: u32 = 0;
pub const MANIFOLD_SURFACE_CYLINDER: u32 = 1;
pub const MANIFOLD_SURFACE_TORUS: u32 = 2;

// ── Distance types (enums.h) ────────────────────────────────────────
pub const DISTANCE_L2: u32 = 0;
pub const DISTANCE_COSINE: u32 = 1;
pub const DISTANCE_SURFACE_CYLINDER: u32 = 2;
pub const DISTANCE_SURFACE_TORUS: u32 = 3;

// ── Filtration types (enums.h) ──────────────────────────────────────
pub const FILTRATION_ALPHA: u8 = 0;
pub const FILTRATION_RIPS: u8 = 1;

// ── Cohomology algorithm types (enums.h) ────────────────────────────
pub const COHOMOLOGY_FIELD: u32 = 0;
pub const COHOMOLOGY_INTEGRAL: u32 = 1;

extern "C" {
    // ── Global lifecycle ────────────────────────────────────────────
    pub fn s_init();
    pub fn s_free();

    // ── Point cloud management ──────────────────────────────────────
    pub fn _get_point_cloud_index() -> u64;

    pub fn _add_point_cloud_raw(
        k: u64,
        manifold: u32,
        distance: u32,
        n: u32,
        d: u32,
        points: *mut f64,
        nb: u32,
        border: *mut f64,
    ) -> i32;

    pub fn _unload_point_cloud(k: u64) -> i32;
    pub fn _unload_point_cloud_all() -> i32;

    // ── Complex management ──────────────────────────────────────────
    pub fn _add_complex(k: u64, fp: u32, fd: u32) -> i32;
    pub fn _add_complex_ladic(k: u64, l: u32) -> i32;
    pub fn _unload_complex(k: u64) -> i32;
    pub fn _unload_complex_all() -> i32;

    // ── Triangulation + filtration ──────────────────────────────────
    pub fn _generate_delaunay(k: u64) -> i32;
    pub fn _build_complex(k: u64, dim_sk: u32, rmax: f64, filtration_type: u8) -> i32;

    // ── Cohomology ──────────────────────────────────────────────────
    /// Returns a malloc'd JSON string. Caller must free with `libc::free()`.
    pub fn _cohomology(k: u64, measures_only: u8) -> *mut c_char;

    // ── Export ───────────────────────────────────────────────────────
    /// Returns a malloc'd JSON string. Caller must free with `libc::free()`.
    pub fn _export_cocycles(k: u64) -> *mut c_char;
    /// Returns a malloc'd JSON string. Caller must free with `libc::free()`.
    pub fn _export_complex(k: u64) -> *mut c_char;
    /// Returns a malloc'd JSON string. Caller must free with `libc::free()`.
    pub fn _export_coboundary(k: u64, d: u32) -> *mut c_char;

    // ── Hodge Laplacian ─────────────────────────────────────────────
    pub fn _Laplace_eigen(
        k: u64,
        n: u32,
        nev: u32,
        na: *mut *mut u32,
        a: *mut *mut *mut f64,
        v: *mut *mut *mut f64,
    ) -> i32;

    // ── Filtration scales ───────────────────────────────────────────
    pub fn _get_scales(
        k: u64,
        nr: *mut u64,
        r: *mut *mut u64,
        precision: *mut f64,
    ) -> i32;

    // ── Time series ─────────────────────────────────────────────────
    pub fn _add_time_series_raw(k: u64, d: u32, n: u32, data: *mut f64) -> i32;
    pub fn _unload_time_series(k: u64) -> i32;
    pub fn _unload_time_series_all() -> i32;
    pub fn _time_delay_embeddings(k: u64, tau: u32, d: u32, k_p: *mut u64) -> i32;

    // ── Simplices (raw) ─────────────────────────────────────────────
    pub fn _add_simplices_raw(
        k: u64,
        type_s: u8,
        n_simplices: u32,
        verts_per_simplex: u32,
        indices: *mut u64,
    ) -> i32;

    // ── Graphs ──────────────────────────────────────────────────────
    pub fn _add_graph_raw(
        k: u64,
        n_edges: u32,
        endpoints: *mut u64,
        weights: *mut f64,
    ) -> i32;

    // ── Local systems + sheaf cohomology ────────────────────────────
    pub fn _set_local_system(
        k: u64,
        n_edges: u32,
        endpoints: *mut u64,
        weights: *mut i64,
    ) -> i32;

    pub fn _set_sheaf(
        k: u64,
        rank: u32,
        n_edges: u32,
        endpoints: *mut u64,
        matrices: *mut i64,
    ) -> i32;

    /// Returns a malloc'd JSON string. Caller must free with `libc::free()`.
    pub fn _sheaf_cohomology(k: u64, measures_only: u8) -> *mut c_char;

    pub fn _sheaf_Laplace_eigen(
        k: u64,
        n: u32,
        nev: u32,
        na: *mut *mut u32,
        a: *mut *mut *mut f64,
        v: *mut *mut *mut f64,
    ) -> i32;

    // ── Sequences + information theory ──────────────────────────────
    pub fn _nsrps(ns: u64, s: *mut u32) -> i32;
    pub fn _entropy_dist(np: u64, p: *mut f64, alpha: f64, e: *mut f64) -> i32;
    pub fn _dist_words_seq_bin(
        nx: u64,
        x: *mut u8,
        w: u64,
        np: *mut u64,
        p: *mut *mut f64,
    ) -> i32;
    pub fn _nrps_seq_bin(nx: u64, x: *mut u8) -> i32;

    // ── Dynamical systems ───────────────────────────────────────────
    pub fn ds_logistic(lambda_: f64, np: u64, p: *mut *mut f64) -> i32;
}
