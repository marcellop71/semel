//! Safe RAII wrapper around the semel C library.
//!
//! Semel uses global state and is NOT thread-safe. All calls must happen
//! on a single thread. The global context is initialised exactly once
//! via `std::sync::Once`.
//!
//! # Example
//!
//! ```no_run
//! use semel::{SemelSlot, Filtration};
//!
//! let mut slot = SemelSlot::new();
//! let mut points = vec![1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0];
//! slot.add_point_cloud(4, 2, &mut points).unwrap();
//! slot.add_complex(5, 1).unwrap();
//! slot.generate_delaunay().unwrap();
//! slot.build_complex(2, 2.0, Filtration::Alpha).unwrap();
//! let result = slot.cohomology(false).unwrap();
//! println!("{}", result);
//! ```

use std::ffi::CStr;
use std::sync::Once;

use anyhow::{bail, Result};
use serde_json::Value;

pub use semel_sys;

static INIT: Once = Once::new();

// ── Enums ────────────────────────────────────────────────────────────

/// Manifold type for the ambient space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum Manifold {
    FlatUnbounded = 0,
    SurfaceCylinder = 1,
    SurfaceTorus = 2,
}

/// Distance metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum Distance {
    L2 = 0,
    Cosine = 1,
    SurfaceCylinder = 2,
    SurfaceTorus = 3,
}

/// Filtration type for complex construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Filtration {
    Alpha = 0,
    Rips = 1,
}

// ── Result types ─────────────────────────────────────────────────────

/// Eigendecomposition result from the Hodge Laplacian.
#[derive(Debug, Clone)]
pub struct EigenDecomposition {
    /// Eigenvalues for each dimension 0..n-2.
    pub eigenvalues: Vec<Vec<f64>>,
    /// Eigenvectors for each dimension. `None` for sparse mode (nev > 0).
    /// When present, row-major matrix of shape `(na[i], na[i])`.
    pub eigenvectors: Vec<Option<Vec<f64>>>,
}

/// Entropy measures for a probability distribution.
#[derive(Debug, Clone)]
pub struct EntropyMeasures {
    pub shannon: f64,
    pub collision: f64,
    pub min: f64,
    pub renyi: f64,
}

/// Filtration scales.
#[derive(Debug, Clone)]
pub struct FiltrationScales {
    pub scales: Vec<u64>,
    pub precision: f64,
}

// ── Global init ──────────────────────────────────────────────────────

/// Ensure semel global context is initialised (idempotent).
pub fn ensure_init() {
    INIT.call_once(|| unsafe {
        semel_sys::s_init();
    });
}

// ── SemelSlot ────────────────────────────────────────────────────────

/// RAII handle for a single semel point-cloud / complex slot.
///
/// On drop, the complex and point cloud are unloaded in order.
pub struct SemelSlot {
    k: u64,
    has_complex: bool,
    has_point_cloud: bool,
}

impl SemelSlot {
    /// Allocate a new slot index from semel's internal counter.
    pub fn new() -> Self {
        ensure_init();
        let k = unsafe { semel_sys::_get_point_cloud_index() };
        SemelSlot {
            k,
            has_complex: false,
            has_point_cloud: false,
        }
    }

    /// Wrap an existing point-cloud slot index (e.g. from time-delay embeddings).
    pub fn from_point_cloud_index(k: u64) -> Self {
        ensure_init();
        SemelSlot {
            k,
            has_complex: false,
            has_point_cloud: true,
        }
    }

    /// Return the internal slot index.
    pub fn index(&self) -> u64 {
        self.k
    }

    // ── Point cloud ──────────────────────────────────────────────

    /// Load a point cloud into this slot (flat unbounded, L2 distance).
    ///
    /// `points` is row-major: n points of dimension d, so `points.len() == n * d`.
    pub fn add_point_cloud(&mut self, n: u32, d: u32, points: &mut [f64]) -> Result<()> {
        self.add_point_cloud_with(n, d, points, Manifold::FlatUnbounded, Distance::L2, None)
    }

    /// Load a point cloud with custom manifold and distance metric.
    pub fn add_point_cloud_with(
        &mut self,
        n: u32,
        d: u32,
        points: &mut [f64],
        manifold: Manifold,
        distance: Distance,
        border: Option<&mut [f64]>,
    ) -> Result<()> {
        assert_eq!(points.len(), (n as usize) * (d as usize));
        let (nb, border_ptr) = match border {
            Some(b) => (b.len() as u32, b.as_mut_ptr()),
            None => (0, std::ptr::null_mut()),
        };
        let rc = unsafe {
            semel_sys::_add_point_cloud_raw(
                self.k,
                manifold as u32,
                distance as u32,
                n,
                d,
                points.as_mut_ptr(),
                nb,
                border_ptr,
            )
        };
        if rc != 0 {
            bail!("_add_point_cloud_raw failed (rc={})", rc);
        }
        self.has_point_cloud = true;
        Ok(())
    }

    // ── Complex management ───────────────────────────────────────

    /// Attach a cochain complex over GF(fp) with field descriptor fd.
    pub fn add_complex(&mut self, fp: u32, fd: u32) -> Result<()> {
        let rc = unsafe { semel_sys::_add_complex(self.k, fp, fd) };
        if rc != 0 {
            bail!("_add_complex failed (rc={})", rc);
        }
        self.has_complex = true;
        Ok(())
    }

    /// Attach a cochain complex with l-adic (integral) coefficients.
    ///
    /// `l = 0` gives full integral cohomology; `l > 0` filters by the prime l.
    pub fn add_complex_ladic(&mut self, l: u32) -> Result<()> {
        let rc = unsafe { semel_sys::_add_complex_ladic(self.k, l) };
        if rc != 0 {
            bail!("_add_complex_ladic failed (rc={})", rc);
        }
        self.has_complex = true;
        Ok(())
    }

    // ── Triangulation + filtration ───────────────────────────────

    /// Compute Delaunay triangulation of the loaded point cloud.
    pub fn generate_delaunay(&self) -> Result<()> {
        let rc = unsafe { semel_sys::_generate_delaunay(self.k) };
        if rc != 0 {
            bail!("_generate_delaunay failed (rc={})", rc);
        }
        Ok(())
    }

    /// Build the filtered simplicial complex up to dimension `dim_sk`
    /// with maximum filtration radius `rmax`.
    pub fn build_complex(&self, dim_sk: u32, rmax: f64, filtration: Filtration) -> Result<()> {
        let rc =
            unsafe { semel_sys::_build_complex(self.k, dim_sk, rmax, filtration as u8) };
        if rc != 0 {
            bail!("_build_complex failed (rc={})", rc);
        }
        Ok(())
    }

    /// Add simplices directly (bypassing triangulation).
    ///
    /// `indices` is a flat array of vertex indices, row-major:
    /// `n_simplices * verts_per_simplex` entries.
    pub fn add_simplices(
        &self,
        type_s: u8,
        n_simplices: u32,
        verts_per_simplex: u32,
        indices: &mut [u64],
    ) -> Result<()> {
        assert_eq!(
            indices.len(),
            (n_simplices as usize) * (verts_per_simplex as usize)
        );
        let rc = unsafe {
            semel_sys::_add_simplices_raw(
                self.k,
                type_s,
                n_simplices,
                verts_per_simplex,
                indices.as_mut_ptr(),
            )
        };
        if rc != 0 {
            bail!("_add_simplices_raw failed (rc={})", rc);
        }
        Ok(())
    }

    // ── Cohomology ───────────────────────────────────────────────

    /// Run persistent cohomology and return the result as parsed JSON.
    ///
    /// If `measures_only` is true, only summary measures are returned
    /// (no individual cocycle data).
    pub fn cohomology(&self, measures_only: bool) -> Result<Value> {
        parse_and_free(
            unsafe { semel_sys::_cohomology(self.k, measures_only as u8) },
            "_cohomology",
        )
    }

    // ── Graphs ───────────────────────────────────────────────────

    /// Add a graph (edge list with weights) to this slot.
    ///
    /// `endpoints` is a flat array of vertex pairs: `[v0, v1, v2, v3, ...]`.
    /// `weights` has one entry per edge.
    pub fn add_graph(&mut self, endpoints: &mut [u64], weights: &mut [f64]) -> Result<()> {
        let n_edges = weights.len() as u32;
        assert_eq!(endpoints.len(), 2 * n_edges as usize);
        let rc = unsafe {
            semel_sys::_add_graph_raw(
                self.k,
                n_edges,
                endpoints.as_mut_ptr(),
                weights.as_mut_ptr(),
            )
        };
        if rc != 0 {
            bail!("_add_graph_raw failed (rc={})", rc);
        }
        self.has_point_cloud = true; // graphs use the point-cloud slot
        Ok(())
    }

    // ── Local systems + sheaf ────────────────────────────────────

    /// Set a rank-1 local system (transport weights on edges).
    ///
    /// Must be called after `add_complex`/`add_complex_ladic`, before `build_complex`.
    /// `endpoints` is a flat array of vertex pairs: `[v0, v1, v2, v3, ...]`.
    /// `weights` has one entry per edge.
    pub fn set_local_system(
        &self,
        endpoints: &mut [u64],
        weights: &mut [i64],
    ) -> Result<()> {
        let n_edges = weights.len() as u32;
        assert_eq!(endpoints.len(), 2 * n_edges as usize);
        let rc = unsafe {
            semel_sys::_set_local_system(
                self.k,
                n_edges,
                endpoints.as_mut_ptr(),
                weights.as_mut_ptr(),
            )
        };
        if rc != 0 {
            bail!("_set_local_system failed (rc={})", rc);
        }
        Ok(())
    }

    /// Set rank-r sheaf transport matrices on edges.
    ///
    /// `endpoints` is a flat array of vertex pairs: `[v0, v1, v2, v3, ...]`.
    /// `matrices` is a flat array of `n_edges * rank * rank` int64 values
    /// (row-major transport matrices).
    pub fn set_sheaf(
        &self,
        rank: u32,
        endpoints: &mut [u64],
        matrices: &mut [i64],
    ) -> Result<()> {
        let n_edges = (endpoints.len() / 2) as u32;
        assert_eq!(endpoints.len(), 2 * n_edges as usize);
        assert_eq!(
            matrices.len(),
            (n_edges as usize) * (rank as usize) * (rank as usize)
        );
        let rc = unsafe {
            semel_sys::_set_sheaf(
                self.k,
                rank,
                n_edges,
                endpoints.as_mut_ptr(),
                matrices.as_mut_ptr(),
            )
        };
        if rc != 0 {
            bail!("_set_sheaf failed (rc={})", rc);
        }
        Ok(())
    }

    /// Run sheaf cohomology and return the result as parsed JSON.
    pub fn sheaf_cohomology(&self, measures_only: bool) -> Result<Value> {
        parse_and_free(
            unsafe { semel_sys::_sheaf_cohomology(self.k, measures_only as u8) },
            "_sheaf_cohomology",
        )
    }

    /// Sheaf Laplacian eigendecomposition.
    ///
    /// `n`: number of dimensions.
    /// `nev`: number of eigenvalues to compute (0 = all, with eigenvectors).
    pub fn sheaf_laplace_eigen(&self, n: u32, nev: u32) -> Result<EigenDecomposition> {
        laplace_eigen_impl(self.k, n, nev, true)
    }

    // ── Export ────────────────────────────────────────────────────

    /// Export cocycles as parsed JSON.
    pub fn export_cocycles(&self) -> Result<Value> {
        parse_and_free(
            unsafe { semel_sys::_export_cocycles(self.k) },
            "_export_cocycles",
        )
    }

    /// Export the simplicial complex as parsed JSON.
    pub fn export_complex(&self) -> Result<Value> {
        parse_and_free(
            unsafe { semel_sys::_export_complex(self.k) },
            "_export_complex",
        )
    }

    /// Export the coboundary operator for dimension `d` as parsed JSON.
    pub fn export_coboundary(&self, d: u32) -> Result<Value> {
        parse_and_free(
            unsafe { semel_sys::_export_coboundary(self.k, d) },
            "_export_coboundary",
        )
    }

    // ── Hodge Laplacian ──────────────────────────────────────────

    /// Hodge Laplacian eigendecomposition.
    ///
    /// `n`: number of dimensions.
    /// `nev`: number of eigenvalues to compute (0 = all, with eigenvectors).
    pub fn laplace_eigen(&self, n: u32, nev: u32) -> Result<EigenDecomposition> {
        laplace_eigen_impl(self.k, n, nev, false)
    }

    // ── Filtration scales ────────────────────────────────────────

    /// Get the filtration scales for this slot.
    pub fn get_scales(&self) -> Result<FiltrationScales> {
        let mut nr: u64 = 0;
        let mut r: *mut u64 = std::ptr::null_mut();
        let mut precision: f64 = 0.0;
        let rc = unsafe { semel_sys::_get_scales(self.k, &mut nr, &mut r, &mut precision) };
        if rc != 0 {
            bail!("_get_scales failed (rc={})", rc);
        }
        let scales = if nr > 0 && !r.is_null() {
            unsafe { std::slice::from_raw_parts(r, nr as usize) }.to_vec()
        } else {
            Vec::new()
        };
        Ok(FiltrationScales { scales, precision })
    }
}

impl Drop for SemelSlot {
    fn drop(&mut self) {
        if self.has_complex {
            unsafe { semel_sys::_unload_complex(self.k) };
        }
        if self.has_point_cloud {
            unsafe { semel_sys::_unload_point_cloud(self.k) };
        }
    }
}

impl Default for SemelSlot {
    fn default() -> Self {
        Self::new()
    }
}

// ── TimeSeriesSlot ───────────────────────────────────────────────────

/// RAII handle for a time-series slot.
///
/// On drop, the time series is unloaded.
pub struct TimeSeriesSlot {
    k: u64,
    loaded: bool,
}

impl TimeSeriesSlot {
    /// Create a new time-series slot with the given index.
    pub fn new(k: u64) -> Self {
        ensure_init();
        TimeSeriesSlot { k, loaded: false }
    }

    /// Return the slot index.
    pub fn index(&self) -> u64 {
        self.k
    }

    /// Load a time series into this slot.
    ///
    /// `data` is row-major: n samples of dimension d, so `data.len() == n * d`.
    pub fn add_time_series(&mut self, d: u32, n: u32, data: &mut [f64]) -> Result<()> {
        assert_eq!(data.len(), (n as usize) * (d as usize));
        let rc =
            unsafe { semel_sys::_add_time_series_raw(self.k, d, n, data.as_mut_ptr()) };
        if rc != 0 {
            bail!("_add_time_series_raw failed (rc={})", rc);
        }
        self.loaded = true;
        Ok(())
    }

    /// Compute time-delay embeddings and return a [`SemelSlot`] for the
    /// resulting point cloud.
    ///
    /// `tau`: delay parameter. `d`: embedding dimension.
    pub fn time_delay_embeddings(&self, tau: u32, d: u32) -> Result<SemelSlot> {
        let mut k_p: u64 = 0;
        let rc =
            unsafe { semel_sys::_time_delay_embeddings(self.k, tau, d, &mut k_p) };
        if rc != 0 {
            bail!("_time_delay_embeddings failed (rc={})", rc);
        }
        Ok(SemelSlot::from_point_cloud_index(k_p))
    }
}

impl Drop for TimeSeriesSlot {
    fn drop(&mut self) {
        if self.loaded {
            unsafe { semel_sys::_unload_time_series(self.k) };
        }
    }
}

// ── Free functions ───────────────────────────────────────────────────

/// Compute entropy measures for a probability distribution.
///
/// Returns Shannon, collision, min, and Renyi (of order `alpha`) entropies.
pub fn entropy_dist(p: &mut [f64], alpha: f64) -> Result<EntropyMeasures> {
    ensure_init();
    let mut e = [0.0f64; 4];
    let rc = unsafe {
        semel_sys::_entropy_dist(p.len() as u64, p.as_mut_ptr(), alpha, e.as_mut_ptr())
    };
    if rc != 0 {
        bail!("_entropy_dist failed (rc={})", rc);
    }
    Ok(EntropyMeasures {
        shannon: e[0],
        collision: e[1],
        min: e[2],
        renyi: e[3],
    })
}

/// Compute word-length distribution for a binary sequence.
///
/// `x`: binary sequence (0/1 values). `w`: word length.
/// Returns the distribution as a vector of probabilities.
pub fn dist_words_seq_bin(x: &mut [u8], w: u64) -> Result<Vec<f64>> {
    ensure_init();
    let mut np: u64 = 0;
    let mut p: *mut f64 = std::ptr::null_mut();
    let rc = unsafe {
        semel_sys::_dist_words_seq_bin(x.len() as u64, x.as_mut_ptr(), w, &mut np, &mut p)
    };
    if rc != 0 {
        bail!("_dist_words_seq_bin failed (rc={})", rc);
    }
    let result = if np > 0 && !p.is_null() {
        unsafe { std::slice::from_raw_parts(p, np as usize) }.to_vec()
    } else {
        Vec::new()
    };
    Ok(result)
}

/// Compute the NRPS value for a binary sequence.
pub fn nrps_seq_bin(x: &mut [u8]) -> Result<i32> {
    ensure_init();
    let rc = unsafe { semel_sys::_nrps_seq_bin(x.len() as u64, x.as_mut_ptr()) };
    Ok(rc)
}

/// Compute the NSRPS (Non-Squarefree Radical Permutation Sequence).
///
/// `s` is modified in place.
pub fn nsrps(s: &mut [u32]) -> Result<()> {
    ensure_init();
    let rc = unsafe { semel_sys::_nsrps(s.len() as u64, s.as_mut_ptr()) };
    if rc != 0 {
        bail!("_nsrps failed (rc={})", rc);
    }
    Ok(())
}

/// Generate a logistic map trajectory.
///
/// `lambda`: parameter value. `np`: number of points.
/// Returns the trajectory as a vector.
pub fn ds_logistic(lambda: f64, np: u64) -> Result<Vec<f64>> {
    ensure_init();
    let mut p: *mut f64 = std::ptr::null_mut();
    let rc = unsafe { semel_sys::ds_logistic(lambda, np, &mut p) };
    if rc != 0 {
        bail!("ds_logistic failed (rc={})", rc);
    }
    let result = if np > 0 && !p.is_null() {
        unsafe { std::slice::from_raw_parts(p, np as usize) }.to_vec()
    } else {
        Vec::new()
    };
    Ok(result)
}

/// Unload all point clouds (global cleanup).
pub fn unload_point_cloud_all() -> Result<()> {
    ensure_init();
    let rc = unsafe { semel_sys::_unload_point_cloud_all() };
    if rc != 0 {
        bail!("_unload_point_cloud_all failed (rc={})", rc);
    }
    Ok(())
}

/// Unload all complexes (global cleanup).
pub fn unload_complex_all() -> Result<()> {
    ensure_init();
    let rc = unsafe { semel_sys::_unload_complex_all() };
    if rc != 0 {
        bail!("_unload_complex_all failed (rc={})", rc);
    }
    Ok(())
}

/// Unload all time series (global cleanup).
pub fn unload_time_series_all() -> Result<()> {
    ensure_init();
    let rc = unsafe { semel_sys::_unload_time_series_all() };
    if rc != 0 {
        bail!("_unload_time_series_all failed (rc={})", rc);
    }
    Ok(())
}

// ── Internal helpers ─────────────────────────────────────────────────

/// Shared implementation for `_Laplace_eigen` and `_sheaf_Laplace_eigen`.
fn laplace_eigen_impl(
    k: u64,
    n: u32,
    nev: u32,
    sheaf: bool,
) -> Result<EigenDecomposition> {
    let mut na_ptr: *mut u32 = std::ptr::null_mut();
    let mut a_ptr: *mut *mut f64 = std::ptr::null_mut();
    let mut v_ptr: *mut *mut f64 = std::ptr::null_mut();

    let rc = if sheaf {
        unsafe {
            semel_sys::_sheaf_Laplace_eigen(
                k, n, nev, &mut na_ptr, &mut a_ptr, &mut v_ptr,
            )
        }
    } else {
        unsafe {
            semel_sys::_Laplace_eigen(k, n, nev, &mut na_ptr, &mut a_ptr, &mut v_ptr)
        }
    };

    if rc != 0 {
        let name = if sheaf {
            "_sheaf_Laplace_eigen"
        } else {
            "_Laplace_eigen"
        };
        bail!("{} failed (rc={})", name, rc);
    }

    let ndims = (n - 1) as usize;
    let na_slice = unsafe { std::slice::from_raw_parts(na_ptr, ndims) };

    let mut eigenvalues = Vec::with_capacity(ndims);
    let mut eigenvectors = Vec::with_capacity(ndims);

    for i in 0..ndims {
        let nk = na_slice[i] as usize;
        if nk == 0 {
            eigenvalues.push(Vec::new());
            eigenvectors.push(None);
            continue;
        }

        let a_i = unsafe { *a_ptr.add(i) };
        eigenvalues.push(unsafe { std::slice::from_raw_parts(a_i, nk) }.to_vec());

        if nev == 0 {
            let v_i = unsafe { *v_ptr.add(i) };
            eigenvectors
                .push(Some(unsafe { std::slice::from_raw_parts(v_i, nk * nk) }.to_vec()));
        } else {
            eigenvectors.push(None);
        }
    }

    Ok(EigenDecomposition {
        eigenvalues,
        eigenvectors,
    })
}

/// Parse a malloc'd C JSON string into `serde_json::Value`, then free it.
fn parse_and_free(ptr: *mut std::os::raw::c_char, fn_name: &str) -> Result<Value> {
    if ptr.is_null() {
        bail!("{} returned null", fn_name);
    }
    let c_str = unsafe { CStr::from_ptr(ptr) };
    let json_str = c_str
        .to_str()
        .map_err(|e| anyhow::anyhow!("invalid UTF-8 from {}: {}", fn_name, e))?;
    let value: Value = serde_json::from_str(json_str)?;
    unsafe { libc::free(ptr as *mut libc::c_void) };
    Ok(value)
}
