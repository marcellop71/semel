import json

import numpy as np
from ctypes import *
import ctypes
import gc

from .utils import lib_load, lib_unload, dict_keys_str_to_int, numpy_to_pointer, pointer_to_numpy

class Semel:

    MANIFOLD_FLAT_UNBOUNDED = 0
    MANIFOLD_SURFACE_CYLINDER = 1
    MANIFOLD_SURFACE_TORUS = 2

    DISTANCE_L2 = 0
    DISTANCE_COSINE = 1
    DISTANCE_SURFACE_CYLINDER = 2
    DISTANCE_SURFACE_TORUS = 3

    COHOMOLOGY_FIELD = 0
    COHOMOLOGY_INTEGRAL = 1

    FILTRATION_ALPHA = 0
    FILTRATION_RIPS = 1

    def __init__(self):
        self._semel = lib_load("semel", "0.1.0")

    def __del__(self):
        lib_unload(self._semel)

    # high-level apis

    def cohomology(self, k, measures_only=True):
        q = self._cohomology(k, measures_only)

        q['stats']['mass'] = dict_keys_str_to_int(q['stats']['mass'])
        q['stats']['entropy'] = dict_keys_str_to_int(q['stats']['entropy'])
        q_ = { '1d': q['stats'] }

        if not measures_only:
            q_['cocycles'] = dict_keys_str_to_int(q['cocycles'])
            q_['intervals'] = dict_keys_str_to_int(q['intervals'])

            # integral cohomology extra fields
            if 'betti' in q:
                q_['betti'] = dict_keys_str_to_int(q['betti'])
            if 'torsion' in q:
                q_['torsion'] = dict_keys_str_to_int(q['torsion'])
            if 'l' in q:
                q_['l'] = q['l']

        return q_

    def HodgeLaplace(self, k, d, nev=0, eigenvectors=False):
        [a, v] = self._Laplace_eigen(k, d, nev)
        q = { 'eigenvalues': a, 'eigenvectors': v }
        return q

    # low-level apis wrapping semel C apis

    def add_point_cloud(self, k, manifold, P):
        """Add point cloud via raw pointer."""
        if not isinstance(P, np.ndarray):
            P = np.asarray(P, dtype=np.float64)
        P = np.ascontiguousarray(P, dtype=np.float64)
        n, d = P.shape
        m_type = manifold.get('type', 0)
        m_dist = manifold.get('distance', 0)
        border = manifold.get('border', [])
        nb = len(border)
        if nb > 0:
            b = np.ascontiguousarray(border, dtype=np.float64)
            b_ptr = b.ctypes.data_as(POINTER(c_double))
        else:
            b_ptr = None
        self._semel._add_point_cloud_raw(
            c_uint64(k),
            c_uint(m_type), c_uint(m_dist),
            c_uint(n), c_uint(d),
            P.ctypes.data_as(POINTER(c_double)),
            c_uint(nb), b_ptr)

    # keep add_point_cloud_raw as alias
    add_point_cloud_raw = add_point_cloud

    def unload_point_cloud(self, k):
        self._semel._unload_point_cloud(k)

    def unload_point_cloud_all(self):
        self._semel._unload_point_cloud_all()

    def add_complex(self, k, fp=5, fd=1):
        self._semel._add_complex(k, c_uint(fp), c_uint(fd))

    def add_complex_ladic(self, k, l=0):
        """Create complex for integral/l-adic cohomology.
        l: prime for l-adic filtering (0 = full integral cohomology).
        """
        self._semel._add_complex_ladic(c_uint64(k), c_uint(l))

    def unload_complex(self, k):
        self._semel._unload_complex(k)

    def unload_complex_all(self):
        self._semel._unload_complex_all()

    def add_simplices(self, k, typeS=0, S=None):
        """Add simplices via raw pointer.
        S: list of lists or ndarray of shape (n_simplices, verts_per_simplex).
        """
        if S is None or (hasattr(S, '__len__') and len(S) == 0):
            return
        S = np.ascontiguousarray(S, dtype=np.uint64)
        n_simplices, verts_per_simplex = S.shape
        self._semel._add_simplices_raw(
            c_uint64(k), c_ubyte(typeS),
            c_uint(n_simplices), c_uint(verts_per_simplex),
            S.ctypes.data_as(POINTER(c_uint64)))

    # keep add_simplices_raw as alias
    add_simplices_raw = add_simplices

    def generate_delaunay(self, k):
        """Run Delaunay triangulation in C (via qhull) and add simplices directly."""
        rc = self._semel._generate_delaunay(c_uint64(k))
        if rc != 0:
            raise RuntimeError("_generate_delaunay failed (exit code %d)" % rc)

    def add_graph(self, k, G=None):
        """Add weighted graph edges via raw pointer.
        G: list of [i, j, weight] triples, or ndarray of shape (n_edges, 3).
        """
        if G is None or (hasattr(G, '__len__') and len(G) == 0):
            return
        G = np.asarray(G)
        endpoints = np.ascontiguousarray(G[:, :2], dtype=np.uint64)
        weights = np.ascontiguousarray(G[:, 2], dtype=np.float64)
        n_edges = len(G)
        self._semel._add_graph_raw(
            c_uint64(k),
            c_uint(n_edges),
            endpoints.ctypes.data_as(POINTER(c_uint64)),
            weights.ctypes.data_as(POINTER(c_double)))

    def set_local_system(self, k, edges, weights):
        """Set rank-1 local system (transport weights on edges).
        Must be called after add_complex/add_complex_ladic, before build_complex.
        edges: array of shape (n_edges, 2) - vertex pairs
        weights: array of shape (n_edges,) - integer weights
        """
        edges = np.ascontiguousarray(edges, dtype=np.uint64)
        weights = np.ascontiguousarray(weights, dtype=np.int64)
        n_edges = len(weights)
        self._semel._set_local_system(
            c_uint64(k),
            c_uint(n_edges),
            edges.ctypes.data_as(POINTER(c_uint64)),
            weights.ctypes.data_as(POINTER(c_int64)))

    def build_complex(self, k, dim_max, r_max, filtration=0):
        self._semel._build_complex(k, c_uint(dim_max), c_double(r_max), c_ubyte(filtration))

    def export_cocyles(self, k):
        self._semel._export_cocycles.restype = c_char_p
        q = self._semel._export_cocycles(k)
        return json.loads(str(q.decode("utf-8")))

    def export_complex(self, k):
        self._semel._export_complex.restype = c_char_p
        q = self._semel._export_complex(k)
        return json.loads(str(q.decode("utf-8")))

    def export_coboundary(self, k, d):
        self._semel._export_coboundary.restype = c_char_p
        q = self._semel._export_coboundary(k, c_uint(d))
        return json.loads(str(q.decode("utf-8")))

    def _Laplace_eigen(self, k, n, nev=0):
        pna = pointer(c_uint())
        ppa = pointer(pointer(c_double()))
        ppv = pointer(pointer(c_double()))
        self._semel._Laplace_eigen(
            k, n, c_uint(nev),
            byref(pna), byref(ppa), byref(ppv))
        na = pointer_to_numpy(pna, [n-1,], np.int32)
        a = []; v = []
        for i in range(0, n-1):
            nk = na[i]
            if nk == 0:
                a.append(np.array([])); v.append(np.array([]))
                continue
            atmp = pointer_to_numpy(ppa[i], [nk,], np.float64)
            if nev == 0:
                vtmp = pointer_to_numpy(ppv[i], [nk, nk], np.float64)
            else:
                vtmp = None
            a.append(atmp); v.append(vtmp)
        return [a, v]

    def _cohomology(self, k, measures_only=True):
        measures_only_ = 0 if not measures_only else 1
        self._semel._cohomology.restype = c_char_p
        q = self._semel._cohomology(k, measures_only_)
        q = json.loads(str(q.decode("utf-8")))
        ctypes._reset_cache()
        gc.collect()
        return q

    def add_time_series(self, k, S):
        """Add time series via raw pointer.
        S: list of lists or ndarray of shape (n, d).
        """
        S = np.ascontiguousarray(S, dtype=np.float64)
        n, d = S.shape
        self._semel._add_time_series_raw(
            c_uint64(k), c_uint(d), c_uint(n),
            S.ctypes.data_as(POINTER(c_double)))

    def unload_time_series(self, k):
        self._semel._unload_time_series(k)

    def unload_time_series_all(self):
        self._semel._unload_time_series_all()

    def set_sheaf(self, k, rank, edges, matrices):
        """Set rank-r sheaf transport matrices on edges.
        edges: (n_edges, 2) uint64 — vertex pairs
        matrices: (n_edges, rank, rank) int64 — transport matrices
        """
        edges = np.ascontiguousarray(edges, dtype=np.uint64)
        matrices = np.ascontiguousarray(matrices, dtype=np.int64).reshape(-1)
        n_edges = len(edges)
        self._semel._set_sheaf(
            c_uint64(k), c_uint(rank), c_uint(n_edges),
            edges.ctypes.data_as(POINTER(c_uint64)),
            matrices.ctypes.data_as(POINTER(c_int64)))

    def sheaf_cohomology(self, k, measures_only=True):
        """Persistent sheaf cohomology via unrolling."""
        measures_only_ = 0 if not measures_only else 1
        self._semel._sheaf_cohomology.restype = c_char_p
        q = self._semel._sheaf_cohomology(c_uint64(k), measures_only_)
        if q is None:
            return None
        q = json.loads(str(q.decode("utf-8")))
        ctypes._reset_cache()
        gc.collect()

        q['stats']['mass'] = dict_keys_str_to_int(q['stats']['mass'])
        q['stats']['entropy'] = dict_keys_str_to_int(q['stats']['entropy'])
        q_ = { '1d': q['stats'] }

        if not measures_only:
            q_['cocycles'] = dict_keys_str_to_int(q['cocycles'])
            q_['intervals'] = dict_keys_str_to_int(q['intervals'])

            if 'betti' in q:
                q_['betti'] = dict_keys_str_to_int(q['betti'])
            if 'torsion' in q:
                q_['torsion'] = dict_keys_str_to_int(q['torsion'])
            if 'l' in q:
                q_['l'] = q['l']

        return q_

    def sheaf_Laplacian(self, k, d, nev=0):
        """Sheaf Laplacian eigendecomposition."""
        pna = pointer(c_uint())
        ppa = pointer(pointer(c_double()))
        ppv = pointer(pointer(c_double()))
        self._semel._sheaf_Laplace_eigen(
            c_uint64(k), c_uint(d), c_uint(nev),
            byref(pna), byref(ppa), byref(ppv))
        na = pointer_to_numpy(pna, [d-1,], np.int32)
        a = []; v = []
        for i in range(0, d-1):
            nk = na[i]
            if nk == 0:
                a.append(np.array([])); v.append(np.array([]))
                continue
            atmp = pointer_to_numpy(ppa[i], [nk,], np.float64)
            if nev == 0:
                vtmp = pointer_to_numpy(ppv[i], [nk, nk], np.float64)
            else:
                vtmp = None
            a.append(atmp); v.append(vtmp)
        return { 'eigenvalues': a, 'eigenvectors': v }

    def sliding_window_embedding(self, k, d=2, tau=1):
        self._semel._sliding_window_embedding(k, c_uint(d), c_uint(tau))

    def entropy_dist(self, p, alpha=2.0):
        p = np.ascontiguousarray(p, dtype=np.float64)
        e = np.zeros(4, dtype=np.float64)
        self._semel._entropy_dist(
            c_uint64(len(p)),
            p.ctypes.data_as(POINTER(c_double)),
            c_double(alpha),
            e.ctypes.data_as(POINTER(c_double)))
        return {'shannon': e[0], 'collision': e[1], 'min': e[2], 'renyi': e[3]}

    def dist_words_seq_bin(self, x, w):
        x = np.ascontiguousarray(x, dtype=np.uint8)
        np_ = pointer(c_uint64())
        pp = pointer(pointer(c_double()))
        self._semel._dist_words_seq_bin(
            c_uint64(len(x)),
            x.ctypes.data_as(POINTER(c_ubyte)),
            c_uint64(w),
            byref(np_), byref(pp))
        n = np_.contents.value
        return pointer_to_numpy(pp, [n,], np.float64)

    def nrps_seq_bin(self, x):
        x = np.ascontiguousarray(x, dtype=np.uint8)
        self._semel._nrps_seq_bin(
            c_uint64(len(x)),
            x.ctypes.data_as(POINTER(c_ubyte)))

