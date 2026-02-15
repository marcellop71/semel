#include "common.h"
#include "semel.h"

#include <libqhull_r/libqhull_r.h>

semel_ctx_t * ctx = NULL;

int32_t s_init(
    )
{
  if (ctx == NULL)
  {
    ctx = (semel_ctx_t *) malloc(sizeof(semel_ctx_t));
    ctx_log_init(ctx);
    ctx_init(ctx);
  }
  return 0;
}

int32_t s_free(
    )
{
  if (ctx != NULL)
  {
    _unload_point_cloud_all();
    _unload_complex_all();
    _unload_time_series_all();
    ctx_free(ctx);
    ctx_log_free(ctx);
    free(ctx);
    ctx = NULL;
  }
  return 0;
}

// to check
uint64_t _get_point_cloud_index(
    )
{
  uint64_t k = 0;
  uint64_t * p;
  JLL(p, ctx->P, k);
  return (k+1);
}

int32_t _unload_point_cloud(
    uint64_t k)
{
  int rc;
  uint64_t * p;
  JLG(p, ctx->P, k);
  if (p != NULL)
  {
    point_cloud_t * P = (point_cloud_t *) *p;
    point_cloud_free(P);
    free(P);
    JLD(rc, ctx->P, k);
  }
  return 0;
}

int32_t _unload_point_cloud_all(
    )
{
  uint64_t * p;
  uint64_t k;
  k = 0; JLF(p, ctx->P, k);
  while (p != NULL)
  {
    _unload_point_cloud(k);
    JLN(p, ctx->P, k);
  }
  return 0;
}

int32_t _add_point_cloud_raw(
    uint64_t k,
    uint32_t manifold, uint32_t distance,
    uint32_t n, uint32_t d,
    double * points,
    uint32_t nb, double * border)
{
  point_cloud_t * P =
    (point_cloud_t *) malloc(sizeof(point_cloud_t));
  point_cloud_init(P, d, n, manifold, distance, nb, border);

  memcpy(P->p, points, (uint64_t) d * n * sizeof(double));

  uint8_t v [C_1_KB];
  uint64_t * p;
  for (uint32_t i = 0; i < n; ++i)
  {
    uint64_to_256((uint64_t) i, &v[0]);
    JSLI(p, P->pn, v); *p = i;
  }
  P->n_ = n;

  _unload_point_cloud(k);

  uint64_t * pp;
  JLI(pp, ctx->P, k); *pp = (uint64_t) P;

  return 0;
}

int32_t _add_complex(
    uint64_t k, uint32_t fp, uint32_t fd)
{
  complex_cochain_t * Cc =
    (complex_cochain_t *) malloc(sizeof(complex_cochain_t));
  complex_cochain_init(Cc, fp, fd);

  _unload_complex(k);

  uint64_t * p;
  JLI(p, ctx->Cc, k); *p = (uint64_t) Cc;

  return 0;
}

int32_t _add_complex_ladic(
    uint64_t k, uint32_t l)
{
  complex_cochain_t * Cc =
    (complex_cochain_t *) malloc(sizeof(complex_cochain_t));
  complex_cochain_init_z(Cc, l);

  _unload_complex(k);

  uint64_t * p;
  JLI(p, ctx->Cc, k); *p = (uint64_t) Cc;

  return 0;
}

int32_t _unload_complex(
    uint64_t k)
{
  int rc;
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    complex_cochain_free(Cc);
    free(Cc);
    JLD(rc, ctx->Cc, k);
  }
  return 0;
}

int32_t _unload_complex_all(
    )
{
  uint64_t * p;
  uint64_t k;
  k = 0; JLF(p, ctx->Cc, k);
  while (p != NULL)
  {
    _unload_complex(k);
    JLN(p, ctx->Cc, k);
  }
  return 0;
}

int32_t _add_simplices_raw(
    uint64_t k, uint8_t typeS,
    uint32_t n_simplices, uint32_t verts_per_simplex,
    uint64_t * indices)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    complex_simplicial_t * C = &Cc->C;
    C->typeS = typeS;

    uint64_t c;
    JLC(c, C->S, 0, -1);

    uint8_t v [C_1_KB];
    uint64_t * q;

    for (uint32_t i = 0; i < n_simplices; ++i)
    {
      void * s = NULL;
      for (uint32_t j = 0; j < verts_per_simplex; ++j)
      {
        uint64_to_256(indices[(uint64_t) i * verts_per_simplex + j], &v[0]);
        JSLI(q, s, v); *q = 1;
      }
      JLI(q, C->S, i + c); *q = (uint64_t) s;
    }
  }

  return 0;
}

int32_t _generate_delaunay(
    uint64_t k)
{
  uint64_t * p;
  JLG(p, ctx->P, k);
  if (p == NULL)
    return -1;
  point_cloud_t * P = (point_cloud_t *) *p;

  JLG(p, ctx->Cc, k);
  if (p == NULL)
    return -1;
  complex_cochain_t * Cc = (complex_cochain_t *) *p;
  complex_simplicial_t * C = &Cc->C;
  C->typeS = 0;

  // run qhull Delaunay triangulation (reentrant API)
  qhT qh_qh;
  qhT *qh = &qh_qh;
  QHULL_LIB_CHECK
  qh_zero(qh, stderr);

  int dim = (int) P->d;
  int numpoints = (int) P->n;

  // qhull modifies point data in-place, so copy it
  coordT * points = (coordT *) malloc((uint64_t) dim * numpoints * sizeof(coordT));
  memcpy(points, P->p, (uint64_t) dim * numpoints * sizeof(coordT));

  char flags[] = "qhull d Qt Qbb Qc";
  int exitcode = qh_new_qhull(qh, dim, numpoints, points, False, flags, NULL, stderr);

  if (exitcode != 0)
  {
    qh_freeqhull(qh, !qh_ALL);
    int curlong, totlong;
    qh_memfreeshort(qh, &curlong, &totlong);
    free(points);
    return -1;
  }

  // extract Delaunay simplices (skip upper Delaunay facets)
  uint64_t c;
  JLC(c, C->S, 0, -1);

  facetT *facet;
  vertexT *vertex, **vertexp;
  uint64_t si = 0;
  uint8_t v [C_1_KB];
  uint64_t * q;

  FORALLfacets
  {
    if (!facet->upperdelaunay)
    {
      void * s = NULL;
      FOREACHvertex_(facet->vertices)
      {
        uint64_to_256((uint64_t) qh_pointid(qh, vertex->point), &v[0]);
        JSLI(q, s, v); *q = 1;
      }
      JLI(q, C->S, si + c); *q = (uint64_t) s;
      si++;
    }
  }

  qh_freeqhull(qh, !qh_ALL);
  int curlong, totlong;
  qh_memfreeshort(qh, &curlong, &totlong);
  free(points);

  return 0;
}

int32_t _add_graph_raw(
    uint64_t k,
    uint32_t n_edges, uint64_t * endpoints, double * weights)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    complex_simplicial_t * C = &Cc->C;

    uint64_t * q;
    for (uint32_t e = 0; e < n_edges; ++e)
    {
      uint64_t i = endpoints[2 * e];
      uint64_t j = endpoints[2 * e + 1];
      double v = weights[e];
      if (i != j)
      {
        uint64_t ij = (i < j) ? (i << 32) + ((uint32_t) j) : (j << 32) + ((uint32_t) i);
        uint64_t v_;
        memcpy(&v_, &v, sizeof(double));
        JLI(q, C->G, ij); *q = v_;
      }
    }
  }

  return 0;
}

int32_t _set_local_system(
    uint64_t k,
    uint32_t n_edges, uint64_t * endpoints, int64_t * weights)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    int rc;
    JLFA(rc, Cc->local_system);

    uint64_t * q;
    for (uint32_t e = 0; e < n_edges; ++e)
    {
      uint64_t i = endpoints[2 * e];
      uint64_t j = endpoints[2 * e + 1];
      if (i != j)
      {
        uint64_t ij = (i < j) ? (i << 32) + ((uint32_t) j)
                               : (j << 32) + ((uint32_t) i);
        JLI(q, Cc->local_system, ij); *q = (uint64_t) weights[e];
      }
    }
  }
  return 0;
}

int32_t _build_complex(
    uint64_t k, uint32_t dim_sk, double rmax, uint8_t filtration_type)
{
  uint64_t * p, * q;
  JLG(q, ctx->Cc, k);
  if (q != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *q;
    Cc->C.filtration_type = filtration_type;

    point_cloud_t * P = NULL;
    JLG(p, ctx->P, k);
    if (p != NULL)
      { P = (point_cloud_t *) *p; }

    complex_simplicial_build_auto(
      &Cc->C, P, dim_sk, rmax, filtration_type, &ctx->log);

    build_cochain_complex(
      Cc, &ctx->log);
  }

  return 0;
}

char * _export_cocycles(
    uint64_t k)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    return export_cocycles(Cc);
  }
  return NULL;
}

char * _export_complex(
    uint64_t k)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    return export_complex_simplicial(&Cc->C);
  }
  return NULL;
}

char *  _export_coboundary(
    uint64_t k, uint32_t d)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    return export_coboundary_operator(d, Cc);
  }
  return NULL;
}

int32_t _Laplace_eigen(
    uint64_t k, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    if (nev > 0)
      Laplace_eigen_sparse(Cc, n, nev, na, a, v);
    else
    {
      if (Cc->D == NULL)
        coboundary_operator(Cc);
      Laplace_eigen(Cc, n, na, a, v);
    }
  }
  return 0;
}

char * _cohomology(
    uint64_t k, uint8_t measures_only)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    if (Cc->algorithm == COHOMOLOGY_INTEGRAL)
    {
      persistent_cohomology_z(Cc);
      return export_cohomology_z(Cc, measures_only);
    }
    else
    {
      persistent_cohomology(Cc);
      return export_cohomology(Cc, measures_only);
    }
  }
  return NULL;
}

int32_t _get_scales(
    uint64_t k,
    uint64_t * nr, uint64_t ** r, double * precision)
{
  uint64_t * p;
  JLG(p, ctx->P, k);
  if (p != NULL)
  {
    point_cloud_t * P = (point_cloud_t *) *p;
    uint32_t nr32 = 0;
    scales_from_distance_matrix(P->nA, P->A, &nr32, r);
    *nr = (uint64_t) nr32;
    *precision = (double) PRECISION;
  }
  return 0;
}

int32_t _add_time_series_raw(
    uint64_t k, uint32_t d, uint32_t n, double * data)
{
  time_series_t * S =
    (time_series_t *) malloc(sizeof(time_series_t));
  time_series_init(S, d, n);

  memcpy(S->x, data, (uint64_t) d * n * sizeof(double));

  uint64_t * p;
  JLI(p, ctx->S, k); *p = (uint64_t) S;

  return 0;
}

int32_t _unload_time_series(
    uint64_t k)
{
  int rc;
  uint64_t * p;
  JLG(p, ctx->S, k);
  if (p != NULL)
  {
    time_series_free((time_series_t *) *p);
    free((time_series_t *) *p);
    JLD(rc, ctx->S, k);
  }
  return 0;
}

int32_t _unload_time_series_all(
    )
{
  uint64_t * p;
  uint64_t k;
  k = 0; JLF(p, ctx->S, k);
  while (p != NULL)
  {
    _unload_time_series(k);
    JLN(p, ctx->S, k);
  }
  return 0;
}

int32_t _nsrps(
    uint64_t ns, uint32_t * s)
{
  uint64_t t [8];
  t[0] = semel_lap();

  uint64_t c = 0;
  nsrps(ns, s, &c, &ctx->log);

  t[1] = semel_lap();
  //zlog_info(ctx->log.data,
  //  "nsrps: \033[0;36m %.3f\033[0m sec %lu iterations",
  //  ((double) (t[1] - t[0])) / C_1_BLN, c);

  return c;
}

int32_t _entropy_dist(
    uint64_t np, double * p, double alpha, double * e)
{
  return entropy_dist(np, p, alpha, e);
}

int32_t _dist_words_seq_bin(
    uint64_t nx, uint8_t * x, uint64_t w,
    uint64_t * np, double ** p)
{
  return dist_words_seq_bin(nx, x, w, np, p);
}

int32_t _nrps_seq_bin(
    uint64_t nx, uint8_t * x)
{
  return nrps_seq_bin(nx, x);
}

int32_t _set_sheaf(
    uint64_t k, uint32_t rank,
    uint32_t n_edges, uint64_t * endpoints, int64_t * matrices)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;

    // free any existing sheaf transport
    if (Cc->sheaf_transport != NULL)
    {
      int rc;
      uint64_t * sp;
      uint64_t sk = 0;
      JLF(sp, Cc->sheaf_transport, sk);
      while (sp != NULL)
      {
        if (*sp != 0)
          free((int64_t *) *sp);
        JLN(sp, Cc->sheaf_transport, sk);
      }
      JLFA(rc, Cc->sheaf_transport);
    }

    Cc->sheaf_rank = rank;
    uint32_t r2 = rank * rank;

    uint64_t * q;
    for (uint32_t e = 0; e < n_edges; ++e)
    {
      uint64_t i = endpoints[2 * e];
      uint64_t j = endpoints[2 * e + 1];
      if (i != j)
      {
        uint64_t ij = (i < j) ? (i << 32) + ((uint32_t) j)
                               : (j << 32) + ((uint32_t) i);
        int64_t * mat = (int64_t *) malloc(r2 * sizeof(int64_t));
        memcpy(mat, &matrices[(uint64_t) e * r2], r2 * sizeof(int64_t));
        JLI(q, Cc->sheaf_transport, ij); *q = (uint64_t) mat;
      }
    }
  }
  return 0;
}

char * _sheaf_cohomology(
    uint64_t k, uint8_t measures_only)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    if (Cc->sheaf_rank == 0)
      return NULL;

    complex_cochain_t Cc_u;
    if (sheaf_unroll(Cc, &Cc_u) != 0)
      return NULL;

    char * result;
    if (Cc_u.algorithm == COHOMOLOGY_INTEGRAL)
    {
      persistent_cohomology_z(&Cc_u);
      result = export_cohomology_z(&Cc_u, measures_only);
    }
    else
    {
      persistent_cohomology(&Cc_u);
      result = export_cohomology(&Cc_u, measures_only);
    }

    sheaf_unroll_free(&Cc_u);
    return result;
  }
  return NULL;
}

int32_t _sheaf_Laplace_eigen(
    uint64_t k, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v)
{
  uint64_t * p;
  JLG(p, ctx->Cc, k);
  if (p != NULL)
  {
    complex_cochain_t * Cc = (complex_cochain_t *) *p;
    if (Cc->sheaf_rank == 0)
      return -1;

    complex_cochain_t Cc_u;
    if (sheaf_unroll(Cc, &Cc_u) != 0)
      return -1;

    if (nev > 0)
      Laplace_eigen_sparse(&Cc_u, n, nev, na, a, v);
    else
    {
      if (Cc_u.D == NULL)
        coboundary_operator(&Cc_u);
      Laplace_eigen(&Cc_u, n, na, a, v);
    }

    sheaf_unroll_free(&Cc_u);
  }
  return 0;
}
