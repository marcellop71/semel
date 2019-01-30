#include "common.h"
#include "semel.h"

int32_t add_faces(
    void * s, complex_simplicial_t * C, point_cloud_t * P,  uint64_t rmax,
    void ** filtration)
{
  int rc;
  uint64_t * p, * q;

  if (simplex_dimension(s) > 0)
  {
    void * faces = NULL;
    simplex_faces(s, &faces);

    uint64_t v = 0; JLF(p, faces, v);
    while (p !=  NULL)
    {
      add_faces(*p, C, P, rmax, filtration);
      JLN(p, faces, v);
    }
    JLFA(rc, faces);
  }

  uint64_t rd = rmax + 1;
  if (C->G != NULL)
    { rd = (uint64_t) (0.50 * simplex_max_pairwise_dist_G(s, C)); }
  else
  {
    if (P != NULL)
    {
      if (C->filtration_type == FILTRATION_ALPHA)
        { rd = (uint64_t) (PRECISION * simplex_circumradius(s, P)); }
      else
        { rd = (uint64_t) (PRECISION * 0.50 * simplex_max_pairwise_dist(s, P)); }
    }
  }

  if ((rmax == (uint64_t) -1) || (rd <= rmax))
  {
    uint8_t sx [MAX_LEN_SIMPLEX_SERIALIZATION];
    simplex_serialize(s, sx);

    JLG(p, *filtration, rd);
    if (p == NULL)
    {
      void * Vtmp = NULL;
      JSLI(q, Vtmp, sx); *q = (uint64_t) s;
      JLI(q, *filtration, rd); *q = (uint64_t) Vtmp;
    }
    else
    {
      JSLG(q, *p, sx);
      if (q == NULL)
        { JSLI(q, *p, sx); *q = (uint64_t) s; }
      else
        { JSLFA(rc, s); }
    }
  }
  else
    { JSLFA(rc, s); }

  return 0;
}

int32_t complex_simplicial_topdown(
    complex_simplicial_t * C, point_cloud_t * P, uint64_t rmax)
{
  uint64_t * p;
  uint64_t v;

  v = 0; JLF(p, C->S, v);
  while (p != NULL)
  {
    add_faces(*p, C, P, rmax, &C->filtration);
    JLN(p, C->S, v);
  }

  return 0;
}

int32_t add_cofaces(
    void * s, void * N, uint64_t rmax, void ** L, point_cloud_t * P,
    complex_simplicial_t * C)
{
  int32_t rc;
  uint64_t * p, * q, * p_, * q_;
  uint64_t v, w, rd;
  uint8_t vv [C_1_KB];
  uint8_t sx [MAX_LEN_SIMPLEX_SERIALIZATION];

  if (C->filtration_type == FILTRATION_ALPHA)
    { rd = (uint64_t) (PRECISION * simplex_circumradius(s, P)); }
  else
    { rd = (uint64_t) (PRECISION * 0.50 * simplex_max_pairwise_dist(s, P)); }

  if ((rmax == (uint64_t) -1) || (rd <= rmax))
  {
    simplex_serialize(s, &sx[0]);

    JLG(p, C->filtration, rd);
    if (p == NULL)
    {
      void * Vtmp_ = NULL;
      JSLI(q, Vtmp_, sx); *q = (uint64_t) s;
      JLI(q, C->filtration, rd); *q = (uint64_t) Vtmp_;
    }
    else
    {
      JSLG(q, *p, sx);
      if (q == NULL)
        { JSLI(q, *p, sx); *q = (uint64_t) s; }
    }

    if (simplex_dimension(s) < C->dim_sk)
    {
      // iterate over {v} in the lower neighborhood N
      v = 0; JLF(p, N, v);
      while (p != NULL)
      {
        // s_ = s + {v}
        void * s_ = NULL;
        vv[0] = 0; JSLF(p_, s, vv);
        while (p_ != NULL)
        {
          JSLI(q_, s_, vv); *q_ = *p_;
          JSLN(p_, s, vv);
        }
        uint64_to_256(v, &vv[0]);
        JSLI(q_, s_, vv); *q_ = 1;

        // intersection of lower neighborhoods
        void * N_ = NULL;
        w = 0; JLF(q, N, w);
        while (q != NULL)
        {
          JLG(p_, L[v], w);
          if (p_ != NULL)
            { JLI(q_, N_, w); *q_ = *q; }
          JLN(q, N, w);
        }

        add_cofaces(
          s_, N_, rmax, L, P, C);

        JLN(p, N, v);
      }
    }
  }
  else
    { JSLFA(rc, s); JLFA(rc, N); }

  return 0;
}

int32_t expand_neighborhood_graph(
    uint32_t np, uint64_t rmax, void ** L, point_cloud_t * P,
    complex_simplicial_t * C)
{
  uint64_t * p, * q;
  uint64_t w = 0;
  uint64_t rd;
  uint8_t v [C_1_KB];

  for(uint32_t u = 0; u < np; ++u)
  {
    // start with {u}
    void * s = NULL;
    uint64_to_256((uint64_t) u, &v[0]);
    JSLI(p, s, v); *p = 1;

    // get lower neighborhood N of {u} for radius r
    void * N = NULL;
    w = 0; JLF(p, L[u], w);
    while (p != NULL)
    {
      rd = *p;
      if (rd <= rmax)
        { JLI(q, N, w); *q = *p; }
      JLN(p, L[u], w);
    }

    // add cofaces from this lower neighborhood
    add_cofaces(
      s, N, rmax, L, P, C);
  }

  return 0;
}

int32_t annotated_lower_neighborhoods(
    uint64_t rmax, uint32_t nA, double * A,
    void *** L)
{
  uint64_t * q;
  uint32_t i, j;
  uint64_t rd;

  i = 1; j = 0;
  for(uint32_t k = 0; k < nA; ++k)
  {
    rd = (uint64_t) (PRECISION * 0.50 * A[k]);
    if ((rmax == (uint64_t) -1) || (rd <= rmax))
      { JLI(q, (*L)[i], (uint64_t) j); *q = rd; }

    ++j;
    if (j == i)
      { j = 0; ++i; }
  }

  return 0;
}

int32_t get_array_of_indices(
    void * a,
    uint32_t * nidx, uint64_t ** idx)
{
  uint64_t * p;
  uint64_t nr, v, k;
  JLC(nr, a, 0, -1);

  *nidx = (uint32_t) nr;
  *idx = (uint64_t *) malloc((*nidx) * sizeof(uint64_t));

  v = 0; JLF(p, a, v); k = 0;
  while (p !=  NULL)
  {
    (*idx)[k] = v;
    JLN(p, a, v); k++;
  }

  return 0;
}

int32_t scales_from_distance_matrix(
    uint32_t nA, double * A,
    uint32_t * nr, uint64_t ** r)
{
  int rc;
  void * rr = NULL;
  uint64_t * q;
  uint64_t rtmp;

  for(uint32_t k = 0; k < nA; ++k)
  {
    rtmp = (uint64_t) (PRECISION * 0.50 * A[k]);
    JLI(q, rr, (uint64_t) rtmp); *q = 0;
  }

  get_array_of_indices(rr, nr, r);
  JLFA(rc, rr);

  return 0;
}

int32_t complex_simplicial_bottomup(
    complex_simplicial_t * C, point_cloud_t * P, uint64_t rmax)
{
  int rc;
  point_cloud_init_distance_matrix(P);

  void ** L = (void **) malloc(P->n * sizeof(void *));
  for(uint32_t k = 0; k < P->n; ++k)
    { L[k] = NULL; }
  annotated_lower_neighborhoods(
    rmax, P->nA, P->A, &L);

  expand_neighborhood_graph(
    P->n, rmax, L, P, C);

  for(uint32_t k = 0; k < P->n; ++k)
    { JLFA(rc, L[k]); }
  free(L);

  return 0;
}

// gives ids to the simplices in the filtration
// order by: radius, lexicographic in the serialization
int32_t complex_simplicial_build(
    complex_simplicial_t * C)
{
  uint64_t * p, * q, * p_;
  uint64_t r, ir, id;
  uint8_t sx [MAX_LEN_SIMPLEX_SERIALIZATION];
  uint32_t dim;

  //get_array_of_indices(
  //  C->filtration,
  //  &C->nr, &C->r);

  r = 0; JLF(p, C->filtration, r); ir = 0; id = 0;
  while (p != NULL)
  {
    sx[0] = 0; JSLF(q, *p, sx);
    while (q != NULL)
    {
      dim = simplex_dimension(*q);

      JLI(p_, C->simplices, id); *p_ = *q;
      JLI(p_, C->dimension, id); *p_ = dim;
      JLI(p_, C->age, id); *p_ = ir;
      JLI(p_, C->radius, id); *p_ = r;
      JSLI(p_, C->sx2id, sx); *p_ = id;
      JLI(p_, C->scd[dim], C->scc[dim]); *p_ = id;
      JLI(p_, C->sci[dim], id); *p_ = C->scc[dim];

      JSLN(q, *p, sx); id++; C->scc[dim]++; C->scc_++;
    }

    JLN(p, C->filtration, r); ir++;
  }

  simplices_faces(C);
  simplices_cofaces(C);

  return 0;
}

int32_t complex_simplicial_build_auto(
    complex_simplicial_t * C,
    point_cloud_t * P,
    uint32_t dim_sk, double rmax,
    uint8_t filtration_type,
    semel_ctx_log_t * log)
{
  C->filtration_type = filtration_type;
  C->dim_sk = dim_sk;
  uint64_t rmax_ = (uint64_t) (PRECISION * rmax);

  if (C->S == NULL)
    {complex_simplicial_bottomup(C, P, rmax_);}
  else
  {
    if (P != NULL)
      { point_cloud_init_distance_matrix(P); }
    complex_simplicial_topdown(C, P, rmax_);
  }

  complex_simplicial_build(C);

  return 0;
}

int32_t complex_simplicial_init_array(
    uint32_t nC,
    complex_simplicial_t ** C)
{
  *C = (complex_simplicial_t *) malloc(nC * sizeof(complex_simplicial_t));
  for(uint32_t k = 0; k < nC; ++k)
    { complex_simplicial_init(&(*C)[k]); }

  return 0;
}
