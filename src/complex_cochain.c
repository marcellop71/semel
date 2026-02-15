#include "common.h"
#include "semel.h"

uint32_t elementary_dimension(
    void * alpha, complex_cochain_t * Cc)
{
  uint32_t d = 0;

  uint32_t dtmp;
  uint64_t * p, * q;
  uint64_t k;

  k = 0; JLF(p, alpha, k);
  while (p != NULL)
  {
    JLG(q, Cc->C.dimension, k); dtmp = *q;
    if (dtmp > d)
      { d = dtmp; }
    JLN(p, alpha, k);
  }

  return d;
}

fq_t * sign_parity(
    uint64_t k, fq_ctx_t fqctx)
{
  fq_t * c = (fq_t *) malloc(sizeof(fq_t));
  fq_init(*c, fqctx);

  if ((k % 2) == 0)
    { fq_set_si(*c, 1, fqctx); }
  else
    { fq_set_si(*c, -1, fqctx);}

  return c;
}

fq_t * sign_in_coface(
    uint64_t v0, uint64_t v1, fq_ctx_t fqctx, complex_simplicial_t * C)
{
  uint64_t * v0_s, * v1_s;

  fq_t * c = (fq_t *) malloc(sizeof(fq_t));
  fq_init(*c, fqctx);

  JLG(v0_s, C->simplices, v0); JLG(v1_s, C->simplices, v1);
  if ((v0_s != NULL) && (v1_s != NULL))
    { fq_set_si(*c, simplex_sign_in_coface(*v0_s, *v1_s), fqctx); }

  return c;
}

int32_t elementary_simplices_boundary(
    complex_cochain_t * Cc)
{
  uint64_t v, k;
  uint64_t * v_faces, * p, * q;

  v = 0; JLF(v_faces, Cc->C.faces, v);
  while (v_faces != NULL)
  {
    void * v_boundary = NULL;

    k = 0; JLF(p, *v_faces, k);
    while (p != NULL)
    {
      JLI(q, v_boundary, k); *q = (uint64_t) sign_parity(k, Cc->fqctx);
      JLN(p, *v_faces, k);
    }

    JLI(q, Cc->elementary_boundaries, v); *q = (uint64_t) v_boundary;
    JLN(v_faces, Cc->C.faces, v);
  }

  return 0;
}

static int64_t get_local_transport(
    uint64_t vert_0, uint64_t vert_1, void * local_system)
{
  if (local_system == NULL)
    return 1;
  uint64_t edge_key = (vert_0 < vert_1)
    ? (vert_0 << 32) + ((uint32_t) vert_1)
    : (vert_1 << 32) + ((uint32_t) vert_0);
  uint64_t * p;
  JLG(p, local_system, edge_key);
  if (p != NULL)
    return (int64_t) *p;
  return 1;
}

int32_t elementary_simplices_coboundary(
    complex_cochain_t * Cc)
{
  uint64_t v, k;
  uint64_t * v_cofaces, * p, * q;

  v = 0; JLF(v_cofaces, Cc->C.cofaces, v);
  while (v_cofaces != NULL)
  {
    void * v_coboundary = NULL;

    k = 0; JLF(p, *v_cofaces, k);
    while (p != NULL)
    {
      fq_t * c = sign_in_coface(v, k, Cc->fqctx, &Cc->C);

      if (Cc->local_system != NULL)
      {
        uint64_t * v_s, * k_s;
        JLG(v_s, Cc->C.simplices, v);
        JLG(k_s, Cc->C.simplices, k);
        if ((v_s != NULL) && (k_s != NULL))
        {
          uint64_t pos = simplex_index_in_coface(*v_s, *k_s);
          if (pos == 0)
          {
            uint8_t vv [C_1_KB];
            uint64_t * pv;
            vv[0] = 0; JSLF(pv, *k_s, vv);
            uint64_t vert_0 = uint256_to_64(&vv[0]);
            JSLN(pv, *k_s, vv);
            uint64_t vert_1 = uint256_to_64(&vv[0]);
            int64_t w = get_local_transport(vert_0, vert_1, Cc->local_system);
            if (w != 1)
            {
              fq_t wt;
              fq_init(wt, Cc->fqctx);
              fq_set_si(wt, w, Cc->fqctx);
              fq_mul(*c, *c, wt, Cc->fqctx);
              fq_clear(wt, Cc->fqctx);
            }
          }
        }
      }

      JLI(q, v_coboundary, k); *q = (uint64_t) c;
      JLN(p, *v_cofaces, k);
    }

    JLI(q, Cc->elementary_coboundaries, v); *q = (uint64_t) v_coboundary;
    JLN(v_cofaces, Cc->C.cofaces, v);
  }

  return 0;
}

// ── integer coboundary (Z coefficients) ──────────────────────────────

static fmpz_t * sign_in_coface_z(
    uint64_t v0, uint64_t v1, complex_simplicial_t * C)
{
  uint64_t * v0_s, * v1_s;

  fmpz_t * c = (fmpz_t *) malloc(sizeof(fmpz_t));
  fmpz_init(*c);

  JLG(v0_s, C->simplices, v0); JLG(v1_s, C->simplices, v1);
  if ((v0_s != NULL) && (v1_s != NULL))
    { fmpz_set_si(*c, simplex_sign_in_coface(*v0_s, *v1_s)); }

  return c;
}

int32_t elementary_simplices_coboundary_z(
    complex_cochain_t * Cc)
{
  uint64_t v, k;
  uint64_t * v_cofaces, * p, * q;

  v = 0; JLF(v_cofaces, Cc->C.cofaces, v);
  while (v_cofaces != NULL)
  {
    void * v_coboundary = NULL;

    k = 0; JLF(p, *v_cofaces, k);
    while (p != NULL)
    {
      fmpz_t * c = sign_in_coface_z(v, k, &Cc->C);

      if (Cc->local_system != NULL)
      {
        uint64_t * v_s, * k_s;
        JLG(v_s, Cc->C.simplices, v);
        JLG(k_s, Cc->C.simplices, k);
        if ((v_s != NULL) && (k_s != NULL))
        {
          uint64_t pos = simplex_index_in_coface(*v_s, *k_s);
          if (pos == 0)
          {
            uint8_t vv [C_1_KB];
            uint64_t * pv;
            vv[0] = 0; JSLF(pv, *k_s, vv);
            uint64_t vert_0 = uint256_to_64(&vv[0]);
            JSLN(pv, *k_s, vv);
            uint64_t vert_1 = uint256_to_64(&vv[0]);
            int64_t w = get_local_transport(vert_0, vert_1, Cc->local_system);
            if (w != 1)
            {
              fmpz_t wt;
              fmpz_init(wt);
              fmpz_set_si(wt, w);
              fmpz_mul(*c, *c, wt);
              fmpz_clear(wt);
            }
          }
        }
      }

      JLI(q, v_coboundary, k); *q = (uint64_t) c;
      JLN(p, *v_cofaces, k);
    }

    JLI(q, Cc->elementary_coboundaries_z, v); *q = (uint64_t) v_coboundary;
    JLN(v_cofaces, Cc->C.cofaces, v);
  }

  return 0;
}

uint32_t elementary_delta_z(
    void * alpha, complex_cochain_t * Cc,
    void ** d_alpha)
{
  *d_alpha = NULL;

  void * e = Cc->elementary_coboundaries_z;

  int rc;
  uint64_t * p, * q, * p_, * d_sigma_k;
  uint64_t k, j;

  k = 0; JLF(p, alpha, k);
  while (p != NULL)
  {
    JLG(d_sigma_k, e, k);
    if (d_sigma_k != NULL)
    {
      j = 0; JLF(q, *d_sigma_k, j);
      while(q != NULL)
      {
        fmpz_t * c = (fmpz_t *) malloc(sizeof(fmpz_t));
        fmpz_init(*c);
        fmpz_set(*c, *((fmpz_t *) *q));
        fmpz_mul(*c, *c, *((fmpz_t *) *p));

        JLG(p_, *d_alpha, j);
        if (p_ == NULL)
          { JLI(p_, *d_alpha, j); *p_ = (uint64_t) c; }
        else
        {
          fmpz_t * p__ = (fmpz_t *) *p_;
          fmpz_add(*p__, *c, *p__);
          fmpz_clear(*c); free(c);
          if (fmpz_is_zero(*p__))
            { fmpz_clear(*p__); free(p__); JLD(rc, *d_alpha, j); }
        }

        JLN(q, *d_sigma_k, j);
      }
    }

    JLN(p, alpha, k);
  }

  return 0;
}

// ── field coboundary (F_p coefficients) ──────────────────────────────

uint32_t elementary_delta(
    void * alpha, uint32_t method, complex_cochain_t * Cc,
    void ** d_alpha)
{
  *d_alpha = NULL;

  void * e = (method == 0) ?
    Cc->elementary_boundaries : Cc->elementary_coboundaries;

  int rc;
  uint64_t * p, * q, * p_, * d_sigma_k;
  uint64_t k, j;

  // alpha = ... + (*p) [sigma_k] + ...

  k = 0; JLF(p, alpha, k);
  while (p != NULL)
  {
    JLG(d_sigma_k, e, k);
    if (d_sigma_k != NULL)
    {
      // d_sigma_k = ... + (*q) [sigma_j] + ...
      j = 0; JLF(q, *d_sigma_k, j);
      while(q != NULL)
      {
        fq_t * c = (fq_t *) malloc(sizeof(fq_t));
        fq_init(*c, Cc->fqctx);
        fq_set(*c, *((fq_t *) *q), Cc->fqctx);
        fq_mul(*c, *c, *((fq_t *) *p), Cc->fqctx);

        JLG(p_, *d_alpha, j);
        if (p_ == NULL)
          { JLI(p_, *d_alpha, j); *p_ = (uint64_t) c; }
        else
        {
          fq_t * p__ = (fq_t *) *p_;
          fq_add(*p__, *c, *p__, Cc->fqctx);
          fq_clear(*c, Cc->fqctx); free(c);
          if (fq_is_zero(*p__, Cc->fqctx))
            { fq_clear(*p__, Cc->fqctx); free(p__); JLD(rc, *d_alpha, j); }
        }

        JLN(q, *d_sigma_k, j);
      }
    }

    JLN(p, alpha, k);
  }

  return 0;
}

int32_t coboundary_operator(
    complex_cochain_t * Cc)
{
  uint64_t v, w;
  uint64_t * p, * q, * p_, * q_, * p__, * q__;

  for(uint32_t k = 0; k < Cc->C.dim_sk; k++)
  {
    void * Dk = NULL;

    v = 0; JLF(p, Cc->C.scd[k], v);
    while (p != NULL)
    {
      JLG(q, Cc->elementary_coboundaries, *p);
      if (q != NULL)
      {
        void * Dk_v = NULL;

        w = 0; JLF(p_, *q, w);
        while (p_ != NULL)
        {
          JLG(q_, Cc->C.dimension, w);
          if (q_ != NULL)
          {
            JLG(p__, Cc->C.sci[*q_], w);
            if (p__ != NULL)
            {
              char * tmp = fq_get_str_pretty(*((fq_t *) *p_), Cc->fqctx);
              int32_t a = atoi(tmp);
              flint_free(tmp);
              fmpz_t f; fq_ctx_order(f, Cc->fqctx);
              uint32_t o = fmpz_get_ui(f);
              fmpz_clear(f);
              a = (a > (int32_t) (o >> 1)) ? (a - o) : a;

              JLI(q__, Dk_v, *p__); *q__ = a;
            }
          }

          JLN(p_, *q, w);
        }

        JLI(p_, Dk, v); *p_ = (uint64_t) Dk_v;
      }

      JLN(p, Cc->C.scd[k], v);
    }

    JLI(p, Cc->D, k); *p = (uint64_t) Dk;
  }

  return 0;
}

int32_t coboundary_operator_as_gsl_matrix(
    uint32_t k, complex_cochain_t * Cc,
    uint32_t * mD, uint32_t * nD, gsl_matrix ** D)
{
  uint64_t v, w;
  uint64_t * p, * q, * p_;

  *mD = Cc->C.scc[k+1]; *nD = Cc->C.scc[k];
  *D = gsl_matrix_alloc(*mD, *nD);

  gsl_matrix_set_zero(*D);

  JLG(p, Cc->D, k);
  v = 0; JLF(q, *p, v);
  while (q != NULL)
  {
    w = 0; JLF(p_, *q, w);
    while (p_ != NULL)
    {
      gsl_matrix_set(
        *D, w, v, (double) ((int32_t) (*p_)));
      JLN(p_, *q, w);
    }
    JLN(q, *p, v);
  }

  return 0;
}

int32_t build_cochain_complex(
    complex_cochain_t * Cc,
    semel_ctx_log_t * log)
{
  if (Cc->algorithm == COHOMOLOGY_INTEGRAL)
    { elementary_simplices_coboundary_z(Cc); }
  else
    { elementary_simplices_coboundary(Cc); }

  return 0;
}
