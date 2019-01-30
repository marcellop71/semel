#include "common.h"
#include "semel.h"

// ── inverted index helpers ──────────────────────────────────────────
// cob_idx: JudyL( simplex_id -> JudyL( cocycle_id -> 1 ) )
// Tracks which cocycles have a given simplex in their coboundary.

static void cob_idx_add(
    void ** cob_idx, uint64_t simplex_id, uint64_t cocycle_id)
{
  uint64_t * p;
  JLG(p, *cob_idx, simplex_id);
  if (p == NULL)
  {
    void * inner = NULL;
    uint64_t * q;
    JLI(q, inner, cocycle_id); *q = 1;
    JLI(p, *cob_idx, simplex_id); *p = (uint64_t) inner;
  }
  else
  {
    uint64_t * q;
    JLI(q, *p, cocycle_id); *q = 1;
  }
}

static void cob_idx_remove(
    void ** cob_idx, uint64_t simplex_id, uint64_t cocycle_id)
{
  int rc;
  uint64_t * p;
  JLG(p, *cob_idx, simplex_id);
  if (p != NULL)
  {
    JLD(rc, *p, cocycle_id);
    uint64_t cnt; JLC(cnt, *p, 0, -1);
    if (cnt == 0)
      { JLFA(rc, *p); JLD(rc, *cob_idx, simplex_id); }
  }
}

static void cob_idx_add_all(
    void ** cob_idx, void * d_alpha, uint64_t cocycle_id)
{
  uint64_t s, * p;
  s = 0; JLF(p, d_alpha, s);
  while (p != NULL)
    { cob_idx_add(cob_idx, s, cocycle_id); JLN(p, d_alpha, s); }
}

static void cob_idx_remove_all(
    void ** cob_idx, void * d_alpha, uint64_t cocycle_id)
{
  uint64_t s, * p;
  s = 0; JLF(p, d_alpha, s);
  while (p != NULL)
    { cob_idx_remove(cob_idx, s, cocycle_id); JLN(p, d_alpha, s); }
}

// ── persistent cohomology ───────────────────────────────────────────

int32_t persistent_cohomology(
    complex_cochain_t * Cc)
{
  int rc, rc0;
  uint64_t * alpha_i, * alpha_j, * d_alpha_i, * sigma_k, * c_i, * c_j;
  uint64_t * p, * q;
  uint64_t k, i, j, v, ni;
  uint64_t * dim_k, * dim_i, * j_age, * k_age;

  void * V_alpha_dim = NULL;
  void * V_alpha_coboundary = NULL;
  void * cob_idx = NULL;

  fq_t c_i_, c_j_; fq_init(c_i_, Cc->fqctx); fq_init(c_j_, Cc->fqctx);
  fq_t ytmp1, ytmp2; fq_init(ytmp1, Cc->fqctx); fq_init(ytmp2, Cc->fqctx);

  JLFA(rc, Cc->cocycles); JLFA(rc, Cc->intervals);

  // scan simplices in filtration order (k is id)

  k = 0; JLF(sigma_k, Cc->C.simplices, k);
  while (sigma_k != NULL)
  {
    JLG(dim_k, Cc->C.dimension, k);

    // Use inverted index: which cocycles have simplex k in their coboundary?

    void * I = NULL;
    JLG(p, cob_idx, k);
    if (p != NULL)
    {
      uint64_t ci = 0; uint64_t * q_ci;
      JLF(q_ci, *p, ci);
      while (q_ci != NULL)
      {
        JLG(alpha_i, Cc->cocycles, ci);
        if (alpha_i != NULL)
        {
          JLG(dim_i, V_alpha_dim, ci);
          if (dim_i != NULL && (*dim_i) == (*dim_k) - 1)
          {
            JLG(d_alpha_i, V_alpha_coboundary, ci);
            if (d_alpha_i != NULL)
            {
              JLG(q, *d_alpha_i, k);
              if (q != NULL)
              {
                fq_t * c_tmp = (fq_t *) malloc(sizeof(fq_t));
                fq_init(*c_tmp, Cc->fqctx);
                fq_set(*c_tmp, *((fq_t *) *q), Cc->fqctx);
                uint64_t * p_tmp;
                JLI(p_tmp, I, ci); *p_tmp = (uint64_t) c_tmp;
              }
            }
          }
        }

        JLN(q_ci, *p, ci);
      }
    }

    JLC(ni, I, 0, -1);
    if (ni == 0)
    {
      // BIRTH
      void * alpha_k = NULL;
      JLI(p, alpha_k, k); *p = (uint64_t) fq_one_(&Cc->fqctx);
      JLI(p, Cc->cocycles, k); *p = (uint64_t) alpha_k;

      JLI(p, V_alpha_dim, k); *p = (uint64_t) elementary_dimension(alpha_k, Cc);

      void * d_alpha_k = NULL;
      elementary_delta(alpha_k, 1, Cc, &d_alpha_k);
      JLI(p, V_alpha_coboundary, k); *p = (uint64_t) d_alpha_k;

      cob_idx_add_all(&cob_idx, d_alpha_k, k);
    }
    else
    {
      // DEATH: latest-born cocycle j dies
      j = -1; JLL(c_j, I, j);
      fq_set(c_j_, *((fq_t *) *c_j), Cc->fqctx); JLG(alpha_j, Cc->cocycles, j);

      i = 0; JLF(c_i, I, i);
      while ((c_i != NULL) && (i < j))
      { // ── row reduction ──

        fq_set(c_i_, *((fq_t *) *c_i), Cc->fqctx); JLG(alpha_i, Cc->cocycles, i);
        fq_div(ytmp1, c_i_, c_j_, Cc->fqctx); fq_neg(ytmp1, ytmp1, Cc->fqctx);

        // alpha_i += ytmp1 * alpha_j
        v = 0; JLF(p, *alpha_j, v);
        while (p != NULL)
        {
          fq_mul(ytmp2, ytmp1, *((fq_t *) *p), Cc->fqctx);

          fq_t * c = (fq_t *) malloc(sizeof(fq_t));
          fq_init(*c, Cc->fqctx);
          fq_set(*c, ytmp2, Cc->fqctx);

          JLG(q, *alpha_i, v);
          if (q == NULL)
            { JLI(q, *alpha_i, v); *q = (uint64_t) c; }
          else
          {
            fq_t * q_ = (fq_t *) *q;
            fq_add(*q_, *c, *q_, Cc->fqctx);
            fq_clear(*c, Cc->fqctx); free(c);
            if (fq_is_zero(*q_, Cc->fqctx))
              { fq_clear(*q_, Cc->fqctx); free(q_); JLD(rc0, *alpha_i, v); }
          }

          JLN(p, *alpha_j, v);
        }

        // incremental coboundary: d(alpha_i) += ytmp1 * d(alpha_j)
        uint64_t * p_d_j, * p_d_i;
        JLG(p_d_j, V_alpha_coboundary, j);
        JLG(p_d_i, V_alpha_coboundary, i);

        uint64_t s = 0; uint64_t * coeff_j;
        JLF(coeff_j, *p_d_j, s);
        while (coeff_j != NULL)
        {
          fq_mul(ytmp2, ytmp1, *((fq_t *) *coeff_j), Cc->fqctx);

          uint64_t * coeff_i;
          JLG(coeff_i, *p_d_i, s);
          if (coeff_i == NULL)
          {
            fq_t * c_new = (fq_t *) malloc(sizeof(fq_t));
            fq_init(*c_new, Cc->fqctx);
            fq_set(*c_new, ytmp2, Cc->fqctx);
            JLI(coeff_i, *p_d_i, s); *coeff_i = (uint64_t) c_new;
            cob_idx_add(&cob_idx, s, i);
          }
          else
          {
            fq_t * existing = (fq_t *) *coeff_i;
            fq_add(*existing, *existing, ytmp2, Cc->fqctx);
            if (fq_is_zero(*existing, Cc->fqctx))
            {
              fq_clear(*existing, Cc->fqctx); free(existing);
              JLD(rc0, *p_d_i, s);
              cob_idx_remove(&cob_idx, s, i);
            }
          }

          JLN(coeff_j, *p_d_j, s);
        }

        JLN(c_i, I, i);
      } // ── end row reduction ──

      Judy_free_1_b(I, &Cc->fqctx);

      // remove dying cocycle j from inverted index, then free
      JLG(p, V_alpha_coboundary, j);
      if (p != NULL)
      {
        cob_idx_remove_all(&cob_idx, *p, j);
        Judy_free_1_b(*p, &Cc->fqctx);
      }
      JLD(rc0, V_alpha_coboundary, j);

      Judy_free_1_b(*alpha_j, &Cc->fqctx);
      JLD(rc0, Cc->cocycles, j);

      JLG(j_age, Cc->C.age, j); JLG(k_age, Cc->C.age, k);
      if (*k_age > *j_age)
        { JLI(p, Cc->intervals, j); *p = (uint64_t) k; }
    }

    JLN(sigma_k, Cc->C.simplices, k);
  }

  fq_clear(c_i_, Cc->fqctx); fq_clear(c_j_, Cc->fqctx);
  fq_clear(ytmp1, Cc->fqctx); fq_clear(ytmp2, Cc->fqctx);

  Judy_free_2_b(V_alpha_coboundary, &Cc->fqctx);
  JLFA(rc0, V_alpha_dim);
  Judy_free_1_a(cob_idx);

  return 0;
}
