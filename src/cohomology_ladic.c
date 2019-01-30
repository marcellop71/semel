#include "common.h"
#include "semel.h"

// ── inverted index helpers ──────────────────────────────────────────
// cob_idx: JudyL( simplex_id -> JudyL( cocycle_id -> 1 ) )

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

// ── linear combination over Z ───────────────────────────────────────
// target ← a * target + b * source  (applied to JL of fmpz_t *)
// Entries that become zero are removed.

static void lincomb_z(
    void ** target, fmpz_t a, fmpz_t b, void * source)
{
  int rc;
  uint64_t v, * p, * q;

  // first: scale existing target entries by a
  v = 0; JLF(p, *target, v);
  while (p != NULL)
  {
    fmpz_t * pv = (fmpz_t *) *p;
    fmpz_mul(*pv, *pv, a);
    JLN(p, *target, v);
  }

  // second: add b * source entries
  v = 0; JLF(p, source, v);
  while (p != NULL)
  {
    fmpz_t tmp; fmpz_init(tmp);
    fmpz_mul(tmp, b, *((fmpz_t *) *p));

    JLG(q, *target, v);
    if (q == NULL)
    {
      fmpz_t * c = (fmpz_t *) malloc(sizeof(fmpz_t));
      fmpz_init(*c);
      fmpz_set(*c, tmp);
      JLI(q, *target, v); *q = (uint64_t) c;
    }
    else
    {
      fmpz_t * existing = (fmpz_t *) *q;
      fmpz_add(*existing, *existing, tmp);
      if (fmpz_is_zero(*existing))
        { fmpz_clear(*existing); free(existing); JLD(rc, *target, v); }
    }

    fmpz_clear(tmp);
    JLN(p, source, v);
  }

  // third: clean up any target entries that are now zero (from scaling by a)
  v = 0; JLF(p, *target, v);
  while (p != NULL)
  {
    uint64_t next_v = v;
    JLN(q, *target, next_v);  // peek next before possible delete

    fmpz_t * pv = (fmpz_t *) *p;
    if (fmpz_is_zero(*pv))
      { fmpz_clear(*pv); free(pv); JLD(rc, *target, v); }

    v = next_v;
    p = q;
  }
}

// ── content reduction ────────────────────────────────────────────────
// Divide all coefficients by their GCD.  Does not change the
// cohomology class.  Keeps coefficient bit-lengths bounded.

static void normalize_z(void ** ja)
{
  uint64_t v, * p;
  fmpz_t g; fmpz_init(g);

  // compute GCD of all entries
  v = 0; JLF(p, *ja, v);
  while (p != NULL)
  {
    fmpz_gcd(g, g, *((fmpz_t *) *p));
    if (fmpz_is_one(g)) { fmpz_clear(g); return; }   // early exit
    JLN(p, *ja, v);
  }

  if (fmpz_is_zero(g) || fmpz_is_one(g))
    { fmpz_clear(g); return; }

  // divide all entries by g
  v = 0; JLF(p, *ja, v);
  while (p != NULL)
  {
    fmpz_divexact(*((fmpz_t *) *p), *((fmpz_t *) *p), g);
    JLN(p, *ja, v);
  }

  fmpz_clear(g);
}

// ── persistent cohomology over Z ────────────────────────────────────

int32_t persistent_cohomology_z(
    complex_cochain_t * Cc)
{
  int rc, rc0;
  uint64_t * alpha_i, * alpha_j, * d_alpha_i, * sigma_k, * c_i, * c_j;
  uint64_t * p, * q;
  uint64_t k, i, j, ni;
  uint64_t * dim_k, * dim_i, * j_age, * k_age;

  void * V_alpha_dim = NULL;
  void * V_alpha_coboundary = NULL;
  void * cob_idx = NULL;

  JLFA(rc, Cc->cocycles_z); JLFA(rc, Cc->intervals); JLFA(rc, Cc->intervals_torsion);

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
        JLG(alpha_i, Cc->cocycles_z, ci);
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
                fmpz_t * c_tmp = (fmpz_t *) malloc(sizeof(fmpz_t));
                fmpz_init(*c_tmp);
                fmpz_set(*c_tmp, *((fmpz_t *) *q));
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
      JLI(p, alpha_k, k); *p = (uint64_t) fmpz_one_();
      JLI(p, Cc->cocycles_z, k); *p = (uint64_t) alpha_k;

      JLI(p, V_alpha_dim, k); *p = (uint64_t) elementary_dimension(alpha_k, Cc);

      void * d_alpha_k = NULL;
      elementary_delta_z(alpha_k, Cc, &d_alpha_k);
      JLI(p, V_alpha_coboundary, k); *p = (uint64_t) d_alpha_k;

      cob_idx_add_all(&cob_idx, d_alpha_k, k);
    }
    else
    {
      // DEATH: latest-born cocycle j dies
      j = -1; JLL(c_j, I, j);
      fmpz_t c_j_val; fmpz_init(c_j_val);
      fmpz_set(c_j_val, *((fmpz_t *) *c_j));
      JLG(alpha_j, Cc->cocycles_z, j);

      i = 0; JLF(c_i, I, i);
      while ((c_i != NULL) && (i < j))
      { // ── GCD-based row reduction ──

        fmpz_t c_i_val; fmpz_init(c_i_val);
        fmpz_set(c_i_val, *((fmpz_t *) *c_i));
        JLG(alpha_i, Cc->cocycles_z, i);

        fmpz_t g, cj_over_g, ci_over_g;
        fmpz_init(g); fmpz_init(cj_over_g); fmpz_init(ci_over_g);
        fmpz_gcd(g, c_i_val, c_j_val);
        fmpz_divexact(cj_over_g, c_j_val, g);
        fmpz_divexact(ci_over_g, c_i_val, g);
        fmpz_neg(ci_over_g, ci_over_g);

        // α_i ← (c_j/g)*α_i + (-c_i/g)*α_j
        lincomb_z(alpha_i, cj_over_g, ci_over_g, *alpha_j);

        // δ(α_i) ← (c_j/g)*δ(α_i) + (-c_i/g)*δ(α_j)
        uint64_t * p_d_j, * p_d_i;
        JLG(p_d_j, V_alpha_coboundary, j);
        JLG(p_d_i, V_alpha_coboundary, i);

        // remove i from cob_idx before updating
        cob_idx_remove_all(&cob_idx, *p_d_i, i);

        lincomb_z(p_d_i, cj_over_g, ci_over_g, *p_d_j);

        // re-add i to cob_idx with updated coboundary
        cob_idx_add_all(&cob_idx, *p_d_i, i);

        // content reduction: keep coefficients small
        normalize_z(alpha_i);
        normalize_z(p_d_i);

        fmpz_clear(g); fmpz_clear(cj_over_g); fmpz_clear(ci_over_g);
        fmpz_clear(c_i_val);

        JLN(c_i, I, i);
      } // ── end row reduction ──

      Judy_free_1_z(I);

      // store torsion order: |c_j_val| is the pivot value at death
      fmpz_t abs_pivot; fmpz_init(abs_pivot);
      fmpz_abs(abs_pivot, c_j_val);
      uint64_t torsion = fmpz_get_ui(abs_pivot);
      if (torsion == 1) { torsion = 0; } // free: no torsion
      fmpz_clear(abs_pivot);

      // remove dying cocycle j from inverted index, then free
      JLG(p, V_alpha_coboundary, j);
      if (p != NULL)
      {
        cob_idx_remove_all(&cob_idx, *p, j);
        Judy_free_1_z(*p);
      }
      JLD(rc0, V_alpha_coboundary, j);

      Judy_free_1_z(*alpha_j);
      JLD(rc0, Cc->cocycles_z, j);

      JLG(j_age, Cc->C.age, j); JLG(k_age, Cc->C.age, k);
      if (*k_age > *j_age)
      {
        JLI(p, Cc->intervals, j); *p = (uint64_t) k;
        JLI(p, Cc->intervals_torsion, j); *p = torsion;
      }

      fmpz_clear(c_j_val);
    }

    JLN(sigma_k, Cc->C.simplices, k);
  }

  Judy_free_2_z(V_alpha_coboundary);
  JLFA(rc0, V_alpha_dim);
  Judy_free_1_a(cob_idx);

  return 0;
}
