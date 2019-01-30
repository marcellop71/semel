#include "common.h"
#include "semel.h"

// ── sheaf transport lookup ─────────────────────────────────────────

static int64_t * get_sheaf_transport(
    uint64_t vert_0, uint64_t vert_1, void * sheaf_transport)
{
  if (sheaf_transport == NULL)
    return NULL;
  uint64_t edge_key = (vert_0 < vert_1)
    ? (vert_0 << 32) + ((uint32_t) vert_1)
    : (vert_1 << 32) + ((uint32_t) vert_0);
  uint64_t * p;
  JLG(p, sheaf_transport, edge_key);
  if (p != NULL)
    return (int64_t *) *p;
  return NULL;
}

// ── extract first two vertices from a simplex (JSL sorted order) ───

static void simplex_first_two_verts(
    void * s, uint64_t * v0, uint64_t * v1)
{
  uint8_t vv [C_1_KB];
  uint64_t * pv;
  vv[0] = 0; JSLF(pv, s, vv);
  *v0 = uint256_to_64(&vv[0]);
  JSLN(pv, s, vv);
  *v1 = uint256_to_64(&vv[0]);
}

// ── F_p coboundary builder ─────────────────────────────────────────

static int32_t sheaf_build_coboundary_fq(
    complex_cochain_t * Cc_orig, complex_cochain_t * Cc_u)
{
  uint32_t r = Cc_orig->sheaf_rank;
  uint64_t sid, cid;
  uint64_t * p_cofaces, * p_cf;

  sid = 0; JLF(p_cofaces, Cc_orig->C.cofaces, sid);
  while (p_cofaces != NULL)
  {
    uint64_t * sid_s, * cid_s;
    JLG(sid_s, Cc_orig->C.simplices, sid);

    cid = 0; JLF(p_cf, *p_cofaces, cid);
    while (p_cf != NULL)
    {
      JLG(cid_s, Cc_orig->C.simplices, cid);
      if (sid_s != NULL && cid_s != NULL)
      {
        int64_t sign = simplex_sign_in_coface(*sid_s, *cid_s);
        uint64_t pos = simplex_index_in_coface(*sid_s, *cid_s);

        if (pos == 0)
        {
          // position-0 face: use transport matrix T
          uint64_t v0, v1;
          simplex_first_two_verts(*cid_s, &v0, &v1);
          int64_t * T = get_sheaf_transport(v0, v1, Cc_orig->sheaf_transport);

          for (uint32_t i = 0; i < r; i++)
          {
            uint64_t vsid = sid * r + i;

            // look up or create coboundary for vsid
            uint64_t * q_cb;
            JLG(q_cb, Cc_u->elementary_coboundaries, vsid);
            void * v_coboundary;
            if (q_cb == NULL)
            {
              v_coboundary = NULL;
            }
            else
            {
              v_coboundary = (void *) *q_cb;
            }

            for (uint32_t j = 0; j < r; j++)
            {
              int64_t Tji;
              if (T != NULL)
                Tji = T[j * r + i]; // T[j][i] in row-major T
              else
                Tji = (i == j) ? 1 : 0; // identity

              if (Tji == 0)
                continue;

              int64_t coeff = sign * Tji;
              uint64_t vcid = cid * r + j;

              uint64_t * q_existing;
              JLG(q_existing, v_coboundary, vcid);
              if (q_existing == NULL)
              {
                fq_t * c = (fq_t *) malloc(sizeof(fq_t));
                fq_init(*c, Cc_u->fqctx);
                fq_set_si(*c, coeff, Cc_u->fqctx);
                uint64_t * q_new;
                JLI(q_new, v_coboundary, vcid);
                *q_new = (uint64_t) c;
              }
              else
              {
                fq_t tmp;
                fq_init(tmp, Cc_u->fqctx);
                fq_set_si(tmp, coeff, Cc_u->fqctx);
                fq_add(*((fq_t *) *q_existing), *((fq_t *) *q_existing), tmp, Cc_u->fqctx);
                fq_clear(tmp, Cc_u->fqctx);
                if (fq_is_zero(*((fq_t *) *q_existing), Cc_u->fqctx))
                {
                  fq_clear(*((fq_t *) *q_existing), Cc_u->fqctx);
                  free((fq_t *) *q_existing);
                  int rc;
                  JLD(rc, v_coboundary, vcid);
                }
              }
            }

            // store back
            uint64_t * q_store;
            JLI(q_store, Cc_u->elementary_coboundaries, vsid);
            *q_store = (uint64_t) v_coboundary;
          }
        }
        else
        {
          // non-position-0 face: identity block
          for (uint32_t i = 0; i < r; i++)
          {
            uint64_t vsid = sid * r + i;
            uint64_t vcid = cid * r + i;

            uint64_t * q_cb;
            JLG(q_cb, Cc_u->elementary_coboundaries, vsid);
            void * v_coboundary;
            if (q_cb == NULL)
              v_coboundary = NULL;
            else
              v_coboundary = (void *) *q_cb;

            uint64_t * q_existing;
            JLG(q_existing, v_coboundary, vcid);
            if (q_existing == NULL)
            {
              fq_t * c = (fq_t *) malloc(sizeof(fq_t));
              fq_init(*c, Cc_u->fqctx);
              fq_set_si(*c, sign, Cc_u->fqctx);
              uint64_t * q_new;
              JLI(q_new, v_coboundary, vcid);
              *q_new = (uint64_t) c;
            }
            else
            {
              fq_t tmp;
              fq_init(tmp, Cc_u->fqctx);
              fq_set_si(tmp, sign, Cc_u->fqctx);
              fq_add(*((fq_t *) *q_existing), *((fq_t *) *q_existing), tmp, Cc_u->fqctx);
              fq_clear(tmp, Cc_u->fqctx);
              if (fq_is_zero(*((fq_t *) *q_existing), Cc_u->fqctx))
              {
                fq_clear(*((fq_t *) *q_existing), Cc_u->fqctx);
                free((fq_t *) *q_existing);
                int rc;
                JLD(rc, v_coboundary, vcid);
              }
            }

            uint64_t * q_store;
            JLI(q_store, Cc_u->elementary_coboundaries, vsid);
            *q_store = (uint64_t) v_coboundary;
          }
        }
      }

      JLN(p_cf, *p_cofaces, cid);
    }

    JLN(p_cofaces, Cc_orig->C.cofaces, sid);
  }

  return 0;
}

// ── Z coboundary builder ───────────────────────────────────────────

static int32_t sheaf_build_coboundary_z(
    complex_cochain_t * Cc_orig, complex_cochain_t * Cc_u)
{
  uint32_t r = Cc_orig->sheaf_rank;
  uint64_t sid, cid;
  uint64_t * p_cofaces, * p_cf;

  sid = 0; JLF(p_cofaces, Cc_orig->C.cofaces, sid);
  while (p_cofaces != NULL)
  {
    uint64_t * sid_s, * cid_s;
    JLG(sid_s, Cc_orig->C.simplices, sid);

    cid = 0; JLF(p_cf, *p_cofaces, cid);
    while (p_cf != NULL)
    {
      JLG(cid_s, Cc_orig->C.simplices, cid);
      if (sid_s != NULL && cid_s != NULL)
      {
        int64_t sign = simplex_sign_in_coface(*sid_s, *cid_s);
        uint64_t pos = simplex_index_in_coface(*sid_s, *cid_s);

        if (pos == 0)
        {
          uint64_t v0, v1;
          simplex_first_two_verts(*cid_s, &v0, &v1);
          int64_t * T = get_sheaf_transport(v0, v1, Cc_orig->sheaf_transport);

          for (uint32_t i = 0; i < r; i++)
          {
            uint64_t vsid = sid * r + i;

            uint64_t * q_cb;
            JLG(q_cb, Cc_u->elementary_coboundaries_z, vsid);
            void * v_coboundary;
            if (q_cb == NULL)
              v_coboundary = NULL;
            else
              v_coboundary = (void *) *q_cb;

            for (uint32_t j = 0; j < r; j++)
            {
              int64_t Tji;
              if (T != NULL)
                Tji = T[j * r + i];
              else
                Tji = (i == j) ? 1 : 0;

              if (Tji == 0)
                continue;

              int64_t coeff = sign * Tji;
              uint64_t vcid = cid * r + j;

              uint64_t * q_existing;
              JLG(q_existing, v_coboundary, vcid);
              if (q_existing == NULL)
              {
                fmpz_t * c = (fmpz_t *) malloc(sizeof(fmpz_t));
                fmpz_init(*c);
                fmpz_set_si(*c, coeff);
                uint64_t * q_new;
                JLI(q_new, v_coboundary, vcid);
                *q_new = (uint64_t) c;
              }
              else
              {
                fmpz_t tmp;
                fmpz_init(tmp);
                fmpz_set_si(tmp, coeff);
                fmpz_add(*((fmpz_t *) *q_existing), *((fmpz_t *) *q_existing), tmp);
                fmpz_clear(tmp);
                if (fmpz_is_zero(*((fmpz_t *) *q_existing)))
                {
                  fmpz_clear(*((fmpz_t *) *q_existing));
                  free((fmpz_t *) *q_existing);
                  int rc;
                  JLD(rc, v_coboundary, vcid);
                }
              }
            }

            uint64_t * q_store;
            JLI(q_store, Cc_u->elementary_coboundaries_z, vsid);
            *q_store = (uint64_t) v_coboundary;
          }
        }
        else
        {
          for (uint32_t i = 0; i < r; i++)
          {
            uint64_t vsid = sid * r + i;
            uint64_t vcid = cid * r + i;

            uint64_t * q_cb;
            JLG(q_cb, Cc_u->elementary_coboundaries_z, vsid);
            void * v_coboundary;
            if (q_cb == NULL)
              v_coboundary = NULL;
            else
              v_coboundary = (void *) *q_cb;

            uint64_t * q_existing;
            JLG(q_existing, v_coboundary, vcid);
            if (q_existing == NULL)
            {
              fmpz_t * c = (fmpz_t *) malloc(sizeof(fmpz_t));
              fmpz_init(*c);
              fmpz_set_si(*c, sign);
              uint64_t * q_new;
              JLI(q_new, v_coboundary, vcid);
              *q_new = (uint64_t) c;
            }
            else
            {
              fmpz_t tmp;
              fmpz_init(tmp);
              fmpz_set_si(tmp, sign);
              fmpz_add(*((fmpz_t *) *q_existing), *((fmpz_t *) *q_existing), tmp);
              fmpz_clear(tmp);
              if (fmpz_is_zero(*((fmpz_t *) *q_existing)))
              {
                fmpz_clear(*((fmpz_t *) *q_existing));
                free((fmpz_t *) *q_existing);
                int rc;
                JLD(rc, v_coboundary, vcid);
              }
            }

            uint64_t * q_store;
            JLI(q_store, Cc_u->elementary_coboundaries_z, vsid);
            *q_store = (uint64_t) v_coboundary;
          }
        }
      }

      JLN(p_cf, *p_cofaces, cid);
    }

    JLN(p_cofaces, Cc_orig->C.cofaces, sid);
  }

  return 0;
}

// ── sheaf_unroll ───────────────────────────────────────────────────

int32_t sheaf_unroll(
    complex_cochain_t * Cc_orig, complex_cochain_t * Cc_u)
{
  uint32_t r = Cc_orig->sheaf_rank;
  if (r == 0)
    return -1;

  // 1. Initialize Cc_u with same field parameters
  if (Cc_orig->algorithm == COHOMOLOGY_INTEGRAL)
  {
    complex_cochain_init_z(Cc_u, Cc_orig->l_prime);
  }
  else
  {
    const fmpz * fp = fq_ctx_prime(Cc_orig->fqctx);
    uint32_t p = fmpz_get_ui(fp);
    uint32_t d = fq_ctx_degree(Cc_orig->fqctx);
    complex_cochain_init(Cc_u, p, d);
  }
  Cc_u->algorithm = Cc_orig->algorithm;
  Cc_u->l_prime = Cc_orig->l_prime;

  // 2. Populate virtual simplices
  uint64_t sid;
  uint64_t * p_s;
  uint64_t * q;

  uint32_t max_dim = 0;

  sid = 0; JLF(p_s, Cc_orig->C.simplices, sid);
  while (p_s != NULL)
  {
    uint64_t * p_dim, * p_age, * p_rad;
    JLG(p_dim, Cc_orig->C.dimension, sid);
    JLG(p_age, Cc_orig->C.age, sid);
    JLG(p_rad, Cc_orig->C.radius, sid);

    uint32_t dim = (p_dim != NULL) ? (uint32_t) *p_dim : 0;
    uint64_t age = (p_age != NULL) ? *p_age : 0;
    uint64_t rad = (p_rad != NULL) ? *p_rad : 0;

    if (dim > max_dim)
      max_dim = dim;

    for (uint32_t c = 0; c < r; c++)
    {
      uint64_t vid = sid * r + c;

      JLI(q, Cc_u->C.simplices, vid); *q = 1; // dummy simplex marker
      JLI(q, Cc_u->C.dimension, vid); *q = dim;
      JLI(q, Cc_u->C.age, vid); *q = age;
      JLI(q, Cc_u->C.radius, vid); *q = rad;
    }

    JLN(p_s, Cc_orig->C.simplices, sid);
  }

  Cc_u->C.dim_sk = max_dim;

  // Build scd (counter -> id), sci (id -> counter), scc per dimension
  for (uint32_t d = 0; d <= max_dim; d++)
  {
    uint64_t counter = 0;
    uint64_t vid = 0;
    uint64_t * p_dim;

    JLF(p_dim, Cc_u->C.dimension, vid);
    while (p_dim != NULL)
    {
      if (*p_dim == d)
      {
        JLI(q, Cc_u->C.scd[d], counter); *q = vid;
        JLI(q, Cc_u->C.sci[d], vid); *q = counter;
        counter++;
      }
      JLN(p_dim, Cc_u->C.dimension, vid);
    }
    Cc_u->C.scc[d] = counter;
  }

  // 3. Build sheaf coboundary
  if (Cc_orig->algorithm == COHOMOLOGY_INTEGRAL)
    sheaf_build_coboundary_z(Cc_orig, Cc_u);
  else
    sheaf_build_coboundary_fq(Cc_orig, Cc_u);

  return 0;
}

// ── sheaf_unroll_free ──────────────────────────────────────────────

int32_t sheaf_unroll_free(
    complex_cochain_t * Cc_u)
{
  complex_cochain_free(Cc_u);
  return 0;
}
