#include "common.h"
#include "semel.h"

// simple dynamic string buffer for building JSON without jansson

typedef struct {
  char * buf;
  size_t len;
  size_t cap;
} strbuf_t;

static void strbuf_init(strbuf_t * sb, size_t initial_cap)
{
  sb->cap = (initial_cap > 64) ? initial_cap : 64;
  sb->buf = (char *) malloc(sb->cap);
  sb->buf[0] = '\0';
  sb->len = 0;
}

static void strbuf_ensure(strbuf_t * sb, size_t extra)
{
  while (sb->len + extra + 1 > sb->cap)
  {
    sb->cap *= 2;
    sb->buf = (char *) realloc(sb->buf, sb->cap);
  }
}

static void strbuf_append(strbuf_t * sb, const char * s)
{
  size_t slen = strlen(s);
  strbuf_ensure(sb, slen);
  memcpy(sb->buf + sb->len, s, slen + 1);
  sb->len += slen;
}

static void strbuf_append_int(strbuf_t * sb, int64_t v)
{
  char tmp[32];
  snprintf(tmp, sizeof(tmp), "%" PRId64, v);
  strbuf_append(sb, tmp);
}

static void strbuf_append_double(strbuf_t * sb, double v)
{
  char tmp[64];
  snprintf(tmp, sizeof(tmp), "%.8g", v);
  strbuf_append(sb, tmp);
}

static char * strbuf_detach(strbuf_t * sb)
{
  return sb->buf;
}

// export functions

char * export_cocycles(
    complex_cochain_t * Cc)
{
  uint64_t v;
  uint64_t * p;

  strbuf_t sb;
  strbuf_init(&sb, 1024);
  strbuf_append(&sb, "[");

  int first_cocycle = 1;
  v = 0; JLF(p, Cc->cocycles, v);
  while (p != NULL)
  {
    if (!first_cocycle) { strbuf_append(&sb, ","); }
    first_cocycle = 0;

    void * alpha = (void *) *p;
    uint64_t dim = elementary_dimension(alpha, Cc);

    strbuf_append(&sb, "{\"dimension\":");
    strbuf_append_int(&sb, dim);
    strbuf_append(&sb, ",\"expression\":[");

    uint64_t k;
    uint64_t * q;
    int first_term = 1;
    k = 0; JLF(q, alpha, k);
    while (q != NULL)
    {
      char * tmp = fq_get_str_pretty(*((fq_t *) *q), Cc->fqctx);
      int32_t c = atoi(tmp);
      flint_free(tmp);
      fmpz_t f; fq_ctx_order(f, Cc->fqctx);
      uint32_t o = fmpz_get_ui(f);
      fmpz_clear(f);
      c = (c > (int32_t) (o >> 1)) ? (c - o) : c;

      if (!first_term) { strbuf_append(&sb, ","); }
      first_term = 0;

      strbuf_append(&sb, "[");
      strbuf_append_int(&sb, c);
      strbuf_append(&sb, ",");
      strbuf_append_int(&sb, k);
      strbuf_append(&sb, "]");

      JLN(q, alpha, k);
    }

    strbuf_append(&sb, "]}");
    JLN(p, Cc->cocycles, v);
  }

  strbuf_append(&sb, "]");
  return strbuf_detach(&sb);
}

char * export_complex_simplicial(
    complex_simplicial_t * C)
{
  uint64_t * p, * q;
  uint64_t k;
  uint8_t v [C_1_KB];

  strbuf_t sb;
  strbuf_init(&sb, 4096);

  uint64_t n;
  JLC(n, C->simplices, 0, -1);

  strbuf_append(&sb, "{\"number-of-simplices\":");
  strbuf_append_int(&sb, n);
  strbuf_append(&sb, ",\"complex\":[");

  int first_simplex = 1;
  k = 0; JLF(p, C->simplices, k);
  while (p != NULL)
  {
    if (!first_simplex) { strbuf_append(&sb, ","); }
    first_simplex = 0;

    JLG(q, C->age, k);
    strbuf_append(&sb, "{\"age\":");
    strbuf_append_int(&sb, *q);

    JLG(q, C->dimension, k);
    strbuf_append(&sb, ",\"dimension\":");
    strbuf_append_int(&sb, *q);

    JLG(q, C->sci[simplex_dimension(*p)], k);
    strbuf_append(&sb, ",\"index\":");
    strbuf_append_int(&sb, *q);

    strbuf_append(&sb, ",\"simplex\":[");
    int first_vert = 1;
    v[0] = 0; JSLF(q, *p, v);
    while (q != NULL)
    {
      if (!first_vert) { strbuf_append(&sb, ","); }
      first_vert = 0;

      if (C->typeS == 0)
      {
        strbuf_append_int(&sb, uint256_to_64(v));
      }
      if (C->typeS == 1)
      {
        strbuf_append(&sb, "\"");
        strbuf_append(&sb, (char *) &v[0]);
        strbuf_append(&sb, "\"");
      }

      JSLN(q, *p, v);
    }

    strbuf_append(&sb, "]}");
    JLN(p, C->simplices, k);
  }

  strbuf_append(&sb, "]}");
  return strbuf_detach(&sb);
}

char * export_cohomology(
    complex_cochain_t * Cc, uint8_t measures_only)
{
  uint64_t * alpha, * sigma;
  uint64_t alpha_dim, * p, * j_dim, * j_radius, * k_radius;
  uint64_t i, j;

  // collect cocycles by dimension: dim -> list of birth radii
  // collect intervals by dimension: dim -> list of (birth, death) pairs
  // use MAX_DIM-sized arrays

  uint32_t max_dim_seen = 0;

  // first pass: count intervals per dimension
  uint64_t interval_count[MAX_DIM];
  memset(interval_count, 0, sizeof(interval_count));

  j = 0; JLF(p, Cc->intervals, j);
  while (p != NULL)
  {
    JLG(j_dim, Cc->C.dimension, j);
    uint32_t d = (uint32_t) *j_dim;
    interval_count[d]++;
    if (d > max_dim_seen) { max_dim_seen = d; }
    JLN(p, Cc->intervals, j);
  }

  // allocate interval arrays
  double * interval_birth[MAX_DIM];
  double * interval_death[MAX_DIM];
  uint64_t interval_idx[MAX_DIM];
  memset(interval_birth, 0, sizeof(interval_birth));
  memset(interval_death, 0, sizeof(interval_death));
  memset(interval_idx, 0, sizeof(interval_idx));

  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (interval_count[d] > 0)
    {
      interval_birth[d] = (double *) malloc(interval_count[d] * sizeof(double));
      interval_death[d] = (double *) malloc(interval_count[d] * sizeof(double));
    }
  }

  // second pass: fill interval arrays
  j = 0; JLF(p, Cc->intervals, j);
  while (p != NULL)
  {
    uint64_t k = *p;
    JLG(j_dim, Cc->C.dimension, j);
    JLG(j_radius, Cc->C.radius, j);
    JLG(k_radius, Cc->C.radius, k);

    uint32_t d = (uint32_t) *j_dim;
    uint64_t idx = interval_idx[d]++;
    interval_birth[d][idx] = ((double) *j_radius) / PRECISION;
    interval_death[d][idx] = ((double) *k_radius) / PRECISION;

    JLN(p, Cc->intervals, j);
  }

  // count cocycles per dimension
  uint64_t cocycle_count[MAX_DIM];
  memset(cocycle_count, 0, sizeof(cocycle_count));

  i = 0; JLF(alpha, Cc->cocycles, i);
  while (alpha != NULL)
  {
    alpha_dim = elementary_dimension(*alpha, Cc);
    cocycle_count[alpha_dim]++;
    if (alpha_dim > max_dim_seen) { max_dim_seen = (uint32_t) alpha_dim; }
    JLN(alpha, Cc->cocycles, i);
  }

  // allocate cocycle birth arrays
  double * cocycle_birth[MAX_DIM];
  uint64_t cocycle_idx[MAX_DIM];
  memset(cocycle_birth, 0, sizeof(cocycle_birth));
  memset(cocycle_idx, 0, sizeof(cocycle_idx));

  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (cocycle_count[d] > 0)
      { cocycle_birth[d] = (double *) malloc(cocycle_count[d] * sizeof(double)); }
  }

  // fill cocycle birth arrays
  i = 0; JLF(alpha, Cc->cocycles, i);
  while (alpha != NULL)
  {
    alpha_dim = elementary_dimension(*alpha, Cc);
    j = 0; JLF(sigma, *alpha, j);
    JLG(j_radius, Cc->C.radius, j);

    uint64_t idx = cocycle_idx[alpha_dim]++;
    cocycle_birth[alpha_dim][idx] = ((double) *j_radius) / PRECISION;

    JLN(alpha, Cc->cocycles, i);
  }

  // compute mass and entropy per dimension
  double mass[MAX_DIM];
  double entropy[MAX_DIM];
  memset(mass, 0, sizeof(mass));
  memset(entropy, 0, sizeof(entropy));

  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    for (uint64_t k = 0; k < interval_count[d]; ++k)
      { mass[d] += interval_death[d][k] - interval_birth[d][k]; }

    if (mass[d] > 0.0)
    {
      for (uint64_t k = 0; k < interval_count[d]; ++k)
      {
        double tmp = (interval_death[d][k] - interval_birth[d][k]) / mass[d];
        if (tmp > 0.0)
          entropy[d] += - tmp * log(tmp);
      }
    }
  }

  // build JSON output
  strbuf_t sb;
  strbuf_init(&sb, 4096);

  strbuf_append(&sb, "{\"stats\":{\"mass\":{");
  int first = 1;
  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (interval_count[d] > 0)
    {
      if (!first) { strbuf_append(&sb, ","); }
      first = 0;
      strbuf_append(&sb, "\"");
      strbuf_append_int(&sb, d);
      strbuf_append(&sb, "\":");
      strbuf_append_double(&sb, mass[d]);
    }
  }

  strbuf_append(&sb, "},\"entropy\":{");
  first = 1;
  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (interval_count[d] > 0)
    {
      if (!first) { strbuf_append(&sb, ","); }
      first = 0;
      strbuf_append(&sb, "\"");
      strbuf_append_int(&sb, d);
      strbuf_append(&sb, "\":");
      strbuf_append_double(&sb, entropy[d]);
    }
  }
  strbuf_append(&sb, "}}");

  if (measures_only == 0)
  {
    // cocycles
    strbuf_append(&sb, ",\"cocycles\":{");
    first = 1;
    for (uint32_t d = 0; d <= max_dim_seen; ++d)
    {
      if (cocycle_count[d] > 0)
      {
        if (!first) { strbuf_append(&sb, ","); }
        first = 0;
        strbuf_append(&sb, "\"");
        strbuf_append_int(&sb, d);
        strbuf_append(&sb, "\":[");
        for (uint64_t k = 0; k < cocycle_count[d]; ++k)
        {
          if (k > 0) { strbuf_append(&sb, ","); }
          strbuf_append(&sb, "[");
          strbuf_append_double(&sb, cocycle_birth[d][k]);
          strbuf_append(&sb, "]");
        }
        strbuf_append(&sb, "]");
      }
    }
    strbuf_append(&sb, "}");

    // intervals
    strbuf_append(&sb, ",\"intervals\":{");
    first = 1;
    for (uint32_t d = 0; d <= max_dim_seen; ++d)
    {
      if (interval_count[d] > 0)
      {
        if (!first) { strbuf_append(&sb, ","); }
        first = 0;
        strbuf_append(&sb, "\"");
        strbuf_append_int(&sb, d);
        strbuf_append(&sb, "\":[");
        for (uint64_t k = 0; k < interval_count[d]; ++k)
        {
          if (k > 0) { strbuf_append(&sb, ","); }
          strbuf_append(&sb, "[");
          strbuf_append_double(&sb, interval_birth[d][k]);
          strbuf_append(&sb, ",");
          strbuf_append_double(&sb, interval_death[d][k]);
          strbuf_append(&sb, "]");
        }
        strbuf_append(&sb, "]");
      }
    }
    strbuf_append(&sb, "}");
  }

  strbuf_append(&sb, "}");

  // cleanup
  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    free(interval_birth[d]);
    free(interval_death[d]);
    free(cocycle_birth[d]);
  }

  return strbuf_detach(&sb);
}

// ── integer cohomology export ────────────────────────────────────────

static int is_l_power(uint64_t n, uint32_t l)
{
  if (l < 2 || n == 0) return 0;
  while (n > 1)
  {
    if (n % l != 0) return 0;
    n /= l;
  }
  return 1;
}

char * export_cocycles_z(
    complex_cochain_t * Cc)
{
  uint64_t v;
  uint64_t * p;

  strbuf_t sb;
  strbuf_init(&sb, 1024);
  strbuf_append(&sb, "[");

  int first_cocycle = 1;
  v = 0; JLF(p, Cc->cocycles_z, v);
  while (p != NULL)
  {
    if (!first_cocycle) { strbuf_append(&sb, ","); }
    first_cocycle = 0;

    void * alpha = (void *) *p;
    uint64_t dim = elementary_dimension(alpha, Cc);

    strbuf_append(&sb, "{\"dimension\":");
    strbuf_append_int(&sb, dim);
    strbuf_append(&sb, ",\"expression\":[");

    uint64_t k;
    uint64_t * q;
    int first_term = 1;
    k = 0; JLF(q, alpha, k);
    while (q != NULL)
    {
      fmpz_t * fval = (fmpz_t *) *q;
      int64_t c = fmpz_get_si(*fval);

      if (!first_term) { strbuf_append(&sb, ","); }
      first_term = 0;

      strbuf_append(&sb, "[");
      strbuf_append_int(&sb, c);
      strbuf_append(&sb, ",");
      strbuf_append_int(&sb, k);
      strbuf_append(&sb, "]");

      JLN(q, alpha, k);
    }

    strbuf_append(&sb, "]}");
    JLN(p, Cc->cocycles_z, v);
  }

  strbuf_append(&sb, "]");
  return strbuf_detach(&sb);
}

char * export_cohomology_z(
    complex_cochain_t * Cc, uint8_t measures_only)
{
  uint64_t * alpha, * sigma;
  uint64_t alpha_dim, * p, * j_dim, * j_radius, * k_radius;
  uint64_t i, j;

  uint32_t max_dim_seen = 0;

  // first pass: count intervals per dimension and collect torsion
  uint64_t interval_count[MAX_DIM];
  memset(interval_count, 0, sizeof(interval_count));

  j = 0; JLF(p, Cc->intervals, j);
  while (p != NULL)
  {
    JLG(j_dim, Cc->C.dimension, j);
    uint32_t d = (uint32_t) *j_dim;
    interval_count[d]++;
    if (d > max_dim_seen) { max_dim_seen = d; }
    JLN(p, Cc->intervals, j);
  }

  // allocate interval arrays (birth, death, torsion_order)
  double * interval_birth[MAX_DIM];
  double * interval_death[MAX_DIM];
  uint64_t * interval_torsion[MAX_DIM];
  uint64_t interval_idx[MAX_DIM];
  memset(interval_birth, 0, sizeof(interval_birth));
  memset(interval_death, 0, sizeof(interval_death));
  memset(interval_torsion, 0, sizeof(interval_torsion));
  memset(interval_idx, 0, sizeof(interval_idx));

  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (interval_count[d] > 0)
    {
      interval_birth[d] = (double *) malloc(interval_count[d] * sizeof(double));
      interval_death[d] = (double *) malloc(interval_count[d] * sizeof(double));
      interval_torsion[d] = (uint64_t *) malloc(interval_count[d] * sizeof(uint64_t));
    }
  }

  // second pass: fill interval arrays
  j = 0; JLF(p, Cc->intervals, j);
  while (p != NULL)
  {
    uint64_t k = *p;
    JLG(j_dim, Cc->C.dimension, j);
    JLG(j_radius, Cc->C.radius, j);
    JLG(k_radius, Cc->C.radius, k);

    uint32_t d = (uint32_t) *j_dim;
    uint64_t idx = interval_idx[d]++;
    interval_birth[d][idx] = ((double) *j_radius) / PRECISION;
    interval_death[d][idx] = ((double) *k_radius) / PRECISION;

    // get torsion order
    uint64_t * t_p;
    JLG(t_p, Cc->intervals_torsion, j);
    interval_torsion[d][idx] = (t_p != NULL) ? *t_p : 0;

    JLN(p, Cc->intervals, j);
  }

  // count cocycles per dimension (surviving = infinite bars)
  uint64_t cocycle_count[MAX_DIM];
  memset(cocycle_count, 0, sizeof(cocycle_count));

  i = 0; JLF(alpha, Cc->cocycles_z, i);
  while (alpha != NULL)
  {
    alpha_dim = elementary_dimension(*alpha, Cc);
    cocycle_count[alpha_dim]++;
    if (alpha_dim > max_dim_seen) { max_dim_seen = (uint32_t) alpha_dim; }
    JLN(alpha, Cc->cocycles_z, i);
  }

  // allocate cocycle birth arrays
  double * cocycle_birth[MAX_DIM];
  uint64_t cocycle_idx[MAX_DIM];
  memset(cocycle_birth, 0, sizeof(cocycle_birth));
  memset(cocycle_idx, 0, sizeof(cocycle_idx));

  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (cocycle_count[d] > 0)
      { cocycle_birth[d] = (double *) malloc(cocycle_count[d] * sizeof(double)); }
  }

  // fill cocycle birth arrays
  i = 0; JLF(alpha, Cc->cocycles_z, i);
  while (alpha != NULL)
  {
    alpha_dim = elementary_dimension(*alpha, Cc);
    j = 0; JLF(sigma, *alpha, j);
    JLG(j_radius, Cc->C.radius, j);

    uint64_t idx = cocycle_idx[alpha_dim]++;
    cocycle_birth[alpha_dim][idx] = ((double) *j_radius) / PRECISION;

    JLN(alpha, Cc->cocycles_z, i);
  }

  // compute betti numbers and torsion invariant factors per dimension
  // betti[d] = number of free generators (cocycles with torsion_order=0 + infinite bars)
  // torsion_factors[d] = list of invariant factors > 1

  // compute mass and entropy per dimension (same as field version)
  double mass[MAX_DIM];
  double entropy[MAX_DIM];
  memset(mass, 0, sizeof(mass));
  memset(entropy, 0, sizeof(entropy));

  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    for (uint64_t k = 0; k < interval_count[d]; ++k)
      { mass[d] += interval_death[d][k] - interval_birth[d][k]; }

    if (mass[d] > 0.0)
    {
      for (uint64_t k = 0; k < interval_count[d]; ++k)
      {
        double tmp = (interval_death[d][k] - interval_birth[d][k]) / mass[d];
        if (tmp > 0.0)
          entropy[d] += - tmp * log(tmp);
      }
    }
  }

  // build JSON output
  strbuf_t sb;
  strbuf_init(&sb, 4096);

  // stats
  strbuf_append(&sb, "{\"stats\":{\"mass\":{");
  int first = 1;
  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (interval_count[d] > 0)
    {
      if (!first) { strbuf_append(&sb, ","); }
      first = 0;
      strbuf_append(&sb, "\"");
      strbuf_append_int(&sb, d);
      strbuf_append(&sb, "\":");
      strbuf_append_double(&sb, mass[d]);
    }
  }

  strbuf_append(&sb, "},\"entropy\":{");
  first = 1;
  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    if (interval_count[d] > 0)
    {
      if (!first) { strbuf_append(&sb, ","); }
      first = 0;
      strbuf_append(&sb, "\"");
      strbuf_append_int(&sb, d);
      strbuf_append(&sb, "\":");
      strbuf_append_double(&sb, entropy[d]);
    }
  }
  strbuf_append(&sb, "}}");

  if (measures_only == 0)
  {
    // l field (if l-adic filtering)
    if (Cc->l_prime > 0)
    {
      strbuf_append(&sb, ",\"l\":");
      strbuf_append_int(&sb, Cc->l_prime);
    }

    // betti numbers
    strbuf_append(&sb, ",\"betti\":{");
    first = 1;
    for (uint32_t d = 0; d <= max_dim_seen; ++d)
    {
      // betti = surviving cocycles (all are free in the integral sense here)
      uint64_t betti_d = cocycle_count[d];
      if (betti_d > 0 || interval_count[d] > 0)
      {
        if (!first) { strbuf_append(&sb, ","); }
        first = 0;
        strbuf_append(&sb, "\"");
        strbuf_append_int(&sb, d);
        strbuf_append(&sb, "\":");
        strbuf_append_int(&sb, betti_d);
      }
    }
    strbuf_append(&sb, "}");

    // torsion: invariant factors > 1 per dimension
    strbuf_append(&sb, ",\"torsion\":{");
    first = 1;
    for (uint32_t d = 0; d <= max_dim_seen; ++d)
    {
      // collect torsion factors for this dimension
      int has_torsion = 0;
      for (uint64_t k = 0; k < interval_count[d]; ++k)
      {
        uint64_t t = interval_torsion[d][k];
        if (t > 0)
        {
          if (Cc->l_prime > 0 && !is_l_power(t, Cc->l_prime))
            continue;
          has_torsion = 1;
          break;
        }
      }

      if (has_torsion)
      {
        if (!first) { strbuf_append(&sb, ","); }
        first = 0;
        strbuf_append(&sb, "\"");
        strbuf_append_int(&sb, d);
        strbuf_append(&sb, "\":[");
        int first_t = 1;
        for (uint64_t k = 0; k < interval_count[d]; ++k)
        {
          uint64_t t = interval_torsion[d][k];
          if (t > 0)
          {
            if (Cc->l_prime > 0 && !is_l_power(t, Cc->l_prime))
              continue;
            if (!first_t) { strbuf_append(&sb, ","); }
            first_t = 0;
            strbuf_append_int(&sb, t);
          }
        }
        strbuf_append(&sb, "]");
      }
    }
    strbuf_append(&sb, "}");

    // cocycles (infinite bars)
    strbuf_append(&sb, ",\"cocycles\":{");
    first = 1;
    for (uint32_t d = 0; d <= max_dim_seen; ++d)
    {
      if (cocycle_count[d] > 0)
      {
        if (!first) { strbuf_append(&sb, ","); }
        first = 0;
        strbuf_append(&sb, "\"");
        strbuf_append_int(&sb, d);
        strbuf_append(&sb, "\":[");
        for (uint64_t k = 0; k < cocycle_count[d]; ++k)
        {
          if (k > 0) { strbuf_append(&sb, ","); }
          strbuf_append(&sb, "[");
          strbuf_append_double(&sb, cocycle_birth[d][k]);
          strbuf_append(&sb, ",0]");
        }
        strbuf_append(&sb, "]");
      }
    }
    strbuf_append(&sb, "}");

    // intervals with torsion annotation
    strbuf_append(&sb, ",\"intervals\":{");
    first = 1;
    for (uint32_t d = 0; d <= max_dim_seen; ++d)
    {
      if (interval_count[d] > 0)
      {
        if (!first) { strbuf_append(&sb, ","); }
        first = 0;
        strbuf_append(&sb, "\"");
        strbuf_append_int(&sb, d);
        strbuf_append(&sb, "\":[");
        for (uint64_t k = 0; k < interval_count[d]; ++k)
        {
          if (k > 0) { strbuf_append(&sb, ","); }
          strbuf_append(&sb, "[");
          strbuf_append_double(&sb, interval_birth[d][k]);
          strbuf_append(&sb, ",");
          strbuf_append_double(&sb, interval_death[d][k]);
          strbuf_append(&sb, ",");
          strbuf_append_int(&sb, interval_torsion[d][k]);
          strbuf_append(&sb, "]");
        }
        strbuf_append(&sb, "]");
      }
    }
    strbuf_append(&sb, "}");
  }

  strbuf_append(&sb, "}");

  // cleanup
  for (uint32_t d = 0; d <= max_dim_seen; ++d)
  {
    free(interval_birth[d]);
    free(interval_death[d]);
    free(interval_torsion[d]);
    free(cocycle_birth[d]);
  }

  return strbuf_detach(&sb);
}

char * export_coboundary_operator(
    uint32_t k, complex_cochain_t * Cc)
{
  uint64_t v, w;
  uint64_t * p, * q, * p_;

  coboundary_operator(Cc);

  strbuf_t sb;
  strbuf_init(&sb, 1024);

  strbuf_append(&sb, "{\"m\":");
  strbuf_append_int(&sb, Cc->C.scc[k+1]);
  strbuf_append(&sb, ",\"n\":");
  strbuf_append_int(&sb, Cc->C.scc[k]);
  strbuf_append(&sb, ",\"matrix_sparse\":[");

  int first = 1;
  JLG(p, Cc->D, k);
  v = 0; JLF(q, *p, v);
  while (q != NULL)
  {
    w = 0; JLF(p_, *q, w);
    while (p_ != NULL)
    {
      if (!first) { strbuf_append(&sb, ","); }
      first = 0;

      strbuf_append(&sb, "[");
      strbuf_append_int(&sb, w);
      strbuf_append(&sb, ",");
      strbuf_append_int(&sb, v);
      strbuf_append(&sb, ",");
      strbuf_append_int(&sb, *p_);
      strbuf_append(&sb, "]");

      JLN(p_, *q, w);
    }
    JLN(q, *p, v);
  }

  strbuf_append(&sb, "]}");
  return strbuf_detach(&sb);
}
