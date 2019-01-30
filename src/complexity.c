#include "common.h"
#include "semel.h"

int32_t pairs_most_common(
    uint64_t nx, uint32_t * x,
    uint32_t * x0, uint32_t * x1, uint64_t * x01, uint32_t * x2)
{
  uint32_t xmax = 0;
  for(uint64_t i=0;i<nx;++i)
    { if (x[i] > xmax) { xmax = x[i]; } }
  *x2 = xmax + 1;

  void * p = NULL;
  int32_t rc;

  uint64_t * q;
  uint64_t tmp, tmpmax, k, kmax;
  for(uint64_t i=0;i<nx;i+=2)
  {
    tmp = ((uint64_t) x[i]) | (((uint64_t) x[i+1]) << 32);
    k = 0;
    JLG(q, p, tmp);
    if (q != NULL) { k = *q; }
    JLI(q, p, tmp); *q = k + 1;
  }

  kmax = 0; tmpmax = 0;
  tmp = 0; JLF(q, p, tmp);
  while (q != NULL)
  {
    k = *q;
    if (k > kmax)
      { kmax = k; tmpmax = tmp; }
    JLN(q, p, tmp);
  }

  *x0 = (uint32_t) tmpmax;
  *x1 = (uint32_t) (tmpmax >> 32);
  *x01 = kmax;

  JLFA(rc, p);

  return 0;
}

int32_t pairs_substitution(
    uint64_t nx, uint32_t * x,
    uint32_t x0, uint32_t x1, uint32_t x2)
{
  uint64_t i, j;
  for(i=0, j=0;i<nx;i+=2)
  {
    if ((x[i] == x0) && (x[i+1] == x1))
      { x[j] = x2; ++j; }
    else
      { x[j] = x[i]; x[j+1] = x[i+1]; j+= 2; }
  }

  if ((nx % 2) == 1)
    { x[j] = x[i]; ++j; }

  for(i=j;i<nx;++i)
    { x[i] = 0; }

  return 0;
}

int32_t nsrps(
    uint64_t nx, uint32_t * x,
    uint64_t * c, semel_ctx_log_t * log)
{
  uint32_t x0, x1, x2;
  uint64_t x01;
  if (nx > 1)
  {
    pairs_most_common(nx, x, &x0, &x1, &x01, &x2);
    pairs_substitution(nx, x, x0, x1, x2);
    nx -= x01; ++(*c);

    nsrps(nx, x, c, log);
  }

  return 0;
}

static int32_t seq_iminmax_16(
  uint64_t nx, uint16_t * x,
  uint64_t * imin, uint64_t * imax)
{
  *imin = 0; *imax = 0;
  uint16_t min = x[*imin];
  uint16_t max = x[*imax];

  uint64_t i;
  for(i=1;i<nx;i++)
  {
    if (x[i] < min)
      { *imin = i; min = x[i]; }
    if (x[i] > max)
      { *imax = i; max = x[i]; }
  }

  return 0;
}

static int32_t seq_emp_dist_16(
  uint64_t nx, uint16_t * x, uint16_t * y)
{
  uint64_t ny = (1 << 16);
  uint64_t i;
  for(i=0;i<ny;i++)
    { y[i] = 0; }
  for(i=0;i<nx;i++)
    { y[x[i]]++; }
  return 0;
}

static int32_t seq_group_8to16(
    uint64_t nx, uint8_t * x, uint16_t * y)
{
  uint64_t i, j;
  j = 0;
  for(i=0;i<nx-2;i+=2)
    { y[j++] = (x[i+1] << 8) | x[i]; }
  return 0;
}

int32_t entropy_dist(
    uint64_t np, double * p, double alpha, double * e)
{
  double e_shannon, e_collision, e_min, e_renyi;
  double psum, pmin, pmax;

  psum = 0.0; pmax = 0.0; pmin = 1.0;
  e_shannon = 0.0; e_collision = 0.0; e_min = 0.0; e_renyi = 0.0;
  for(uint64_t i=0;i<np;i++)
  {
    psum += p[i];
    if (p[i] < pmin) { pmin = p[i]; }
    if (p[i] > pmax) { pmax = p[i]; }
    if (p[i] > 0.0) { e_shannon += p[i] * log2(p[i]); }
    e_collision += p[i]*p[i];
    e_renyi += pow(p[i], alpha);
  }
  e_shannon = -e_shannon;
  e_collision = -log2(e_collision);
  e_min = -log2(pmax);
  e_renyi = (1.0 / (1.0 - alpha)) * log2(e_renyi);

  e[0] = e_shannon;
  e[1] = e_collision;
  e[2] = e_min;
  e[3] = e_renyi;

  return 0;
}

int32_t dist_words_seq_bin(
    uint64_t nx, uint8_t * x, uint64_t w,
    uint64_t * np, double ** p)
{
  uint64_t i, j;
  uint64_t itmp;

  if (w < 32)
  {
    *np = 1 << w;
    *p = (double *) malloc((*np) * sizeof(double));

    for(i=0;i<(*np);i++)
      { (*p)[i] = 0; }

    for(i=0;i<nx-w;i++)
    {
      itmp = 0;
      for(j=0;j<w;j++)
        { itmp |= (x[i+j] << j); }
      (*p)[itmp]++;
    }

    for(i=0;i<(*np);i++)
      { (*p)[i] /= ((double) nx); }
  }

  return 0;
}

int32_t nrps_seq_bin(
    uint64_t nx, uint8_t * x)
{
  uint64_t ny = nx / 2;
  uint16_t * y = (uint16_t *) malloc(ny * sizeof(uint16_t));
  seq_group_8to16(nx, x, y);
  uint64_t nyd = (1 << 16);
  uint16_t * yd = (uint16_t *) malloc(nyd * sizeof(uint16_t));
  seq_emp_dist_16(ny, y, yd);
  uint64_t imin, imax;
  seq_iminmax_16(nyd, yd, &imin, &imax);
  free(y);
  free(yd);
  return 0;
}
