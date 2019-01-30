#include "common.h"
#include "semel.h"

int32_t ds_logistic(
    double lambda_,
    uint64_t np, double ** p)
{
  mpfr_t x, lambda, y0, y1, one;
  const int prec = 200;

  double x0 = ((double) rand()) / ((double) RAND_MAX);
  mpfr_init2(x, prec);
  mpfr_set_d(x, x0, MPFR_RNDD);

  mpfr_init2(lambda, prec);
  mpfr_set_d(lambda, lambda_, MPFR_RNDD);

  mpfr_init2(y0, prec);
  mpfr_init2(y1, prec);

  mpfr_init2(one, prec);
  mpfr_set_d(one, 1.0, MPFR_RNDD);

  *p = (double *) malloc(np * sizeof(double));

  for(uint64_t k=0; k<np; ++k)
  {
    mpfr_sub(y0, one, x, MPFR_RNDD);
    mpfr_mul(y1, x, y0, MPFR_RNDD);
    mpfr_mul(x, lambda, y1, MPFR_RNDD);
    (*p)[k] = mpfr_get_d(x, MPFR_RNDD);
  }

  mpfr_clear(x);
  mpfr_clear(lambda);
  mpfr_clear(y0);
  mpfr_clear(y1);
  mpfr_clear(one);
  mpfr_free_cache();

  return 0;
}
