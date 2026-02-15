#include "common.h"
#include "semel.h"

// serialization for uint64_t so that
// lexicographic order is natural order on integers
int32_t uint64_to_256(
    uint64_t x, uint8_t * sx)
{
  mpz_t x_;
  mpz_init2(x_, NBITS);
  sx[0] = 0;
  mpz_set_ui(x_, x);
  gmp_sprintf((char *) sx, "<%064Zx>", x_);
  mpz_clear(x_);
  return 0;
}

uint64_t uint256_to_64(
    uint8_t * sx)
{
  mpz_t x_;
  mpz_init2(x_, NBITS);
  gmp_sscanf((char *) sx, "<%064Zx>", x_);
  uint64_t x = mpz_get_ui(x_);
  mpz_clear(x_);
  return x;
}

// serialization for simplices so that
// ordering by dimension and veritices
int32_t simplex_serialize(
    void * s, uint8_t * sx)
{
  uint8_t v [C_1_KB];
  char stmp0 [C_1_KB], stmp1 [C_1_KB];
  stmp0[0]= 0; stmp1[0]= 0;

  uint64_t k;
  uint64_t * p;

  v[0]= 0; JSLF(p, s, v); k = 0;
  while (p != NULL)
  {
    strcat(stmp0, (char *) &v[0]);
    JSLN(p, s, v); k++;
  }

  uint64_to_256(k, (uint8_t *) &stmp1[0]);
  sprintf((char *) sx, "%s[%s]", stmp1, stmp0);

  return 0;
}

