#include "common.h"
#include "semel.h"

fq_t * fq_one_(
    fq_ctx_t * fqctx)
{
  fq_t * x = (fq_t *) malloc(sizeof(fq_t));
  fq_init(*x, *fqctx);
  fq_set_si(*x, 1, *fqctx);
  return x;
}

int32_t Judy_free_1_a(
    void * ja)
{
  int rc;
  uint64_t id;
  uint64_t * p;

  id = 0; JLF(p, ja, id);
  while (p != NULL)
  {
    JLFA(rc, *p);
    JLN(p, ja, id);
  }
  JLFA(rc, ja);

  return 0;
}

int32_t Judy_free_2_a(
    void * ja)
{
  int rc;
  uint64_t id;
  uint64_t * p;

  id = 0; JLF(p, ja, id);
  while (p != NULL)
  {
    Judy_free_1_a(*p);
    JLN(p, ja, id);
  }
  JLFA(rc, ja);

  return 0;
}

int32_t Judy_free_1_b(
    void * ja, fq_ctx_t * fqctx)
{
  int rc;
  uint64_t i;
  uint64_t * p;

  i = 0; JLF(p, ja, i);
  while (p != NULL)
  {
    fq_t * p_ = (fq_t *) *p;
    fq_clear(*p_, *fqctx); free(p_);
    JLN(p, ja, i);
  }
  JLFA(rc, ja);

  return 0;
}

int32_t Judy_free_2_b(
    void * ja, fq_ctx_t * fqctx)
{
  int rc;
  uint64_t i, j;
  uint64_t * p, * q;

  i = 0; JLF(p, ja, i);
  while (p != NULL)
  {
    j = 0; JLF(q, *p, j);
    while (q != NULL)
    {
      fq_t * q_ = (fq_t *) *q;
      fq_clear(*q_, *fqctx); free(q_);
      JLN(q, *p, j);
    }
    JLFA(rc, *p);

    JLN(p, ja, i);
  }
  JLFA(rc, ja);

  return 0;
}

fmpz_t * fmpz_one_(void)
{
  fmpz_t * x = (fmpz_t *) malloc(sizeof(fmpz_t));
  fmpz_init(*x);
  fmpz_set_si(*x, 1);
  return x;
}

int32_t Judy_free_1_z(
    void * ja)
{
  int rc;
  uint64_t i;
  uint64_t * p;

  i = 0; JLF(p, ja, i);
  while (p != NULL)
  {
    fmpz_t * f = (fmpz_t *) *p;
    fmpz_clear(*f); free(f);
    JLN(p, ja, i);
  }
  JLFA(rc, ja);

  return 0;
}

int32_t Judy_free_2_z(
    void * ja)
{
  int rc;
  uint64_t i;
  uint64_t * p;

  i = 0; JLF(p, ja, i);
  while (p != NULL)
  {
    Judy_free_1_z(*p);
    JLN(p, ja, i);
  }
  JLFA(rc, ja);

  return 0;
}

