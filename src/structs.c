#include "common.h"

int32_t point_cloud_init(
    point_cloud_t * P,
    uint32_t d, uint32_t n,
    uint32_t manifold, uint32_t distance,
    uint32_t nb, double * b)
{
  P->d = d;
  P->n = n;
  P->p = NULL;
  if ((d * n) > 0)
    { P->p = (double *) malloc((P->d * P->n) * sizeof(double)); }
  P->pn = NULL;

  P->n_ = 0;

  P->manifold = manifold;
  P->distance = distance;
  P->nb = nb;
  P->b = NULL;
  if (nb > 0)
  {
    P->b = (double *) malloc(P->nb * sizeof(double));
    for(uint32_t i=0;i<nb;++i)
      { P->b[i] = b[i]; }
  }

  P->nA = 0;
  P->A = NULL;
  P->Amax = 0.0;

  return 0;
}

int32_t point_cloud_free(
    point_cloud_t * P)
{
  int rc;

  if (P->p != NULL) { free(P->p); }
  JSLFA(rc, P->pn);

  if (P->b != NULL) { free(P->b); }

  if (P->A != NULL) { free(P->A); }

  return 0;
}

int32_t complex_simplicial_init(
    complex_simplicial_t * C)
{
  C->S = NULL;
  C->typeS = 0;

  C->G = NULL;

  C->filtration = NULL;

  C->simplices = NULL;
  C->sx2id = NULL;
  C->faces = NULL;
  C->cofaces = NULL;

  C->dimension = NULL;
  C->age = NULL;
  C->radius = NULL;

  //C->nr = 0;
  //C->r = NULL;

  C->dim_sk = 0;
  for(uint32_t i = 0; i < MAX_DIM; i++)
    { C->scd[i] = NULL; C->sci[i] = NULL; C->scc[i] = 0; }
  C->scc_ = 0;

  return 0;
}

int32_t complex_simplicial_free(
    complex_simplicial_t * C)
{
  int rc;

  uint64_t * p, * q;
  uint64_t r;
  uint8_t sx [MAX_LEN_SIMPLEX_SERIALIZATION];

  JLFA(rc, C->S);
  JLFA(rc, C->G);

  r = 0; JLF(p, C->filtration, r);
  while (p != NULL)
  {
    sx[0] = 0; JSLF(q, *p, sx);
    while (q != NULL)
    {
      JSLFA(rc, *q);
      JSLN(q, *p, sx);
    }
    JSLFA(rc, *p);

    JLN(p, C->filtration, r);
  }
  JLFA(rc, C->filtration);
  
  JLFA(rc, C->simplices);
  JLFA(rc, C->sx2id);

  Judy_free_1_a(C->faces);
  Judy_free_1_a(C->cofaces);

  JLFA(rc, C->dimension);
  JLFA(rc, C->age);
  JLFA(rc, C->radius);

  //if (C->r != NULL) { free(C->r); C->r = NULL; }

  for(uint32_t i = 0; i < MAX_DIM; i++)
    { JLFA(rc, C->scd[i]); JLFA(rc, C->sci[i]); }

  return 0;
}

int32_t complex_cochain_init(
    complex_cochain_t * Cc, uint32_t fp, uint32_t fd)
{
  complex_simplicial_init(&Cc->C);

  fmpz_t fp_;
  fmpz_init(fp_);
  fmpz_init_set_ui(fp_, fp);
  slong fd_ = (slong) fd;
  fq_ctx_init(Cc->fqctx, fp_, fd_, "F");

  Cc->D = NULL;

  Cc->elementary_boundaries = NULL;
  Cc->elementary_coboundaries = NULL;

  Cc->cocycles = NULL;
  //Cc->cocycles_vanished = NULL;
  Cc->intervals = NULL;

  Cc->cocycles_harmonic = NULL;

  Cc->algorithm = COHOMOLOGY_FIELD;
  Cc->l_prime = 0;

  Cc->elementary_coboundaries_z = NULL;
  Cc->cocycles_z = NULL;
  Cc->intervals_torsion = NULL;

  Cc->local_system = NULL;

  Cc->sheaf_rank = 0;
  Cc->sheaf_transport = NULL;

  return 0;
}

int32_t complex_cochain_init_z(
    complex_cochain_t * Cc, uint32_t l)
{
  complex_simplicial_init(&Cc->C);

  // init a dummy fq context (p=2, d=1) for safety — not used in integral path
  fmpz_t fp_;
  fmpz_init(fp_);
  fmpz_init_set_ui(fp_, 2);
  fq_ctx_init(Cc->fqctx, fp_, 1, "F");

  Cc->D = NULL;

  Cc->elementary_boundaries = NULL;
  Cc->elementary_coboundaries = NULL;

  Cc->cocycles = NULL;
  Cc->intervals = NULL;

  Cc->cocycles_harmonic = NULL;

  Cc->algorithm = COHOMOLOGY_INTEGRAL;
  Cc->l_prime = l;

  Cc->elementary_coboundaries_z = NULL;
  Cc->cocycles_z = NULL;
  Cc->intervals_torsion = NULL;

  Cc->local_system = NULL;

  Cc->sheaf_rank = 0;
  Cc->sheaf_transport = NULL;

  return 0;
}

int32_t complex_cochain_free(
    complex_cochain_t * Cc)
{
  complex_simplicial_free(&Cc->C);

  fq_ctx_clear(Cc->fqctx);

  int rc;

  Judy_free_2_a(Cc->D);
  Judy_free_2_b(Cc->elementary_boundaries, &Cc->fqctx);
  Judy_free_2_b(Cc->elementary_coboundaries, &Cc->fqctx);
  Judy_free_2_b(Cc->cocycles, &Cc->fqctx);

  //JLFA(rc, Cc->cocycles_vanished);

  JLFA(rc, Cc->intervals);

  JLFA(rc, Cc->cocycles_harmonic);

  Judy_free_2_z(Cc->elementary_coboundaries_z);
  Judy_free_2_z(Cc->cocycles_z);
  JLFA(rc, Cc->intervals_torsion);

  JLFA(rc, Cc->local_system);

  // free sheaf transport matrices
  if (Cc->sheaf_transport != NULL)
  {
    uint64_t * sp;
    uint64_t sk = 0;
    JLF(sp, Cc->sheaf_transport, sk);
    while (sp != NULL)
    {
      if (*sp != 0)
        free((int64_t *) *sp);
      JLN(sp, Cc->sheaf_transport, sk);
    }
    JLFA(rc, Cc->sheaf_transport);
  }

  return 0;
}

int32_t time_series_init(
    time_series_t * S,
    uint32_t d, uint32_t n)
{
  S->d = d;
  S->n = n;
  if ((d * n) > 0)
    { S->x = (double *) malloc((S->d * S->n) * sizeof(double)); }

  return 0;
}

int32_t time_series_free(
    time_series_t * S)
{
  if (S->x != NULL) { free(S->x); }

  return 0;
}
