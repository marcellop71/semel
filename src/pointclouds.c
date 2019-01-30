#include "common.h"
#include "semel.h"

int32_t point_cloud_init_distance_matrix(
    point_cloud_t * P)
{
  if (P->A == NULL)
  {
    if (P->n >= 2)
    {
      P->nA = (P->n * (P->n - 1)) >> 1;
      P->A = (double *) malloc(P->nA * sizeof(double));
      distance_matrix(
        P->d, P->n, P->p,
        P->manifold, P->distance, P->nb, P->b,
        &P->A, &P->Amax);
    }
  }

  return 0;
}
