#include "common.h"
#include "semel.h"

double vector_dist(
    uint32_t d, double * p, double * q,
    uint32_t manifold, uint32_t distance, uint32_t nb, double * b)
{
  double z = 0.0;

  if (distance == DISTANCE_L2)
  {
    for(uint32_t i = 0; i < d; ++i)
     { z += pow((double) (p[i] - q[i]), 2.0); }
    z = sqrt(z);
    z = round(z * PRECISION) / PRECISION;
  }

  if (distance == DISTANCE_COSINE)
  {
    double zpq, zp, zq;
    zpq = 0.0; zp = 0.0; zq = 0.0;
    for(uint32_t i = 0; i < d; ++i)
    {
      zpq += ((double) p[i]) * ((double) q[i]);
      zp += pow((double) p[i], 2.0);
      zq += pow((double) q[i], 2.0);
    }
    zp = sqrt(zp); zq = sqrt(zq);
    z = zpq / (zp * zq);
    z = round(z * PRECISION) / PRECISION;
  }

  if (distance == DISTANCE_SURFACE_TORUS)
  {
    double tmp_x_1 = fabs(p[0] - q[0]);
    double tmp_y_1 = fabs(p[1] - q[1]);
    double tmp_x_2 = b[0] - tmp_x_1;
    double tmp_y_2 = b[1] - tmp_y_1;
    double tmp_x = (tmp_x_1 < tmp_x_2) ? tmp_x_1 : tmp_x_2;
    double tmp_y = (tmp_y_1 < tmp_y_2) ? tmp_y_1 : tmp_y_2;
    z = pow(tmp_x, 2.0) + pow(tmp_y, 2.0);
    z = sqrt(z);
    z = round(z * PRECISION) / PRECISION;
  }

  return z;
}

double vector_norm_squared(
    uint32_t d, double * p)
{
  double z = 0;
  for(uint32_t i = 0; i < d; ++i)
   { z += pow(p[i], 2.0); }
  return z;
}

double vector_dot(
    uint32_t d, double * p, double * q)
{
  double z = 0;
  for(uint32_t i = 0; i < d; ++i)
   { z += p[i] * q[i]; }
  return z;
}

int32_t vector_add(
    uint32_t d, double * p, double * q, double * x)
{
  for(uint32_t i = 0; i < d; ++i)
   { x[i] = p[i] + q[i]; }
  return 0;
}

int32_t vector_sub(
    uint32_t d, double * p, double * q, double * x)
{
  for(uint32_t i = 0; i < d; ++i)
   { x[i] = p[i] - q[i]; }
  return 0;
}

int32_t vector_cross(
    uint32_t d, double * p, double * q, double * cp)
{
  if (d == 3)
  {
    cp[0] = p[1] * q[2] - p[2] * q[1];
    cp[1] = p[2] * q[0] - p[0] * q[2];
    cp[2] = p[0] * q[1] - p[1] * q[0];
  }
  return 0;
}

typedef struct {
  uint32_t d;
  double * p;
  uint32_t manifold, distance, nb;
  double * b;
  double * A;
  uint32_t i_start, i_end;
} dist_thread_arg_t;

static void * distance_matrix_worker(void * arg)
{
  dist_thread_arg_t * a = (dist_thread_arg_t *) arg;
  for(uint32_t i = a->i_start; i < a->i_end; ++i)
  {
    uint32_t k = (i * (i - 1)) >> 1;
    for(uint32_t j = 0; j < i; ++j)
    {
      a->A[k++] = vector_dist(
        a->d, &a->p[a->d * i], &a->p[a->d * j],
        a->manifold, a->distance, a->nb, a->b);
    }
  }
  return NULL;
}

int32_t distance_matrix(
    uint32_t d, uint32_t n, double * p,
    uint32_t manifold, uint32_t distance, uint32_t nb, double * b,
    double ** A, double * Amax)
{
  long nthreads = sysconf(_SC_NPROCESSORS_ONLN);
  if (nthreads < 1) { nthreads = 1; }
  if ((uint32_t) nthreads > n) { nthreads = n; }

  pthread_t * threads = (pthread_t *) malloc(nthreads * sizeof(pthread_t));
  dist_thread_arg_t * args = (dist_thread_arg_t *) malloc(nthreads * sizeof(dist_thread_arg_t));

  uint32_t rows_per_thread = n / nthreads;
  uint32_t remainder = n % nthreads;
  uint32_t row = 0;

  for(long t = 0; t < nthreads; ++t)
  {
    args[t].d = d;
    args[t].p = p;
    args[t].manifold = manifold;
    args[t].distance = distance;
    args[t].nb = nb;
    args[t].b = b;
    args[t].A = *A;
    args[t].i_start = row;
    args[t].i_end = row + rows_per_thread + (t < (long) remainder ? 1 : 0);
    row = args[t].i_end;
    pthread_create(&threads[t], NULL, distance_matrix_worker, &args[t]);
  }

  for(long t = 0; t < nthreads; ++t)
    { pthread_join(threads[t], NULL); }

  double Amax_ = 0.0;
  uint32_t nA = (n * (n - 1)) >> 1;
  for(uint32_t k = 0; k < nA; ++k)
  {
    if ((*A)[k] > Amax_) { Amax_ = (*A)[k]; }
  }
  *Amax = Amax_;

  free(threads);
  free(args);
  return 0;
}

uint32_t num_points_from_card_of_distance_matrix(
    uint32_t a)
{
  uint32_t n = (1 + (uint32_t) sqrt(1 + (8 * a))) / 2;
  return n;
}

double get_dist_G(
    uint32_t i, uint32_t j, complex_simplicial_t * C)
{
  if (i == j) { return 0.0; }

  uint64_t * p;
  uint64_t ij = (i < j) ? (((uint64_t) i) << 32) + j : (((uint64_t) j) << 32) + i;
  JLG(p, C->G, ij);
  if (p != NULL)
  {
    double d;
    memcpy(&d, p, sizeof(double));
    return d;
  }

  return 1E+100;
}

double get_dist(
    uint32_t i, uint32_t j, point_cloud_t * P)
{
  if (i > j)
    { return get_dist(j, i, P); }
  if (i == j)
    { return 0.0; }
  if (i < j)
  {
    if (P->A != NULL)
      { return P->A[((j * (j - 1)) >> 1) + i]; }
    else
      { return vector_dist(
          P->d, &P->p[P->d * i], &P->p[P->d * j],
          P->manifold, P->distance, P->nb, P->b); }
  }
  return 0.0;
}

double max_pairwise_dist_G(
    uint32_t d, uint32_t * pi, complex_simplicial_t * C)
{
  double dmax = 0.0;
  uint32_t n = d + 1;
  double d_;
  for(uint32_t i = 0; i < n; ++i)
  {
    for(uint32_t j = i+1; j < n; ++j)
    {
      d_ = get_dist_G(pi[i], pi[j], C);
      if (d_ > dmax)
        { dmax = d_; }
    }
  }
  return dmax;
}

double max_pairwise_dist(
    uint32_t d, uint32_t * pi, point_cloud_t * P)
{
  double dmax = 0.0;
  uint32_t n = d + 1;
  double d_;
  for(uint32_t i = 0; i < n; ++i)
  {
    for(uint32_t j = i+1; j < n; ++j)
    {
      d_ = get_dist(pi[i], pi[j], P);
      if (d_ > dmax)
        { dmax = d_; }
    }
  }
  return dmax;
}

double circumradius_jung(
    uint32_t d, uint32_t * pi, point_cloud_t * P)
{
  double dmax = max_pairwise_dist(d, pi, P);
  double r = sqrt(((double) d) / (2.0 * (((double) d) + 1.0))) * dmax;
  return r;
}

double circumradius(
    uint32_t d, uint32_t * pi, point_cloud_t * P)
{
  double r;

  if (d <= 0)
    { r = 0.0; return r; }

  if (d == 1)
  {
    double a = get_dist(pi[0], pi[1], P);
    r = 0.5 * a;
    return r;
  }

  if (d == 2)
  {
    double a = get_dist(pi[0], pi[1], P);
    double b = get_dist(pi[1], pi[2], P);
    double c = get_dist(pi[0], pi[2], P);
    double s = 0.5 * (a + b + c);
    double tmp = s * (a + b - s) * (a + c - s) * (b + c - s);
    r = (tmp > 0) ? (a * b * c) / (4.0 * sqrt(tmp)) : 0.0;
    return r;
  }

  uint32_t n = d + 1;
  gsl_permutation * p;
  gsl_matrix * X;
  int signum;
  double d_;

  X = gsl_matrix_alloc(n + 1, n + 1);
  gsl_matrix_set(X, 0, 0, 0.0);
  for(uint32_t i = 0; i < n; ++i)
  {
    gsl_matrix_set(X, i + 1, 0, 1.0);
    gsl_matrix_set(X, 0, i + 1, 1.0);
    gsl_matrix_set(X, i + 1, i + 1, 0.0);
    for(uint32_t j = i + 1; j < n; ++j)
    {
      d_ = get_dist(pi[i], pi[j], P); d_ = d_ * d_;
      gsl_matrix_set(X, i + 1, j + 1, d_);
      gsl_matrix_set(X, j + 1, i + 1, d_);
    }
  }

  p = gsl_permutation_alloc(n + 1); signum = 0;
  gsl_linalg_LU_decomp(X, p, &signum);
  double det1 = gsl_linalg_LU_det(X, signum);
  gsl_permutation_free(p);
  gsl_matrix_free(X);

  X = gsl_matrix_alloc(n, n);
  for(uint32_t i = 0; i < n; ++i)
  {
    gsl_matrix_set(X, i, i, 0.0);
    for(uint32_t j = i+1; j < n; ++j)
    {
      d_ = get_dist(pi[i], pi[j], P); d_ = d_ * d_;
      gsl_matrix_set(X, i, j, d_);
      gsl_matrix_set(X, j, i, d_);
    }
  }

  p = gsl_permutation_alloc(n); signum = 0;
  gsl_linalg_LU_decomp(X, p, &signum);
  double det2 = gsl_linalg_LU_det(X, signum);
  gsl_permutation_free(p);
  gsl_matrix_free(X);

  r = (det1 > 0) ? sqrt(- 0.50 * (det2 / det1)) : 0.0;

  double rmin = circumradius_jung(d, pi, P);
  if (rmin < r)
    { r = rmin; }

  return r;
}

int32_t circumsphere(
    uint32_t d, uint32_t n, double * p,
    double * cc, double * cr)
{
  *cr = 0.0;
  uint32_t i, j, k;

  if (n <= (d + 1))
  {
    if (n == 1)
      { memcpy(cc, p, d * sizeof(double)); *cr = 0.0; }
    if (n == 2)
    {
      for(k = 0; k < d; ++k)
        { cc[k] = 0.50 * p[k + (d * 0)] + 0.50 * p[k + (d * 1)]; }
      *cr = 0.50 * vector_dist(d, &p[d * 0], &p[d * 1],
        MANIFOLD_FLAT_UNBOUNDED, DISTANCE_L2, 0, NULL);
    }
    if (n == 3)
    {
      if (d == 2)
      {
        double * d2 = (double *) malloc(n * sizeof(double));
        for(i = 0; i < n; ++i)
        {
          double tmp = 0.0;
          for(k = 0; k < d; ++k)
            { tmp += pow(p[k + (d * i)], 2.0); }
          d2[i] = tmp;
        }

        gsl_matrix * Dx, * Dy, * a, * c;

        Dx = gsl_matrix_alloc(n, d + 1);
        Dy = gsl_matrix_alloc(n, d + 1);
        a = gsl_matrix_alloc(n, d + 1);
        c = gsl_matrix_alloc(n, d + 1);

        for(i = 0; i < n; ++i)
        {
          for(j = 0; j < d; ++j)
            { gsl_matrix_set(a, i, j, p[j + (d * i)]); }
          gsl_matrix_set(a, i, d, 1.0);

          gsl_matrix_set(Dx, i, 0, d2[i]);
          gsl_matrix_set(Dx, i, 1, p[1 + (d * i)]);
          gsl_matrix_set(Dx, i, d, 1.0);

          gsl_matrix_set(Dy, i, 0, d2[i]);
          gsl_matrix_set(Dy, i, 1, p[0 + (d * i)]);
          gsl_matrix_set(Dy, i, d, 1.0);

          gsl_matrix_set(c, i, 0, d2[i]);
          for(j = 0; j < d; ++j)
            { gsl_matrix_set(c, i, j+1, p[j + (d * i)]); }
        }

        gsl_permutation * perm = gsl_permutation_alloc(d + 1);
        int signum = 0;

        double Dx_, Dy_, a_, c_;
        gsl_linalg_LU_decomp(Dx, perm, &signum); Dx_ = gsl_linalg_LU_det(Dx, signum);
        gsl_linalg_LU_decomp(Dy, perm, &signum); Dy_ = gsl_linalg_LU_det(Dy, signum);
        gsl_linalg_LU_decomp(a, perm, &signum); a_ = gsl_linalg_LU_det(a, signum);
        gsl_linalg_LU_decomp(c, perm, &signum); c_ = gsl_linalg_LU_det(c, signum);

        cc[0] = (a_ > 0) ? Dx_ / (2.0 * a_) : 0.0;
        cc[1] = (a_ > 0) ? - Dy_ / (2.0 * a_) : 0.0;

        *cr = (fabs(a_) > 0) ? sqrt((Dx_ * Dx_) + (Dy_ * Dy_) - 4.0 * (a_ * c_)) / (2.0 * fabs(a_)) : 0.0;

        gsl_permutation_free(perm);
        gsl_matrix_free(Dx); gsl_matrix_free(Dy);
        gsl_matrix_free(a); gsl_matrix_free(c);
        free(d2);
      }
      if (d == 3)
      {
        double * d01 = (double *) malloc(d * sizeof(double));
        double * d02 = (double *) malloc(d * sizeof(double));
        double * d12 = (double *) malloc(d * sizeof(double));
        vector_sub(d, &p[(d * 0)], &p[(d * 1)], d01);
        vector_sub(d, &p[(d * 0)], &p[(d * 2)], d02);
        vector_sub(d, &p[(d * 1)], &p[(d * 2)], d12);

        double * cp0 = (double *) malloc(d * sizeof(double));
        vector_cross(d, d01, d12, cp0);
        double y = vector_norm_squared(d, cp0);
        free(cp0);

        double d01_ = vector_dist(d, &p[d * 0], &p[d * 1],
          MANIFOLD_FLAT_UNBOUNDED, DISTANCE_L2, 0, NULL);
        double d02_ = vector_dist(d, &p[d * 0], &p[d * 2],
          MANIFOLD_FLAT_UNBOUNDED, DISTANCE_L2, 0, NULL);
        double d12_ = vector_dist(d, &p[d * 1], &p[d * 2],
          MANIFOLD_FLAT_UNBOUNDED, DISTANCE_L2, 0, NULL);

        if (y > 0)
        {
          double a0 = (y > 0) ? (pow(d12_, 2.0) * vector_dot(d, d01, d02)) / (2.0 * y) : 0.0;
          double a1 = (y > 0) ? (pow(d02_, 2.0) * vector_dot(d, d01, d12)) / (2.0 * y) : 0.0;
          double a2 = (y > 0) ? (pow(d01_, 2.0) * vector_dot(d, d02, d12)) / (2.0 * y) : 0.0;

          for(k = 0; k < d; ++k)
            { cc[k] = a0 * p[k + (d * 0)] - a1 * p[k + (d * 1)] + a2 * p[k + (d * 2)]; }

          *cr = (y > 0) ? (d01_ * d02_ * d12_) / (2.0 * sqrt(y)) : 0.0;
        }
        else
        {
          double dmax_;

          dmax_ = d01_;
          for(k = 0; k < d; ++k)
            { cc[k] = 0.50 * p[k + (d * 0)] + 0.50 * p[k + (d * 1)]; }

          if (d02_ > dmax_)
          {
            dmax_ = d02_;
            for(k = 0; k < d; ++k)
              { cc[k] = 0.50 * p[k + (d * 0)] + 0.50 * p[k + (d * 2)]; }
          }

          if (d12_ > dmax_)
          {
            dmax_ = d12_;
            for(k = 0; k < d; ++k)
              { cc[k] = 0.50 * p[k + (d * 1)] + 0.50 * p[k + (d * 2)]; }
          }

          *cr = 0.50 * dmax_;
        }

        free(d01); free(d02); free(d12);
      }
    }
    if (n == 4)
    {
      if (d == 3)
      {
        double * d2 = (double *) malloc(n * sizeof(double));
        for(i = 0; i < n; ++i)
        {
          double tmp = 0.0;
          for(k = 0; k < d; ++k)
            { tmp += pow(p[k + (d * i)], 2.0); }
          d2[i] = tmp;
        }

        gsl_matrix * Dx, * Dy, * Dz, * a, * c;

        Dx = gsl_matrix_alloc(n, d + 1);
        Dy = gsl_matrix_alloc(n, d + 1);
        Dz = gsl_matrix_alloc(n, d + 1);
        a = gsl_matrix_alloc(n, d + 1);
        c = gsl_matrix_alloc(n, d + 1);

        for(i = 0; i < n; ++i)
        {
          for(j = 0; j < d; ++j)
            { gsl_matrix_set(a, i, j, p[j + (d * i)]); }
          gsl_matrix_set(a, i, d, 1.0);

          gsl_matrix_set(Dx, i, 0, d2[i]);
          gsl_matrix_set(Dx, i, 1, p[1 + (d * i)]);
          gsl_matrix_set(Dx, i, 2, p[2 + (d * i)]);
          gsl_matrix_set(Dx, i, d, 1.0);

          gsl_matrix_set(Dy, i, 0, d2[i]);
          gsl_matrix_set(Dy, i, 1, p[0 + (d * i)]);
          gsl_matrix_set(Dy, i, 2, p[2 + (d * i)]);
          gsl_matrix_set(Dy, i, d, 1.0);

          gsl_matrix_set(Dz, i, 0, d2[i]);
          gsl_matrix_set(Dz, i, 1, p[0 + (d * i)]);
          gsl_matrix_set(Dz, i, 2, p[1 + (d * i)]);
          gsl_matrix_set(Dz, i, d, 1.0);

          gsl_matrix_set(c, i, 0, d2[i]);
          for(j = 0; j < d; ++j)
            { gsl_matrix_set(c, i, j+1, p[j + (d * i)]); }
        }

        gsl_permutation * perm = gsl_permutation_alloc(d + 1);
        int signum = 0;

        double Dx_, Dy_, Dz_, a_, c_;
        gsl_linalg_LU_decomp(Dx, perm, &signum); Dx_ = gsl_linalg_LU_det(Dx, signum);
        gsl_linalg_LU_decomp(Dy, perm, &signum); Dy_ =  gsl_linalg_LU_det(Dy, signum);
        gsl_linalg_LU_decomp(Dz, perm, &signum); Dz_ = gsl_linalg_LU_det(Dz, signum);
        gsl_linalg_LU_decomp(a, perm, &signum); a_ = gsl_linalg_LU_det(a, signum);
        gsl_linalg_LU_decomp(c, perm, &signum); c_ = gsl_linalg_LU_det(c, signum);

        cc[0] = (a_ > 0) ? Dx_ / (2.0 * a_) : 0.0;
        cc[1] = (a_ > 0) ? - Dy_ / (2.0 * a_) : 0.0;
        cc[2] = (a_ > 0) ? Dz_ / (2.0 * a_) : 0.0;

        *cr = (fabs(a_) > 0) ? sqrt((Dx_ * Dx_) + (Dy_ * Dy_) + (Dz_ * Dz_) - 4.0 * (a_ * c_)) / (2.0 * fabs(a_)) : 0.0;

        gsl_permutation_free(perm);
        gsl_matrix_free(Dx); gsl_matrix_free(Dy); gsl_matrix_free(Dz);
        gsl_matrix_free(a); gsl_matrix_free(c);
        free(d2);
      }
    }
  }

  return 0;
}

int32_t miniball(
    uint32_t d,
    uint32_t n, double * p,
    void ** tau, void ** v,
    double * cc, double * cr)
{
  int32_t rc;
  uint64_t ntau;
  J1C(ntau, *tau, 0, -1);
  if (ntau == 0)
  {
    uint32_t np_;
    J1C(np_, *v, 0, -1);
    if (np_ == 0)
      { *cr = 0.0; }
    else
    {
      double * p_ = (double *) malloc((d * np_) * sizeof(double));

      uint64_t i, j;
      i = 0; J1F(rc, *v, i); j = 0;
      while (rc == 1)
      {
        memcpy(&p_[d * j], &p[d * i], d * sizeof(double));
        J1N(rc, *v, i); j++;
      }

      circumsphere(d, np_, p_, cc, cr);
      free(p_);
    }
  }
  else
  {
    uint64_t ip0;
    int32_t k = 1 + (rand() % ntau);
    J1BC(rc, *tau, k, ip0);
    J1U(rc, *tau, ip0);

    miniball(d, n, p, tau, v, cc, cr);
    if ((*cr == 0.0) || (vector_dist(d, &p[d * ip0], cc,
      MANIFOLD_FLAT_UNBOUNDED, DISTANCE_L2, 0, NULL) > *cr))
    {
      J1S(rc, *v, ip0);
      miniball(d, n, p, tau, v, cc, cr);
    }
  }

  return 0;
}

int32_t smallest_ball(
    uint32_t d, uint32_t n, double * p,
    double * cc, double * cr)
{
  int32_t rc;
  void * tau = NULL;
  void * v = NULL;
  for(uint32_t k = 0; k < n; ++k)
    { J1S(rc, tau, k); }

  miniball(
    d, n, p, &tau, &v, cc, cr);

  return 0;
}
