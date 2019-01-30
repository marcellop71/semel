#include "common.h"
#include "semel.h"

int32_t Laplace_eigen(
    complex_cochain_t * Cc, uint32_t n,
    uint32_t ** na, double *** a, double *** v)
{
  uint32_t * mD, * nD;
  gsl_matrix ** D;

  mD = (uint32_t *) malloc(n * sizeof(uint32_t));
  nD = (uint32_t *) malloc(n * sizeof(uint32_t));
  D = (gsl_matrix **) malloc(n * sizeof(gsl_matrix *));

  *na = (uint32_t *) malloc((n-1) * sizeof(uint32_t));
  *a = (double **) malloc((n-1) * sizeof(double *));
  *v = (double **) malloc((n-1) * sizeof(double *));

  for(uint32_t k = 0; k < n; k++)
    { coboundary_operator_as_gsl_matrix(
        k, Cc, &mD[k], &nD[k], &D[k]); }

  for(uint32_t k = 0; k < (n-1); k++)
  {
    uint32_t nL = mD[k];
    double alpha = 1.0, beta = 0.0;

    gsl_matrix * D0, * D1, * L;

    D0 = gsl_matrix_alloc(nL, nL);
    D1 = gsl_matrix_alloc(nL, nL);
    gsl_blas_dgemm(
      CblasNoTrans, CblasTrans, alpha, D[k], D[k], beta, D0);
    gsl_blas_dgemm(
      CblasTrans, CblasNoTrans, alpha, D[k+1], D[k+1], beta, D1);
    gsl_matrix_add(D0, D1);

    L = gsl_matrix_alloc(nL, nL);
    gsl_matrix_memcpy(L, D0);

    gsl_eigen_symm_workspace * w =
      //gsl_eigen_symm_alloc(nL);
      gsl_eigen_symmv_alloc(nL);

    gsl_vector * eig = gsl_vector_alloc(nL);
    gsl_matrix * eiv = gsl_matrix_alloc(nL, nL);
    //gsl_eigen_symm(L, eig, w);
    gsl_eigen_symmv(L, eig, eiv, w);
    //gsl_sort_vector(eig);
    gsl_eigen_symmv_sort(eig, eiv, GSL_EIGEN_SORT_VAL_ASC);

    (*na)[k] = nL;
    (*a)[k] = (double *) malloc(nL * sizeof(double));
    (*v)[k] = (double *) malloc((nL * nL) * sizeof(double));
    for(uint32_t j = 0; j < nL; j++)
    {
      (*a)[k][j] = gsl_vector_get(eig, j);
      gsl_vector_view colj = gsl_matrix_column(eiv, j);;
      for(uint32_t i = 0; i < nL; i++)
        { (*v)[k][j + (nL * i)] = gsl_vector_get(&colj, i); }
    }

    gsl_vector_free(eig);
    gsl_matrix_free(eiv);
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(D0);
    gsl_matrix_free(D1);
    gsl_matrix_free(L);
  }

  for(uint32_t k = 0; k < n; k++)
    { gsl_matrix_free(D[k]); }
  free(mD); free(nD); free(D);

  return 0;
}

// ── sparse Hodge Laplacian via Lanczos ───────────────────────────────

typedef struct {
  uint32_t nrow, ncol, nnz;
  double * val;
  uint32_t * col_idx;
  uint32_t * row_ptr;
} csr_t;

static void csr_free(csr_t * A)
{
  free(A->val);
  free(A->col_idx);
  free(A->row_ptr);
  free(A);
}

// Build CSR of D_k (coboundary operator from k-cochains to (k+1)-cochains)
// directly from Cc->elementary_coboundaries and the simplicial index maps.
// D_k has mD = scc[k+1] rows, nD = scc[k] columns.
static csr_t * csr_from_coboundary(uint32_t k, complex_cochain_t * Cc)
{
  uint32_t mD = Cc->C.scc[k+1];
  uint32_t nD = Cc->C.scc[k];

  // first pass: count nonzeros per row
  uint32_t * row_count = (uint32_t *) calloc(mD > 0 ? mD : 1, sizeof(uint32_t));

  fmpz_t f;
  fmpz_init(f);
  fq_ctx_order(f, Cc->fqctx);
  uint32_t o = fmpz_get_ui(f);
  fmpz_clear(f);

  uint64_t v;
  uint64_t * p, * q, * p_, * q_;

  v = 0; JLF(p, Cc->C.scd[k], v);
  while (p != NULL)
  {
    JLG(q, Cc->elementary_coboundaries, *p);
    if (q != NULL)
    {
      uint64_t w = 0; JLF(p_, (void *) *q, w);
      while (p_ != NULL)
      {
        JLG(q_, Cc->C.dimension, w);
        if (q_ != NULL && *q_ == (k + 1))
        {
          uint64_t * ri;
          JLG(ri, Cc->C.sci[*q_], w);
          if (ri != NULL && *ri < mD)
            row_count[*ri]++;
        }
        JLN(p_, (void *) *q, w);
      }
    }
    JLN(p, Cc->C.scd[k], v);
  }

  // build row_ptr
  csr_t * A = (csr_t *) malloc(sizeof(csr_t));
  A->nrow = mD; A->ncol = nD;
  A->row_ptr = (uint32_t *) malloc((mD + 1) * sizeof(uint32_t));
  A->row_ptr[0] = 0;
  for (uint32_t i = 0; i < mD; i++)
    A->row_ptr[i+1] = A->row_ptr[i] + row_count[i];
  A->nnz = A->row_ptr[mD];

  A->val = (double *) malloc(A->nnz * sizeof(double));
  A->col_idx = (uint32_t *) malloc(A->nnz * sizeof(uint32_t));

  // reset row_count as insertion cursor
  memset(row_count, 0, mD * sizeof(uint32_t));

  // second pass: fill values
  v = 0; JLF(p, Cc->C.scd[k], v);
  while (p != NULL)
  {
    JLG(q, Cc->elementary_coboundaries, *p);
    if (q != NULL)
    {
      uint64_t w = 0; JLF(p_, (void *) *q, w);
      while (p_ != NULL)
      {
        JLG(q_, Cc->C.dimension, w);
        if (q_ != NULL && *q_ == (k + 1))
        {
          uint64_t * ri;
          JLG(ri, Cc->C.sci[*q_], w);
          if (ri != NULL && *ri < mD)
          {
            char * tmp = fq_get_str_pretty(*((fq_t *) *p_), Cc->fqctx);
            int32_t a = atoi(tmp);
            free(tmp);
            a = (a > (int32_t) (o >> 1)) ? (a - o) : a;

            uint32_t row = (uint32_t) *ri;
            uint32_t pos = A->row_ptr[row] + row_count[row];
            if (pos < A->nnz)
            {
              A->val[pos] = (double) a;
              A->col_idx[pos] = (uint32_t) v;
            }
            row_count[row]++;
          }
        }
        JLN(p_, (void *) *q, w);
      }
    }
    JLN(p, Cc->C.scd[k], v);
  }

  free(row_count);
  return A;
}

// Build CSR of A^T from CSR of A
static csr_t * csr_transpose(csr_t * A)
{
  csr_t * AT = (csr_t *) malloc(sizeof(csr_t));
  AT->nrow = A->ncol; AT->ncol = A->nrow; AT->nnz = A->nnz;

  // count nonzeros per row of AT (= per column of A)
  uint32_t * count = (uint32_t *) calloc(AT->nrow, sizeof(uint32_t));
  for (uint32_t i = 0; i < A->nnz; i++)
    count[A->col_idx[i]]++;

  AT->row_ptr = (uint32_t *) malloc((AT->nrow + 1) * sizeof(uint32_t));
  AT->row_ptr[0] = 0;
  for (uint32_t i = 0; i < AT->nrow; i++)
    AT->row_ptr[i+1] = AT->row_ptr[i] + count[i];

  AT->val = (double *) malloc(AT->nnz * sizeof(double));
  AT->col_idx = (uint32_t *) malloc(AT->nnz * sizeof(uint32_t));

  memset(count, 0, AT->nrow * sizeof(uint32_t));
  for (uint32_t r = 0; r < A->nrow; r++)
  {
    for (uint32_t j = A->row_ptr[r]; j < A->row_ptr[r+1]; j++)
    {
      uint32_t c = A->col_idx[j];
      uint32_t pos = AT->row_ptr[c] + count[c];
      AT->val[pos] = A->val[j];
      AT->col_idx[pos] = r;
      count[c]++;
    }
  }

  free(count);
  return AT;
}

// y = A * x  (CSR SpMV)
static void csr_spmv(csr_t * A, double * x, double * y)
{
  for (uint32_t i = 0; i < A->nrow; i++)
  {
    double s = 0.0;
    for (uint32_t j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
      s += A->val[j] * x[A->col_idx[j]];
    y[i] = s;
  }
}

// y = L * x = Dk * DkT * x + Dkp1T * Dkp1 * x
static void hodge_matvec(
    csr_t * Dk, csr_t * DkT,
    csr_t * Dkp1, csr_t * Dkp1T,
    double * x, double * y,
    double * w1, double * w2, uint32_t n)
{
  // down part: w1 = DkT * x (size Dk->ncol), y = Dk * w1 (size n)
  csr_spmv(DkT, x, w1);
  csr_spmv(Dk, w1, y);

  // up part: w2 = Dkp1 * x (size Dkp1->nrow), w1 reuse = Dkp1T * w2 (size n)
  csr_spmv(Dkp1, x, w2);
  csr_spmv(Dkp1T, w2, w1);

  for (uint32_t i = 0; i < n; i++)
    y[i] += w1[i];
}

int32_t Laplace_eigen_sparse(
    complex_cochain_t * Cc, uint32_t n, uint32_t nev,
    uint32_t ** na, double *** a, double *** v)
{
  *na = (uint32_t *) malloc((n-1) * sizeof(uint32_t));
  *a = (double **) malloc((n-1) * sizeof(double *));
  *v = (double **) malloc((n-1) * sizeof(double *));

  for (uint32_t k = 0; k < (n-1); k++)
  {
    uint32_t nL = Cc->C.scc[k+1];

    // clamp nev to matrix size
    uint32_t nev_k = (nev < nL) ? nev : nL;

    if (nL == 0 || nev_k == 0)
    {
      (*na)[k] = 0;
      (*a)[k] = NULL;
      (*v)[k] = NULL;
      continue;
    }

    // if nev_k >= nL, fall back to dense solve for full spectrum
    if (nev_k >= nL)
    {
      // build dense Laplacian and solve with GSL
      csr_t * Dk  = csr_from_coboundary(k, Cc);
      csr_t * Dkp1 = csr_from_coboundary(k+1, Cc);
      csr_t * DkT = csr_transpose(Dk);
      csr_t * Dkp1T = csr_transpose(Dkp1);

      uint32_t ws1 = nL;
      if (Dk->ncol > ws1) ws1 = Dk->ncol;
      double * work1 = (double *) calloc(ws1, sizeof(double));
      double * work2 = (double *) calloc(Dkp1->nrow > 1 ? Dkp1->nrow : 1, sizeof(double));

      gsl_matrix * L = gsl_matrix_alloc(nL, nL);
      double * col = (double *) calloc(nL, sizeof(double));
      double * out = (double *) calloc(nL, sizeof(double));
      for (uint32_t j = 0; j < nL; j++)
      {
        memset(col, 0, nL * sizeof(double));
        col[j] = 1.0;
        hodge_matvec(Dk, DkT, Dkp1, Dkp1T, col, out, work1, work2, nL);
        for (uint32_t i = 0; i < nL; i++)
          gsl_matrix_set(L, i, j, out[i]);
      }
      free(col); free(out); free(work1); free(work2);
      csr_free(Dk); csr_free(DkT); csr_free(Dkp1); csr_free(Dkp1T);

      gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(nL);
      gsl_vector * eig = gsl_vector_alloc(nL);
      gsl_matrix * eiv = gsl_matrix_alloc(nL, nL);
      gsl_eigen_symmv(L, eig, eiv, w);
      gsl_eigen_symmv_sort(eig, eiv, GSL_EIGEN_SORT_VAL_ASC);

      (*na)[k] = nL;
      (*a)[k] = (double *) malloc(nL * sizeof(double));
      (*v)[k] = (double *) malloc((nL * nL) * sizeof(double));
      for (uint32_t j = 0; j < nL; j++)
      {
        (*a)[k][j] = gsl_vector_get(eig, j);
        gsl_vector_view colj = gsl_matrix_column(eiv, j);
        for (uint32_t i = 0; i < nL; i++)
          (*v)[k][j + (nL * i)] = gsl_vector_get(&colj.vector, i);
      }

      gsl_vector_free(eig); gsl_matrix_free(eiv);
      gsl_eigen_symmv_free(w); gsl_matrix_free(L);
      continue;
    }

    // build CSR matrices
    csr_t * Dk   = csr_from_coboundary(k, Cc);
    csr_t * Dkp1 = csr_from_coboundary(k+1, Cc);
    csr_t * DkT   = csr_transpose(Dk);
    csr_t * Dkp1T = csr_transpose(Dkp1);

    // work buffers for hodge_matvec
    uint32_t wsize1 = nL;
    if (Dk->ncol > wsize1) wsize1 = Dk->ncol;
    uint32_t wsize2 = Dkp1->nrow > 1 ? Dkp1->nrow : 1;
    double * work1 = (double *) calloc(wsize1, sizeof(double));
    double * work2 = (double *) calloc(wsize2, sizeof(double));

    // Lanczos iteration with full reorthogonalization
    // to find the nev_k smallest eigenvalues of L
    uint32_t m = nL / 2;
    if (m < nev_k * 20) m = nev_k * 20;
    if (m > nL) m = nL;

    // Lanczos vectors Q (nL x m), stored column-major
    double * Q = (double *) calloc((size_t) nL * m, sizeof(double));
    double * alpha_l = (double *) calloc(m, sizeof(double)); // diagonal
    double * beta_l  = (double *) calloc(m, sizeof(double)); // subdiagonal
    double * w_vec = (double *) malloc(nL * sizeof(double));

    // initial vector: q_0 = (1,1,...,1) / ||...||
    double nrm = sqrt((double) nL);
    for (uint32_t i = 0; i < nL; i++)
      Q[i] = 1.0 / nrm;

    // Lanczos steps
    uint32_t m_actual = m;
    for (uint32_t j = 0; j < m; j++)
    {
      double * qj = &Q[(size_t) j * nL];

      // w = L * q_j
      hodge_matvec(Dk, DkT, Dkp1, Dkp1T, qj, w_vec, work1, work2, nL);

      // alpha_j = q_j . w
      double aj = 0.0;
      for (uint32_t i = 0; i < nL; i++)
        aj += qj[i] * w_vec[i];
      alpha_l[j] = aj;

      // w = w - alpha_j * q_j
      for (uint32_t i = 0; i < nL; i++)
        w_vec[i] -= aj * qj[i];

      // w = w - beta_j * q_{j-1}
      if (j > 0)
      {
        double * qjm1 = &Q[(size_t) (j-1) * nL];
        for (uint32_t i = 0; i < nL; i++)
          w_vec[i] -= beta_l[j] * qjm1[i];
      }

      // full reorthogonalization against all previous vectors
      for (uint32_t p = 0; p <= j; p++)
      {
        double * qp = &Q[(size_t) p * nL];
        double dot = 0.0;
        for (uint32_t i = 0; i < nL; i++)
          dot += w_vec[i] * qp[i];
        for (uint32_t i = 0; i < nL; i++)
          w_vec[i] -= dot * qp[i];
      }

      // beta_{j+1} = ||w||
      double bj1 = 0.0;
      for (uint32_t i = 0; i < nL; i++)
        bj1 += w_vec[i] * w_vec[i];
      bj1 = sqrt(bj1);

      if (j + 1 < m)
      {
        if (bj1 < 1e-14)
        {
          m_actual = j + 1;
          break;
        }
        beta_l[j+1] = bj1;
        double * qj1 = &Q[(size_t) (j+1) * nL];
        for (uint32_t i = 0; i < nL; i++)
          qj1[i] = w_vec[i] / bj1;
      }
    }

    free(w_vec); free(work1); free(work2);

    // eigendecompose the m_actual x m_actual tridiagonal matrix T
    // using GSL
    gsl_matrix * T = gsl_matrix_calloc(m_actual, m_actual);
    for (uint32_t i = 0; i < m_actual; i++)
    {
      gsl_matrix_set(T, i, i, alpha_l[i]);
      if (i + 1 < m_actual)
      {
        gsl_matrix_set(T, i, i+1, beta_l[i+1]);
        gsl_matrix_set(T, i+1, i, beta_l[i+1]);
      }
    }

    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(m_actual);
    gsl_vector * eig = gsl_vector_alloc(m_actual);
    gsl_matrix * eiv = gsl_matrix_alloc(m_actual, m_actual);
    gsl_eigen_symmv(T, eig, eiv, w);
    gsl_eigen_symmv_sort(eig, eiv, GSL_EIGEN_SORT_VAL_ASC);

    // take the nev_k smallest eigenvalues (Ritz values)
    uint32_t nout = (nev_k < m_actual) ? nev_k : m_actual;
    (*na)[k] = nout;
    (*a)[k] = (double *) malloc(nout * sizeof(double));
    (*v)[k] = (double *) calloc(nout, sizeof(double));

    for (uint32_t j = 0; j < nout; j++)
      (*a)[k][j] = gsl_vector_get(eig, j);

    gsl_vector_free(eig); gsl_matrix_free(eiv);
    gsl_eigen_symmv_free(w); gsl_matrix_free(T);
    free(Q); free(alpha_l); free(beta_l);
    csr_free(Dk); csr_free(DkT); csr_free(Dkp1); csr_free(Dkp1T);
  }

  return 0;
}
