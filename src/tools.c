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

int32_t get_status(
    semel_ctx_log_t * log)
{
  FILE * fs;
  char * buf;
  uint64_t n;
  fs = fopen("/sys/fs/cgroup/memory/memory.usage_in_bytes", "r");
  if (fs == NULL)
    { return -1; }
  fseek(fs, 0L, SEEK_END);
  n = ftell(fs);
  fseek(fs, 0L, SEEK_SET);
  buf = (char *) calloc(n, sizeof(char));
  fread(buf, sizeof(char), n, fs);
  fclose(fs);
  zlog_info(log->data, "%s", buf);
  free(buf);
  return 0;
}

static inline uint64_t swapbits(uint64_t p, uint64_t m, int k) {
  uint64_t q = ((p >> k)^p) & m;
  return p ^ q ^ (q << k);
}

uint64_t bitreverse(
    uint64_t n)
{
  static const uint64_t m1 = ((1UL << 63) - 1) / (1 + (1 << 1) + (1 << 2));
  static const uint64_t m2 = ((1UL << 63) - 1) / (1 + (1 << 3) + (1 << 6));
  static const uint64_t m3 = ((1UL << 9) - 1) + (((1UL << 9)-1)<<36);
  static const uint64_t m4 = (1UL << 27) - 1;
  n = swapbits(n, m1, 2);
  n = swapbits(n, m2, 6);
  n = swapbits(n, m3, 18);
  n = swapbits(n, m4, 36);
  n = (n >> 63) | (n << 1);
  return n;
}

// n + ((-n) & 63) rounds up n to the next multiple of 64
static inline uint32_t roundup64(
    uint32_t n)
{
  return n + ((-n) & 63);
}

static inline uint32_t numblock64(
    uint32_t n)
{
  return roundup64(n) >> 6;
}

// (n * (n + 1)) >> 1 is the number boxes in an (non-strict) upper triangular matrix n x n
static inline uint32_t triangular(
    uint32_t n)
{
  return (n * (n + 1)) >> 1;
}

// (n * (n + 1)) >> 1 is the number boxes in a strict upper triangular matrix n x n
//static inline uint32_t triangular_strict(
//    uint32_t n)
//{
//  return (n * (n - 1)) >> 1;
//}

uint32_t low(
    uint32_t nx, uint64_t * x)
{
  uint32_t z = 0;

  uint64_t tmp = 0;
  int32_t i = nx-1;
  // x & (-x) isolate the rightmost 1-bit of x, producing 0 if none
  while ((i >= 0) && ((tmp = x[i] & (-x[i])) == 0))
    { --i; }
  if (i >= 0)
  {
    while (tmp > 0)
      { tmp <<= 1; z++; }
    z += 64 * i;
  }

  return z;
}

int32_t identity(
    uint32_t n, uint64_t * V)
{
  uint32_t n_ = numblock64(n);
  uint32_t k;
  for(k = 0; k < (n_ - 1); ++k)
  {
    for(uint32_t i = 0; i < 64; ++i)
      { V[triangular(k) + i] = (1UL << i); }
  }
  for(uint32_t i = 0; i < (n % 64); ++i)
    { V[triangular(k) + i] = (1UL << i); }

  return 0;
}

// uint64_t * A = (uint64_t *) malloc((n * n) * sizeof(uint64_t))
// uint64_t * B = (uint64_t *) malloc((n_ * n) * sizeof(uint64_t))
int32_t dense_compress(
    uint32_t n, uint64_t * A, uint64_t * B)
{
  uint32_t n_ = numblock64(n);
  for(uint32_t j = 0; j < n; ++j)
  {
    for(uint32_t i = 0; i < n_; ++i)
    {
      uint64_t x = 0;
      for(uint32_t k = 0; k < 64; ++k)
      {
        if (A[(i + 64 - 1 - k) + (j * n)] == 1)
          { x |= (1UL << k); }
      }
      B[i + (j * n_)] = x;
    }
  }

  return 0;
}

// uint64_t * A = (uint64_t *) malloc((n_ * n) * sizeof(uint64_t))
// uint64_t * B = (uint64_t *) malloc((n * n) * sizeof(uint64_t))
int32_t dense_decompress(
    uint32_t n, uint64_t * A, uint64_t * B)
{
  uint32_t n_ = numblock64(n);
  for(uint32_t j = 0; j < n; ++j)
  {
    for(uint32_t i = 0; i < n_; ++i)
    {
      for(uint32_t k = 0; k < 64; ++k)
        { B[(i + 64 - 1 - k) + (j * n_)] = (A[i + (j * n)] & (1UL << k)) >> k; }
    }
  }

  return 0;
}

int32_t dense_anti_transpose(
    uint32_t n, uint64_t * A, uint64_t * B)
{
  for(uint32_t j = 0; j < n; ++j)
  {
    for(uint32_t i = 0; i < n; ++i)
      { B[i + (j * n)] = A[(n - 1 - i) + ((n - 1 - j) * n)]; }
  }

  return 0;
}

int32_t dense_compressed_to_upper_triangle(
    uint32_t n, uint64_t * A, uint64_t * B)
{
  uint32_t n_ = numblock64(n);
  uint32_t c = 0;
  for(uint32_t j = 0; j < n; ++j)
  {
    uint32_t nj = numblock64(j);
    for(uint32_t i = 0; i < nj; ++i)
      { B[c] = A[i + (j * n_)]; c++; }
  }

  return 0;
}
