#ifndef __common__
#define __common__

#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <inttypes.h>

#include <math.h>

#include <pthread.h>
#include "zlog.h"
#include "Judy.h"
//#include "hiredis/hiredis.h"

#include "gmp.h"
#include "mpfr.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/fq.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

#include "defs.h"
#include "enums.h"
#include "fwddecls.h"
#include "structs.h"
#include "ctx.h"

#endif
