/* Dimension-adaptive sparse-grid cubature based on nested Clenshaw-Curtis
   rules, by Steven G. Johnson.  */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#include "cubature.h"

/* error return codes */
#define SUCCESS 0
#define FAILURE 1

/***************************************************************************/
/* CLENSHAW-CURTIS WEIGHTS */

/* A level-j Clenshaw-Curtis rule corresponds to 2*n+1=1+2^j points for j > 0,
   where the points for an interval [-1,1] are x = +/- cos(i*pi/(2*n)) for 
   i = 0..n-1 and one point at x=0.  We define j=0 to be the 1-point rule
   with a point at x=0, and j=-1 to be the 0-point rule that returns 0.

   We allow j to go from 0 to 63.

   For sparse grids (and adaptation in general), we instead store the
   *difference* (CCD) rules that compute the difference between the j
   and j-1 CC rules.

   For convenience, we store the weights and points in order of
   increasing j, so that a point shared with j-1 comes before a point
   only for levels >= j.  Also, we only store the nonnegative x points,
   since the negative x points have the same weights and are mirror images.

   So, x[0] is 0, x[1] is cos(0)=1, x[2] is cos(pi/4), 
       x[3:4] are cos([1,3]/8), x[5:8] are cos([1,3,5,7]/16), etcetera.

   And w[j][i] are the level-j CCD weights corresponding to x[i].
*/

typedef struct ccd_rules_s {
     int jmax; /* j <= jmax < 64 has been allocated */
     double *x; /* points x[i] in increasing j order; see above */
     double *w[64]; /* weights w[j][i] for x[i] of level-j CCD rule */
} ccd_rules;

/* Compute the CC weights via FFTW's DCT, for 2*n+1 quadrature points.
   The weights are for integration of functions on [-1,1].  w[i] is the
   weight of the points x = +/- cos(i*pi/(2*n)) for 0 <= i < n, and w[n]
   is the weight of the point x=0. */
static void clencurt_weights(int n, double *w)
{
     fftw_plan p;
     int j;
     double scale = 1.0 / n;

     if (n == 0) { w[0] = 2; return; } /* trivial 1-point rule */

/*
   The general principle is this: in Fejer and Clenshaw-Curtis quadrature,
   we take our function f(x) evaluated at f(cos(theta)), discretize
   it in theta to a vector f of points, compute the DCT, then multiply
   by some coefficients c (the integrals of cos(theta) * sin(theta)).
   In matrix form, given the DCT matrix D, this is:

             c' * D * f = (D' * c)' * f = w' * f

   where w is the vector of weights for each function value.  It
   is obviously much nicer to precompute w if we are going to be
   integrating many functions.   Since w = D' * c, and the transpose
   D' of a DCT is another DCT (e.g. the transpose of DCT-I is a DCT-I,
   modulo some rescaling in FFTW's normalization), to compute w is just
   another DCT.

   There is an additional wrinkle, because the odd-frequency DCT components
   of f integrate to zero, so every other entry in c is zero.  We can
   take advantage of this in computing w, because we can essentially do
   a radix-2 step in the DCT where one of the two subtransforms is zero.
   Therefore, for 2n+1 inputs, we only need to do a DCT of size n+1, and
   the weights w are a nice symmetric function.
*/
     p = fftw_plan_r2r_1d(n+1, w, w, FFTW_REDFT00, FFTW_ESTIMATE);
     for (j = 0; j <= n; ++j) w[j] = scale / (1 - 4*j*j);
     fftw_execute(p);
     w[0] *= 0.5;
     fftw_destroy_plan(p);
}

static void init_ccd_rules(ccd_rules *ccd)
{
     int j;
     ccd->jmax = -1;
     ccd->x = NULL;
     for (j = 0; j < 64; ++j) ccd->w[j] = NULL;
}

static void destroy_ccd_rules(ccd_rules *ccd)
{
     if (ccd) {
	  int j;
	  free(ccd->x);
	  for (j = 0; j < 64; ++j) free(ccd->w[j]);
	  init_ccd_rules(ccd);
     }
}

#define CCD_N(j) ((j)>0 ? (1 << ((j)-1)) : 0)

#define K_PI (3.1415926535897932384626433832795028841971)

/* expand ccd up to rules for the new jmax */
static int grow_ccd_rules(ccd_rules *ccd, int jmax)
{
     int j;
     double *wj, *wj1; /* scratch space for j and j-1 weights */

     if (!ccd || jmax >= 64) return FAILURE;
     if (jmax < 0 || jmax <= ccd->jmax) return SUCCESS;

     ccd->x = (double *) realloc(ccd->x, (1+CCD_N(jmax)) * sizeof(double));
     if (!ccd->x) return FAILURE;

     wj = (double *) malloc((2+CCD_N(jmax)+CCD_N(jmax-1)) * sizeof(double));
     if (!wj) return FAILURE;
     wj1 = wj + 1 + CCD_N(jmax);

     for (j = ccd->jmax + 1; j <= jmax; ++j) {
	  ccd->w[j] = (double *) malloc((1+CCD_N(j)) * sizeof(double));
	  if (!ccd->w[j]) return FAILURE;
     }

     for (j = ccd->jmax + 1; j <= jmax; ++j) {
	  if (j == 0) {
	       ccd->x[0] = 0;
	       ccd->w[j][0] = 2; /* integral of 1 on [-1,1] */
	  }
	  else if (j == 1) {
	       clencurt_weights(CCD_N(j), wj);
	       ccd->x[1] = 1;
	       ccd->w[1][0] = wj[1] - 2;
	       ccd->w[1][1] = wj[0];
	  }
	  else /* (j > 1) */ {
	       int k, i;
	       double *x, *w = ccd->w[j];
	       int nj = CCD_N(j), nj1 = CCD_N(j-1);
	       clencurt_weights(nj, wj);
	       clencurt_weights(nj1, wj1);
	       w[0] = wj[nj] - wj1[nj1];
	       w[1] = wj[0] - wj1[0];
	       for (k = 2; k < j; ++k) {
		    double *wk = w + 1 + CCD_N(k-1);
		    int sk = 1 << (j - k), sk1 = 1 << (j - k - 1);
		    for (i = 0; sk*(2*i+1) < nj; ++i)
			 wk[i] = wj[sk*(2*i+1)] - wj1[sk1*(2*i+1)];
	       }
	       x = ccd->x + 1 + nj1;
	       w += 1 + nj1;
	       for (i = 0; 2*i+1 < nj; ++i) {
		    x[i] = cos((2*i+1) * K_PI / (2 * nj));
		    w[i] = wj[2*i+1];
	       }
	  }
     }

     ccd->jmax = jmax;
     free(wj);
     return SUCCESS;
}

/* total number of distinct weights for a given level j */
static size_t ccd_nf(unsigned j) { return 1 + (j > 0 ? (1U << (j-1)) : 0); }

/* total number of new weights for level j compared to level j-1 */
static size_t ccd_dnf(unsigned j) { 
     return j == 0 ? 1U : (j == 1 ? 1U : (1U << (j-2)));
}

/* actual # function evaluations for a given level compared to previous one,
   counting pairs as 2 */
static size_t ccd_dnf2(unsigned j) {
     return 2 * ccd_dnf(j) - (j == 0);
}

/***************************************************************************/
/* STORAGE OF LEVEL VECTOR */

/* the sparse-grid quadrature is a sum of difference (CCD) rules in
   each dimension with levels J = (j_0, j_1, ...., j_{d-1}) in d
   dimensions for J <= (m_0, m_1, ...) = M and |J|_1 < |M|_\infty,
   where | |_1 is the L_1 norm and |M|_\infty is the L_\infty norm.

   Since J is used as a key in a red-black tree (below), we want to
   make comparisons of two J values quickly.  We do this by storing
   each j_i as one byte (unsigned char), but making the J array equivalent
   to an array of size_t (8 bytes on a 64-bit machine), so that we
   can compare size_t values instead of every j individually. */

typedef unsigned char *Jarr;

/* length of J array in sizeof(size_t)'s, given dim, packing 1 j per
   byte, but padding to an integer multiple of sizeof(size_t) bytes */
static unsigned J_length(unsigned dim) {
     return ((dim + sizeof(size_t) - 1) / sizeof(size_t));
}

/* given J array, return j_i */
static unsigned J_get(const Jarr J, unsigned i) { return J[i]; }

/* set j_i to j (which must be < 256) */
static void J_set(Jarr J, unsigned i, unsigned j) { J[i] = j; }

/* return J array, uninitialized */
static Jarr J_alloc(unsigned dim) {
     return (Jarr) malloc(J_length(dim) * sizeof(size_t));
}

/* set J array to zero */
static void J_zero(Jarr J, unsigned dim) {
     memset(J, 0, J_length(dim) * sizeof(size_t));
}

/* return maximum component of J */
static unsigned J_max(const Jarr J, unsigned dim)
{
     unsigned i, jmax = 0;
     for (i = 0; i < dim; ++i) {
	  unsigned j = J_get(J, i);
	  if (j > jmax) jmax = j;
     }
     return jmax;
}

/* return sum of components of J */
static unsigned J_sum(const Jarr J, unsigned dim)
{
     unsigned i, jsum = 0;
     for (i = 0; i < dim; ++i)
	  jsum += J_get(J, i);
     return jsum;
}

/* return number of dimensions in which j_i == m_i or j_i == L*/
static unsigned J_equal_count(const Jarr J, const Jarr M, unsigned L,
			      unsigned dim) {
     unsigned i, count = 0;
     for (i = 0; i < dim; ++i) {
	  unsigned j = J_get(J, i);
	  count += j == J_get(M, i) || j == L;
     }
     return count;
}

/* return a newly allocated copy of J, or NULL on failure */
static Jarr J_dup(const Jarr J, unsigned dim) {
     unsigned len = J_length(dim);
     Jarr dup = (Jarr) malloc(len * sizeof(size_t));
     if (dup) memcpy(dup, J, len * sizeof(size_t));
     return dup;
}

/* return the number of function evaluations (or pair sums) that need
   to be stored for a given J, not counting evaluations stored for
   smaller J (in any dimension) */
static size_t J_dnf(const Jarr J, unsigned dim) {
     size_t dnf = 1;
     unsigned i;
     for (i = 0; i < dim; ++i)
	  dnf *= ccd_dnf(J_get(J, i));
     return dnf;
}

/* increment the J array so that each j_i runs from 0 to m_i,
   where m_i are the components of the packed array M.
   Returns 0 when incrementing is done. */
static int inc_J(Jarr J, const Jarr M, unsigned dim)
{
     unsigned i, j;
     /* (this could be made more efficient if needed) */
     for (i = 0; i < dim && J_get(J, i) == J_get(M, i); ++i) ;
     if (i == dim) return 0;
     J_set(J, i, J_get(J, i) + 1);
     for (j = 0; j < i; ++j) J_set(J, j, 0);
     return 1;
}

static unsigned imin2(unsigned i, unsigned j) { return(i<j ? i : j); }

/* like inc_J, but also requires L1 norm of J to be <= N */
static int inc_JN(Jarr J, const Jarr M, unsigned N, unsigned dim) {
     unsigned i, j, n = J_sum(J, dim);
     for (i = 0; i < dim; ++i) {
	  unsigned j = J_get(J, i);
	  n -= j;
	  if (j < imin2(J_get(M, i), N - n))
	       break;
     }
     if (i == dim) return 0;
     J_set(J, i, J_get(J, i) + 1);
     for (j = 0; j < i; ++j) J_set(J, j, 0);
     return 1;
}

/* return total # function evals for a given M and max L1 norm = N,
   counting evaluations for all J <= M && |J| <= N, and counting
   actual function evaluations not pairs.  Jt is a scratch array. */
static size_t MN_nf(const Jarr M, unsigned N, unsigned dim,
		     Jarr J) {
     size_t nf = 0;
     J_zero(J, dim);
     do {
	  size_t nf_cur = 1;
	  unsigned i;
	  for (i = 0; i < dim; ++i)
	       nf_cur *= ccd_dnf2(J_get(J, i));
	  nf += nf_cur;
     } while (inc_JN(J, M, N, dim));
     return nf;
}

/***************************************************************************/
/* RED-BLACK TREE OF FUNCTION VALUES */

/* 
   For every J, we store the function evaluations which first appear
   for that J in a red-black tree, indexed by a compressed
   representation of J described above.  By "first appear" we mean the
   function evaluations at points that occur for this J but not if we
   decrease *any* of the levels j_i in J.  Since function evaluations
   always occur in pairs, with equal weights, in a CCD rule, we only
   store the sum of the function evaluations.
*/

typedef struct J_data_s {
     Jarr J; /* packed J vector */
     double *f; /* function evaluations */
} J_data;

static void J_data_destroy(J_data *d) {
     if (d) {
	  free(d->f);
	  free(d->J);
	  free(d);
     }
}

/* create J_data with a copy of J, assuming the functino has dimensionality
   fdim (i.e. each function evaluation returns a vector of fdim>0 values) */
static J_data *J_data_create(const Jarr J, unsigned dim, unsigned fdim) {
     J_data *d = (J_data *) malloc(sizeof(J_data));
     if (!d) return NULL;
     d->J = NULL; d->f = NULL;
     d->J = J_dup(J, dim);
     if (!d->J) { J_data_destroy(d); return NULL; }
     d->f = (double *) malloc(sizeof(double) * (fdim * J_dnf(J, dim)));
     if (!d->f) { J_data_destroy(d); return NULL; }
     return d;
}

static int J_data_compare(const J_data *d1, const J_data *d2, unsigned len) {
     unsigned i;
     size_t *J1 = (size_t *) d1->J, *J2 = (size_t *) d2->J;
     for (i = 0; i < len; ++i)
	  if (J1[i] != J2[i]) 
	       return J1[i] < J2[i] ? -1 : 1;
     return 0;
}

/* red-black tree of J_data */
typedef J_data *rb_key;
typedef unsigned rb_key_data; /* the J length, passed to compare func */
#define rb_destroy_key J_data_destroy
#define rb_compare J_data_compare
#include "redblack.h"

/***************************************************************************/
/* POINT ENUMERATION */

/* the lowest-level integrand takes an array of J_data's to evaluate
   and assigns all of the function values for data->f.   Here, we provide
   some helper functions to enumerate the points for a given J_data. */

/* set nps (length dim) array to the number of (positive-coordinate)
   coordinates for J in each dimension. */
static void J_get_nps(size_t *nps, const Jarr J, unsigned dim) {
     unsigned i;
     for (i = 0; i < dim; ++i)
	  nps[i] = ccd_dnf(J_get(J, i));
}

/* Increment the ks array, so that calling this repeatedly starting
   from ks[i]=0 will cause every ks[i] to range from 0 to nps[i]-1.
   Returns 0 when ks can no longer be incremented. */
static int inc_ks(size_t *ks, const size_t *nps, unsigned dim) {
     unsigned i;
     for (i = 0; i < dim && ks[i] == nps[i] - 1; ++i) ;
     if (i == dim) return 0;
     ks[i] += 1;
     memset(ks, 0, sizeof(size_t) * i);
     return 1;
}

/* like inc_ks, except neg[i] goes from 0 to 1, and we flip the sign
   of the corresponding component of x.  However, if x[i] == 0,
   then that dimension is skipped. */
static int inc_negative(unsigned *neg, double *x, unsigned dim) {
     unsigned i, j;
     for (i = 0; i < dim && (neg[i] || x[i] == 0); ++i) ;
     if (i == dim) return 0;
     neg[i] = 1;
     x[i] = -x[i];
     for (j = 0; j < i; ++j) if (x[j] != 0) { neg[j] = 0; x[j] = -x[j]; }
     return 1;
}

/* set array x[dim] to the ks-th point in J, where 0 <= ks[i] <
   nps[i] with nps set by J_get_nps.  All coordinates are set to
   nonnegative values; the actual function value should be summed over
   positive + negative values in each nonzero coordinate. */
static void J_point(double *x, 
		     const Jarr J, const size_t *ks,
		     unsigned dim,
		     const ccd_rules *ccd) {
     const double *ccdx = ccd->x;
     unsigned i;
     for (i = 0; i < dim; ++i) {
	  unsigned j = J_get(J, i);
	  x[i] = ccdx[j > 0 ? ccd_nf(j-1) + ks[i] : 0];
     }
}

/***************************************************************************/
/* J EVALUATION: evaluate one term J in the integral */

/* evaluate the J-th term in the integrand, where J is d->J, assuming
   that d->f and all J' <= J terms in the tree t have been evaluated.
   Jt points to a scratch array of length J_length(dim), ks is a
   scratch array of length dim, and nps is a scratch array of length
   dim.

   Return results in sums[fdim] array, where fdim is the number of
   integrands. 
*/
static void J_eval(unsigned fdim, double *sums, 
		   const J_data *d, 
		   rb_tree *t, const ccd_rules *ccd, unsigned dim,
		   Jarr Jt, size_t *ks, size_t *nps)
{
     const Jarr J = d->J;
     unsigned fi, i;
     J_data key;

     for (fi = 0; fi < fdim; ++fi) sums[fi] = 0;
     J_zero(Jt, dim);
     key.J = Jt; key.f = NULL;
     do { /* evaluate Jt-specific contributions for all Jt <= J */
	  double *f = rb_tree_find(t, &key)->k->f; /* lookup Jt func. evals */
	  size_t k = 0; /* scalar index corresponding to ks */
	  for (i = 0; i < dim; ++i) ks[i] = 0;
	  J_get_nps(nps, Jt, dim);
	  do { /* loop over ks from 0 to nps */
	       double w = 1;
	       for (i = 0; i < dim; ++i) {
		    unsigned j0 = J_get(J, i);
		    unsigned j = J_get(Jt, i);
		    w *= ccd->w[j0][j > 0 ? ccd_nf(j-1) + ks[i] : 0];
	       }
	       for (fi = 0; fi < fdim; ++fi)
		    sums[fi] += w * f[fdim * k + fi];
	       ++k;
	  } while (inc_ks(ks, nps, dim));
     } while (inc_J(Jt, J, dim));
}

/***************************************************************************/
/* INTEGRATION (low-level internal routine) */

/* internal, low-level integrand: evaluates J[i]->f[...] for
   i = 0 to nJ-1.  */
typedef void (*integrand_)(unsigned nJ, J_data **Je,
			   unsigned dim, unsigned fdim,
			   const ccd_rules *ccd,
			   void *fdata);

/* dimensions and errors, which we sort in descending by errmax as we
   pick dimensions to refine and track the per-dimension errors */
typedef struct {
     unsigned i; /* index of dimension = 0 ... dim-1 */
     double errmax; /* max. integrand error in that dimension */
} dim_error;

/* sort in descending order by errmax */
static int dim_error_compare(const void *d1_, const void *d2_) {
     dim_error *d1 = (dim_error *) d1_;
     dim_error *d2 = (dim_error *) d2_;
     double diff = d1->errmax - d2->errmax;
     return (diff > 0 ? -1 : (diff < 0 ? +1 : 0));
}
static void dim_error_sort(dim_error *d, unsigned dim) {
     qsort(d, dim, sizeof(dim_error), dim_error_compare);
}

/* ensure that J array (length *n_alloc) contains at least n elements,
   growing it as needed */
static J_data **grow_J(J_data **Jd, size_t n, size_t *n_alloc)
{
     if (n > *n_alloc) {
	  *n_alloc = 2*n;
	  Jd = (J_data **) realloc(Jd, sizeof(J_data *) * *n_alloc);
     }
     return Jd;
}

static int converged(double err, double val, double abstol, double reltol) {
     return(err <= abstol || (val != 0 && fabs(err / val) <= reltol));
}

/* Set val and err (arrays of length fdim) to the estimated integral of f
   and corresponding errors, where f is fdim integrands.

   On input, M is an initial maximum J in each dimension.  This are
   incremented adaptively until the specified error tolerances have
   been met for every integrand, or until maxEval is exceeded.

   Returns FAILURE if a memory allocation error occurs, otherwise SUCCESS. */
static int integrate(unsigned dim, unsigned fdim, integrand_ f, void *fdata,
		     Jarr M,
		     size_t maxEval, double reqAbsError, double reqRelError,
		     double *val, double *err)
{
     J_data **Je = NULL; /* array of points to evaluate */
     size_t ie, ne, ne_alloc = 0; /* length of Je array & allocated length */
     rb_tree t; /* red-black tree of cubature terms J and function values */
     Jarr J = NULL;
     unsigned i, fi, di, N;
     ccd_rules ccd;
     int ret = FAILURE;
     size_t numEval;
     double *Jsum = NULL; /* length fdim array of per-J integral contrib. */
     double *derrs = NULL; /* dim x fdim array of per-dimension err. est. */
     double *rem_err = NULL; /* fdim array of remaining errors after pops */
     size_t *scratch = NULL; /* 2*dim scratch array */
     dim_error *dims = NULL; /* dims[dim] array for picking worst dims */
     
     if (fdim == 0) return SUCCESS; /* no integrands */

     init_ccd_rules(&ccd);
     rb_tree_init(&t, J_length(dim));

     if (dim == 0) { /* trivial 0-dimensional "integral" = 1 f evaluation */
	  J_data J, *pJ;
	  J.J = NULL;
	  J.f = val;
	  pJ = &J;
	  f(1, &pJ, dim, fdim, &ccd, fdata);
	  for (fi = 0; fi < fdim; ++fi) err[fi] = 0;
	  return SUCCESS;
     }

     dims = (dim_error *) malloc(dim * sizeof(dim_error));
     if (!dims) goto done;

     J = J_alloc(dim);
     if (!J) goto done;
     derrs = (double *) malloc(sizeof(double) * (dim * fdim));
     if (!derrs) goto done;
     Jsum = (double *) malloc(sizeof(double) * fdim);
     if (!Jsum) goto done;
     scratch = (size_t *) malloc(sizeof(size_t) * 2 * dim);
     if (!scratch) goto done;
     rem_err = Jsum; /* don't need at same time as Jsum, can re-use */
     
     memset(derrs, 0, sizeof(double) * (dim * fdim));
     memset(val, 0, sizeof(double) * fdim);
     memset(err, 0, sizeof(double) * fdim);

     if (FAILURE == grow_ccd_rules(&ccd, J_max(M, dim))) goto done;
     
     ne = 0;
     J_zero(J, dim);
     N = J_max(M, dim);
     do {
	  if (!(Je = grow_J(Je, ++ne, &ne_alloc))) goto done;
	  Je[ne-1] = J_data_create(J, dim, fdim);
	  if (!Je[ne-1]) goto done;
	  if (!rb_tree_insert(&t, Je[ne-1])) {
	       J_data_destroy(Je[ne-1]);
	       goto done;
	  }
     } while (inc_JN(J, M, N, dim));
     numEval = MN_nf(M, N, dim, J);
     f(ne, Je, dim, fdim, &ccd, fdata);

     for (ie = 0; ie < ne; ++ie) {
	  unsigned jsum = J_sum(Je[ie]->J, dim);
	  unsigned jmax = jsum == N ? J_max(Je[ie]->J, dim) : N+1;
	  unsigned count = J_equal_count(Je[ie]->J, M, jmax, dim);
	  J_eval(fdim, Jsum, Je[ie], &t, &ccd, dim, J,
		 scratch, scratch + dim);
	  for (fi = 0; fi < fdim; ++fi) val[fi] += Jsum[fi];
	  if (count > 0)
	       for (i = 0; i < dim; ++i) {
		    unsigned j = J_get(Je[ie]->J, i);
		    if (j == J_get(M, i) || j == jmax) {
			 for (fi = 0; fi < fdim; ++fi) {
			      double erri = fabs(Jsum[fi]) / count;
			      err[fi] += erri;
			      derrs[i*fdim + fi] += erri;
			 }
		    }
	       }
     }

     /* set up array of dimensions and errors for refinement */
     for (i = 0; i < dim; ++i) {
	  double errmax = 0;
	  for (fi = 0; fi < fdim; ++fi)
	       if (derrs[i*fdim + fi] > errmax)
		    errmax = derrs[i*fdim + fi];
	  dims[i].i = i;
	  dims[i].errmax = errmax;
     }

     while (numEval < maxEval || !maxEval) {
	  unsigned Nprev = N;

	  /* Convergence check: */
	  for (fi=0; fi < fdim && converged(err[fi], val[fi],
					    reqAbsError, reqRelError); ++fi) ;
	  if (fi == fdim)
	       break; /* convergence */

	  /* Refine dimensions: */

	  dim_error_sort(dims, dim);
	  for (fi = 0; fi < fdim; ++fi) rem_err[fi] = err[fi];
	  /* increment M in all dimensions to be refined */
	  di = 0;
	  do {
	       i = dims[di].i;
	       J_set(M, i, J_get(M, i) + 1);
	       N = J_max(M, dim);
	       numEval = MN_nf(M, N, dim, J);
	       for (fi = 0; fi < fdim; ++fi) 
		    rem_err[fi] -= derrs[i*fdim + fi];
	       memset(derrs + i*fdim, 0, sizeof(double) * fdim);
	       for (fi=0; fi < fdim && converged(rem_err[fi], val[fi],
						 reqAbsError, reqRelError);
		    ++fi) ;
	       ++di;
	       if (fi == fdim) break; /* other regions have small errs */
	  } while (di < dim && (numEval < maxEval || !maxEval));
	  ne = 0;
	  J_zero(J, dim);
	  /* add all new J's (note: generic inc_JN is inefficient here) */
	  while (inc_JN(J, M, N, dim)) { 
	       for (i = 0; i < di /* first di dims were incremented */
			 && J_get(J,dims[i].i) < J_get(M,dims[i].i); ++i) ;
	       if (i == di && (N == Nprev || J_sum(J,dim) < N))
		   continue; /* not a new point */

	       if (!(Je = grow_J(Je, ++ne, &ne_alloc))) goto done;
	       Je[ne-1] = J_data_create(J, dim, fdim);
	       if (!Je[ne-1]) goto done;
	       if (!rb_tree_insert(&t, Je[ne-1])) {
		    J_data_destroy(Je[ne-1]);
		    goto done;
	       }
	  }
	  if (FAILURE == grow_ccd_rules(&ccd, J_max(M, dim))) goto done;
	  /* evaulate integrand at new points */
	  f(ne, Je, dim, fdim, &ccd, fdata);
	  /* accumulate new terms in integrand and errors */
	  for (ie = 0; ie < ne; ++ie) {
	       unsigned jsum = J_sum(Je[ie]->J, dim);
	       unsigned jmax = jsum == N ? J_max(Je[ie]->J, dim) : N+1;
	       unsigned count = J_equal_count(Je[ie]->J, M, jmax, dim);
	       J_eval(fdim, Jsum, Je[ie], &t, &ccd, dim, J,
		      scratch, scratch + dim);
	       for (fi = 0; fi < fdim; ++fi) val[fi] += Jsum[fi];
	       if (count > 0)
		    for (i = 0; i < dim; ++i) {
			 unsigned j = J_get(Je[ie]->J, i);
			 if (j == J_get(M, i) || j == jmax)
			      for (fi = 0; fi < fdim; ++fi)
				   derrs[i*fdim + fi] += fabs(Jsum[fi])/count;
		    }
	  }
	  for (fi = 0; fi < fdim; ++fi) {
	       err[fi] = 0;
	       for (i = 0; i < dim; ++i) err[fi] += derrs[i*fdim + fi];
	  }
	  for (di = 0; di < dim; ++di) {
	       i = dims[di].i;
	       dims[di].errmax = 0;
	       for (fi = 0; fi < fdim; ++fi)
		    if (derrs[i*fdim + fi] > dims[di].errmax)
			 dims[di].errmax = derrs[i*fdim + fi];
	  }
     }

     ret = SUCCESS;
done:
     free(Je);
     free(scratch);
     free(Jsum);
     free(derrs);
     free(J);
     free(dims);
     destroy_ccd_rules(&ccd);
     rb_tree_destroy(&t);
     return ret;
}
		     

/***************************************************************************/
/* SERIAL API: user integrand is evaluated one pt at a time */

typedef struct {
     integrand f;
     void *fdata;
     const double *xmin;
     const double *xmax;
     double *x0, *x, *fval;
     size_t *ks, *nps;
     unsigned *negative;
} sintegrand_data;

static void sintegrand(unsigned nJ, J_data **Je,
		       unsigned dim, unsigned fdim,
		       const ccd_rules *ccd,
		       void *d_) {
     sintegrand_data *d = (sintegrand_data *) d_;
     integrand f = d->f;
     void *fdata = d->fdata;
     double *x0 = d->x0, *x = d->x, *fval = d->fval;
     const double *xmin = d->xmin;
     const double *xmax = d->xmax;
     size_t *ks = d->ks, *nps = d->nps;
     unsigned *negative = d->negative;
     unsigned iJ, i;

     for (iJ = 0; iJ < nJ; ++iJ) { /* J points */
	  Jarr J = Je[iJ]->J;
	  unsigned k = 0;
	  memset(ks, 0, sizeof(size_t) * dim);
	  J_get_nps(nps, J, dim);
	  do { /* x points "owned" by this J */
	       double *fJ = Je[iJ]->f + fdim * k;
	       J_point(x0, J, ks, dim, ccd);
	       memset(negative, 0, sizeof(int) * dim);
	       memset(fJ, 0, sizeof(double) * fdim);
	       do { /* accumulate all possible sign flips */
		    for (i = 0; i < dim; ++i)
			 x[i] = xmin[i] + (x0[i]+1) * 0.5 * (xmax[i]-xmin[i]);
		    f(dim, x, fdata, fdim, fval);
		    for (i = 0; i < fdim; ++i)
			 fJ[i] += fval[i];
	       } while (inc_negative(negative, x0, dim));
	       ++k;
	  } while (inc_ks(ks, nps, dim));
     }
}

int sadapt_integrate(unsigned fdim, integrand f, void *fdata,
		     unsigned dim, const double *xmin, const double *xmax, 
		     unsigned maxEval, double reqAbsError, double reqRelError, 
		     double *val, double *err)
{
     sintegrand_data d;
     double *scratch = NULL;
     size_t *iscratch = NULL;
     Jarr M = NULL;
     int ret = FAILURE;
     unsigned i;
     double volscale = 1; /* scale factor of integration volume */
     
     /* trivial cases: */
     if (fdim == 0) return SUCCESS;
     if (dim == 0) {
	  f(dim, xmin, fdata, fdim, val);
	  memset(err, 0, sizeof(double) * fdim);
	  return SUCCESS;
     }

     for (i = 0; i < dim; ++i)
	  volscale *= fabs(xmax[i] - xmin[i]) * 0.5;
     if (volscale == 0) { /* empty integration domain */
	  memset(val, 0, sizeof(double) * fdim);
	  memset(err, 0, sizeof(double) * fdim);
	  return SUCCESS;
     }
     reqAbsError /= volscale;

     d.negative = NULL;
     scratch = (double *) malloc(sizeof(double) * (2*dim + fdim));
     if (!scratch) goto done;
     iscratch = (size_t *) malloc(sizeof(size_t) * (dim * 2));
     if (!iscratch) goto done;
     d.negative = (unsigned *) malloc(sizeof(int) * dim);
     if (!d.negative) goto done;

     M = J_alloc(dim);
     if (!M) goto done;
     J_zero(M, dim); /* initial M == 0, so integration starts with 1 pt */

     d.f = f;
     d.fdata = fdata;
     d.xmin = xmin;
     d.xmax = xmax;
     d.x0 = scratch;
     d.x = d.x0 + dim;
     d.fval = d.x + dim;
     d.ks = iscratch;
     d.nps = d.ks + dim;

     ret = integrate(dim, fdim, sintegrand, &d, M,
		     maxEval, reqAbsError, reqRelError, val, err);

     for (i = 0; i < fdim; ++i) {
	  val[i] *= volscale;
	  err[i] *= volscale;
     }

done:
     free(M);
     free(d.negative);
     free(iscratch);
     free(scratch);
     return ret;
}

/***************************************************************************/
