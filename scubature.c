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

   We only allow j to go from 0 to 31, which allows us to use a 4-bit
   packed format for storing j, mainly for comparison efficiency in
   high-dimensional cubature more than storage efficiency (16 j's per
   64-bit integer).  If you need more than 1+2^31 points (and much
   more than that in multidimensions), you probably need a different
   integration algorithm anyway.

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
     int jmax; /* j <= jmax < 32 has been allocated */
     double *x; /* points x[i] in increasing j order; see above */
     double *w[32]; /* weights w[j][i] for x[i] of level-j CCD rule */
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
     for (j = 0; j < 32; ++j) ccd->w[j] = NULL;
}

static void destroy_ccd_rules(ccd_rules *ccd)
{
     if (ccd) {
	  int j;
	  free(ccd->x);
	  for (j = 0; j < 32; ++j) free(ccd->w[j]);
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

     if (!ccd || jmax >= 32) return FAILURE;
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

/* actual # function evaluations for a given level, counting pairs as 2 */
static size_t ccd_nf2(unsigned j) {
     return 1 + (j > 0 ? (1U << j) : 0);
}

/***************************************************************************/
/* COMPRESSED STORAGE OF LEVEL VECTOR */

/* the sparse-grid quadrature is a sum of difference (CCD) rules in
   each dimension with levels J = (j_0, j_1, ...., j_{d-1}) in d
   dimensions for J <= (m_0, m_1, ...) = M and |J|_1 < |M|_\infty,
   where | |_1 is the L_1 norm and |M|_\infty is the L_\infty norm.

   Since J is used as a key in a red-black tree (below), we want to
   make comparisons of two J values quickly.  We do this by packing
   the levels j_i into arrays Jp of unsigned integers at least the
   same size as pointers (type size_t), 4 bits each (allowing 0 <= j_i
   < 32).  On a 64-bit machine, this packs 16 dimensions per integer,
   and we use the ordinary integer comparison function for the tree
   (in dictionary order if we have more than 16 dimensions).
*/

/* length of Jp array, given dim, packing 2 j's per byte */
static unsigned Jp_length(unsigned dim) {
     return (dim + (sizeof(size_t) * 2 - 1)) / (sizeof(size_t) * 2);
}

/* given Jp array, return j_i */
static unsigned Jp_get(const size_t *Jp, unsigned i) {
     unsigned i0 = i % (sizeof(size_t) * 2);
     return (Jp[i / (sizeof(size_t) * 2)] >> (i0 * 4)) & 0xf;
}

/* set j_i to j (which must be < 32) */
static void Jp_set(size_t *Jp, unsigned i, unsigned j) {
     unsigned i0 = i % (sizeof(size_t) * 2);
     unsigned i1 = i / (sizeof(size_t) * 2);
     size_t Jp1 = Jp[i1];
     unsigned j0 = (Jp1 >> (i0 * 4)) & 0xf; /* Jp_get(Jp, i) */
     Jp[i1] = Jp1 ^ ((j0 ^ j) << (i0 * 4));
}

/* return Jp array, uninitialized */
static size_t *Jp_alloc(unsigned dim) {
     return (size_t *) malloc(Jp_length(dim) * sizeof(size_t));
}

/* set Jp array to zero */
static void Jp_zero(size_t *Jp, unsigned dim) {
     memset(Jp, 0, Jp_length(dim) * sizeof(size_t));
}

/* return maximum component of Jp */
static unsigned Jp_max(const size_t *Jp, unsigned dim)
{
     unsigned i, jmax = 0;
     for (i = 0; i < dim; ++i) {
	  unsigned j = Jp_get(Jp, i);
	  if (j > jmax) jmax = j;
     }
     return jmax;
}

/* return sum of components of Jp */
static unsigned Jp_sum(const size_t *Jp, unsigned dim)
{
     unsigned i, jsum = 0;
     for (i = 0; i < dim; ++i)
	  jsum += Jp_get(Jp, i);
     return jsum;
}

/* return number of dimensions in which j_i == m_i */
static unsigned Jp_equal_count(const size_t *Jp, const size_t *Mp,
			       unsigned dim) {
     unsigned i, count = 0;
     for (i = 0; i < dim; ++i)
	  count += Jp_get(Jp, i) == Jp_get(Mp, i);
     return count;
}

/* return a newly allocated copy of Jp, or NULL on failure */
static size_t *Jp_dup(const size_t *Jp, unsigned dim) {
     unsigned len = Jp_length(dim);
     size_t *dup = (size_t *) malloc(len * sizeof(size_t));
     if (dup) memcpy(dup, Jp, len * sizeof(size_t));
     return dup;
}

/* return the number of function evaluations (or pair sums) that need
   to be stored for a given J, not counting evaluations stored for
   smaller J (in any dimension) */
static size_t Jp_dnf(const size_t *Jp, unsigned dim) {
     size_t dnf = 1;
     unsigned i;
     for (i = 0; i < dim; ++i)
	  dnf *= ccd_dnf(Jp_get(Jp, i));
     return dnf;
}

/* increment the Jp array so that each j_i runs from 0 to m_i,
   where m_i are the components of the packed array Mp.
   Returns 0 when incrementing is done. */
static int inc_Jp(size_t *Jp, const size_t *Mp, unsigned dim)
{
     unsigned i;
     /* (this could be made more efficient if needed) */
     for (i = 0; i < dim && Jp_get(Jp, i) == Jp_get(Mp, i); ++i) ;
     if (i == dim) return 0;
     Jp_set(Jp, i, Jp_get(Jp, i) + 1);
     for (i = i - 1; i >= 0; --i) Jp_set(Jp, i, 0);
     return 1;
}

static unsigned imin2(unsigned i, unsigned j) { return(i<j ? i : j); }

/* like inc_Jp, but also requires L1 norm of J to be <= N */
static int inc_JpN(size_t *Jp, const size_t *Mp, unsigned N, unsigned dim) {
     unsigned i, n = Jp_sum(Jp, dim);
     for (i = 0; i < dim; ++i) {
	  unsigned j = Jp_get(Jp, i);
	  n -= j;
	  if (j < imin2(Jp_get(Mp, i), N - n))
	       break;
     }
     if (i == dim) return 0;
     Jp_set(Jp, i, Jp_get(Jp, i) + 1);
     for (i = i - 1; i >= 0; --i) Jp_set(Jp, i, 0);
     return 1;
}

/* return total # function evals for a given Mp and max L1 norm = N,
   counting evaluations for all J <= M && |J| <= N, and counting
   actual function evaluations not pairs.  Jt is a scratch array. */
static size_t MpN_nf(const size_t *Mp, unsigned N, unsigned dim,
		     size_t *Jp) {
     size_t nf = 0, nf_cur;
     Jp_zero(Jp, dim);
     N = imin2(N, Jp_sum(Mp, dim)); /* ensure |J|==N is achieved for J <= M */
     do { /* inefficient way to loop over all |J|==N */
	  unsigned i;
	  if (Jp_sum(Jp, dim) != N) continue;
	  nf_cur = 1;
	  for (i = 0; i < dim; ++i)
	       nf_cur *= ccd_nf2(Jp_get(Jp, i));
	  nf += nf_cur;
     } while (inc_JpN(Jp, Mp, N, dim));
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
     size_t *Jp; /* packed J vector */
     double *f; /* function evaluations */
} J_data;

static void J_data_destroy(J_data *d) {
     if (d) {
	  free(d->f);
	  free(d->Jp);
	  free(d);
     }
}

/* create J_data with a copy of Jp, assuming the functino has dimensionality
   fdim (i.e. each function evaluation returns a vector of fdim>0 values) */
static J_data *J_data_create(const size_t *Jp, unsigned dim, unsigned fdim) {
     J_data *d = (J_data *) malloc(sizeof(J_data));
     if (!d) return NULL;
     d->Jp = NULL; d->f = NULL;
     d->Jp = Jp_dup(Jp, dim);
     if (!d->Jp) { J_data_destroy(d); return NULL; }
     d->f = (double *) malloc(sizeof(double) * (fdim * Jp_dnf(Jp, dim)));
     if (!d->f) { J_data_destroy(d); return NULL; }
     return d;
}

static int J_data_compare(const J_data *d1, const J_data *d2, unsigned dim) {
     unsigned i;
     size_t *Jp1 = d1->Jp, *Jp2 = d2->Jp;
     for (i = 0; i < dim; ++i)
	  if (Jp1[i] != Jp2[i]) 
	       return Jp1[i] < Jp2[i] ? -1 : 1;
     return 0;
}

/* red-black tree of J_data */
typedef J_data *rb_key;
typedef unsigned rb_key_data; /* the dimension, passed to compare func */
#define rb_destroy_key J_data_destroy
#define rb_compare J_data_compare
#include "redblack.h"

/***************************************************************************/
/* POINT ENUMERATION */

/* the lowest-level integrand takes an array of J_data's to evaluate
   and assigns all of the function values for data->f.   Here, we provide
   some helper functions to enumerate the points for a given J_data. */

/* set nps (length dim) array to the number of (positive-coordinate)
   coordinates for Jp in each dimension. */
static void Jp_get_nps(size_t *nps, const size_t *Jp, unsigned dim) {
     unsigned i;
     for (i = 0; i < dim; ++i)
	  nps[i] = ccd_dnf(Jp_get(Jp, i));
}

/* Increment the ks array, so that calling this repeatedly starting
   from ks[i]=0 will cause every ks[i] to range from 0 to nps[i]-1.
   Returns 0 when ks can no longer be incremented. */
static int inc_ks(size_t *ks, const size_t *nps, unsigned dim) {
     unsigned i;
     for (i = 0; i < dim && ks[i] == nps[i] - 1; ++i) ;
     if (i == dim) return 0;
     ks[i] += 1;
     for (i = i - 1; i >= 0; --i) ks[i] = 0;
     return 1;
}

/* like inc_ks, except neg[i] goes from 0 to 1, and we flip the sign
   of the corresponding component of x.  However, if x[i] == 0,
   then that dimension is skipped. */
static int inc_negative(unsigned *neg, double *x, unsigned dim) {
     unsigned i;
     for (i = 0; i < dim && (neg[i] || x[i] == 0); ++i) ;
     if (i == dim) return 0;
     neg[i] = 1;
     x[i] = -x[i];
     for (i = i - 1; i >= 0; --i) if (x[i] != 0) { neg[i] = 0; x[i] = -x[i]; }
     return 1;
}

/* set array x[dim] to the ks-th point in Jp, where 0 <= ks[i] <
   nps[i] with nps set by Jp_get_nps.  All coordinates are set to
   nonnegative values; the actual function value should be summed over
   positive + negative values in each nonzero coordinate. */
static void Jp_point(double *x, 
		     const size_t *Jp, const size_t *ks,
		     unsigned dim,
		     const ccd_rules *ccd) {
     const double *ccdx = ccd->x;
     unsigned i;
     for (i = 0; i < dim; ++i) {
	  unsigned j = Jp_get(Jp, i);
	  x[i] = ccdx[j > 0 ? ccd_nf(j-1) + ks[i] : 0];
     }
}

/***************************************************************************/
/* J EVALUATION: evaluate one term J in the integral */

/* evaluate the J-th term in the integrand, where J is d->Jp, assuming
   that d->f and all J' <= J terms in the tree t have been evaluated.
   Jt points to a scratch array of length Jp_length(dim), ks is a
   scratch array of length dim, and nps is a scratch array of length
   dim.

   Return results in sums[fdim] array, where fdim is the number of
   integrands. 
*/
static void J_eval(unsigned fdim, double *sums, 
		   const J_data *d, 
		   rb_tree *t, const ccd_rules *ccd, unsigned dim,
		   size_t *Jt, size_t *ks, size_t *nps)
{
     const size_t *Jp = d->Jp;
     unsigned fi, i;
     J_data key;

     for (fi = 0; fi < fdim; ++fi) sums[fi] = 0;
     Jp_zero(Jt, dim);
     key.Jp = Jt; key.f = NULL;
     do { /* evaluate Jt-specific contributions for all Jt <= Jp */
	  double *f = rb_tree_find(t, &key)->k->f; /* lookup Jt func. evals */
	  size_t k = 0; /* scalar index corresponding to ks */
	  for (i = 0; i < dim; ++i) ks[i] = 0;
	  Jp_get_nps(nps, Jt, dim);
	  do { /* loop over ks from 0 to nps */
	       double w = 1;
	       for (i = 0; i < dim; ++i) {
		    unsigned j = Jp_get(Jt, i);
		    w *= ccd->w[j][j > 0 ? ccd_nf(j-1) + ks[i] : 0];
	       }
	       for (fi = 0; fi < fdim; ++fi)
		    sums[fi] += w * f[fdim * k + fi];
	       ++k;
	  } while (inc_ks(ks, nps, dim));
     } while (inc_Jp(Jt, Jp, dim));
}

/***************************************************************************/
/* INTEGRATION (low-level internal routine) */

/* internal, low-level integrand: evaluates J[i]->f[...] for
   i = 0 to nJ-1.  */
typedef void (*integrand_)(unsigned nJ, J_data **J,
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
static J_data **grow_J(J_data **J, size_t n, size_t *n_alloc)
{
     if (n > *n_alloc) {
	  *n_alloc = 2*n;
	  J = (J_data **) realloc(J, sizeof(J_data *) * *n_alloc);
     }
     return J;
}

static int converged(double err, double val, double abstol, double reltol) {
     return(err <= abstol || (val != 0 && fabs(err / val) <= reltol));
}

/* Set val and err (arrays of length fdim) to the estimated integral of f
   and corresponding errors, where f is fdim integrands.

   On input, Mp is an initial maximum J in each dimension.  This are
   incremented adaptively until the specified error tolerances have
   been met for every integrand, or until maxEval is exceeded.

   Returns FAILURE if a memory allocation error occurs, otherwise SUCCESS. */
static int integrate(unsigned dim, unsigned fdim, integrand_ f, void *fdata,
		     size_t *Mp,
		     size_t maxEval, double reqAbsError, double reqRelError,
		     double *val, double *err)
{
     J_data **Je = NULL; /* array of points to evaluate */
     size_t ie, ne, ne_alloc = 0; /* length of Je array & allocated length */
     rb_tree t; /* red-black tree of cubature terms J and function values */
     size_t *Jp = NULL;
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
     rb_tree_init(&t, dim);

     if (dim == 0) { /* trivial 0-dimensional "integral" = 1 f evaluation */
	  J_data J, *pJ;
	  J.Jp = NULL;
	  J.f = val;
	  pJ = &J;
	  f(1, &pJ, dim, fdim, &ccd, fdata);
	  for (fi = 0; fi < fdim; ++fi) err[fi] = 0;
	  return SUCCESS;
     }

     dims = (dim_error *) malloc(dim * sizeof(dim_error));
     if (!dims) goto done;

     Jp = Jp_alloc(dim);
     if (!Jp) goto done;
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

     if (FAILURE == grow_ccd_rules(&ccd, Jp_max(Mp, dim))) goto done;
     
     ne = 0;
     Jp_zero(Jp, dim);
     N = Jp_max(Mp, dim);
     do {
	  if (!(Je = grow_J(Je, ++ne, &ne_alloc))) goto done;
	  Je[ne-1] = J_data_create(Jp, dim, fdim);
	  if (!Je[ne-1]) goto done;
	  if (!rb_tree_insert(&t, Je[ne-1])) {
	       J_data_destroy(Je[ne-1]);
	       goto done;
	  }
     } while (inc_JpN(Jp, Mp, N, dim));
     numEval = MpN_nf(Mp, N, dim, Jp);
     f(ne, Je, dim, fdim, &ccd, fdata);

     for (ie = 0; ie < ne; ++ie) {
	  unsigned count = Jp_equal_count(Je[ie]->Jp, Mp, dim);
	  J_eval(fdim, Jsum, Je[ie], &t, &ccd, dim, Jp,
		 scratch, scratch + dim);
	  for (fi = 0; fi < fdim; ++fi) val[fi] += Jsum[fi];
	  if (count > 0)
	       for (i = 0; i < dim; ++i)
		    if (Jp_get(Je[ie]->Jp, i) == Jp_get(Mp, i)) {
			 for (fi = 0; fi < fdim; ++fi) {
			      double erri = Jsum[fi] / count;
			      err[fi] += erri;
			      derrs[i*fdim + fi] += erri;
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
	       Jp_set(Mp, i, Jp_get(Mp, i) + 1);
	       N = Jp_max(Mp, dim);
	       numEval = MpN_nf(Mp, N, dim, Jp);
	       for (fi = 0; fi < fdim; ++fi) 
		    rem_err[fi] -= derrs[i*fdim + fi];
	       memset(derrs + i*fdim, 0, sizeof(double) * fdim);
	       for (fi=0; fi < fdim && converged(rem_err[fi], val[fi],
						 reqAbsError, reqRelError);
		    ++fi) ;
	       if (fi == fdim) break; /* other regions have small errs */
	       ++di;
	  } while (di < dim && (numEval < maxEval || !maxEval));
	  ne = 0;
	  Jp_zero(Jp, dim);
	  /* add all new J's (note: generic inc_JpN is inefficient here) */
	  while (inc_JpN(Jp, Mp, N, dim)) { 
	       for (i = 0; i < di /* first di dims were incremented */
			 && Jp_get(Jp,dims[i].i) < Jp_get(Mp,dims[i].i); ++i) ;
	       if (i == di) continue; /* not a new point */
	       if (!(Je = grow_J(Je, ++ne, &ne_alloc))) goto done;
	       Je[ne-1] = J_data_create(Jp, dim, fdim);
	       if (!Je[ne-1]) goto done;
	       if (!rb_tree_insert(&t, Je[ne-1])) {
		    J_data_destroy(Je[ne-1]);
		    goto done;
	       }
	  }
	  if (FAILURE == grow_ccd_rules(&ccd, Jp_max(Mp, dim))) goto done;
	  /* evaulate integrand at new points */
	  f(ne, Je, dim, fdim, &ccd, fdata);
	  /* accumulate new terms in integrand and errors */
	  for (ie = 0; ie < ne; ++ie) {
	       unsigned count = Jp_equal_count(Je[ie]->Jp, Mp, dim);
	       J_eval(fdim, Jsum, Je[ie], &t, &ccd, dim, Jp,
		      scratch, scratch + dim);
	       for (fi = 0; fi < fdim; ++fi) val[fi] += Jsum[fi];
	       if (count > 0)
		    for (i = 0; i < dim; ++i)
			 if (Jp_get(Je[ie]->Jp, i) == Jp_get(Mp, i))
			      for (fi = 0; fi < fdim; ++fi)
				   derrs[i*fdim + fi] += Jsum[fi] / count;
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
     free(Jp);
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

static void sintegrand(unsigned nJ, J_data **J,
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
	  size_t *Jp = J[iJ]->Jp;
	  unsigned k = 0;
	  memset(ks, 0, sizeof(size_t) * dim);
	  Jp_get_nps(nps, Jp, dim);
	  do { /* x points "owned" by this J */
	       double *fJ = J[iJ]->f + fdim * k;
	       Jp_point(x0, Jp, ks, dim, ccd);
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
     size_t *iscratch = NULL, *Mp = NULL;
     int ret = FAILURE;
     
     /* trivial cases: */
     if (fdim == 0) return SUCCESS;
     if (dim == 0) {
	  f(dim, xmin, fdata, fdim, val);
	  memset(err, 0, sizeof(double) * fdim);
     }

     d.negative = NULL;
     scratch = (double *) malloc(sizeof(double) * (2*dim + fdim));
     if (!scratch) goto done;
     iscratch = (size_t *) malloc(sizeof(size_t) * (dim * 2));
     if (!iscratch) goto done;
     d.negative = (unsigned *) malloc(sizeof(int) * dim);
     if (!d.negative) goto done;

     Mp = Jp_alloc(dim);
     if (!Mp) goto done;
     Jp_zero(Mp, dim); /* initial M == 0, so integration starts with 1 pt */
     
     d.f = f;
     d.fdata = fdata;
     d.xmin = xmin;
     d.xmax = xmax;
     d.x0 = scratch;
     d.x = d.x0 + dim;
     d.fval = d.x + dim;
     d.ks = iscratch;
     d.nps = d.ks + dim;

     ret = integrate(dim, fdim, sintegrand, &d, Mp,
		     maxEval, reqAbsError, reqRelError, val, err);

done:
     free(Mp);
     free(d.negative);
     free(iscratch);
     free(scratch);
     return ret;
}

/***************************************************************************/
