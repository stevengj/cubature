/* Adaptive multidimensional integration of a vector of integrands.
 *
 * Copyright (c) 2005-2009 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* Adaptive multidimensional integration on hypercubes (or, really,
   hyper-rectangles) using cubature rules.

   A cubature rule takes a function and a hypercube and evaluates
   the function at a small number of points, returning an estimate
   of the integral as well as an estimate of the error, and also
   a suggested dimension of the hypercube to subdivide.

   Given such a rule, the adaptive integration is simple:

   1) Evaluate the cubature rule on the hypercube(s).
      Stop if converged.

   2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.

   3) Goto (1).

 The basic algorithm is based on the adaptive cubature described in
 
     A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric
     integration over an N-dimensional rectangular region,"
     J. Comput. Appl. Math. 6 (4), 295-302 (1980).

 and subsequently extended to integrating a vector of integrands in

     J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm
     for the approximate calculation of multiple integrals,"
     ACM Trans. Math. Soft. 17 (4), 437-451 (1991).

 Note, however, that we do not use any of code from the above authors
 (in part because their code is Fortran 77, but mostly because it is
 under the restrictive ACM copyright license).  I did make use of some
 GPL code from Rudolf Schuerer's HIntLib and from the GNU Scientific
 Library as listed in the copyright notice above, on the other hand.

 TODO:

   * Putting these routines into the GNU GSL library would be nice.

   * It would be very nice to implement some kind of parallelization
     scheme.  e.g. at the very least, within a given cubature rule
     one could call the integrand with an array of points to evaluate
     (which the integrand could then farm out to separate threads etc.).

   * For high-dimensional integrals, it would be nice to implement
     a sparse-grid cubature scheme using Clenshaw-Curtis quadrature.
     Currently, for dimensions > 7 or so, quasi Monte Carlo methods win.

   * Berntsen et. al also describe a "two-level" error estimation scheme
     that they claim makes the algorithm more robust.  It might be
     nice to implement this, at least as an option (although I seem
     to remember trying it once and it made the number of evaluations
     substantially worse for my test integrands).

   * In the cubature-rule evaluation for each hypercube, I end up
     allocating and then deallocating a bunch of temporary arrays
     to hold the integrand vectors.  One could consolidate these
     by allocating a single scratch workspace in adapt_integrate,
     if for some reason the overhead of malloc is a concern.
*/

/* USAGE: Call adapt_integrate with your function as described in cubature.h.

	  To compile a test program, compile cubature.c with
	  -DTEST_INTEGRATOR as described at the end. */

#include "cubature.h"

/***************************************************************************/
/* Basic datatypes */

typedef struct {
     double val, err;
} esterr;

static double relError(esterr ee)
{
     return (ee.val == 0.0 ? HUGE_VAL : fabs(ee.err / ee.val));
}

typedef struct {
     unsigned dim;
     double *data;	/* length 2*dim = center followed by half-widths */
     double vol;	/* cache volume = product of widths */
} hypercube;

static double compute_vol(const hypercube *h)
{
     unsigned i;
     double vol = 1;
     for (i = 0; i < h->dim; ++i)
	  vol *= 2 * h->data[i + h->dim];
     return vol;
}

static hypercube make_hypercube(unsigned dim, const double *center, const double *halfwidth)
{
     unsigned i;
     hypercube h;
     h.dim = dim;
     h.data = (double *) malloc(sizeof(double) * dim * 2);
     for (i = 0; i < dim; ++i) {
	  h.data[i] = center[i];
	  h.data[i + dim] = halfwidth[i];
     }
     h.vol = compute_vol(&h);
     return h;
}

static hypercube make_hypercube_range(unsigned dim, const double *xmin, const double *xmax)
{
     hypercube h = make_hypercube(dim, xmin, xmax);
     unsigned i;
     for (i = 0; i < dim; ++i) {
	  h.data[i] = 0.5 * (xmin[i] + xmax[i]);
	  h.data[i + dim] = 0.5 * (xmax[i] - xmin[i]);
     }
     h.vol = compute_vol(&h);
     return h;
}

static void destroy_hypercube(hypercube *h)
{
     free(h->data);
     h->dim = 0;
}

typedef struct {
     hypercube h;
     unsigned splitDim;
     unsigned fdim; /* dimeinsionality of vector integrand */
     esterr *ee; /* array of length fdim */
     double errmax; /* max ee[k].err */
} region;

static region make_region(const hypercube *h, unsigned fdim)
{
     region R;
     R.h = make_hypercube(h->dim, h->data, h->data + h->dim);
     R.splitDim = 0;
     R.fdim = fdim;
     R.ee = (esterr *) malloc(sizeof(esterr) * fdim);
     return R;
}

static void destroy_region(region *R)
{
     destroy_hypercube(&R->h);
     free(R->ee);
     R->ee = 0;
}

static void cut_region(region *R, region *R2)
{
     unsigned d = R->splitDim, dim = R->h.dim;
     *R2 = *R;
     R->h.data[d + dim] *= 0.5;
     R->h.vol *= 0.5;
     R2->h = make_hypercube(dim, R->h.data, R->h.data + dim);
     R->h.data[d] -= R->h.data[d + dim];
     R2->h.data[d] += R->h.data[d + dim];
     R2->ee = (esterr *) malloc(sizeof(esterr) * R2->fdim);
}

typedef struct rule_s {
     unsigned dim;              /* the dimensionality */
     unsigned num_points;       /* number of evaluation points */
     unsigned (*evalError)(struct rule_s *r,
			   unsigned fdim, integrand f, void *fdata,
			   const hypercube *h, esterr *ee);
     void (*destroy)(struct rule_s *r);
} rule;

static void destroy_rule(rule *r)
{
     if (r->destroy) r->destroy(r);
     free(r);
}

static region eval_region(region R, integrand f, void *fdata, rule *r)
{
     int k;
     R.splitDim = r->evalError(r, R.fdim, f, fdata, &R.h, R.ee);
     R.errmax = R.ee[0].err;
     for (k = 1; k < R.fdim; ++k)
	  if (R.ee[k].err > R.errmax) R.errmax = R.ee[k].err;
     return R;
}

/***************************************************************************/
/* Functions to loop over points in a hypercube. */

/* Based on orbitrule.cpp in HIntLib-0.0.10 */

/* ls0 returns the least-significant 0 bit of n (e.g. it returns
   0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etcetera). */

static unsigned ls0(unsigned n)
{
#if defined(__GNUC__) && \
    ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ > 3)
     return __builtin_ctz(~n); /* gcc builtin for version >= 3.4 */
#else
     const unsigned bits[256] = {
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	  0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8,
     };
     unsigned bit = 0;
     while ((n & 0xff) == 0xff) {
	  n >>= 8;
	  bit += 8;
     }
     return bit + bits[n & 0xff];
#endif
}

/**
 *  Evaluate the integral on all 2^n points (+/-r,...+/-r)
 *
 *  A Gray-code ordering is used to minimize the number of coordinate updates
 *  in p.
 */
static void evalR_Rfs(double *sum, unsigned fdim, integrand f, void *fdata, unsigned dim, double *p, const double *c, const double *r)
{
     double *val;
     unsigned i,j;
     unsigned signs = 0; /* 0/1 bit = +/- for corresponding element of r[] */

     /* We start with the point where r is ADDed in every coordinate
        (this implies signs=0). */
     for (i = 0; i < dim; ++i)
	  p[i] = c[i] + r[i];

     val = (double *) malloc(sizeof(double) * fdim);
     for (j = 0; j < fdim; ++j) sum[j] = 0.0;

     /* Loop through the points in Gray-code ordering */
     for (i = 0;; ++i) {
	  unsigned mask, d;

	  f(dim, p, fdata, fdim, val);
	  for (j = 0; j < fdim; ++j) sum[j] += val[j];

	  d = ls0(i);	/* which coordinate to flip */
	  if (d >= dim)
	       break;

	  /* flip the d-th bit and add/subtract r[d] */
	  mask = 1U << d;
	  signs ^= mask;
	  p[d] = (signs & mask) ? c[d] - r[d] : c[d] + r[d];
     }
     free(val);
}

static void evalRR0_0fs(double *sum, unsigned fdim, integrand f, void *fdata, unsigned dim, double *p, const double *c, const double *r)
{
     unsigned i, j, k;
     double *val;

     val = (double *) malloc(sizeof(double) * fdim);
     for (k = 0; k < fdim; ++k) sum[k] = 0.0;

     for (i = 0; i < dim - 1; ++i) {
	  p[i] = c[i] - r[i];
	  for (j = i + 1; j < dim; ++j) {
	       p[j] = c[j] - r[j];
	       f(dim, p, fdata, fdim, val);
	       for (k = 0; k < fdim; ++k) sum[k] += val[k];
	       p[i] = c[i] + r[i];
	       f(dim, p, fdata, fdim, val);
	       for (k = 0; k < fdim; ++k) sum[k] += val[k];
	       p[j] = c[j] + r[j];
	       f(dim, p, fdata, fdim, val);
	       for (k = 0; k < fdim; ++k) sum[k] += val[k];
	       p[i] = c[i] - r[i];
	       f(dim, p, fdata, fdim, val);
	       for (k = 0; k < fdim; ++k) sum[k] += val[k];

	       p[j] = c[j];	/* Done with j -> Restore p[j] */
	  }
	  p[i] = c[i];		/* Done with i -> Restore p[i] */
     }
     free(val);
}

static unsigned evalR0_0fs4d(unsigned fdim, integrand f, void *fdata, unsigned dim, double *p, const double *c, double *sum0, const double *r1, double *sum1, const double *r2, double *sum2)
{
     double maxdiff = 0;
     unsigned i, j, dimDiffMax = 0;
     double *val, *val0, *D;
     double ratio = r1[0] / r2[0];

     val = (double *) malloc(sizeof(double) * fdim * 3);
     val0 = val + fdim; D = val0 + fdim;

     ratio *= ratio;
     f(dim, p, fdata, fdim, val0);
     for (j = 0; j < fdim; ++j) sum0[j] += val0[j];

     for (i = 0; i < dim; i++) {
	  double diff;

	  p[i] = c[i] - r1[i];
	  f(dim, p, fdata, fdim, val);
	  for (j = 0; j < fdim; ++j) {
	       sum1[j] += val[j];
	       D[j] = val[j];
	  }

	  p[i] = c[i] + r1[i];
	  f(dim, p, fdata, fdim, val);
	  for (j = 0; j < fdim; ++j) {
	       sum1[j] += val[j];
	       D[j] += val[j] - 2*val0[j];
	  }

	  p[i] = c[i] - r2[i];
	  f(dim, p, fdata, fdim, val);
	  for (j = 0; j < fdim; ++j) {
	       sum2[j] += val[j];
	       D[j] -= ratio * val[j];
	  }

	  p[i] = c[i] + r2[i];
	  f(dim, p, fdata, fdim, val);
	  diff = 0;
	  for (j = 0; j < fdim; ++j) {
	       sum2[j] += val[j];
	       diff += fabs(D[j] - ratio * (val[j] - 2*val0[j]));
	  }

	  p[i] = c[i];

	  if (diff > maxdiff) {
	       maxdiff = diff;
	       dimDiffMax = i;
	  }
     }
     free(val);
     return dimDiffMax;
}

#define num0_0(dim) (1U)
#define numR0_0fs(dim) (2 * (dim))
#define numRR0_0fs(dim) (2 * (dim) * (dim-1))
#define numR_Rfs(dim) (1U << (dim))

/***************************************************************************/
/* Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
   cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
   and A. A. Malik.  See:

         A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
         symmetric numerical integration rules," SIAM
         J. Numer. Anal. 20 (3), 580-588 (1983).
*/

typedef struct {
     rule parent;

     /* temporary arrays of length dim */
     double *widthLambda, *widthLambda2, *p;

     /* dimension-dependent constants */
     double weight1, weight3, weight5;
     double weightE1, weightE3;
} rule75genzmalik;

#define real(x) ((double)(x))
#define to_int(n) ((int)(n))

static int isqr(int x)
{
     return x * x;
}

static void destroy_rule75genzmalik(rule *r_)
{
     rule75genzmalik *r = (rule75genzmalik *) r_;
     free(r->p);
}

static unsigned rule75genzmalik_evalError(rule *r_, unsigned fdim, integrand f, void *fdata, const hypercube *h, esterr *ee)
{
     /* lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19) */
     const double lambda2 = 0.3585685828003180919906451539079374954541;
     const double lambda4 = 0.9486832980505137995996680633298155601160;
     const double lambda5 = 0.6882472016116852977216287342936235251269;
     const double weight2 = 980. / 6561.;
     const double weight4 = 200. / 19683.;
     const double weightE2 = 245. / 486.;
     const double weightE4 = 25. / 729.;

     rule75genzmalik *r = (rule75genzmalik *) r_;
     unsigned i, j, dimDiffMax, dim = r_->dim;
     double *sums, *sum1, *sum2, *sum3, *sum4, *sum5;
     const double *center = h->data;
     const double *halfwidth = h->data + dim;

     sums = (double *) malloc(sizeof(double) * fdim * 5);
     sum1 = sums; sum2 = sum1 + fdim; sum3 = sum2 + fdim; sum4 = sum3 + fdim;
     sum5 = sum4 + fdim;
     for (j = 0; j < fdim; ++j) sum1[j] = sum2[j] = sum3[j] = 0.0;

     for (i = 0; i < dim; ++i)
	  r->p[i] = center[i];

     for (i = 0; i < dim; ++i)
	  r->widthLambda2[i] = halfwidth[i] * lambda2;
     for (i = 0; i < dim; ++i)
	  r->widthLambda[i] = halfwidth[i] * lambda4;

     /* Evaluate function in the center, in f(lambda2,0,...,0) and
        f(lambda3=lambda4, 0,...,0).  Estimate dimension with largest error */
     dimDiffMax = evalR0_0fs4d(fdim, f, fdata, dim, r->p, center, sum1, r->widthLambda2, sum2, r->widthLambda, sum3);

     /* Calculate sum4 for f(lambda4, lambda4, 0, ...,0) */
     evalRR0_0fs(sum4, fdim, f, fdata, dim, r->p, center, r->widthLambda);

     /* Calculate sum5 for f(lambda5, lambda5, ..., lambda5) */
     for (i = 0; i < dim; ++i)
	  r->widthLambda[i] = halfwidth[i] * lambda5;
     evalR_Rfs(sum5, fdim, f, fdata, dim, r->p, center, r->widthLambda);

     /* Calculate fifth and seventh order results */
     
     for (j = 0; j < fdim; ++j) {
	  double result, res5th;
	  result = h->vol * (r->weight1 * sum1[j] + weight2 * sum2[j] + r->weight3 * sum3[j] + weight4 * sum4[j] + r->weight5 * sum5[j]);
	  res5th = h->vol * (r->weightE1 * sum1[j] + weightE2 * sum2[j] + r->weightE3 * sum3[j] + weightE4 * sum4[j]);
	  
	  ee[j].val = result;
	  ee[j].err = fabs(res5th - result);
     }

     free(sums);
     return dimDiffMax;
}

static rule *make_rule75genzmalik(unsigned dim)
{
     rule75genzmalik *r;

     if (dim < 2) return 0; /* this rule does not support 1d integrals */

     /* Because of the use of a bit-field in evalR_Rfs, we are limited
	to be < 32 dimensions (or however many bits are in unsigned).
	This is not a practical limitation...long before you reach
	32 dimensions, the Genz-Malik cubature becomes excruciatingly
	slow and is superseded by other methods (e.g. Monte-Carlo). */
     if (dim >= sizeof(unsigned) * 8) return 0;

     r = (rule75genzmalik *) malloc(sizeof(rule75genzmalik));
     if (!r) return 0;
     r->parent.dim = dim;

     r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim)))
		   / real(19683));
     r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
     r->weight5 = real(6859) / real(19683) / real(1U << dim);
     r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim)))
		    / real(729));
     r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);

     r->p = (double *) malloc(sizeof(double) * dim * 3);
     if (!r->p) { free(r); return 0; }
     r->widthLambda = r->p + dim;
     r->widthLambda2 = r->p + 2 * dim;

     r->parent.num_points = num0_0(dim) + 2 * numR0_0fs(dim)
	  + numRR0_0fs(dim) + numR_Rfs(dim);

     r->parent.evalError = rule75genzmalik_evalError;
     r->parent.destroy = destroy_rule75genzmalik;

     return (rule *) r;
}

/***************************************************************************/
/* 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
   GNU GSL (which in turn is based on QUADPACK). */

static unsigned rule15gauss_evalError(rule *r,
				      unsigned fdim, integrand f, void *fdata,
				      const hypercube *h, esterr *ee)
{
     /* Gauss quadrature weights and kronrod quadrature abscissae and
	weights as evaluated with 80 decimal digit arithmetic by
	L. W. Fullerton, Bell Labs, Nov. 1981. */
     const unsigned n = 8;
     const double xgk[8] = {  /* abscissae of the 15-point kronrod rule */
	  0.991455371120812639206854697526329,
	  0.949107912342758524526189684047851,
	  0.864864423359769072789712788640926,
	  0.741531185599394439863864773280788,
	  0.586087235467691130294144838258730,
	  0.405845151377397166906606412076961,
	  0.207784955007898467600689403773245,
	  0.000000000000000000000000000000000
	  /* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
	     xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
     };
     static const double wg[4] = {  /* weights of the 7-point gauss rule */
	  0.129484966168869693270611432679082,
	  0.279705391489276667901467771423780,
	  0.381830050505118944950369775488975,
	  0.417959183673469387755102040816327
     };
     static const double wgk[8] = { /* weights of the 15-point kronrod rule */
	  0.022935322010529224963732008058970,
	  0.063092092629978553290700663189204,
	  0.104790010322250183839876322541518,
	  0.140653259715525918745189590510238,
	  0.169004726639267902826583426598550,
	  0.190350578064785409913256402421014,
	  0.204432940075298892414161999234649,
	  0.209482141084727828012999174891714
     };

     const double center = h->data[0];
     const double halfwidth = h->data[1];

     double *f_center, *val, *fv1, *fv2;
     double *result_gauss, *result_kronrod, *result_abs;
     unsigned j, k;
     double *scratch = (double *) malloc(sizeof(double) * fdim * 19);
     f_center = scratch; val = scratch + fdim;
     result_gauss = val + fdim; result_kronrod = result_gauss + fdim;
     result_abs = result_kronrod + fdim;
     fv1 = result_abs + fdim; fv2 = fv1 + 7*fdim;

     f(1, &center, fdata, fdim, f_center);
     for (k = 0; k < fdim; ++k) {
	  result_gauss[k] = f_center[k] * wg[n/2 - 1];
	  result_kronrod[k] = f_center[k] * wgk[n - 1];
	  result_abs[k] = fabs(result_kronrod[k]);
     }
     
     for (j = 0; j < (n - 1) / 2; ++j) {
	  int j2 = 2*j + 1;
	  double x, w = halfwidth * xgk[j2];

	  x = center - w; 
	  f(1, &x, fdata, fdim, val);
	  for (k = 0; k < fdim; ++k) {
	       fv1[j2*fdim + k] = val[k];
	       result_gauss[k] += wg[j] * val[k];
	       result_kronrod[k] += wgk[j2] * val[k];
	       result_abs[k] += wgk[j2] * fabs(val[k]);
	  }

	  x = center + w; 
	  f(1, &x, fdata, fdim, val);
	  for (k = 0; k < fdim; ++k) {
	       fv2[j2*fdim + k] = val[k];
	       result_gauss[k] += wg[j] * val[k];
	       result_kronrod[k] += wgk[j2] * val[k];
	       result_abs[k] += wgk[j2] * fabs(val[k]);
	  }
     }

     for (j = 0; j < n/2; ++j) {
	  int j2 = 2*j;
	  double x, w = halfwidth * xgk[j2];

	  x = center - w; 
	  f(1, &x, fdata, fdim, val);
	  for (k = 0; k < fdim; ++k) {
	       fv1[j2*fdim + k] = val[k];
	       result_kronrod[k] += wgk[j2] * val[k];
	       result_abs[k] += wgk[j2] * fabs(val[k]);
	  }

	  x = center + w; 
	  f(1, &x, fdata, fdim, val);
	  for (k = 0; k < fdim; ++k) {
	       fv2[j2*fdim + k] = val[k];
	       result_kronrod[k] += wgk[j2] * val[k];
	       result_abs[k] += wgk[j2] * fabs(val[k]);
	  }
     }

     /* compute result and error estimate: */
     for (k = 0; k < fdim; ++k) {
	  double result_asc, mean, err;

	  ee[k].val = result_kronrod[k] * halfwidth;

	  mean = result_kronrod[k] * 0.5;
	  result_asc = wgk[n - 1] * fabs(f_center[k] - mean);
	  for (j = 0; j < n - 1; ++j)
	       result_asc += wgk[j] * (fabs(fv1[j*fdim+k]-mean) 
				       + fabs(fv2[j*fdim+k]-mean));
	  err = fabs(result_kronrod[k] - result_gauss[k]) * halfwidth;
	  result_abs[k] *= halfwidth;
	  result_asc *= halfwidth;
	  if (result_asc != 0 && err != 0) {
	       double scale = pow((200 * err / result_asc), 1.5);
	       if (scale < 1)
		    err = result_asc * scale;
	       else
		    err = result_asc;
	  }
	  if (result_abs[k] > DBL_MIN / (50 * DBL_EPSILON)) {
	       double min_err = 50 * DBL_EPSILON * result_abs[k];
	       if (min_err > err)
		    err = min_err;
	  }
	  ee[k].err = err;
     }

     free(scratch);
     return 0; /* no choice but to divide 0th dimension */
}

static rule *make_rule15gauss(unsigned dim)
{
     rule *r;
     if (dim != 1) return 0; /* this rule is only for 1d integrals */
     r = (rule *) malloc(sizeof(rule));
     if (!r) return 0;
     r->dim = dim;
     r->num_points = 15;
     r->evalError = rule15gauss_evalError;
     r->destroy = 0;
     return r;
}

/***************************************************************************/
/* binary heap implementation (ala _Introduction to Algorithms_ by
   Cormen, Leiserson, and Rivest), for use as a priority queue of
   regions to integrate. */

typedef region heap_item;
#define KEY(hi) ((hi).errmax)

typedef struct {
     unsigned n, nalloc;
     heap_item *items;
     unsigned fdim;
     esterr *ee; /* array of length fdim of the total integrand & error */
} heap;

static void heap_resize(heap *h, unsigned nalloc)
{
     h->nalloc = nalloc;
     h->items = (heap_item *) realloc(h->items, sizeof(heap_item) * nalloc);
}

static heap heap_alloc(unsigned nalloc, unsigned fdim)
{
     heap h;
     unsigned i;
     h.n = 0;
     h.nalloc = 0;
     h.items = 0;
     h.fdim = fdim;
     h.ee = (esterr *) malloc(sizeof(esterr) * fdim);
     for (i = 0; i < fdim; ++i) h.ee[i].val = h.ee[i].err = 0;
     heap_resize(&h, nalloc);
     return h;
}

/* note that heap_free does not deallocate anything referenced by the items */
static void heap_free(heap *h)
{
     h->n = 0;
     heap_resize(h, 0);
     h->fdim = 0;
     free(h->ee);
}

static void heap_push(heap *h, heap_item hi)
{
     int insert;
     unsigned i, fdim = h->fdim;

     for (i = 0; i < fdim; ++i) {
	  h->ee[i].val += hi.ee[i].val;
	  h->ee[i].err += hi.ee[i].err;
     }
     insert = h->n;
     if (++(h->n) > h->nalloc)
	  heap_resize(h, h->n * 2);

     while (insert) {
	  int parent = (insert - 1) / 2;
	  if (KEY(hi) <= KEY(h->items[parent]))
	       break;
	  h->items[insert] = h->items[parent];
	  insert = parent;
     }
     h->items[insert] = hi;
}

static heap_item heap_pop(heap *h)
{
     heap_item ret;
     int i, n, child;

     if (!(h->n)) {
	  fprintf(stderr, "attempted to pop an empty heap\n");
	  exit(EXIT_FAILURE);
     }

     ret = h->items[0];
     h->items[i = 0] = h->items[n = --(h->n)];
     while ((child = i * 2 + 1) < n) {
	  int largest;
	  heap_item swap;

	  if (KEY(h->items[child]) <= KEY(h->items[i]))
	       largest = i;
	  else
	       largest = child;
	  if (++child < n && KEY(h->items[largest]) < KEY(h->items[child]))
	       largest = child;
	  if (largest == i)
	       break;
	  swap = h->items[i];
	  h->items[i] = h->items[largest];
	  h->items[i = largest] = swap;
     }

     {
	  unsigned i, fdim = h->fdim;
	  for (i = 0; i < fdim; ++i) {
	       h->ee[i].val -= ret.ee[i].val;
	       h->ee[i].err -= ret.ee[i].err;
	  }
     }
     return ret;
}

/***************************************************************************/

/* adaptive integration, analogous to adaptintegrator.cpp in HIntLib */

static int ruleadapt_integrate(rule *r, unsigned fdim, integrand f, void *fdata, const hypercube *h, unsigned maxEval, double reqAbsError, double reqRelError, double *val, double *err)
{
     unsigned maxIter;		/* maximum number of adaptive subdivisions */
     heap regions;
     unsigned i, j;
     int status = -1; /* = ERROR */

     if (maxEval) {
	  if (r->num_points > maxEval)
	       return status; /* ERROR */
	  maxIter = (maxEval - r->num_points) / (2 * r->num_points);
     }
     else
	  maxIter = UINT_MAX;

     regions = heap_alloc(1, fdim);

     heap_push(&regions, eval_region(make_region(h, fdim), f, fdata, r));
     /* another possibility is to specify some non-adaptive subdivisions: 
	if (initialRegions != 1)
	   partition(h, initialRegions, EQUIDISTANT, &regions, f,fdata, r); */

     for (i = 0; i < maxIter; ++i) {
	  region R, R2;
	  for (j = 0; j < fdim && (regions.ee[j].err <= reqAbsError
				   || relError(regions.ee[j]) <= reqRelError);
	       ++j) ;
	  if (j == fdim) {
	       status = 0; /* converged! */
	       break;
	  }
	  R = heap_pop(&regions); /* get worst region */
	  cut_region(&R, &R2);
	  heap_push(&regions, eval_region(R, f, fdata, r));
	  heap_push(&regions, eval_region(R2, f, fdata, r));
     }

     /* re-sum integral and errors */
     for (j = 0; j < fdim; ++j) val[j] = err[j] = 0;  
     for (i = 0; i < regions.n; ++i) {
	  for (j = 0; j < fdim; ++j) { 
	       val[j] += regions.items[i].ee[j].val;
	       err[j] += regions.items[i].ee[j].err;
	  }
	  destroy_region(&regions.items[i]);
     }

     /* printf("regions.nalloc = %d\n", regions.nalloc); */
     heap_free(&regions);

     return status;
}

int adapt_integrate(unsigned fdim, integrand f, void *fdata, 
		    unsigned dim, const double *xmin, const double *xmax, 
		    unsigned maxEval, double reqAbsError, double reqRelError, 
		    double *val, double *err)
{
     rule *r;
     hypercube h;
     int status;
     
     if (fdim == 0) /* nothing to do */ return 0;
     if (dim == 0) { /* trivial integration */
	  f(0, xmin, fdata, fdim, val);
	  *err = 0;
	  return 0;
     }
     r = dim == 1 ? make_rule15gauss(dim) : make_rule75genzmalik(dim);
     if (!r) { 
	  unsigned i;
	  for (i = 0; i < fdim; ++i) {
	       val[i] = 0;
	       err[i] = HUGE_VAL; 
	  }
	  return -2; /* ERROR */
     }
     h = make_hypercube_range(dim, xmin, xmax);
     status = ruleadapt_integrate(r, fdim, f, fdata, &h,
				  maxEval, reqAbsError, reqRelError,
				  val, err);
     destroy_hypercube(&h);
     destroy_rule(r);
     return status;
}

/***************************************************************************/

/* Compile with -DTEST_INTEGRATOR for a self-contained test program.
   
   Usage: ./integrator <dim> <tol> <integrand> <maxeval>

   where <dim> = # dimensions, <tol> = relative tolerance,
   <integrand> is either 0/1/2 for the three test integrands (see below),
   and <maxeval> is the maximum # function evaluations (0 for none).
*/
   
#ifdef TEST_INTEGRATOR

int count = 0;
unsigned integrand_fdim = 0;
int *which_integrand = NULL;
const double radius = 0.50124145262344534123412; /* random */

/* Simple constant function */
double
fconst (double x[], size_t dim, void *params)
{
  return 1;
}

/*** f0, f1, f2, and f3 are test functions from the Monte-Carlo
     integration routines in GSL 1.6 (monte/test.c).  Copyright (c)
     1996-2000 Michael Booth, GNU GPL. ****/

/* Simple product function */
double f0 (unsigned dim, const double *x, void *params)
{
     double prod = 1.0;
     unsigned int i;
     for (i = 0; i < dim; ++i)
	  prod *= 2.0 * x[i];
     return prod;
}

/* Gaussian centered at 1/2. */
double f1 (unsigned dim, const double *x, void *params)
{
     double a = *(double *)params;
     double sum = 0.;
     unsigned int i;
     for (i = 0; i < dim; i++) {
	  double dx = x[i] - 0.5;
	  sum += dx * dx;
     }
     return (pow (M_2_SQRTPI / (2. * a), (double) dim) *
	     exp (-sum / (a * a)));
}

/* double gaussian */
double f2 (unsigned dim, const double *x, void *params)
{
     double a = *(double *)params;
     double sum1 = 0.;
     double sum2 = 0.;
     unsigned int i;
     for (i = 0; i < dim; i++) {
	  double dx1 = x[i] - 1. / 3.;
	  double dx2 = x[i] - 2. / 3.;
	  sum1 += dx1 * dx1;
	  sum2 += dx2 * dx2;
     }
     return 0.5 * pow (M_2_SQRTPI / (2. * a), dim) 
	  * (exp (-sum1 / (a * a)) + exp (-sum2 / (a * a)));
}

/* Tsuda's example */
double f3 (unsigned dim, const double *x, void *params)
{
     double c = *(double *)params;
     double prod = 1.;
     unsigned int i;
     for (i = 0; i < dim; i++)
	  prod *= c / (c + 1) * pow((c + 1) / (c + x[i]), 2.0);
     return prod;
}

/* test integrand from W. J. Morokoff and R. E. Caflisch, "Quasi=
   Monte Carlo integration," J. Comput. Phys 122, 218-230 (1995).
   Designed for integration on [0,1]^dim, integral = 1. */
static double morokoff(unsigned dim, const double *x, void *params)
{
     double p = 1.0 / dim;
     double prod = pow(1 + p, dim);
     unsigned int i;
     for (i = 0; i < dim; i++)
	  prod *= pow(x[i], p);
     return prod;
}

/*** end of GSL test functions ***/

void f_test(unsigned dim, const double *x, void *data_,
	    unsigned fdim, double *retval)
{
     double val;
     unsigned i, j;
     ++count;
     (void) data_; /* not used */
     for (j = 0; j < fdim; ++j) {
     double fdata = which_integrand[j] == 6 ? (1.0+sqrt (10.0))/9.0 : 0.1;
     switch (which_integrand[j]) {
	 case 0: /* simple smooth (separable) objective: prod. cos(x[i]). */
	      val = 1;
	      for (i = 0; i < dim; ++i)
		   val *= cos(x[i]);
	      break;
	 case 1: { /* integral of exp(-x^2), rescaled to (0,infinity) limits */
	      double scale = 1.0;
	      val = 0;
	      for (i = 0; i < dim; ++i) {
		   double z = (1 - x[i]) / x[i];
		   val += z * z;
		   scale *= M_2_SQRTPI / (x[i] * x[i]);
	      }
	      val = exp(-val) * scale;
	      break;
	 }
	 case 2: /* discontinuous objective: volume of hypersphere */
	      val = 0;
	      for (i = 0; i < dim; ++i)
		   val += x[i] * x[i];
	      val = val < radius * radius;
	      break;
	 case 3:
	      val = f0(dim, x, &fdata);
	      break;
	 case 4:
	      val = f1(dim, x, &fdata);
	      break;
	 case 5:
	      val = f2(dim, x, &fdata);
	      break;
	 case 6:
	      val = f3(dim, x, &fdata);
	      break;
	 case 7:
	      val = morokoff(dim, x, &fdata);
	      break;
	 default:
	      fprintf(stderr, "unknown integrand %d\n", which_integrand[j]);
	      exit(EXIT_FAILURE);
     }
     /* if (count < 100) printf("%d: f(%g, ...) = %g\n", count, x[0], val); */
     retval[j] = val;
     }
}

/* surface area of n-dimensional unit hypersphere */
static double S(unsigned n)
{
     double val;
     int fact = 1;
     if (n % 2 == 0) { /* n even */
	  val = 2 * pow(M_PI, n * 0.5);
	  n = n / 2;
	  while (n > 1) fact *= (n -= 1);
	  val /= fact;
     }
     else { /* n odd */
	  val = (1 << (n/2 + 1)) * pow(M_PI, n/2);
	  while (n > 2) fact *= (n -= 2);
	  val /= fact;
     }
     return val;
}

static double exact_integral(int which, unsigned dim, const double *xmax) {
     unsigned i;
     double val;
     switch(which) {
	 case 0:
	      val = 1;
	      for (i = 0; i < dim; ++i)
		   val *= sin(xmax[i]);
	      break;
	 case 2:
	      val = dim == 0 ? 1 : S(dim) * pow(radius * 0.5, dim) / dim;
	      break;
	 default:
	      val = 1.0;
     }
     return val;
}

#include <ctype.h>
int main(int argc, char **argv)
{
     double *xmin, *xmax;
     double tol, *val, *err;
     unsigned i, dim, maxEval;

     dim = argc > 1 ? atoi(argv[1]) : 2;
     tol = argc > 2 ? atof(argv[2]) : 1e-2;
     maxEval = argc > 4 ? atoi(argv[4]) : 0;
     
     /* parse: e.g. "x/y/z" is treated as fdim = 3, which_integrand={x,y,z} */
     if (argc <= 3) {
	  integrand_fdim = 1;
	  which_integrand = (int *) malloc(sizeof(int) * integrand_fdim);
	  which_integrand[0] = 0; /* default */
     }
     else {
	  unsigned j = 0;
	  integrand_fdim = 1;
	  for (i = 0; argv[3][i]; ++i) if (argv[3][i] == '/') ++integrand_fdim;
	  if (!integrand_fdim) {
	       fprintf(stderr, "invalid which_integrand \"%s\"", argv[3]);
	       exit(EXIT_FAILURE);
	  }
	  which_integrand = (int *) malloc(sizeof(int) * integrand_fdim);
	  which_integrand[0] = 0;
	  for (i = 0; argv[3][i]; ++i) {
	       if (argv[3][i] == '/')
		    which_integrand[++j] = 0;
	       else if (isdigit(argv[3][i]))
		    which_integrand[j] = 
			 which_integrand[j]*10 + argv[3][i] - '0';
	       else {
		    fprintf(stderr, "invalid which_integrand \"%s\"", argv[3]);
		    exit(EXIT_FAILURE);
	       }
	  }
     }
     val = (double *) malloc(sizeof(double) * integrand_fdim);
     err = (double *) malloc(sizeof(double) * integrand_fdim);

     xmin = (double *) malloc(dim * sizeof(double));
     xmax = (double *) malloc(dim * sizeof(double));
     for (i = 0; i < dim; ++i) {
	  xmin[i] = 0;
	  xmax[i] = 1;
     }

     printf("%u-dim integral, tolerance = %g\n", dim, tol);
     adapt_integrate(integrand_fdim, f_test, NULL, 
		     dim, xmin, xmax, 
		     maxEval, 0, tol, val, err);
     for (i = 0; i < integrand_fdim; ++i) {
	  printf("integrand %d: integral = %g, est err = %g, true err = %g\n", 
		 which_integrand[i], val[i], err[i], 
		 fabs(val[i] - exact_integral(which_integrand[i], dim, xmax)));
     }
     printf("#evals = %d\n", count);

     free(xmax);
     free(xmin);
     free(err);
     free(val);
     free(which_integrand);

     return 0;
}

#endif
