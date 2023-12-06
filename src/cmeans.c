/* C code for (weighted) fuzzy c-means, rewritten from scratch by KH. */

#include <stdlib.h>
#include <math.h>
#include <R.h>

/* Enhance readability of matrix-subscripting for matrices stored in
   row-major order. */
#define MSUB(x, i, j, n)	x[(i) + (n) * (j)]

static double *d;
static double *dwrk, *dwrk_x, *dwrk_w;
static int *iwrk;

static void
cmeans_setup(int nr_x, int nr_p, int dist)
{
    int len_u_d = nr_x * nr_p;
    
    d = (double *) R_alloc(len_u_d, sizeof(double));
    if(dist == 1) {
	/* Needed for weighted medians. */
	dwrk_x = (double *) R_alloc(nr_x, sizeof(double));
	dwrk_w = (double *) R_alloc(nr_x, sizeof(double));
	dwrk = (double *) R_alloc(nr_x, sizeof(double));
	iwrk = (int *) R_alloc(nr_x, sizeof(int));
    }
}

/*
static void
cmeans_copy_vector(double *from, double *to, int len)
{
    int i;
    for(i = 0; i < len; i++)
	to[i] = from[i];
}

static double
cmeans_delta_old_new(double *old, double *new, int len)
{
    int i;
    double sum = 0;
    for(i = 0; i < len; i++)
	sum += fabs(new[i] - old[i]);
    return(sum / len);
}
*/

static int
cmeans_sign(double x)
{
    if(x == 0) return(0);
    return((x > 0) ? 1 : -1);
}

static double
cmeans_weighted_median(double *x, double *w, int len)
{
    int i;
    double sum, val, marg, mval, cumsum_w, cumsum_w_x;

    /* Sort x. */
    for(i = 0; i < len; i++)
	iwrk[i] = i;
    rsort_with_index(x, iwrk, len);
    
    /* Permute w using iwrk, and normalize. */
    sum = 0;    
    for(i = 0; i < len; i++) {
	dwrk[i] = w[iwrk[i]];
	sum += dwrk[i];
    }
    for(i = 0; i < len; i++) {
	w[i] = dwrk[i] / sum;
    }

    cumsum_w = cumsum_w_x = 0;
    mval = R_PosInf;
    marg = *x;			/* -Wall */
    for(i = 0; i < len; i++) {
	cumsum_w += w[i];
	cumsum_w_x += w[i] * x[i];
	val = x[i] * (cumsum_w - .5) - cumsum_w_x;
	if(val < mval) {
	    marg = x[i];
	    mval = val;
	}
    }

    return(marg);
}

/* Update the dissimilarities (between objects and prototypes) for a
 * single object (i.e., a single row of the dissimilarity matrix. */
static void
ufcl_dissimilarities(double *x, double *p,
		     int nr_x, int nc, int nr_p,
		     int dist, int ix, double *d)
{
    int ip, j;
    double sum, v;
    
    for(ip = 0; ip < nr_p; ip++) {
	sum = 0;
	for(j = 0; j < nc; j++) {
	    v = MSUB(x, ix, j, nr_x) - MSUB(p, ip, j, nr_p);
	    if(dist == 0)
		sum += v * v;
	    else if(dist == 1)
		sum += fabs(v);
	}
	MSUB(d, ix, ip, nr_x) = sum;
    }
}

static void
cmeans_dissimilarities(double *x, double *p,
		       int nr_x, int nc, int nr_p,
		       int dist, double *d)
{
    int ix;
    
    for(ix = 0; ix < nr_x; ix++) {
	/* Loop over all objects ... */
	ufcl_dissimilarities(x, p, nr_x, nc, nr_p, dist, ix, d);
    }
}

/* Update the memberships for a single object (i.e., a single row of the
 * membership matrix.) */
static void
ufcl_memberships(double *d, int nr_x, int nr_p,
		 double exponent, int ix,
		 double *u)
{
    int ip, n_of_zeroes;
    double sum, v;

    n_of_zeroes = 0;
    for(ip = 0; ip < nr_p; ip++) {
	if(MSUB(d, ix, ip, nr_x) == 0)
	    n_of_zeroes++;
    }
    if(n_of_zeroes > 0) {
	v = 1 / n_of_zeroes;
	for(ip = 0; ip < nr_p; ip++)
	    MSUB(u, ix, ip, nr_x) =
		((MSUB(d, ix, ip, nr_x) == 0) ? v : 0);
    }
    else {
	/* Use the assumption that in general, pow() is more
	 * expensive than subscripting. */
	sum = 0;
	for(ip = 0; ip < nr_p; ip++) {
	    v = 1 / pow(MSUB(d, ix, ip, nr_x), exponent);
	    sum += v;
	    MSUB(u, ix, ip, nr_x) = v;
	}
	for(ip = 0; ip < nr_p; ip++)
	    MSUB(u, ix, ip, nr_x) /= sum;
    }
}

static void
cmeans_memberships(double *d,
		   int nr_x, int nr_p,
		   double exponent, double *u)
{
    int ix;

    for(ix = 0; ix < nr_x; ix++) {
	/* Loop over all objects ... */
	ufcl_memberships(d, nr_x, nr_p, exponent, ix, u);
    }
}

static void
cmeans_prototypes(double *x, double *u, double *w,
		  int nr_x, int nc, int nr_p,
		  double f, int dist, double *p)
{
    int ix, ip, j;
    double sum, v;

    if(dist == 0) {
	/* Euclidean: weighted means. */
	for(ip = 0; ip < nr_p; ip++) {
	    for(j = 0; j < nc; j++)
		MSUB(p, ip, j, nr_p) = 0;
	    sum = 0;
	    for(ix = 0; ix < nr_x; ix++) {
		v = w[ix] * pow(MSUB(u, ix, ip, nr_x), f);
		sum += v;
		for(j = 0; j < nc; j++)
		    MSUB(p, ip, j, nr_p) += v * MSUB(x, ix, j, nr_x);
	    }
	    for(j = 0; j < nc; j++)
		MSUB(p, ip, j, nr_p) /= sum;
	}
    }
    else {
	/* Manhattan: weighted medians. */
	for(ip = 0; ip < nr_p; ip++)
	    for(j = 0; j < nc; j++) {
		for(ix = 0; ix < nr_x; ix++) {
		    dwrk_x[ix] = MSUB(x, ix, j, nr_x);
		    dwrk_w[ix] = w[ix] * pow(MSUB(u, ix, ip, nr_x), f);
		}
		MSUB(p, ip, j, nr_p) =
		    cmeans_weighted_median(dwrk_x, dwrk_w, nr_x);
	    }
    }
}

static double
cmeans_error_fn(double *u, double *d, double *w,
		int nr_x, int nr_p, double f)
{
    int ix, ip;    
    double sum;

    sum = 0;
    for(ix = 0; ix < nr_x; ix++)    
	for(ip = 0; ip < nr_p; ip++)
	    sum += w[ix] * pow(MSUB(u, ix, ip, nr_x), f)
		* MSUB(d, ix, ip, nr_x);
    return(sum);
}

void
cmeans(double *x, int *nr_x, int *nc, double *p, int *nr_p, double *w,
       double *f, int *dist, int *itermax, double *reltol, int *verbose,
       double *u, double *ermin, int *iter)
{
    double exponent = 1 / (*f - 1);
    double old_value, new_value;

    cmeans_setup(*nr_x, *nr_p, *dist);
    
    cmeans_dissimilarities(x, p, *nr_x, *nc, *nr_p, *dist, d);
    cmeans_memberships(d, *nr_x, *nr_p, exponent, u);
    old_value = new_value = cmeans_error_fn(u, d, w, *nr_x, *nr_p, *f);
    
    *iter = 0;
    while((*iter)++ < *itermax) {
	cmeans_prototypes(x, u, w, *nr_x, *nc, *nr_p, *f, *dist, p);
	cmeans_dissimilarities(x, p, *nr_x, *nc, *nr_p, *dist, d);
	cmeans_memberships(d, *nr_x, *nr_p, exponent, u);
	new_value = cmeans_error_fn(u, d, w, *nr_x, *nr_p, *f);
	if(fabs(old_value - new_value) < *reltol * (old_value + *reltol)) {
	    if(*verbose)
		Rprintf("Iteration: %3d converged, Error: %13.10f\n",
			*iter, new_value);
	    break;
	}
	else {
	    if(*verbose) {
		*ermin = cmeans_error_fn(u, d, w, *nr_x, *nr_p, *f);
		Rprintf("Iteration: %3d, Error: %13.10f\n",
			*iter, new_value);
	    }
	    old_value = new_value;
	}
    }

    *ermin = new_value;
}

/* Update prototypes based on a single object. */
static void
ufcl_prototypes(double *x, double *u, double *w,
		int nr_x, int nc, int nr_p,
		double f, int dist, double lrate, int ix, double *p)
{
    int ip, j;
    double grad;

    for(ip = 0; ip < nr_p; ip++) {
	for(j = 0; j < nc; j++) {
	    grad = MSUB(x, ix, j, nr_x) - MSUB(p, ip, j, nr_p);
	    if(dist == 1)
		grad = cmeans_sign(grad);
	    MSUB(p, ip, j, nr_p) +=
		lrate * w[ix] * pow(MSUB(u, ix, ip, nr_x), f) * grad;
	}
    }
}

void
ufcl(double *x, int *nr_x, int *nc, double *p, int *nr_p, double *w,
     double *f, int *dist, int *itermax, double *reltol, int *verbose,
     double *rate_par,
     double *u, double *ermin, int *iter)
{
    double exponent = 1 / (*f - 1);
    double old_value, new_value;

    int ix;
    double lrate;

    cmeans_setup(*nr_x, *nr_p, 0);
    
    /* Need some starting values ... */
    cmeans_dissimilarities(x, p, *nr_x, *nc, *nr_p, *dist, d);
    cmeans_memberships(d, *nr_x, *nr_p, exponent, u);
    old_value = new_value = cmeans_error_fn(u, d, w, *nr_x, *nr_p, *f);
    
    *iter = 0;
    while((*iter)++ < *itermax) {
	/* Turns out that sampling the objects is a bad idea ... */
	lrate = *rate_par * (1 - (double) *iter / *itermax);
	for(ix = 0; ix < *nr_x; ix++) {
	    ufcl_dissimilarities(x, p, *nr_x, *nc, *nr_p, *dist, ix, d);
	    ufcl_memberships(d, *nr_x, *nr_p, exponent, ix, u);
	    ufcl_prototypes(x, u, w, *nr_x, *nc, *nr_p, *f, *dist,
			    lrate, ix, p);
	}
	new_value = cmeans_error_fn(u, d, w, *nr_x, *nr_p, *f);
	if(fabs(old_value - new_value) < *reltol * (old_value + *reltol)) {
	    if(*verbose)
		Rprintf("Iteration: %3d converged, Error: %13.10f\n",
			*iter, new_value);
	    break;
	}
	else {
	    if(*verbose) {
		*ermin = cmeans_error_fn(u, d, w, *nr_x, *nr_p, *f);
		Rprintf("Iteration: %3d, Error: %13.10f\n",
			*iter, new_value);
	    }
	    old_value = new_value;
	}
    }

    *ermin = new_value;
}	
