#include <stdio.h>
#include <stdlib.h>
#include <svm.h>

/*
 * svm_model
 */
struct svm_model
{
	int n;			 /* number of SVs */
	double *sv_coef;	 /* sv_coef[i] is the coefficient of SV[i] */
	struct svm_node ** SV;	 /* SVs */
	double rho;		 /* the constant in the decision function */

	struct svm_parameter param;	 /* parameter */

	int free_sv;		 /* XXX: 1 if svm_model is created by svm_load_model */
				      /* 0 if svm_model is created by svm_train */
};

struct svm_node ** sparsify (double *x, int r, int c)
{
    struct svm_node** sparse;
    int         i, ii, count;
    
    sparse = (struct svm_node **) malloc (r * sizeof(struct svm_node *));
    for (i = 0; i < r; i++) {
	/* determine nr. of non-zero elements */
	for (count = ii = 0; ii < c; ii++)
	    if (x[i * c + ii] != 0) count++;

	/* allocate memory for column elements */
	sparse[i] = (struct svm_node *) malloc ((count + 1) * sizeof(struct svm_node));

	/* set column elements */
	for (count = ii = 0; ii < c; ii++)
	    if (x[i * c + ii] != 0) {
		sparse[i][count].index = ii;
		sparse[i][count].value = x[i * c + ii];
		count++;
	    }

	/* set termination element */
	sparse[i][count].index = -1;
    }

    return sparse;
}

void svmtrain (double *x, int *r, int *c,
	       double *y,
	       int    *svm_type,
	       int    *kernel_type,
	       double *degree,
	       double *gamma,
	       double *coef0,
	       double *cost,
	       double *nu,
	       double *cache,
	       double *tolerance,
	       double *epsilon,
	       int    *nr,
	       int    *index,
	       double *coefs,
	       double *rho)
{
    struct svm_parameter par;
    struct svm_problem   prob;
    struct svm_model    *model;
    int                  i, ii;
    
    /* 1. set parameter */
    par.svm_type    = *svm_type;
    par.kernel_type = *kernel_type;
    par.degree      = *degree;
    par.gamma       = *gamma;
    par.coef0       = *coef0;
    par.cache_size  = *cache;
    par.eps         = *tolerance;
    par.C           = *cost;
    par.nu          = *nu;
    par.p           = *epsilon;
    
    /* 2. set problem */
    prob.l = *r;
    prob.y = y;
    
    prob.x = sparsify (x, *r, *c);
    
    /* 3. call svm_train */
    model = svm_train (&prob, &par);
    
    /* 4. set up return values */
    for (i = ii = 0;
	(i < *r) && (ii < model->n);
	i++)
	if (prob.x[i] == model->SV[ii]) {
	    /* copy coef. */
	    coefs [i] = model->sv_coef[ii];

	    /* set index = true */
	    index [i] = 1;

	    ii++;
	}
    
    *rho = model->rho;
    *nr  = model->n;

    /* 5. clean up memory */
    svm_destroy_model (model);
    for (i = 0; i < *r; i++) free (prob.x[i]);
    free (prob.x);
}
	     
void svmclassify (double *v, int *r, int *c,
		  double *coefs,
		  double *rho,

		  int    *svm_type,
		  int    *kernel_type,
		  double *degree,
		  double *gamma,
		  double *coef0,

		  double *x, int *xr,
		  double *ret)
{
    struct svm_model m;
    struct svm_node ** train;
    double label = 0;
    int i;
    
    /* set up model */
    m.n       = *r;
    m.sv_coef = coefs;
    m.SV      = sparsify (v, *r, *c);
    m.rho     = *rho;

    /* set up parameter */
    m.param.svm_type    = *svm_type;
    m.param.kernel_type = *kernel_type;
    m.param.degree      = *degree;
    m.param.gamma       = *gamma;
    m.param.coef0       = *coef0;

    m.free_sv           = 1;

    /* create sparse training matrix */
    train = sparsify (x, *xr, *c);

    /* call svm-function for each x-row */
    for (i = 0; i < *xr; i++)
	svm_classify (&m, train[i], label, ret + i);

    /* clean up memory */
    for (i = 0; i < *xr; i++)
	free (train[i]);
    free (train);

    for (i = 0; i < *r; i++)
	free (m.SV[i]);
    free (m.SV);
    
}	     
		


