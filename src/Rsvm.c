#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svm.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

/*
 * svm_model
 */
struct svm_model
{
    struct svm_parameter param; /* parameter */
    int nr_class;		/* number of classes, = 2 in
				   regression/one class svm */
    int l;			/* total #SV */
    struct svm_node **SV;	/* SVs (SV[l]) */
    double **sv_coef;	        /* coefficients for SVs in decision functions
			           (sv_coef[n-1][l]) */
    double *rho;		/* constants in decision functions
				   (rho[n*(n-1)/2]) */
    
    /* for classification only */

    int *label;		        /* label of each class (label[n]) */
    int *nSV;		        /* number of SVs for each class (nSV[n]) */
				/* nSV[0] + nSV[1] + ... + nSV[n-1] = l */
    /* XXX */
    int free_sv;		/* 1 if svm_model is created by
				   svm_load_model */
				/* 0 if svm_model is created by
				   svm_train */
};

/*
 * results from cross-validation
 */

struct crossresults
{
    double* results;
    double  total1;
    double  total2;
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

struct svm_node ** transsparse (double *x, int r, int *rowindex, int *colindex)
{
    struct svm_node** sparse;
    int i, ii, count = 0, nnz = 0;

    sparse = (struct svm_node **) malloc (r * sizeof(struct svm_node*));
    for (i = 0; i < r; i++) {
	/* allocate memory for column elements */
	nnz = rowindex[i+1] - rowindex[i];
	sparse[i] = (struct svm_node *) malloc ((nnz + 1) * sizeof(struct svm_node));

	/* set column elements */
	for (ii = 0; ii < nnz; ii++) {
	    sparse[i][ii].index = colindex[count] - 1;
	    sparse[i][ii].value = x[count];
	    count++;
	}

	/* set termination element */
	sparse[i][ii].index = -1;
    }    

    return sparse;
    
}    


/* Cross-Validation-routine from svm-train */
void do_cross_validation(struct svm_problem *prob,
			 struct svm_parameter *param,
			 int nr_fold,
			 double* cresults,
			 double* ctotal1,
			 double* ctotal2)
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;

	/* random shuffle */
	for(i=0; i<prob->l; i++)
	{
		int j = rand()%(prob->l-i);
		struct svm_node *tx;
		double ty;
			
		tx = prob->x[i];
		prob->x[i] = prob->x[j];
		prob->x[j] = tx;

		ty = prob->y[i];
		prob->y[i] = prob->y[j];
		prob->y[j] = ty;
	}

	for(i=0; i<nr_fold; i++)
	{
		int begin = i*prob->l/nr_fold;
		int end = (i+1)*prob->l/nr_fold;
		int j,k;
		struct svm_problem subprob;

		subprob.l = prob->l-(end-begin);
		subprob.x = Malloc(struct svm_node*,subprob.l);
		subprob.y = Malloc(double,subprob.l);
			
		k=0;
		for(j = 0; j < begin; j++)
		{
			subprob.x[k] = prob->x[j];
			subprob.y[k] = prob->y[j];
			++k;
		}
		for(j = end; j<prob->l; j++)
		{
			subprob.x[k] = prob->x[j];
			subprob.y[k] = prob->y[j];
			++k;
		}

		if(param->svm_type == EPSILON_SVR ||
		   param->svm_type == NU_SVR)
		{
			struct svm_model *submodel = svm_train(&subprob,param);
			double error = 0;
			for(j=begin;j<end;j++)
			{
				double v = svm_predict(submodel,prob->x[j]);
				double y = prob->y[j];
				error += (v-y)*(v-y);
				sumv += v;
				sumy += y;
				sumvv += v*v;
				sumyy += y*y;
				sumvy += v*y;
			}
			svm_destroy_model(submodel);
			/* printf("Mean squared error = %g\n",
			   error/(end-begin)); */
			cresults[i] = error/(end-begin);
			total_error += error;			
		}
		else
		{
			struct svm_model *submodel = svm_train(&subprob,param);
			int correct = 0;
			for(j=begin;j<end;j++)
			{
				double v = svm_predict(submodel,prob->x[j]);
				if(v == prob->y[j])
					++correct;
			}
			svm_destroy_model(submodel);
			/* printf("Accuracy = %g%% (%d/%d)\n", */
			/* 100.0*correct/(end-begin),correct,(end-begin)); */
			cresults[i] = 100.0*correct/(end-begin);
			total_correct += correct;
		}

		free(subprob.x);
		free(subprob.y);
	}
	
	if(param->svm_type == EPSILON_SVR || param->svm_type == NU_SVR)
	{
	    /* printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
	        printf("Cross Validation Squared correlation coefficient = %g\n",
	    	((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
	    	((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
	    	); */
	    *ctotal1 = total_error/prob->l;
	    *ctotal2 = ((prob->l * sumvy - sumv * sumy) *
			(prob->l * sumvy - sumv*sumy))  /
		       ((prob->l * sumvv - sumv * sumv) *
		        (prob->l * sumyy - sumy * sumy));
	}
	else
	    /* printf("Cross Validation Accuracy =
	       %g%%\n",100.0*total_correct/prob.l); */
	    *ctotal1 = 100.0 * total_correct / prob->l;
}


void svmtrain (double *x, int *r, int *c, 
	       double *y,
	       int    *rowindex, int *colindex,
	       int    *svm_type,
	       int    *kernel_type,
	       double *degree,
	       double *gamma,
	       double *coef0,
	       double *cost,
	       double *nu,
	       int    *weightlabels,
	       double *weights,
	       int    *nweights,
	       double *cache,
	       double *tolerance,
	       double *epsilon,
	       int    *shrinking,
	       int    *cross,
	       int    *sparse,
	       
	       int    *nclasses,
	       int    *nr,
	       int    *index,
	       int    *labels,
	       int    *nSV,
	       double *rho,
	       double *coefs,

	       double *cresults,
	       double *ctotal1,
	       double *ctotal2,
	       char   **error)
{
    struct svm_parameter par;
    struct svm_problem   prob;
    struct svm_model    *model = NULL;
    int i, ii;
    const char* s;
    
    /* set parameters */
    par.svm_type    = *svm_type;
    par.kernel_type = *kernel_type;
    par.degree      = *degree;
    par.gamma       = *gamma;
    par.coef0       = *coef0;
    par.cache_size  = *cache;
    par.eps         = *tolerance;
    par.C           = *cost;
    par.nu          = *nu;
    par.nr_weight   = *nweights;
    if (par.nr_weight > 0) {
	par.weight      = (double *) malloc (sizeof(double) * par.nr_weight);
	memcpy (par.weight, weights, par.nr_weight * sizeof(double));
	par.weight_label = (int *) malloc (sizeof(int) * par.nr_weight);
	memcpy (par.weight_label, weightlabels, par.nr_weight * sizeof(int));
    }
    par.p           = *epsilon;
    par.shrinking   = *shrinking;

    /* set problem */
    prob.l = *r;
    prob.y = y;
    
    if (*sparse > 0)
	prob.x = transsparse(x, *r, rowindex, colindex);
    else
	prob.x = sparsify(x, *r, *c);
    
    /* check parameters & copy error message */
    s = svm_check_parameter(&prob, &par);
    if (s) {
	strcpy(*error, s);
    } else {
	/* call svm_train */
	model = svm_train(&prob, &par);
    
	/* set up return values */
	for (ii = 0; ii < model->l; ii++)
	    for (i = 0; i < *r;	i++)
		if (prob.x[i] == model->SV[ii]) index[ii] = i+1;
	
	*nr  = model->l;
	*nclasses = model->nr_class;
	memcpy (rho, model->rho, *nclasses * (*nclasses - 1)/2 * sizeof(double));
	for (i = 0; i < *nclasses-1; i++)
	    memcpy (coefs + i * *nr, model->sv_coef[i],  *nr * sizeof (double));
	
	if (*svm_type < 2) {
	    memcpy (labels, model->label, *nclasses * sizeof(int));
	    memcpy (nSV, model->nSV, *nclasses * sizeof(int));
	}
	
	/* Perform cross-validation, if requested */
	if (*cross > 0)
	    do_cross_validation (&prob, &par, *cross, cresults,
				 ctotal1, ctotal2);
	
	/* clean up memory */
	svm_destroy_model (model);
    }
    
    /* clean up memory */
    if (par.nr_weight > 0) {
	free(par.weight);
	free(par.weight_label);
    }
    
    for (i = 0; i < *r; i++) free (prob.x[i]);
    free (prob.x);
}
	     
void svmpredict  (double *v, int *r, int *c,
		  int    *rowindex,
		  int    *colindex,
		  double *coefs,
		  double *rho,
		  int    *nclasses,
		  int    *totnSV,
		  int    *labels,
		  int    *nSV,
		  int    *sparsemodel,

		  int    *svm_type,
		  int    *kernel_type,
		  double *degree,
		  double *gamma,
		  double *coef0,

		  double *x, int *xr,
		  int    *xrowindex,
		  int    *xcolindex,
		  int    *sparsex,
		  
		  double *ret)
{
    struct svm_model m;
    struct svm_node ** train;
    int i;
    
    /* set up model */
    m.l        = *totnSV;
    m.nr_class = *nclasses;
    m.sv_coef  = (double **) malloc (m.nr_class * sizeof(double));
    for (i = 0; i < m.nr_class - 1; i++) {
      m.sv_coef[i] = (double *) malloc (m.l * sizeof (double));
      memcpy (m.sv_coef[i], coefs + i*m.l, m.l * sizeof (double));
    }
    
    if (*sparsemodel > 0)
	m.SV   = transsparse(v, *r, rowindex, colindex);
    else
	m.SV   = sparsify(v, *r, *c);
    
    m.rho      = rho;
    m.label    = labels;
    m.nSV      = nSV;

    /* set up parameter */
    m.param.svm_type    = *svm_type;
    m.param.kernel_type = *kernel_type;
    m.param.degree      = *degree;
    m.param.gamma       = *gamma;
    m.param.coef0       = *coef0;

    m.free_sv           = 1;

    /* create sparse training matrix */
    if (*sparsex > 0)
	train = transsparse(x, *xr, xrowindex, xcolindex);
    else
	train = sparsify(x, *xr, *c);

    /* call svm-function for each x-row */
    for (i = 0; i < *xr; i++)
	ret[i] = svm_predict (&m, train[i]);

    /* clean up memory */
    for (i = 0; i < *xr; i++)
	free (train[i]);
    free (train);

    for (i = 0; i < *r; i++)
	free (m.SV[i]);
    free (m.SV);
    
    for (i = 0; i < m.nr_class - 1; i++)
      free(m.sv_coef[i]);
    free(m.sv_coef);
}	     
		



