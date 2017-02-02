
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void
cmeans(double *x, int *nr_x, int *nc, double *p, int *nr_p, double *w,
       double *f, int *dist, int *itermax, double *reltol, int *verbose,
       double *u, double *ermin, int *iter);

int
cshell(int *xrows, int *xcols, double *x, int *ncenters,
       double *centers, int *itermax, int *iter, 
       int *verbose, int *dist, double *U, double *UANT, double
       *f, double *ermin, double *radius, int *flag);

int
e1071_floyd(int *n, double *A, double *C, int *P);

void
ufcl(double *x, int *nr_x, int *nc, double *p, int *nr_p, double *w,
     double *f, int *dist, int *itermax, double *reltol, int *verbose,
     double *rate_par,
     double *u, double *ermin, int *iter);

void
svmtrain (double *x, int *r, int *c, 
          double *y,
          int    *rowindex, int *colindex,
          int    *svm_type,
          int    *kernel_type,
          int    *degree,
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
          int    *probability,
	       
          int    *nclasses,
          int    *nr,
          int    *index,
          int    *labels,
          int    *nSV,
          double *rho,
          double *coefs,
          double *sigma,
          double *probA,
          double *probB,

          double *cresults,
          double *ctotal1,
          double *ctotal2,
          char   **error);

void
svmpredict  (int    *decisionvalues,
    	     int    *probability,

             double *v, int *r, int *c,
	     int    *rowindex,
	     int    *colindex,
	     double *coefs,
	     double *rho,
	     int    *compprob,
	     double *probA,
	     double *probB,
	     int    *nclasses,
	     int    *totnSV,
	     int    *labels,
	     int    *nSV,
	     int    *sparsemodel,

	     int    *svm_type,
	     int    *kernel_type,
	     int    *degree,
	     double *gamma,
	     double *coef0,

	     double *x, int *xr,
	     int    *xrowindex,
	     int    *xcolindex,
	     int    *sparsex,
		  
	     double *ret,
	     double *dec,
	     double *prob);

void
svmwrite (double *v, int *r, int *c,
	  int    *rowindex,
	  int    *colindex,
	  double *coefs,
	  double *rho,
	  int    *compprob,
          double *probA,
          double *probB,
	  int    *nclasses,
	  int    *totnSV,
	  int    *labels,
	  int    *nSV,
	  int    *sparsemodel,

	  int    *svm_type,
	  int    *kernel_type,
	  int    *degree,
	  double *gamma,
	  double *coef0,

	  char **filename);


static const R_CMethodDef CEntries[] = {
    {"cmeans", (DL_FUNC) &cmeans, 14},
    {"cshell", (DL_FUNC) &cshell, 15},
    {"e1071_floyd", (DL_FUNC) &e1071_floyd, 4},
    {"svmpredict", (DL_FUNC) &svmpredict, 30},
    {"svmtrain", (DL_FUNC) &svmtrain, 37},
    {"svmwrite", (DL_FUNC) &svmwrite, 21},
    {"ufcl", (DL_FUNC) &ufcl, 15},
    {NULL, NULL, 0}
};

void R_init_e1071(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
