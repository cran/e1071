#ifdef __cplusplus
extern "C" {
#endif

struct svm_node
{
	int index;
	double value;
};

struct svm_problem
{
	int l;
	double *y;
	struct svm_node **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, C_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID };	/* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	double degree;	// for poly
	double gamma;	// for poly/rbf/sigmoid
	double coef0;	// for poly/sigmoid

	// these are for training only
	double cache_size; // in MB
	double eps;	// stopping criteria
	double C;	// for C_SVC and C_SVR
	double nu;	// for NU_SVC and ONE_CLASS
	double p;	// for C_SVR
};

struct svm_model;

struct svm_model *svm_train(const struct svm_problem *prob,
			    const struct svm_parameter *param);

int svm_save_model(const char *model_file_name, const struct svm_model *model);

struct svm_model *svm_load_model(const char *model_file_name);

int svm_classify(const struct svm_model *model, const struct svm_node *x,
		 double label, double *decision_value);

void svm_destroy_model(struct svm_model *model);

#ifdef __cplusplus
}
#endif
