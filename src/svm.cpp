#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "svm.h"
template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
#define EPS_A 1e-12
#define INF DBL_MAX

//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size
//
class Cache
{
public:
	Cache(int l,int size);
	~Cache();

	// return 1 if returned data is valid (cache hit)
	// return 0 if returned data needs to be filled (cache miss)
	int get_data(const int index, double **data);
private:
	int l;
	int size;
	struct head_t
	{
		head_t *prev, *next;	// a cicular list
		int index;
		double *data;
	};

	head_t* head;
	head_t* lru_head;
	head_t** index_to_head;
	void move_to_last(head_t *h);
};

Cache::Cache(int _l,int _size):l(_l),size(_size)
{
	head = new head_t[size];
	int i;
	for(i=0;i<size;i++)
	{
		head[i].next = &head[i+1];
		head[i].prev = &head[i-1];
		head[i].index = -1;
		head[i].data = new double[l];
	}

	head[0].prev = &head[size-1];
	head[size-1].next = &head[0];
	lru_head = &head[0];

	index_to_head = new head_t *[l];
	for(i=0;i<l;i++)
		index_to_head[i] = 0;
}

Cache::~Cache()
{
	for(int i=0;i<size;i++)
		delete[] head[i].data;
	delete[] head;
	delete[] index_to_head;
}

void Cache::move_to_last(head_t *h)
{
	if(lru_head == h)
		lru_head = lru_head->next;
	else
	{
		// delete from current location
		h->prev->next = h->next;
		h->next->prev = h->prev;

		// insert to last position
		h->next = lru_head;
		h->prev = lru_head->prev;
		h->prev->next = h;
		h->next->prev = h;
	}
}

int Cache::get_data(const int index,double **data)
{
	head_t *h=index_to_head[index];
	if(h)
	{
		move_to_last(h);
		*data = h->data;
		return 1;
	}		
	else
	{
		// get one from lru_head
		h = lru_head;
		lru_head = lru_head->next;
		if(h->index!=-1)
			index_to_head[h->index] = 0;
		h->index = index;
		index_to_head[index] = h;
		*data = h->data;
		return 0;
	}
}

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//
class Kernel {
public:
	Kernel(int l, const svm_node * const * x, const svm_parameter& param);
	virtual ~Kernel();

	static double k_function(const svm_node *x, const svm_node *y,
				 const svm_parameter& param);
	virtual double *get_Q(int column) const = 0;

protected:

	double (Kernel::*kernel_function)(int i, int j) const;

private:
	const svm_node * const * const x;

	// svm_parameter
	const int kernel_type;
	const double degree;
	const double gamma;
	const double coef0;

	double *x_square;
	static double dot(const svm_node *px, const svm_node *py);
	double kernel_linear(int i, int j) const
	{
		return dot(x[i],x[j]);
	}
	double kernel_poly(int i, int j) const
	{
		return pow(gamma*dot(x[i],x[j])+coef0,degree);
	}
	double kernel_rbf(int i, int j) const
	{
		return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
	}
	double kernel_sigmoid(int i, int j) const
	{
		return tanh(gamma*dot(x[i],x[j])+coef0);
	}
};

Kernel::Kernel(int l, const svm_node * const * _x, const svm_parameter& param)
:x(_x), kernel_type(param.kernel_type), degree(param.degree),
 gamma(param.gamma), coef0(param.coef0)
{

	switch(kernel_type)
	{
		case LINEAR:
			kernel_function = &Kernel::kernel_linear;
			break;
		case POLY:
			kernel_function = &Kernel::kernel_poly;
			break;
		case RBF:
			kernel_function = &Kernel::kernel_rbf;
			break;
		case SIGMOID:
			kernel_function = &Kernel::kernel_sigmoid;
			break;
		default:
			fprintf(stderr,"unknown kernel function.\n");
			exit(1);
	}

	if(kernel_type == RBF)
	{
		x_square = new double[l];
		for(int i=0;i<l;i++)
			x_square[i] = dot(x[i],x[i]);
	}
	else
		x_square = 0;
}

Kernel::~Kernel()
{
	delete x_square;
}

double Kernel::dot(const svm_node *px, const svm_node *py)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			sum += px->value * py->value;
			++px;
			++py;
		}
		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}			
	}
	return sum;
}

double Kernel::k_function(const svm_node *x, const svm_node *y,
			  const svm_parameter& param)
{
	switch(param.kernel_type)
	{
		case LINEAR:
			return dot(x,y);
		case POLY:
			return pow(param.gamma*dot(x,y)+param.coef0,param.degree);
		case RBF:
		{
			double sum = 0;
			while(x->index != -1 && y->index !=-1)
			{
				if(x->index == y->index)
				{
					double d = x->value - y->value;
					sum += d*d;
					++x;
					++y;
				}
				else
				{
					if(x->index > y->index)
					{	
						sum += y->value * y->value;
						++y;
					}
					else
					{
						sum += x->value * x->value;
						++x;
					}
				}
			}

			while(x->index != -1)
			{
				sum += x->value * x->value;
				++x;
			}

			while(y->index != -1)
			{
				sum += y->value * y->value;
				++y;
			}
			
			return exp(-param.gamma*sum);
		}
		case SIGMOID:
			return tanh(param.gamma*dot(x,y)+param.coef0);
		default:
			fprintf(stderr,"unknown kernel function.\n");
			exit(1);
	}
}

// Generalized SMO+SVMlight algorithm
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + b^T \alpha
//
//		0 <= alpha_i <= C
//		y^T \alpha = \delta
//
// Given:
//
//	Q, b, y, C, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping criterion
//
// solution will be put in \alpha, objective value will be put in obj
//
class Solver {
public:
	Solver( int _l, const Kernel& Q, const double *b, const double *_y,
		double *_alpha, double _C, double _eps,
		double& obj, double& rho);

	~Solver();

private:
	const int l;
	const double *y;
	double * const alpha;
	const double C;
	const double eps;

	double * G;		// gradient of objective function
	enum { LOWER_BOUND, UPPER_BOUND, FREE };
	char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	void update_alpha_status(int i)
	{
		if(alpha[i] >= C-EPS_A)
			alpha_status[i] = UPPER_BOUND;
		else if(alpha[i] <= EPS_A)
			alpha_status[i] = LOWER_BOUND;
		else alpha_status[i] = FREE;
	}
	bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
	bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
	int select_working_set(int &i, int &j);
};

Solver::Solver( int _l, const Kernel& Q, const double *b, const double *_y,
		double *_alpha, double _C, double _eps,
		double& obj, double& rho)
:l(_l),y(_y),alpha(_alpha),C(_C),eps(_eps)
{
	// initialize alpha_status
	{
		alpha_status = new char[l];
		for(int i=0;i<l;i++)
			update_alpha_status(i);
	}

	// initialize gradient
	{
		G = new double[l];
		int i;
		for(i=0;i<l;i++)
			G[i] = b[i];
		for(i=0;i<l;i++)
			if(!is_lower_bound(i))
			{
				double *Q_i = Q.get_Q(i);
				double alpha_i = alpha[i];
				for(int j=0;j<l;j++)
					G[j] += alpha_i*Q_i[j];
			}
	}

	// optimization step

	int iter = 0;
	int counter = 0;

	while(1)
	{
		int i,j;
		if(select_working_set(i,j)!=0)
			break;

		++iter;

		// show progress

		if(++counter == 100)
		{
			fprintf(stderr,".");
			counter = 0;
		}

		// update alpha[i] and alpha[j]
		
		const double *Q_i = Q.get_Q(i);
		const double *Q_j = Q.get_Q(j);

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];
		double s = y[i]*old_alpha_i + y[j]*old_alpha_j;

		double H = s/y[j];
		double L = (s-C*y[i])/y[j];

		if(y[i]*y[j] < 0) { double t = H; H = L; L = t;}

		H = min(C,H);
		L = max(0.0,L);

		alpha[j] += y[i] * (y[j]*G[i] - y[i]*G[j]) /
			(y[i]*(y[i]*Q_j[j] - 2*Q_i[j]*y[j]) + y[j]*y[j]*Q_i[i]);

		if(alpha[j] > H) alpha[j] = H;
		else if(alpha[j] < L) alpha[j] = L;

		alpha[i] = (s - y[j]*alpha[j])/y[i];

		// update alpha_status

		update_alpha_status(i);
		update_alpha_status(j);

		// update G

		double delta_alpha_i = alpha[i] - old_alpha_i;
		double delta_alpha_j = alpha[j] - old_alpha_j;
		
		for(int k=0;k<l;k++)
		{
			G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
		}
	}
	
	// calculate rho

	{
		double r;
		int nr_free = 0;
		double ub = INF, lb = -INF, sum_free = 0;
		for(int i=0;i<l;i++)
		{
			double yG = y[i]*G[i];

			if(is_lower_bound(i))
			{
				if(y[i] > 0)
					ub = min(ub,yG);
				else
					lb = max(lb,yG);
			}
			else if(is_upper_bound(i))
			{
				if(y[i] < 0)
					ub = min(ub,yG);
				else
					lb = max(lb,yG);
			}
			else
			{
				++nr_free;
				sum_free += yG;
			}
		}

		if(nr_free>0)
			r = sum_free/nr_free;
		else
			r = (ub+lb)/2;

		rho = r;
	}

	// calculate objective value

	{
		double v = 0;
		int i;
		for(i=0;i<l;i++)
			v += alpha[i] * (G[i]/2 + b[i]);

		obj = v;
	}

	printf("\noptimization finished, #iter = %d\n",iter);
}

Solver::~Solver()
{
	delete[] alpha_status;
	delete[] G;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set(int &out_i, int &out_j)
{
	// return i,j which maximize -grad(f)^T d , under constraint
	// if alpha_i == C, d != +1
	// if alpha_i == 0, d != -1

	double Gmax1 = -INF;		// max { -grad(f)_i * d | y_i*d = +1 }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;		// max { -grad(f)_i * d | y_i*d = -1 }
	int Gmax2_idx = -1;

	for(int i=0;i<l;i++)
	{
		if(y[i]>0)	// y > 0
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax1)
				{
					Gmax1 = -G[i];
					Gmax1_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax2)
				{
					Gmax2 = G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y < 0
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax2)
				{
					Gmax2 = -G[i];
					Gmax2_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax1)
				{
					Gmax1 = G[i];
					Gmax1_idx = i;
				}
			}
		}
	}

	if(Gmax1+Gmax2 < eps)
 		return 1;

	out_i = Gmax1_idx;
	out_j = Gmax2_idx;
	return 0;
}

//
// Q matrices for different formulations
//
class C_SVC_Q: public Kernel
{ 
public:
	C_SVC_Q(const svm_problem& prob, const svm_parameter& param, const double *_y)
	:Kernel(prob.l, prob.x, param), l(prob.l), y(_y)
	{
		cache = new Cache(l,(int)min((double)l,(param.cache_size*(1<<20))/(sizeof(double)*l)));
	}
	
	double *get_Q(int i) const
	{
		double *data;

		if(cache->get_data(i,&data) == 0)
		{
			for(int j=0;j<l;j++)
				data[j] = y[i]*y[j]*(this->*kernel_function)(i,j);
		}
		return data;
	}

	~C_SVC_Q()
	{
		delete cache;
	}
private:
	const int l;
	const double *y;
	Cache *cache;
};

class NU_SVC_Q: public Kernel
{
public:
	NU_SVC_Q(const svm_problem& prob, const svm_parameter& param, const double *_y)
	:Kernel(prob.l, prob.x, param), l(prob.l), y(_y)
	{
		cache = new Cache(l,(int)min((double)l,(param.cache_size*(1<<20))/(sizeof(double)*l)));
	}
	
	double *get_Q(int i) const
	{
		double *data;

		if(cache->get_data(i,&data) == 0)
		{
			for(int j=0;j<l;j++)
				data[j] = y[i]*y[j]*(1+(this->*kernel_function)(i,j));
		}
		return data;
	}

	~NU_SVC_Q()
	{
		delete cache;
	}
private:
	const int l;
	const double *y;
	Cache *cache;
};

class ONE_CLASS_Q: public Kernel
{
public:
	ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
	:Kernel(prob.l, prob.x, param), l(prob.l)
	{
		cache = new Cache(l,(int)min((double)l,(param.cache_size*(1<<20))/(sizeof(double)*l)));
	}
	
	double *get_Q(int i) const
	{
		double *data;

		if(cache->get_data(i,&data) == 0)
		{
			for(int j=0;j<l;j++)
				data[j] = (this->*kernel_function)(i,j);
		}
		return data;
	}

	~ONE_CLASS_Q()
	{
		delete cache;
	}
private:
	const int l;
	Cache *cache;
};

class C_SVR_Q: public Kernel
{ 
public:
	C_SVR_Q(const svm_problem& prob, const svm_parameter& param)
	:Kernel(prob.l, prob.x, param), l(prob.l)
	{
		cache = new Cache(2*l,(int)min((double)(2*l),(param.cache_size*(1<<20))/(sizeof(double)*(2*l))));
	}
	
	double *get_Q(int i) const
	{
		double *data;

		if(cache->get_data(i,&data) == 0)
		{
			if(i<l)
				for(int j=0;j<l;j++)
				{
					data[j] = (this->*kernel_function)(i,j);
					data[j+l] = -data[j];
				}
			else
				for(int j=0;j<l;j++)
				{
					data[j+l] = (this->*kernel_function)(i-l,j);
					data[j] = -data[j+l];
				}
		}
		return data;
	}

	~C_SVR_Q()
	{
		delete cache;
	}
private:
	const int l;
	const double *y;
	Cache *cache;
};

//
// svm_model
//
struct svm_model
{
	int n;			// number of SVs
	double *sv_coef;	// sv_coef[i] is the coefficient of SV[i]
	svm_node const ** SV;	// SVs
	double rho;		// the constant in the decision function

	svm_parameter param;	// parameter

	int free_sv;		// XXX: 1 if svm_model is created by svm_load_model
				//      0 if svm_model is created by svm_train
};

static void solve_c_svc(
	const svm_problem *prob, const svm_parameter* param,
	double *alpha, double& obj, double& rho)
{
	int l = prob->l;
	double *minus_ones = new double[l];
	int i;

	for(i=0;i<l;i++)
	{
		alpha[i] = 0;
		minus_ones[i] = -1;
	}

	Solver s(l, C_SVC_Q(*prob,*param,prob->y), minus_ones, prob->y,
		 alpha, param->C, param->eps, obj, rho);

	double sum_alpha=0;
	for(i=0;i<l;i++)
		sum_alpha += alpha[i];

	printf("nu = %f\n", sum_alpha/(param->C*prob->l));

	delete[] minus_ones;

	for(i=0;i<l;i++)
		alpha[i] *= prob->y[i];
}

static void solve_nu_svc(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, double& obj, double& rho)
{
	int l = prob->l;
	double *zeros = new double[l];
	double *ones = new double[l];
	int i;

	int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound
	if(n>=prob->l)
	{
		fprintf(stderr,"nu must be in (0,1)\n");
		exit(1);
	}
	for(i=0;i<n;i++)
		alpha[i] = 1;
	alpha[n] = param->nu * prob->l - n;
	for(i=n+1;i<l;i++)
		alpha[i] = 0;

	for(i=0;i<l;i++)
	{	
		zeros[i] = 0;
		ones[i] = 1;
	}

	Solver s(l, NU_SVC_Q(*prob,*param,prob->y), zeros, ones,
		 alpha, 1.0, param->eps, obj, rho);

	printf("C = %f\n",1/rho);

	delete[] zeros;
	delete[] ones;

	for(i=0;i<l;i++)
		alpha[i] *= prob->y[i]/rho;
}

static void solve_one_class(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, double& obj, double& rho)
{
	int l = prob->l;
	double *zeros = new double[l];
	double *ones = new double[l];
	int i;

	int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound
	if(n>=prob->l)
	{
		fprintf(stderr,"nu must be in (0,1)\n");
		exit(1);
	}
	for(i=0;i<n;i++)
		alpha[i] = 1;
	alpha[n] = param->nu * prob->l - n;
	for(i=n+1;i<l;i++)
		alpha[i] = 0;

	for(i=0;i<l;i++)
	{
		zeros[i] = 0;
		ones[i] = 1;
	}

	Solver s(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
		 alpha, 1.0, param->eps, obj, rho);

	delete[] zeros;
	delete[] ones;
}

static void solve_c_svr(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, double& obj, double& rho)
{
	int l = prob->l;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	double *y = new double[2*l];
	int i;

	for(i=0;i<l;i++)
	{
		alpha2[i] = 0;
		linear_term[i] = param->p - prob->y[i];
		y[i] = 1;

		alpha2[i+l] = 0;
		linear_term[i+l] = param->p + prob->y[i];
		y[i+l] = -1;
	}

	Solver s(2*l, C_SVR_Q(*prob,*param), linear_term, y,
		 alpha2, param->C, param->eps, obj, rho);

	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];

	delete[] alpha2;
	delete[] linear_term;
	delete[] y;
}

//
// Interface functions
//
svm_model *svm_train(const svm_problem *prob, const svm_parameter *param)
{
	svm_model *model = (svm_model *)malloc(sizeof(svm_model));
	model->param = *param;
	double *alpha = new double[prob->l];
	double obj, rho;
	switch(param->svm_type)
	{
		case C_SVC:
			solve_c_svc(prob,param,alpha,obj,rho);
			break;
		case NU_SVC:
			solve_nu_svc(prob,param,alpha,obj,rho);
			break;
		case ONE_CLASS:
			solve_one_class(prob,param,alpha,obj,rho);
			break;
		case C_SVR:
			solve_c_svr(prob,param,alpha,obj,rho);
			break;
	}

	model->rho = rho;
	printf("obj = %f, rho = %f\n",obj,rho);

	// output SVs

	int nSV = 0;
	int nBSV = 0;
	for(int i=0;i<prob->l;i++)
	{
		if(fabs(alpha[i]) >= EPS_A)
		{
			++nSV;
			if(fabs(alpha[i]) >= param->C - EPS_A)
				++nBSV;
		}
	}

	printf("nSV = %d, nBSV = %d\n",nSV,nBSV);

	model->n = nSV;
	model->sv_coef = (double *)malloc(sizeof(double) * nSV);
	model->SV = (const svm_node **)malloc(sizeof(svm_node*) * nSV);

	{
		int j = 0;
		for(int i=0;i<prob->l;i++)
		{
			if(fabs(alpha[i]) >= EPS_A)
			{
				model->sv_coef[j] = alpha[i];
				model->SV[j] = prob->x[i];
				++j;
			}
		}
	}

	delete[] alpha;
	model->free_sv = 0;	// XXX
	return model;
}

int svm_classify(const svm_model *model, const svm_node *x,
		 double label, double *decision_value)
{
	const int n = model->n;
	const double *sv_coef = model->sv_coef;

	double sum = 0;
	if(model->param.svm_type == NU_SVC)
	{
		for(int i=0;i<n;i++)
			sum += sv_coef[i] * (1+Kernel::k_function(x,model->SV[i],model->param));
	}
	else
	{
		for(int i=0;i<n;i++)
			sum += sv_coef[i] * Kernel::k_function(x,model->SV[i],model->param);
		sum-=(model->rho);
	}

	*decision_value = sum;

	if(model->param.svm_type == ONE_CLASS)
		return (sum > 0);
	else
		return (sum*label)>0;
}

const char *svm_type_table[] =
{
	"c_svc","nu_svc","one_class","c_svr",NULL
};

const char *kernel_type_table[]=
{
	"linear","polynomial","rbf","sigmoid",NULL
};

int svm_save_model(const char *model_file_name, const svm_model *model)
{
	FILE *fp = fopen(model_file_name,"w");
	if(fp==NULL) return -1;

	const svm_parameter& param = model->param;

	fprintf(fp,"svm_type %s\n", svm_type_table[param.svm_type]);
	fprintf(fp,"kernel_type %s\n", kernel_type_table[param.kernel_type]);

	if(param.kernel_type == POLY)
		fprintf(fp,"degree %g\n", param.degree);

	if(param.kernel_type == POLY || param.kernel_type == RBF || param.kernel_type == SIGMOID)
		fprintf(fp,"gamma %g\n", param.gamma);

	if(param.kernel_type == POLY || param.kernel_type == SIGMOID)
		fprintf(fp,"coef0 %g\n", param.coef0);

	if(param.svm_type == C_SVC || param.svm_type == ONE_CLASS ||
	   param.svm_type == C_SVR)
		fprintf(fp, "rho %g\n", model->rho);

	fprintf(fp, "SV %d\n", model->n);

	const int n = model->n;
	const double * const sv_coef = model->sv_coef;
	const svm_node **SV = model->SV;

	for(int i=0;i<n;i++)
	{
		fprintf(fp, "%.16g ",sv_coef[i]);
		const svm_node *p = SV[i];
		while(p->index != -1)
		{
			fprintf(fp,"%d:%.8g ",p->index,p->value);
			p++;
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	return 0;
}

svm_model *svm_load_model(const char *model_file_name)
{
	FILE *fp = fopen(model_file_name,"rb");
	if(fp==NULL) return NULL;
	
	// read n,b,param

	svm_model *model = (svm_model *)malloc(sizeof(svm_model));
	svm_parameter& param = model->param;

	char cmd[81];
	while(1)
	{
		fscanf(fp,"%80s",cmd);

		if(strcmp(cmd,"svm_type")==0)
		{
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;svm_type_table[i];i++)
			{
				if(strcmp(svm_type_table[i],cmd)==0)
				{
					param.svm_type=i;
					break;
				}
			}
			if(svm_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown svm type.\n");
				exit(1);
			}
		}
		else if(strcmp(cmd,"kernel_type")==0)
		{		
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;kernel_type_table[i];i++)
			{
				if(strcmp(kernel_type_table[i],cmd)==0)
				{
					param.kernel_type=i;
					break;
				}
			}
			if(kernel_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown kernel function.\n");
				exit(1);
			}
		}
		else if(strcmp(cmd,"degree")==0)
			fscanf(fp,"%lf",&param.degree);
		else if(strcmp(cmd,"gamma")==0)
			fscanf(fp,"%lf",&param.gamma);
		else if(strcmp(cmd,"coef0")==0)
			fscanf(fp,"%lf",&param.coef0);
		else if(strcmp(cmd,"rho")==0)
			fscanf(fp,"%lf",&model->rho);
		else if(strcmp(cmd,"SV")==0)
		{
			fscanf(fp,"%d",&model->n);
			while(1)
			{
				int c = getc(fp);
				if(c==EOF || c=='\n') break;	
			}
			break;
		}
		else
		{
			fprintf(stderr,"unknown text in model file\n");
			exit(1);
		}
	}

	// read sv_coef and SV

	int elements = 0;
	long pos = ftell(fp);

	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
			case '\n':
				// count the '-1' element
			case ':':
				++elements;
				break;
			case EOF:
				goto out;
			default:
				;
		}
	}
out:
	fseek(fp,pos,SEEK_SET);

	const int n = model->n;
	model->sv_coef = (double*)malloc(sizeof(double)*n);
	model->SV = (const svm_node**)malloc(sizeof(svm_node*)*n);
	svm_node *x_space = (svm_node*)malloc(sizeof(svm_node)*elements);

	int j=0;
	for(int i=0;i<n;i++)
	{
		model->SV[i] = &x_space[j];
		fscanf(fp,"%lf",&model->sv_coef[i]);
		while(1)
		{
			int c;
			do {
				c = getc(fp);
				if(c=='\n') goto out2;
			} while(isspace(c));
			ungetc(c,fp);
			fscanf(fp,"%d:%lf",&(x_space[j].index),&(x_space[j].value));
			++j;
		}	
out2:
		x_space[j++].index = -1;
	}

	fclose(fp);

	model->free_sv = 1;	// XXX
	return model;
}

void svm_destroy_model(svm_model* model)
{
	if(model->free_sv)
		free((void *)(model->SV[0]));
	free(model->sv_coef);
	free(model->SV);
	free(model);
}
