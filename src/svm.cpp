#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "svm.h"
typedef signed char schar;
template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
template <class T> inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
	dst = new T[n];
	memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
#define EPS_A 1e-12
#define INF HUGE_VAL

//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size limit in bytes
//
class Cache
{
public:
	Cache(int l,int size);
	~Cache();

	// request data [0,len)
	// return some position p where [p,len) need to be filled
	// (p >= len if nothing needs to be filled)
	int get_data(const int index, double **data, int len);
	void swap_index(int i, int j);	// future_option
private:
	int l;
	int size;
	struct head_t
	{
		head_t *prev, *next;	// a cicular list
		double *data;
		int len;		// data[0,len) is cached in this entry
	};

	head_t* head;
	head_t lru_head;
	void lru_delete(head_t *h);
	void lru_insert(head_t *h);
};

Cache::Cache(int l_,int size_):l(l_),size(size_)
{
	head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
	size /= sizeof(double);
	size -= l * sizeof(head_t) / sizeof(double);
	lru_head.next = lru_head.prev = &lru_head;
}

Cache::~Cache()
{
	for(head_t *h = lru_head.next; h != &lru_head; h=h->next)
		free(h->data);
	free(head);
}

void Cache::lru_delete(head_t *h)
{
	// delete from current location
	h->prev->next = h->next;
	h->next->prev = h->prev;
}

void Cache::lru_insert(head_t *h)
{
	// insert to last position
	h->next = &lru_head;
	h->prev = lru_head.prev;
	h->prev->next = h;
	h->next->prev = h;
}

int Cache::get_data(const int index, double **data, int len)
{
	head_t *h = &head[index];
	if(h->len) lru_delete(h);
	int more = len - h->len;

	if(more > 0)
	{
		// free old space
		while(size < more)
		{
			head_t *old = lru_head.next;
			lru_delete(old);
			free(old->data);
			size += old->len;
			old->data = 0;
			old->len = 0;
		}

		// allocate new space
		h->data = (double *)realloc(h->data,sizeof(double)*len);
		size -= more;
		swap(h->len,len);
	}

	lru_insert(h);
	*data = h->data;
	return len;
}

void Cache::swap_index(int i, int j)
{
	if(i==j) return;

	if(head[i].len) lru_delete(&head[i]);
	if(head[j].len) lru_delete(&head[j]);
	swap(head[i].data,head[j].data);
	swap(head[i].len,head[j].len);
	if(head[i].len) lru_insert(&head[i]);
	if(head[j].len) lru_insert(&head[j]);

	if(i>j) swap(i,j);
	for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
	{
		if(h->len > i)
		{
			if(h->len > j)
				swap(h->data[i],h->data[j]);
			else
			{
				// give up
				lru_delete(h);
				free(h->data);
				size += h->len;
				h->data = 0;
				h->len = 0;
			}
		}
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
	virtual double *get_Q(int column, int len) const = 0;
	virtual void swap_index(int i, int j) const	// no so const...
	{
		swap(x[i],x[j]);
		if(x_square) swap(x_square[i],x_square[j]);
	}
protected:

	double (Kernel::*kernel_function)(int i, int j) const;

private:
	const svm_node **x;
	double *x_square;

	// svm_parameter
	const int kernel_type;
	const double degree;
	const double gamma;
	const double coef0;

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

Kernel::Kernel(int l, const svm_node * const * x_, const svm_parameter& param)
:kernel_type(param.kernel_type), degree(param.degree),
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

	clone(x,x_,l);

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
	delete[] x;
	delete[] x_square;
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
//		y_i = +1 or -1
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
	Solver() {};
	virtual ~Solver() {};
	void Solve(int l, const Kernel& Q, const double *b_, const schar *y_,
		   double *alpha_, double C, double eps,
		   double& obj, double& rho, int shrinking);
protected:
	int active_size;
	schar *y;
	double *G;		// gradient of objective function
	enum { LOWER_BOUND, UPPER_BOUND, FREE };
	char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	double *alpha;
	const Kernel *Q;
	double eps;
	double C;
	double *b;
	int *active_set;
	double *G_bar;		// gradient, if we treat free variables as 0
	int l;
	bool unshrinked;		// XXX

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
	bool is_free(int i) { return alpha_status[i] == FREE; }
	void swap_index(int i, int j);
	void reconstruct_gradient();
	virtual int select_working_set(int &i, int &j);
	virtual double calculate_rho();
	virtual void do_shrinking();
};

void Solver::swap_index(int i, int j)
{
	Q->swap_index(i,j);
	swap(y[i],y[j]);
	swap(G[i],G[j]);
	swap(alpha_status[i],alpha_status[j]);
	swap(alpha[i],alpha[j]);
	swap(b[i],b[j]);
	swap(active_set[i],active_set[j]);
	swap(G_bar[i],G_bar[j]);
}

void Solver::reconstruct_gradient()
{
	// reconstruct inactive elements of G from G_bar and free variables

	if(active_size == l) return;

	int i;
	for(i=active_size;i<l;i++)
		G[i] = C * G_bar[i] + b[i];
	
	for(i=0;i<active_size;i++)
		if(is_free(i))
		{
			const double *Q_i = Q->get_Q(i,l);
			double alpha_i = alpha[i];
			for(int j=active_size;j<l;j++)
				G[j] += alpha_i * Q_i[j];
		}
}

void Solver::Solve(int l, const Kernel& Q, const double *b_, const schar *y_,
		   double *alpha_, double C, double eps,
		   double& obj, double& rho, int shrinking)
{
	this->l = l;
	this->Q = &Q;
	clone(b, b_,l);
	clone(y, y_,l);
	clone(alpha,alpha_,l);
	this->C = C;
	this->eps = eps;
	unshrinked = false;

	// initialize alpha_status
	{
		alpha_status = new char[l];
		for(int i=0;i<l;i++)
			update_alpha_status(i);
	}

	// initialize active set (for shrinking)
	{
		active_set = new int[l];
		for(int i=0;i<l;i++)
			active_set[i] = i;
		active_size = l;
	}

	// initialize gradient
	{
		G = new double[l];
		G_bar = new double[l];
		int i;
		for(i=0;i<l;i++)
		{
			G[i] = b[i];
			G_bar[i] = 0;
		}
		for(i=0;i<l;i++)
			if(!is_lower_bound(i))
			{
				double *Q_i = Q.get_Q(i,l);
				double alpha_i = alpha[i];
				int j;
				for(j=0;j<l;j++)
					G[j] += alpha_i*Q_i[j];
				if(is_upper_bound(i))
					for(j=0;j<l;j++)
						G_bar[j] += Q_i[j];
			}
	}

	// optimization step

	int iter = 0;
	int counter = min(l,1000)+1;

	while(1)
	{
		// show progress and do shrinking

		if(--counter == 0)
		{
			counter = min(l,1000);
			if(shrinking) do_shrinking();
			fprintf(stderr,".");
		}

		int i,j;
		if(select_working_set(i,j)!=0)
		{
			// reconstruct the whole gradient
			reconstruct_gradient();
			// reset active set size and check
			active_size = l;
			fprintf(stderr,"*");
			if(select_working_set(i,j)!=0)
				break;
			else
				counter = 1;	// do shrinking next iteration
		}
		
		++iter;

		// update alpha[i] and alpha[j]
		
		const double *Q_i = Q.get_Q(i,active_size);
		const double *Q_j = Q.get_Q(j,active_size);

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];

		if(y[i]!=y[j])
		{
			double L = max(0.0,alpha[j]-alpha[i]);
			double H = min(C,C+alpha[j]-alpha[i]);
			alpha[j] += (-G[i]-G[j])/(Q_i[i]+Q_j[j]+2*Q_i[j]);
			if(alpha[j] >= H) alpha[j] = H;
			else if(alpha[j] <= L) alpha[j] = L;
			alpha[i] += (alpha[j] - old_alpha_j);
		}
		else
		{
			double L = max(0.0,alpha[i]+alpha[j]-C);
			double H = min(C,alpha[i]+alpha[j]);
			alpha[j] += (G[i]-G[j])/(Q_i[i]+Q_j[j]-2*Q_i[j]);
			if(alpha[j] >= H) alpha[j] = H;
			else if(alpha[j] <= L) alpha[j] = L;			
			alpha[i] -= (alpha[j] - old_alpha_j);
		}

		// update alpha_status and G_bar

		{
			bool ui = is_upper_bound(i);
			bool uj = is_upper_bound(j);
			update_alpha_status(i);
			update_alpha_status(j);
			int k;
			if(ui != is_upper_bound(i))
			{
				Q_i = Q.get_Q(i,l);
				if(ui)
					for(k=0;k<l;k++)
						G_bar[k] -= Q_i[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += Q_i[k];
			}

			if(uj != is_upper_bound(j))
			{
				Q_j = Q.get_Q(j,l);
				if(uj)
					for(k=0;k<l;k++)
						G_bar[k] -= Q_j[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += Q_j[k];
			}
		}

		// update G

		double delta_alpha_i = alpha[i] - old_alpha_i;
		double delta_alpha_j = alpha[j] - old_alpha_j;
		
		for(int k=0;k<active_size;k++)
		{
			G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
		}
	}

	// calculate rho

	rho = calculate_rho();

	// calculate objective value
	{
		double v = 0;
		int i;
		for(i=0;i<l;i++)
			v += alpha[i] * (G[i] + b[i]);

		obj = v/2;
	}

	// put back the solution
	{
		for(int i=0;i<l;i++)
			alpha_[active_set[i]] = alpha[i];
	}

	// juggle everything back
	/*{
		for(int i=0;i<l;i++)
			while(active_set[i] != i)
				swap_index(i,active_set[i]);
				// or Q.swap_index(i,active_set[i]);
	}*/

	printf("\noptimization finished, #iter = %d\n",iter);

	delete[] b;
	delete[] y;
	delete[] alpha;
	delete[] alpha_status;
	delete[] active_set;
	delete[] G;
	delete[] G_bar;
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

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)	// y = +1
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
		else		// y = -1
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

void Solver::do_shrinking()
{
	int i,j,k;
	if(select_working_set(i,j)!=0) return;
	double Gm1 = -y[j]*G[j];
	double Gm2 = y[i]*G[i];

	// shrink
	
	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm1) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || -(Gm1 + Gm2) > eps*10) return;
	
	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm1) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double Solver::calculate_rho()
{
	double r;
	int nr_free = 0;
	double ub = INF, lb = -INF, sum_free = 0;
	for(int i=0;i<active_size;i++)
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

	return r;
}

//
// Solver for nu-svm classification
//

class Solver_NU_SVC : public Solver
{
public:
	Solver_NU_SVC() {}
	void Solve(int l, const Kernel& Q, const double *b, const schar *y_,
		   double *alpha_, double C_, double eps_,
		   double& obj, double& rho, const schar *y1_, double& r_,
		   int shrinking)
	{
		clone(y1,y1_,l);
		r = &r_;
		Solver::Solve(l,Q,b,y_,alpha_,C_,eps_,obj,rho,shrinking);
	}
	~Solver_NU_SVC()
	{
		delete[] y1;
	}
private:
	int select_working_set(int &i, int &j);
	double calculate_rho();
	void do_shrinking();
	schar *y1;
	double* r;
};

int Solver_NU_SVC::select_working_set(int &out_i, int &out_j)
{
	// return i,j which maximize -grad(f)^T d , under constraint
	// if alpha_i == C, d != +1
	// if alpha_i == 0, d != -1

	double Gmax1 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = +1 }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = -1 }
	int Gmax2_idx = -1;

	double Gmax3 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = +1 }
	int Gmax3_idx = -1;

	double Gmax4 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = -1 }
	int Gmax4_idx = -1;

	for(int i=0;i<active_size;i++)
	{
		if(y1[i]==+1)	// y == +1
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
		else		// y == -1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax3)
				{
					Gmax3 = -G[i];
					Gmax3_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax4)
				{
					Gmax4 = G[i];
					Gmax4_idx = i;
				}
			}
		}
	}

	if(max(Gmax1+Gmax2,Gmax3+Gmax4) < eps)
 		return 1;

	if(Gmax1+Gmax2 > Gmax3+Gmax4)
	{
		out_i = Gmax1_idx;
		out_j = Gmax2_idx;
	}
	else
	{
		out_i = Gmax3_idx;
		out_j = Gmax4_idx;
	}
	return 0;
}

void Solver_NU_SVC::do_shrinking()
{
	double Gmax1 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = +1 }
	double Gmax2 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = -1 }
	double Gmax3 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = +1 }
	double Gmax4 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = -1 }

	int k;
	for(k=0;k<active_size;k++)
	{
		if(!is_upper_bound(k))
		{
			if(y1[k]==+1)
			{
				if(-G[k] > Gmax1) Gmax1 = -G[k];
			}
			else	if(-G[k] > Gmax3) Gmax3 = -G[k];
		}
		if(!is_lower_bound(k))
		{
			if(y1[k]==+1)
			{	
				if(G[k] > Gmax2) Gmax2 = G[k];
			}
			else	if(G[k] > Gmax4) Gmax4 = G[k];
		}
	}

	double Gm1 = -Gmax2;
	double Gm2 = -Gmax1;
	double Gm3 = -Gmax4;
	double Gm4 = -Gmax3;

	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y1[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y1[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm4) continue;
		}
		else continue;

		--active_size;
		swap(y1[k],y1[active_size]);
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || max(-(Gm1+Gm2),-(Gm3+Gm4)) > eps*10) return;
	
	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm4) continue;
		}
		else continue;

		swap(y1[k],y1[active_size]);
		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double Solver_NU_SVC::calculate_rho()
{
	int nr_free1 = 0,nr_free2 = 0;
	double ub1 = INF, ub2 = INF;
	double lb1 = -INF, lb2 = -INF;
	double sum_free1 = 0, sum_free2 = 0;

	for(int i=0;i<active_size;i++)
	{
		if(y1[i]==+1)
		{
			if(is_lower_bound(i))
				ub1 = min(ub1,G[i]);
			else if(is_upper_bound(i))
				lb1 = max(lb1,G[i]);
			else
			{
				++nr_free1;
				sum_free1 += G[i];
			}
		}
		else
		{
			if(is_lower_bound(i))
				ub2 = min(ub2,G[i]);
			else if(is_upper_bound(i))
				lb2 = max(lb2,G[i]);
			else
			{
				++nr_free2;
				sum_free2 += G[i];
			}
		}
	}

	double r1,r2;
	if(nr_free1 > 0)
		r1 = sum_free1/nr_free1;
	else
		r1 = (ub1+lb1)/2;
	
	if(nr_free2 > 0)
		r2 = sum_free2/nr_free2;
	else
		r2 = (ub2+lb2)/2;
	
	*r = (r1+r2)/2;
	return (r1-r2)/2;
}

//
// Q matrices for different formulations
//
class C_SVC_Q: public Kernel
{ 
public:
	C_SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y_)
	:Kernel(prob.l, prob.x, param)
	{
		clone(y,y_,prob.l);
		cache = new Cache(prob.l,(int)(param.cache_size*(1<<20)));
	}
	
	double *get_Q(int i, int len) const
	{
		double *data;
		int start;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			for(int j=start;j<len;j++)
				data[j] = y[i]*y[j]*(this->*kernel_function)(i,j);
		}
		return data;
	}

	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		Kernel::swap_index(i,j);
		swap(y[i],y[j]);
	}

	~C_SVC_Q()
	{
		delete[] y;
		delete cache;
	}
private:
	schar *y;
	Cache *cache;
};

#define NU_SVC_Q C_SVC_Q	/* same stuff */

class ONE_CLASS_Q: public Kernel
{
public:
	ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
	:Kernel(prob.l, prob.x, param)
	{
		cache = new Cache(prob.l,(int)(param.cache_size*(1<<20)));
	}
	
	double *get_Q(int i, int len) const
	{
		double *data;
		int start;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			for(int j=start;j<len;j++)
				data[j] = (this->*kernel_function)(i,j);
		}
		return data;
	}

	~ONE_CLASS_Q()
	{
		delete cache;
	}
private:
	Cache *cache;
};

class C_SVR_Q: public Kernel
{ 
public:
	C_SVR_Q(const svm_problem& prob, const svm_parameter& param)
	:Kernel(prob.l, prob.x, param)
	{
		int l = prob.l;
		cache = new Cache(2*l,(int)(param.cache_size*(1<<20)));
		sign = new schar[2*l];
		index = new int[2*l];
		for(int k=0;k<l;k++)
		{
			sign[k] = 1;
			sign[k+l] = -1;
			index[k] = k;
			index[k+l] = k;
		}
	}

	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		swap(sign[i],sign[j]);
		swap(index[i],index[j]);
	}
	
	double *get_Q(int i, int len) const
	{
		double *data;
		int start;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			schar si = sign[i];
			for(int j=start;j<len;j++)
			{
				double v = (this->*kernel_function)(index[i],index[j]);
				if(si == sign[j]) data[j] = v;
				else data[j] = -v;
			}
		}
		return data;
	}

	~C_SVR_Q()
	{
		delete cache;
		delete[] sign;
		delete[] index;
	}
private:
	Cache *cache;
	schar *sign;
	int *index;
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
	double *alpha, double& obj, double& rho, double& upper_bound)
{
	int l = prob->l;
	double *minus_ones = new double[l];
	schar *y = new schar[l];

	int i;

	for(i=0;i<l;i++)
	{
		alpha[i] = 0;
		minus_ones[i] = -1;
		if(prob->y[i] > 0) y[i] = +1; else y[i]=-1;
	}

	Solver s;
	s.Solve(l, C_SVC_Q(*prob,*param,y), minus_ones, y,
		alpha, param->C, param->eps, obj, rho, param->shrinking);

	double sum_alpha=0;
	for(i=0;i<l;i++)
		sum_alpha += alpha[i];

	printf("nu = %f\n", sum_alpha/(param->C*prob->l));

	for(i=0;i<l;i++)
		alpha[i] *= y[i];

	upper_bound = param->C;

	delete[] minus_ones;
	delete[] y;
}

static void solve_nu_svc(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, double& obj, double& rho, double& upper_bound)
{
	int i;
	int l = prob->l;
	double nu = param->nu;

	int y_pos = 0;
	int y_neg = 0;
	schar *y = new schar[l];

	for(i=0;i<l;i++)
		if(prob->y[i]>0)
		{
			y[i] = +1;
			++y_pos;
		}
		else
		{
			y[i] = -1;
			++y_neg;
		}

	if(nu < 0 || nu*l/2 > min(y_pos,y_neg))
	{
		fprintf(stderr,"specified nu is infeasible\n");
		exit(1);
	}

	double sum_pos = nu*l/2;
	double sum_neg = nu*l/2;

	for(i=0;i<l;i++)
		if(y[i] == +1)
		{
			alpha[i] = min(1.0,sum_pos);
			sum_pos -= alpha[i];
		}
		else
		{
			alpha[i] = min(1.0,sum_neg);
			sum_neg -= alpha[i];
		}

	double *zeros = new double[l];
	schar *ones = new schar[l];

	for(i=0;i<l;i++)
	{	
		zeros[i] = 0;
		ones[i] = 1;
	}

	double r;
	Solver_NU_SVC s;
	s.Solve(l, NU_SVC_Q(*prob,*param,y), zeros, ones,
		alpha, 1.0, param->eps, obj, rho, y, r, param->shrinking);

	printf("C = %f\n",1/r);

	for(i=0;i<l;i++)
		alpha[i] *= y[i]/r;

	rho /= r;
	obj /= (r*r);
	upper_bound = 1/r;

	delete[] y;
	delete[] zeros;
	delete[] ones;
}

static void solve_one_class(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, double& obj, double& rho, double& upper_bound)
{
	int l = prob->l;
	double *zeros = new double[l];
	schar *ones = new schar[l];
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

	Solver s;
	s.Solve(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
		alpha, 1.0, param->eps, obj, rho, param->shrinking);

	upper_bound = 1;

	delete[] zeros;
	delete[] ones;
}

static void solve_c_svr(
	const svm_problem *prob, const svm_parameter *param,
	double *alpha, double& obj, double& rho, double & upper_bound)
{
	int l = prob->l;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	schar *y = new schar[2*l];
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

	Solver s;
	s.Solve(2*l, C_SVR_Q(*prob,*param), linear_term, y,
		alpha2, param->C, param->eps, obj, rho, param->shrinking);

	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];

	upper_bound = param->C;

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
	double obj, rho, upper_bound;
	switch(param->svm_type)
	{
		case C_SVC:
			solve_c_svc(prob,param,alpha,obj,rho,upper_bound);
			break;
		case NU_SVC:
			solve_nu_svc(prob,param,alpha,obj,rho,upper_bound);
			break;
		case ONE_CLASS:
			solve_one_class(prob,param,alpha,obj,rho,upper_bound);
			break;
		case C_SVR:
			solve_c_svr(prob,param,alpha,obj,rho,upper_bound);
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
			if(fabs(alpha[i]) >= upper_bound - EPS_A)
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

double svm_classify(const svm_model *model, const svm_node *x)
{
	const int n = model->n;
	const double *sv_coef = model->sv_coef;

	double sum = 0;

	for(int i=0;i<n;i++)
		sum += sv_coef[i] * Kernel::k_function(x,model->SV[i],model->param);
	sum-=(model->rho);
	return sum;
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
