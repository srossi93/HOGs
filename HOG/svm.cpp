 /***************************************************************************
 *   Copyright (C) 2009 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI)    *
 *   by modifying the LIBSVM 2.89 package by C.-C. Chang and C.-J. Lin.    *
 *   Copyright information for LIBSVM 2.89 is below.                       *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *					          EXTENDED LIBSVM                              *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * Copyright (c) 2000-2009 Chih-Chung Chang and Chih-Jen Lin               *
 * All rights reserved.                                                    *
 *                                                                         *
 * Redistribution and use in source and binary forms, with or without      *
 * modification, are permitted provided that the following conditions      *
 * are met:                                                                *
 *                                                                         *
 * 1. Redistributions of source code must retain the above copyright       *
 * notice, this list of conditions and the following disclaimer.           *
 *                                                                         *
 * 2. Redistributions in binary form must reproduce the above copyright    *
 * notice, this list of conditions and the following disclaimer in the     *
 * documentation and/or other materials provided with the distribution.    *
 *                                                                         *
 * 3. Neither name of copyright holders nor the names of its contributors  *
 * may be used to endorse or promote products derived from this software   *
 * without specific prior written permission.                              *
 *                                                                         *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     *
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   *
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR  *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,   *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,     *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR      *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            *
 *                                                                         *
 ***************************************************************************/


#include "image.h"
#include "library.h"
#include "svm.h"


/***************************************************************************************************************************************/
int libsvm_version = LIBSVM_VERSION;


template <class T> inline T MIN(T x,T y) { return (x<y)?x:y; }


template <class T> inline T MAX(T x,T y) { return (x>y)?x:y; }


template <class T> inline void SWAP(T& x, T& y) { T t=x; x=y; y=t; }


template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
	dst = new T[n];
	memcpy((void *)dst,(void *)src,sizeof(T)*n);
}


inline double powi(double base, int times)
{
	double tmp = base, ret = 1.0;
	int t;
	for(t=times; t>0; t/=2)
	{
		if(t%2==1) ret*=tmp;
		tmp = tmp * tmp;
	}
	return ret;
}


#define INF HUGE_VAL
#define TAU 1e-12
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define FMAX(a,b) ((a)>(b) ? (a) : (b))


const char *svm_type_table[] = {"c_svc","one_class","epsilon_svr",NULL};
const char *kernel_type_table[] = {"linear","polynomial","rbf","sigmoid",NULL};


/***************************************************************************************************************************************/
class Cache
{
public:
	Cache(int l,long int size);
	~Cache();
	int get_data(const int index, float **data, int len);
	void swap_index(int i, int j);	
private:
	int l;
	long int size;
	struct head_t
	{
		head_t *prev, *next;	// a circular list
		float *data;
		int len;				// data[0,len) is cached in this entry
	};
	head_t *head;
	head_t lru_head;
	void lru_delete(head_t *h);
	void lru_insert(head_t *h);
};


Cache::Cache(int l_,long int size_):l(l_),size(size_)
{
	head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
	size /= sizeof(float);
	size -= l * sizeof(head_t) / sizeof(float);
	size = MAX(size, 2 * (long int) l);	// cache must be large enough for two columns
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


int Cache::get_data(const int index, float **data, int len)
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
		h->data = (float *)realloc(h->data,sizeof(float)*len);
		size -= more;
		SWAP(h->len,len);
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
	SWAP(head[i].data,head[j].data);
	SWAP(head[i].len,head[j].len);
	if(head[i].len) lru_insert(&head[i]);
	if(head[j].len) lru_insert(&head[j]);
	if(i>j) SWAP(i,j);
	for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
	{
		if(h->len > i)
		{
			if(h->len > j)
				SWAP(h->data[i],h->data[j]);
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


/***************************************************************************************************************************************/
class QMatrix
{
public:
	virtual float *get_Q(int column, int len) const = 0;
	virtual float *get_QD() const = 0;
	virtual void swap_index(int i, int j) const = 0;
	virtual ~QMatrix() {}
};


class Kernel: public QMatrix
{
public:
	Kernel(int l, svm_node * x, const svm_parameter& param);
	virtual ~Kernel();
	static double k_function(const svm_node *x, const svm_node *y, const svm_parameter& param);
	virtual float *get_Q(int column, int len) const = 0;
	virtual float *get_QD() const = 0;
	virtual void swap_index(int i, int j) const	// no so const...
	{
		SWAP(x[i],x[j]);
		if(x_square) SWAP(x_square[i],x_square[j]);
	}
protected:
	double (Kernel::*kernel_function)(int i, int j) const;
private:
	svm_node *x;
	double *x_square;
	const int kernel_type;
	const int degree;
	const double gamma;
	const double coef0;
	const double lambda;
	static double dot(const svm_node *px, const svm_node *py);
	static double dot(const svm_node &px, const svm_node &py);
	double kernel_linear(int i, int j) const
	{
		return dot(x[i],x[j]);
	}
	double kernel_poly(int i, int j) const
	{
		return powi(gamma*dot(x[i],x[j])+coef0,degree);
	}
	double kernel_rbf(int i, int j) const
	{
		return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
	}
	double kernel_sigmoid(int i, int j) const
	{
		return tanh(gamma*dot(x[i],x[j])+coef0);
	}
	double kernel_mrf(int i, int j) const
	{
		double sum=0;
		int r,last=x[i].dim-1;
		for(r=0;r<last;r++)
			sum+=(x[i].values[r]-x[j].values[r])*(x[i].values[r]-x[j].values[r]);
		return exp(-gamma*sum)+lambda*x[i].values[last]*x[j].values[last];
	}
	double kernel_gustavo(int i, int j) const
	{
		double spectral=0,spatial=0;
		int r,dim=x[i].dim/3;
		for(r=0;r<dim;r++)
			spectral+=(x[i].values[r]-x[j].values[r])*(x[i].values[r]-x[j].values[r]);
		for(r=dim;r<3*dim;r++)
			spatial+=(x[i].values[r]-x[j].values[r])*(x[i].values[r]-x[j].values[r]);
		return lambda*exp(-coef0*spatial)+(1-lambda)*exp(-gamma*spectral);
	}
};


Kernel::Kernel(int l, svm_node * x_, const svm_parameter& param):kernel_type(param.kernel_type), degree(param.degree),gamma(param.gamma), coef0(param.coef0), lambda(param.lambda)
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
		case MRF:
			kernel_function = &Kernel::kernel_mrf;
			break;
		case GUSTAVO:
			kernel_function = &Kernel::kernel_gustavo;
			break;
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
	int dim = MIN(px->dim, py->dim);
	for (int i = 0; i < dim; i++)
		sum += (px->values)[i] * (py->values)[i];
	return sum;
}


double Kernel::dot(const svm_node &px, const svm_node &py)
{
	double sum = 0;
	int dim = MIN(px.dim, py.dim);
	for (int i = 0; i < dim; i++)
		sum += px.values[i] * py.values[i];
	return sum;
}


double Kernel::k_function(const svm_node *x, const svm_node *y, const svm_parameter& param)
{
	switch(param.kernel_type)
	{
		case LINEAR:
			return dot(x,y);
		case POLY:
			return powi(param.gamma*dot(x,y)+param.coef0,param.degree);
		case RBF:
		{
			double sum = 0;
			int dim = MIN(x->dim, y->dim), i;
			for (i = 0; i < dim; i++)
			{
				double d = x->values[i] - y->values[i];
				sum += d*d;
			}
			for (; i < x->dim; i++)
				sum += x->values[i] * x->values[i];
			for (; i < y->dim; i++)
				sum += y->values[i] * y->values[i];
			return exp(-param.gamma*sum);
		}
		case SIGMOID:
			return tanh(param.gamma*dot(x,y)+param.coef0);
		case MRF:
		{
			double sum=0;
			word last=x->dim-1,i; // x and y are assumed to have the same dimensions
			for(i=0;i<last;i++)
				sum+=(x->values[i]-y->values[i])*(x->values[i]-y->values[i]);
			return exp(-param.gamma*sum)+param.lambda*x->values[last]*y->values[last];
		}
		case GUSTAVO:
		{
			double spectral=0,spatial=0;
			word dim=x->dim/3,i; // x and y are assumed to have the same dimensions
			for(i=0;i<dim;i++)
				spectral+=(x->values[i]-y->values[i])*(x->values[i]-y->values[i]);
			for(i=dim;i<3*dim;i++)
				spatial+=(x->values[i]-y->values[i])*(x->values[i]-y->values[i]);
			return param.lambda*exp(-param.coef0*spatial)+(1-param.lambda)*exp(-param.gamma*spectral);
		}
		default:
			return 0;  // Unreachable 
	}
}


/***************************************************************************************************************************************/
// An SMO algorithm in Fan et al., JMLR 6(2005), p. 1889--1918
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + p^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//	Q, p, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping tolerance
//
// solution will be put in \alpha, objective value will be put in obj

class Solver
{
public:
	Solver() {};
	virtual ~Solver() {};

	struct SolutionInfo
	{
		double obj;
		double rho;
		double upper_bound_p;
		double upper_bound_n;
		double r;	// for Solver_NU
		longword nSV;
	};
	void Solve(int l, const QMatrix& Q, const double *p_, const signed char *y_, double *alpha_, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking);
protected:
	int active_size;
	signed char *y;
	double *G;		// gradient of objective function
	enum { LOWER_BOUND, UPPER_BOUND, FREE };
	char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	double *alpha;
	const QMatrix *Q;
	const float *QD;
	double eps;
	double Cp,Cn;
	double *p;
	int *active_set;
	double *G_bar;		// gradient, if we treat free variables as 0
	int l;
	bool unshrink;	// XXX
	double get_C(int i)
	{
		return (y[i] > 0)? Cp : Cn;
	}
	void update_alpha_status(int i)
	{
		if(alpha[i] >= get_C(i))
			alpha_status[i] = UPPER_BOUND;
		else if(alpha[i] <= 0)
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
private:
	bool be_shrunk(int i, double Gmax1, double Gmax2);	
};


void Solver::swap_index(int i, int j)
{
	Q->swap_index(i,j);
	SWAP(y[i],y[j]);
	SWAP(G[i],G[j]);
	SWAP(alpha_status[i],alpha_status[j]);
	SWAP(alpha[i],alpha[j]);
	SWAP(p[i],p[j]);
	SWAP(active_set[i],active_set[j]);
	SWAP(G_bar[i],G_bar[j]);
}


void Solver::reconstruct_gradient()
{
	// reconstruct inactive elements of G from G_bar and free variables
	if(active_size == l) return;

	int i,j;
	int nr_free = 0;

	for(j=active_size;j<l;j++)
		G[j] = G_bar[j] + p[j];

	for(j=0;j<active_size;j++)
		if(is_free(j))
			nr_free++;
	if (nr_free*l > 2*active_size*(l-active_size))
	{
		for(i=active_size;i<l;i++)
		{
			const float *Q_i = Q->get_Q(i,active_size);
			for(j=0;j<active_size;j++)
				if(is_free(j))
					G[i] += alpha[j] * Q_i[j];
		}
	}
	else
	{
		for(i=0;i<active_size;i++)
			if(is_free(i))
			{
				const float *Q_i = Q->get_Q(i,l);
				double alpha_i = alpha[i];
				for(j=active_size;j<l;j++)
					G[j] += alpha_i * Q_i[j];
			}
	}
}


void Solver::Solve(int l, const QMatrix& Q, const double *p_, const signed char *y_, double *alpha_, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking)
{
	this->l = l;
	this->Q = &Q;
	QD=Q.get_QD();
	clone(p, p_,l);
	clone(y, y_,l);
	clone(alpha,alpha_,l);
	this->Cp = Cp;
	this->Cn = Cn;
	this->eps = eps;
	unshrink = false;
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
			G[i] = p[i];
			G_bar[i] = 0;
		}
		for(i=0;i<l;i++)
			if(!is_lower_bound(i))
			{
				const float *Q_i = Q.get_Q(i,l);
				double alpha_i = alpha[i];
				int j;
				for(j=0;j<l;j++)
					G[j] += alpha_i*Q_i[j];
				if(is_upper_bound(i))
					for(j=0;j<l;j++)
						G_bar[j] += get_C(i) * Q_i[j];
			}
	}
	// optimization step
	int iter = 0;
	int counter = MIN(l,1000)+1;
	while(1)
	{
		// show progress and do shrinking
		if(--counter == 0)
		{
			counter = MIN(l,1000);
			if(shrinking) do_shrinking();
			cout << ".";
			flush(cout);
		}
		int i,j;
		if(select_working_set(i,j)!=0)
		{
			// reconstruct the whole gradient
			reconstruct_gradient();
			// reset active set size and check
			active_size = l;
			cout << "*";
			flush(cout);
			if(select_working_set(i,j)!=0)
				break;
			else
				counter = 1;	// do shrinking next iteration
		}
		++iter;
		// update alpha[i] and alpha[j], handle bounds carefully
		const float *Q_i = Q.get_Q(i,active_size);
		const float *Q_j = Q.get_Q(j,active_size);
		double C_i = get_C(i);
		double C_j = get_C(j);
		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];
		if(y[i]!=y[j])
		{
			double quad_coef = Q_i[i]+Q_j[j]+2*Q_i[j];
			if (quad_coef <= 0)
				quad_coef = TAU;
			double delta = (-G[i]-G[j])/quad_coef;
			double diff = alpha[i] - alpha[j];
			alpha[i] += delta;
			alpha[j] += delta;
			if(diff > 0)
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = diff;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = -diff;
				}
			}
			if(diff > C_i - C_j)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = C_i - diff;
				}
			}
			else
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = C_j + diff;
				}
			}
		}
		else
		{
			double quad_coef = Q_i[i]+Q_j[j]-2*Q_i[j];
			if (quad_coef <= 0)
				quad_coef = TAU;
			double delta = (G[i]-G[j])/quad_coef;
			double sum = alpha[i] + alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;

			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum - C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum - C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
		}
		// update G
		double delta_alpha_i = alpha[i] - old_alpha_i;
		double delta_alpha_j = alpha[j] - old_alpha_j;
		for(int k=0;k<active_size;k++)
		{
			G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
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
						G_bar[k] -= C_i * Q_i[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_i * Q_i[k];
			}

			if(uj != is_upper_bound(j))
			{
				Q_j = Q.get_Q(j,l);
				if(uj)
					for(k=0;k<l;k++)
						G_bar[k] -= C_j * Q_j[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_j * Q_j[k];
			}
		}
	}
	// calculate rho
	si->rho = calculate_rho();
	// calculate objective value
	{
		double v = 0;
		int i;
		for(i=0;i<l;i++)
			v += alpha[i] * (G[i] + p[i]);

		si->obj = v/2;
	}
	// put back the solution
	{
		for(int i=0;i<l;i++)
			alpha_[active_set[i]] = alpha[i];
	}
	// juggle everything back
	si->upper_bound_p = Cp;
	si->upper_bound_n = Cn;
	cout << endl << "optimization finished, #iter = " << iter << endl;
	flush(cout);
	delete[] p;
	delete[] y;
	delete[] alpha;
	delete[] alpha_status;
	delete[] active_set;
	delete[] G;
	delete[] G_bar;
}


int Solver::select_working_set(int &out_i, int &out_j)				// return 1 if already optimal, return 0 otherwise
{
	// return i,j such that
	// i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
	// j: minimizes the decrease of obj value
	//    (if quadratic coefficeint <= 0, replace it with tau)
	//    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
	double Gmax = -INF;
	double Gmax2 = -INF;
	int Gmax_idx = -1;
	int Gmin_idx = -1;
	double obj_diff_min = INF;
	for(int t=0;t<active_size;t++)
		if(y[t]==+1)	
		{
			if(!is_upper_bound(t))
				if(-G[t] >= Gmax)
				{
					Gmax = -G[t];
					Gmax_idx = t;
				}
		}
		else
		{
			if(!is_lower_bound(t))
				if(G[t] >= Gmax)
				{
					Gmax = G[t];
					Gmax_idx = t;
				}
		}
	int i = Gmax_idx;
	const float *Q_i = NULL;
	if(i != -1) // NULL Q_i not accessed: Gmax=-INF if i=-1
		Q_i = Q->get_Q(i,active_size);
	for(int j=0;j<active_size;j++)
	{
		if(y[j]==+1)
		{
			if (!is_lower_bound(j))
			{
				double grad_diff=Gmax+G[j];
				if (G[j] >= Gmax2)
					Gmax2 = G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef=Q_i[i]+QD[j]-2.0*y[i]*Q_i[j];
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						Gmin_idx=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
		else
		{
			if (!is_upper_bound(j))
			{
				double grad_diff= Gmax-G[j];
				if (-G[j] >= Gmax2)
					Gmax2 = -G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef=Q_i[i]+QD[j]+2.0*y[i]*Q_i[j];
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						Gmin_idx=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
	}
	if(Gmax+Gmax2 < eps)
		return 1;

	out_i = Gmax_idx;
	out_j = Gmin_idx;
	return 0;
}


bool Solver::be_shrunk(int i, double Gmax1, double Gmax2)
{
	if(is_upper_bound(i))
	{
		if(y[i]==+1)
			return(-G[i] > Gmax1);
		else
			return(-G[i] > Gmax2);
	}
	else if(is_lower_bound(i))
	{
		if(y[i]==+1)
			return(G[i] > Gmax2);
		else	
			return(G[i] > Gmax1);
	}
	else
		return(false);
}


void Solver::do_shrinking()
{
	int i;
	double Gmax1 = -INF;		// MAX { -y_i * grad(f)_i | i in I_up(\alpha) }
	double Gmax2 = -INF;		// MAX { y_i * grad(f)_i | i in I_low(\alpha) }
	// find maximal violating pair first
	for(i=0;i<active_size;i++)
	{
		if(y[i]==+1)	
		{
			if(!is_upper_bound(i))	
			{
				if(-G[i] >= Gmax1)
					Gmax1 = -G[i];
			}
			if(!is_lower_bound(i))	
			{
				if(G[i] >= Gmax2)
					Gmax2 = G[i];
			}
		}
		else	
		{
			if(!is_upper_bound(i))	
			{
				if(-G[i] >= Gmax2)
					Gmax2 = -G[i];
			}
			if(!is_lower_bound(i))	
			{
				if(G[i] >= Gmax1)
					Gmax1 = G[i];
			}
		}
	}
	if(unshrink == false && Gmax1 + Gmax2 <= eps*10) 
	{
		unshrink = true;
		reconstruct_gradient();
		active_size = l;
		cout << "*";
		flush(cout);
	}
	for(i=0;i<active_size;i++)
		if (be_shrunk(i, Gmax1, Gmax2))
		{
			active_size--;
			while (active_size > i)
			{
				if (!be_shrunk(active_size, Gmax1, Gmax2))
				{
					swap_index(i,active_size);
					break;
				}
				active_size--;
			}
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

		if(is_upper_bound(i))
		{
			if(y[i]==-1)
				ub = MIN(ub,yG);
			else
				lb = MAX(lb,yG);
		}
		else if(is_lower_bound(i))
		{
			if(y[i]==+1)
				ub = MIN(ub,yG);
			else
				lb = MAX(lb,yG);
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


/***************************************************************************************************************************************/
class SVC_Q: public Kernel
{ 
public:
	SVC_Q(const svm_problem& prob, const svm_parameter& param, const signed char *y_):Kernel(prob.l, prob.x, param)
	{
		clone(y,y_,prob.l);
		cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
		QD = new float[prob.l];
		for(unsigned long i=0;i<prob.l;i++)
			QD[i]= (float)(this->*kernel_function)(i,i);
	}	
	float *get_Q(int i, int len) const
	{
		float *data;
		int start, j;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			for(j=start;j<len;j++)
				data[j] = (float)(y[i]*y[j]*(this->*kernel_function)(i,j));
		}
		return data;
	}
	float *get_QD() const
	{
		return QD;
	}
	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		Kernel::swap_index(i,j);
		SWAP(y[i],y[j]);
		SWAP(QD[i],QD[j]);
	}
	~SVC_Q()
	{
		delete[] y;
		delete cache;
		delete[] QD;
	}
private:
	signed char *y;
	Cache *cache;
	float *QD;
};


class ONE_CLASS_Q: public Kernel
{
public:
	ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param):Kernel(prob.l, prob.x, param)
	{
		cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
		QD = new float[prob.l];
		for(unsigned long i=0;i<prob.l;i++)
			QD[i]= (float)(this->*kernel_function)(i,i);
	}	
	float *get_Q(int i, int len) const
	{
		float *data;
		int start, j;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			for(j=start;j<len;j++)
				data[j] = (float)(this->*kernel_function)(i,j);
		}
		return data;
	}
	float *get_QD() const
	{
		return QD;
	}
	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		Kernel::swap_index(i,j);
		SWAP(QD[i],QD[j]);
	}
	~ONE_CLASS_Q()
	{
		delete cache;
		delete[] QD;
	}
private:
	Cache *cache;
	float *QD;
};


class SVR_Q: public Kernel
{ 
public:
	SVR_Q(const svm_problem& prob, const svm_parameter& param):Kernel(prob.l, prob.x, param)
	{
		l = prob.l;
		cache = new Cache(l,(long int)(param.cache_size*(1<<20)));
		QD = new float[2*l];
		sign = new signed char[2*l];
		index = new int[2*l];
		for(int k=0;k<l;k++)
		{
			sign[k] = 1;
			sign[k+l] = -1;
			index[k] = k;
			index[k+l] = k;
			QD[k]= (float)(this->*kernel_function)(k,k);
			QD[k+l]=QD[k];
		}
		buffer[0] = new float[2*l];
		buffer[1] = new float[2*l];
		next_buffer = 0;
	}
	void swap_index(int i, int j) const
	{
		SWAP(sign[i],sign[j]);
		SWAP(index[i],index[j]);
		SWAP(QD[i],QD[j]);
	}	
	float *get_Q(int i, int len) const
	{
		float *data;
		int j, real_i = index[i];
		if(cache->get_data(real_i,&data,l) < l)
		{
			for(j=0;j<l;j++)
				data[j] = (float)(this->*kernel_function)(real_i,j);
		}
		// reorder and copy
		float *buf = buffer[next_buffer];
		next_buffer = 1 - next_buffer;
		signed char si = sign[i];
		for(j=0;j<len;j++)
			buf[j] = (float) si * (float) sign[j] * data[index[j]];
		return buf;
	}
	float *get_QD() const
	{
		return QD;
	}
	~SVR_Q()
	{
		delete cache;
		delete[] sign;
		delete[] index;
		delete[] buffer[0];
		delete[] buffer[1];
		delete[] QD;
	}
private:
	int l;
	Cache *cache;
	signed char *sign;
	int *index;
	mutable int next_buffer;
	float *buffer[2];
	float *QD;
};


/***************************************************************************************************************************************/
static void solve_c_svc(const svm_problem *prob, const svm_parameter* param, double *alpha, Solver::SolutionInfo* si, double Cp, double Cn)
{
	unsigned long l = prob->l,i;
	double *minus_ones = new double[l];
	signed char *y = new signed char[l];
	for(i=0;i<l;i++)
	{
		alpha[i] = 0;
		minus_ones[i] = -1;
		y[i] = (prob->y[i]>0 ? 1 : -1);
	}
	Solver s;
	s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y, alpha, Cp, Cn, param->eps, si, param->shrinking);
	si->nSV=0;
	for(i=0;i<l;i++)
	{
		alpha[i] *= y[i];
		if(fabs(alpha[i])>0)
			++si->nSV;
	}
	delete[] minus_ones,y;
}


static void solve_epsilon_svr(const svm_problem *prob, const svm_parameter *param,	double *alpha, Solver::SolutionInfo* si)
{
	unsigned long l = prob->l,i;
	double *alpha2 = new double[2*l], *linear_term = new double[2*l];
	signed char *y = new signed char[2*l];
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
	s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y, alpha2, param->C, param->C, param->eps, si, param->shrinking);
	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];
	delete[] alpha2,linear_term,y;
}


static void solve_radius_L1(const svm_problem *prob, const svm_parameter *param, double *beta, Solver::SolutionInfo* si)
{
	unsigned long l = prob->l,i;
	double *zeros = new double[l];
	signed char *ones = new signed char[l];
	beta[0] = 1;		// primo beta a 1, 
	for(i=1;i<l;i++)
		beta[i] = 0;	// gli altri a 0
	for(i=0;i<l;i++)
	{
		zeros[i] = 0;
		ones[i] = 1;
	}
	Solver s;
	s.Solve(l,ONE_CLASS_Q(*prob,*param),zeros,ones,beta,INF,INF,param->eps,si,param->shrinking);
	delete[] zeros,ones;
}


/***************************************************************************************************************************************/
static double span_bound_compute_svc(const svm_problem *prob, const svm_parameter* param, double *alpha, Solver::SolutionInfo* si)
{
	longword i;
	vector<longword> usv,bsv;
	for(i=0;i<prob->l;i++)
		if(((alpha[i]>0)&&(fabs(alpha[i])<si->upper_bound_p)) || ((alpha[i]<0)&&(fabs(alpha[i])<si->upper_bound_n)))
			usv.push_back(i);
		else if(fabs(alpha[i])>0)
			bsv.push_back(i);
	longword nUSV=usv.size(),nBSV=bsv.size(),b,u,spancounter=0;
	if(nUSV>0)
	{
		Matrix QU = Matrix(nUSV,nUSV,1.1), QUB = Matrix(nUSV,nBSV);
		Vector c = Vector(nUSV,0), q = Vector(nUSV);
		for(u=0;u<nUSV;u++)
		{
			for(i=u+1;i<nUSV;i++)
			{
				QU[u][i]=Kernel::k_function(&prob->x[usv[u]],&prob->x[usv[i]],*param);
				QU[i][u]=QU[u][i];
			}
			for(b=0;b<nBSV;b++)
				QUB[u][b]=Kernel::k_function(&prob->x[usv[u]],&prob->x[bsv[b]],*param);
		}
		Matrix L=InvLowerCholesky(QU);
		double a=0;
		for(u=0;u<nUSV;u++)
		{
			for(i=0;i<=u;i++)
				c[u]+=L[u][i];
			a+=c[u]*c[u];
			QU[u][u]=1;
		}
		if(a!=0)
			a=-1/a;
		else
			errorlog("span_bound_compute_svc. Singular lower diagonal matrix.");
		for(u=0;u<nUSV;u++)
		{
			double alphaspan=0,d=0,discriminant=si->rho;
			for(i=u;i<nUSV;i++)
			{
				alphaspan+=L[i][u]*L[i][u];
				d+=L[i][u]*c[i];
				discriminant+=alpha[usv[i]]*QU[i][u];
			}
			alphaspan = (alphaspan+a*d*d>0 ? fabs(alpha[usv[u]])/(alphaspan+a*d*d) : 0);
			for(b=0;b<nBSV;b++)
				discriminant+=alpha[bsv[b]]*QUB[u][b];
			if(alphaspan>prob->y[usv[u]]*discriminant)
				spancounter++;
		}
		for(b=0;b<nBSV;b++)
		{
			double alphaspan=1-a,f=0,discriminant=si->rho;
			for(u=0;u<nUSV;u++)
			{
				double e=0;
				for(i=0;i<=u;i++)
					e+=L[u][i]*QUB[i][b];
				alphaspan+=e*e;
				f+=e*c[u];
				discriminant+=alpha[usv[u]]*QUB[u][b];
			}
			alphaspan = (alphaspan+a*f*(2+f)>0 ? fabs(alpha[bsv[b]])*(alphaspan+a*f*(2+f)) : 0);
			for(i=0;i<nBSV;i++)
				discriminant+=alpha[bsv[i]]*Kernel::k_function(&prob->x[bsv[i]],&prob->x[bsv[b]],*param);
			if(alphaspan>prob->y[bsv[b]]*discriminant)
				spancounter++;
		}
		return float(spancounter)/prob->l;
	}
	else if(si->nSV>0)
		return float(si->nSV)/prob->l;
	else 
		return 1;
}


static double span_bound_compute_svr(const svm_problem *prob, const svm_parameter* param, double *beta, Solver::SolutionInfo* si)
{
	flush(cout);
	longword i,j;
	vector<longword> usv,bsv;		// determino gli indici degli usv e bsv
	for(i=0;i<prob->l;i++)
		if(fabs(beta[i])==param->C)
			bsv.push_back(i);
		else if(fabs(beta[i])>0)
			usv.push_back(i);
	longword nUSV=usv.size(),nBSV=bsv.size(),b,u;
	cout << nUSV << " " << nBSV << " ";
	if(nUSV>0)
	{
		double spanbound=0;
		for(i=0;i<prob->l;i++)	
			if(fabs(beta[i])>0)
			{
				spanbound+=(prob->y[i]*beta[i]-param->p*fabs(beta[i])-Kernel::k_function(&prob->x[i],&prob->x[i],*param)*beta[i]*beta[i]);
				for(j=i+1;j<prob->l;j++)
					if(fabs(beta[j])>0)
						spanbound-=2*Kernel::k_function(&prob->x[i],&prob->x[j],*param)*beta[i]*beta[j];
			}	
		spanbound=param->p+spanbound/(param->C*prob->l);
		Matrix QU = Matrix(nUSV,nUSV,1.1), QUB = Matrix(nUSV,nBSV);
		Vector c = Vector(nUSV,0), q = Vector(nUSV);
		for(u=0;u<nUSV;u++)
		{
			for(i=u+1;i<nUSV;i++)
			{
				QU[u][i]=Kernel::k_function(&prob->x[usv[u]],&prob->x[usv[i]],*param);
				QU[i][u]=QU[u][i];
			}
			for(b=0;b<nBSV;b++)
				QUB[u][b]=Kernel::k_function(&prob->x[usv[u]],&prob->x[bsv[b]],*param);
		}
		Matrix L=InvLowerCholesky(QU);
		double a=0;
		for(u=0;u<nUSV;u++)
		{
			for(i=0;i<=u;i++)
				c[u]+=L[u][i];
			a+=c[u]*c[u];
			QU[u][u]=1;
		}
		if(a!=0)
			a=-1/a;
		else
			errorlog("span_bound_compute. Singular lower diagonal matrix.");
		for(u=0;u<nUSV;u++)
		{
			double betaspan=0,d=0;
			for(i=u;i<nUSV;i++)
			{
				betaspan+=L[i][u]*L[i][u];
				d+=L[i][u]*c[i];
			}
			betaspan = (betaspan+a*d*d>0 ? fabs(beta[usv[u]])/(betaspan+a*d*d) : 0);
			spanbound+=betaspan/prob->l;
		}
		for(b=0;b<nBSV;b++)
		{
			double betaspan=1-a,f=0;
			for(u=0;u<nUSV;u++)
			{
				double e=0;
				for(i=0;i<=u;i++)
					e+=L[u][i]*QUB[i][b];
				betaspan+=e*e;
				f+=e*c[u];
			}
			betaspan = (betaspan+a*f*(2+f)>0 ? fabs(beta[bsv[b]])*(betaspan+a*f*(2+f)) : 0);
			spanbound+=betaspan/prob->l;
		}
		return spanbound;
	}
	else if(powell_svr_tanh)
		return (prob->l+1)*param_max[0]+2*param_max[2]+1;
	else
		return (prob->l+1)*exp(theta_max[0]+1.0)+2*exp(theta_max[2]+1.0)+1;
}


/***************************************************************************************************************************************/
static void var_train_compute(const svm_problem *prob, const svm_parameter* param, double *beta, Vector &var_train, Matrix &L, vector<longword> &usv)
{
	longword i,j,u;
	for(i=0;i<prob->l;i++)
		if((fabs(beta[i])>0)&&(fabs(beta[i])<param->C))
			usv.push_back(i);
	longword nUSV=usv.size();
	cout << "#USV = " << nUSV << endl;
	var_train = Vector(prob->l);
	for(i=0;i<prob->l;i++)
		var_train[i]=Kernel::k_function(&prob->x[i],&prob->x[i],*param);
	if(nUSV>0)
	{
		//ofstream out("prova.xls");
		Matrix QU = Matrix(nUSV,nUSV,1.00001);	// era 1.00001
		for(u=0;u<nUSV;u++)
			for(i=u+1;i<nUSV;i++)
			{
				QU[u][i]=Kernel::k_function(&prob->x[usv[u]],&prob->x[usv[i]],*param);
				QU[i][u]=QU[u][i];
			}
		L=InvLowerCholesky(QU);
		for(j=0;j<prob->l;j++)
		{
			if(nUSV>2)
			{
				for(u=0;u<nUSV;u++)
				{
					double LQ=0;
					for(i=0;i<u+1;i++)
						LQ+=L[u][i]*Kernel::k_function(&prob->x[usv[i]],&prob->x[j],*param);
					var_train[j]-=LQ*LQ;
				}
				var_train[j]=(var_train[j]<0 ? 0 : var_train[j]);
			}
			else if(nUSV==2)
			{
				double q[2];
				for(u=0;u<nUSV;u++)
					q[u]=Kernel::k_function(&prob->x[usv[u]],&prob->x[i],*param);
				var_train[i]-=(QU[1][1]*q[0]*q[0]+QU[0][0]*q[1]*q[1]-2*QU[0][1]*q[0]*q[1])/(QU[0][0]*QU[1][1]-2*QU[0][1]*QU[0][1]);
			}
			else if(nUSV==1)
			{
				double q=Kernel::k_function(&prob->x[usv[0]],&prob->x[i],*param);
				var_train[i]-=q*q/QU[0][0];
			}
			//var_train[j]=fabs(var_train[j]);
		}
	}
}


/*static void var_train_compute(const svm_problem *prob, const svm_parameter* param, double *beta, Vector &var_train, Matrix &QU, vector<longword> &usv)
{
	longword i,j,u;
	for(i=0;i<prob->l;i++)
		if((fabs(beta[i])>0)&&(fabs(beta[i])<param->C))
			usv.push_back(i);
	longword nUSV=usv.size();
	cout << "#USV = " << nUSV << endl;
	var_train = Vector(prob->l,1); // così K(x,x) è già scritto in ogni componente
	if(nUSV>0)						// se ho usv
	{
		QU = Matrix(nUSV,nUSV,1);
		Vector q = Vector(nUSV);
		Matrix QQ = Matrix(nUSV,nUSV);
		for(u=0;u<nUSV;u++)
			for(j=u+1;j<nUSV;j++)
			{
				QU[u][j]=Kernel::k_function(&prob->x[usv[u]],&prob->x[usv[j]],*param);
				QU[j][u]=QU[u][j];
			}
		for(i=0;i<prob->l;i++)
			if(nUSV>2)
			{
				double qstar=Kernel::k_function(&prob->x[i],&prob->x[i],*param);
				for(u=0;u<nUSV;u++)
				{
					q[u]=Kernel::k_function(&prob->x[usv[u]],&prob->x[i],*param);
					QQ[u][u]=QU[u][u]-q[u]*q[u]/qstar;
				}
				for(u=0;u<nUSV;u++)
					for(j=u+1;j<nUSV;j++)
					{
						QQ[u][j]=QU[u][j]-q[u]*q[j]/qstar;
						QQ[j][u]=QQ[u][j];
					}
				var_train[i]=1/(1/qstar-2*ConjugateGrad(QQ,q)/(qstar*qstar));
				var_train[i]=(var_train[i]<0 ? 0 : var_train[i]);
			}
			else if(nUSV==2)
			{
				for(u=0;u<nUSV;u++)
					q[u]=Kernel::k_function(&prob->x[usv[u]],&prob->x[i],*param);
				var_train[i]=Kernel::k_function(&prob->x[i],&prob->x[i],*param)-(QU[1][1]*q[0]*q[0]+QU[0][0]*q[1]*q[1]-2*QU[0][1]*q[0]*q[1])/(QU[0][0]*QU[1][1]-2*QU[0][1]*QU[0][1]);
			}
			else if(nUSV==1)
			{
				q[0]=Kernel::k_function(&prob->x[usv[0]],&prob->x[i],*param);
				var_train[i]=Kernel::k_function(&prob->x[i],&prob->x[i],*param)-q[0]*q[0]/QU[0][0];
			}
	}
}*/


static void var_noise_compute_ML(const struct svm_problem *prob, struct svm_model *model, Vector var_train)
{
	longword i;
	word iter,maxiter=100;
	Vector eta = Vector(prob->l);
	for(i=0;i<prob->l;i++)
	{
		double error=svm_predict(model,&prob->x[i])-prob->y[i];
		eta[i]=error*error;
	}
	ofstream out("prova2.xls");
	for(double sigma=0;sigma<=1;sigma+=0.01)
	{
		double Phi=0;
		for(i=0;i<prob->l;i++)
			Phi+=eta[i]/(var_train[i]+sigma)+log(var_train[i]+sigma);
		out << sigma << "	" << Phi/prob->l << endl;
	}
	out.close();
	double alpha=0,v=1,eps_stop=0.0001;
	bool go_on=true;
	for(iter=0;(iter<maxiter)&&go_on;iter++)
	{
		double Phi1=0,Phi2=0;
		for(i=0;i<prob->l;i++)
		{
			double u=var_train[i]+v;
			Phi1+=v*(u-eta[i])/(u*u*prob->l);
			Phi2+=v*v*(2*eta[i]-u)/(u*u*u*prob->l);
		}
		Phi2+=Phi1;
		if(Phi2!=0)
		{
			double delta=Phi1/Phi2;
			go_on=(sqrt(v)*fabs(exp(-0.5*delta)-1)>eps_stop);
			alpha-=delta;
			v=exp(alpha);
		}
		else go_on=false;
	}
	model->noise_var=v;
}


static double ErrorCI(Vector var_train, Vector eta, double std)
{
	longword i;
	const double step=0.01;
	double J=0,p;
	for(p=0.5*step;p<1;p+=step)
	{
		double rate=0;
		for(i=0;i<eta.size();i++)
			if(1-2*Queue(float(fabs(eta[i])/sqrt(var_train[i]+std*std)))<p) rate++;
		rate/=eta.size();
		J+=fabs(p-rate)*step;
	}
	return J;
}


static void var_noise_compute_CI(const struct svm_problem *prob, struct svm_model *model, Vector var_train)
{
	Vector eta = Vector(prob->l);
	double rmse=0;
	for(longword i=0;i<prob->l;i++)
	{
		eta[i]=svm_predict(model,&prob->x[i])-prob->y[i];
		rmse+=eta[i]*eta[i];
	}
	const double tol=2.0e-4,zeps=1.0e-10,cgold=0.3819660;
	unsigned short iter;
	double left=1e-10,right=2*sqrt(rmse/prob->l),middle=sqrt(rmse/prob->l),a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0,rho;
	a=(left < right ? left : right);
	b=(left > right ? left : right);
	x=w=v=middle;
	fw=fv=fx=ErrorCI(var_train,eta,x);
	for(iter=1;iter<=brentiter;iter++)
	{
		xm=0.5*(a+b);
		tol2=2*(tol1=tol*fabs(x)+zeps);
		if(fabs(x-xm)<=(tol2-0.5*(b-a)))
		{
			rho=x;
			model->noise_var=rho*rho;
			return;
		}
		if(fabs(e)>tol1)
		{
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2*(q-r);
			if(q>0)
				p=-p;
			q=fabs(q);
			etemp=e;
			e=d;
			if(fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
				d=cgold*(e=(x>=xm ? a-x : b-x));
			else
			{
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
		else
		{
			d=cgold*(e=(x>=xm ? a-x : b-x));
		}
		u=(fabs(d)>=tol1 ? x+d : x+SIGN(tol1,d));
		fu=ErrorCI(var_train,eta,u);
		if(fu<=fx)
		{
			if(u>=x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		}
		else
		{
			if(u<x)
				a=u;
			else
				b=u;
			if(fu<=fw || w == x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if(fu<=fv || v==x || v==w)
			{
				v=u;
				fv=fu;
			}
		}
	}
	rho=x;
	model->noise_var=rho*rho;
}


/***************************************************************************************************************************************/
struct decision_function
{
	double *alpha;
	double rho;
	double spanbound;
	Matrix L;
	Vector var_train;
};


decision_function svm_train_one(const svm_problem *prob, const svm_parameter *param, double Cp, double Cn)
{
	double *alpha = new double[prob->l],SB=1;
	Solver::SolutionInfo si;
	Vector var_train;
	Matrix L;
	vector<longword> usv;
	switch(param->svm_type)
	{
		case C_SVC:
			solve_c_svc(prob,param,alpha,&si,Cp,Cn);
			if(param->PSB)
				SB=span_bound_compute_svc(prob,param,alpha,&si);
			break;
		case EPSILON_SVR:
			solve_epsilon_svr(prob,param,alpha,&si);
			if(param->PSB)
				SB=span_bound_compute_svr(prob,param,alpha,&si);
			if(param->GGGH)
				var_train_compute(prob,param,alpha,var_train,L,usv);
			break;
	}
	decision_function f;
	f.alpha = alpha;
	f.rho = si.rho;
	f.spanbound = SB;
	if(param->GGGH)
	{
		f.var_train = var_train;
		f.L=L;
	}
	return f;
}


/***************************************************************************************************************************************/
// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
void svm_group_classes(const svm_problem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm)
{
	int l = prob->l;
	int max_nr_class = 16;
	int nr_class = 0;
	int *label = (int *)malloc((max_nr_class)*sizeof(int));
	int *count = (int *)malloc((max_nr_class)*sizeof(int));
	int *data_label = new int[l];	
	int i;
	for(i=0;i<l;i++)
	{
		int this_label = (int)prob->y[i];
		int j;
		for(j=0;j<nr_class;j++)
		{
			if(this_label == label[j])
			{
				++count[j];
				break;
			}
		}
		data_label[i] = j;
		if(j == nr_class)
		{
			if(nr_class == max_nr_class)
			{
				max_nr_class *= 2;
				label = (int *)realloc(label,max_nr_class*sizeof(int));
				count = (int *)realloc(count,max_nr_class*sizeof(int));
			}
			label[nr_class] = this_label;
			count[nr_class] = 1;
			++nr_class;
		}
	}
	int *start = new int[nr_class];
	start[0] = 0;
	for(i=1;i<nr_class;i++)
		start[i] = start[i-1]+count[i-1];
	for(i=0;i<l;i++)
	{
		perm[start[data_label[i]]] = i;
		++start[data_label[i]];
	}
	start[0] = 0;
	for(i=1;i<nr_class;i++)
		start[i] = start[i-1]+count[i-1];

	*nr_class_ret = nr_class;
	*label_ret = label;
	*start_ret = start;
	*count_ret = count;
	delete[] data_label;
}


svm_model *svm_train(const svm_problem *prob, const svm_parameter *param)
{
	svm_model *model = new svm_model[1];
	model->param = *param;
	model->free_sv = 0;
	if(param->svm_type == EPSILON_SVR)
	{
		model->nr_class = 2;
		model->label = NULL;
		model->nSV = NULL;
		model->sv_coef = new double*[1];
		decision_function f = svm_train_one(prob,param,0,0);
		model->rho = new double[1];
		model->rho[0] = f.rho;
		int nSV = 0;
		unsigned long i;
		for(i=0;i<prob->l;i++)
			if(fabs(f.alpha[i]) > 0) ++nSV;
		model->l = nSV;
		model->SV = new svm_node[nSV];
		model->sv_coef[0] = new double[nSV];
		model->SB = f.spanbound;
		vector<longword> usv;
		cout << "#SV = " << model->l << endl << "SB = " << model->SB << endl;
		int j = 0;
		for(i=0;i<prob->l;i++)
			if(fabs(f.alpha[i]) > 0)
			{
				model->SV[j] = prob->x[i];
				model->sv_coef[0][j] = f.alpha[i];
				if(fabs(f.alpha[i])<param->C)
					usv.push_back(j);
				++j;
			}
		delete[] f.alpha;
		if(param->GGGH)
		{
			model->usv=vector<longword>(usv);
			model->L=Matrix(f.L);
			var_noise_compute_CI(prob,model,f.var_train);
		}
	}
	else
	{
		// classification
		int l = prob->l;
		int nr_class;
		int *label = NULL;
		int *start = NULL;
		int *count = NULL;
		int *perm = new int[l];
		svm_group_classes(prob,&nr_class,&label,&start,&count,perm);		
		svm_node *x = new svm_node[l];
		int i;
		for(i=0;i<l;i++)
			x[i] = prob->x[perm[i]];
		// train k*(k-1)/2 models
		bool *nonzero = new bool[l];
		for(i=0;i<l;i++)
			nonzero[i] = false;
		decision_function *f = new decision_function[nr_class*(nr_class-1)/2];
		model->SB=0;
		int p = 0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				svm_problem sub_prob;
				int si = start[i], sj = start[j];
				int ci = count[i], cj = count[j];
				sub_prob.l = ci+cj;
				sub_prob.x = new svm_node[sub_prob.l];
				sub_prob.y = new double[sub_prob.l];
				int k;
				for(k=0;k<ci;k++)
				{
					sub_prob.x[k] = x[si+k];
					sub_prob.y[k] = +1;
				}
				for(k=0;k<cj;k++)
				{
					sub_prob.x[ci+k] = x[sj+k];
					sub_prob.y[ci+k] = -1;
				}
				f[p] = svm_train_one(&sub_prob,param,param->C,param->C);
				model->SB+=f[p].spanbound*ci*cj/(prob->l*prob->l);
				for(k=0;k<ci;k++)
					if(!nonzero[si+k] && fabs(f[p].alpha[k]) > 0)
						nonzero[si+k] = true;
				for(k=0;k<cj;k++)
					if(!nonzero[sj+k] && fabs(f[p].alpha[ci+k]) > 0)
						nonzero[sj+k] = true;
				delete[] sub_prob.x,sub_prob.y;
				++p;
			}
		// build output
		model->nr_class = nr_class;
		model->label = new int[nr_class];
		for(i=0;i<nr_class;i++)
			model->label[i] = label[i];
		model->rho = new double[nr_class*(nr_class-1)/2];
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			model->rho[i] = f[i].rho;
		int total_sv = 0;
		int *nz_count = new int[nr_class];
		model->nSV = new int[nr_class];
		for(i=0;i<nr_class;i++)
		{
			int nSV = 0;
			for(int j=0;j<count[i];j++)
				if(nonzero[start[i]+j])
				{	
					++nSV;
					++total_sv;
				}
			model->nSV[i] = nSV;
			nz_count[i] = nSV;
			cout << "#SV per class " << i+1 << " = " << nSV << endl;
		}
		if(nr_class>2)
			cout << "Total #SV = " << total_sv << endl;
		if(param->PSB)
			cout << "SB = " << model->SB << endl;
		flush(cout);

		/*ofstream outsb("SB.xls",ios::app);
		outsb << model->SB << "	";
		outsb.close();*/

		model->l = total_sv;
		model->SV = new svm_node[total_sv];
		p = 0;
		for(i=0;i<l;i++)
			if(nonzero[i]) model->SV[p++] = x[i];
		int *nz_start = new int[nr_class];
		nz_start[0] = 0;
		for(i=1;i<nr_class;i++)
			nz_start[i] = nz_start[i-1]+nz_count[i-1];
		model->sv_coef = new double*[nr_class-1];
		for(i=0;i<nr_class-1;i++)
			model->sv_coef[i] = new double[total_sv];
		p = 0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				// classifier (i,j): coefficients with
				// i are in sv_coef[j-1][nz_start[i]...],
				// j are in sv_coef[i][nz_start[j]...]
				int si = start[i];
				int sj = start[j];
				int ci = count[i];
				int cj = count[j];
				int q = nz_start[i];
				int k;
				for(k=0;k<ci;k++)
					if(nonzero[si+k])
						model->sv_coef[j-1][q++] = f[p].alpha[k];
				q = nz_start[j];
				for(k=0;k<cj;k++)
					if(nonzero[sj+k])
						model->sv_coef[i][q++] = f[p].alpha[ci+k];
				++p;
			}
		free(label);
		free(count);
		delete[] perm,start,x,nonzero;
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			delete[] f[i].alpha;
		delete[] f,nz_count,nz_start;
	}
	return model;
}


/***************************************************************************************************************************************/
double svm_cross_validation(const svm_problem *prob, const svm_parameter *param, int nr_fold)		// stratified CV
{
	int i;
	int *fold_start = new int[nr_fold+1];
	int l = prob->l;
	int *perm = new int[l];
	int nr_class;
	// stratified cv may not give leave-one-out rate
	// Each class to l folds -> some folds may have zero elements
	if((param->svm_type == C_SVC)&&(nr_fold < l))
	{
		int *start = NULL;
		int *label = NULL;
		int *count = NULL;
		svm_group_classes(prob,&nr_class,&label,&start,&count,perm);
		int *fold_count = new int[nr_fold];
		int c;
		int *index = new int[l];
		for(i=0;i<l;i++)
			index[i]=perm[i];
		for (c=0; c<nr_class; c++) 
			for(i=0;i<count[c];i++)
			{
				int j = i+rand()%(count[c]-i);
				SWAP(index[start[c]+j],index[start[c]+i]);
			}
		for(i=0;i<nr_fold;i++)
		{
			fold_count[i] = 0;
			for (c=0; c<nr_class;c++)
				fold_count[i]+=(i+1)*count[c]/nr_fold-i*count[c]/nr_fold;
		}
		fold_start[0]=0;
		for (i=1;i<=nr_fold;i++)
			fold_start[i] = fold_start[i-1]+fold_count[i-1];
		for (c=0; c<nr_class;c++)
			for(i=0;i<nr_fold;i++)
			{
				int begin = start[c]+i*count[c]/nr_fold;
				int end = start[c]+(i+1)*count[c]/nr_fold;
				for(int j=begin;j<end;j++)
				{
					perm[fold_start[i]] = index[j];
					fold_start[i]++;
				}
			}
		fold_start[0]=0;
		for (i=1;i<=nr_fold;i++)
			fold_start[i] = fold_start[i-1]+fold_count[i-1];
		free(label);
		free(count);	
		delete[] start,index,fold_count;
	}
	else
	{
		for(i=0;i<l;i++) perm[i]=i;
		for(i=0;i<l;i++)
		{
			int j = i+rand()%(l-i);
			SWAP(perm[i],perm[j]);
		}
		for(i=0;i<=nr_fold;i++)
			fold_start[i]=i*l/nr_fold;
	}
	double *target = new double[prob->l],error=0;
	for(i=0;i<nr_fold;i++)
	{
		int begin = fold_start[i];
		int end = fold_start[i+1];
		int j,k;
		struct svm_problem subprob;
		subprob.l = l-(end-begin);
		subprob.x = new struct svm_node[subprob.l];
		subprob.y = new double[subprob.l];
		k=0;
		for(j=0;j<begin;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			++k;
		}
		for(j=end;j<l;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			++k;
		}
		struct svm_model *submodel = svm_train(&subprob,param);
		for(j=begin;j<end;j++)
			target[perm[j]] = svm_predict(submodel,prob->x+perm[j]);
		svm_destroy_model(submodel);
		delete[] subprob.x,subprob.y;
	}
	delete[] fold_start,perm;	
	for(i=0;i<long(prob->l);i++)
		if((param->svm_type==C_SVC)&&(prob->y[i]!=target[i]))			// C_SVC: return CV error rate
			error++;
		else if(param->svm_type==EPSILON_SVR)							// EPSILON_SVR: return CV MAE
			error+=fabs(prob->y[i]-target[i]);
	cout << "CV error = " << error/prob->l << endl;
	return error/prob->l;
}


double svm_hold_out_valid(const svm_problem *prob, const svm_problem *valprob, const svm_parameter *param)		// stratified CV
{
	unsigned long i;
	struct svm_model *model = svm_train(prob,param);
	double error=0;
	for(i=0;i<valprob->l;i++)
	{
		double target=svm_predict(model,&valprob->x[i]);
		if((param->svm_type==C_SVC)&&(valprob->y[i]!=target))			// C_SVC: return hold-out error rate
			error++;
		else if(param->svm_type==EPSILON_SVR)							// EPSILON_SVR: return hold-out MAE
			error+=fabs(valprob->y[i]-target);
	}
	svm_destroy_model(model);
	cout << "hold-out error on validation set = " << error/prob->l << endl;
	return error/prob->l;
}


/***************************************************************************************************************************************/
double ErrorSB(struct svm_problem *prob, struct svm_parameter *param, Vector &theta)
{
	if(param->svm_type==C_SVC)
	{
		if((param->svm_type==C_SVC)&&((fabs(theta[0])>10)||(theta[1]<-9)||(theta[1]>5)))			// outside allowed search region for powell
			return 1;
		else
		{
			param->C=exp(theta[0]);										// theta[0] = ln(C)
			param->gamma=0.5*exp(-theta[1]);							// theta[1] = ln(sigmasquare)
			if(param->svm_type==EPSILON_SVR) param->p=exp(theta[2]);
			struct svm_model *model = svm_train(prob,param);
			double SB = (model->nSV>0 ? model->SB : 1);
			svm_destroy_model(model);
			return SB;
		}
	}
	else if(powell_svr_tanh)
	{
		param->C=param_min[0]+0.5*(tanh(theta[0])+1)*(param_max[0]-param_min[0]);
		param->gamma=0.5/(param_min[1]+0.5*(tanh(theta[1])+1)*(param_max[1]-param_min[1]));
		param->p=param_min[2]+0.5*(tanh(theta[2])+1)*(param_max[2]-param_min[2]);
		struct svm_model *model = svm_train(prob,param);
		double SB = model->SB;
		svm_destroy_model(model);
		return SB;
	}
	else
	{
		if((param->svm_type==EPSILON_SVR)&&((theta[0]<theta_min[0])||(theta[0]>theta_max[0])||(theta[1]<theta_min[1])||(theta[1]>theta_max[1])||(theta[2]<theta_min[2])||(theta[2]>theta_max[2])))
			return (prob->l+1)*exp(theta_max[0]+1.0)+2*exp(theta_max[2]+1.0)+1;
		else
		{
			param->C=exp(theta[0]);										// theta[0] = ln(C)
			param->gamma=0.5*exp(-theta[1]);							// theta[1] = ln(sigmasquare)
			param->p=exp(theta[2]);										// theta[2] = ln(epsilon)
			struct svm_model *model = svm_train(prob,param);
			double SB = model->SB;
			svm_destroy_model(model);
			return SB;
		}
	}
}


double ErrorLine(struct svm_problem *prob, struct svm_parameter *param, double rho, Vector &theta, Vector &searchdir)
{
	unsigned short j;
	Vector thetamove = Vector(theta);
	for(j=0;j<theta.size();j++)
		thetamove[j]+=rho*searchdir[j];
	return ErrorSB(prob,param,thetamove);
}


void MinBracket(struct svm_problem *prob, struct svm_parameter *param, double &left, double &middle, double &right, Vector &theta, Vector &searchdir)
{
	const double gold=1.618034,glimit=100.0,tiny=1.0e-20;
	double parablim,parab,r,q,Jparab,dummy,Jleft=ErrorLine(prob,param,left,theta,searchdir),Jmiddle=ErrorLine(prob,param,middle,theta,searchdir);
	if(Jmiddle>Jleft)
	{
		SHFT(dummy,left,middle,dummy)
		SHFT(dummy,Jmiddle,Jleft,dummy)
	}
	right=middle+gold*(middle-left);
	double Jright=ErrorLine(prob,param,right,theta,searchdir);
	while(Jmiddle>Jright)
	{
		r=(middle-left)*(Jmiddle-Jright);
		q=(middle-right)*(Jmiddle-Jleft);
		parab=middle-((middle-right)*q-(middle-left)*r)/(2.0*SIGN(FMAX(fabs(q-r),tiny),q-r));
		parablim=middle+glimit*(right-middle);
		if((middle-parab)*(parab-right)>0)
		{
			Jparab=ErrorLine(prob,param,parab,theta,searchdir);
			if(Jparab<Jright)
			{
				left=middle;
				middle=parab;
				Jleft=Jmiddle;
				Jmiddle=Jparab;
				return;
			}
			else if(Jparab>Jmiddle)
			{
				right=parab;
				Jright=Jparab;
				return;
			}
			parab=right+gold*(right-middle);
			Jparab=ErrorLine(prob,param,parab,theta,searchdir);
		}
		else if((right-parab)*(parab-parablim)>0)
		{
			Jparab=ErrorLine(prob,param,parab,theta,searchdir);
			if (Jparab < Jright)
			{
				SHFT(middle,right,parab,right+gold*(right-middle))
				SHFT(Jmiddle,Jright,Jparab,ErrorLine(prob,param,parab,theta,searchdir))
			}
		}
		else if((parab-parablim)*(parablim-right)>=0)
		{
			parab=parablim;
			Jparab=ErrorLine(prob,param,parab,theta,searchdir);
		}
		else
		{
			parab=right+gold*(right-middle);
			Jparab=ErrorLine(prob,param,parab,theta,searchdir);
		}
		SHFT(left,middle,right,parab)
		SHFT(Jleft,Jmiddle,Jright,Jparab)
	}
}


double Brent(struct svm_problem *prob, struct svm_parameter *param, double left, double middle, double right, double &rho, Vector &theta, Vector &searchdir)
{
	const double tol=2.0e-4,zeps=1.0e-10,cgold=0.3819660;
	unsigned short iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0;
	a=(left < right ? left : right);
	b=(left > right ? left : right);
	x=w=v=middle;
	fw=fv=fx=ErrorLine(prob,param,x,theta,searchdir);
	for(iter=1;iter<=brentiter;iter++)
	{
		xm=0.5*(a+b);
		tol2=2*(tol1=tol*fabs(x)+zeps);
		if(fabs(x-xm)<=(tol2-0.5*(b-a)))
		{
			rho=x;
			return fx;
		}
		if(fabs(e)>tol1)
		{
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2*(q-r);
			if(q>0)
				p=-p;
			q=fabs(q);
			etemp=e;
			e=d;
			if(fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
				d=cgold*(e=(x>=xm ? a-x : b-x));
			else
			{
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
		else
		{
			d=cgold*(e=(x>=xm ? a-x : b-x));
		}
		u=(fabs(d)>=tol1 ? x+d : x+SIGN(tol1,d));
		fu=ErrorLine(prob,param,u,theta,searchdir);
		if(fu<=fx)
		{
			if(u>=x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		}
		else
		{
			if(u<x)
				a=u;
			else
				b=u;
			if(fu<=fw || w == x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if(fu<=fv || v==x || v==w)
			{
				v=u;
				fv=fu;
			}
		}
	}
	rho=x;
	return fx;
}


double LineSearch(struct svm_problem *prob, struct svm_parameter *param, Vector &theta, Vector &searchdir)
{
	unsigned short j;
	double left=0,middle=1,right,rho;
	Vector Theta = Vector(theta), SearchDir = Vector(searchdir);
	MinBracket(prob,param,left,middle,right,Theta,SearchDir);
	double J=Brent(prob,param,left,middle,right,rho,Theta,SearchDir);
	for(j=0;j<theta.size();j++)
	{
		searchdir[j]*=rho;
		theta[j]+=searchdir[j];
	}
	return J;
}


void Powell(struct svm_problem *prob, struct svm_parameter *param, Vector &theta)
{
	const double tiny=1.0e-25,eps_stop=2*DBL_MIN;
	unsigned short i,I,j,iter;
	Matrix pseudoconjug = Matrix(theta.size(),theta.size(),0);
	for(i=0;i<theta.size();i++) pseudoconjug[i][i]=1;
	double delta,J=DBL_MAX,Jold,Jextra,t;
	Vector thetaold = Vector(theta.size(),0), thetaextra = Vector(theta.size()), searchdir = Vector(theta.size());
	theta[0]=theta[1]=0;
	if(param->svm_type==C_SVC)
		theta[0]=theta[1]=0;
	else if(powell_svr_tanh)
		for(i=0;i<3;i++)
			theta[i]=0.5*log(param_init[i]-param_min[i])-0.5*log(param_max[i]-param_init[i]);
	else
	{
		theta[0]=theta[1]=0;
		theta[2]=-5;
	}
	for(iter=1;iter<powelliter;iter++)
	{
		Jold=J;
		I=0;
		delta=0;
		for(i=0;i<theta.size();i++)
		{
			for(j=0;j<theta.size();j++)
				searchdir[j]=pseudoconjug[j][i];
			Jextra=J;
			J=LineSearch(prob,param,theta,searchdir);
			if(Jextra-J>delta)
			{
				delta=Jextra-J;
				I=i;
			}
		}
		if(2*(Jold-J)<=eps_stop*(fabs(Jold)+fabs(J))+tiny)
			return;
		for(j=0;j<theta.size();j++)
		{
			thetaextra[j]=2*theta[j]-thetaold[j];
			searchdir[j]=theta[j]-thetaold[j];
			thetaold[j]=theta[j];
		}
		Jextra=ErrorSB(prob,param,thetaextra);
		if(Jextra < Jold)
		{
			t=2*(Jold-2*J+Jextra)*(Jold-J-delta)*(Jold-J-delta)-delta*(Jold-Jextra)*(Jold-Jextra);
			if (t<0)
			{
				J=LineSearch(prob,param,theta,searchdir);
				for(j=0;j<theta.size();j++)
				{
					pseudoconjug[j][I]=pseudoconjug[j][1];
					pseudoconjug[j][1]=searchdir[j];
				}
			}
		}
	}
}


void Powell(struct svm_problem *prob, struct svm_parameter *param)
{
	Vector theta = Vector((param->svm_type==C_SVC ? 2 : 3));
	Powell(prob,param,theta);
	param->C=exp(theta[0]);
	param->gamma=0.5*exp(-theta[1]);
	if(param->svm_type==C_SVC)
	{
		param->C=exp(theta[0]);
		param->gamma=0.5*exp(-theta[1]);
	}
	else if(powell_svr_tanh)
	{
		param->C=param_min[0]+0.5*(tanh(theta[0])+1)*(param_max[0]-param_min[0]);
		param->gamma=0.5/(param_min[1]+0.5*(tanh(theta[1])+1)*(param_max[1]-param_min[1]));
		param->p=param_min[2]+0.5*(tanh(theta[2])+1)*(param_max[2]-param_min[2]);
	}
	else
	{
		param->C=exp(theta[0]);
		param->gamma=0.5*exp(-theta[1]);
		param->p=exp(theta[2]);
	}
}


/***************************************************************************************************************************************/
vector<bool> PDDP(Matrix X, longword target)
{
	// target = n. cluster voluti
	// X.rows() = n. campioni da clusterizzare
	// X.cols() = n. feature
	// X[i][j] = feature j-ma del campione i-mo
	// ritorna vettore di booleani (tanti quanti i campioni) con true se un campione è un centro-cluster e false altrimenti,
	// scegliendo in ciascun cluster il campione più vicino al baricentro (mean[]...)
	// opzioni: così o ritornare direttamente il baricentro?
	longword row=X.rows(),k,valid;
	if(row<=target)
		errorlog("PDDP. Fewer target clusters than samples.");
	word col=X.cols(),i,j;
	Matrix mean = Matrix(target,col,0);
	Tensor square = Tensor(target,col,col,0);
	vector<longword> size = vector<longword>(target,0);
	vector<longword> label = vector<longword>(row,0);
	size[0]=row;
	for(k=0;k<row;k++)
		for(i=0;i<col;i++)
		{
			mean[0][i]+=X[k][i]/size[0];
			square[0][i][i]+=X[k][i]*X[k][i]/size[0];
			for(j=i+1;j<col;j++)
				square[0][i][j]=X[k][i]*X[k][j]/size[0];
		}
	for(valid=1;valid<target;valid++)
	{
		longword sizemax=1,split=0;
		for(k=0;k<valid;k++)
			if(size[k]>sizemax)
			{
				sizemax=size[k];
				split=k;
			}
		Matrix cov = Matrix(square[split]);
		for(i=0;i<col;i++)
		{
			cov[i][i]-=mean[split][i]*mean[split][i];
			for(j=i+1;j<col;j++)
			{
				cov[i][j]-=mean[split][i]*mean[split][j];
				cov[j][i]=cov[i][j];
			}
		}
		Vector eigen=EigenPower(cov);
		bool singular=(eigen[0]==2);
		vector<bool> projmove = vector<bool>(row,false);
		longword check=0;
		if(!singular)
		{
			for(k=0;k<row;k++)
				if(label[k]==split)
				{
					double proj=0;
					for(i=0;i<col;i++)
						proj+=eigen[i]*(X[k][i]-mean[split][i]);
					if(proj>0)
					{
						projmove[k]=true;
						check++;
					}
				}
			singular = (check==0)||(check==size[split]);
		}
		for(k=0;k<row;k++)
			if((label[k]==split)&&((!singular&&projmove[k])||(singular&&(size[valid]<size[split]))))
			{
				label[k]=valid;
				for(i=0;i<col;i++)
				{
					mean[split][i]=(mean[split][i]*size[split]-X[k][i])/(size[split]-1);
					mean[valid][i]=(mean[valid][i]*size[valid]+X[k][i])/(size[valid]+1);
					square[split][i][i]=(square[split][i][i]*size[split]-X[k][i]*X[k][i])/(size[split]-1);
					square[valid][i][i]=(square[valid][i][i]*size[valid]+X[k][i]*X[k][i])/(size[valid]+1);
					for(j=i+1;j<col;j++)
					{
						square[split][i][j]=(square[split][i][j]*size[split]-X[k][i]*X[k][j])/(size[split]-1);
						square[valid][i][j]=(square[valid][i][j]*size[valid]+X[k][i]*X[k][j])/(size[valid]+1);
					}
				}
				size[split]--;
				size[valid]++;
			}
	}
	vector<bool> selected = vector<bool>(row,false);
	for(valid=0;valid<target;valid++)
	{
		longword closest=0;
		double mindist=DBL_MAX;
		for(k=0;k<row;k++)
			if(label[k]==valid)
			{
				double dist=0;
				for(i=0;i<col;i++)
					dist+=(mean[valid][i]-X[k][i])*(mean[valid][i]-X[k][i]);
				if(dist<mindist)
				{
					mindist=dist;
					closest=k;
				}
			}
		selected[closest]=true;
	}
	return selected;
}


/***************************************************************************************************************************************/
void svm_predict_values(const svm_model *model, const svm_node *x, double* dec_values)
{
	if(model->param.svm_type == EPSILON_SVR)
	{
		double *sv_coef = model->sv_coef[0];
		double sum = 0;
		
		for(unsigned long i=0;i<model->l;i++)
			sum += sv_coef[i] * Kernel::k_function(x,model->SV+i,model->param);
		sum -= model->rho[0];
		*dec_values = sum;
	}
	else
	{
		int i;
		int nr_class = model->nr_class;
		int l = model->l;
//		double *kvalue = new double[l];
		Vector kvalue = Vector(l);
		for(i=0;i<l;i++)
			kvalue[i] = Kernel::k_function(x,model->SV+i,model->param);
//		int *start = new int[nr_class];
		vector<int> start = vector<int>(nr_class);
		start[0] = 0;
		for(i=1;i<nr_class;i++)
			start[i] = start[i-1]+model->nSV[i-1];
		int p=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				double sum = 0;
				int si = start[i];
				int sj = start[j];
				int ci = model->nSV[i];
				int cj = model->nSV[j];
				int k;
				double *coef1 = model->sv_coef[j-1];
				double *coef2 = model->sv_coef[i];
				for(k=0;k<ci;k++)
					if(kvalue[si+k]!=0)
						sum += coef1[si+k] * kvalue[si+k];
				for(k=0;k<cj;k++)
					if(kvalue[sj+k]!=0)
						sum += coef2[sj+k] * kvalue[sj+k];
				sum -= model->rho[p];
				dec_values[p] = sum;
				p++;
			}
//		delete[] kvalue,start;
	}
}


double svm_predict(const svm_model *model, const svm_node *x)
{
	if(model->param.svm_type == EPSILON_SVR)
	{
		double res;
		svm_predict_values(model, x, &res);
		return res;
	}
	else
	{
		int i;
		int nr_class = model->nr_class;
		double *dec_values = new double[nr_class*(nr_class-1)/2];
		svm_predict_values(model, x, dec_values);
		vector<int> vote = vector<int>(nr_class,0);
//		int *vote = new int[nr_class];
//		for(i=0;i<nr_class;i++)
//			vote[i] = 0;
		int pos=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				if(dec_values[pos++] > 0)
					++vote[i];
				else
					++vote[j];
			}
		int vote_max_idx = 0;
		for(i=1;i<nr_class;i++)
			if(vote[i] > vote[vote_max_idx])
				vote_max_idx = i;
//		delete[] vote,dec_values;
		delete[] dec_values;
		return model->label[vote_max_idx];
	}
}


double svm_predict(const svm_model *model, const svm_node *x, double sigma[3])
{
	if(model->param.GGGH)
	{
		word i,j;
		/*double qstar=Kernel::k_function(x,x,model->param),approxvar;
		Vector q = Vector(model->L.rows());
		Matrix QQ = Matrix(model->L.rows(),model->L.rows());*/
		double approxvar=Kernel::k_function(x,x,model->param);
		if(model->L.rows()>2)
		{
			for(i=0;i<model->L.rows();i++)
			{
				double LQ=0;
				for(j=0;j<i+1;j++)
					LQ+=model->L[i][j]*Kernel::k_function(&model->SV[model->usv[j]],x,model->param);
				approxvar-=LQ*LQ;
			}
			/*for(i=0;i<model->L.rows();i++)
			{
				q[i]=Kernel::k_function(&model->SV[model->usv[i]],x,model->param);
				QQ[i][i]=model->L[i][i]-q[i]*q[i]/qstar;
			}
			for(i=0;i<model->L.rows();i++)
				for(j=i+1;j<model->L.rows();j++)
				{
					QQ[i][j]=model->L[i][j]-q[i]*q[j]/qstar;
					QQ[j][i]=QQ[i][j];
				}
			approxvar=1/(1/qstar-2*ConjugateGrad(QQ,q)/(qstar*qstar));*/
		}
		else if(model->L.rows()==2)
		{
			double q[2],QU[2][2];
			for(i=0;i<model->L.rows();i++)
			{
				q[i]=Kernel::k_function(&model->SV[model->usv[i]],x,model->param);
				QU[i][i]=Kernel::k_function(&model->SV[model->usv[i]],&model->SV[model->usv[i]],model->param);
			}
			QU[0][1]=QU[1][0]=Kernel::k_function(&model->SV[model->usv[0]],&model->SV[model->usv[1]],model->param);
			approxvar-=(QU[1][1]*q[0]*q[0]+QU[0][0]*q[1]*q[1]-2*QU[0][1]*q[0]*q[1])/(QU[0][0]*QU[1][1]-2*QU[0][1]*QU[0][1]);
		}
		else if(model->L.rows()==1)
		{
			double q=Kernel::k_function(&model->SV[model->usv[0]],x,model->param);
			approxvar-=q*q/Kernel::k_function(&model->SV[model->usv[0]],&model->SV[model->usv[0]],model->param);
		}
		approxvar=(approxvar<0 ? 0 :approxvar);
		sigma[0]=sqrt(approxvar);
		sigma[1]=sqrt(model->noise_var);
		sigma[2]=sqrt(approxvar+model->noise_var);
	}
	return svm_predict(model,x);
}


void svm_regression_bias_training(struct svm_problem *prob, struct svm_model *model)
{
	if(model->param.svm_type!=EPSILON_SVR)
		exit(0);
	longword i;
	double bias=0;
	for(i=0;i<prob->l;i++)
		bias+=svm_predict(model,&prob->x[i])-prob->y[i];
	model->rho[0]+=bias/prob->l;
}


/***************************************************************************************************************************************/
void init_default_param(svm_parameter *param, int type)
{
	param->svm_type = type;
	param->GGGH = false;
	param->degree = 3;			// unused
	param->coef0 = 0;			// unused
	//param->p = 0.01;			// unused
	param->cache_size = 100;	// default = 100
	param->eps = 1e-3;			// default = 1e-3
	param->shrinking = 1;		// default = 1
}


void svm_destroy_model(svm_model* model)
{
	if(model->free_sv && model->l > 0)
	for (unsigned long i = 0; i < model->l; i++)
		delete[] model->SV[i].values;
	for(int i=0;i<model->nr_class-1;i++)
		delete[] model->sv_coef[i];
	delete[] model->SV,model->sv_coef,model->rho,model->label,model->nSV,model;
}


void svm_destroy_prob(svm_problem *prob)
{
	unsigned long i;
	for(i=0;i<prob->l;i++)
		delete[] prob->x[i].values;
	delete[] prob->x,prob->y;
}


const char *svm_check_parameter(const svm_problem *prob, const svm_parameter *param)
{
	// svm_type
	int svm_type = param->svm_type;
	if(svm_type != C_SVC &&
	   svm_type != EPSILON_SVR)
		return "unknown svm type";
	// kernel_type, degree
	int kernel_type = param->kernel_type;
	if(kernel_type != LINEAR &&
	   kernel_type != POLY &&
	   kernel_type != RBF &&
	   kernel_type != SIGMOID)
		return "unknown kernel type";
	if(param->degree < 0)
		return "degree of polynomial kernel < 0";
	// cache_size,eps,C,p,shrinking
	if(param->cache_size <= 0)
		return "cache_size <= 0";
	if(param->eps <= 0)
		return "eps <= 0";
	if(svm_type == C_SVC ||
	   svm_type == EPSILON_SVR)
		if(param->C <= 0)
			return "C <= 0";
	if(svm_type == EPSILON_SVR)
		if(param->p < 0)
			return "p < 0";
	if(param->shrinking != 0 &&
	   param->shrinking != 1)
		return "shrinking != 0 and shrinking != 1";
	return NULL;
}


/***************************************************************************************************************************************/
void svm_regression_read_ground(struct svm_problem *prob, Vector gain, Vector offset, char name[80])
{
	ifstream in(name);
	longword i=0;
	word r;
	double dummy;
	while(!in.eof())
	{
		in >> dummy;
		i++;
	}
	in.close();
	prob->l=i/gain.size();
	prob->y = new double[prob->l];
	prob->x = new svm_node[prob->l];
	ifstream in1(name);
	for(i=0;i<prob->l;i++)
	{
		prob->x[i].dim=gain.size()-1;
		prob->x[i].values = new double[prob->x[i].dim];
		for(r=0;r<prob->x[i].dim;r++)
		{
			in1 >> prob->x[i].values[r];
			prob->x[i].values[r]=gain[r]*prob->x[i].values[r]+offset[r];
		}
		in1 >> prob->y[i];
		prob->y[i]=gain[gain.size()-1]*prob->y[i]+offset[gain.size()-1];
	}
	in1.close();
}


/*void svm_regression_write_test(struct svm_model *model, char testname[80], char outname[80])
{
	longword i,j,nUSV=model->L.rows();
	ofstream out(outname);
	ifstream in(testname);
	svm_node *test = new svm_node;
	test->dim=model->SV[0].dim;
	test->values = new double[test->dim];
	double target,estimate;
	while(!in.eof())
	{
		for(j=0;j<test->dim;j++)
		{
			in >> test->values[j];
			out << test->values[j] << "	";
		}
		in >> target;
		estimate=svm_predict(model,test);
		if(model->param.GGGH)
		{
			double qstar=Kernel::k_function(test,test,model->param);
			Vector q = Vector(nUSV);
			Matrix QQ = Matrix(nUSV,nUSV);
			for(i=0;i<nUSV;i++)
			{
				q[i]=Kernel::k_function(&model->SV[model->usv[i]],test,model->param);
				QQ[i][i]=model->L[i][i]-q[i]*q[i]/qstar;
			}
			for(i=0;i<nUSV;i++)
				for(j=i+1;j<nUSV;j++)
				{
					QQ[i][j]=model->L[i][j]-q[i]*q[j]/qstar;
					QQ[j][i]=QQ[i][j];
				}
			double errorvar=1/(1/qstar-2*ConjugateGrad(QQ,q)/(qstar*qstar));
			//double errorvar=1+2*ConjugateGrad(model->L,q);
			out << target << "	" << estimate << "	" << errorvar+model->noise_var << endl;
		}
		else
			out << target << "	" << estimate << endl;
	}
	delete[] test->values;
	delete test;
	in.close();
	out.close();
}*/


void svm_regression_write_test(struct svm_model *model, char testname[80], char outname[80])
{
	longword i,j;
	word r;
	ofstream out(outname);
	ifstream in(testname);
	svm_node *test = new svm_node;
	test->dim=model->SV[0].dim;
	test->values = new double[test->dim];
	double target,estimate;
	while(!in.eof())
	{
		for(r=0;r<test->dim;r++)
		{
			in >> test->values[r];
			out << test->values[r] << "	";
		}
		in >> target;
		estimate=svm_predict(model,test);
		if(model->param.GGGH)
		{
			double errorvar=1;
			for(i=0;i<model->L.rows();i++)
			{
				double LQ=0;
				for(j=0;j<i+1;j++)
					LQ+=model->L[i][j]*Kernel::k_function(&model->SV[model->usv[j]],test,model->param);
				errorvar-=LQ*LQ;
			}
			out << target << "	" << estimate << "	" << sqrt(fabs(errorvar)+model->noise_var) << endl;
		}
		else
			out << target << "	" << estimate << endl;
	}
	delete[] test->values;
	delete test;
	in.close();
	out.close();
}


Vector svm_regression_holdout(struct svm_model *model, struct svm_problem *prob, double gain, double offset, char outname[80])
{
	longword i,j,testsize=(model->param.GGGH ? 0 : prob->l);
	short r;
	double rmse=0,mae=0,bias=0,normbias=0,normstd=0,estimate,sigma[3],rate[5]={0,0,0,0,0};
	ofstream out(outname,ios::app);
	out << "x";
	for(r=0;r<prob->x[0].dim;r++) out << ";";
	out << "y;AST;estimate;";
	if(model->param.GGGH)
		out << 	"sigmas;sigman;sigma;norm;area" << endl;
	else
		out << endl;
	for(i=0;i<prob->l;i++)
	{
		if(model->param.GGGH)
		{
			estimate=svm_predict(model,&prob->x[i],sigma);
			for(j=0;j<prob->x[i].dim;j++)
				out << prob->x[i].values[j] << ";";
			double norm=(estimate-prob->y[i])/sigma[2],area=erf(fabs(norm)*sqrt(0.5));
			for(r=4;r>(10*area-5>0 ? 10*area-5 : -1);r--)
				rate[r]++;
			normbias+=norm;
			normstd+=norm*norm;
			if(gain*sigma[2]<SigmaThr)
			{
				out << prob->y[i] << ";"  << gain*prob->y[i]+offset << ";" << gain*estimate+offset << ";" << gain*sigma[0] << ";" << gain*sigma[1] << ";" << gain*sigma[2] << ";" << norm << ";" << area << endl;
				mae+=fabs(prob->y[i]-estimate);
				rmse+=(prob->y[i]-estimate)*(prob->y[i]-estimate);
				bias+=estimate-prob->y[i];
				testsize++;
			}
		}
		else
		{
			estimate=svm_predict(model,&prob->x[i]);
			out << prob->y[i] << ";" << gain*prob->y[i]+offset << ";" << gain*estimate+offset << endl;
			mae+=fabs(prob->y[i]-estimate);
			rmse+=(prob->y[i]-estimate)*(prob->y[i]-estimate);
			bias+=estimate-prob->y[i];
		}
	}
	out << endl << "C =;" << model->param.C << endl << "sigma =;" << sqrt(0.5/model->param.gamma) << endl << "epsilon =;" << model->param.p 
		<< endl << "MAE;" << gain*mae/testsize << endl << "RMSE;" << gain*sqrt(rmse/testsize) << endl << "bias =;" << gain*bias/testsize
		<< endl << "#SV =;" << model->l << endl;
	if(model->param.GGGH)
	{
		out << "#USV =;" << model->L.rows() << endl << "normbias =;" << normbias/testsize
			<< endl << "normstd =;" << sqrt(normstd/testsize-normbias*normbias/(testsize*testsize)) << endl << endl << "USV Cholesky matrix:" << endl;
		for(i=0;i<model->L.rows();i++)
		{
			for(j=0;j<model->L.cols();j++)
				out << model->L[i][j] << ";";
			out << endl;
		}
		out << endl << "p;rate" << endl;
		for(r=0;r<5;r++)
			out << 0.5+0.1*r << ";" << rate[r]/testsize << endl;
	}
	out.close();
	Vector measures;
	measures.push_back(gain*mae/prob->l);
	measures.push_back(gain*sqrt(rmse/prob->l));
	measures.push_back(gain*bias/prob->l);
	if(model->param.GGGH) for(r=0;r<5;r++) measures.push_back(rate[r]/testsize);
	return measures;
}


Matrix PCA(struct svm_problem *prob)
{
	word d=prob->x[0].dim,i,j,k;
	Vector mean = Vector(d,0), buffer = Vector(d);
	Matrix cov = Matrix(d,d,0), proj = Matrix(d,d,0);
	for(i=0;i<prob->l;i++)
		for(j=0;j<d;j++)
		{
			mean[j]+=prob->x[i].values[j]/prob->l;
			for(k=j;k<d;k++)
				cov[j][k]+=prob->x[i].values[j]*prob->x[i].values[k]/prob->l;
		}
	for(i=0;i<d;i++)
	{
		cov[i][i]-=mean[i]*mean[i];
		for(j=i+1;j<d;j++)
		{
			cov[i][j]-=mean[i]*mean[j];
			cov[j][i]=cov[i][j];
		}
		proj[i][i]=1;
	}
	Matrix diag = Matrix(cov);
	double od;
	do
	{
		od=0;
		word p,q;
		double maxod=DBL_MIN;
		for(i=0;i<d-1;i++)
			for(j=i+1;j<d;j++)
			{
				od+=2*diag[i][j]*diag[i][j];
				if(fabs(diag[i][j])>maxod)
				{
					maxod=fabs(diag[i][j]);
					p=i;
					q=j;
				}
			}
		double x=2*diag[p][q], y=diag[p][p]-diag[q][q];
		double t=sqrt(x*x+y*y);
		double c=sqrt(0.5*(1+fabs(y)/t));
		double s=sqrt(0.5*(1-fabs(y)/t));
		if(x*y<0) s=-s;
		for(j=0;j<d;j++)
		{
			double temp=c*diag[p][j]+s*diag[q][j];
			diag[q][j]=-s*diag[p][j]+c*diag[q][j];
			diag[p][j]=temp;
		}
		for(j=0;j<d;j++)
		{
			double temp=c*diag[j][p]+s*diag[j][q];
			diag[j][q]=-s*diag[j][p]+c*diag[j][q];
			diag[j][p]=temp;
		}
		diag[p][q]=0;
		diag[q][p]=0;
		for(j=0;j<d;j++)
		{
			t=c*proj[p][j]+s*proj[q][j];
			proj[q][j]=-s*proj[p][j]+c*proj[q][j];
			proj[p][j]=t;
		}
	}
	while(od>0);
	for(i=0;i<prob->l;i++)
	{
		for(j=0;j<d;j++)
			buffer[j]=prob->x[i].values[j];
		for(j=0;j<d;j++)
		{
			prob->x[i].values[j]=0;
			for(k=0;k<d;k++)
				prob->x[i].values[j]+=proj[j][k]*buffer[k];
		}
	}
	return proj;
}


Vector QPparamMRF(Matrix XTR)
{
	const double rho=100;
	unsigned int i,j,k;
	Matrix XX = Matrix(XTR.rows(),XTR.rows(),0);
	Matrix invXX = Matrix(XTR.rows(),XTR.rows(),0);
	Matrix XPI = Matrix(XTR.rows(),XTR.cols(),0);
	for(i=0;i<XTR.rows();i++)
		for(j=0;j<XTR.rows();j++)
		{
			for(k=0;k<XTR.cols();k++)
				XX[i][j]+=XTR[i][k]*XTR[j][k];
		}
	InvDetCholesky(XX,invXX);
	for(i=0;i<XTR.rows();i++)
		for(j=0;j<XTR.cols();j++)
			for(k=0;k<XTR.rows();k++)
				XPI[i][j]+=invXX[i][k]*XTR[k][j];
	double XPIabsmax=-1,tolerance=1/(rho+XTR.rows()-1),normlambda=0;
	for(i=0;i<XTR.cols();i++) for(j=0;j<XTR.rows();j++) if(fabs(XPI[j][i])>XPIabsmax) XPIabsmax=fabs(XPI[j][i]);

	svm_problem *prob = new svm_problem;
	prob->l=XTR.rows()+XTR.cols();
	prob->x = new svm_node[prob->l];
	prob->y = new double[prob->l];
	signed char *y = new signed char[prob->l];
	double *alpha = new double[prob->l], *p = new double[prob->l];
	for(i=0;i<prob->l;i++)
	{
		p[i]=0;
		y[i]=1;
		prob->y[i]=1;
		prob->x[i].dim=XTR.rows();
		prob->x[i].values = new double[prob->x[i].dim];
		for(j=0;j<prob->x[i].dim;j++)
		{
			prob->x[i].values[j]=(i<XTR.rows() ? (i==j ? 1 : 0) : -XPI[j][i-XTR.rows()]/XPIabsmax);
			p[i]-=tolerance*prob->x[i].values[j];
		}
		alpha[i]=(i<XTR.rows() ? 1/double(XTR.rows())-tolerance : 0);
	}
	svm_parameter *param = new svm_parameter;
	param->svm_type = EPSILON_SVR;
	param->kernel_type=LINEAR;
	param->GGGH=param->PSB=param->HK=false;
	param->C=param->gamma=param->lambda=param->coef0=param->p=0;
	param->degree=0;
	param->cache_size = 100;	// default = 100
	param->eps = 1e-3;			// default = 1e-3
	param->shrinking = 1;		// default = 1
	Solver s;
	Solver::SolutionInfo si;
	s.Solve(prob->l,SVC_Q(*prob,*param,y),p,y,alpha,INF,INF,param->eps,&si,param->shrinking);
	Vector lambda;
	ofstream out("prova.csv");
	for(i=0;i<prob->l;i++)
	{
		if(alpha[i]>0) out << i+1 << ";" << alpha[i] << endl;
		delete[] prob->x[i].values;
		if(i<XTR.rows()) normlambda+=alpha[i];
		//if(i<XTR.rows()) lambda.push_back(alpha[i]+tolerance);
	}
	out << endl;
	for(i=0;i<XTR.rows();i++)
	{
		lambda.push_back((alpha[i]+tolerance)/(normlambda+XTR.rows()*tolerance));
		out << lambda[i] << ";";
	}
	out.close();
	delete[] prob->y,prob->x,p,y,alpha;
	delete prob,param;
	return lambda;
}


//Vector QPparamMRF(Matrix XTR)
//{
//	unsigned int i,j;
//	double Xabsmax=-1;
//	for(i=0;i<XTR.cols();i++) for(j=0;j<XTR.rows();j++) if(fabs(XTR[j][i])>Xabsmax) Xabsmax=fabs(XTR[j][i]);
//
//	svm_problem *prob = new svm_problem;
//	prob->l=XTR.rows()+XTR.cols();
//	prob->x = new svm_node[prob->l];
//	prob->y = new double[prob->l];
//	signed char *ones = new signed char[prob->l];
//	double *alpha = new double[prob->l];
//	double *zeros = new double[prob->l];
//	for(i=0;i<prob->l;i++)
//	{
//		prob->x[i].dim=XTR.cols();
//		prob->x[i].values = new double[prob->x[i].dim];
//		for(j=0;j<prob->x[i].dim;j++) prob->x[i].values[j]=(i<XTR.rows() ? XTR[i][j]/Xabsmax : (i==j+XTR.rows() ? -1 : 0));
//		prob->y[i]=1;
//		zeros[i]=0;
//		ones[i]=1;
//		alpha[i]=(i<XTR.rows() ? 1 : 0);
//	}
//	svm_parameter *param = new svm_parameter;
//	param->svm_type = EPSILON_SVR;
//	param->kernel_type=LINEAR;
//	param->GGGH=param->PSB=param->HK=false;
//	param->C=param->gamma=param->lambda=param->coef0=param->p=0;
//	param->degree=0;
//	param->cache_size = 100;	// default = 100
//	param->eps = 1e-3;			// default = 1e-3
//	param->shrinking = 1;		// default = 1
//	Solver s;
//	Solver::SolutionInfo si;
//	s.Solve(prob->l,SVC_Q(*prob,*param,ones),zeros,ones,alpha,INF,INF,param->eps,&si,param->shrinking);
//	Vector lambda;
//	ofstream out("prova.csv");
//	for(i=0;i<prob->l;i++)
//	{
//		if(alpha[i]>0) out << i+1 << ";" << alpha[i] << endl;
//		delete[] prob->x[i].values;
//		if(i<XTR.rows()) lambda.push_back(alpha[i]);
//	}
//	out << endl;
//	for(i=0;i<XTR.rows();i++) out << lambda[i] << ";";
//	out.close();
//	delete[] prob->y,prob->x,zeros,ones,alpha;
//	delete prob,param;
//	return lambda;
//}


//Vector QPparamMRF(Matrix XTR)
//{
//	unsigned int i,j;
//	svm_problem *prob = new svm_problem;
//	prob->l=XTR.rows()+XTR.cols();
//	double Xabsmax=-1;
//	for(i=0;i<XTR.cols();i++) for(j=0;j<XTR.rows();j++) if(fabs(XTR[j][i])>Xabsmax) Xabsmax=fabs(XTR[j][i]);
//	prob->x = new svm_node[prob->l];
//	prob->y = new double[prob->l];
//	signed char *zero = new signed char[prob->l];
//	double *alpha = new double[prob->l];
//	const double tolerance=1;
//	for(i=0;i<prob->l;i++)
//	{
//		prob->x[i].dim=XTR.rows();
//		prob->x[i].values = new double[prob->x[i].dim];
//		for(j=0;j<XTR.rows();j++) prob->x[i].values[j]=(i<XTR.cols() ? XTR[j][i]/Xabsmax : (i==XTR.cols()+j ? 1 : 0));
//		prob->y[i]=(i<XTR.cols() ? -1 : -tolerance);
//		zero[i]=0;
//	}
//	svm_parameter *param = new svm_parameter;
//	param->svm_type = EPSILON_SVR;
//	param->kernel_type=LINEAR;
//	param->GGGH=param->PSB=param->HK=false;
//	param->C=param->gamma=param->lambda=param->coef0=param->p=0;
//	param->degree=0;
//	param->cache_size = 100;	// default = 100
//	param->eps = 1e-3;			// default = 1e-3
//	param->shrinking = 1;		// default = 1
//	Solver s;
//	Solver::SolutionInfo *si;
//	s.Solve(prob->l,ONE_CLASS_Q(*prob,*param),prob->y,zero,alpha,INF,INF,param->eps,si,param->shrinking);
//	Vector lambda = Vector(XTR.rows(),0);
//	ofstream out("prova.csv");
//	for(i=0;i<prob->l;i++)
//	{
//		if(alpha[i]>0)
//		{
//			for(j=0;j<lambda.size();j++) lambda[j]+=alpha[i]*prob->x[i].values[j];
//			out << i+1 << ";" << alpha[i] << endl;
//		}
//		delete[] prob->x[i].values;
//	}
//	out.close();
//	delete[] prob->y,prob->x,zero,alpha;
//	delete prob,param;
//	return lambda;
//}



/*void svm_regression_resubstitution(struct svm_model *model, struct svm_problem *prob, double gain, char outname[80])
{
	longword i,j,k;
	ofstream out(outname);
	double rmse=0,mae=0;
	for(i=0;i<prob->l;i++)
	{
		double estimate=svm_predict(model,&prob->x[i]);
		if(model->param.GGGH)
		{
			double qstar=Kernel::k_function(&prob->x[i],&prob->x[i],model->param);
			Vector q = Vector(model->L.rows());
			Matrix QQ = Matrix(model->L.rows(),model->L.rows());
			for(j=0;j<model->L.rows();j++)
			{
				q[j]=Kernel::k_function(&model->SV[model->usv[j]],&prob->x[i],model->param);
				QQ[j][j]=model->L[j][j]-q[j]*q[j]/qstar;
			}
			for(j=0;j<model->L.rows();j++)
				for(k=j+1;k<model->L.rows();k++)
				{
					QQ[j][k]=model->L[j][k]-q[j]*q[k]/qstar;
					QQ[k][j]=QQ[j][k];
				}
			double approxvar=1/(1/qstar-2*ConjugateGrad(QQ,q)/(qstar*qstar));
			approxvar=(approxvar<0 ? 0 :approxvar);
			out << gain*prob->y[i] << "	" << gain*estimate << "	" << gain*sqrt(approxvar) << "	" << gain*sqrt(approxvar+model->noise_var) << endl;
		}
		else
			out << gain*prob->y[i] << "	" << gain*estimate << endl;
		mae+=fabs(prob->y[i]-estimate);
		rmse+=(prob->y[i]-estimate)*(prob->y[i]-estimate);
	}
	out << endl << "MAE	" << gain*mae/prob->l << endl << "RMSE	" << gain*sqrt(rmse/prob->l) << endl;
	out.close();
}*/

