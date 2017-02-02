/***************************************************************************
 *   Copyright (C) 2008 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI).   *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *         LIBRARY OF SUPPORT AND NUMERICAL FUNCTIONS AND CLASSES          *
 *                                                                         *
 ***************************************************************************/
#ifndef LIBRARY 
#define LIBRARY

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <limits>

using namespace std;

enum mrfoption { MANUAL,HK};
const word RatioGain=8;
const double PI=3.14159265359;
const double EuleroMascheroni=0.57721566490153286060;
const word icmiter=10;
const longword hkiter=1000000;
const double gammacoef[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
const double EdgePenalty=100;
const double nan=sqrt(-1.0);
const bool GraphCutSave=true;


/***************************************************************************************************************************************/
static void errorlog(char message[256])
{
	struct tm *newtime;
	time_t aclock;
	time(&aclock);
	newtime=localtime(&aclock);
	ofstream out("log.txt",ios::app);
	out << asctime(newtime) << message << endl;
	out.close();
	exit(0);
}


static longword longround(double x)
{
	if(x-floor(x)<0.5)
		return longword(floor(x));
	else
		return longword(ceil(x));
}


/***************************************************************************************************************************************/
typedef vector<double> Vector;


class Matrix: public vector<Vector>
{
	public:
		Matrix(size_t rows=0, size_t cols=1)
	: vector<Vector>(rows, Vector(cols)) {}
		Matrix(size_t rows, size_t cols, const double &value)
	: vector<Vector>(rows, Vector(cols,value)) {}
		Matrix(const Matrix &orig)
	: vector<Vector>(orig) {}
		size_t rows() const { return this[0].size();}
		size_t cols() const { return (*this)[0].size();}
};


class Tensor: public vector<Matrix>
{
	public:
		Tensor(size_t rows=0, size_t cols=0, size_t depth=1)
	: vector<Matrix>(rows,Matrix(cols,depth)) {}
		Tensor(size_t rows, size_t cols, size_t depth, const double &value)
	: vector<Matrix>(rows,Matrix(cols,depth,value)) {}
		Tensor(const Tensor &orig)
	: vector<Matrix> (orig) {}
		size_t rows() const { return this[0].size();}
		size_t cols() const { return (*this)[0].size();}
		size_t depth() const { return (*this)[0][0].size();}
};


/***************************************************************************************************************************************/
static double InvDetCholesky(Matrix A, Matrix &AI)
{
	short N=A.cols(),i,j,k;
	if(N!=A.rows())
		errorlog("InvDetCholesky. Rectangular matrix.");
	Vector p = Vector(N,0);
	Matrix L = Matrix(A);
	double sum,det=1;
	for(i=0;i<N;i++)
		for(j=i;j<N;j++)
		{
			sum=L[i][j];
			for(k=0;k<i;k++)
				sum-=L[i][k]*L[j][k];
			if(i==j)
			{
				if(sum<=0)
					errorlog("InvDetCholesky. Nonpositive definite matrix.");
				p[i]=sqrt(sum);
			}
			else
				L[j][i]=sum/p[i];
		}
	for(i=0;i<N;i++)					// invert the lower triangular matrix L
	{
		det*=p[i]*p[i];
		L[i][i]=1.0/p[i];
		for(j=i+1;j<N;j++)
		{
			sum=0;
			for(k=i;k<j;k++)
				sum-=L[j][k]*L[k][i];
			L[j][i]=sum/p[j];
			L[i][j]=0;
		}
	}
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
		{
			AI[i][j]=0;
			for(k=0;k<N;k++)
				AI[i][j]+=L[k][i]*L[k][j];
		}
	return det;
}


static double InvDetCholeskySecure(Matrix A, Matrix &AI)
{
	short N=A.cols(),i,j,k;
	if(N!=A.rows())
		errorlog("InvDetCholesky. Rectangular matrix.");
	Vector p = Vector(N,0);
	Matrix L = Matrix(A);
	double sum,det=1;
	bool valid=true;
	for(i=0;(i<N)&&valid;i++)
		for(j=i;(j<N)&&valid;j++)
		{
			sum=L[i][j];
			for(k=0;k<i;k++)
				sum-=L[i][k]*L[j][k];
			if(i==j)
			{
				if(sum<=0)
					valid=false;
				p[i]=sqrt(sum);
			}
			else
				L[j][i]=sum/p[i];
		}
	if(valid)
	{
		for(i=0;i<N;i++)					// invert the lower triangular matrix L
		{
			det*=p[i]*p[i];
			L[i][i]=1.0/p[i];
			for(j=i+1;j<N;j++)
			{
				sum=0;
				for(k=i;k<j;k++)
					sum-=L[j][k]*L[k][i];
				L[j][i]=sum/p[j];
				L[i][j]=0;
			}
		}
		for(i=0;i<N;i++)
			for(j=0;j<N;j++)
			{
				AI[i][j]=0;
				for(k=0;k<N;k++)
					AI[i][j]+=L[k][i]*L[k][j];
			}
	}
	else
	{
		det=1;
		for(i=0;i<N;i++)
		{
			det*=(A[i][i]>0 ? A[i][i] : 1);
			AI[i][i]=(A[i][i]>0 ? 1/A[i][i] : 1);
			for(j=i+1;j<N;j++)
			{
				AI[i][j]=0;
				AI[j][i]=0;
			}
		}
	}
	return det;
}


static Matrix InvLowerCholesky(Matrix A)
{
	short N=A.cols(),i,j,k;
	if(N!=A.rows())
		errorlog("InvLowerCholesky. Rectangular matrix.");
	Vector p = Vector(N,0);
	Matrix L = Matrix(A);
	double sum;
	for(i=0;i<N;i++)
		for(j=i;j<N;j++)
		{
			sum=L[i][j];
			for(k=0;k<i;k++)
				sum-=L[i][k]*L[j][k];
			if(i==j)
			{
				if(sum<=0)
					errorlog("InvLowerCholesky. Nonpositive definite matrix.");
				p[i]=sqrt(sum);
			}
			else
				L[j][i]=sum/p[i];
		}
	for(i=0;i<N;i++)					// invert the lower triangular matrix L
	{
		L[i][i]=1.0/p[i];
		for(j=i+1;j<N;j++)
		{
			sum=0;
			for(k=i;k<j;k++)
				sum-=L[j][k]*L[k][i];
			L[j][i]=sum/p[j];
			L[i][j]=0;
		}
	}
	return L;
}


static double DetCholesky(Matrix A)
{
	short N=A.cols(),i,j,k;
	if(N!=A.rows())
		errorlog("DetCholesky. Rectangular matrix.");
	Vector p = Vector(N,0);
	Matrix L = Matrix(A);
	double sum,det=1;
	for(i=0;i<N;i++)
		for(j=i;j<N;j++)
		{
			sum=L[i][j];
			for(k=0;k<i;k++)
				sum-=L[i][k]*L[j][k];
			if(i==j)
			{
				if(sum<=0)
					errorlog("DetCholesky. Nonpositive definite matrix.");
				p[i]=sqrt(sum);
			}
			else
				L[j][i]=sum/p[i];
		}
	for(i=0;i<N;i++)
		det*=p[i]*p[i];
	return det;
}


/***************************************************************************************************************************************/
static Vector EigenPower(Matrix A)
{
	word N=A.cols(),i,j,iter,maxiter=100;
	if(N!=A.rows())
		errorlog("EigenPower. Rectangular matrix.");
	Vector Y = Vector(N,1);
	double delta=DBL_MAX,eps_stop=0.0001;
	for(iter=0;(iter<maxiter)&&(delta>eps_stop);iter++)
	{
		Vector W = Vector(N,0);
		double Winfty=0;	// ha un segno per tener conto di convergenza a meno del segno
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
				W[i]+=A[i][j]*Y[j];
			if(fabs(W[i])>fabs(Winfty))
				Winfty=W[i];
		}
		delta=0;
		if(Winfty!=0)
			for(i=0;i<N;i++)
			{
				double Ynew=W[i]/Winfty;
				if(fabs(Y[i]-Ynew)>delta)
					delta=fabs(Y[i]-Ynew);
				Y[i]=Ynew;
			}
		else
		{
			Y = Vector(N,2);
			delta=0;
		}
	}
	return Y;
}


/***************************************************************************************************************************************/
static double ConjugateGrad(Matrix Q, Vector L)
// data una matrice quadrata simmetrica Q ed un vettore L, risolve min(X'QX/2-L'X) (occhio ai segni!),
// scrive la soluzione in X e restituisce in uscita il valore minimo della funzione-obiettivo
// (v. Stoer, Bulirsch, pp. 606-sgg.).
{
	word N=Q.cols(),i,j;
	Vector X = Vector(N,0), P = Vector(L), R = Vector(L), Rnew = Vector(N);
	double delta=DBL_MAX,F=0;
	const double eps_stop=0.00001;
	for(word iter=0;(iter<N+1)&&(delta>eps_stop);iter++)
	{
		delta=0;
		double RR=0,PQP=0,RRnew=0;				// calcolo della "learning rate" alpha
		for(i=0;i<N;i++)
		{
			RR+=R[i]*R[i];
			PQP+=P[i]*P[i]*Q[i][i];
			for(j=i+1;j<N;j++)
				PQP+=2*P[i]*P[j]*Q[i][j];
		}
		double alpha=RR/PQP;
		for(i=0;i<N;i++)
			X[i]+=alpha*P[i];					// passo di discesa
		for(i=0;i<N;i++)
		{
			double QP=0;
			for(j=0;j<N;j++)
				QP+=Q[i][j]*P[j];
			Rnew[i]=R[i]-alpha*QP;
			RRnew+=Rnew[i]*Rnew[i];
		}
		double beta=RRnew/RR;					// calcolo del coefficiente beta
		for(i=0;i<N;i++)						// aggiornamento delle direzioni P e Q
		{
			P[i]=Rnew[i]+beta*P[i];
			R[i]=Rnew[i];
			if(fabs(P[i])>delta) delta=fabs(P[i]);
		}
	}
	for(i=0;i<N;i++)							// calcolo del valore finale (minimo) della funzione-obiettivo
	{
		F+=0.5*Q[i][i]*X[i]*X[i]-L[i]*X[i];
		for(j=i+1;j<N;j++)
			F+=Q[i][j]*X[i]*X[j];
	}
	return F;
}


/***************************************************************************************************************************************/
static Vector HoKashyap(Matrix XTR, longword maxiter)
{
	word col=XTR.rows(),row=XTR.cols(),i,j;
	longword r,iter;
	const double rho=1,eps=0.0001;
	Vector w = Vector(col,1);
	Vector wold = Vector(col,1);
	Vector b = Vector(row,0);
	Matrix XX = Matrix(col,col,0);
	Matrix invXX = Matrix(col,col,0);
	double delta=DBL_MAX;
	for(i=0;i<col;i++)
		for(j=0;j<col;j++)
		{
			for(r=0;r<row;r++)
				XX[i][j]+=XTR[i][r]*XTR[j][r];
		}
	InvDetCholesky(XX,invXX);
	for(r=0;r<row;r++)
	{
		for(i=0;i<col;i++)
			b[r]+=w[i]*XTR[i][r];
		if(b[r]<0)
			b[r]=1;
	}
	for(iter=0;(iter<maxiter)&&(delta>eps);iter++)
	{
		cout << iter+1 << " ";
		flush(cout);
		for(i=0;i<col;i++)
		{
			wold[i]=w[i];
			w[i]=0;
			for(r=0;r<row;r++)
				for(j=0;j<col;j++)
					w[i]+=invXX[i][j]*XTR[j][r]*b[r];
		}
		delta=0;
		for(r=0;r<row;r++)
		{
			double e=-b[r];
			for(i=0;i<col;i++)
				e+=w[i]*XTR[i][r];
			if(e>0)
			{
				b[r]+=rho*e;
				if(rho*e>delta) delta=rho*e;
			}
		}
		for(i=0;i<col;i++)
			if(fabs(w[i]-wold[i])>delta) delta=fabs(w[i]-wold[i]);
	}
	return w;
}


/***************************************************************************************************************************************/
static double LnGamma(double x)
{
	if(x<=0)
		errorlog("LnGamma. Unsupported nonpositive argument.");
	double y=x,t=x+5.5,s=1.000000000190015;
	t-=(x+0.5)*log(t);
	word j;
	for(j=0;j<6;j++)
		s+=gammacoef[j]/++y;
	return -t+log(2.5066282746310005*s/x);
}


static double Gamma(double x)
{
	return exp(LnGamma(x));
}


static double PolyGamma(word n, double x)
{
	if(x<=0)
		errorlog("PolyGamma. Unsupported nonpositive argument.");
	word j;
	x--;
	double theta0=1.000000000190015,theta1=0,theta2=0,theta3=0,theta4=0,y=x,Z=x+5.5;
	for(j=0;j<6;j++)
	{
		theta0+=gammacoef[j]/++y;
		theta1-=gammacoef[j]/(y*y);
		theta2+=2*gammacoef[j]/(y*y*y);
		theta3-=6*gammacoef[j]/(y*y*y*y);
		theta4+=24*gammacoef[j]/(y*y*y*y*y);
	}
	switch(n)
	{
		case 0:
			return log(Z)+(x+0.5)/Z+theta1/theta0-1;
			break;
		case 1:
			return (x+10.5)/(Z*Z)+(theta2*theta0-theta1*theta1)/(theta0*theta0);
			break;
		case 2:
			return -(x+15.5)/(Z*Z*Z)+(theta3*theta0*theta0-3*theta2*theta1*theta0+2*theta1*theta1*theta1)/(theta0*theta0*theta0);
			break;
		case 3:
			return (2*x+41)/(Z*Z*Z*Z)+(theta4*theta0*theta0*theta0-4*theta3*theta1*theta0*theta0+3*theta2*theta0*(4*theta1*theta1-theta2*theta0)-6*theta1*theta1*theta1*theta1)/(theta0*theta0*theta0*theta0);
			break;
		default:
			errorlog("PolyGamma. Unsupported order.");
			return 0;
			break;
	}
}


static double InvPolyGamma(word n, double psi)
{
	if(n>2)
		errorlog("InvPolyGamma. Unsupported order.");
	else if((n==1)&&(psi<=0))
		errorlog("InvPolyGamma. Unsupported nonpositive argument.");
	else if((n==2)&&(psi>=0))
		errorlog("InvPolyGamma. Unsupported nonnegative argument.");
	const double eps=0.0001;
	longword iter=0,maxiter=100;
	double delta=1,x=-10;
	for(iter=0;(iter<maxiter)&&(fabs(delta)>eps);iter++)
	{
		delta=PolyGamma(n,exp(x))-psi;
		x-=delta*exp(-x)/PolyGamma(n+1,exp(x));
	}
	return exp(x);
}


static double Hermite(double x, word order)
{
	switch(order)
	{
	case 0:  return 1;													break;
	case 1:  return x;													break;
	case 2:  return x*x-1;												break; 
	case 3:	 return x*x*x-3*x;											break;
	case 4:  return x*x*x*x-6*x*x+3;									break;
	default: return x*Hermite(x,order-1)-(order-1)*Hermite(x,order-2);	break;
	}
}


static void gser(double *gamser, double a, double x, double *gln)
{
	int n,ITMAX=100;
	const double EPS1=3.0e-7;
	double sum,del,ap;
	*gln=LnGamma(a);
	if (x <= 0.0)
	{
		if (x < 0.0) exit(1);
		*gamser=0.0;
		return;
	}
	else
	{
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++)
		{
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS1)
			{
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
	}
	return;
	}
}


static void gcf(double *gammcf, double a, double x, double *gln)
{
	int i,ITMAX=100;
	const double FPMIN=1.0e-30,EPS=1.0e-16;
	double an,b,c,d,del,h;
	*gln=LnGamma(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++)
	{
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	//if (i > ITMAX) exit(1);
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}


static double IncGammaP(double a, double x)
{
	double gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) return 0;
	else if (x < (a+1.0))
	{
		gser(&gamser,a,x,&gln);
		return gamser;
	}
	else
	{
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}


static double InvGammaP(double a, double P)
{
	const double eps=0.0001;
	longword iter=0,maxiter=100;
	double delta=1,x=-10;
	for(iter=0;(iter<maxiter)&&(fabs(delta)>eps);iter++)
	{
		delta=IncGammaP(a,exp(x))-P;
		x-=delta*exp(LnGamma(a)-a*x-exp(x));
	}
	return exp(x);
}


static double erf(double x)
{
	if(x>0) return IncGammaP(0.5,x*x);
	else if(x<0) return -IncGammaP(0.5,x*x);
	else return 0;
}


static double Queue(double x)
{
	return 0.5*(1-erf(x*sqrt(0.5)));
}


static double NormalCDF(double x)
{
	return 0.5+(x<0 ? -0.5 : 0.5)*IncGammaP(0.5,x*x*0.5);
}



/***************************************************************************************************************************************/
static double PottsBesag(Image *energymap)
{
	/*if(energymap->get_bands()>2)
		errorlog("PottsBesag. Unsupported number of channels.");*/
	word w=energymap->get_width(),h=energymap->get_height(),b=energymap->get_bands(),m,n,r,iter,maxiter=25;
	const double eps=0.01;
	double beta=0,deltabeta=DBL_MAX;
	for(iter=0;(iter<maxiter)&&(deltabeta>eps);iter++)
	{
		double phiprime=0,phisecond=0;
		for(m=1;m<w-1;m++)
			for(n=1;n<h-1;n++)
			{
				double A=0,B=1,C=0,lambda=exp(beta);
				for(r=0;r<b;r++)
				{
					double eta=energymap->get(m,n,r),e=exp(lambda*eta);
					A+=eta*e;
					B+=e;
					C+=eta*eta*e;
				}
				if(B>0)
				{
					phiprime+=lambda*A/B;
					phisecond+=lambda*A/B+lambda*lambda*(C*B-A*A)/(B*B);
				}
			}
		if(phisecond!=0)
		{
			deltabeta=fabs(phiprime/phisecond);
			beta-=phiprime/phisecond;
		}
		else
			deltabeta=0;
	}
	return exp(beta);
}


static bool Finite(double x)
{
	return (_finite(x)!=0);
}


static int DoubleCompare(const void *a, const void *b)
{
	double A=*(double*)a,B=*(double*)b;
	if(A>B) return 1;
	else if(A<B) return -1;
	else return 0;
}


static void QuickSort(double *v, word count)
{
	qsort(v,count,sizeof(double),DoubleCompare);
}
#endif