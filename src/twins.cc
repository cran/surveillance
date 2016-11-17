/*******************************************************************
 * Authors:
 *   Mathias Hofmann
 *   Michael Hoehle <hoehle@stat.uni-muenchen.de>
 *   Volker Schmid
 * Contributors:
 *   Michaela Paul
 *   Daniel Sabanes Bove <daniel.sabanesbove@ifspm.uzh.ch>
 *   Sebastian Meyer <sebastian.meyer@ifspm.uzh.ch>
 * History:
 *   July 2016 (SM) -- dropped deprecated "register" storage class specifier
 *   April 2012 (SM) -- replaced exit() calls by Rf_error()
 *   March 2012 (DSB) -- changed long types to int to be in accordance with R
 *                       (we observed bad allocations in 64 bit machines)
 *   May 2010 (DSB) -- modified from Oct 2008
 *
 * Markov Chain Monte Carlo (MCMC) estimation in the Branching Process
 * like Epidemic Model. Instead of a slow R solution this code
 * provides a faster C++ solution. Can be invoked through R or be
 * programmed as a librrary. This code uses the Gnu Scientific Library
 * (GSL) available from http://sources.redhat.com/gsl/
 *
 * For now this code is quick & dirty. A more OO framework would be nice
 * to enable better programming, but this will probably be speedwise slower.
 *******************************************************************/

#include <iostream>
#include <fstream>

/*New C++ uses header iostream (without the .h) followed by a namespace*/
using namespace std;

#include <math.h>

/* Replaced calls to GSL with functions from the R API */
#include <R.h>
#include <Rmath.h>
/*wrappers to what used to be GSL functions*/
#include "gsl_wrappers.h"

//  Dynamic_2d_array class by David Maisonave (609-345-1007) (www.axter.com)
//  Description:
//	The dynamic array class listed below is more efficient then other 
//	similar classes that use temporary objects as return types, or use 
//	an std::vector as a return type.
//
//	It's also more compatible with a C style 2D array, in that the 
//	array is in one continuous memory block. This makes it possible 
//	to pass this array object to a C Function that has a C-Style 
//	2D array for a parameter.
//  Example usage:
/*

Dynamic_2d_array<int> MyIntArray(12, 34);
MyIntArray[0][1] = 123;
cout << MyIntArray[0][1] << endl;

*/

template < class T >
class Dynamic_2d_array
{
public:
  // constructor
  Dynamic_2d_array(size_t row, size_t col) :
    m_row(row),
    m_col(col), 
    m_data((row!=0 && col!=0) ? new T[row*col] : NULL)
  {}

  // copy ctr
  Dynamic_2d_array(const Dynamic_2d_array& src) :
    m_row(src.m_row),
    m_col(src.m_col), 
    m_data((src.m_row!=0 && src.m_col!=0) ? new T[src.m_row*src.m_col] : NULL)
  {
    for(size_t r=0; r<m_row; ++r)
      for(size_t c=0; c<m_col; ++c) 
	(*this)[r][c] = src[r][c];
  }
  
  // destructor
  ~Dynamic_2d_array()
  {
    delete[] m_data;
  }
  
  // non-const access
  inline T* operator[](size_t i) 
  {
    return (m_data + (m_col*i));
  }
  
  // const access
  inline T const*const operator[](size_t i) const 
  {
    return (m_data + (m_col*i));
  }

private:
  const size_t m_row;
  const size_t m_col;
  T* m_data; 
};


// Note that the class Dynamic_2d_array automatically allocates and
// deallocates the memory.
typedef Dynamic_2d_array<long> LongMatrix;
typedef Dynamic_2d_array<double> DoubleMatrix;
typedef Dynamic_2d_array<int> IntMatrix;



// Analogous class for vectors (== 1D arrays)
template < class T >
class Dynamic_1d_array
{
public:
  // constructor
  Dynamic_1d_array(size_t length) :
    m_length(length),
    m_data((length !=0) ? new T[length] : NULL)
  {}

  // copy ctr
  Dynamic_1d_array(const Dynamic_1d_array& src) :
    m_length(src.m_length),
    m_data((src.m_length!=0) ? new T[src.m_length] : NULL)
  {
    for(size_t i=0; i<m_length; ++i)
      (*this)[i] = src[i];
  }
  
  // destructor
  ~Dynamic_1d_array()
  {
    delete[] m_data;
  }
  
  // non-const access
  inline T& operator[](size_t i) 
  {
    return m_data[i];
  }
  
  // const access
  inline T const& operator[](size_t i) const 
  {
    return m_data[i];
  }

private:
  const size_t m_length;
  T* m_data; 
};


// Note that the class Dynamic_1d_array automatically allocates and
// deallocates the memory.
typedef Dynamic_1d_array<long> LongVector;
typedef Dynamic_1d_array<double> DoubleVector;
typedef Dynamic_1d_array<int> IntVector;





/************************************
  Globals
*************************************/

/*Setup params*/
int overdispersion;
int varnu;
int la_rev;
int K_geom;
int la_estim;
int nu_trend;
int theta_pred_estim;
int xi_estim;
int delta_rev;
int xi_estim_delta;
int epsilon_rev;
int xi_estim_epsilon;
int xi_estim_psi;
double psiRWSigma = 0.25;
double xRWSigma = 0.25;
double taubetaRWSigma = 0.25;

/*Priors*/
double alpha_lambda = 1.0;
double beta_lambda = 1.0;

double alpha_xi = 1.0;
double beta_xi = 1.0;

double p_K = 1.0;

double alpha_nu = 1.0;
double beta_nu = 1.0;

double alpha_psi = 1.0;
double beta_psi = 10.0;


double alpha_a=1;
double alpha_b=0.001;
double beta_a=1.0;
double beta_b=.00001;
double gamma_a=1;
double gamma_b=0.001;
double delta_a=1;
double delta_b=0.001;
double epsilon_a=1;
double epsilon_b=0.001;



/*********************************************************************
 * Compute sum from 1 to I and 1 to n of a vektor with indices 0,...,I 
 * of a vektor with indices 0,...,n
 * Parameters:
 *
 * X a vector with indices 0,..,I of a vector with indices 0,...,n
 * I "length" of vector (true length due to zero indice is I+1)
 *********************************************************************/
double sumIn(const LongMatrix& X, int I, int n) {
  double res = 0;
  for (int i=1; i<=I; i++){
    for (int t=1; t<=n; t++) {
      res += X[i][t]; 
    }
  }
  return(res);
}

/*********************************************************************
 * Compute sum from 1 to I and 1 to n of a vektor with indices 0,...,I 
 * of a vektor with indices 0,...,n
 * This is the double version
 * Parameters:
 *
 * X a vector with indices 0,..,I of a vector with indices 0,...,n
 * I "length" of vector (true length due to zero indice is I+1)
 *********************************************************************/
double sumIn(const DoubleMatrix& X, int I, int n) {
  double res = 0;
  for (int i=1; i<=I; i++){
    for (int t=1; t<=n; t++) {
      res += X[i][t]; 
    }
  }
  return(res);
}

/*********************************************************************
 * Compute sum from 1 to I and 1 to n of a vektor with indices 0,...,I 
 * of a vektor with indices 0,...,n
 * Parameters:
 *
 * X a vector with indices 0,..,I of a vector with indices 0,...,n
 * I "length" of vector (true length due to zero indice is I+1)
 *********************************************************************/
double sumIn2(const LongMatrix& X, int I, int n) {
  double res = 0;
  for (int i=1; i<=I; i++){
    for (int t=2; t<=n; t++) {
      res += X[i][t]; 
    }
  }
  return(res);
}


/*********************************************************************
 * Compute sum from 1 to I and 1 to n of a vektor with indices 0,...,I 
 * of a vektor with indices 0,...,n
 * This is the double version
 * Parameters:
 *
 * X a vector with indices 0,..,I of a vector with indices 0,...,n
 * I "length" of vector (true length due to zero indice is I+1)
 *********************************************************************/
double sumIn2(const DoubleMatrix& X, int I, int n) {
  double res = 0;
  for (int i=1; i<=I; i++){
    for (int t=2; t<=n; t++) {
      res += X[i][t]; 
    }
  }
  return(res);
}


/*********************************************************************
 * Compute sum from 1 to I of a vektor with indices 0,...,I 
 * of a vektor with indices 0,...,n
 * Parameters:
 *
 * X a vector with indices 0,..,I of a vector with indices 0,...,n
 * I "length" of vector (true length due to zero indice is I+1)
 *********************************************************************/
double sumI1(const LongMatrix& X, int I, int t) {
  double res = 0;
  for (int i=1; i<=I; i++) { res += X[i][t]; }
  return(res);
}

/*********************************************************************
 * Compute sum from 1 to I of a vektor with indices 0,...,I 
 * of a vektor with indices 0,...,n
 * This is the double version
 * Parameters:
 *
 * X a vector with indices 0,..,I of a vector with indices 0,...,n
 * I "length" of vector (true length due to zero indice is I+1)
 *********************************************************************/
double sumI1(const DoubleMatrix& X, int I, int t) {
  double res = 0;
  for (int i=1; i<=I; i++) { res += X[i][t]; }
  return(res);
}

/*********************************************************************
 * factorial function
 *********************************************************************/
long factorial(long x){
  long fac=1;
  if(x<0){ Rf_error("negative value passed to factorial function\n");}
  else{
    if(x==0){fac=1;}
    else{
      for(int i=1;i<=x;i++){
        fac*=i;
      }
    }
  }
  return(fac);
}

/*********************************************************************
 * logit function
 *********************************************************************/
double logit(double y){
  if(y <= 0 || y >= 1){
    Rf_error("y <= 0 or y >= 1 in logit function.\n");
  }
  double logit;
  logit = log(y/(1-y));
    return(logit);
}




/*********************************************************************
 * inverse logit function
 *********************************************************************/
double invlogit(double y){
  double invlogit;
  invlogit = 1/(1 + exp(-y));
    return(invlogit);
}


/*********************************************************************
 * inverse logit function diff.
 *********************************************************************/
double invlogitd(double y){
  double invlogitd;
  invlogitd = exp(-y)/pow((1.0 + exp(-y)),2);
    return(invlogitd);
}




/*********************************************************************
 * Makes one Metropolis-Hastings update step, log-scale
 *********************************************************************/
double updateMHlog(double &par, double parStar, double logFpar, double logFparStar, double &acceptedpar) {
  double accpar = exp(logFparStar - logFpar); 
  if (gsl_rng_uniform() <= accpar) {par = parStar; acceptedpar++;}
  return(0);
}

/*********************************************************************
 * Makes one Metropolis-Hastings update step
 *********************************************************************/
double updateMH(double &par, double parStar, double Fpar, double FparStar, double &acceptedpar) {
  double accpar = FparStar/Fpar; 
  if (gsl_rng_uniform() <= accpar) {par = parStar; acceptedpar++;}
  return(0);
}






/*********************************************************************
 * Tunes a parameter
 *********************************************************************/
double tune(double& parameter, double accepted, double samples, double& tunepar, double a=0.3, double b=0.4){
      tunepar=1;     
      if ((accepted/samples>a) && (accepted/samples<b)) {
        tunepar=0;
      } else if (accepted/samples>b) {
        parameter *= 1.5;
      } else if (accepted/samples<a) {
        parameter /= 2.0;
      }
      return(0);
}

double ABS(double x)
{
  if (x>0){return x;}else{return -x;}
}
double MIN(double a, double b)
{
  if (a<b){return a;}else{return b;}
}

double xMx(double* Q, double* x, int noa, int b)
{
  double erg=0.0;
  
  for (int i=0; i<noa; i++)
    {
      for (int j=0; j<noa; j++)
	{
	  if (ABS(i-j) < b)
	    {
	      erg += x[i]*x[j]*Q[(int)(MIN(i,j)*b+ABS(i-j))];
	      if (i==j)
		{
		  erg -= 0.5*x[i]*x[j]*Q[(int)(MIN(i,j)*b+ABS(i-j))];
		}
	    }
	}
    }

  return erg;
}

double xMx2(double* Q, double* x, int n, int b)
{
  double erg=0.0;
  
  for (int i=0; i<=n; i++)
    {
      for (int j=0; j<=n; j++)
	{
	  if (ABS(i-j) < b)
	    {
	      erg += x[i]*x[j]*Q[(int)(MIN(i,j)*b+ABS(i-j))];
	    }
	}
    }

  return erg;
}


/* BERECHNET A-1, k ist matrixlaenge*/
void invers(double* A,int k)
{
  DoubleVector ergebnis(k * k);
  
  if (k==1)
    {
      ergebnis[0] = 1.0/A[0];
    }

  if (k==2)
    {
      double det = A[0]*A[3]-A[1]*A[2];
      ergebnis[0] = A[3]/det;
      ergebnis[1] = 0.0-A[1]/det;
      ergebnis[2] = 0.0-A[2]/det;
      ergebnis[3] = A[0]/det;
    }

  if (k>2)
    {
      REprintf("Error in the twins.cc function invers()\n");
    }

  for (int i=0; i< k*k; i++)
    {
      A[i]=ergebnis[i];
    }
  
  return;
}

void mxschreibe(double* A, int a, int b)
{
  for (int i=0; i<a; i++)
    {
      for (int j=0; j<b; j++)
	{
	  Rprintf("%f ", A[i*b+j]);
	}
      Rprintf("\n");
    }
  Rprintf("\n");
  return;
}


int mxcheck(int n, const IntMatrix& matrix)
{
  int zs=0;
  for (int i=0; i<n; i++)
    {
      zs=0;
      for (int j=0; j<n; j++)
	{
	  zs+=matrix[i][j];
	  if (matrix[i][j]!=matrix[j][i])
	    {
	      REprintf("Error: Matrix is not symmetric! (Row: %d, Column %d\n",i,j); 
	      return 1;
	    }
	}
      if (zs != 0)
	{
	  REprintf("Error: Row sum not zero in row %d",i,"\n");
	  return 1;
	}
    }
  return 0;
}





/*updatealphabeta
  Erzeugt Normalverteilten Zufallsvektor der Laenge noa*/
void gausssample(double* temp, int noa)
{
  
  for (int i=0; i< noa; i++)
    {
      temp[i]=gsl_ran_gaussian(1);
    }
  
  return;
}

double hyper(int rw, double* theta, double k_a, double k_b, int n)
{

  double aa;
  double bb;
  double kappa_neu=0;
  
  if (rw==1)
    {
      double summe = 0.0;
      
      aa = k_a + 0.5 * double(n-1-rw);
      
      for (int i=3; i <= n; i++)
	{
	  summe = summe + (theta[i] - theta[i-1]) * (theta[i] - theta[i-1]);
	}
      bb = k_b + 0.5 * summe;
      
      kappa_neu = gsl_ran_gamma(aa,1/bb);
    }
  if (rw==2)
    {
      
      double dopp_diff;
      double summe = 0.0;
      
      aa = k_a + 0.5 * double(n-1-rw);
      
      for (int i=3; i < n; i++)
	{
	  dopp_diff = theta[i-1] - 2*theta[i] + theta[i+1];
	  summe = summe + dopp_diff * dopp_diff;
	}
      bb = k_b + 0.5 * summe;
      
      
      
      kappa_neu = gsl_ran_gamma( aa, 1/bb);
    }
  if (rw==0)
    {
      double summe = 0.0;
      
      aa = k_a + 0.5 * double(n-1-rw);
      
      for (int i=2; i <= n; i++){
	summe = summe + (theta[i] * theta[i]);
      }
      
      bb = k_b + 0.5 * summe;
      
      kappa_neu = gsl_ran_gamma(aa,1/bb);
    }
  
  return kappa_neu;
}


double update_tau_alpha(const DoubleVector& alpha, int I, double aa, double bb, DoubleVector& xreg)
{
  
  aa += double(I);
  for (int i=1; i<=I; i++)
    {
      bb += (alpha[i]-xreg[i])*(alpha[i]-xreg[i]);
    }
  double neu=gsl_ran_gamma( aa, 1/bb);
  return neu;
}


double update_tau_gamma(const DoubleVector& alpha, int ncov, double aa, double bb)
{
  
  aa += double(ncov);
  for (int i=0; i<ncov; i++)
    {
      bb += alpha[i]*alpha[i];
    }
  double neu=gsl_ran_gamma( aa, 1/bb);
  return neu;
}

void berechneQ(double* temp, int age_block, double kappa, int noa,int nop, double delta)
{

  if (age_block==1)
    {
      int index=0;
      temp[index]=kappa+delta*nop;
      index++;
      temp[index]=-kappa;
      index++;
      for (int i=1; i<noa-1; i++)
	{
	  temp[index]=2*kappa+delta*nop;
	  index++;
	  temp[index]=-kappa;
	  index++;
	}
      temp[index]=kappa+delta*nop;
    }
  if (age_block==2)
    {
      int index=0;
      temp[index]=kappa+delta*nop;
      index++;
      temp[index]=-2*kappa;
      index++;
      temp[index]=kappa;
      index++;
      temp[index]=5*kappa+delta*nop;
      index++;
      temp[index]=-4*kappa;
      index++;
      temp[index]=kappa;
      index++;
      for (int i=2; i<noa-2; i++)
	{
	  temp[index]=6*kappa+delta*nop;
	  index++;
	  temp[index]=-4*kappa;
	  index++;
	  temp[index]=kappa;
	  index++;
	}
      temp[index]=5*kappa+delta*nop;
      index++;
      temp[index]=-2*kappa;
      index++;
      index++;
      temp[index]=kappa+delta*nop;
    }

  return;
}



double sumg(int ncov, const DoubleMatrix& xcov, DoubleVector& gamma, int t, int scov)
{
  double sum=0;
  for (int i=scov; i<ncov; i++)
    {
      sum += xcov[i][t]*gamma[i];
    }
  return sum;
}




void alphaupdate(DoubleVector& gamma, DoubleVector& alpha, DoubleVector& beta, DoubleVector& delta, const DoubleMatrix& lambda, double p, int I, int n, const LongMatrix& Y, const LongMatrix& X, long& acc_alpha, double taualpha, int ncov, const DoubleMatrix& xcov, DoubleVector& xreg, const DoubleMatrix& omega, const DoubleMatrix& omegaX, int scov, int mode){

  for (int i=1; i<=I; i++)
    {
  
      double tau=taualpha;
      double my=0;
      for (int t=2; t<=n; t++)
	{
	  tau+=omegaX[i][t]*delta[t]*exp(sumg(ncov,xcov,gamma,t,scov)+alpha[i]+beta[t]);
	  my+=(X[i][t]);
	  my-=(1-alpha[i])*omegaX[i][t]*delta[t]*exp(sumg(ncov,xcov,gamma,t,scov)+alpha[i]+beta[t]);
	}
      my=my+xreg[i]*taualpha;
      my=my/tau;
      double alphaneu=gsl_ran_gaussian(sqrt(1/tau))+my;
      double tauneu=taualpha;
      double myneu=0;
      for (int t=2; t<=n; t++)
	{
	  tauneu+=omegaX[i][t]*delta[t]*exp(sumg(ncov,xcov,gamma,t,scov)+alphaneu+beta[t]);
	  myneu+=(X[i][t]);
	  myneu-=(1-alphaneu)*omegaX[i][t]*delta[t]*exp(sumg(ncov,xcov,gamma,t,scov)+alphaneu+beta[t]);
	}
      myneu=myneu+xreg[i]*taualpha;
      myneu=myneu/tauneu;
      double akzw=0.5*log(tauneu/(2*PI))-0.5*tauneu*(alphaneu-myneu)*(alphaneu-myneu); /* log Proposalw. alt|neu*/
      akzw -= ((0.5*log(tau/(2*PI))-0.5*tau*(alpha[i]-my)*(alpha[i]-my))); /* log Proposalw. neu|alt*/
      akzw += (-0.5*taualpha*(alpha[i]-xreg[i])*(alpha[i]-xreg[i]));
      akzw -= (-0.5*taualpha*(alphaneu-xreg[i])*(alphaneu-xreg[i]));
      for (int t=2; t<=n; t++)
	{
	  akzw += (alpha[i]*X[i][t]-omegaX[i][t]*delta[t]*exp(sumg(ncov,xcov,gamma,t,scov)+alpha[i]+beta[t]));
	  akzw -= (alphaneu*X[i][t]-omegaX[i][t]*delta[t]*exp(sumg(ncov,xcov,gamma,t,scov)+alphaneu+beta[t]));
	}

      if (exp(akzw) >= gsl_rng_uniform())
	{
	  alpha[i]=alphaneu;
	  acc_alpha += 1;
	}
    }

  return;
}



void erzeuge_b_Q(DoubleVector& gamma  , double* my, double* Q, const DoubleVector& alpha, DoubleVector& delta, DoubleVector& beta, const LongMatrix& X, const LongMatrix& Z, const LongMatrix& Y,
		 int n, int I, double taubeta, int rw, const DoubleMatrix& lambda, double p, const DoubleMatrix& xcov, int ncov, const DoubleMatrix& omega, const DoubleMatrix& omegaX,int scov, int mode)
{
  if (mode==1)
    {
      /* b-vektor des Proposals*/
      for (int t=0;t<n; t++)
	{
	  my[t]=0.0;
	  for (int i=1; i<=I; i++)
	    {
	      my[t]+=X[i][t+2];
	      my[t]-=(1-beta[t+2])*omegaX[i][t+2]*delta[t+2]*(exp(sumg(ncov,xcov,gamma,t+2,scov)+alpha[i]+beta[t+2]));
	    }
	}
    }
  if (mode==2)
    {
      /* b-vektor des Proposals*/
      for (int t=2;t<=n; t++)
	{
	  my[t-2]=0.0;
	  for (int i=1; i<=I; i++)
	    {
	      my[t-2]+=Y[i][t];
	      my[t-2]-=(1-beta[t])*(omega[i][t]*Z[i][t-1]*exp(sumg(ncov,xcov,gamma,t,scov)+alpha[i]+beta[t]));
	    }
	}
      
      
    }
  
  /* Praezisionsmatrix*/
  berechneQ(Q, rw, taubeta, n, 1, 0.0);
  
  if (mode==1)
    {
      for (int i=1; i<=I; i++)
	{
	  for (int t=0; t<n; t++)
	    {
	      Q[(rw+1)*(t)]+= omegaX[i][t+2]*delta[t+2]*exp(sumg(ncov,xcov,gamma,t+2,scov)+alpha[i]+beta[t+2]);
	    }	
	}
    }
  
  if (mode==2)
    {
      for (int i=1; i<=I; i++)
	{
	  for (int t=2; t<=n; t++)
	    {
	      Q[(rw+1)*(t-2)]+= omega[i][t]*Z[i][t-1]*exp(sumg(ncov,xcov,gamma,t,scov)+alpha[i]+beta[t]);
	    }
	}	
    }
  
  return;
}


void erzeuge_b_Q_2(double* my, double* Q, const DoubleVector& alpha, DoubleVector& beta, DoubleVector& gamma, DoubleVector& delta, const LongMatrix& X,
		   int n, int I, double taubeta, int rw, const DoubleMatrix& xcov, int ncov, int scov, const DoubleMatrix& omega)
{

  /* b-vektor des Proposals*/
  for (int t=0;t<=n; t++)
    {
      my[t]=0.0;
      for (int i=1; i<=I; i++)
	{
	  my[t]+=X[i][t+2];
	  my[t]-=(1-beta[t])*omega[i][t+2]*delta[t+2]*exp(sumg(ncov,xcov,gamma,t+2,scov)+alpha[i]+beta[t]);
	}
    }


  /* Praezisionsmatrix*/
  berechneQ(Q, rw, taubeta, n+1, 1, 0.0);

  for (int i=1; i<=I; i++)
    {
      for (int t=0; t<=n; t++)
	{
	  Q[(rw+1)*(t)]+= omega[i][t+2]*delta[t+2]*exp(sumg(ncov,xcov,gamma,t+2,scov)+alpha[i]+beta[t]);
	}	
    }



  return;
}



void machnu(DoubleVector& mu, const DoubleVector& alpha, DoubleVector& beta, DoubleVector& delta, DoubleMatrix& nu, int I, int n, int ncov, const DoubleMatrix& xcov, int scov)
{
  for (int i=1; i<=I; i++)
    { 
      for (int t=2; t<=n; t++)
	{
	  nu[i][t]=delta[t]*exp(sumg(ncov,xcov,mu,t,scov)+alpha[i]+beta[t]);
	}
    }
  return;
}


void update_gamma_j(int j, const DoubleVector& alpha, DoubleVector& beta, DoubleVector& gamma, DoubleVector& delta, int ncov, const DoubleMatrix& xcov, const LongMatrix& X, int n, int I, double taugamma, DoubleVector& gammaneu, long& acc_gamma, const DoubleMatrix& omega, int scov)
{
  double g = 0;
  double gd = 0;
  double gdd = 0;
  double c = 0;
  for(int i=1;i<=I;i++){
    for(int t=2;t<=n;t++){
      g -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gamma,t,scov)); /* g ist g(gamma[j]^0)*/
      gd -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gamma,t,scov))*xcov[j][t];
      gdd -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gamma,t,scov))*xcov[j][t]*xcov[j][t];
      c += xcov[j][t]*X[i][t];
    }
  }

  double s = sqrt(1/(taugamma-gdd)); /* s ist s*/
  double b = c + gd - gamma[j]*gdd;
  double m = b*s*s;

  double gammajStar = gsl_ran_gaussian(s) + m;

  /* Debug stuff deleted */

  for(int k=0;k<ncov;k++){
    gammaneu[k] = gamma[k];
  }

  gammaneu[j] = gammajStar;

  double g2 = 0;
  double gd2 = 0;
  double gdd2 = 0;
  for(int i=1;i<=I;i++){
    for(int t=2;t<=n;t++){
      g2 -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gammaneu,t,scov)); /* g2 ist g(gamma[j])*/
      gd2 -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gammaneu,t,scov))*xcov[j][t];
      gdd2 -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gammaneu,t,scov))*xcov[j][t]*xcov[j][t];

    }
  }

  double s2 = sqrt(1/(taugamma-gdd2)); /*s2 ist s^0*/
  double b2 = c + gd2 - gammajStar*gdd2; 
  double m2 = b2*s2*s2;

  double a = 0;

  a += gammajStar*c;
  a -= gamma[j]*c;
  a -= 0.5*taugamma*gammajStar*gammajStar;
  a += 0.5*taugamma*gamma[j]*gamma[j];
  a += g2;
  a -= g;
  a += log(s);
  a -= log(s2);
  a += 0.5*((gammajStar-m)/s)*((gammajStar-m)/s);
  a -= 0.5*((gamma[j]-m2)/s2)*((gamma[j]-m2)/s2);

  if(exp(a)>gsl_rng_uniform()){
    gamma[j] = gammajStar;
    acc_gamma += 1;
  }

  return;
}



void update_beta_t(int t, const DoubleVector& alpha, DoubleVector& beta, DoubleVector& gamma, DoubleVector& delta, int ncov, const DoubleMatrix& xcov, const LongMatrix& X, int n, int I, double taubeta, long& acc_beta, const DoubleMatrix& omega, int scov)
{
  double h = 0;
  double c = 0;
  double d = 0;
  for(int i=1;i<=I;i++){
    h -= omega[i][t]*delta[t]*exp(alpha[i] + beta[t] + sumg(ncov,xcov,gamma,t,scov)); /* h ist h(beta[t]^0), beta ist \beta^0, betatStar ist \beta*/
    c += X[i][t];
  }
  if(t==2){
    c -= taubeta*(beta[t+2]-2*beta[t+1]);
    d = taubeta;
  }
  if(t==3){
    c -= taubeta*((beta[t+2]-2*beta[t+1]) + (-2*beta[t+1] - 2*beta[t-1]));
    d = 5*taubeta;
  }
  if((t>=4)&&(t<=(n-2))){
    c -= taubeta*((beta[t+2]-2*beta[t+1]) + (-2*beta[t+1] - 2*beta[t-1]) + (beta[t-2] - 2*beta[t-1]));
    d = 6*taubeta;
  }
  if(t==(n-1)){
    c -= taubeta*((-2*beta[t+1] - 2*beta[t-1]) + (beta[t-2] - 2*beta[t-1]));
    d = 5*taubeta;
  }
  if(t==n){
    c -= taubeta*(beta[t-2] - 2*beta[t-1]);
    d = taubeta;
  }



  double s = sqrt(1/(d - h)); /* s ist s*/
  double b = c + (1 - beta[t])*h;
  double m = b*s*s;

  double betatStar = gsl_ran_gaussian(s) + m;

  double h2 = 0;
  for(int i=1;i<=I;i++){
    h2 -= omega[i][t]*delta[t]*exp(alpha[i] + betatStar + sumg(ncov,xcov,gamma,t,scov)); /* h2 ist h(beta[t])*/
  }

  double s2 = sqrt(1/(d - h2)); /* s2 ist s^0*/
  double b2 = c + (1 - betatStar)*h2;
  double m2 = b2*s2*s2;


  double a = 0;

  a += betatStar*c;
  a -= beta[t]*c;
  a -= 0.5*d*betatStar*betatStar;
  a += 0.5*d*beta[t]*beta[t];
  a += h2;
  a -= h;
  a += log(s);
  a -= log(s2);
  a += 0.5*((betatStar-m)/s)*((betatStar-m)/s);
  a -= 0.5*((beta[t]-m2)/s2)*((beta[t]-m2)/s2);

  if(exp(a)>gsl_rng_uniform()){
    beta[t] = betatStar;
    acc_beta += 1;
  }

  return;
}






void update_lambda_br(DoubleMatrix& lambda, DoubleMatrix& lambda_br,DoubleVector& xi_lambda, IntMatrix& breakpoints, IntMatrix& breakpointsStar, IntVector& K, IntVector& KStar, IntVector& Km1, double alpha_lambda, double beta_lambda, const LongMatrix& Y, const LongMatrix& Z, int n, int I, double& acceptedbr, const DoubleMatrix& omega, int theta_pred_estim, int xi_estim, int K_geom, double p_K, double alpha_xi, double beta_xi)
{
  /*update breakpoints of lambda using reversible jump MCMC*/


  int newbreakpoint =0;
  int removebreakpoint=0;
  int newbreakpointnumber=0;
  int u;
  double v=1;
  double a;
  double alpha_la;
  double beta_la;
  for(int i=1;i<=I;i++){
    if(!theta_pred_estim){
      a=gsl_rng_uniform();
      if(a<0.5){u=1;}else{u=2;}
 
      if(K[i]==1){u=2;v=.5;}  /*K[i] is number of segments of lambda*/
      if(K[i]==(n-1)){u=1;v=.5;} /*if(!theta_pred_estim) max of K[i] is n-1*/


	
      /*decide if new brreakpoint or remove breakpoint*/
      if(u==1){/*remove breakpoint*/
	if(K[i]==2){v=2;} 
	KStar[i]=K[i]-1;
	a=gsl_rng_uniform();
	removebreakpoint=(int)floor(a*(double)(K[i]-1))+1; 
	/*generate breakpointsStar*/
	for(int k=1;k<removebreakpoint;k++){
	  breakpointsStar[i][k]=breakpoints[i][k];
	}
	for(int k=(removebreakpoint+1);k<=K[i];k++){
	  breakpointsStar[i][k-1]=breakpoints[i][k];
	}
      }/*if(u==1)*/
      if(u==2){/*new breakpoint*/
	if(K[i]==(n-2)){v=2;} 
	KStar[i]=K[i]+1;
	/*sample new breakpoint*/
	int need=1;
	while(need==1){
	  need=0;
	  a=gsl_rng_uniform();
	  newbreakpoint=int(floor(a*double(n)))+1; 
	  if (newbreakpoint<=2){need=1;}
	  if (newbreakpoint>n){need=1;}
	  for(int k=1;k<=K[i];k++){
	    if(newbreakpoint==breakpoints[i][k]){
	      need=1;
	    }
	  }
	}/*while(need==1)*/
	/*generate breakpointsStar*/
	for(int k=1;k<=K[i];k++){
	  if((newbreakpoint>breakpoints[i][k-1])&&(newbreakpoint<breakpoints[i][k])){
	    newbreakpointnumber = k;
	  }
	}
	for(int k=1;k<newbreakpointnumber;k++){
	  breakpointsStar[i][k]=breakpoints[i][k];
	  need=0;

	}
	breakpointsStar[i][newbreakpointnumber]=newbreakpoint;
	for(int k=newbreakpointnumber;k<=K[i];k++){ /*changed K[I]*/
	  breakpointsStar[i][k+1]=breakpoints[i][k];
	}
      }/*if(u==2)*/
    }/*if(!theta_pred_estim)*/

    if(theta_pred_estim){

      a=gsl_rng_uniform();
      if(a<0.5){u=1;}else{u=2;}
 
      if(K[i]==1){u=2;v=.5;}
      if(K[i]==(n)){u=1;v=.5;} 


	
      /*decide if new brreakpoint or remove breakpoint*/
      if(u==1){/*remove breakpoint*/
	if(K[i]==2){v=2;} 
	KStar[i]=K[i]-1;
	a=gsl_rng_uniform();
	removebreakpoint=(int)floor(a*(double)(K[i]-1))+1; 
	/*generate breakpointsStar*/
	for(int k=1;k<removebreakpoint;k++){
	  breakpointsStar[i][k]=breakpoints[i][k];
	}
	for(int k=(removebreakpoint+1);k<=K[i];k++){ /*changed K[I]*/
	  breakpointsStar[i][k-1]=breakpoints[i][k];
	}
      }/*if(u==1)*/
      if(u==2){/*new breakpoint*/
	if(K[i]==(n-1)){v=2;} 
	KStar[i]=K[i]+1;
	/*sample new breakpoint*/
	int need=1;
	while(need==1){
	  need=0;
	  a=gsl_rng_uniform();
	  newbreakpoint=int(floor(a*double(n+1)))+1; 
	  if (newbreakpoint<=2){need=1;}
	  if (newbreakpoint>(n+1)){need=1;}
	  for(int k=1;k<=K[i];k++){
	    if(newbreakpoint==breakpoints[i][k]){
	      need=1;
	    }
	  }
	}/*while(need==1)*/
	/*generate breakpointsStar*/
	for(int k=1;k<=K[i];k++){
	  if((newbreakpoint>breakpoints[i][k-1])&&(newbreakpoint<breakpoints[i][k])){
	    newbreakpointnumber = k;
	  }
	}
	for(int k=1;k<newbreakpointnumber;k++){
	  breakpointsStar[i][k]=breakpoints[i][k];
	  need=0;

	}
	breakpointsStar[i][newbreakpointnumber]=newbreakpoint;
	for(int k=newbreakpointnumber;k<=K[i];k++){ 
	  breakpointsStar[i][k+1]=breakpoints[i][k];
	}
      }/*if(u==2)*/
    }/*if(theta_pred_estim)    */
      

    double sumY1=0.0;
    double sumoZ1=0.0;
    double sumY2=0.0;
    double sumoZ2=0.0;
    double sumY3=0.0;
    double sumoZ3=0.0;
    double sumY4=0.0;
    double sumoZ4=0.0;
    if(newbreakpointnumber!=(KStar[i]-1)){
      if (u==2)
	{
	  for (int t=breakpointsStar[i][newbreakpointnumber-1]; t<breakpointsStar[i][newbreakpointnumber];t++)
	    {
	      sumY1+=Y[i][t];
	      sumoZ1+=(omega[i][t]*Z[i][t-1]);
	    }
	  for (int t=breakpointsStar[i][newbreakpointnumber]; t<breakpointsStar[i][newbreakpointnumber+1];t++)
	    {
	      sumY2+=Y[i][t];
	      sumoZ2+=(omega[i][t]*Z[i][t-1]);
	    }
	  sumY3=sumY1+sumY2;
	  sumoZ3=sumoZ1+sumoZ2;
	}
    }
    if(removebreakpoint!=(KStar[i])){
      if (u==1)
	{
	  for (int t=breakpoints[i][removebreakpoint-1]; t<breakpoints[i][removebreakpoint];t++)
	    {
	      sumY3+=Y[i][t];
	      sumoZ3+=(omega[i][t]*Z[i][t-1]);
	    }
	  for (int t=breakpoints[i][removebreakpoint]; t<breakpoints[i][removebreakpoint+1];t++)
	    {
	      sumY4+=Y[i][t];
	      sumoZ4+=(omega[i][t]*Z[i][t-1]);
	    }
	  sumY1=sumY3+sumY4;
	  sumoZ1=sumoZ3+sumoZ4;
	}
    }
    if(newbreakpointnumber==(KStar[i]-1)){
      if (u==2){
	for (int t=breakpointsStar[i][newbreakpointnumber-1]; t<breakpointsStar[i][newbreakpointnumber];t++)
	  {
	    sumY1+=Y[i][t];
	    sumoZ1+=(omega[i][t]*Z[i][t-1]);
	  }
        if(breakpoints[i][newbreakpointnumber]!=(n+1)){
	  for (int t=breakpointsStar[i][newbreakpointnumber]; t<(breakpointsStar[i][newbreakpointnumber+1]-1);t++)
	    {
	      sumY2+=Y[i][t];
	      sumoZ2+=(omega[i][t]*Z[i][t-1]);
	    }
	}
	sumY3=sumY1+sumY2;
	sumoZ3=sumoZ1+sumoZ2;
      }
    }
    if(removebreakpoint==(KStar[i])){
      if (u==1){
	for (int t=breakpoints[i][removebreakpoint-1]; t<breakpoints[i][removebreakpoint];t++)
	  {
	    sumY3+=Y[i][t];
	    sumoZ3+=(omega[i][t]*Z[i][t-1]);
	  }
        if(breakpoints[i][removebreakpoint]!=(n+1)){
	  for (int t=breakpoints[i][removebreakpoint]; t<(breakpoints[i][removebreakpoint+1]-1);t++)
	    {
	      sumY4+=Y[i][t];
	      sumoZ4+=(omega[i][t]*Z[i][t-1]);
	    }
        }
	sumY1=sumY3+sumY4;
	sumoZ1=sumoZ3+sumoZ4;

      }
    }

    alpha_la = alpha_lambda;
    beta_la = beta_lambda;
    if(xi_estim){
      beta_la = xi_lambda[i];
    }
    sumY1+=alpha_la;
    sumY2+=alpha_la;
    sumY3+=alpha_la;
    sumY4+=alpha_la;
    sumoZ1+=beta_la;
    sumoZ2+=beta_la;
    sumoZ3+=beta_la;
    sumoZ4+=beta_la;

	
    double accbr = alpha_la*log(beta_la)-gsl_sf_lngamma(alpha_la); /* alpha_lambda beta_lambda*/
    if(K_geom){
      accbr += log(1-p_K);
    }
    accbr=accbr*pow((double)(-1.0),u)+log(v);

    if (u==2)
      {
	accbr=accbr+gsl_sf_lngamma(sumY1)+gsl_sf_lngamma(sumY2)-gsl_sf_lngamma(sumY3);
        accbr=accbr-sumY1*log(sumoZ1)-sumY2*log(sumoZ2)+sumY3*log(sumoZ3);
      }
    if (u==1)
      {
	accbr=accbr+gsl_sf_lngamma(sumY1)-gsl_sf_lngamma(sumY3)-gsl_sf_lngamma(sumY4);
	accbr=accbr-sumY1*log(sumoZ1)+sumY3*log(sumoZ3)+sumY4*log(sumoZ4);
      }
    if (gsl_rng_uniform() <= exp(accbr)) {
      /* for(int i=1;i<=I;i++){ changed*/
      K[i] = KStar[i]; 
      for(int k=0;k<=KStar[i];k++){
	breakpoints[i][k]=breakpointsStar[i][k];
      }
      acceptedbr++;
      /*} changed*/
    }

    Km1[i]=K[i]-1;


    /*update xi_lambda*/
    if(xi_estim){
      double a = alpha_xi + alpha_lambda*K[i];
      double b = beta_xi;
      for(int k=1;k<=K[i];k++){
	b += lambda_br[i][k];
      }
      xi_lambda[i] = gsl_ran_gamma (a, 1/b);
    }else{
      xi_lambda[i] = beta_lambda;
    }
  } //for i

  //if theta_pred_estim
  //update timevariing lambda with breakpoints via reversible jump MCMC
  //breakpoints = 2,...,n+1, number of breakpoins in region i is K[i]+1, number of time regions is K[i] (breakpoints[i][0]=2, ..., breakpoints[i][K[i]]=n+1)
  //if !theta_pred_estim
  //update timevariing lambda with breakpoints via reversible jump MCMC
  //breakpoints = 2,...,n+2, number of breakpoins in region i is K[i]+1, number of time regions is K[i] (breakpoints[i][0]=2, ..., breakpoints[i][K[i]]=n+2)

  for(int i=1;i<=I;i++){
    for(int k=1;k<=(K[i]-1);k++){
      double a = alpha_la;
      double b = xi_lambda[i];
      for(int t=breakpoints[i][k-1];t<breakpoints[i][k];t++){ //breakpoints = 2,...,n+2
	a += Y[i][t];
	b += omega[i][t]*Z[i][t-1];
      }
      lambda_br[i][k] = gsl_ran_gamma (a, 1/b);
      for(int t=breakpoints[i][k-1];t<breakpoints[i][k];t++){
	lambda[i][t] = lambda_br[i][k]; //calculate the lambda[i][t] from the lambda_br[i][k]
      }
    }


    double a = alpha_la;
    double b = xi_lambda[i];
    if(breakpoints[i][K[i]-1]<(n+1)){ 
      for(int t=breakpoints[i][K[i]-1];t<(breakpoints[i][K[i]]-1);t++){ //breakpoints = 2,...,n+2
	a += Y[i][t];
	b += omega[i][t]*Z[i][t-1];
      }
    }
    lambda_br[i][K[i]] = gsl_ran_gamma (a, 1/b);
    for(int t=breakpoints[i][K[i]-1];t<breakpoints[i][K[i]];t++){
      lambda[i][t] = lambda_br[i][K[i]]; //calculate the lambda[i][t] from the lambda_br[i][k]
    }
  } //for i
  return;
}







void update_delta_br(DoubleVector& delta, DoubleVector& delta_br,double &xi_delta, IntVector& breakpoints_delta, IntVector& breakpointsStar_delta, int& K_delta, int& KStar_delta, int& Km1_delta, double delta_a, double delta_b, const LongMatrix& X, const DoubleMatrix& nu, int n, int I, double& acceptedbr_delta, const DoubleMatrix& omega, int xi_estim_delta, int K_geom, double p_K, double alpha_xi, double beta_xi)
{
  //update breakpoints of lambda using reversible jump MCMC


  int newbreakpoint = 0;
  int removebreakpoint = 0;
  int newbreakpointnumber = 0;
  int u;
  double v=1;
  double a;
  double alpha_de;
  double beta_de;

  a=gsl_rng_uniform();
  if(a<0.5){u=1;}else{u=2;}
 
  if(K_delta==1){u=2;v=.5;}
  if(K_delta==(n-1)){u=1;v=.5;} 


	
  //decide if new brreakpoint or remove breakpoint
  if(u==1){//remove breakpoint
    if(K_delta==2){v=2;} 
    KStar_delta=K_delta-1;
    a=gsl_rng_uniform();
    removebreakpoint=(int)floor(a*(double)(K_delta-1))+1; 
    //generate breakpointsStar_delta
    for(int k=1;k<removebreakpoint;k++){
      breakpointsStar_delta[k]=breakpoints_delta[k];
    }
    for(int k=(removebreakpoint+1);k<=K_delta;k++){
      breakpointsStar_delta[k-1]=breakpoints_delta[k];
    }
  }//if(u==1)



  if(u==2){//new breakpoint
    if(K_delta==(n-2)){v=2;} 
    KStar_delta=K_delta+1;
    //sample new breakpoint
    int need=1;
    while(need==1){
      need=0;
      a=gsl_rng_uniform();
      newbreakpoint=int(floor(a*double(n)))+1; 
      if (newbreakpoint<=2){need=1;}
      if (newbreakpoint>n){need=1;}
      for(int k=1;k<=K_delta;k++){
	if(newbreakpoint==breakpoints_delta[k]){
	  need=1;
	}
      }
    }//while(need==1)
    //generate breakpointsStar_delta
    for(int k=1;k<=K_delta;k++){
      if((newbreakpoint>breakpoints_delta[k-1])&&(newbreakpoint<breakpoints_delta[k])){
	newbreakpointnumber = k;
      }
    }

    for(int k=1;k<newbreakpointnumber;k++){
      breakpointsStar_delta[k]=breakpoints_delta[k];
      need=0;

    }
    breakpointsStar_delta[newbreakpointnumber]=newbreakpoint;

    for(int k=newbreakpointnumber;k<=K_delta;k++){ //changed K
      breakpointsStar_delta[k+1]=breakpoints_delta[k];
    }
  }//if(u==2)


  double sumX1=0.0;
  double sumoZ1=0.0;
  double sumX2=0.0;
  double sumoZ2=0.0;
  double sumX3=0.0;
  double sumoZ3=0.0;
  double sumX4=0.0;
  double sumoZ4=0.0;
  if(newbreakpointnumber!=(KStar_delta-1)){
    if (u==2)
      {
	for(int i=1;i<=I;i++){
	  for (int t=breakpointsStar_delta[newbreakpointnumber-1]; t<breakpointsStar_delta[newbreakpointnumber];t++)
	    {
	      sumX1+=X[i][t];
	      sumoZ1+=(omega[i][t]*nu[i][t]/delta[t]);
	    }
	  for (int t=breakpointsStar_delta[newbreakpointnumber]; t<breakpointsStar_delta[newbreakpointnumber+1];t++)
	    {
	      sumX2+=X[i][t];
	      sumoZ2+=(omega[i][t]*nu[i][t]/delta[t]);
	    }
	}
	sumX3=sumX1+sumX2;
	sumoZ3=sumoZ1+sumoZ2;
      }
  }
  if(removebreakpoint!=(KStar_delta)){
    if (u==1)
      {
	for(int i=1;i<=I;i++){
	  for (int t=breakpoints_delta[removebreakpoint-1]; t<breakpoints_delta[removebreakpoint];t++)
	    {
	      sumX3+=X[i][t];
	      sumoZ3+=(omega[i][t]*nu[i][t]/delta[t]);
	    }
	  for (int t=breakpoints_delta[removebreakpoint]; t<breakpoints_delta[removebreakpoint+1];t++)
	    {
	      sumX4+=X[i][t];
	      sumoZ4+=(omega[i][t]*nu[i][t]/delta[t]);
	    }
	}
	sumX1=sumX3+sumX4;
	sumoZ1=sumoZ3+sumoZ4;
      }
  }
  if(newbreakpointnumber==(KStar_delta-1)){
    if (u==2){
      for(int i=1;i<=I;i++){
	for (int t=breakpointsStar_delta[newbreakpointnumber-1]; t<breakpointsStar_delta[newbreakpointnumber];t++)
	  {
	    sumX1+=X[i][t];
	    sumoZ1+=(omega[i][t]*nu[i][t]/delta[t]);
	  }
        if(breakpoints_delta[newbreakpointnumber]!=(n+1)){
	  for (int t=breakpointsStar_delta[newbreakpointnumber]; t<(breakpointsStar_delta[newbreakpointnumber+1]-1);t++)
	    {
	      sumX2+=X[i][t];
	      sumoZ2+=(omega[i][t]*nu[i][t]/delta[t]);
	    }
	}
      }
      sumX3=sumX1+sumX2;
      sumoZ3=sumoZ1+sumoZ2;
    }
  }
  if(removebreakpoint==(KStar_delta)){
    if (u==1){
      for(int i=1;i<=I;i++){
	for (int t=breakpoints_delta[removebreakpoint-1]; t<breakpoints_delta[removebreakpoint];t++)
	  {
	    sumX3+=X[i][t];
	    sumoZ3+=(omega[i][t]*nu[i][t]/delta[t]);
	  }
        if(breakpoints_delta[removebreakpoint]!=(n+1)){
	  for (int t=breakpoints_delta[removebreakpoint]; t<(breakpoints_delta[removebreakpoint+1]-1);t++)
	    {
	      sumX4+=X[i][t];
	      sumoZ4+=(omega[i][t]*nu[i][t]/delta[t]);
	    }
        }
      }
      sumX1=sumX3+sumX4;
      sumoZ1=sumoZ3+sumoZ4;

    }
  }

  // cout << sumX1 << sumX2 << sumX3 << sumX4 << sumoZ1 << sumoZ2 << sumoZ3 << sumoZ4 << endl;




  alpha_de = delta_a;
  beta_de = delta_b;
  if(xi_estim_delta){
    beta_de = xi_delta;

  }
  sumX1+=alpha_de;
  sumX2+=alpha_de;
  sumX3+=alpha_de;
  sumX4+=alpha_de;
  sumoZ1+=beta_de;
  sumoZ2+=beta_de;
  sumoZ3+=beta_de;
  sumoZ4+=beta_de;

	
  double accbr = alpha_de*log(beta_de)-gsl_sf_lngamma(alpha_de); // delta_a delta_b
  if(K_geom){
    accbr += log(1-p_K);
  }
  accbr=accbr*pow((double)(-1.0),u)+log(v);

  if (u==2)
    {
      accbr=accbr+gsl_sf_lngamma(sumX1)+gsl_sf_lngamma(sumX2)-gsl_sf_lngamma(sumX3);
      accbr=accbr-sumX1*log(sumoZ1)-sumX2*log(sumoZ2)+sumX3*log(sumoZ3);
    }
  if (u==1)
    {
      accbr=accbr+gsl_sf_lngamma(sumX1)-gsl_sf_lngamma(sumX3)-gsl_sf_lngamma(sumX4);
      accbr=accbr-sumX1*log(sumoZ1)+sumX3*log(sumoZ3)+sumX4*log(sumoZ4);
    }
  if (gsl_rng_uniform() <= exp(accbr)) {
    K_delta = KStar_delta; 
    for(int k=0;k<=KStar_delta;k++){
      breakpoints_delta[k]=breakpointsStar_delta[k];
    }
    acceptedbr_delta++;

  }


  Km1_delta=K_delta-1;

  // cout << "accbr" << accbr << endl;




  //update xi_delta
  if(xi_estim_delta){
    double a = alpha_xi + delta_a*K_delta;
    double b = beta_xi;
    for(int k=1;k<=K_delta;k++){
      b += delta_br[k];
    }
    xi_delta = gsl_ran_gamma (a, 1/b);
  }else{
    xi_delta = delta_b;
  }

  //  cout << "xi_delta  " << xi_delta << endl;
  // cout << "delta_br[1]  " << delta_br[1] << endl;






  //update timevariing delta with breakpoints_delta via reversible jump MCMC


  for(int k=1;k<=(K_delta-1);k++){
    double a = alpha_de;
    double b = xi_delta;
    for(int i=1;i<=I;i++){
      for(int t=breakpoints_delta[k-1];t<breakpoints_delta[k];t++){ //breakpoints_delta = 2,...,n+2
	a += X[i][t];
	b += omega[i][t]*nu[i][t]/delta[t];
      }
    }
    delta_br[k] = gsl_ran_gamma (a, 1/b);
    for(int t=breakpoints_delta[k-1];t<breakpoints_delta[k];t++){
      delta[t] = delta_br[k]; //calculate the delta[t] from the delta_br[k]
    }

  }


  //cout << "nu[1][2]   " << nu[1][2] << endl;
  //cout << "delta[2]   " << delta[2] << endl;

  a = alpha_de;
  double b = xi_delta;
  //	cout << a << "           " << b << endl;
  if(breakpoints_delta[K_delta-1]<(n+1)){ 
    for(int i=1;i<=I;i++){
      for(int t=breakpoints_delta[K_delta-1];t<(breakpoints_delta[K_delta]-1);t++){ //breakpoints_delta = 2,...,n+2
	a += X[i][t];
	b += omega[i][t]*nu[i][t]/delta[t];
      }
    }
  }
  delta_br[K_delta] = gsl_ran_gamma (a, 1/b);
  for(int t=breakpoints_delta[K_delta-1];t<breakpoints_delta[K_delta];t++){
    delta[t] = delta_br[K_delta]; //calculate the delta[t] from the delta_br[k]
  }
  //	cout << a << "           " << b << endl;

  // 	  cout << delta_br[1] << endl;





  return;
}



void update_epsilon_br(DoubleVector& epsilon, DoubleVector& epsilon_br,double& xi_epsilon, IntVector& breakpoints_epsilon, IntVector& breakpointsStar_epsilon, int& K_epsilon, int& KStar_epsilon, int& Km1_epsilon, double epsilon_a, double epsilon_b, const LongMatrix& S, int n, int I, double& acceptedbr_epsilon, const DoubleMatrix& omega, int xi_estim_epsilon, int K_geom, double p_K, double alpha_xi, double beta_xi)
{
  /*update breakpoints of lambda using reversible jump MCMC*/
  
  /*define and initialize to avoid compiler warnings */
  int newbreakpoint = 0;
  int removebreakpoint = 0;
  int newbreakpointnumber = 0;
  int u;
  double v=1;
  double a;
  double alpha_ep;
  double beta_ep;

  a=gsl_rng_uniform();
  if(a<0.5){u=1;}else{u=2;}
 
  if(K_epsilon==1){u=2;v=.5;}
  if(K_epsilon==(n-1)){u=1;v=.5;} 


	
  //decide if new brreakpoint or remove breakpoint
  if(u==1){//remove breakpoint
    if(K_epsilon==2){v=2;} 
    KStar_epsilon=K_epsilon-1;
    a=gsl_rng_uniform();
    removebreakpoint=(int)floor(a*(double)(K_epsilon-1))+1; 
    //generate breakpointsStar_epsilon
    for(int k=1;k<removebreakpoint;k++){
      breakpointsStar_epsilon[k]=breakpoints_epsilon[k];
    }
    for(int k=(removebreakpoint+1);k<=K_epsilon;k++){
      breakpointsStar_epsilon[k-1]=breakpoints_epsilon[k];
    }
  }//if(u==1)



  if(u==2){//new breakpoint
    if(K_epsilon==(n-2)){v=2;} 
    KStar_epsilon=K_epsilon+1;
    //sample new breakpoint
    int need=1;
    while(need==1){
      need=0;
      a=gsl_rng_uniform();
      newbreakpoint=int(floor(a*double(n)))+1; 
      if (newbreakpoint<=2){need=1;}
      if (newbreakpoint>n){need=1;}
      for(int k=1;k<=K_epsilon;k++){
	if(newbreakpoint==breakpoints_epsilon[k]){
	  need=1;
	}
      }
    }//while(need==1)
    //generate breakpointsStar_epsilon
    for(int k=1;k<=K_epsilon;k++){
      if((newbreakpoint>breakpoints_epsilon[k-1])&&(newbreakpoint<breakpoints_epsilon[k])){
	newbreakpointnumber = k;
      }
    }

    for(int k=1;k<newbreakpointnumber;k++){
      breakpointsStar_epsilon[k]=breakpoints_epsilon[k];
      need=0;

    }
    breakpointsStar_epsilon[newbreakpointnumber]=newbreakpoint;

    for(int k=newbreakpointnumber;k<=K_epsilon;k++){ //changed K
      breakpointsStar_epsilon[k+1]=breakpoints_epsilon[k];
    }
  }//if(u==2)


  double sumX1=0.0;
  double sumoZ1=0.0;
  double sumX2=0.0;
  double sumoZ2=0.0;
  double sumX3=0.0;
  double sumoZ3=0.0;
  double sumX4=0.0;
  double sumoZ4=0.0;
  if(newbreakpointnumber!=(KStar_epsilon-1)){
    if (u==2)
      {
	for(int i=1;i<=I;i++){
	  for (int t=breakpointsStar_epsilon[newbreakpointnumber-1]; t<breakpointsStar_epsilon[newbreakpointnumber];t++)
	    {
	      sumX1+=S[i][t];
	      sumoZ1+=(omega[i][t]);
	    }
	  for (int t=breakpointsStar_epsilon[newbreakpointnumber]; t<breakpointsStar_epsilon[newbreakpointnumber+1];t++)
	    {
	      sumX2+=S[i][t];
	      sumoZ2+=(omega[i][t]);
	    }
	}
	sumX3=sumX1+sumX2;
	sumoZ3=sumoZ1+sumoZ2;
      }
  }
  if(removebreakpoint!=(KStar_epsilon)){
    if (u==1)
      {
	for(int i=1;i<=I;i++){
	  for (int t=breakpoints_epsilon[removebreakpoint-1]; t<breakpoints_epsilon[removebreakpoint];t++)
	    {
	      sumX3+=S[i][t];
	      sumoZ3+=(omega[i][t]);
	    }
	  for (int t=breakpoints_epsilon[removebreakpoint]; t<breakpoints_epsilon[removebreakpoint+1];t++)
	    {
	      sumX4+=S[i][t];
	      sumoZ4+=(omega[i][t]);
	    }
	}
	sumX1=sumX3+sumX4;
	sumoZ1=sumoZ3+sumoZ4;
      }
  }
  if(newbreakpointnumber==(KStar_epsilon-1)){
    if (u==2){
      for(int i=1;i<=I;i++){
	for (int t=breakpointsStar_epsilon[newbreakpointnumber-1]; t<breakpointsStar_epsilon[newbreakpointnumber];t++)
	  {
	    sumX1+=S[i][t];
	    sumoZ1+=(omega[i][t]);
	  }
        if(breakpoints_epsilon[newbreakpointnumber]!=(n+1)){
	  for (int t=breakpointsStar_epsilon[newbreakpointnumber]; t<(breakpointsStar_epsilon[newbreakpointnumber+1]-1);t++)
	    {
	      sumX2+=S[i][t];
	      sumoZ2+=(omega[i][t]);
	    }
	}
      }
      sumX3=sumX1+sumX2;
      sumoZ3=sumoZ1+sumoZ2;
    }
  }
  if(removebreakpoint==(KStar_epsilon)){
    if (u==1){
      for(int i=1;i<=I;i++){
	for (int t=breakpoints_epsilon[removebreakpoint-1]; t<breakpoints_epsilon[removebreakpoint];t++)
	  {
	    sumX3+=S[i][t];
	    sumoZ3+=(omega[i][t]);
	  }
        if(breakpoints_epsilon[removebreakpoint]!=(n+1)){
	  for (int t=breakpoints_epsilon[removebreakpoint]; t<(breakpoints_epsilon[removebreakpoint+1]-1);t++)
	    {
	      sumX4+=S[i][t];
	      sumoZ4+=(omega[i][t]);
	    }
        }
      }
      sumX1=sumX3+sumX4;
      sumoZ1=sumoZ3+sumoZ4;

    }
  }

  // cout << sumX1 << sumX2 << sumX3 << sumX4 << sumoZ1 << sumoZ2 << sumoZ3 << sumoZ4 << endl;




  alpha_ep = epsilon_a;
  beta_ep = epsilon_b;
  if(xi_estim_epsilon){
    beta_ep = xi_epsilon;

  }
  sumX1+=alpha_ep;
  sumX2+=alpha_ep;
  sumX3+=alpha_ep;
  sumX4+=alpha_ep;
  sumoZ1+=beta_ep;
  sumoZ2+=beta_ep;
  sumoZ3+=beta_ep;
  sumoZ4+=beta_ep;

  // cout << "xi_epsilon  " << xi_epsilon <<endl;
	
  double accbr = alpha_ep*log(beta_ep)-gsl_sf_lngamma(alpha_ep); // epsilon_a epsilon_b
  if(K_geom){
    accbr += log(1-p_K);
  }
  accbr=accbr*pow((double)(-1.0),u)+log(v);

  if (u==2)
    {
      accbr=accbr+gsl_sf_lngamma(sumX1)+gsl_sf_lngamma(sumX2)-gsl_sf_lngamma(sumX3);
      accbr=accbr-sumX1*log(sumoZ1)-sumX2*log(sumoZ2)+sumX3*log(sumoZ3);
    }
  if (u==1)
    {
      accbr=accbr+gsl_sf_lngamma(sumX1)-gsl_sf_lngamma(sumX3)-gsl_sf_lngamma(sumX4);
      accbr=accbr-sumX1*log(sumoZ1)+sumX3*log(sumoZ3)+sumX4*log(sumoZ4);
    }
  if (gsl_rng_uniform() <= exp(accbr)) {
    K_epsilon = KStar_epsilon; 
    for(int k=0;k<=KStar_epsilon;k++){
      breakpoints_epsilon[k]=breakpointsStar_epsilon[k];
    }
    acceptedbr_epsilon++;

  }


  Km1_epsilon=K_epsilon-1;

  // cout << "accbr" << accbr << endl;




  //update xi_epsilon
  if(xi_estim_epsilon){
    double a = alpha_xi + epsilon_a*K_epsilon;
    double b = beta_xi;
    for(int k=1;k<=K_epsilon;k++){
      b += epsilon_br[k];
    }
    xi_epsilon = gsl_ran_gamma (a, 1/b);
  }else{
    xi_epsilon = epsilon_b;
  }

  //  cout << "xi_epsilon  " << xi_epsilon << endl;
  // cout << "epsilon_br[1]  " << epsilon_br[1] << endl;






  //update timevariing epsilon with breakpoints_epsilon via reversible jump MCMC


  for(int k=1;k<=(K_epsilon-1);k++){
    double a = alpha_ep;
    double b = xi_epsilon;
    for(int i=1;i<=I;i++){
      for(int t=breakpoints_epsilon[k-1];t<breakpoints_epsilon[k];t++){ //breakpoints_epsilon = 2,...,n+2
	a += S[i][t];
	b += omega[i][t];
      }
    }
    epsilon_br[k] = gsl_ran_gamma (a, 1/b);
    for(int t=breakpoints_epsilon[k-1];t<breakpoints_epsilon[k];t++){
      epsilon[t] = epsilon_br[k]; //calculate the epsilon[t] from the epsilon_br[k]
    }

  }


  //cout << "nu[1][2]   " << nu[1][2] << endl;
  //cout << "epsilon[2]   " << epsilon[2] << endl;

  a = alpha_ep;
  double b = xi_epsilon;
  //	cout << a << "           " << b << endl;
  if(breakpoints_epsilon[K_epsilon-1]<(n+1)){ 
    for(int i=1;i<=I;i++){
      for(int t=breakpoints_epsilon[K_epsilon-1];t<(breakpoints_epsilon[K_epsilon]-1);t++){ //breakpoints_epsilon = 2,...,n+2
	a += S[i][t];
	b += omega[i][t];
      }
    }
  }
  epsilon_br[K_epsilon] = gsl_ran_gamma (a, 1/b);
  for(int t=breakpoints_epsilon[K_epsilon-1];t<breakpoints_epsilon[K_epsilon];t++){
    epsilon[t] = epsilon_br[K_epsilon]; //calculate the epsilon[t] from the epsilon_br[k]
  }
  //	cout << a << "           " << b << endl;

  // 	  cout << epsilon_br[1] << endl;





  return;
}









/***********************************************************************
 Convert surveillance data object to a format suitable for the
 twin software (which has a notation similar to the one in the paper)
 That is: nothing is at index zero (neither column nor row)

 Params:
  x - the vector with the data (can become a matrix later on)
  n - length of the vector
  I - number of columns in the data (currently this is always zero)

 Returns:
  matrix with the data
*/
LongMatrix surveillancedata2twin(int* x, int n, int I) {

  //Allocate data structure 
  LongMatrix Z(I+1, n+1);

  /* Fill with zeros at all 0 index (rows & columns) */
  for (int t=0; t<=n; t++){ Z[0][t]=0; }
  for (int i=0; i<=I; i++){ Z[i][0]=0; }

  //Start @ index 1. (Z[0] is not defined)
  for (int t=1; t<=n; t++) {
    for (int i=1; i<=I; i++) { 
      Z[i][t]=x[t-1]; 
    } 
  }
  /* Done */
  return(Z);
}



/*********************************************************************
 * Read integer data from file. The format is: n Z_1 Z_2 ... Z_n
 *********************************************************************/

// this could also be changed to return LongMatrix, but it is not used...
// long** readData(char* fileName, long *size, long *size2) {
//   //Read Data from file and check if everything went allright

//   //cout << "\"" << fileName  << "\"" <<  endl;

//   ifstream fin(fileName);
  
//   if (! fin) {
//     cerr << "Error: File \"" << fileName << "\" not found!" << endl;
//     return(NULL);
//   }
  
//   // cout << fin << endl;
//   //Count the number of data entries
//   int n=0;
//   fin >> n;
//   Rprintf("n=%d\n",n);
//   int I=1; 
//   //fin >> I;
//   //cout << "I=" << I << endl;

//   long **Z = new long*[I+1];
//   for (long i=0; i<=I; i++){
//     Z[i] = new long[n+1];
//   }

	
//   for (long t=0; t<=n; t++){
//     Z[0][t]=0;
//   }
	
//   for (long i=0; i<=I; i++){
//     Z[i][0]=0;
//   }
	
//   //Start @ index 1. (Z[0] is not defined)
//   int t=1;
//   while (!fin.eof() && (t<=n)) { 	
//     int i=1;
//     while (!fin.eof() && (i<=I)) { 
//       fin >> Z[i][t];
//       i++;
//     }
//     t++;
//   }

//   fin.close();

//   //Return the result consisting of Z and n
//   *size = n;
//   *size2 = I;

//   return(Z);
// }


/* Calculate the deviance of the data we use that the data, Z, is a
 *  sum of Poisson distributed variables, i.e. it is Poisson
 *  distributed.
 *
 *    Z_t = S_t + X_t + Y_t, i.e. 
 *    Z_t ~ Po(nu*p + nu*(1-p) + lambda*W_{t-1})
 *
 *    D = -2log p(Z|theta) + 2 log p(Z|\mu(theta)=Z)
 */
double satdevalt(int n, int I, const LongMatrix& X, const LongMatrix& Y, const LongMatrix& Z, 
		 const DoubleMatrix& omega, const DoubleMatrix& lambda, const DoubleMatrix& nu, double *xi, DoubleMatrix& eta, DoubleMatrix& eta2, DoubleMatrix& varr, double psi, int overdispersion) {
  double res = 0;
  //Loop over all data
  for (int i=1; i<=I; i++) {
    for (int t=2; t<=n; t++) {
      //Use the equation derived for the saturated deviance in the paper
      //calculate the mean and variance of Z[i][t]
      eta[i][t] = (nu[i][t]*xi[i]+lambda[i][t]*Z[i][t-1]);
      eta2[i][t] = eta[i][t];
      if(overdispersion){
        varr[i][t] = eta2[i][t]*(1+eta2[i][t]/psi);
      }else{
	varr[i][t] = eta2[i][t];
      }
      //calculate the Deviance in the Poisson and NegBin case
      if(!overdispersion){
        if (Z[i][t] == 0) {
          res += 2 * eta[i][t];
        } else {
          res += 2 * ( Z[i][t] * log(Z[i][t]/eta[i][t]) - Z[i][t] + eta[i][t]);
        }
      }
      if(overdispersion){
        if (Z[i][t] == 0) {
          res += 2 * ( - (Z[i][t]+psi) * log((Z[i][t]+psi)/(eta[i][t]+psi)));
        } else {
          res += 2 * ( - (Z[i][t]+psi) * log((Z[i][t]+psi)/(eta[i][t]+psi)) + Z[i][t] * log(Z[i][t]/eta[i][t]));
        }
      }
    }
  }
  return(res);
}


/* Calculate the deviance of the data we use that the data, Z, is a
 *  sum of Poisson distributed variables, i.e. it is Poisson
 *  distributed.
 *
 *    Z_t = X_t + Y_t, i.e. 
 *    Z_t ~ Po(nu_t + lambda_t*Z_{t-1})
 *
 *    D = -2log p(Z|theta)
 */
double satdev(int n, int I, const LongMatrix& Z, 
	      const DoubleMatrix& lambda, const DoubleMatrix& nu, double *xi, DoubleVector& epsilon, DoubleMatrix& eta, double psi, int overdispersion) {
  double res = 0;
  //Loop over all data
  for (int i=1; i<=I; i++) {
    for (int t=2; t<=n; t++) {
      //Use the equation derived for the saturated deviance in the paper
      //calculate the mean and variance of Z[i][t]
      eta[i][t] = (epsilon[t] + nu[i][t]*xi[i]+lambda[i][t]*Z[i][t-1]);
      //calculate the Deviance in the Poisson and NegBin case
      if(!overdispersion){
        res -= 2 * ( Z[i][t] * log(eta[i][t]) - gsl_sf_lngamma(Z[i][t]+1) - eta[i][t]);
      }
      if(overdispersion){
	res -= 2 * ( gsl_sf_lngamma(Z[i][t]+psi) - gsl_sf_lngamma(Z[i][t]+1) - gsl_sf_lngamma(psi) - (Z[i][t]+psi)*log(eta[i][t]+psi) + psi*log(psi) + Z[i][t]*log(eta[i][t]));
      }
    }
  }
  return(res);
}



// Calculate chi square the sum of the qudratic pearson residuals (z-mean)/sd


double chisq(int n, int I, const LongMatrix& Z, 
	     const DoubleMatrix& lambda, const DoubleMatrix& nu, double *xi, DoubleVector& epsilon, DoubleMatrix& eta, DoubleMatrix& varr, DoubleMatrix& rpearson, double psi, int overdispersion) {
  double res = 0;
  //Loop over all data
  for (int i=1; i<=I; i++) {
    for (int t=2; t<=n; t++) {
      //calculate the mean and variance of Z[i][t]
      eta[i][t] = (epsilon[t] + nu[i][t]*xi[i]+lambda[i][t]*Z[i][t-1]);
      if(overdispersion){
        varr[i][t] = eta[i][t]*(1+eta[i][t]/psi);
      }else{
	varr[i][t] = eta[i][t];
      }
      rpearson[i][t] = (Z[i][t]-eta[i][t])/sqrt(varr[i][t]);
      //calculate chisq in the Poisson and NegBin case
      res += rpearson[i][t]*rpearson[i][t];
    }
  }
  return(res);
}
  







/**********************************************************************
 * Estimation in the basic epidemic model
 *
 */
void bplem_estimate(int verbose, ofstream &logfile, ofstream &logfile2, ofstream &acclog, const LongMatrix& Z, double* xi, int n, int I, int T, int nfreq, int burnin, int filter,
		    int samples, int rw) {
  //Model parameters - start values
  double nu_const = alpha_nu/beta_nu;
  double lambda_const = 0.5;
  double psi    = alpha_psi / beta_psi;
  double x      = logit(lambda_const);
  
  if(!verbose) {
    Rprintf("------------------------------------------------\n");
    if (!la_rev){
      Rprintf("lambda:  Ga(%f, %f)-->\t%f\n", alpha_lambda,  beta_lambda, lambda_const);
    }
    if(!varnu){
      Rprintf("nu:      Ga(%f, %f)-->\t%f\n", alpha_nu, beta_nu, nu_const);
    }
    if(overdispersion){
      Rprintf("psi:     Ga(%f, %f)-->\t%f\n", alpha_psi, beta_psi, psi);
    }
    Rprintf("------------------------------------------------\n");
  }
  
  
  //Allocate arrays for all latent variables and initialize them

  // first all 2D arrays (matrices)

  LongMatrix X(I+1, n+1);
  LongMatrix Y(I+1, n+1);
  LongMatrix S(I+1, n+1);
  DoubleMatrix omega(I+1, n+1);
  DoubleMatrix sumX(I+1, n+1);
  DoubleMatrix sumY(I+1, n+1);
  DoubleMatrix sumS(I+1, n+1);
  DoubleMatrix sumomega(I+1, n+1);  
  DoubleMatrix nu(I+1, n+1);  
  DoubleMatrix lambda(I+1, n+2);  
  DoubleMatrix lambda_br(I+1, n+2);  
  DoubleMatrix eta(I+1, n+1);  
  DoubleMatrix eta2(I+1, n+1);  
  DoubleMatrix varr(I+1, n+1);  
  DoubleMatrix rpearson(I+1, n+1);
  DoubleMatrix Sumeta(I+1, n+1);
  DoubleMatrix Sumvarr(I+1, n+1);
  DoubleMatrix Sumrpearson(I+1, n+1);
  IntMatrix breakpoints(I+1, n+2);
  IntMatrix breakpointsStar(I+1, n+2);
  LongMatrix bp(I+1, n+2);

  // long** X = new long*[I+1];
  // long** Y = new long*[I+1];
  // long** S = new long*[I+1];
  // double **omega= new double*[I+1];
  // double** sumX = new double*[I+1];
  // double** sumY = new double*[I+1];
  // double** sumS = new double*[I+1];
  // double **sumomega= new double*[I+1];
  // double **nu= new double*[I+1];
  // double *alpha=new double[I+1];
  // double* beta= new double[n+1];
  // double **lambda=new double*[I+1];
  // double **lambda_br=new double*[I+1];
  // double **eta=new double*[I+1];
  // double **eta2=new double*[I+1];
  // double **varr=new double*[I+1];
  // double **rpearson=new double*[I+1];
  // double **Sumeta=new double*[I+1];
  // double **Sumvarr=new double*[I+1];
  // double **Sumrpearson=new double*[I+1];
  // int **breakpoints=new int*[I+1];
  // int **breakpointsStar=new int*[I+1];
  // long **bp=new long*[I+1];

  // We would have to delete the pointers manually at the end of the routine
  // in order not to corrupt the memory!!!

  // for (long i=0; i<=I; i++){
  //   X[i]=new long[n+1];
  //   Y[i]=new long[n+1];
  //   S[i]=new long[n+1];
  //   omega[i]=new double[n+1];
  //   sumX[i]=new double[n+1];
  //   sumY[i]=new double[n+1];
  //   sumS[i]=new double[n+1];
  //   sumomega[i]=new double[n+1];
  //   nu[i]=new double[n+1];
  //   lambda[i]=new double[n+2];
  //   lambda_br[i]=new double[n+2];
  //   breakpoints[i]=new int[n+2];
  //   breakpointsStar[i]=new int[n+2];
  //   bp[i]=new long[n+2];
  //   eta[i]=new double[n+1];
  //   eta2[i]=new double[n+1];
  //   varr[i]=new double[n+1];
  //   rpearson[i]=new double[n+1];
  //   Sumeta[i]=new double[n+1];
  //   Sumvarr[i]=new double[n+1];
  //   Sumrpearson[i]=new double[n+1];
  // }

  // then the rest (1D arrays and numbers)

  DoubleVector alpha(I + 1);
  DoubleVector beta(n + 1);
  DoubleVector delta(n + 2);
  DoubleVector delta_br(n + 2);
  double xi_delta = 1;
  DoubleVector epsilon(n + 2);
  DoubleVector epsilon_br(n + 2);
  double xi_epsilon = 1;
  double xi_psi = 1;
  
  IntVector K(I + 1);
  IntVector Km1(I + 1);
  IntVector KStar(I + 1);
  DoubleVector xi_lambda(I + 1);
  IntVector breakpoints_delta(n+2);
  IntVector breakpointsStar_delta(n+2);
  LongVector bp_delta(n+2);
  int K_delta = 0; 
  int Km1_delta = 0;
  int KStar_delta = 0;
  IntVector breakpoints_epsilon(n+2);
  IntVector breakpointsStar_epsilon(n+2);
  LongVector bp_epsilon(n+2);
  int K_epsilon = 0;
  int Km1_epsilon = 0;
  int KStar_epsilon = 0;

  LongVector   Xnp1(I + 1);
  LongVector   Snp1(I + 1);
  LongVector   Ynp1(I + 1);
  LongVector   Znp1(I + 1);	
  DoubleVector omeganp1(I + 1);
  DoubleVector nunp1(I + 1);
  
  
  if(!varnu){
    for (int i=0; i<=I; i++) {
      for (int t=0; t<=n; t++) {
        nu[i][t] = alpha_nu/beta_nu;
      }
    }
  }
  
  for (int i=0; i<=I; i++) {
    for (int t=0; t<=n; t++) {
      lambda[i][t] = lambda_const;
    }
  }
  
  for (int i=0; i<=I; i++) {
    for (int t=0; t<=n; t++) {
      X[i][t] = 0; S[i][t] = 0; Y[i][t] = Z[i][t]; omega[i][t] = 1; eta[i][t] = 0; bp[i][t] = 0; bp_delta[t] = 0; bp_epsilon[t] = 0;
      sumX[i][t] = 0;  sumY[i][t] = 0;  sumS[i][t] = 0; sumomega[i][t] = 0; Sumeta[i][t] = 0; Sumrpearson[i][t] = 0;
    }
    bp[i][n+1] = 0; xi_lambda[i] = 1;
    bp_delta[n+1] = 0; bp_epsilon[n+1] = 0;
  }
  
  
  /* Fuer Saisonkomponenente */
  int ncov;
  int scov = 0;
  if(delta_rev){
    scov = 1;
  }

  // determine the number of covariates and allocate then
  // the vectors and design matrix.
  ncov = nu_trend ? (nfreq * 2 + 2) : (nfreq * 2 + 1);
  DoubleVector gamma(ncov);
  DoubleVector gammaneu(ncov);
  DoubleMatrix xcov(ncov, n+2);

  // bad, do not do that:
  // double* gamma;
  // double* gammaneu = NULL;
  // double** xcov;

  if(!nu_trend){
    // ncov=nfreq*2+1;
    // gamma = new double[ncov];
    // gammaneu = new double[ncov];
    // xcov = new double*[ncov];
    // for (int i=0; i<ncov; i++)
    //   {
    // 	xcov[i]=new double[n+2];
    //   }
    if(varnu){
      for (int t=2; t<=(n+1); t++)
	{
	  xcov[0][t]=1.0;
	}
      
      for (int i=1; i<=nfreq; i++)
	{
	  for (int t=2; t<=(n+1); t++)
	    {
	      xcov[i*2-1][t]=sin(2*PI*(t-1)*i/T); //schwingung um einen Zeitpunkt nach hinten verschoben. beginnt bei t=2
	      xcov[i*2][t]=cos(2*PI*(t-1)*i/T);
	      
	    } 
	  // cout << endl;
	}
    } //if varnu
  } //if !nu_trend
  // Saisonkomponente mit linearem trend
  else{
    // ncov=nfreq*2+2;
    // gamma = new double[ncov];
    // gammaneu = new double[ncov];
    // xcov = new double*[ncov];
    // for (int i=0; i<ncov; i++)
    //   {
    // 	xcov[i]=new double[n+2];
    //   }
    if(varnu){
      for (int t=2; t<=(n+1); t++)
	{
	  xcov[0][t]=1.0;
	  xcov[ncov-1][t]=t-n/2;
	}
      
      for (int i=1; i<=nfreq; i++)
	{
	  for (int t=2; t<=(n+1); t++)
	    {
	      xcov[i*2-1][t]=sin(2*PI*(t-1)*i/T); //schwingung um einen Zeitpunkt nach hinten verschoben. beginnt bei t=2
	      xcov[i*2][t]=cos(2*PI*(t-1)*i/T);
	      
	    } 
	  //cout << endl;
	}
    } //if varnu
  } // if nu_trend
  
  

  // Regionenanteil
  DoubleVector xreg(I + 1);

  if(varnu){
    for (int i=1; i<=I; i++)
      {
	xreg[i]=log(xi[i]);
	xi[i]=1;
	//	    cout << xreg[i]<<endl;
      }
  }
  double taualpha=alpha_a/alpha_b;
  double taubeta=beta_a/beta_b;
  /*double taubetaStar;*/
  /*double taubeta=beta_a/beta_b;*/
  double taugamma=gamma_a/gamma_b;
  
  
  //  if(la_rev){
  for(int i=0; i<=I; i++){
    K[i]=1;
    Km1[i]=0;
    KStar[i]=1; 
    breakpoints[i][0]=2;
    breakpoints[i][1]=n+2;
    breakpointsStar[i][0]=2; 
    breakpointsStar[i][1]=n+2; 
    for(int j=0; j<=(n+1); j++){
      lambda_br[i][j]= 0;
    }
    for(int t=2; t<=(n+1); t++){ 
      lambda[i][t]= 0; 
    } 
  }
  
  //  }
  //  if(delta_rev){
  K_delta=1;
  Km1_delta=0;
  KStar_delta=1; 
  breakpoints_delta[0]=2;
  breakpoints_delta[1]=n+2;
  breakpointsStar_delta[0]=2; 
  breakpointsStar_delta[1]=n+2; 
  for(int j=0; j<=(n+1); j++){
    delta_br[j]= 1;
  }
  for(int t=0; t<=(n+1); t++){ 
    delta[t]= 1; 
  } 
  //  }
  

  //  if(epsilon_rev){
  K_epsilon=1;
  Km1_epsilon=0;
  KStar_epsilon=1; 
  breakpoints_epsilon[0]=2;
  breakpoints_epsilon[1]=n+2;
  breakpointsStar_epsilon[0]=2; 
  breakpointsStar_epsilon[1]=n+2; 
  for(int j=0; j<=(n+1); j++){
    epsilon_br[j]= 0;
  }
  for(int t=0; t<=(n+1); t++){ 
    epsilon[t]= 0; 
  } 
  // }

  
  if(varnu){ 
    for (int i=0; i<=I; i++) {
      alpha[i] = log(xreg[i]);
    }
    
    for (int t=0; t<=n; t++) {
      beta[t] = 0.0;
    }
    for (int j=0; j<ncov; j++) {
      gamma[j] = 0.0;
    }
    
    double Zmin = Z[1][2];
    for(int i=1;i<=I;i++){
      for(int t=2;t<=n;t++){
	if(Z[i][t]<Zmin){
	  Zmin = Z[i][t];
	}
      }
    }
    
    // gamma[0] =  log(alpha_nu/beta_nu);
    if(scov==0){
      gamma[0] = log(Zmin+1); //  //
    }
    machnu(gamma, alpha, beta,delta, nu, I, n, ncov, xcov,scov);
    
  }//if varnu  
  
  
  for (int i=1;i<=I; i++) {
    X[i][2] = (long)floor(nu[i][2]);
    Y[i][2] = (long)floor(lambda[i][2]*nu[i][2]/(1 - lambda[i][2]));
    omeganp1[i] = 1;
  }
  
  
	
  //Variables for statistics
  double acceptedPsi = 0;
  double acceptedlambda = 0;
  double acceptedbr = 0;
  double acceptedbr_delta = 0;
  double acceptedbr_epsilon = 0;
  long acc_beta=0;
  long acc_alpha=0;
  long acc_gamma=0;


  /*hoehle: min/max is deprecated double tuneSampleSize = 1000<?burnin; */
  double tuneSampleSize = MIN(1000,burnin);

  double need = 0; 
  double tunex = 0; 
  double tunepsi = 0; 
  double tunetaubeta = 0; 

  
  //Write the header to the logfile
  logfile << "i\t";
  if (!la_rev&&la_estim)
    {
      logfile << "lambda\t";
    }
  logfile << "psi\t";
  logfile << "xipsi\t";

  if(!varnu){
    logfile << "nu\t";
  }

  if(varnu){
    for (int j=0;j<ncov;j++) {
      logfile << "gamma[" << j << "]\t";
    }
    if(delta_rev){
      logfile << "Kdelta\t" << "xidelta\t";
      for (int j=2; j<=n; j++) {
        logfile << "delta[" << j << "]\t";
      }
    }
  }//if varnu

  if(epsilon_rev){
    logfile << "Kepsilon\t" << "xiepsilon\t";
    for (int j=2; j<=n; j++) {
      logfile << "epsilon[" << j << "]\t";
    }
  }


  if (la_rev&&la_estim)
    {
      logfile << "K\t";
      logfile << "xilambda\t";
      for (int t=2; t<=n; t++)
	{
	  logfile << "lambda[" << t <<"]\t";
	}
    }
  logfile << "Znp1\t";
  logfile << "D\t" << endl;


  //Write the header to the logfile2
  for (int t=1;t<=n;t++) {
    logfile2 << "X[" << t << "]\tY[" << t << "]\tomega["
	     << t << "]\tbp[" << t <<"]\t";
  }
  logfile2 << "bp[" << n+1 <<"]\t";
  logfile2 << endl;


  //Calculate the necessary number of samples.
  int sampleSize = filter*samples + burnin;
  if (!verbose) {
    Rprintf("Total number of samples = %d + %d * %d = %d\n", burnin, filter, samples, sampleSize);
    //if (overdispersion) {
    //cout << "(overdispersion)" << endl;
    //}else {
    //cout << "(no overdispersion)" << endl;
    //}
  }


  /*Loop over samples - start at 1*/
  int sampleCounter=1;
  while ( sampleCounter<=sampleSize) {
 
    if ((!verbose) && ((sampleCounter % 10) == 0)) {
      Rprintf(".");
    }
    //Progress bar: 0,..,100% is shown.
    if (sampleCounter > tuneSampleSize && (!verbose) && (sampleCounter % (int)floor(sampleSize/100.0) == 0)) {
      Rprintf("%d%%", sampleCounter*100 / sampleSize);
    }
    if(0){ 
      if(varnu){     
	if ((sampleCounter % 100 == 0)) {
	  Rprintf("alpha\t%f  beta\t%f  %f  gamma[0]\t%f  gamma[1]\t%f  gamma[2]\t%f  %f  lambda\t%f\n", (double)acc_alpha/I, beta[2], (double)acc_beta, gamma[0], gamma[1], gamma[2],(double)acc_gamma, lambda[1][2]);
	  
	  /*      cout<< "alpha\t" << (double)acc_alpha/I<<"  "
		  << "beta\t" <<" "<< beta[2] <<" "<< (double)acc_beta<<"  "
		  << "gamma[0]\t" <<" "<< gamma[0] <<" "<< "gamma[1]\t" <<" "
		  << gamma[1] <<" "<< "gamma[2]\t" <<" "<< gamma[2] <<" "
		  << (double)acc_gamma<<"  "
		  << "lambda\t" << lambda[1][2] << endl;*/
	}
      } 
      if(la_rev){     
	if ((sampleCounter % 100 == 0)) {
	  Rprintf("K\t%d\n", K[1]);
	}
      }
      
      if(delta_rev){     
	if ((sampleCounter % 100 == 0)) {
	  Rprintf("K_delta\t%f  delta[2]\t%f\n", K_delta, delta[2]);
	}
      }
      if(epsilon_rev){     
	if ((sampleCounter % 100 == 0)) {
	  Rprintf("K_epsilon\t%f  epsilon[2]\t%f\n", K_epsilon, epsilon[2]);
	}
      }
    }

    //    cout << ":"<<flush;

    /* Temporary variables, some which are not really used */
    double a,b;  /* c, binp*/
    
    /*Calculate sums */
    double XSum = sumIn2(X,I,n);
 
    double xiSum = 0;
    for (int i=1;i<=I; i++) {
      xiSum = xiSum + xi[i];
    }

    if(!la_estim){
      for (int i=1;i<=I; i++) {
	for (int t=2;t<=n; t++) { 
	  lambda[i][t] =  0;
	}
      }
      if(sampleCounter%3==0){
	acceptedlambda++;
      }
    }
    
    if(la_estim){
      if (!la_rev)
	{
	  //Update x=logit(lambda) mit random walk proposal fur variables nu
	  double xStar;
	  xStar = gsl_ran_gaussian(xRWSigma) + x;
	  double lambdaStar;
	  lambdaStar = invlogit(xStar);
	  double logFx = 0;
	  logFx = log(gsl_ran_beta_pdf(lambda_const,alpha_lambda,beta_lambda)) + log(invlogitd(x));
	  double logFxStar = 0;
	  logFxStar = log(gsl_ran_beta_pdf(lambdaStar,alpha_lambda,beta_lambda)) + log(invlogitd(xStar));
	  for (int i=1;i<=I; i++) {
	    for (int t=2;t<=n; t++) {
	      logFx = logFx + gsl_ran_poisson_log_pdf(Y[i][t], lambda_const*omega[i][t]*Z[i][t-1]);
	      logFxStar = logFxStar + gsl_ran_poisson_log_pdf(Y[i][t], lambdaStar*omega[i][t]*Z[i][t-1]);
	    }
	  }
	  double accx = exp(logFxStar - logFx);
   
	  if (gsl_rng_uniform() <= accx) {x = xStar; lambda_const = invlogit(xStar); acceptedlambda++;}
	  for (int i=1;i<=I; i++) {
	    for (int t=2;t<=n; t++) { 
	      lambda[i][t] =  lambda_const;
	    }
	  }
	}
    }

    if(!varnu){

      double omegaxiSum = 0;
      for (int i=1;i<=I; i++) {
	for (int t=2;t<=n; t++) { 
	  omegaxiSum += omega[i][t]*xi[i];
	}
      }
    
      //Update nu
      a  = alpha_nu + XSum;
      b  = beta_nu + omegaxiSum;
      nu_const =  gsl_ran_gamma (a, 1/b);
      for (int i=1;i<=I; i++) {
	for (int t=2;t<=n; t++) { 
	  nu[i][t] =  nu_const;
	}
      }
    }
    
  



    ///////////////////////////////////////////////////////////////////////////
    //In case of over-dispersion log(psi) has to be updated using a Random Walk
    //Metropolis-Hastings Step.
		
    if (overdispersion) {
      double a_psi = alpha_psi;
      double b_psi = beta_psi;
      if(xi_estim_psi){
	a_psi = 1;
	b_psi = xi_psi;
      }
      //loglik * prior of the old value
      double logFPsi = gsl_ran_gamma_log_pdf(psi,a_psi,1/b_psi) + log(psi);
      for (int i=1; i<=I; i++) {
	for (int t=2; t<=n; t++){
	  logFPsi += gsl_ran_gamma_log_pdf(omega[i][t],psi,1/psi);
	}
      }
      //Generate the new value from Gaussian distrib N(psi, psiRWSigma^2)
      double logpsiStar = gsl_ran_gaussian(psiRWSigma) + log(psi);
      double psiStar = exp(logpsiStar);
      //loglik * prior of the new value
      double logFPsiStar = gsl_ran_gamma_log_pdf(psiStar,a_psi,1/b_psi) + logpsiStar;
      for (int i=1; i<=I; i++) {
	for (int t=2; t<=n; t++) {
	  logFPsiStar += gsl_ran_gamma_log_pdf(omega[i][t],psiStar,1/psiStar);
	}
      }
      //Acceptance prob - fmin(1, <>) superflous.
      double accpsi = exp(logFPsiStar-logFPsi);
      
      //Do we accept?
      if ((psi>0) && (gsl_rng_uniform() <= accpsi)) {psi = psiStar; acceptedPsi++;}
    }

    //update xi_psi
    if(xi_estim_psi){
      double a = alpha_psi + 1;
      double b = beta_psi + psi;
      xi_psi = gsl_ran_gamma (a, 1/b);
    }
    //////////////////////////////////////////////////////////////////////////

    

    //State information to file if we are on an filter'th sample
   
    if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
      logfile << sampleCounter << "\t";
      if (!la_rev){
	logfile << lambda_const  << "\t";
      }
      logfile << psi << "\t";
      logfile << xi_psi << "\t";
      if(!varnu){
	logfile << nu_const << "\t";
      }
    }

    if(varnu){
      // Unterprogramme fuer den Update von alpha und beta
      if (I>=2)
	{
	  alphaupdate(gamma, alpha, beta, delta, lambda, 1, I, n, Y, X, acc_alpha, taualpha, ncov, xcov, xreg, omega, omega, scov,1);
	  taualpha=update_tau_alpha(alpha, I, alpha_a, alpha_b, xreg);

	  if (sampleCounter%3==0)
	    {
	      if(scov==0){
		double asum=0;
		for (int i=1; i<=I; i++)
		  {
		    asum+=(alpha[i]-xreg[i]);
		  }
		for (int i=1; i<=I; i++)
		  {
		    alpha[i]-=(asum/I);
		  }
		gamma[0]=gamma[0]+(asum/I);
	      }
	    }
	}
      else 
	{
	  alpha[1]=0.0;
	}


      //Update fuer zeitlichen effekt mit RW

      if (rw>0)
	{

	  //  update_beta_nurrw(gamma, alpha, beta, delta, X, Z, Y, n, I,  taubeta,  rw, 1, lambda, acc_beta, sampleCounter, my, my2, temp, z, theta, Q, Q2, L, L2, xcov, ncov, scov, omega, omega, 1);
	  //update_beta_block(alpha, beta, gamma, delta, X, n, I, taubeta, rw, acc_beta, sampleCounter, n1, n2, my, my2, z, theta, beta0, Q, Q2, L, L2, xcov, ncov, scov, omega);

	  /*hofmann - no fortran
	    update_beta_tau_block(alpha, beta, gamma, delta, beta_a, beta_b, X, n, I, taubeta, rw, acc_beta, taubetaRWSigma, taubetaStar, sampleCounter, n1, n2, my, my2, z, theta, beta0, Q, Q2, L, L2, xcov, ncov, scov, omega);
	  */
	  //taubeta=beta_a/beta_b;
	  // taubeta=hyper(rw, beta, beta_a, beta_b, n);
	  //taubeta=720;

	  //if(sampleCounter%500==1){cout << taubeta << endl << endl;}
	  //  for(int t=2;t<=n;t++){
	  //  update_beta_t(t, alpha, beta, gamma, delta, ncov, xcov, X, n, I, taubeta, acc_beta, omega, scov);
	  //   }
  
	 
	  if(scov==0){
	    //  if (sampleCounter%1==0)
	    //   {
   	    double bsum=0;
	    for (int t=2; t<=n; t++)
	      {
	        bsum+=(beta[t]);
	      }
	    for (int t=2; t<=n; t++)
	      {
	        beta[t]-=(bsum/(n-1));
	      }
	    gamma[0]=gamma[0]+(bsum/(n-1));
	    //    }
	  }
	 
	} //if (rw>0)


      //update saison

      //update_gamma( alpha,  beta, gamma,ncov,  xcov,  X,  Z, Y, n, I, taugamma, 1, lambda, acc_gamma, P, P2, gammaalt, z2, L, Q, omega, omega,1);
      taugamma=gamma_b;
      // cout << gamma[0]<<"  "  << gamma[1] << endl;
      for(int j=scov;j<ncov;j++){
	update_gamma_j(j, alpha, beta, gamma, delta, ncov, xcov, X, n, I, taugamma, gammaneu, acc_gamma, omega, scov);
	// gamma[j]=0;
      }
      //    cout << gamma[0]<<"  "  << gamma[1] << endl;
      //cout << "lambda [1][2]   " << lambda [1][2] << endl;
      //cout << "X [1][2]   " << X [1][2] << endl;
      //cout << "omega [1][2]   " << omega [1][2] << endl;
      // update delta_br
      if(delta_rev){
	update_delta_br(delta, delta_br, xi_delta, breakpoints_delta, breakpointsStar_delta, K_delta, KStar_delta, Km1_delta, delta_a, delta_b, X, nu, n, I, acceptedbr_delta, omega, xi_estim_delta, K_geom, p_K, alpha_xi, beta_xi);
      }
      //cout << delta[2]<<"  " << K_delta << endl;

      // Berechnet nu_it=log(alpha_i+beta_t)
      machnu(gamma, alpha, beta,  delta, nu, I, n, ncov, xcov, scov);

      //cout << nu[1][2] << endl << endl;
      //cout << beta[2] << endl << endl;
      //       cout << "test" << endl;
      //       cout << "test2" << endl;


      if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
	//     for (int i=1;i<=I; i++) {
	//       for (int t=1; t<=n; t++) {
	// 	 logfile << nu[i][t] << "\t";
	//       }
	//     }
	//    logfile << mu << "\t";
	for (int j=0; j<ncov; j++) {
	  logfile << gamma[j] << "\t";
	}
	if(delta_rev){
	  if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
	    logfile << Km1_delta<<"\t"<< xi_delta<<"\t";
	    for (int j=2; j<=n; j++) {
	      logfile << delta[j] << "\t";
	    }
	  }

	  if (sampleCounter>burnin) {
	    for (int k=1; k<=K_delta; k++) {
	      for (int j=2; j<=n; j++) {
		if (breakpoints_delta[k]==j){
		  bp_delta[j]+=1;
		}
	      }
	    }
	  }
	}//if(delta_rev)

      }//if

    }//if varnu

    if(epsilon_rev){
      update_epsilon_br(epsilon, epsilon_br, xi_epsilon, breakpoints_epsilon, breakpointsStar_epsilon, K_epsilon, KStar_epsilon, Km1_epsilon, epsilon_a, epsilon_b, S, n, I, acceptedbr_epsilon, omega, xi_estim_epsilon, K_geom, p_K, alpha_xi, beta_xi);

      if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
        logfile << Km1_epsilon<<"\t"<< xi_epsilon<<"\t";
	for (int j=2; j<=n; j++) {
	  logfile << epsilon[j] << "\t";
        }
      }

      if (sampleCounter>burnin) {
        for (int k=1; k<=K_epsilon; k++) {
          for (int j=2; j<=n; j++) {
            if (breakpoints_epsilon[k]==j){
  	      bp_epsilon[j]+=1;
            }
          }
        }
      }
    }//if(epsilon_rev)


    if(la_estim){
      if (la_rev)
	{
	  update_lambda_br(lambda, lambda_br, xi_lambda, breakpoints, breakpointsStar, K, KStar, Km1, alpha_lambda, beta_lambda, Y, Z, n, I, acceptedbr, omega, theta_pred_estim, xi_estim, K_geom, p_K, alpha_xi, beta_xi);



	  if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
	    logfile << Km1[1]<<"\t"<< xi_lambda[1]<<"\t";
	    for (int j=2; j<=n; j++) {
	      logfile << lambda[1][j] << "\t";
	    }
	  }
	  for (int i=1;i<=I; i++) {
	    if (sampleCounter>burnin) {
	      for (int k=1; k<=K[i]; k++) {
		for (int j=2; j<=n; j++) {
		  if (breakpoints[i][k]==j){
		    bp[i][j]+=1;
		  }
		}
	      }
	    }
	  }
	}//if(la_rev)
    } // if(la_estim)
    // cout << S[1][106] << endl;
    //   cout << "test" << endl;
    //Loop over the individual X[t], Y[t], S[t], and omega[t]
    for (int i=1;i<=I; i++) {
      for (int t=2; t<=n; t++) {
        //Update X
        double binp = nu[i][t]*xi[i] / (epsilon[t] + nu[i][t]*xi[i] + lambda[i][t] * Z[i][t-1]);
	X[i][t] =  gsl_ran_binomial( binp, Z[i][t]);

        //Update S
        binp =  epsilon[t] / (epsilon[t] + lambda[i][t] * Z[i][t-1]);
	//hoehle 9 Apr 2009 -- protection against Z[i][t-1]==0 case, leading to binp = nan
	if (Z[i][t-1] == 0) {binp = 1;}
	S[i][t] =  gsl_ran_binomial( binp, (Z[i][t] - X[i][t]));
  
    
        //Update Y
        Y[i][t] = Z[i][t] - X[i][t] - S[i][t];
        
	//Debug
	//cout << "i=" << i << "\tt=" << t << "\tX=" << X[i][t] << "\tY=" << Y[i][t] << "\tZ=" << Z[i][t] << "\tS=" << S[i][t] << "\tepsilon=" << epsilon[t] << "\tbinp=" << binp << endl;

        //Update omega[t] in case of overdispersion
	if(overdispersion){
	  double a = psi + Z[i][t];
	  double b = psi + epsilon[t] + nu[i][t] + lambda[i][t]*Z[i][t-1];
	  omega[i][t] = gsl_ran_gamma(a,1/b);
	}
        //Write state to log-file.
        if (sampleCounter>burnin) {
	  sumX[i][t] += X[i][t]; sumY[i][t] += Y[i][t]; sumS[i][t] += S[i][t]; sumomega[i][t] += omega[i][t]; Sumeta[i][t] += eta[i][t]; Sumvarr[i][t] += varr[i][t]; Sumrpearson[i][t] += rpearson[i][t];
	}
      }//for t
    }//for i
    //   cout << "test2" << endl;
    // cout << Z[1][2] << endl;
    // cout << X[1][2] << endl;
    // cout << Y[1][2] << endl;
    // cout << S[1][2] << endl;


    //Praediktive Verteilung fuer variables nu
    for (int i=1;i<=I;i++) {
      if(!theta_pred_estim){
	double p_thetanp1 = ((double(K[i]))/double(n)); //(1+double(K[i]))
	if(K_geom){
	  p_thetanp1 = (double(K[i])*(1.0-p_K)*(1.0-pow((double)1.0-p_K,double(n-1))))/((double(n)-1.0)*(1.0-pow((double)1.0-p_K,double(n))));
	}
        if(gsl_rng_uniform()<=p_thetanp1){
	  if (sampleCounter>burnin) {
	    bp[i][n+1] += 1;
	  }
	  double alpha_la = alpha_lambda;
	  double beta_la = beta_lambda;
	  if(xi_estim){
	    beta_la = xi_lambda[i];
	  }
          lambda[i][n+1]=gsl_ran_gamma(alpha_la,1/beta_la);
        }
      }
      if(overdispersion){
	omeganp1[i] =  gsl_ran_gamma(psi,1/psi);
      }else{
	omeganp1[i] = 1;
      }
      if(varnu){
	a = 0;
	for(int j=scov;j<ncov;j++){
	  a += gamma[j]*xcov[j][n+1];
	}
	if(rw>0){
	  a += gsl_ran_gaussian(sqrt(1/taubeta)) + (2*beta[n-1]-beta[n]);
	}
	if(delta_rev){
	  double p_thetanp1 = ((double(K[i]))/double(n)); //(1+double(K[i]))
	  if(K_geom){
	    p_thetanp1 = ((double(K[i]))*(1.0-p_K)*(1.0-pow((double)1.0-p_K,double(n-1))))/((double(n)-1.0)*(1.0-pow((double)1.0-p_K,double(n))));
	  }
	  if(gsl_rng_uniform()<=p_thetanp1){
	    if (sampleCounter>burnin) {
	      bp_delta[n+1] += 1;
	    }
	    double alpha_de = delta_a;
	    double beta_de = delta_b;
	    if(xi_estim){
	      beta_de = xi_delta;
	    }
	    delta[n+1]=gsl_ran_gamma(alpha_de,1/beta_de);
	  }
	  a += log(delta[n+1]);
	}	  
	nunp1[i] = exp(a);
      }else{
	nunp1[i]=nu[i][n];
      }
      if(epsilon_rev){
	double p_thetanp1 = ((double(K[i]))/double(n)); //(1+double(K[i]))
	if(K_geom){
	  p_thetanp1 = ((double(K[i]))*(1.0-p_K)*(1.0-pow((double)1.0-p_K,double(n-1))))/((double(n)-1.0)*(1.0-pow((double)1.0-p_K,double(n))));
	}
	if(gsl_rng_uniform()<=p_thetanp1){
	  if (sampleCounter>burnin) {
	    bp_epsilon[n+1] += 1;
	  }
	  double alpha_ep = epsilon_a;
	  double beta_ep = epsilon_b;
	  if(xi_estim){
	    beta_ep = xi_epsilon;
	  }
	  epsilon[n+1]=gsl_ran_gamma(alpha_ep,1/beta_ep);
	}
      }	  

      Xnp1[i] = gsl_ran_poisson(omeganp1[i]*nunp1[i]*xi[i]);
      Ynp1[i] = gsl_ran_poisson(lambda[i][n+1]*omeganp1[i]*(Z[i][n]));
      Snp1[i] = gsl_ran_poisson(omeganp1[i]*epsilon[n+1]);
      Znp1[i] = Xnp1[i] + Ynp1[i] + Snp1[i];
      if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
	logfile << Znp1[1] << "\t";
      }
    }
 
     

    if ((sampleCounter>burnin) && ((sampleCounter-burnin) % filter == 0)) {
      logfile << satdev(n,I,Z,lambda,nu,xi,epsilon,eta,psi,overdispersion) << endl;
    }

	
    logfile.flush();

    //Tuning	
    if(sampleCounter == tuneSampleSize){
      if (!la_rev)
	{
	  Rprintf("Current xRWSigma= %f --> acc rate= %f\n", xRWSigma, acceptedlambda/tuneSampleSize);
	  tune(xRWSigma, acceptedlambda, tuneSampleSize,tunex);
	  Rprintf("Corrected xRWSigma= %f\n", xRWSigma);
	}
      if(overdispersion){
	Rprintf("\nCurrent psiRWSigma= %f --> acc rate = %f\n", psiRWSigma, acceptedPsi/tuneSampleSize);
        tune(psiRWSigma, acceptedPsi, tuneSampleSize,tunepsi);
        Rprintf("Corrected psiRWSigma= %f\n", psiRWSigma);
      }
      
      if(varnu&&(rw>0)){
        Rprintf("Current taubetaRWSigma= %f --> acc rate %f\n", taubetaRWSigma, acc_beta/tuneSampleSize);
        tune(taubetaRWSigma, acc_beta, tuneSampleSize,tunetaubeta,0.1,0.4);
        Rprintf("Corrected taubetaRWSigma= %f\n", taubetaRWSigma);
      }
      
      //tunetaubeta = 0;

      need=tunex + tunepsi + tunetaubeta;
      if(need > 0){
        acceptedlambda = 0;
        acceptedbr = 0;
        acceptedbr_delta = 0;
        acceptedbr_epsilon = 0;
        acceptedPsi = 0;
        sampleCounter = 0;
        if(varnu){
          acc_beta=0;
          acc_alpha=0;
          acc_gamma=0;
        }
        //Fix seed of generator to reproduce results.
	//        gsl_rng_set(r,seed);
      }//if
    }//if

    sampleCounter++;
  }//while counter
      

  //Write means to logfile2

  for (int t=1;t<=n;t++) {
    logfile2 << (double)sumX[1][t]/((double)samples*(double)filter) << "\t" << (double)sumY[1][t]/((double)samples*(double)filter)<< "\t" << (double)sumomega[1][t]/((double)samples*(double)filter) << "\t"<< (double)bp[1][t]/((double)samples*(double)filter) << "\t";
  }
  logfile2 << (double)bp[1][n+1]/((double)samples*(double)filter) << "\t";
  logfile2 << endl;


  //Write accepted status to file 
  if(overdispersion){acclog << "psi\t" <<  psiRWSigma << "\t" << (double)acceptedPsi/(double)sampleSize << endl;}
  if (!la_rev){acclog << "lambda\t" << xRWSigma  << "\t" << (double)acceptedlambda/(double)sampleSize << endl;}
  if (la_rev){acclog << "br\t" << 0  << "\t" << (double)acceptedbr/(double)sampleSize << endl;}
  if(I>1){acclog << "alpha\t" << 0 <<"\t" <<(double)acc_alpha/((double)sampleSize*I)<<endl;}
  if(varnu&&(rw>0)){acclog <<"beta\t"<<0 <<"\t"<< (double)acc_beta/((double)sampleSize*(double)(n-1.0))<<endl;}
  if(varnu){acclog <<"gamma\t"<<0 <<"\t"<< (double)acc_gamma/((double)sampleSize*(double)ncov)<<endl;}
  if(delta_rev){acclog <<"brdelta\t"<<0 <<"\t"<< (double)acceptedbr_delta/((double)sampleSize)<<endl;}
  if(epsilon_rev){acclog <<"brepsilon\t"<<0 <<"\t"<< (double)acceptedbr_epsilon/((double)sampleSize)<<endl;}
  //cout << "done"<<endl;
  return;
    



}//funktion

/* hoehle: interface for calling twins from R */

extern "C" { 

  void twins (int *x_ptr, int *n_ptr, int *I_ptr,
	      char **logFile_ptr, char **logFile2_ptr, 
	      int *burnin_ptr, int *filter_ptr,  int *sampleSize_ptr, 
	      double *alpha_xi_ptr,  double *beta_xi_ptr, 
	      int *T_ptr, int *nfreq_ptr, 
	      double *psiRWSigma_ptr, double *alpha_psi_ptr, double *beta_psi_ptr,		int *nu_trend_ptr) {

  //Splash screen
  Rprintf("MCMC Estimation in BPLE Model v1.0.1 (using R API).\n");

  /* Datafile and Logfile variables */
  char *logFile  = *logFile_ptr;
  char *logFile2 = *logFile2_ptr;

  /* Fix a lot of tuning parameters */
  overdispersion = 1;
  alpha_lambda = 1;
  beta_lambda = 1;
  alpha_nu = 1;
  beta_nu = 1;
  xRWSigma = 1;
  varnu = 1;
  alpha_a = 1;
  alpha_b = 1;
  beta_a = 1;
  beta_b = 1;
  int rw = 0;
  la_rev = 1;
  /* Daniel: now nu_trend is set through the R option "nu_trend" 
     nu_trend = 0; */
  theta_pred_estim = 0;
  xi_estim = 1;
  delta_rev = 0;
  xi_estim_delta = 0;
  delta_a = 1;
  delta_b = 1;
  epsilon_rev = 0;
  xi_estim_epsilon = 0;
  epsilon_a = 1;
  epsilon_b = 1;
  la_estim = 1;
  xi_estim_psi = 0;
  K_geom = 0;
  p_K = 0;
  gamma_a = 1;
  gamma_b = 0.000001;

  /* hoehle -- new code to fetch params  if called through interface */
  int burnin = *burnin_ptr;
  int filter = *filter_ptr;
  int sampleSize = *sampleSize_ptr;
  int T = *T_ptr;
  int nfreq = *nfreq_ptr;

  alpha_xi = *alpha_xi_ptr;
  beta_xi = *beta_xi_ptr;
  psiRWSigma = *psiRWSigma_ptr;
  alpha_psi = *alpha_psi_ptr;
  beta_psi = *beta_psi_ptr;
  nu_trend = *nu_trend_ptr;

  /* Status information */
  Rprintf("dim(x) = %d\t%d\n", *n_ptr, *I_ptr);
  Rprintf("logfile is in \"%s\".\n",logFile);
  Rprintf("logfile2 is in \"%s\".\n", logFile2);
  Rprintf("burnin = %d (%d)\n", burnin, *burnin_ptr);
  Rprintf("filter = %d (%d)\n", filter, *filter_ptr);
  Rprintf("sampleSize = %d (%d)\n", sampleSize, *sampleSize_ptr);
  Rprintf("T = %d\n", T);
  Rprintf("nfreq = %d\n",nfreq);
  Rprintf("alpha_xi = %f\n", alpha_xi);
  Rprintf("beta_xi = %f\n", beta_xi);
  Rprintf("psiRWSigma = %f\n", psiRWSigma);
  Rprintf("alpha_psi = %f\n", alpha_psi);
  Rprintf("beta_psi = %f\n", beta_psi);
  Rprintf("nu_trend = %d\n", nu_trend);


  /*********************************************************************** 
   * Open the log file for reading 
   ***********************************************************************/
  ofstream logfile,logfile2,accfile;

  char accFile[200];
  sprintf(accFile, "%s.acc", logFile);

  logfile.open(logFile);
  logfile2.open(logFile2);
  accfile.open(accFile);

  if (!logfile) { Rf_error("Error opening the log file.\n");}
  if (!accfile) { Rf_error("Error opening the acc file.\n");}

  /* Allocate a random number generator -- this is now the R RNG and
     Fix seed of generator to reproduce results (i.e. fetch current seed
     value from R.
  */
  GetRNGstate();

  //Read Data
  int I = *I_ptr;
  int n = *n_ptr;

  LongMatrix Z = surveillancedata2twin(x_ptr,n,I);
  /*  Z = readData(dataFile,&n,&I);*/

  double xi[2] = {0, 1};
  // instead of
  // xi[1] = 1;


  Rprintf(" ====== The data =======\n");
  for (int t=0; t<=n; t++) {
    for (int i=0; i<=I; i++) {
      Rprintf("%d\t", Z[i][t]);
    }
    Rprintf("\n");
  }

    

  /**********************************************************************
   * Do the MCMC estimation, results are written to log.txt
   **********************************************************************/
  bplem_estimate(0,logfile,logfile2,accfile,Z,xi,n,I,T,nfreq,burnin,filter,sampleSize, rw);

  logfile.close();
  logfile2.close();
  accfile.close();

  // Now we can go back to R...
  Rprintf("\nDone with twins -- going back to R.\n");

  //Done - save the current seed value for use in R
  PutRNGstate();
} /* end of twins */

} /* end of extern "C" */
