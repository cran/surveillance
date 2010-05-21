/*******************************************************************
 * Author: Michael Höhle <hoehle@stat.uni-muenchen.de>
 * Date:   Aug 2008 *
 *
 * Header file containing wrappers for GSL related calls
 * to R calls using the R API. This code is used in twins.cc
 *******************************************************************/

/* new definitions to replace GSL code */

// Remove the dead RNG variable (DSB 04/05/2010):
// int r;


double gsl_rng_uniform () {
  //  GetRNGstate();
  double res = runif(0,1);
  //PutRNGstate();
  return(res);
}

double gsl_ran_gaussian(double sigma) {
  //GetRNGstate();
  double res = rnorm(0.0,sigma);
  //PutRNGstate();
  return(res);
}

double gsl_ran_gamma(double a, double b) {
  //GetRNGstate();
  double res = rgamma(a,b);
  //PutRNGstate();
  return(res);
}

unsigned int gsl_ran_poisson(double lambda) {
  //GetRNGstate();
  unsigned int res = rpois(lambda);
  //PutRNGstate();
  return(res);
}


unsigned int gsl_ran_binomial(double p, unsigned int n) {
  //GetRNGstate();
  unsigned int res = rbinom(n,p);
  //PutRNGstate();
  return(res);
}

//hoehle: The original function assumes mu>0, which needs not be the case!
//This version handles that part. This is the log version.

double
gsl_ran_poisson_log_pdf (const unsigned int k, const double mu)
{
  double p;
  if (mu==0) {
    return(log((double)(k == 0)));
  } else {
    double lf = lgammafn(k+1); /*gsl2R: gsl_sf_lnfact(k) */

    p = k*log(mu) - lf - mu;
    return p;
  }
}

double gsl_sf_lngamma(double x) {
  return(lgammafn(x));
}

double gsl_ran_beta_pdf (double x, double a, double b) {
  return(dbeta(x,a,b,0));
}

/**********************************************************************
 * Log version of the Gamma pdf with mean a*b and variance a*b^2.
 *
 **********************************************************************/

double
gsl_ran_gamma_log_pdf (const double x, const double a, const double b)
{
  if (x < 0)
    {
      //This is problematic!
      return log((double)0) ;
    }
  else if (x == 0)
    {
      if (a == 1)
        return log(1/b) ;
      else
        return log((double)0) ;
    }
  else if (a == 1)
    {
      return -x/b - log(b) ;
    }
  else 
    {
      double p;
      /*gsl2R:      double lngamma = gsl_sf_lngamma (a);*/
      double lngamma = lgammafn(a);
      p = (a-1)*log(x) - x/b - lngamma - a*log(b);
      return p;
    }
}

/* Seed random number generator */
//void gsl_rng_set(int r, long seed) {
//  set.seed(seed);
//}
