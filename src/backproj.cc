#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP eq3a(SEXP rlambdaOld, SEXP  ry, SEXP rincuPmf) {
BEGIN_RCPP

    // get arguments
  NumericVector lambdaOld(rlambdaOld);
  int T = lambdaOld.length();
  NumericVector y(ry);
  NumericVector incuPmf(rincuPmf);

  // Create long enough vectors for queries about dincu and pincu
  NumericVector dincu(T);
  NumericVector pincu(T);
  pincu[0] = dincu[0];    
  for (int i=1; i<incuPmf.length(); i++) {
    dincu[i] = incuPmf[i];
    pincu[i] = pincu[i-1] + dincu[i];
  }
  for (int i=incuPmf.length(); i<T; i++) {
    dincu[i] = 0.0;
    pincu[i] = 1.0;
  }  
  
      
  // Allocate result vector
  NumericVector phiNew(T);
  
  // many loops
  for(int t = 0; t < T; ++t)
    {
      double sum3a = 0.0;
      
      for(int d = 0; d <= T - (t + 1); ++d)
	{
	  double tmp = 0.0;
	  
	  for(int i = 0; i < t + d; ++i)
	    {
	      tmp += lambdaOld[i] * dincu[t + d - i];
	    }
	  
	  tmp = dincu[d] / tmp;
	  
	  if(R_IsNaN(tmp) || (! R_finite(tmp)))
	    {
	      tmp = 0.0;
	    }
	  
	  sum3a += y[t + d] * tmp;	    
	}
      
      
      //printf("Querying index %d\\n", T-(t+1));
      
      phiNew[t] = lambdaOld[t] / pincu[T - (t + 1)] * sum3a;
      
      if(R_IsNaN(phiNew[t]) || (! R_finite(phiNew[t])))
	{
	  phiNew[t] = 0.0;
	}
    }
  
  //Show values
  //for (int i=0; i<T; ++i) {
  //  printf( "%f\\n", phiNew[i]);
  //}
  
  //return List::create(Named( "phi" = phiNew));
  return(wrap(phiNew));

END_RCPP
}
