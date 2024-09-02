#include <R.h>
#include <Rinternals.h>

//.Call(C_eq3a, as.numeric(lambda.old), as.numeric(Y), as.numeric(incu.pmf))
SEXP eq3a(SEXP rlambdaOld, SEXP  ry, SEXP rincuPmf) {

  int T, i, t, d;

  // get arguments
  double *lambdaOld = REAL(rlambdaOld);
  T = LENGTH(rlambdaOld);
  double *y = REAL(ry);              // also of length T
  double *incuPmf = REAL(rincuPmf);  // for times 0, ..., d_max

  // Create long enough vectors for queries about dincu and pincu
  double *dincu = (double*)R_alloc(T, sizeof(double));
  double *pincu = (double*)R_alloc(T, sizeof(double));
  dincu[0] = incuPmf[0];
  pincu[0] = dincu[0];
  for (i=1; i<LENGTH(rincuPmf); i++) {
    dincu[i] = incuPmf[i];
    pincu[i] = pincu[i-1] + dincu[i];
  }
  for (i=LENGTH(rincuPmf); i<T; i++) {
    dincu[i] = 0.0;
    pincu[i] = 1.0;
  }

  // Allocate result vector
  SEXP rphiNew = PROTECT(allocVector(REALSXP, T));
  double *phiNew = REAL(rphiNew);

  // many loops
  for(t = 0; t < T; ++t)
    {
      double sum3a = 0.0;
      
      for(d = 0; d <= T - (t + 1); ++d)
	{
	  double tmp = 0.0;
	  
	  for(i = 0; i < t + d; ++i)
	    {
	      tmp += lambdaOld[i] * dincu[t + d - i];
	    }
	  
	  tmp = dincu[d] / tmp;
	  
	  if(!R_FINITE(tmp))
	    {
	      tmp = 0.0;
	    }
	  
	  sum3a += y[t + d] * tmp;	    
	}
      
      //printf("Querying index %d\\n", T-(t+1));
      
      phiNew[t] = lambdaOld[t] / pincu[T - (t + 1)] * sum3a;
      
      if(!R_FINITE(phiNew[t]))
	{
	  phiNew[t] = 0.0;
	}
    }

  //Show values
  //for (i=0; i<T; ++i) {
  //  printf( "%f\\n", phiNew[i]);
  //}

  UNPROTECT(1);
  return(rphiNew);

}
