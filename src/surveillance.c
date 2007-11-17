/**
   
  C routines for the surveillance package

  Author: Michael Höhle <http://www.stat.uni-muenchen.de/~hoehle>
  Date:   21 Sep 2007

  Atm the only C routines are concerned with the GLR computations
  in the algorithm algo.prc.

  //should check that these really work...
  void lr_cusum  - intercept chart with known kappa
  void glr_cusum - intercept chart with estimated kappa
  void glr_cusum_window -- window limited intercept chart with estimated kappa

  //removedvoid glr_epi
  void glr_epi_window

  //History
  21 Sep 2007 -- modified code to get around error of extreme strict (=pedantic) MacOS compiling on CRAN
  28 Nov 2006 -- file created
*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**********************************************************************
  C implementation of the LR test for the seasonal Poisson chart with
  fixed change in the intercept

  Params:
   x   - array of observed values (pos 0 is first interesting value)
  mu0  - array with the means once in-control (pos 0 is first interesting value)
  lx   - length of the x and mu0 array
  kappa- the change in intercept to detect (here known in advance)

  c_ARL- when to sound alarm threshold
  ret_N- here the return value is stored
  ret_lr- GLR value for each n to be returned
**********************************************************************/

void lr_cusum(int* x,double* mu0, int *lx_R, double *kappa_R, double *c_ARL_R,int *ret_N, double *ret_lr) {
  /* Pointers to something useful */
  int lx = *lx_R;
  double c_ARL = *c_ARL_R;
  double kappa = *kappa_R;
  
  /* Loop variables */
  register int n=0; 
  int stop = 0;
  int N = lx;

  /* Loop over all 0 <= n <= length(x) */
  while ((n < lx)) {
    /*Compute for one n*/
    /*printf("n=%d\n",n);*/

    
    double zn = kappa * x[n] + (1-exp(kappa))*mu0[n];

    /* Add up */
    if (n==0) {
      /*ret_lr[n] = zn;*/
      ret_lr[n] = fmax(0,zn);
    } else {
      /*ret_lr[n] = zn + fmax(0,ret_lr[n-1]);*/
      ret_lr[n] = fmax(0,ret_lr[n-1] + zn);
    }

    /* Find the first time that the GLR increases c_ARL there we stop */
    if ((ret_lr[n] > c_ARL) && !stop) {
      N = n;
      stop = 1;
      break;
    }
  
    /* Advance counter */
    n++;
  }

  /* Return value (add 1 for R/SPlus array compability */
  *ret_N = N+1;
}






/**********************************************************************
  Fast C implementation of the sequential GLR test without windowing
  for Poisson distributed variables

  Params:
   x   - array of observed values (pos 0 is first interesting value)
  mu0  - array with the means once in-control (pos 0 is first interesting value)
  lx   - length of the x and mu0 array
  n0   - number of burn-in values (number of observations, not array index!)
  c_ARL- when to sound alarm threshold
  ret_N- here the return value is stored
  ret_glr- GLR value for each n to be returned
  dir  - direction of testing
**********************************************************************/

void glr_cusum(int* x,double* mu0, int *lx_R, int *n0_R, double *c_ARL_R,int *ret_N, double *ret_glr, int *dir_R) {
  /* Pointers to something useful */
  
  int lx = *lx_R;
  int n0 = *n0_R;
  int dir = *dir_R;
  double c_ARL = *c_ARL_R;
  
  /* Loop variables */
  register int n, k,l; /*n0-1*/
  for (n=0; n<n0-1; n++) { ret_glr[n] = 0; }
  
  int stop = 0;
  int N = lx;


  /* Precalculation of log(mu0) */
  double logmu0[lx];
  for (l=0;l<lx;l++) { logmu0[l] = log(mu0[l]); }

  /* Debug */
  /* for (l=0;l<lx;l++) {printf("logmu0[%d] = log(%f) = %f\n",l,mu0[l],logmu0[l]); } */

  
  /* Loop over all n0 <= n <= length(x) */
  while ((n < lx)) {
    /* Compute for one n */
    /* printf("n=%d\n",n); */

    /* Define max of the GLR stats */
    double maxGLR = -1e99;

    /* For the recursive computation of kappa_ml */
    double sumx = 0;
    double summu0 = 0;
    

    /* Loop over all k */
    for (k=n; k>=0; k--) { /* Backwards loop makes calculations faster */
      /* Recursive update of the kappa.ml quantitities */
      sumx += x[k];
      summu0 += mu0[k];
      /* Calculate MLE of kappa */
      double kappa_ml = dir*fmax(0,dir*log(sumx/summu0));
      
      

      /* Calculate sum of likelihood ratios (See deriv (1) on 20.9.2006 sheet)*/
      /* 
	 double sum = 0;
	 for (l=k; l<=n; l++) {
	 sum += x[l]*kappa_ml + mu0[l]*(1-exp(kappa_ml));
	 }
      */

      /* Debug */
      /* printf("For n=%d, k=%d we have one.k = %f\n",n,k,sum); */

      /* Recursive updating of the likelihood ratios -- See
	 notes on the 21 september printout. This is fast! */
      double sum = kappa_ml * sumx + (1-exp(kappa_ml))*summu0;

      /* if (sum- sum2 > 1e-6) {
	 printf("For n=%d, k=%d we have one.k = %f\n",n,k,sum);
	 printf("fast impl %f\n",sum2);
	 } 
      */
      
      /* Save the max value */
      if (sum > maxGLR) { maxGLR = sum;}
    }

    /* Save the return value */
    ret_glr[n] = maxGLR;

    /* Debug */
    /* printf("For n=%d the best GLR value is %f\n",n,maxGLR); */

    /* Find the first time that the GLR increases c_ARL there we stop */
    if ((maxGLR > c_ARL) && !stop) {
      N = n;
      stop = 1;
      break;
    }
  
    /*Advance counter*/
    n++;
  }

  /* Return value (add 1 for R/SPlus array compability */
  *ret_N = N+1;
}

/**********************************************************************
  Fast C implementation of the sequential GLR test without windowing
  for Poisson distributed variables

  Params:
   x   - array of observed values (pos 0 is first interesting value)
  mu0  - array with the means once in-control (pos 0 is first interesting value)
  lx   - length of the x and mu0 array
  Mtilde  - number of vals we will need to estimate a detection
  M - max time to go back in time from N
  c_ARL- when to sound alarm threshold
**********************************************************************/

void glr_cusum_window(int* x,double* mu0, int *lx_R, int *M_R, int *Mtilde_R, double *c_ARL_R,int *ret_N) {
  /* Pointers to something useful */
  int lx = *lx_R;
  int M = *M_R;
  int Mtilde = *Mtilde_R;
  double c_ARL = *c_ARL_R;

  /* Loop variables (n>Mtilde, so we start with n=Mtilde (due to -1 in index) */
  register int n = Mtilde, k,l;
  int stop = 0;
  int N = lx;

  /* Precalculation of log(mu0) */
  double logmu0[lx];
  for (l=0;l<lx;l++) { logmu0[l] = log(mu0[l]); }

  /* Debug */
  /* for (l=0;l<lx;l++) {printf("logmu0[%d] = log(%f) = %f\n",l,mu0[l],logmu0[l]); } */

  
  /* Loop over all Mtilde <= n <= length(x) */
  while ((n < lx)) {
    /* Compute for one n */
    /* printf("n=%d\n",n);*/

    /* Define max of the GLR stats */
    double maxGLR = -1e99;

    /* For the recursive computation of kappa_ml compute for (n-Mtilde+1):n */
    double sumx = 0;
    double summu0 = 0;
    for (l=n-Mtilde+1; l<=n; l++) {
      sumx += x[l];
      summu0 += mu0[l];
    }

    /* Loop over all max(0,n-M) <= k <= n-Mtilde -- do this backwards */
    /* for (k=max(0,n-M); k<= (n-Mtilde); k++) { */
    for (k=n-Mtilde; k>=fmax(0,n-M); k--) {

      /* Recursive update of the kappa.ml quantitities */
      sumx += x[k];
      summu0 += mu0[k];
      double kappa_ml = fmax(0,log(sumx/summu0));

      /* Calculate sum of likelihood ratios */
      /* 
	 double sum = 0;
	 for (l=k; l<=n; l++) {
 	sum += x[l]*kappa_ml + mu0[l]*(1-exp(kappa_ml));
	}
      */

      /*Calculate sum of likelihood ratios using recursive updating (fast!)*/
      double sum = kappa_ml * sumx + (1-exp(kappa_ml))*summu0;

      /*
	if (sum- sum2 > 1e-6) {
	printf("For n=%d, k=%d we have one.k = %f\n",n,k,sum);
	printf("fast impl %f\n",sum2);
      }
      */
    

      /* Debug */
      /* printf("For n=%d, k=%d we have one.k = %f\n",n,k,sum);*/

      /* Save the max value */
      if (sum > maxGLR) { maxGLR = sum;}
    }

    /* Debug*/
    /* printf("For n=%d the best GLR value is %f\n",n,maxGLR);*/

    /* Find the first time that the GLR increases c_ARL there we stop */
    if ((maxGLR > c_ARL) && !stop) {
      N = n;
      stop = 1;
      break;
    }
  
    /* Advance counter */
    n++;
  }

  /* Return value (add 1 for R/SPlus array compability */
  *ret_N = N+1;
}


void glr_cusum_window_old(int* x,double* mu0, int *lx_R, int *M_R, int *Mtilde_R, double *c_ARL_R,int *ret_N) {
  /* Pointers to something useful */
  int lx = *lx_R;
  int M = *M_R;
  int Mtilde = *Mtilde_R;
  double c_ARL = *c_ARL_R;

  /* Loop variables (n > Mtilde, i.e. n > Mtilde-1+1) */
  register int n = Mtilde, k,l;
  double mu1[lx];
  int stop = 0;
  int N = lx;

  /* Precalculation of log(mu0) */
  double logmu0[lx];
  for (l=0;l<lx;l++) { logmu0[l] = log(mu0[l]); }

  /* Debug */
  /* for (l=0;l<lx;l++) {printf("logmu0[%d] = log(%f) = %f\n",l,mu0[l],logmu0[l]); } */

  
  /* Loop over all Mtilde <= n <= length(x) */
  while ((n < lx)) {
    /* Compute for one n */
    /* printf("n=%d\n",n); */

    /* Define max of the GLR stats */
    double maxGLR = -1e99;

    /* Loop over all max(0,n-M) <= k <= n-Mtilde (+1 to adjust for start?!?!) */
    for (k=(int)fmax(0,n-M); k<= (n-Mtilde); k++) {
      
      /* Calculate kappa.ml */
      double sumx = 0;
      double summu0 = 0;
      for (l=k; l<=n; l++) {
	sumx += x[l];
	summu0 += mu0[l];
      }
      double kappa_ml = fmax(0,log(sumx/summu0));

      /* Value of the mean under H1 */
      for (l=k; l<=n; l++) {
	mu1[l] = mu0[l]*exp(kappa_ml);
      }

      /* Calculate sum of likelihood ratios */
      double sum = 0;
      for (l=k; l<=n; l++) {
	sum += x[l]*(log(mu1[l]) - logmu0[l]) - mu1[l]+mu0[l];
      }

      /* Debug */
      /*printf("For n=%d, k=%d we have one.k = %f\n",n,k,sum); */

      /*Save the max value */
      if (sum > maxGLR) { maxGLR = sum;}
    }

    /*Debug */
    /*printf("For n=%d the best GLR value is %f\n",n,maxGLR);*/

    /*Find the first time that the GLR increases c_ARL there we stop*/
    if ((maxGLR > c_ARL) && !stop) {
      N = n;
      stop = 1;
      break;
    }
  
    /*Advance counter*/
    n++;
  }

  /*Return value (add 1 for R/SPlus array compability*/
  *ret_N = N+1;
}

/*======================================================================
  GLR in the Epidemic Poisson model
  ======================================================================
*/

/*Helper functions*/

/* Score function */
inline double score(double phi, int *x, double *xm1, double *mu0, int k, int n) {
  register int i;
  double sum = 0;
  /*printf("[1] ");*/
  for (i=k; i<=n; i++) {
    sum += (x[i]*xm1[i])/(exp(phi)*xm1[i]+mu0[i]) - xm1[i];
  }
  /*printf("\n");*/
  return(exp(phi)*sum);
}

inline double sqr(double x) {
  return(x*x);
}

/*fisher information*/
inline double fisher(double phi,int *x,double *xm1, double *mu0, int k,int n,double scorephi) {
  register int i;
  double sum = 0;
  for (i=k; i<=n; i++) {
    sum += (x[i]*sqr(xm1[i]))/sqr(exp(phi)*xm1[i]+mu0[i]);
  }
  return(-scorephi + exp(2.0*phi)*sum);
}

/**********************************************************************
 GLR detector for the epidemic Poisson model described in 
 Held et. al (2005).

 Parameters:
  x   -- the data (as array)
  mu0 -- base means under H0
  lx -- length of x
  Mtilde_R -- number of obs needed to get good estimate (typically 1)
  M        -- Mtilde < M
  xm10 -- observed value of x_0 (0 for initialization, but known if >1st round)
  c_ARL_R -- constant determining when to signal alarm
  ret_N -- the return value
  ret_lr --- GLR value for each n to be returned
**********************************************************************/

void glr_epi_window(int* x,double* mu0, int *lx_R, int *Mtilde_R, int *M_R, double *xm10, double *c_ARL_R,int *ret_N, double *ret_glr) {
  /* printf("====> begin glr_epi\n"); */

  /* Pointers to something useful */
  int lx = *lx_R; /* length of x */
  int Mtilde = *Mtilde_R;
  int M = *M_R;
  double c_ARL = *c_ARL_R;

  /* Loop variables */
  register int n, k,i;

  /* Init return values up to the first position */
  int n0 = fmax(Mtilde-1,0); /*hoehle: 25.9: changepoint can happen at position one: fmax(Mtilde-1,1);*/
  for (n=0; n<n0; n++) { ret_glr[n] = 0; }
  n=n0;

  /* Compute x_{t-1} */
  double xm1[lx];
  xm1[0] = *xm10; /* used to be 0 */
  for (i=1; i<lx; i++) { xm1[i] = x[i-1]; }

  /* Declare variables */
  int stop = 0;
  int N = lx;

  /* printf("Length of the data = %d\n",lx); */
  /* printf("starting at n0 %d\n",n0); */

  /*Loop over all n0 <= n <= length(x)*/
  while ((n < lx)) {
    /*Compute for one n */
    /*printf("n=%d\n",n); */

    /*Define max of the GLR stats*/
    double maxGLR = -1e99;

    /*For the computation of lambda_ml (save old value) */
    double lnk = 0;
    double lambda_ml = 0.5;

    /*Loop over all k  */
    int low = (M==(-1)) ? 0 : fmax(0,n-M);
    int up  = n-(Mtilde-1);
    for (k=low; k<=up; k++) {
      /*printf("n = %d, k = %d\n",n,k); */

      /*Init phi at the last MLE*/
      double phi_old = 2.0;
      double phi_new = log(lambda_ml);
      int iter = 0, maxIter = 1e3;
      /*Compute the MLE (move tuning parameters up as arguments)? */
      /*printf("phi_new = %f\n",phi_new); */
      /*printf("diff = %f\n",fabs(exp(phi_new) - exp(phi_old))); */

      while ((phi_new>-18) & (fabs(exp(phi_new) - exp(phi_old)) > 1e-6) & (iter<maxIter)) {
	iter++;

	/*Do one Newton-Raphson step */
        phi_old = phi_new;

	/*Compute score and fisher functions */
	double scorephi = score(phi_old,x,xm1,mu0,k,n);
        phi_new = phi_old + scorephi/fisher(phi_old,x,xm1,mu0,k,n,scorephi);
	/*printf("score(%f) =  = %f\n",phi_old,scorephi);
	  printf("fisher(%f) =  = %f\n",phi_old,fisher(phi_old,x,xm1,mu0,k,n,scorephi)); 
	  printf("phi_new = %f\n",phi_new); */
      }
      /*Compute the MLE */
      lambda_ml = exp(phi_new);

      /*Compute l_{n,k} (we use no recursive thing for the 2nd term -- not worth it) */
      lnk = 0;
      for (i=k;i<=n;i++) {
	lnk += x[i] * log( lambda_ml * xm1[i]/mu0[i] + 1.0) - lambda_ml*xm1[i];
      }

      /*Debug */
      /*printf("For n=%d, k=%d we have lambda_ml = %f and one.k = %f\n",n,k,lambda_ml,lnk); */

      /*Save the max value */
      if (lnk > maxGLR) { maxGLR = lnk;}
    }

    /*Debug */
    /*printf("For n=%d the best GLR value is %f\n",n,maxGLR); */

    /*Save the return value */
    ret_glr[n] = maxGLR;

    /*Find the first time that the GLR increases c_ARL there we stop */
    if ((maxGLR > c_ARL) && !stop) {
      N = n;
      stop = 1;
      break;
    }
  
    /*Advance counter */
    n++;
  }

  /*Set the remaining values to zero */
  for (i=n+1;i<lx;i++) { ret_glr[i] = 0; }

  /*Return value (add 1 for R/SPlus array compability*/
  *ret_N = N+1;
}


/* Test purposes */
int main( int argc, char *argv[] ) {
  int x[] = { 5,10,10,11,11,8,12,8,13,8,7,7,7,6,4,2,4,7,5,7,6,1,3,2,2,2,1,3,1,1,6,3,2,2,1,2,1,2,3,2,2,4,1,3,5,5,3,6,6,9,5,11,12,4,8,3,8,10,14,12,10,5,8,10,12,7,4,6,4,8,4,3,2,6,1,5,1,1,1,2,1,0,1,3,0,2,1,1,1,1,1,1,2,0,4,1,8,2,3,13,15,8,13,21,12,11,12,10,15,16,20,23,14,15,14,13,9,8,20,10,8,8,6,4,3,6,4,2,6,3,5,3,4,2,2,4,2,3,1,2,3,3,4,1,8,1,7,6,5,9,10,17,6,13,13,12,11,10,12,12,8,8,6,14,7,5,4,7,5,8,4,4,3,5,2,0,1,1,1,2,3,1,2,2,3,2,0,4,3,1,4,2,3,9,4,3,3,7,12,7,10,9,14,12,10,10,8,8,10,19,9,4,9,11,8,6,6,5,5,9,6,5,3,3,2,4,4,3,2,5,1,2,3,2,0,2,1,1,6,2,2,6,3,2,9,4,6,8,6,8};

  double mu0[] = { 
    8.740325,9.264172,9.715987,10.07551,10.32563,
    10.45395,10.45395,10.32563,10.07551,9.715987,
    9.264172,8.740325,8.16617,7.563268,6.951627,
    6.348679,5.768659,5.222371,4.717282,4.257862,
    3.846057,3.481838,3.163735,2.889329,2.655669,
    2.459603,2.298031,2.168088,2.067267,1.993501,
    1.945212,1.921334,1.921334,1.945212,1.993501,
    2.067267,2.168088,2.298031,2.459603,2.655669,
    2.889329,3.163735,3.481838,3.846057,4.257862,
    4.717282,5.222371,5.768659,6.348679,6.951627,
    7.563268,8.16617,8.740325,9.264172,9.715987,
    10.07551,10.32563,10.45395,10.45395,10.32563,
    10.07551,9.715987,9.264172,8.740325,8.16617,
    7.563268,6.951627,6.348679,5.768659,5.222371,
    4.717282,4.257862,3.846057,3.481838,3.163735,
    2.889329,2.655669,2.459603,2.298031,2.168088,
    2.067267,1.993501,1.945212,1.921334,1.921334,
    1.945212,1.993501,2.067267,2.168088,2.298031,
    2.459603,2.655669,2.889329,3.163735,3.481838,
    3.846057,4.257862,4.717282,5.222371,5.768659,
    6.348679,6.951627,7.563268,8.16617,8.740325,
    9.264172,9.715987,10.07551,10.32563,10.45395,
    10.45395,10.32563,10.07551,9.715987,9.264172,
    8.740325,8.16617,7.563268,6.951627,6.348679,
    5.768659,5.222371,4.717282,4.257862,3.846057,
    3.481838,3.163735,2.889329,2.655669,2.459603,
    2.298031,2.168088,2.067267,1.993501,1.945212,
    1.921334,1.921334,1.945212,1.993501,2.067267,
    2.168088,2.298031,2.459603,2.655669,2.889329,
    3.163735,3.481838,3.846057,4.257862,4.717282,
    5.222371,5.768659,6.348679,6.951627,7.563268,
    8.16617,8.740325,9.264172,9.715987,10.07551,
    10.32563,10.45395,10.45395,10.32563,10.07551,
    9.715987,9.264172,8.740325,8.16617,7.563268,
    6.951627,6.348679,5.768659,5.222371,4.717282,
    4.257862,3.846057,3.481838,3.163735,2.889329,
    2.655669,2.459603,2.298031,2.168088,2.067267,
    1.993501,1.945212,1.921334,1.921334,1.945212,
    1.993501,2.067267,2.168088,2.298031,2.459603,
    2.655669,2.889329,3.163735,3.481838,3.846057,
    4.257862,4.717282,5.222371,5.768659,6.348679,
    6.951627,7.563268,8.16617,8.740325,9.264172,
    9.715987,10.07551,10.32563,10.45395,10.45395,
    10.32563,10.07551,9.715987,9.264172,8.740325,
    8.16617,7.563268,6.951627,6.348679,5.768659,
    5.222371,4.717282,4.257862,3.846057,3.481838,
    3.163735,2.889329,2.655669,2.459603,2.298031,
    2.168088,2.067267,1.993501,1.945212,1.921334,
    1.921334,1.945212,1.993501,2.067267,2.168088,
    2.298031,2.459603,2.655669,2.889329,3.163735
  };
  
  int N;
  int n0 = 10;
  int lx = 150;
  int i;
  int dir=1;
  double c_ARL = 5.0;
  double val[150];
  glr_cusum(x,mu0,&lx,&n0,&c_ARL,&N,val,&dir);
  for (i=0;i<150;i++) printf("val[%d]=%f\n",i,val[i]);
  int M = 50;
  int Mtilde = 5;
  glr_cusum_window_old(x,mu0,&lx,&M,&Mtilde,&c_ARL,&N);
  printf("N = %d\n",N);

  return(0);
}


/*Stupid tester */
void foo(double *x) {
  *x = log(10);
}
 
