/**

  C routines for the surveillance package

  Author: (C) Michael Höhle <http://www.stat.uni-muenchen.de/~hoehle>
  Date:   8 Jan 2008

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, a copy is available at
  http://www.r-project.org/Licenses/

  Atm the only C routines are concerned with the GLR computations
  in the algorithm algo.prc.

  //should check that these really work...
  void lr_cusum  - intercept chart with known kappa
  void glr_cusum - intercept chart with estimated kappa
  void glr_cusum_window -- window limited intercept chart with estimated kappa

  //removedvoid glr_epi
  void glr_epi_window

  //History
  17 Feb 2009 -- added LR scheme for negative binomial (still experimental)
  08 Jan 2007 -- added the files for the negative binomial computations
  21 Sep 2007 -- modified code to get around error of extreme strict (=pedantic) MacOS compiling on CRAN
  28 Nov 2006 -- file created
*/

/*#define DEBUG*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

/* header */
/* void lr_cusum(int* ,double* , int *, double *, double *,int *, double *) ;
void glr_cusum(int* ,double* , int *, int *, double *,int *, double *, int *, int *, int *) ; */

/* Helper function for x^2 */
static R_INLINE double sqr(double x) {
  return(x*x);
}


/*======================================================================
  Poisson GLR detector
  ======================================================================
*/


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
  ret_cases - The number of cases to be returned
  ret - what should be returned (value of lr-statistic, cases)?
**********************************************************************/

void lr_cusum(int* x,double* mu0, int *lx_R, double *kappa_R, double *c_ARL_R,int *ret_N, double *ret_lr, double *ret_cases, int *ret_R) {
  /* Pointers to something useful */
  int lx = *lx_R;
  double c_ARL = *c_ARL_R;
  double kappa = *kappa_R;
  int ret = *ret_R;

  /* Loop variables */
  register int n=0;
  int stop = 0;
  int N = lx;

  /* Loop over all 0 <= n <= length(x) */
  while ((n < lx)) {
    /*Compute for one n*/
    /*printf("n=%d\n",n);*/


    double zn = kappa * x[n] + (1-exp(kappa))*mu0[n];

#ifdef DEBUG
    printf("For kappa=%f and mu[%d]=%f:\tx[%d]=%f, LR=%f\n",kappa,n,mu0[n],n,x[n],zn);
#endif

    /* Add up */
    if (n==0) {
      ret_lr[n] = fmax(0,zn);
      /*5.11.2009 -- Bug fix. There was a small programming error for the
	computing the cases for n==0.
	if (ret==2) ret_cases[n] = (c_ARL + mu0[n]*(kappa-1))/kappa ; */
      if (ret==2) ret_cases[n] = (c_ARL + mu0[n]*(exp(kappa)-1))/kappa ;
    }
    else {
      ret_lr[n] = fmax(0,ret_lr[n-1] + zn);
      if (ret==2) ret_cases[n] = (c_ARL - ret_lr[n-1] + mu0[n]*(exp(kappa)-1))/kappa ;
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



/***********************************************************************
 Function for the computation of the glr-statistic with time-varying
 in-control value

 Params
 n   - timepoint n where the glr-statistic should be computed
 x   - array with observations
 mu0 - array with estimated in-comtrol parameters
 dir - direction of testing (up (1) or down (-1)

 the function returns max_1<=k<=n sup_theta sum_t=k^n log f_theta(x_t)/f_theta0(x_t)

************************************************************************/

double glr (int n, int x[], double mu0[], int dir){

    /* For the recursive computation of kappa_ml */
    double sumx = 0;
    double summu0 = 0;

    /* Define max of the GLR stats */
    double maxGLR = -1e99;

    /* Loop variable */
    register int k;

    /* For fitting and summation */
    double kappa_ml = 0;
    double sum = 0;

    /* Loop over all k */
    for (k=n; k>=0; k--) { /* Backwards loop makes calculations faster */
      /* Recursive update of the kappa.ml quantitities */

      sumx += x[k];
      summu0 += mu0[k];

      /* Calculate MLE of kappa */
      kappa_ml = dir*fmax(0,dir*log(sumx/summu0));

      /* Recursive updating of the likelihood ratios -- See
	 notes on the 21 september printout. This is fast! */
      sum = kappa_ml * sumx + (1-exp(kappa_ml))*summu0;

      /* save max value */
      if (sum > maxGLR) { maxGLR = sum;}
    }
    return(maxGLR);
}


/***********************************************************************
 Function for the computation of the window-limited glr-statistic with time-varying
 in-control value

 Params
 n   - timepoint n where the glr-statistic should be computed
 x   - array with observations
 mu0 - array with estimated in-comtrol parameters
 dir - direction of testing (up (1) or down (-1)
 M - max time to go back in time from N
 Mtilde  - number of vals we will need to estimate a detection

 the function returns max(0,n-M) <= k <= n-Mtilde sup_theta sum_t=k^n log f_theta(x_t)/f_theta0(x_t)

************************************************************************/

double glr_window (int n, int x[], double mu0[], int dir, int M, int Mtilde){

/* Define max of the GLR stats */
    double maxGLR = -1e99;

    /* Loop variable */
    register int k,l;

    /* For the recursive computation of kappa_ml compute for (n-Mtilde+1):n */
    double sumx = 0;
    double summu0 = 0;
    /* For fitting and summation */
    double sum = 0;
    double kappa_ml = 0;

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
      kappa_ml = dir*fmax(0,dir*log(sumx/summu0));;

      /*Calculate sum of likelihood ratios using recursive updating (fast!)*/
      sum = kappa_ml * sumx + (1-exp(kappa_ml))*summu0;

      /* Save the max value */
      if (sum > maxGLR) { maxGLR = sum;}
    }

    return(maxGLR);
}


/**********************************************************************
  Fast C implementation of the sequential GLR test without windowing
  for Poisson distributed variables, this function can test in both
  directions (up/down) and there is the possibility ( in opposite to old function
  glr_cusum) to return the number of cases at timepoint n to produce an alarm at
  any timepoint 1<=k<=n

  Params:
   x   - array of observed values (pos 0 is first interesting value)
  mu0  - array with the means once in-control (pos 0 is first interesting value)
  lx   - length of the x and mu0 array
  n0   - number of burn-in values (number of observations, not array index!)
  c_ARL- when to sound alarm threshold
  ret_N- here the return value is stored
  ret_glr- GLR value for each n to be returned
  dir  - direction of testing
  ret - what should be returned (value of glr-statistic, cases)?

**********************************************************************/

void glr_cusum(int* x,double* mu0, int *lx_R, int *n0_R, double *c_ARL_R,int *ret_N, double *ret_glr, double *ret_cases, int *dir_R, int *ret_R) {

  /* Pointers to something useful */
  int lx = *lx_R;
  int n0 = *n0_R;
  int dir = *dir_R;
  int ret = *ret_R;
  double c_ARL = *c_ARL_R;

  /* Loop variables */
  register int n; /*l,n0-1*/

  for (n=0; n<n0-1; n++) { ret_glr[n] = 0; }
  for (n=0; n<n0-1; n++) { ret_cases[n] = 0; }

  int stop = 0;
  int N = lx;


  /* Precalculation of log(mu0) -- apparently not used anymore*/
  //double logmu0[lx];
  //for (l=0;l<lx;l++) {logmu0[l] = log(mu0[l]); }

  /* Debug */
  /* for (l=0;l<lx;l++) {printf("logmu0[%d] = log(%f) = %f\n",l,mu0[l],logmu0[l]); } */

  /* Loop over all n0 <= n <= length(x) */
  while ((n < lx)) {
    /* Compute for one n */


      /* to compute the glr-statistic with helper function (see above) */
     ret_glr[n] = glr(n,x,mu0,dir);


      /* to find the number of cases that are necassary to produce an alarm */
      /* optionally, if ret == 2*/
      if (ret == 2){
        /* change the value at timepoint n until an alarm is produced */

        int xnnew = -1;

        /* glr-statistic for the new x value, initialize it so, that the loop starts */
        double glrnew = c_ARL - dir;

        /* save the old value of x */
        int xnold = x[n];

        /* increase/decrease xnnew until the glr-statistic with the new x is >= c_ARL */
        while ((dir*glrnew < c_ARL*dir)){

          /* increase/decrease xnnew */
          xnnew = xnnew + 1;
          /* put this value in vector x at timepoint n */
          x[n] = xnnew;

          /* compute the glr-statistic */
          glrnew = glr(n,x,mu0,dir);
        }

       /* save the value */
       ret_cases[n] = xnnew;
       /* set x[n] back to original value so that we can go to next n*/
       x[n] = xnold;

      }


      /* Find the first time that the GLR increases c_ARL there we stop */
      if ((ret_glr[n] >= c_ARL) && !stop) {
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

void glr_cusum_window(int* x,double* mu0, int *lx_R, int *M_R, int *Mtilde_R, double *c_ARL_R,int *ret_N, double *ret_glr, double *ret_cases, int *dir_R, int *ret_R) {
  /* Pointers to something useful */
  int lx = *lx_R;
  int M = *M_R;
  int Mtilde = *Mtilde_R;
  int dir = *dir_R;
  int ret = *ret_R;
  double c_ARL = *c_ARL_R;

  /* Loop variables (n>Mtilde, so we start with n=Mtilde (due to -1 in index) */
  register int n = Mtilde; /*l*/
  int stop = 0;
  int N = lx;

  /* Precalculation of log(mu0) -- apparently not used anymore */
  //double logmu0[lx];
  //for (l=0;l<lx;l++) { logmu0[l] = log(mu0[l]); }

  /* Debug */
  /* for (l=0;l<lx;l++) {printf("logmu0[%d] = log(%f) = %f\n",l,mu0[l],logmu0[l]); } */


  /* Loop over all n0 <= n <= length(x) */
  while ((n < lx)) {
    /* Compute for one n */


      /* to compute the glr-statistic with helper function (see above) */
     ret_glr[n] = glr_window(n,x,mu0,dir,M,Mtilde);


      /* to find the number of cases that are necassary to produce an alarm */
      /* optionally, if ret == 2*/
      if (ret == 2){
        /* change the value at timepoint n as long as an alarm is produced */

        int xnnew = -1;

        /* glr-statistic for the new x value, initialize it so, that the loop starts */
        double glrnew = c_ARL - dir;

        /* save the old value of x */
        int xnold = x[n];


        /* increase/decrease xnnew as long the glr-statistic with the new x is >= c_ARL */
        while ((dir*glrnew < c_ARL*dir)){

          /* increase/decrease xnnew */
          xnnew = xnnew + 1;
          /* put this value in vector x at timepoint n */
          x[n] = xnnew;

          /* compute the glr-statistic */
          glrnew = glr_window(n,x,mu0,dir,M,Mtilde);
        }

       /* save the value */
       ret_cases[n] = xnnew;
       /* set x[n] back to original value so that we can go to next n*/
       x[n] = xnold;

      }

    /* Debug*/
    /* printf("For n=%d the best GLR value is %f\n",n,maxGLR);*/

    /* Find the first time that the GLR increases c_ARL there we stop */
    if ((ret_glr[n] >= c_ARL) && !stop) {
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



/*======================================================================
  GLR in the Epidemic Poisson model
  ======================================================================
*/

/*Helper functions*/

/* Score function */
static R_INLINE double score(double phi, int *x, double *xm1, double *mu0, int k, int n) {
  register int i;
  double sum = 0;
  /*printf("[1] ");*/
  for (i=k; i<=n; i++) {
    sum += (x[i]*xm1[i])/(exp(phi)*xm1[i]+mu0[i]) - xm1[i];
  }
  /*printf("\n");*/
  return(exp(phi)*sum);
}

/*fisher information*/
static R_INLINE double fisher(double phi,int *x,double *xm1, double *mu0, int k,int n,double scorephi) {
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


/*
  ======================================================================
                              Negative binomial chart


  Comment/ToDo: move to seperate files?
  ======================================================================
*/

/**********************************************************************
  C implementation of the LR test for the negative binomial based
  intercept shift chart with known, but possible time varying, incontrol
  mean. See Hoehle and Paul (2008) paper for details.

  Params:
   x   - array of observed values (pos 0 is first interesting value)
  mu0  - array with the means once in-control (pos 0 is first interesting value)
  alpha -- fixed dispersion parameter of the NegBin distribution (see Lawless87)
  lx   - length of the x and mu0 array
  kappa- the change in intercept to detect (here known in advance)

  c_ARL- when to sound alarm threshold
  ret_N- here the return value is stored
  ret_lr- GLR value for each n to be returned
  ret_cases - The number of cases to be returned
  ret - what should be returned (value of lr-statistic, cases)?
**********************************************************************/

void lr_cusum_nb(int* x, double* mu0, double* alpha_R, int *lx_R, double *kappa_R, double *c_ARL_R,int *ret_N, double *ret_lr, double *ret_cases, int *ret_R) {

#ifdef DEBUG
  printf("====> begin lr_cusum_nb\n");
#endif

  /* Pointers to something useful */
  int lx = *lx_R;
  double c_ARL = *c_ARL_R;
  double kappa = *kappa_R;
  double alpha = *alpha_R;
  int ret = *ret_R;

#ifdef DEBUG
  printf("lx = %d\n",lx);
  printf("alpha = %f\n",alpha);
#endif

  /* Loop variables */
  register int n=0;
  int stop = 0;
  int N = lx;

  /* Loop over all 0 <= n <= length(x) */
  while ((n < lx)) {
    /*Compute for one n*/
#ifdef DEBUG
    printf("n=%d\n",n);
#endif

    /* LR for one NB variable as given in the first equation of Sect 2.1
       in the Hoehle and Paul (2008) paper
    */
    double zn = kappa * x[n] + (x[n]+1/alpha)*log( (1+alpha*mu0[n])/(1+alpha*mu0[n]*exp(kappa)) );

    /* Recursive CUSUM as given in (4) by Hoehle and Paul (2008) */
    if (n==0) {
      /* Statistic */
      ret_lr[n] = fmax(0,zn);
      /* Number of cases it takes to sound an alarm - backcalc'ed by backcalc.mws*/
      if (ret==2) ret_cases[n] = -(log((1+alpha*mu0[n])/(1+alpha*mu0[n]*exp(kappa)))-c_ARL*alpha)/alpha/(kappa+log((1+alpha*mu0[n])/(1+alpha*mu0[n]*exp(kappa))));
    }
    else {
      /* Statistic */
      ret_lr[n] = fmax(0,ret_lr[n-1] + zn);
      /* Number of cases it takes to sound an alarm -- backcalc.mws*/
      if (ret==2) ret_cases[n] = -(ret_lr[n-1]*alpha+log((1+alpha*mu0[n])/(1+alpha*mu0[n]*exp(kappa)))-c_ARL*alpha)/alpha/(kappa+log((1+alpha*mu0[n])/(1+alpha*mu0[n]*exp(kappa))));
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


/* ======================================================================

                    Functions for the intercept chart

   ======================================================================
*/

/* Score function for intercept chart*/
static R_INLINE double nbScore(double kappa, int *x, double *mu0, double alpha, int k, int n) {
  register int i;
  double sum = 0;
  /*printf("[1] ");*/
  for (i=k; i<=n; i++) {
    sum += (x[i]-exp(kappa)*mu0[i])/(1+alpha*exp(kappa)*mu0[i]);
  }
  /*printf("\n");*/
  return(sum);
}


/*fisher information for intercept chart -- its minus the hesse */
static R_INLINE double nbFisher(double kappa,int *x, double *mu0, double alpha, int k,int n) {
  register int i;
  double sum = 0;
  for (i=k; i<=n; i++) {
      sum += mu0[i]*(alpha*x[i]+1)/sqr(1+alpha*exp(kappa)*mu0[i]);
  }
  return( exp(kappa)*sum);
}

/* Formula to compute a single l_{n,k} for the intercept chart */
static R_INLINE double nblnk(double kappa,int *x, double *mu0, double alpha, int k,int n) {
  register int i;
  double lnk = 0;

  for (i=k;i<=n;i++) {
    lnk += kappa * x[i] + (x[i] + 1/alpha) * log( (1+alpha*mu0[i])/(1+alpha*mu0[i]*exp(kappa)));
  }
  return(lnk);
}


/**********************************************************************
 GLR detector for the negative binomial model described in
 Hoehle and Paul (2008).

 Parameters:
  x   -- the data (as array)
  mu0 -- base means under H0
  alpha -- fixed dispersion parameter of the NegBin distribution (see Lawless (1987))
  lx -- length of x
  Mtilde_R -- number of obs needed to get good estimate (typically 1)
  M        -- Mtilde < M
  c_ARL_R -- constant determining when to signal alarm
  ret_N -- the return value
  ret_lr --- GLR value for each n to be returned
**********************************************************************/

void glr_nb_window(int* x,double* mu0, double* alpha_R, int *lx_R, int *Mtilde_R, int *M_R, double *c_ARL_R,int *ret_N, double *ret_glr, int *dir_R) {
#ifdef DEBUG
  printf("====> begin glr_nb_window\n");
#endif
  /* Pointers to something useful */
  int lx = *lx_R; /* length of x */
  int Mtilde = *Mtilde_R;
  int M = *M_R;
  double c_ARL = *c_ARL_R;
  double alpha = *alpha_R;
  int dir = *dir_R;

  /* Loop variables */
  register int n, k,i;

  /*changepoint can happen at position one (ie. index zero in C*/
  int n0 = fmax(Mtilde-1,0);

#ifdef DEBUG
  printf("Length of the data = %d\n",lx);
  printf("starting at n0= %d\n",n0);
#endif
  /* Show the data */
  /*for (n=0; n<lx; n++) { printf("x[%d] = %d\n",n,x[n]); }*/

  /* Init return values up to the first position */
  for (n=0; n<n0; n++) { ret_glr[n] = 0; }
  n=n0;

  /* Declare variables */
  int stop = 0;
  int N = lx;


  /*Loop over all n0 <= n <= length(x)*/
  while ((n < lx)) {
    /*Compute for one n */
    /*printf("n=%d\n",n); */

    /*Define max of the GLR stats*/
    double maxGLR = -1e99;

    /*For the computation of kappa_ml (save old value) */
    double lnk = 0;
    double kappa_ml = 0.5; /*starting value */

    /*Loop over all k  */
    int low = (M==(-1)) ? 0 : fmax(0,n-M);
    int up  = n-(Mtilde-1);
#ifdef DEBUG
    printf("M= %d, Mtilde= %d, low = %d, up = %d\n", M, Mtilde, low, up);
#endif

    for (k=low; k<=up; k++) {
#ifdef DEBUG
      printf("n = %d, k = %d\n",n,k);
#endif
      /*Init kappa at the last MLE*/
      double kappa_old = 0.4;
      double kappa_new = kappa_ml;
      int iter = 0, maxIter = 1e3;
      /*Compute the MLE (move tuning parameters up as arguments)? */
      /*printf("kappa_new = %f\n",kappa_new);
	printf("diff = %f\n",fabs(exp(kappa_new) - exp(kappa_old))); */

      while ((kappa_new>-18) & (fabs(kappa_new - kappa_old) > 1e-6) & (iter<maxIter)) {
	iter++;

	/*Do one Newton-Raphson step */
        kappa_old = kappa_new;
        kappa_new = kappa_old + nbScore(kappa_old,x,mu0,alpha,k,n)/nbFisher(kappa_old,x,mu0,alpha,k,n);

	/* Debug info */
	/*printf("score(%f) =  = %f\n",kappa_old,nbScore(kappa_old,x,mu0,alpha,k,n));
	printf("fisher(%f) =  = %f\n",kappa_old,nbFisher(kappa_old,x,mu0,alpha,k,n));
	printf("kappa_new = %f\n",kappa_new);
	printf("diff = %f\n",fabs(exp(kappa_new) - exp(kappa_old))); */
      }
      /*Compute the MLE */
      kappa_ml = dir*fmax(0,dir*kappa_new);

      /*Compute l_{n,k} */
      lnk = nblnk(kappa_ml, x,mu0,alpha, k, n);

      /*Debug */
#ifdef DEBUG
      printf("For n=%d, k=%d we have kappa_ml = %f and l_{n,k} = %f\n",n,k,kappa_ml,lnk);
#endif
      /*Save the max value */
      if (lnk > maxGLR) { maxGLR = lnk;}
    }

    /*Debug */
#ifdef DEBUG
    printf("For n=%d the highest GLR value is %f\n",n,maxGLR);
#endif
    /*Save the return value */
    ret_glr[n] = maxGLR;

    /*Find the first time that the GLR increases c_ARL there we stop */ /*hoehle: now >= */
    if ((maxGLR >= c_ARL) && !stop) {
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



/* ======================================================================

              Functions for the general negative binomial chart
    	      Current implementation: Epidemic chart

   ======================================================================
*/

/********** Epidemic Chart ***********/

/* alternative \mu_{1,t}(theta) */
static R_INLINE double mu1(int i, double theta, double *mu0, double *xm1) {
  return( mu0[i] + exp(theta) * xm1[i]);
}
/* first derivative */
/*
static R_INLINE double d1mu1(int i, double theta, double *mu0, double *xm1) {
  return( exp(theta) * xm1[i]);
}
*/
/* second derivative */
/*
static R_INLINE double d2mu1(int i, double theta, double *mu0, double *xm1) {
  return( exp(theta) * xm1[i]);
}
*/
/************************************/

/********** Intercept Chart (only upwards) ***********/

/* /\* alternative \mu_{1,t}(theta) *\/ */
/* static R_INLINE double mu1(int i, double theta, double *mu0, double *xm1) { */
/*   return( mu0[i] * exp(exp(theta))) ; */
/* } */
/* /\* first derivative *\/ */
/* static R_INLINE double d1mu1(int i, double theta, double *mu0, double *xm1) { */
/*   return( mu0[i] * exp(theta + exp(theta))); */
/* } */
/* /\* second derivative *\/ */
/* static R_INLINE double d2mu1(int i, double theta, double *mu0, double *xm1) { */
/*   return( mu0[i] * ( exp(theta + exp(theta)) + exp(2*theta + exp(theta)))); */
/* } */

/********** Intercept Chart (only upwards) regular parameterization requiring a fmax***********/

/* alternative \mu_{1,t}(theta) */
/* static R_INLINE double mu1(int i, double theta, double *mu0, double *xm1) { */
/*   return( mu0[i] * exp(theta)) ; */
/* } */
/* /\* first derivative *\/ */
/* static R_INLINE double d1mu1(int i, double theta, double *mu0, double *xm1) { */
/*   return( mu0[i] * exp(theta)); */
/* } */
/* /\* second derivative *\/ */
/* static R_INLINE double d2mu1(int i, double theta, double *mu0, double *xm1) { */
/*   return( mu0[i] * exp(theta)); */
/* } */

/************************************/



/* Score function for the general negative binomial chart*/
/*
static R_INLINE double nbGeneralScore(double theta, int *x, double *xm1, double *mu0, double alpha, int k, int n) {
  register int i;
  double sum = 0;
  double mu1i = 0, d1mu1i = 0;
  for (i=k; i<=n; i++) {
    mu1i = mu1(i, theta, mu0, xm1);
    d1mu1i = d1mu1(i, theta, mu0, xm1);
    sum += d1mu1i * (x[i]-mu1i) / (mu1i * (1+alpha*mu1i));
  }
  return(sum);
}
*/

/*fisher information for the general chart -- its minus the hesse */
/*
static R_INLINE double nbGeneralFisher(double theta,int *x, double *xm1, double *mu0, double alpha, int k,int n) {
  register int i;
  double sum = 0;
  double mu1i = 0, d1mu1i = 0, d2mu1i;
  for (i=k; i<=n; i++) {
    mu1i = mu1(i, theta, mu0, xm1);
    d1mu1i = d1mu1(i, theta, mu0, xm1);
    d2mu1i = d2mu1(i, theta, mu0, xm1);

    sum += ((d2mu1i*(x[i]-mu1i) - sqr(d1mu1i))/ (mu1i * (1+alpha*mu1i))) -
      sqr(d1mu1i)*(x[i]-mu1i)*(1+2*alpha*mu1i)/sqr(mu1i * (1+alpha*mu1i));
  }
  return( sum );
}
*/

/* Formula to compute a single l_{n,k} for the general chart */
static R_INLINE double nbGeneralLnk(double theta,int *x, double *xm1,  double *mu0, double alpha, int k,int n) {
  register int i;
  double lnk = 0, mu1i=0;

  for (i=k;i<=n;i++) {
    mu1i = mu1(i, theta, mu0, xm1);
    lnk += x[i]*( log(mu1i)-log(mu0[i])+log(1+alpha*mu0[i])-log(1+alpha*mu1i)) +
      1/alpha*(log(1+alpha*mu0[i])-log(1+alpha*mu1i));
  }
  return(lnk);
}


/**********************************************************************
 GLR detector for the general negative binomial model described in
 Hoehle and Paul (2008). Currently, the epidemic chart is implemented.
 To obtain other alternatives modify the mu1, d1mu1 and d2mu1 functions.

 Parameters:
  x   -- the data (as array)
  mu0 -- base means under H0
  alpha -- fixed dispersion parameter of the NegBin distribution (see Lawless (1987))
  lx -- length of x
  Mtilde_R -- number of obs needed to get good estimate (typically 1)
  M        -- Mtilde < M
  c_ARL_R -- constant determining when to signal alarm
  ret_N -- the return value
  ret_lr --- GLR value for each n to be returned
**********************************************************************/

void glr_nbgeneral_window(int* x,double* mu0, double* alpha_R, int *lx_R, int *Mtilde_R, int *M_R, double *xm10, double *c_ARL_R,int *ret_N, double *ret_glr, int *dir_R) {
#ifdef DEBUG
  printf("====> begin glr_nbgeneral_window \n");
#endif

  /* Pointers to something useful */
  int lx = *lx_R; /* length of x */
  int Mtilde = *Mtilde_R;
  int M = *M_R;
  double c_ARL = *c_ARL_R;
  double alpha = *alpha_R;
  /* int dir = *dir_R; -- currently direction is not supported?? */

  /* Loop variables */
  register int n, k,i;

  /*changepoint can happen at position one (ie. index zero in C*/
  int n0 = fmax(Mtilde-1,0);

  /* Compute x_{t-1} */
  double xm1[lx];
  xm1[0] = *xm10; /* used to be 0 */
  for (i=1; i<lx; i++) { xm1[i] = x[i-1]; }

#ifdef DEBUG
  printf("Length of the data = %d\n",lx);
  printf("starting at n0= %d\n",n0);
#endif
  /* Show the data */
  /*for (n=0; n<lx; n++) { printf("x[%d] = %d\n",n,x[n]); }*/

  /* Init return values up to the first position */
  for (n=0; n<n0; n++) { ret_glr[n] = 0; }
  n=n0;

  /* Declare variables */
  int stop = 0;
  int N = lx;

  /*Loop over all n0 <= n <= length(x)*/
  while ((n < lx)) {
    /*Compute for one n */
#ifdef DEBUG
    printf("n=%d\n",n);
#endif

    /*Define max of the GLR stats*/
    double maxGLR = -1e99;

    /*For the computation of theta_ml (save old value) */
    double lnk = 0;
    double theta_ml = 0.5; /*starting value */

    /*Loop over all k  */
    int low = (M==(-1)) ? 0 : fmax(0,n-M);
    int up  = n-(Mtilde-1);
#ifdef DEBUG
    printf("M= %d, Mtilde= %d, low = %d, up = %d\n", M, Mtilde, low, up);
#endif

    for (k=low; k<=up; k++) {
#ifdef DEBUG
      printf("n = %d, k = %d\n",n,k);
#endif
      /*Init theta at the last MLE*/
      double theta_old = 0.4;
      double theta_new = theta_ml;
      int iter = 0, maxIter = 1e3;
      /*Compute the MLE (move tuning parameters up as arguments)? */
      /*printf("theta_new = %f\n",theta_new);
	printf("diff = %f\n",fabs(exp(theta_new) - exp(theta_old))); */

      while ((theta_new>-18) & (fabs(theta_new - theta_old) > 1e-6) & (iter<maxIter)) {
	iter++;

	/*Do one Newton-Raphson step */
        theta_old = theta_new;
        theta_new = theta_old + nbScore(theta_old,x,mu0,alpha,k,n)/nbFisher(theta_old,x,mu0,alpha,k,n);

	/* Debug info */
	/*printf("score(%f) =  = %f\n",theta_old,nbScore(theta_old,x,mu0,alpha,k,n));
	printf("fisher(%f) =  = %f\n",theta_old,nbFisher(theta_old,x,mu0,alpha,k,n));
	printf("theta_new = %f\n",theta_new);
	printf("diff = %f\n",fabs(exp(theta_new) - exp(theta_old)));*/
      }
      /*Compute the MLE */
      /*theta_ml = theta_new ; */
      theta_ml = theta_new; /*dir*fmax(0,dir*theta_new);*/

      /*Compute l_{n,k} */
      lnk = nbGeneralLnk(theta_ml, x, xm1, mu0, alpha, k, n);

      /*Debug */
#ifdef DEBUG
      printf("For n=%d, k=%d we have theta_ml = %f and l_{n,k} = %f\n",n,k,theta_ml,lnk);
#endif
      /*Save the max value */
      if (lnk > maxGLR) { maxGLR = lnk;}
    }

    /*Debug */
#ifdef DEBUG
    printf("For n=%d the highest GLR value is %f\n",n,maxGLR);
#endif
    /*Save the return value */
    ret_glr[n] = maxGLR;

    /*Find the first time that the GLR increases c_ARL there we stop */ /*hoehle: now >= */
    if ((maxGLR >= c_ARL) && !stop) {
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
/*
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
  int ret=1;
  double c_ARL = 5.0;
  double val[150];
  double cases[150];

  glr_cusum(x,mu0,&lx,&n0,&c_ARL,&N,val,cases,&dir,&ret);
  for (i=0;i<150;i++) printf("val[%d]=%f\n",i,val[i]);
  int M = 50;
  int Mtilde = 5;
  printf("N = %d\n",N);

  return(0);
}
*/

/*Stupid tester */
void foo(double *x) {
  *x = log(10);
}

