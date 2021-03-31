/*
 *  16-Aug 2012 / (C) Michael Hoehle
 *  This file is a modified version of the code ks.c available
 *  at http://svn.r-project.org/R/trunk/src/library/stats/src/ks.c (r60102)
 *  The file is copyright 1999-2009 by The R Core Team under GPL-2 
 *  (or later) as shown below. As stated in the GPL-2 license 
 *  the present file is again available under GPL-2.
 *
 *  License:
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* ks.c
   Compute the asymptotic distribution of the one- and two-sample
   two-sided Kolmogorov-Smirnov statistics, and the exact distributions
   in the two-sided one-sample and two-sample cases.
*/

#include <R.h>
#include <Rmath.h>		/* constants */

/*#include "ctest.h"*/

static double K(int n, double d);
static void m_multiply(double *A, double *B, double *C, int m);
static void m_power(double *A, int eA, double *V, int *eV, int m, int n);

/* Two-sample two-sided asymptotic distribution */
void
pkstwo(Sint *n, double *x, double *tol)
{
/* x[1:n] is input and output
 *
 * Compute
 *   \sum_{k=-\infty}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = 1 + 2 \sum_{k=1}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^\infty \exp(-(2k-1)^2\pi^2/(8x^2))
 *
 * See e.g. J. Durbin (1973), Distribution Theory for Tests Based on the
 * Sample Distribution Function.  SIAM.
 *
 * The 'standard' series expansion obviously cannot be used close to 0;
 * we use the alternative series for x < 1, and a rather crude estimate
 * of the series remainder term in this case, in particular using that
 * ue^(-lu^2) \le e^(-lu^2 + u) \le e^(-(l-1)u^2 - u^2+u) \le e^(-(l-1))
 * provided that u and l are >= 1.
 *
 * (But note that for reasonable tolerances, one could simply take 0 as
 * the value for x < 0.2, and use the standard expansion otherwise.)
 *
 */
    double new, old, s, w, z;
    Sint i, k, k_max;

    k_max = (Sint) sqrt(2 - log(*tol));

    for(i = 0; i < *n; i++) {
	if(x[i] < 1) {
	    z = - (M_PI_2 * M_PI_4) / (x[i] * x[i]);
	    w = log(x[i]);
	    s = 0;
	    for(k = 1; k < k_max; k += 2) {
		s += exp(k * k * z - w);
	    }
	    x[i] = s / M_1_SQRT_2PI;
	}
	else {
	    z = -2 * x[i] * x[i];
	    s = -1;
	    k = 1;
	    old = 0;
	    new = 1;
	    while(fabs(old - new) > *tol) {
		old = new;
		new += 2 * s * exp(z * k * k);
		s *= -1;
		k++;
	    }
	    x[i] = new;
	}
    }
}

/* The two-sided one-sample 'exact' distribution */
void
pkolmogorov2x(double *x, Sint *n)
{
    /* x is input and output. */

    *x = K(*n, *x);
}

static double
K(int n, double d)
{
    /* Compute Kolmogorov's distribution.
       Code published in
	 George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
	 "Evaluating Kolmogorov's distribution".
	 Journal of Statistical Software, Volume 8, 2003, Issue 18.
	 URL: http://www.jstatsoft.org/v08/i18/.
    */

   int k, m, i, j, g, eH, eQ;
   double h, s, *H, *Q;

   /* 
      The faster right-tail approximation is omitted here.
      s = d*d*n; 
      if(s > 7.24 || (s > 3.76 && n > 99)) 
          return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
   */
   k = (int) (n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   H = (double*) Calloc(m * m, double);
   Q = (double*) Calloc(m * m, double);
   for(i = 0; i < m; i++)
       for(j = 0; j < m; j++)
	   if(i - j + 1 < 0)
	       H[i * m + j] = 0;
	   else
	       H[i * m + j] = 1;
   for(i = 0; i < m; i++) {
       H[i * m] -= pow(h, i + 1);
       H[(m - 1) * m + i] -= pow(h, (m - i));
   }
   H[(m - 1) * m] += ((2 * h - 1 > 0) ? pow(2 * h - 1, m) : 0);
   for(i = 0; i < m; i++)
       for(j=0; j < m; j++)
	   if(i - j + 1 > 0)
	       for(g = 1; g <= i - j + 1; g++)
		   H[i * m + j] /= g;
   eH = 0;
   m_power(H, eH, Q, &eQ, m, n);
   s = Q[(k - 1) * m + k - 1];
   for(i = 1; i <= n; i++) {
       s = s * i / n;
       if(s < 1e-140) {
	   s *= 1e140;
	   eQ -= 140;
       }
   }
   s *= pow(10., eQ);
   Free(H);
   Free(Q);
   return(s);
}

static void
m_multiply(double *A, double *B, double *C, int m)
{
    /* Auxiliary routine used by K().
       Matrix multiplication.
    */
    int i, j, k;
    double s;
    for(i = 0; i < m; i++)
	for(j = 0; j < m; j++) {
	    s = 0.;
	    for(k = 0; k < m; k++)
		s+= A[i * m + k] * B[k * m + j];
	    C[i * m + j] = s;
	}
}

static void
m_power(double *A, int eA, double *V, int *eV, int m, int n)
{
    /* Auxiliary routine used by K().
       Matrix power.
    */
    double *B;
    int eB , i;

    if(n == 1) {
	for(i = 0; i < m * m; i++)
	    V[i] = A[i];
	*eV = eA;
	return;
    }
    m_power(A, eA, V, eV, m, n / 2);
    B = (double*) Calloc(m * m, double);
    m_multiply(V, V, B, m);
    eB = 2 * (*eV);
    if((n % 2) == 0) {
	for(i = 0; i < m * m; i++)
	    V[i] = B[i];
	*eV = eB;
    }
    else {
	m_multiply(A, B, V, m);
	*eV = eA + eB;
    }
    if(V[(m / 2) * m + (m / 2)] > 1e140) {
	for(i = 0; i < m * m; i++)
	    V[i] = V[i] * 1e-140;
	*eV += 140;
    }
    Free(B);
}
