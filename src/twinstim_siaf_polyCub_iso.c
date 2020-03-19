/*******************************************************************************
 * Call polyCub_iso from polyCubAPI.h for a specific intrfr function
 *
 * Copyright (C) 2017,2020 Sebastian Meyer
 *
 * This file is part of the R package "surveillance",
 * free software under the terms of the GNU General Public License, version 2,
 * a copy of which is available at https://www.R-project.org/Licenses/.
 ******************************************************************************/

#include <math.h>
#include <polyCubAPI.h>
#include <R_ext/Error.h>


/*** C-implementation of "intrfr" functions ***/

// power-law kernel
static double intrfr_powerlaw(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double d = exp(logpars[1]);
    double onemd = 1.0 - d;
    double twomd = 2.0 - d;
    if (fabs(onemd) < 1e-7) {
        return R - sigma * log1p(R/sigma);
    } else if (fabs(twomd) < 1e-7) {
        return log1p(R/sigma) - R/(R+sigma);
    } else {
        return (R*pow(R+sigma,onemd) - (pow(R+sigma,twomd) - pow(sigma,twomd))/twomd) / onemd;
    }
}

static double intrfr_powerlaw_dlogsigma(double R, double *logpars)
{
    double newlogpars[2] = {logpars[0], log1p(exp(logpars[1]))};
    // sigma*d = exp(logsigma+logd)
    return -exp(logpars[0]+logpars[1]) * intrfr_powerlaw(R, newlogpars);
}

static double intrfr_powerlaw_dlogd(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double d = exp(logpars[1]);
    double onemd = 1.0 - d;
    double twomd = 2.0 - d;
    if (fabs(onemd) < 1e-7) {
        return sigma * logpars[0] * (1.0-logpars[0]/2.0) - log(R+sigma) * (R+sigma) +
            sigma/2.0 * pow(log(R+sigma),2.0) + R;
    } else if (fabs(twomd) < 1e-7) {
        return (-log(R+sigma) * ((R+sigma)*log(R+sigma) + 2.0*sigma) +
                (R+sigma)*logpars[0]*(logpars[0]+2.0) + 2.0*R) / (R+sigma);
    } else {
        return (pow(sigma,twomd) * (logpars[0]*(-d*d + 3.0*d - 2.0) - 2.0*d + 3.0) +
                pow(R+sigma,onemd) * (log(R+sigma)*onemd*twomd * (sigma - R*onemd) +
                                      R*(d*d+1.0) + 2.0*d*(sigma-R) - 3.0*sigma)
                ) * d/onemd/onemd/twomd/twomd;
    }
}

// student kernel
static double intrfr_student(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double d = exp(logpars[1]);
    double onemd = 1.0 - d;
    if (fabs(onemd) < 1e-7) {
        return log(R*R+sigma*sigma) / 2.0 - logpars[0];
    } else {
        return ( pow(R*R+sigma*sigma,onemd) - pow(sigma*sigma,onemd) )/2/onemd;
    }
}

static double intrfr_student_dlogsigma(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double d = exp(logpars[1]);
    return sigma*sigma * ( pow(R*R+sigma*sigma,-d) - pow(sigma,-2.0*d) );
}

static double intrfr_student_dlogd_primitive(double x, double sigma, double d)
{
    double x2ps2 = x*x + sigma*sigma;
    double dm1 = d - 1.0;
    return (d*dm1*log(x2ps2) + d) / (2.0*dm1*dm1 * pow(x2ps2,dm1));
}
static double intrfr_student_dlogd(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double d = exp(logpars[1]);
    if (fabs(d-1.0) < 1e-7) {
        return pow(logpars[0], 2.0) - pow(log(R*R+sigma*sigma), 2.0) / 4.0;
    } else {
        return intrfr_student_dlogd_primitive(R, sigma, d) -
            intrfr_student_dlogd_primitive(0.0, sigma, d);
    }
}

// lagged power-law kernel
static double intrfr_powerlawL_sigmadxplint(double R, double sigma, double d)
{
    double twomd = 2.0 - d;
    double xplint = (fabs(twomd) < 1e-7) ? log(R/sigma) : (pow(R,twomd)-pow(sigma,twomd))/twomd;
    return pow(sigma,d) * xplint;
}

static double intrfr_powerlawL(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double upper = (R > sigma) ? sigma : R;
    double res = upper*upper / 2.0;  // integral over x*constant part
    if (R <= sigma) {
        return res;
    } else {
        return res + intrfr_powerlawL_sigmadxplint(R, sigma, exp(logpars[1]));
    }
}

static double intrfr_powerlawL_dlogsigma(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    if (R <= sigma) {
        return 0.0;
    }
    double d = exp(logpars[1]);
    return d * intrfr_powerlawL_sigmadxplint(R, sigma, d);
}

static double intrfr_powerlawL_dlogd(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    if (R <= sigma) {
        return 0.0;
    }
    double d = exp(logpars[1]);
    double twomd = 2.0 - d;
    double sigmadRtwomdd = pow(sigma,d) * pow(R,twomd) * d;
    return (fabs(twomd) < 1e-7) ? -pow(sigma*log(R/sigma), 2.0) :
        (sigmadRtwomdd * (-twomd)*log(R/sigma) - d*sigma*sigma + sigmadRtwomdd)/(twomd*twomd);
}

// Gaussian kernel
static double intrfr_gaussian(double R, double *logsigma)
{
    double sigma2 = exp(2*logsigma[0]);
    return sigma2 * (1 - exp(-R*R/2/sigma2));
}

static double intrfr_gaussian_dlogsigma(double R, double *logsigma)
{
    double sigma2 = exp(2*logsigma[0]);
    double R2sigma2 = R*R/2/sigma2;
    return 2*sigma2 * (1 - (1+R2sigma2)/exp(R2sigma2));
}

// Exponential kernel
static double intrfr_exponential(double R, double *logsigma)
{
    double sigma = exp(logsigma[0]);
    return sigma * (sigma - (R+sigma)*exp(-R/sigma));
}

static double intrfr_exponential_dlogsigma(double R, double *logsigma)
{
    double sigma = exp(logsigma[0]);
    return 2*sigma*sigma - ((R+sigma)*(R+sigma) + sigma*sigma)*exp(-R/sigma);
}


/*** function to be called from R ***/

void C_siaf_polyCub1_iso(
    double *x, double *y,  // vertex coordinates (open)
    int *L,                // number of vertices
    int *intrfr_code,      // F(R) identifier
    double *pars,          // parameters for F(R)
    int *subdivisions, double *epsabs, double *epsrel, // Rdqags options
    int *stop_on_error,
    double *value, double *abserr, int *neval) // results
{
    intrfr_fn intrfr;
    switch(*intrfr_code) { // = INTRFR_CODE in ../R/twinstim_siaf_polyCub_iso.R
    case 10: intrfr = intrfr_powerlaw; break;
    case 11: intrfr = intrfr_powerlaw_dlogsigma; break;
    case 12: intrfr = intrfr_powerlaw_dlogd; break;
    case 20: intrfr = intrfr_student; break;
    case 21: intrfr = intrfr_student_dlogsigma; break;
    case 22: intrfr = intrfr_student_dlogd; break;
    case 30: intrfr = intrfr_powerlawL; break;
    case 31: intrfr = intrfr_powerlawL_dlogsigma; break;
    case 32: intrfr = intrfr_powerlawL_dlogd; break;
    case 40: intrfr = intrfr_gaussian; break;
    case 41: intrfr = intrfr_gaussian_dlogsigma; break;
    case 50: intrfr = intrfr_exponential; break;
    case 51: intrfr = intrfr_exponential_dlogsigma; break;        
    default: error("unknown intrfr_code"); break;
    }
    double center_x = 0.0;
    double center_y = 0.0;
    polyCub_iso(x, y, L, intrfr, pars, &center_x, &center_y,
                subdivisions, epsabs, epsrel, stop_on_error,
                value, abserr, neval);
    return;
}
