/*******************************************************************************
 * Call polyCub_iso from polyCubAPI.h for a specific intrfr function
 *
 * Copyright (C) 2017 Sebastian Meyer
 *
 * This file is part of the R package "surveillance",
 * free software under the terms of the GNU General Public License, version 2,
 * a copy of which is available at http://www.R-project.org/Licenses/.
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
    if (d == 1.0) {
        return R - sigma * log1p(R/sigma);
    } else if (d == 2.0) {
        return log1p(R/sigma) - R/(R+sigma);
    } else {
        return (R*pow(R+sigma,1.0-d) - (pow(R+sigma,2.0-d) - pow(sigma,2.0-d))/(2.0-d)) / (1.0-d);
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
    if (d == 1.0) {
        return sigma * logpars[0] * (1.0-logpars[0]/2.0) - log(R+sigma) * (R+sigma) +
            sigma/2.0 * pow(log(R+sigma),2.0) + R;
    } else if (d == 2.0) {
        return (-log(R+sigma) * ((R+sigma)*log(R+sigma) + 2.0*sigma) +
                (R+sigma)*logpars[0]*(logpars[0]+2.0) + 2.0*R) / (R+sigma);
    } else {
        return (pow(sigma,2.0-d) * (logpars[0]*(-d*d + 3.0*d - 2.0) - 2.0*d + 3.0) +
                pow(R+sigma,1.0-d) * (log(R+sigma)*(d-1.0)*(d-2.0) * (R*(d-1.0) + sigma) +
                                      R*(d*d+1.0) + 2.0*d*(sigma-R) - 3.0*sigma)
                ) * d / pow(d-1.0,2.0) / pow(d-2.0,2.0);
    }
}

// student kernel
static double intrfr_student(double R, double *logpars)
{
    double sigma = exp(logpars[0]);
    double d = exp(logpars[1]);
    if (d == 1.0) {
        return log(R*R+sigma*sigma) / 2.0 - logpars[0];
    } else {
        return ( pow(R*R+sigma*sigma,1.0-d) - pow(sigma*sigma,1.0-d) ) / (2.0-2.0*d);
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
    if (d == 1.0) {
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
    double xplint = (twomd == 0.0) ? log(R/sigma) : (pow(R,twomd)-pow(sigma,twomd))/twomd;
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
    return (twomd == 0.0) ? -pow(sigma*log(R/sigma), 2.0) :
        (sigmadRtwomdd * (-twomd)*log(R/sigma) - d*sigma*sigma + sigmadRtwomdd)/(twomd*twomd);
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
    default: error("unknown intrfr_code"); break;
    }
    double center_x = 0.0;
    double center_y = 0.0;
    polyCub_iso(x, y, L, intrfr, pars, &center_x, &center_y,
                subdivisions, epsabs, epsrel, stop_on_error,
                value, abserr, neval);
    return;
}
