\name{meanResponse}
\alias{meanResponse}
\encoding{latin1}
\title{Calculate mean response needed in algo.hhh}
\description{
        Calculates the mean response for the model specified in designRes 
        according to equations (1.2) and (1.1) in Held et al., 2005
        for univariate time series and equations (3.3) and (3.2) (with extensions)
        for multivariate time series. See details.
}
\usage{
meanResponse(theta, designRes)
}
\arguments{
	\item{theta}{vector of parameters 
	
			\eqn{\theta = (\alpha_1,\ldots,\alpha_m,
								\bold{\lambda}, 
								\bold{\phi}, \bold{\beta}, 
								\bold{\gamma}_1, \ldots,  \bold{\gamma}_m, 
								\bold{\psi}),}{\theta = (\alpha_1,\ldots,\alpha_m,
								\lambda, \phi, \beta, 
								\gamma_1, \ldots, \gamma_m, \psi),}
								
			where \eqn{\bold{\lambda}=(\lambda_1,\ldots,\lambda_m)}{\lambda=(\lambda_1,\ldots,\lambda_m)}, 
			\eqn{\bold{\phi}=(\phi_1,\ldots,\phi_m)}{\phi=(\phi_1,\ldots,\phi_m)}, 
			\eqn{\bold{\beta}=(\beta_1,\ldots,\beta_m)}{\beta=(\beta_1,\ldots,\beta_m)}, 
			\eqn{\bold{\gamma_1}=(\gamma_{11},\ldots,\gamma_{1,2S_1})}{\gamma_1=(\gamma_11,\ldots,\gamma_(1,2S_1))}, 
			\eqn{\bold{\gamma_m}=(\gamma_{m1},\ldots,\gamma_{m,2S_m})}{\gamma_m=(\gamma_m1,\ldots,\gamma_(m,2S_m))}, 
			\eqn{\bold{\psi}=(\psi_1,\ldots,\psi_m)}{\psi=(\psi_1,\ldots,\psi_m)}.
			
			If the model specifies less parameters, those components are omitted.} 
	\item{designRes}{Result of a call to \code{make.design} }
}

\details{        
	Calculates the mean response for a Poisson or a negative binomial model
	with mean
	\deqn{\mu_t = \lambda y_{t-lag} + \nu_t }{\mu_t = \lambda y_t-lag + \nu_t }
	where  
	\deqn{\log( \nu_t) = \alpha + \beta t + \sum_{j=1}^{S}(\gamma_{2j-1} \sin(\omega_j t) + \gamma_{2j} \cos(\omega_j t) ) }{
			log(\nu_t) = \alpha + \beta t + \sum_(j=1)^S (\gamma_(2j-1) * \sin(\omega_j * t) + \gamma_2j * \cos(\omega_j * t) )  }
	and \eqn{\omega_j = 2\pi j/period }{\omega_j = 2 * \pi * j / period} are Fourier frequencies with
	known period, e.g. \code{period}=52 for weekly data,
	for a univariate time series.
	
	Per default, the number of cases at time point \eqn{t-1}, i.e. \eqn{lag=1}, enter
	as autoregressive covariates into the model. Other lags can also be considered.
	
	The seasonal terms in the predictor can also be expressed as
	\eqn{\gamma_{s} \sin(\omega_s t) + \delta_{s} \cos(\omega_s t) = A_s \sin(\omega_s t + \epsilon_s)}{
	 \gamma_s sin(\omega_s * t) + \delta_s cos(\omega_s * t) = A_s sin(\omega_s * t + \epsilon_s)}
	with amplitude \eqn{A_s=\sqrt{\gamma_s^2 +delta_s^2}}{A_s=sqrt{\gamma_s^2 +delta_s^2}} 
	and phase difference \eqn{\tan(\epsilon_s) = \delta_s / \gamma_s}. The amplitude and
	phase shift can be obtained from a fitted model by specifying \code{amplitudeShift=TRUE}
	in the \code{coef} method.
	
	For multivariate time series the mean structure is
	\deqn{\mu_{it} = \lambda_i y_{i,t-lag} + \phi_i \sum_{j \sim i} w_{ji} y_{j,t-lag} + n_{it} \nu_{it}}{
			\mu_it = \lambda_i * y_i,t-lag + \phi_i * \sum_(j ~ i) w_ji * y_j,t-lag + n_it * \nu_it }
	where 
	\deqn{\log(\nu_{it}) = \alpha_i + \beta_i t + \sum_{j=1}^{S_i} (\gamma_{i,2j-1} \sin(\omega_j t) + \gamma_{i,2j} \cos(\omega_j t) ) }{
			log(\nu_it) = \alpha_i + \beta_i * t + \sum_(j=1)^S_i (\gamma_(i,2j-1) * \sin(\omega_j * t) + \gamma_(i,2j) * \cos(\omega_j * t) )  }
	and \eqn{n_{it}}{n_it} are standardized population counts. The weights \eqn{w_{ji}}{w_ji} are specified in the columns of
	the neighbourhood matrix \code{disProgObj$neighbourhood}.
	
	Alternatively, the mean can be specified as
	\deqn{\mu_{it} = \lambda_i \pi_i y_{i,t-1} + \sum_{j \sim i} \lambda_j (1-\pi_j)/ |k \sim j| y_{j,t-1} + n_{it} \nu_{it}}{
			\mu_it = \lambda_i *\pi_i * y_i,t-1 + \sum_(j ~ i) \lambda_j *(1-\pi_j)/|k ~ j| * y_j,t-1 + n_it * \nu_it }
	if \code{proportion}="single" ("multiple") in \code{designRes$control}. Note that this model specification is still experimental.
	
	
	
}

\value{
	Returns a \code{list} with elements
	\item{mean}{matrix of dimension \eqn{n \times m}{n x m} with the calculated mean 
			response for each time point and unit, where \eqn{n} is the number of time 
			points and \eqn{m} is the number of units. }
	\item{epidemic}{matrix with the epidemic part 
			\eqn{ \lambda_i y_{i,t-1} + \phi_i \sum_{j \sim i} y_{j,t-1}}{
					\lambda_i * y_i,t-1 + \phi_i * \sum_(j ~ i) y_j,t-1} }
	\item{endemic}{matrix with the endemic part of the mean \eqn{ n_{it} \nu_{it} }{ n_it*nu_it } }
	\item{epi.own}{matrix with \eqn{ \lambda_i y_{i,t-1} }{\lambda_i * y_i,t-1} }
	\item{epi.neighbours}{matrix with \eqn{\phi_i \sum_{j \sim i} y_{j,t-1}}{
					\phi_i * \sum_(j ~ i) y_j,t-1} }
}

\author{M. Paul, L. Held}

\keyword{models}

\source{Held, L., \enc{Höhle}{Hoehle}, M., Hofmann, M. (2005) A statistical framework for the 
    analysis of multivariate infectious disease surveillance counts. 
    \emph{Statistical Modelling}, \bold{5}, p. 187--199.}


% make.design
