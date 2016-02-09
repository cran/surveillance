\name{create.grid}
\alias{create.grid}
\encoding{latin1}
\title{Create a Matrix of Initial Values for \code{algo.hhh.grid}}
\description{
  (An auxiliary function for the \strong{deprecated} HHH estimation
  routine \code{\link{algo.hhh.grid}}; use \code{\link{hhh4}} instead.)
  
       For a given model and a list of parameters specified as \code{param = c(lower,upper,length)},
       \code{create.grid} creates a grid of initial values for \code{algo.hhh.grid}.
       The resulting matrix contains all combinations of the supplied parameters 
       which each are a sequence of length \code{length} from \code{lower} to 
       \code{upper}.  
       Note that the autoregressive parameters \eqn{\lambda, \phi} and the
       overdispersion parameter \eqn{\psi} must be positive.  
       Only one sequence of initial values is considered for the autregressive,
       endemic and overdispersion parameters to create the grid,
       e.g. initial values are the same
       for each one of the seasonal and trend parameters. 
}
\usage{
create.grid(disProgObj, control, params = list(epidemic = c(0.1, 0.9, 5),
   endemic=c(-0.5,0.5,3), negbin = c(0.3, 12, 10)))
}
\arguments{
   \item{disProgObj}{object of class \code{disProg} }
   \item{control}{specified model}
	\item{params}{list of parameters: \code{param=c(lower,upper,length)}
		\itemize{
		\item \code{epidemic} autoregressive parameters \eqn{\lambda} and \eqn{\phi}. 
		\item \code{endemic} trend and seasonal parameters \eqn{\beta, \gamma_j}. 
		\item \code{negbin} overdispersion parameter for negative binomial model \eqn{\psi}. 
		} 
        }
}
\value{
  \item{matrix}{matrix with \code{gridSize} starting values as rows}     
}
\seealso{\code{\link{algo.hhh.grid}}}
\author{M. Paul}
\examples{
# simulate data
set.seed(123)
disProgObj <- simHHH(control = list(coefs = list(alpha =-0.5, gamma = 0.4,
                     delta = 0.6)),length=300)$data

# consider the model specified in a control object for algo.hhh.grid
cntrl1 <- list(lambda=TRUE, neighbours=TRUE, 
               linear=TRUE, nseason=1)
cntrl2 <- list(lambda=TRUE, negbin="single")

# create a grid of initial values for respective parameters
grid1 <- create.grid(disProgObj, cntrl1, 
                     params = list(epidemic=c(0.1,0.9,3), 
                                   endemic=c(-1,1,3)))
grid2 <- create.grid(disProgObj, cntrl2, 
                     params = list(epidemic=c(0.1,0.9,5), 
                                   negbin=c(0.3,12,10)))

}
\keyword{misc}
