\name{cleversearch}
\alias{cleversearch}
\concept{minimization}
\concept{maximization}
\title{Optimization over a parameter grid}
\description{
  This function allows greedy/full grid search in any dimension.
}
\usage{
cleversearch(fn, lower, upper, ngrid, startvalue, logscale = TRUE,
             clever = TRUE, verbose = FALSE) 
}
\arguments{
  \item{fn}{a function to be minimized (or maximized), with the only
    argument being the parameter over which minimization is to take
    place. The return value must be scalar.}   
  \item{lower}{numeric vector containing the lower bounds on the
    parameter grid} 
  \item{upper}{numeric vector containing the upper bounds on the
    parameter grid}
  \item{ngrid}{integer number determining the grid length in every dimension}
  \item{startvalue}{optional initial value for the parameter to be
    optimized over} 
  \item{logscale}{logical, whether to construct the grid on logarithmic
    scale}   
  \item{clever}{logical, whether to perform a greedy grid search with
    lookup-table or a full grid evaluation. The latter is only available
    up to 3d.}  
  \item{verbose}{logical. Should the search process be monitored?}
}
\details{
  Unless \code{startvalue} is specified, the search starts at the lower
  bound of the 1d parameter space or at the middle of the 2d/3d grid.
}
\value{
  A list with components
  \item{par}{optimal parameter value that was found.}
  \item{value}{\code{fn} value corresponding to \code{par}.}
  \item{counts}{number of calls to 'fn'.}
}
\seealso{\code{\link{optim}}}
\examples{
simplefun <- function(vec) { return(sum(vec^2)) }
opt <- cleversearch(simplefun, c(-1, -1), c(1, 1), 51, logscale=FALSE)
}
\keyword{optimize}
\keyword{spatial}
