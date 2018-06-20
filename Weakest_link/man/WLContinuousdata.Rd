\name{WLContinuousdata}
\alias{WLContinuousdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulated data set w/ Gaussian predictors for testing weakest Link model analysis methods}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WLContinuousdata(n = 200, p0 = 0.1, p1 = 0.5, a1 = 1, a2 = 2, b1 = 1, b2 = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p0}{
%%     ~~Describe \code{p0} here~~
}
  \item{p1}{
%%     ~~Describe \code{p1} here~~
}
  \item{a1}{
%%     ~~Describe \code{a1} here~~
}
  \item{a2}{
%%     ~~Describe \code{a2} here~~
}
  \item{b1}{
%%     ~~Describe \code{b1} here~~
}
  \item{b2}{
%%     ~~Describe \code{b2} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n = 200, p0 = 0.1, p1 = 0.5, a1 = 1, a2 = 2, b1 = 1, 
    b2 = 2) 
{
    x1 = rnorm(n = n)
    x2 = rnorm(n = n)
    x1x2 = cbind(x1, x2)
    phi1 = a1 + b1 * x1
    phi2 = a2 + b2 * x2
    phi1phi2 = cbind(phi1, phi2)
    minphi = apply(X = phi1phi2, MARGIN = 1, FUN = min)
    yprob = inverselogit(minphi)
    y = rbinom(n = n, size = 1, prob = yprob)
    return(data.frame(x1, x2, y, minphi))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
