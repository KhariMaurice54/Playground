\name{Likelydata}
\alias{Likelydata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Likelydata(data = temp, a1 = 1, a2 = 2, b1 = 1, b2 = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (data = temp, a1 = 1, a2 = 2, b1 = 1, b2 = 2) 
{
    phij1 = inverselogit(a1 + b1 * data$x1)
    phij2 = inverselogit(a2 + b2 * data$x2)
    phij1phij2 = cbind(phij1, phij2)
    minphij = apply(X = phij1phij2, MARGIN = 1, FUN = min)
    return(sum(data$y * log(minphij) + (1 - data$y) * log(1 - 
        minphij)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
