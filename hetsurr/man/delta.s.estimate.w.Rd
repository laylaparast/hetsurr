\name{delta.s.estimate.w}
\alias{delta.s.estimate.w}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates residual treatment effect when W is discrete.
}
\description{
Estimates residual treatment effect when W is discrete.
}
\usage{
delta.s.estimate.w(sone, szero, yone, yzero, wone, wzero, w.all, weight.perturb = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sone}{
sone
}
  \item{szero}{
szero
}
  \item{yone}{
yone
}
  \item{yzero}{
yzero
}
  \item{wone}{
wone
}
  \item{wzero}{
wzero
}
  \item{w.all}{
w.all
}
  \item{weight.perturb}{
perturbation weights
}
}
\value{Returns a vector of length equal to the length of w.all with the estimated residual treatment effect for each w value}
\keyword{internal}