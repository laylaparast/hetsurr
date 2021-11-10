\name{hetsurr.plot}
\alias{hetsurr.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs plots of the treatment effect, residual treatment effect, and proportion explained as a function of the continuous baseline covariate
}
\description{
Constructs plots of the treatment effect, residual treatment effect, and proportion explained as a function of the continuous baseline covariate
}
\usage{
hetsurr.plot(hetsurr.results, xlab.name = "Baseline Covariate", placement = "bottomleft")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hetsurr.results}{
This is a list produced by the hetsurr.fun function. 
}
  \item{xlab.name}{
This is the label for the x-axis; default is "Baseline Covariate".
}
  \item{placement}{
This is the placement for the legend; default is "bottomleft".
}
}
\value{
Produces a plot with 3 panels. 
}
\references{
Parast, L., Cai, T., Tian L (2021). Testing for Heterogeneity in the Utility of a Surrogate Marker. Biometrics, In press.  
}
\author{
Layla Parast
}
\examples{
	data(example.data)
	names(example.data)
	#computationally intensive
	\donttest{
	het.ob = hetsurr.fun(y1=example.data$y1, y0=example.data$y0, s1=example.data$s1, 
	s0=example.data$s0, w1=example.data$w1, w0=example.data$w0, wf.grd = seq(0.25, 1.75, 
	length=50), var.want =TRUE,  test.want = TRUE, type = "cont")
	hetsurr.plot(het.ob)
	}
	#reducing dimension of example data to provide a computationally faster example
	het.ob = hetsurr.fun(y1=example.data$y1[1:300], y0=example.data$y0[1:300],
	s1=example.data$s1[1:300], s0=example.data$s0[1:300], w1=example.data$w1[1:300],
	w0=example.data$w0[1:300], wf.grd = seq(0.25, 1.75, length=5), var.want =TRUE,
	test.want = TRUE, type = "cont")
	hetsurr.plot(het.ob)

}

