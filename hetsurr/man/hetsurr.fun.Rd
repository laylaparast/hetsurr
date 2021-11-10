\name{hetsurr.fun}
\alias{hetsurr.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assess and test for heterogeneity in the utility of a surrogate
}
\description{
This function assesses and tests for heterogeneity in the utility of a surrogate with respect to a continuous or discrete baseline covariate. The method is nonparametric and focuses on estimation of the proportion of the treatment effect that is explained by the surrogate marker. Inference information is provided including a confidence band for the proportion and results from formal testing for heterogeneity. When the baseline covariate is continuous, two tests are performed - an omnibus test and a trend test; when the baseline covariate is discrete, only one test is performed - an omnibus test. 
}
\usage{
hetsurr.fun(y1, y0, s1, s0, w1, w0, wf.grd = NULL, h0 = NULL, h1 = NULL, h2 = NULL, 
h3 = NULL, h4 = NULL, var.want = FALSE, type = "cont", test.want = FALSE, c.adj = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
the primary outcome in the treatment group
}
  \item{y0}{
the primary outcome in the control group
}
  \item{s1}{
the surrogate marker in the treatment group
}
  \item{s0}{
the surrogate marker in the control group
}
  \item{w1}{
the baseline covariate of interest in the treatment group
}
  \item{w0}{
the baseline covariate of interest in the control group
}
  \item{wf.grd}{
the grid of the baseline covariate where estimates are wanted; if this is not specified, this will bet set equal to 50 equally spaced points between the 10th and 90th percentile of the union of w1 and w0 values (if continuous) or the set of unique categories (if discrete)
}
  \item{h0}{
bandwidth; if not specific, this is set to a value of appropriate order
}
  \item{h1}{
bandwidth; if not specific, this is set to a value of appropriate order
}
  \item{h2}{
bandwidth; if not specific, this is set to a value of appropriate order
}
  \item{h3}{
bandwidth; if not specific, this is set to a value of appropriate order
}
  \item{h4}{
bandwidth; if not specific, this is set to a value of appropriate order
}
  \item{var.want}{
TRUE or FALSE; indicates whether variance estimate are wanted, default is FALSE
}
  \item{type}{
"cont" or "discrete"; indicates the form of the baseline covariate w, default is "cont"
}
  \item{test.want}{
TRUE or FALSE; indicates whether testing results are wanted, default is FALSE
}
  \item{c.adj}{
a multiplicative adjustment for the bandwidth; default is 1
}
}
\details{
This function calculates the treatment effect, Delta, the residual treatment effect, Delta_S, and the proportion of treatment effect explained R_S, as a function of the baseline covariate (if continuous) or stratified by category (if discrete). The estimation approach is nonparametric. Variance estimates can be requested including a confidence band for R_S. Results from formal testing for heterogeneity can also be requested. When the baseline covariate is continuous, two tests are performed - an omnibus test and a trend test; when the baseline covariate is discrete, only one test is performed - an omnibus test.
}
\value{
\item{w.values}{values of the baseline covariate at which estimates are provided} 
\item{delta.w}{estimated treatment effect, Delta, at each w value} 
\item{delta.w.s}{estimated residual treatment effect, Delta_S, at each w value} 
\item{R.w.s}{estimated proportion of treatment effect explained, R_S, at each w value} 
\item{se.delta.w}{estimated standard error of Delta} 
\item{se.delta.w.s}{estimated standard error of Delta.S} 
\item{se.R.w.s}{estimated standard error of R.S} 
\item{conf.delta.w.lower}{lower bound of the confidence intervals for Delta} 
\item{conf.delta.w.upper}{upper bound of the confidence intervals for Delta} 
\item{conf.delta.w.s.lower}{lower bound of the confidence intervals for Delta.S} 
\item{conf.delta.w.s.upper}{upper bound of the confidence intervals for Delta.S} 
\item{conf.R.w.s.lower}{lower bound of the confidence intervals for R.S} 
\item{conf.R.w.s.upper}{upper bound of the confidence intervals for R.S} 
\item{band.R.w.s.lower}{lower bound of the confidence band for R.S} 
\item{band.R.w.s.upper}{upper bound of the confidence band for R.S} 
\item{omnibus.test.statistic}{test statistic from omnibus test; for a continuous baseline covariate only } 
\item{omnibus.p.value}{p-value from omnibus test, testing for evidence of heterogeneity; for a continuous baseline covariate only} 
\item{trend.test.statistic}{test statistic from the trend test; for a continuous baseline covariate only} 
\item{trend.p.value}{p-value from trend test, testing for evidence of heterogeneity; for a continuous baseline covariate only} 
\item{test.statistic}{p-value from the test testing for evidence of heterogeneity; for a discrete baseline covariate only} 
\item{p.value}{p-value from the test testing for evidence of heterogeneity; for a discrete baseline covariate only} 
%% ...
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
	length=50), var.want =FALSE, type = "cont")
	het.ob = hetsurr.fun(y1=example.data$y1, y0=example.data$y0, s1=example.data$s1, 
	s0=example.data$s0, w1=example.data$w1, w0=example.data$w0, wf.grd = seq(0.25, 1.75, 
	length=50), var.want =TRUE,  test.want = TRUE, type = "cont")
	}
	#reducing dimension of example data to provide a computationally faster example
	het.ob = hetsurr.fun(y1=example.data$y1[1:500], y0=example.data$y0[1:500],
	s1=example.data$s1[1:500], s0=example.data$s0[1:500], w1=example.data$w1[1:500],
	w0=example.data$w0[1:500], wf.grd = seq(0.25, 1.75, length=5), var.want =FALSE, 
	type = "cont")

}

