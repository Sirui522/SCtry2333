% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deri_plot.R
\name{deri.plot}
\alias{deri.plot}
\title{Plot of directional derivative}
\usage{
deri.plot(data_x, G, inf, sup)
}
\arguments{
\item{data_x}{sample from GMM}

\item{G}{the estimate of clustering}

\item{inf}{inf bound of the range for searching for parameters}

\item{sup}{sup bound of the range for searching for parameters}
}
\value{
maximum value of gradient function in some range.
}
\description{
Plot the directional derivative from the clustering estimator to dirac measures.
}
\examples{
\dontrun{
data_x = c(rnorm(500,0),rnorm(500,2));
M = CN.R(data_x, it_num = 20, inf = -1.5, sup = 3.5);
M_deri = deri.plot(M);M_deri
}
}
