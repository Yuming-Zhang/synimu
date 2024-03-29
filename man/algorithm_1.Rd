% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VaccaroZaki.R
\name{algorithm_1}
\alias{algorithm_1}
\title{Estimation of White Noise and Random Walk Innovation Variances \eqn{R} and \eqn{Q}}
\usage{
algorithm_1(X, remove_last_scale = FALSE)
}
\arguments{
\item{X}{A \code{matrix} of dimension T by p, where T is the length of the time series and p is the number of processes.}

\item{remove_last_scale}{wether the last scale of the Allan variance should be removed}
}
\value{
A \code{list} with the following structure:
\itemize{
 \item R: A \code{vector} of size p of the estimated white noise variances.
 \item Q: A \code{vector} of size p of the estimated random walk innovation variances.
}
}
\description{
This function provides an estimate based on the Allan variance for the white noise and the random walk innovation variances \eqn{R} and \eqn{Q}.
The implementation reflects Algorithm 1 in Richard J. Vaccaro and Ahmed S. Zaki, "Reduced-Drift Virtual Gyro from an Array of Low-Cost Gyros".
}
\author{
Davide Antonio Cucci
}
