% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{sowan_2_net}
\alias{sowan_2_net}
\title{Find matrices that fit the data}
\usage{
sowan_2_net(node_data, n_samples, nloptr_opts = list(algorithm =
  "NLOPT_LD_MMA", xtol_rel = 1e-04, maxeval = 1000))
}
\arguments{
\item{node_data}{is (coercible to) a data.frame with three columns:
\describe{
  \item{First column}{time period id}
  \item{Second column}{node id}
  \item{Third column}{share of node \eqn{i} in the total weight at time \eqn{t}. For each time period, this column should add up to 1}
}
Note that \code{node_data} will be sorted in ascending order by the first and second columns to enforce consistency. Thus, the solution matrices will be arranged accordingly.}

\item{n_samples}{number of solutions to search}

\item{nloptr_opts}{\code{opts} parameters (as list) for \code{nloptr}}
}
\value{
a \code{tibble} describing the solutions found:
\describe{
  \item{index}{index of the solution}
  \item{status}{Status code of \code{nloptr} (should equal 4 if everything went OK)}
  \item{error}{Value of the objective function at the solution}
  \item{sum_M}{Sum of the entries of the solution matrix (should equal 1 if OK)}
  \item{M}{list-column where each element is a solution matrix}
}
}
\description{
This function explores the solution space by solving the problem multiple times starting from random initial values
}
