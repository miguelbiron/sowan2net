#######################################
# utils
#######################################

#' Funny matrices
#'
#' Building blocks for CS matrices
#'
#' @param k number of columns
#' @return a \eqn{(k+1) x k} matrix
get_funny_matrix = function(k){
  return(rbind(rep(1L, k), diag(k)))
}

#' Build a CS matrix
#'
#' @param n number of nodes
#' @return an \eqn{n x (n*(n-1)/2)} matrix
get_CS_matrix = function(n){
  stopifnot(n>1)
  m = as.integer(n*(n-1)/2)
  CS = matrix(0, n, m)
  j = 1
  for(i in 1:(n-1)){
    k = n-i
    CS[i:n, j:(j+k-1)] = get_funny_matrix(k)
    j = j + k
  }
  return(CS)
}

#' Vector to symmetric matrix
#'
#' @param V vector
#' @return an \eqn{n x n} symmetric matrix with diagonal \code{rep(0, n)}
vec_2_sym_mat = function(v){
  m = length(v)
  n = as.integer(0.5*(sqrt(1+8*m) + 1))
  M = matrix(0, n, n)
  M[lower.tri(M)] = v
  M[upper.tri(M)] = t(M)[upper.tri(M)]
  return(M)
}

#' Symmetric matrix to vector
#'
#' Returns the lower triangle (without the diagonal) of a symmetric matrix as a vector
#'
#' @param M symmetric matrix
#' @return a vector
sym_mat_2_vec = function(M){
  return(M[lower.tri(M)])
}

#######################################
# solver functions
#######################################

#' Objective function of the optimization problem
#'
#' @param x a vector of size \eqn{m} with current solution
#' @param C_bar a matrix of size \eqn{nT x m}
#' @param S_vec a vector of size \eqn{nT}
#' @return the value of the objective function at the current solution
f_ob = function(x, C_bar, S_vec){
  return(0.5*sum((as.vector(C_bar %*% x) - S_vec)^2))
}

#' Gradient of the objective function
#'
#' @param x a vector of size \eqn{m} with current solution
#' @param C_bar a matrix of size \eqn{nT x m}
#' @param S_vec a vector of size \eqn{nT}
#' @return the gradient of the objective function at the current solution
grad_f = function(x, C_bar, S_vec){
  return(as.vector(crossprod(C_bar, C_bar %*% x - S_vec)))
}

##############################################################################
# main function
##############################################################################

#' Find matrices that fit the data
#'
#' This function explores the solution space by solving the problem multiple times starting from random initial values
#'
#' @param node_data is (coercible to) a data.frame with three columns:
#' \describe{
#'   \item{First column}{time stamps}
#'   \item{Second column}{node ids}
#'   \item{Third column}{share of node i in the total weight at time t. For each time period, this column should add up to 1}
#' }
#' Note that \code{node_data} will be sorted in ascending order by the first and second columns to enforce consistency. Thus, the solution matrices will be arranged accordingly.
#' @param n_samples number of solutions to search
#' @param xtol_rel relative tolerance parameter for \code{nloptr}
#' @return a tibble describing the solutions found:
#' \describe{
#'   \item{index}{index of the solution}
#'   \item{status}{Status code of \code{nloptr} (should equal 4 if everything went OK)}
#'   \item{error}{Value of the objective function at the solution}
#'   \item{sum_M}{Sum of the entries of the solution matrix (should equal 1 if OK)}
#'   \item{M}{list-column where each row is a solution matrix}
#' }
#' @export
sowan_2_net = function(node_data, n_samples, nloptr_opts = list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-4, maxeval = 1000)){

  # node_data is (coercible as) a data.frame which contains
  #     - First column : time stamps
  #     - Second column: node ids
  #     - Third column : share of node i in the total weight at time t
  #         - For each time period, this column should add up to 1

  # arrange data
  node_data = as.data.frame(node_data[do.call(order, node_data), ])

  # Setup for solver
  n_nodes = length(unique(node_data[, 2])) # n
  n_params = as.integer(n_nodes * (n_nodes - 1) / 2) # m
  n_per = length(unique(node_data[, 1])) # T
  C_bar = do.call("rbind", rep(list(get_CS_matrix(n_nodes)), n_per))
  S_vec = node_data[, 3]
  rm(node_data)

  # Solve with nloptr: search solution space with random starting points
  results = lapply(1:n_samples, function(i){
    fit_l_nloptr = nloptr::nloptr(x0          = 0.5*runif(n_params),
                                  eval_f      = f_ob,
                                  eval_grad_f = grad_f,
                                  lb          = rep(0, n_params),
                                  ub          = rep(0.5, n_params),
                                  opts        = nloptr_opts,
                                  C_bar = C_bar, S_vec = S_vec)

    tibble::tibble(
      index  = i,
      status = fit_l_nloptr$status,
      error  = fit_l_nloptr$objective,
      sum_M  = 2*sum(fit_l_nloptr$solution),
      M      = list(vec_2_sym_mat(fit_l_nloptr$solution))
    )

  })

  return(do.call("rbind", results))

}
