##############################################################################
# main function
##############################################################################

#' Find matrices that fit the data
#'
#' This function explores the solution space by solving the problem multiple times starting from random initial values
#'
#' @param node_data is (coercible to) a data.frame with three columns:
#' \describe{
#'   \item{First column}{time period id}
#'   \item{Second column}{node id}
#'   \item{Third column}{share of node \eqn{i} in the total weight at time \eqn{t}. For each time period, this column should add up to 1}
#' }
#' Note that \code{node_data} will be sorted in ascending order by the first and second columns to enforce consistency. Thus, the solution matrices will be arranged accordingly.
#' @param n_samples number of solutions to search
#' @param nloptr_opts \code{opts} parameters (as list) for \code{nloptr}
#' @return a \code{tibble} describing the solutions found:
#' \describe{
#'   \item{index}{index of the solution}
#'   \item{status}{Status code of \code{nloptr} (should equal 4 if everything went OK)}
#'   \item{error}{Value of the objective function at the solution}
#'   \item{sum_M}{Sum of the entries of the solution matrix (should equal 1 if OK)}
#'   \item{M}{list-column where each element is a solution matrix}
#' }
#' @export
sowan_2_net = function(node_data, n_samples, nloptr_opts = list(algorithm = "NLOPT_LD_MMA", xtol_rel = 1e-4, maxeval = 1000)){

  # arrange data
  node_data = as.data.frame(node_data[do.call(order, node_data), ])

  # setup parameters
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
