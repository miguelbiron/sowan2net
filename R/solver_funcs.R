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
