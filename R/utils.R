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
#' @export
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
#' @export
sym_mat_2_vec = function(M){
  return(M[lower.tri(M)])
}
