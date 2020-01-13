# Thu Oct 17 15:54:44 2019 ------------------------------
#' Generic solver for minimum-entropy pricing kernels
#'
#' @description Given an objective function, its gradient and hessian (in the \code{entropy_foos} object), this function submits them to a convex solver and returns the vector of optimal portfolio weights plus a normalization weight in the first position. This function is supposed to be used internally in the implementation of generic methods for SDF fitting.
#' @param entropy_foos object of \code{entropy_functions} S4 class
#' @param excess_return_matrix T x N \code{matrix} of excess returns
#' @param theta_vector_init (N+1) x 1 \code{numeric} or \cide{matrix} of initial weights
#' @param solver_trace sets the \code{trace} argument in \code{cccp::ctrl} call which is passed to \code{cccp::cccp}
#' @param ... arguments passed to \link{\code{cccp::ctrl}}
#'
#' @details For N assets, the vector of portfolio weights is of length N+1, because in the first position it contains the Lagrange multiplier for the constraint that the average of the pricing kernel time series is equal to 1.0. Thus, the portfolio weights are in fact in positions \code{2:(N+1)}.
#'
#' @return (N+1) x 1 \code{matrix}
#'

solve_entropy_problem <- function(entropy_foos
                                  , excess_return_matrix
                                  , theta_vector_init = rep(1.0, ncol(excess_returns))/ncol(excess_returns)
                                  , solver_trace = FALSE
                                  , ...){
  
  if(inherits(entropy_foos, "cressie_read_functions")){
    pos_ret_constraints <- cccp::nnoc(G = -excess_return_matrix
                                      , h = matrix(0, nrow = nrow(excess_return_matrix), ncol = 1))
  } else {
    pos_ret_constraints <- list()
  }
  
  optimisation_result <- cccp::cccp(f0 = function(x) entropy_foos$objective(x, excess_return_matrix)
                                    , g0 = function(x) entropy_foos$gradient(x, excess_return_matrix)
                                    , h0 = function(x) entropy_foos$hessian(x, excess_return_matrix)
                                    , x0 = theta_vector_init
                                    , optctrl = cccp::ctrl(trace = solver_trace, ...)
                                    , cList = pos_ret_constraints)
  
  if(optimisation_result$status == "optimal"){
    
    return(cccp::getx(optimisation_result))
    
  } else {
    
    warning("The optimisation did not converge, make sure the results are ok")
    
    return(cccp::getx(optimisation_result))
  }
}