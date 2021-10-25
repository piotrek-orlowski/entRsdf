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

solve_entropy_problem <- function(entropy_foos,
                                  excess_return_matrix,
                                  theta_vector_init = rep(1.0, ncol(excess_returns)) / ncol(excess_returns),
                                  solver_trace = FALSE,
                                  ...) {
  if (inherits(entropy_foos, "cressie_read_functions")) {
    sdf_mean <- entropy_foos$get_sdf_mean()
    entropy_power <- entropy_foos$get_power()
    if (entropy_foos$get_power() < 0) {
      pos_ret_constraints <- list(cccp::nnoc(
        G = excess_return_matrix,
        h = matrix(-sdf_mean^entropy_power / entropy_power, nrow = nrow(excess_return_matrix), ncol = 1)
      ))
    } else {
      pos_ret_constraints <- list(cccp::nnoc(
        G = -excess_return_matrix,
        h = matrix(sdf_mean^entropy_power / entropy_power, nrow = nrow(excess_return_matrix), ncol = 1)
      ))
    }
  } else {
    pos_ret_constraints <- list()
  }

  optimisation_result <- cccp::cccp(
    f0 = function(x) entropy_foos$objective(x, excess_return_matrix),
    g0 = function(x) entropy_foos$gradient(x, excess_return_matrix),
    h0 = function(x) entropy_foos$hessian(x, excess_return_matrix),
    x0 = theta_vector_init,
    optctrl = cccp::ctrl(trace = solver_trace, ...),
    cList = pos_ret_constraints
  )

  if (optimisation_result$status == "optimal") {
    return(cccp::getx(optimisation_result))
  } else {
    warning("The optimisation did not converge, make sure the results are ok")

    return(cccp::getx(optimisation_result))
  }
}

solve_bivariate_entropy_problem <- function(entropy_foos,
                                            home_return_matrix,
                                            foreign_return_matrix,
                                            theta_vector_init = rep(0.0, (ncol(home_return_matrix) + ncol(foreign_return_matrix)) - 2L),
                                            solver_trace = FALSE,
                                            ...) {
  home_excess_return_matrix <- apply(
    home_return_matrix[, -1L, drop = FALSE],
    2L,
    function(ret) {
      ret - home_return_matrix[, 1L]
    }
  )

  foreign_excess_return_matrix <- apply(
    foreign_return_matrix[, -1L, drop = FALSE],
    2L,
    function(ret) {
      ret - foreign_return_matrix[, 1L]
    }
  )

  home_size <- dim(home_excess_return_matrix)
  foreign_size <- dim(foreign_excess_return_matrix)

  home_excess_return_matrix <- cbind(home_excess_return_matrix, matrix(0.0, foreign_size[1L], foreign_size[2L]))
  foreign_excess_return_matrix <- cbind(matrix(0.0, home_size[1L], home_size[2L]), foreign_excess_return_matrix)

  pos_ret_constraints <- list(
    cccp::nnoc(
      G = -home_excess_return_matrix,
      h = matrix(home_return_matrix[, 1L], home_size[1L], 1L)
    ),
    cccp::nnoc(
      G = -foreign_excess_return_matrix,
      h = matrix(foreign_return_matrix[, 1L], foreign_size[1L], 1L)
    )
  )

  optimisation_result <- cccp::cccp(
    f0 = function(x) entropy_foos$objective(x, home_return_matrix, foreign_return_matrix),
    g0 = function(x) entropy_foos$gradient(x, home_return_matrix, foreign_return_matrix),
    h0 = function(x) entropy_foos$hessian(x, home_return_matrix, foreign_return_matrix),
    x0 = theta_vector_init,
    optctrl = cccp::ctrl(
      maxiters = 1e2L,
      trace = solver_trace,
      abstol = 1e-12,
      reltol = 1e-10,
      ...
    ),
    cList = pos_ret_constraints
  )

  if (optimisation_result$status == "optimal") {
    return(optimisation_result)
  } else {
    warning("The optimisation did not converge, make sure the results are ok")

    return(optimisation_result)
  }
}

solve_bivariate_entropy_problem_nloptr <- function(entropy_foos,
                                            home_return_matrix,
                                            foreign_return_matrix,
                                            theta_vector_init = rep(0.0, (ncol(home_return_matrix) + ncol(foreign_return_matrix)) - 2L),
                                            solver_trace = FALSE,
                                            ...) {
  home_excess_return_matrix <- apply(
    home_return_matrix[, -1L, drop = FALSE],
    2L,
    function(ret) {
      ret - home_return_matrix[, 1L]
    }
  )
  
  foreign_excess_return_matrix <- apply(
    foreign_return_matrix[, -1L, drop = FALSE],
    2L,
    function(ret) {
      ret - foreign_return_matrix[, 1L]
    }
  )
  
  home_size <- dim(home_excess_return_matrix)
  foreign_size <- dim(foreign_excess_return_matrix)
  
  home_excess_return_matrix <- cbind(home_excess_return_matrix, matrix(0.0, foreign_size[1L], foreign_size[2L]))
  foreign_excess_return_matrix <- cbind(matrix(0.0, home_size[1L], home_size[2L]), foreign_excess_return_matrix)
  
  pos_ret_constraint <- function(theta) {
    rbind(
      home_excess_return_matrix,
      foreign_excess_return_matrix
    ) %*% 
      theta + 
      rbind(
        home_return_matrix[,1L,drop=F],
        foreign_return_matrix[,1L,drop=F]
      )
  }
  
  pos_ret_constraint_jac <- function(theta) {
    rbind(home_excess_return_matrix, foreign_excess_return_matrix)
  }
  
  optimisation_result <- nloptr::slsqp(
    fn = function(x) entropy_foos$objective(x, home_return_matrix, foreign_return_matrix),
    gr = function(x) entropy_foos$gradient(x, home_return_matrix, foreign_return_matrix),
    x0 = theta_vector_init,
    hin = pos_ret_constraint, 
    hinjac = pos_ret_constraint_jac
  )
  
  if (optimisation_result$convergence != -1) {
    return(optimisation_result)
  } else {
    warning("The optimisation did not converge, make sure the results are ok")
    
    return(optimisation_result)
  }
}