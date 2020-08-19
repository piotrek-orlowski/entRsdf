entropy_functions <- R6::R6Class("entropy_functions",
  public = list(
    objective = function(theta_vector, return_matrix) {
      NA_real_
    },
    gradient = function(theta_vector, return_matrix) {
      matrix(NA_real_, nrow = length(theta_vector), ncol = 1L)
    },
    hessian = function(theta_vector, return_matrix) {
      matrix(NA_real_, nrow = length(theta_vector, ncol = length(theta_vector)))
    },
    sdf_recovery = function(theta_vector, return_matrix) {
      matrix(NA_real_, nrow = nrow(return_matrix), ncol = 1L)
    },
    get_description = function() {
      cat(private$description)
    }
  ),
  private = list(
    cr_power = numeric(1),
    sdf_average = 1.0,
    description = "This is a generic class for entropy SDF functions"
  )
)

et_functions <- R6::R6Class("et_functions",
  inherit = entropy_functions,
  public = list(
    initialize = function(sdf_mean = 1.0) {
      private$cr_power <- 0.0
      private$sdf_average <- sdf_mean
    },
    objective = function(theta_vector, return_matrix) {
      res <- exp(-return_matrix %*% theta_vector - 1.0)
      res <- mean(res) + theta_vector[1L]
      return(res)
    },
    gradient = function(theta_vector, return_matrix) {
      res <- exp(-return_matrix %*% theta_vector - 1.0)
      res <- -return_matrix * matrix(res, nrow = nrow(return_matrix), ncol = ncol(return_matrix))
      res <- apply(res, 2L, mean)
      return(res)
    },
    hessian = function(theta_vector, return_matrix) {
      res <- exp(-return_matrix %*% theta_vector - 1.0)
      res <- t(return_matrix) %*% (return_matrix * matrix(res, nrow = nrow(return_matrix), ncol = ncol(return_matrix)))
      res <- 1.0 / nrow(return_matrix) * res
      return(res)
    },
    sdf_recovery = function(theta_vector, return_matrix) {
      res <- exp(-return_matrix %*% theta_vector - 1.0)
      return(res)
    }
  ),
  private = list(
    description = "Functions for fitting and recovering exponential-tilting SDFs. \nThey correspond to setting power = 0 in the Cressie-Read function family."
  )
)

cressie_read_functions <- R6::R6Class("cressie_read_functions",
  inherit = entropy_functions,
  public = list(
    initialize = function(cressie_read_power,
                          sdf_mean) {
      if (cressie_read_power %in% c(-1, 0)) {
        stop("This function is for powers different from -1 and 0. If you wish to use those, please use et_functions or entropy_functions.")
      }
      private$cr_power <- cressie_read_power
      private$sdf_average <- sdf_mean
    },
    objective = function(theta_vector, return_matrix) {
      sdf_expect <- private$sdf_average
      cr_pow <- private$cr_power
      res <- cr_pow * return_matrix %*% theta_vector
      res <- res + sdf_expect^cr_pow
      res <- -1.0 / (1.0 + cr_pow) * res ^ ((cr_pow + 1.0) / cr_pow)
      res <- mean(res)
      return(-res)
    },
    gradient = function(theta_vector, return_matrix) {
      sdf_expect <- private$sdf_average
      cr_pow <- private$cr_power
      res <- cr_pow * return_matrix %*% theta_vector
      res <- res + sdf_expect^cr_pow
      res <- -res^(1.0 / cr_pow)
      res <- apply(return_matrix, 2, function(x) res * x)
      res <- apply(res, 2, mean)
      return(-res)
    },
    hessian = function(theta_vector, return_matrix) {
      sdf_expect <- private$sdf_average
      cr_pow <- private$cr_power
      # return_matrix <- return_matrix - 1.0 / sdf_expect # because excess returns
      res <- cr_pow * return_matrix %*% theta_vector
      res <- res + sdf_expect^cr_pow
      res <- -res^((1.0 - cr_pow) / cr_pow)
      res <- apply(cbind(res, return_matrix), 1, function(x) x[1] * (x[-1] %*% t(x[-1])))
      res <- array(res, c(length(theta_vector), length(theta_vector), nrow(return_matrix)))
      res <- apply(res, c(1, 2), mean)
      return(-res)
    },
    sdf_recovery = function(theta_vector, return_matrix) {
      sdf_expect <- private$sdf_average
      cr_pow <- private$cr_power
      # return_matrix <- return_matrix - 1.0 / sdf_expect # because excess returns
      res <- cr_pow * return_matrix %*% theta_vector
      res <- res + sdf_expect^cr_pow
      res <- res^(1.0 / cr_pow)
      res <- sdf_expect * res / mean(res)
      return(res)
    },
    get_power = function() {
      return(private$cr_power)
    },
    get_risk_free_return = function() {
      return(1.0 / private$sdf_average)
    },
    get_sdf_mean = function() {
      return(private$sdf_average)
    },
    set_power = function(new_power) {
      private$cr_power <- new_power
      invisible(self)
    },
    set_risk_free_return = function(new_risk_free_return) {
      private$sdf_average <- 1.0 / new_risk_free_return
      invisible(self)
    },
    set_sdf_mean = function(new_sdf_mean) {
      private$sdf_average <- new_sdf_mean
      invisible(self)
    }
  ),
  private = list(
    description = "Functions for fitting and recovering general Cressie-Read SDFs. \nThey correspond to any power different from -1 or 0 in the Cressie-Read function family."
  )
)

distance_et_functions <- R6::R6Class("distance_et_functions",
  inherit = entropy_functions,
  public = list(
    objective = function(theta_vector, return_matrix, penalty_value) {

      # Check if there is a pre-estimated lambda (to get a starting value for the first optimisation)
      if (is.null(private$lambda_opt)) {
        lambda_opt <- rep(0, ncol(return_matrix))
      } else {
        lambda_opt <- private$lambda_opt
      }

      # Optimize the inner problem given theta
      # These inner objective / gradient / Hessian are coded in entropic-distance-estimators.cpp
      def_opts <- nloptr::nl.opts()
      def_opts$algorithm <- "NLOPT_LD_LBFGS"
      optimisation_result <- nloptr::nloptr(
        x0 = lambda_opt,
        eval_f = et_distance_objective,
        opts = def_opts,
        theta_extended = theta_vector,
        return_matrix = return_matrix,
        mu_penalty = penalty_value
      )

      # Recover optimal decision
      lambda_opt <- optimisation_result$solution

      private$lambda_opt <- lambda_opt

      # Calculate gradient wrt theta
      gradient_theta <- et_distance_theta_gradient(
        lambda_opt,
        theta_vector,
        return_matrix,
        penalty_value
      )

      # Return objective value and gradient
      # obj value has to switch sign, because we minimize dist for max dual,
      # but the dual here was minimised (that's how solvers work)
      # and we have to switch signs again
      return(list(
        objective = -optimisation_result$objective,
        gradient = gradient_theta
      ))
    },
    get_lambda_stored = function() {
      private$lambda_opt
    },
    sdf_recovery = function(theta_vector, return_matrix) {
      res <- exp(return_matrix %*% theta_vector)
      return(res)
    }
  ),
  private = list(
    lambda_opt = NULL,
    description = "Functions for fitting and recovering exponential-tilting SDFs. \nThey correspond to setting power = 0 in the Cressie-Read function family."
  )
)

bivariate_mgf_functions <- R6::R6Class("bivariate_mgf_functions",
  public = list(
    initialize = function(sdf_powers) {
      private$sdf_powers <- sdf_powers

      mgf_powers <- sdf_powers / (sum(sdf_powers) - 1)
      private$mgf_powers <- mgf_powers
    },
    objective = function(theta_vector, home_return_matrix, foreign_return_matrix) {
      num_home <- ncol(home_return_matrix)
      num_foreign <- ncol(foreign_return_matrix)
      num_tot <- num_home + num_foreign

      return_matrix <- cbind(home_return_matrix, foreign_return_matrix)

      theta_matrix <- matrix(0.0, nrow = num_home + num_foreign, ncol = 2L)

      # first column of each return matrix is Rf
      theta_matrix[1L, 1] <- 1 - sum(theta_vector[1L:(num_home - 1L)])
      theta_matrix[2L:num_home, 1L] <- theta_vector[1L:(num_home - 1L)]

      theta_matrix[num_home + 1L, 2L] <- 1 - sum(theta_vector[(num_home - 1L + 1L):(num_tot - 2L)])
      theta_matrix[(num_home + 2L):num_tot, 2L] <- theta_vector[(num_home - 1L + 1L):(num_tot - 2L)]


      portfolio_return_matrix <- return_matrix %*% theta_matrix

      if (any(portfolio_return_matrix < 0)) {
        return(1e6)
      }
      mgf_powers <- private$mgf_powers

      portfolio_return_matrix[, 1L] <- portfolio_return_matrix[, 1L]^mgf_powers[1L]
      portfolio_return_matrix[, 2L] <- portfolio_return_matrix[, 2L]^mgf_powers[2L]

      portfolio_return_matrix <- apply(
        portfolio_return_matrix,
        1L,
        prod
      )

      res <- mean(portfolio_return_matrix)

      return(res)
    },
    gradient = function(theta_vector, home_return_matrix, foreign_return_matrix) {
      num_home <- ncol(home_return_matrix)
      num_foreign <- ncol(foreign_return_matrix)
      num_tot <- num_home + num_foreign

      return_matrix <- cbind(home_return_matrix, foreign_return_matrix)

      theta_matrix <- matrix(0.0, nrow = num_home + num_foreign, ncol = 2L)

      # first column of each return matrix is Rf
      theta_matrix[1L, 1] <- 1 - sum(theta_vector[1L:(num_home - 1L)])
      theta_matrix[2L:num_home, 1L] <- theta_vector[1L:(num_home - 1L)]

      theta_matrix[num_home + 1L, 2L] <- 1 - sum(theta_vector[(num_home - 1L + 1L):(num_tot - 2L)])
      theta_matrix[(num_home + 2L):num_tot, 2L] <- theta_vector[(num_home - 1L + 1L):(num_tot - 2L)]


      portfolio_return_matrix <- return_matrix %*% theta_matrix

      if (any(portfolio_return_matrix < 0)) {
        return(1e6)
      }
      mgf_powers <- private$mgf_powers

      portfolio_return_powers_min_one <- portfolio_return_powers <- portfolio_return_matrix

      portfolio_return_powers[, 1L] <- portfolio_return_matrix[, 1L]^mgf_powers[1L]
      portfolio_return_powers[, 2L] <- portfolio_return_matrix[, 2L]^mgf_powers[2L]

      portfolio_return_powers_min_one[, 1L] <- portfolio_return_matrix[, 1L]^(mgf_powers[1L] - 1.0)
      portfolio_return_powers_min_one[, 2L] <- portfolio_return_matrix[, 2L]^(mgf_powers[2L] - 1.0)

      home_excess_return_matrix <- apply(
        home_return_matrix[, -1L, drop = FALSE],
        2L,
        function(ret) ret - home_return_matrix[, 1L]
      )

      foreign_excess_return_matrix <- apply(
        foreign_return_matrix[, -1L, drop = FALSE],
        2L,
        function(ret) ret - foreign_return_matrix[, 1L]
      )

      home_excess_return_matrix <- apply(
        home_excess_return_matrix,
        2L,
        function(ret) ret * portfolio_return_powers_min_one[, 1L] * portfolio_return_powers[, 2L] * mgf_powers[1L]
      )

      foreign_excess_return_matrix <- apply(
        foreign_excess_return_matrix,
        2L,
        function(ret) ret * portfolio_return_powers_min_one[, 2L] * portfolio_return_powers[, 1L] * mgf_powers[2L]
      )

      excess_return_matrix <- cbind(home_excess_return_matrix, foreign_excess_return_matrix)

      # excess_return_matrix <- apply(excess_return_matrix
      #                               , 2L
      #                               , function(ret) ret * portfolio_return_powers)

      gradient_res <- apply(excess_return_matrix, 2L, mean)

      return(gradient_res)
    },
    hessian = function(theta_vector, home_return_matrix, foreign_return_matrix) {
      num_home <- ncol(home_return_matrix)
      num_foreign <- ncol(foreign_return_matrix)
      num_tot <- num_home + num_foreign

      return_matrix <- cbind(home_return_matrix, foreign_return_matrix)

      theta_matrix <- matrix(0.0, nrow = num_home + num_foreign, ncol = 2L)

      # first column of each return matrix is Rf
      theta_matrix[1L, 1] <- 1 - sum(theta_vector[1L:(num_home - 1L)])
      theta_matrix[2L:num_home, 1L] <- theta_vector[1L:(num_home - 1L)]

      theta_matrix[num_home + 1L, 2L] <- 1 - sum(theta_vector[(num_home - 1L + 1L):(num_tot - 2L)])
      theta_matrix[(num_home + 2L):num_tot, 2L] <- theta_vector[(num_home - 1L + 1L):(num_tot - 2L)]


      portfolio_return_matrix <- return_matrix %*% theta_matrix

      if (any(portfolio_return_matrix < 0)) {
        return(1e6)
      }
      mgf_powers <- private$mgf_powers

      portfolio_return_powers_min_two <- portfolio_return_powers_min_one <- portfolio_return_powers <- portfolio_return_matrix

      portfolio_return_powers[, 1L] <- portfolio_return_matrix[, 1L] ^ mgf_powers[1L]
      portfolio_return_powers[, 2L] <- portfolio_return_matrix[, 2L] ^ mgf_powers[2L]

      portfolio_return_powers_min_one[, 1L] <- portfolio_return_matrix[, 1L] ^ (mgf_powers[1L] - 1.0)
      portfolio_return_powers_min_one[, 2L] <- portfolio_return_matrix[, 2L] ^ (mgf_powers[2L] - 1.0)

      portfolio_return_powers_min_two[, 1L] <- portfolio_return_matrix[, 1L] ^ (mgf_powers[1L] - 2.0)
      portfolio_return_powers_min_two[, 2L] <- portfolio_return_matrix[, 2L] ^ (mgf_powers[2L] - 2.0)

      home_excess_return_matrix <- apply(
        home_return_matrix[, -1L, drop = FALSE],
        2L,
        function(ret) ret - home_return_matrix[, 1L]
      )

      foreign_excess_return_matrix <- apply(
        foreign_return_matrix[, -1L, drop = FALSE],
        2L,
        function(ret) ret - foreign_return_matrix[, 1L]
      )
      home_size <- dim(home_excess_return_matrix)
      foreign_size <- dim(foreign_excess_return_matrix)

      home_excess_return_matrix <- cbind(home_excess_return_matrix, matrix(0, foreign_size[1L], foreign_size[2L]))

      foreign_excess_return_matrix <- cbind(matrix(0, home_size[1L], home_size[2L]), foreign_excess_return_matrix)

      hessian_home_foreign <- sapply(1L:nrow(home_return_matrix),
        function(ind) {
          res <- outer(
            home_excess_return_matrix[ind, ],
            foreign_excess_return_matrix[ind, ]
          )
          res <- res * prod(mgf_powers)
          res <- res * prod(portfolio_return_powers_min_one[ind, ])
          res
        },
        simplify = FALSE
      )

      hessian_home_home <- sapply(1L:nrow(home_return_matrix),
        function(ind) {
          res <- outer(
            home_excess_return_matrix[ind, ],
            home_excess_return_matrix[ind, ]
          )
          res <- res * mgf_powers[1L] * (mgf_powers[1L] - 1.0)
          res <- res * portfolio_return_powers[ind, 2L] * portfolio_return_powers_min_two[ind, 1L]
          res
        },
        simplify = FALSE
      )

      hessian_foreign_foreign <- sapply(1L:nrow(home_return_matrix),
        function(ind) {
          res <- outer(
            foreign_excess_return_matrix[ind, ],
            foreign_excess_return_matrix[ind, ]
          )
          res <- res * mgf_powers[2L] * (mgf_powers[2L] - 1.0)
          res <- res * portfolio_return_powers[ind, 1L] * portfolio_return_powers_min_two[ind, 2L]
          res
        },
        simplify = FALSE
      )

      hessian_all <- sapply(1L:nrow(home_return_matrix), function(ind) {
        res <- hessian_home_home[[ind]] + hessian_foreign_foreign[[ind]] + hessian_home_foreign[[ind]] + t(hessian_home_foreign[[ind]])
        res
      },
      simplify = FALSE
      )

      hessian_all <- abind::abind(hessian_all, along = 3L)

      hessian_all <- apply(hessian_all, c(1L, 2L), mean)

      return(hessian_all)
    },
    sdf_recovery = function(theta_vector, home_return_matrix, foreign_return_matrix) {
      num_home <- ncol(home_return_matrix)
      num_foreign <- ncol(foreign_return_matrix)
      num_tot <- num_home + num_foreign

      return_matrix <- cbind(home_return_matrix, foreign_return_matrix)

      theta_matrix <- matrix(0.0, nrow = num_home + num_foreign, ncol = 2L)

      # first column of each return matrix is Rf
      theta_matrix[1L, 1] <- 1 - sum(theta_vector[1L:(num_home - 1L)])
      theta_matrix[2L:num_home, 1L] <- theta_vector[1L:(num_home - 1L)]

      theta_matrix[num_home + 1L, 2L] <- 1 - sum(theta_vector[(num_home - 1L + 1L):(num_tot - 2L)])
      theta_matrix[(num_home + 2L):num_tot, 2L] <- theta_vector[(num_home - 1L + 1L):(num_tot - 2L)]


      portfolio_return_matrix <- return_matrix %*% theta_matrix

      # SDF denominator
      portfolio_return_mgf_powers <- sapply(
        seq_along(private$mgf_powers),
        function(ind) {
          portfolio_return_matrix[, ind] ^ private$mgf_powers[ind]
        }
      )
      portfolio_return_mgf_powers <- apply(portfolio_return_mgf_powers, 1L, prod)
      sdf_denominator <- mean(portfolio_return_mgf_powers)

      # SDF numerators
      home_sdf <- portfolio_return_matrix[, 1L]^(1.0 - private$sdf_powers[2L])
      home_sdf <- home_sdf * portfolio_return_matrix[, 2L]^private$sdf_powers[2L]
      home_sdf <- home_sdf ^ (1.0 / (sum(private$sdf_powers) - 1.0))

      foreign_sdf <- portfolio_return_matrix[, 2L] ^ (1.0 - private$sdf_powers[1L])
      foreign_sdf <- foreign_sdf * portfolio_return_matrix[, 1L] ^ private$sdf_powers[1L]
      foreign_sdf <- foreign_sdf ^ (1.0 / (sum(private$sdf_powers) - 1.0))

      home_sdf <- home_sdf / sdf_denominator
      foreign_sdf <- foreign_sdf / sdf_denominator

      return(cbind(home_sdf, foreign_sdf))
    },
    set_powers = function(new_powers) {
      private$mgf_powers <- new_powers
      invisible(self)
    },
    get_powers = function() {
      private$mgf_powers
    }
  ),
  private = list(
    mgf_powers = NULL,
    sdf_powers = NULL,
    lambda_opt = NULL,
    description = "AFD bivariate functions"
  )
)