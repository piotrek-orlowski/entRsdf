entropy_functions <- R6::R6Class("entropy_functions"
                                 , public = list(
                                   objective = function(theta_vector, return_matrix){
                                     NA_real_
                                   }
                                   , gradient = function(theta_vector, return_matrix){
                                     matrix(NA_real_, nrow = length(theta_vector), ncol = 1L)
                                   }
                                   , hessian = function(theta_vector, return_matrix){
                                     matrix(NA_real_, nrow = length(theta_vector, ncol = length(theta_vector)))
                                   }
                                   , sdf_recovery = function(theta_vector, return_matrix){
                                     matrix(NA_real_, nrow = nrow(return_matrix), ncol = 1L)
                                   }
                                   , get_description = function(){
                                     cat(private$description)
                                   }
                                 )
                                 , private = list(
                                   description = "This is a generic class for entropy SDF functions"
                                 ))

et_functions <- R6::R6Class("et_functions"
                               , inherit = entropy_functions
                               , public = list(
                                 objective = function(theta_vector, return_matrix){
                                    return_matrix <- cbind(1.0, return_matrix)
                                    res <- exp(- return_matrix %*% theta_vector - 1.0)
                                    res <- mean(res) + theta_vector[1L]
                                    return(res)
                                  }
                                 , gradient = function(theta_vector, return_matrix){
                                    return_matrix <- cbind(1.0, return_matrix)
                                    res <- exp(- return_matrix %*% theta_vector - 1.0)
                                    res <- - return_matrix * matrix(res, nrow = nrow(return_matrix), ncol = ncol(return_matrix))
                                    res <- apply(res, 2L, mean)
                                    res <- res + c(1.0, rep(0, ncol(return_matrix)-1.0))
                                    return(res)
                                  }
                                 , hessian = function(theta_vector, return_matrix){
                                    return_matrix <- cbind(1.0, return_matrix)
                                    res <- exp(- return_matrix %*% theta_vector - 1.0)
                                    res <- t(return_matrix) %*% (return_matrix * matrix(res, nrow = nrow(return_matrix), ncol = ncol(return_matrix)))
                                    res <- 1.0 / nrow(return_matrix) * res
                                    return(res)
                                  }
                                 , sdf_recovery = function(theta_vector, return_matrix){
                                    return_matrix <- cbind(1.0, return_matrix)
                                    res <- exp(- return_matrix %*% theta_vector - 1.0)
                                    return(res)
                                    }
                               )
                               , private = list(
                                 description = "Functions for fitting and recovering exponential-tilting SDFs. \nThey correspond to setting power = 0 in the Cressie-Read function family."
                               ))
cressie_read_functions <- R6::R6Class("cressie_read_functions"
                                      , inherit = entropy_functions
                                      , public = list(
                                        initialize = function(cressie_read_power
                                                              , sdf_mean){
                                          if(cressie_read_power %in% c(-1, 0)){
                                            stop("This function is for powers different from -1 and 0. If you wish to use those, please use et_functions or entropy_functions.")
                                          }
                                          private$cr_power <-  cressie_read_power
                                          private$sdf_average <- sdf_mean
                                        }
                                        , objective = function(theta_vector, return_matrix){
                                          sdf_expect <- private$sdf_average
                                          cr_pow <- private$cr_power
                                          # return_matrix <- return_matrix - 1.0 / sdf_expect # because excess returns
                                          res <- cr_pow * return_matrix %*% theta_vector
                                          res <- res + sdf_expect^cr_pow
                                          res <- - 1.0 / (1.0 + cr_pow) * res^((cr_pow + 1.0)/cr_pow)
                                          res <- mean(res)
                                          return(-res)
                                        }
                                        , gradient = function(theta_vector, return_matrix){
                                          sdf_expect <- private$sdf_average
                                          cr_pow <- private$cr_power
                                          # return_matrix <- return_matrix - 1.0 / sdf_expect # because excess returns
                                          res <- cr_pow * return_matrix %*% theta_vector
                                          res <- res + sdf_expect^cr_pow
                                          res <- -res^(1.0 / cr_pow)
                                          res <- apply(return_matrix, 2, function(x) res * x)
                                          res <- apply(res, 2, mean)
                                          return(-res)
                                        }
                                        , hessian = function(theta_vector, return_matrix){
                                          sdf_expect <- private$sdf_average
                                          cr_pow <- private$cr_power
                                          # return_matrix <- return_matrix - 1.0 / sdf_expect # because excess returns
                                          res <- cr_pow * return_matrix %*% theta_vector
                                          res <- res + sdf_expect^cr_pow
                                          res <- -res^((1.0 - cr_pow)/cr_pow)
                                          res <- apply(cbind(res, return_matrix), 1, function(x) x[1] * (x[-1] %*% t(x[-1])))
                                          res <- array(res, c(length(theta_vector), length(theta_vector), nrow(return_matrix)))
                                          res <- apply(res, c(1,2), mean)
                                          return(-res)
                                        }
                                        , sdf_recovery = function(theta_vector, return_matrix){
                                          sdf_expect <- private$sdf_average
                                          cr_pow <- private$cr_power
                                          # return_matrix <- return_matrix - 1.0 / sdf_expect # because excess returns
                                          res <- cr_pow * return_matrix %*% theta_vector
                                          res <- res + sdf_expect^cr_pow
                                          res <- res^(1.0 / cr_pow)
                                          res <- sdf_expect * res / mean(res)
                                          return(res)
                                        }
                                        , get_power = function(){
                                          return(private$cr_power)
                                        }
                                        , get_risk_free_return = function(){
                                          return(1.0/private$sdf_average)
                                        }
                                        , set_power = function(new_power){
                                          private$cr_power <- new_power
                                          invisible(self)
                                        }
                                        , set_risk_free_return = function(new_risk_free_return){
                                          private$sdf_average <- 1.0 / new_risk_free_return
                                          invisible(self)
                                        }
                                        , set_sdf_mean = function(new_sdf_mean){
                                          private$sdf_average <- new_sdf_mean
                                          invisible(self)
                                        }
                                      )
                                      , private = list(
                                        cr_power = -0.5 # Hellinger
                                        , sdf_average = 1.0
                                        , description = "Functions for fitting and recovering general Cressie-Read SDFs. \nThey correspond to any power different from -1 or 0 in the Cressie-Read function family."
                                      )
                                      )

distance_et_functions <- R6::R6Class("distance_et_functions"
                                     , inherit = entropy_functions
                                     , public = list(
                                        objective = function(theta_vector, return_matrix, penalty_value){
                                          
                                          # Check if there is a pre-estimated lambda (to get a starting value for the first optimisation)
                                          if(is.null(private$lambda_opt)){
                                            lambda_opt <- rep(0, ncol(return_matrix))
                                          } else {
                                            lambda_opt <- private$lambda_opt
                                          }
                                          
                                          # Optimize the inner problem given theta
                                          # These inner objective / gradient / Hessian are coded in entropic-distance-estimators.cpp
                                          def_opts <- nloptr::nl.opts()
                                          def_opts$algorithm <- "NLOPT_LD_LBFGS"
                                          optimisation_result <- nloptr::nloptr(x0 = lambda_opt
                                                                                , eval_f = et_distance_objective
                                                                                , opts = def_opts
                                                                                , theta_extended = theta_vector
                                                                                , return_matrix = return_matrix
                                                                                , mu_penalty = penalty_value)
                                          
                                          # Recover optimal decision
                                          lambda_opt <- optimisation_result$solution
                                          
                                          private$lambda_opt <- lambda_opt 
                                          
                                          # Calculate gradient wrt theta
                                          gradient_theta <- et_distance_theta_gradient(lambda_opt
                                                                            , theta_vector
                                                                            , return_matrix
                                                                            , penalty_value)
                                          
                                          # Return objective value and gradient
                                          # obj value has to switch sign, because we minimize dist for max dual,
                                          # but the dual here was minimised (that's how solvers work)
                                          # and we have to switch signs again
                                          return(list(objective = - optimisation_result$objective
                                                 , gradient = gradient_theta))
                                        }
                                        # , gradient = function(theta_vector, return_matrix, penalty_value){
                                        #   
                                        #   if(is.null(private$lambda_opt)){
                                        #     self$objective(theta_vector, return_matrix, penalty_value)$objective
                                        #     lambda_opt <- private$lambda_opt
                                        #   } else {
                                        #     lambda_opt <- private$lambda_opt
                                        #   }
                                        #   
                                        #   res <- et_distance_theta_gradient(lambda_opt
                                        #                                     , theta_vector
                                        #                                     , return_matrix
                                        #                                     , penalty_value)
                                        #   return(res)
                                        # }
                                        , get_lambda_stored = function(){
                                          private$lambda_opt
                                        }
                                        , sdf_recovery = function(theta_vector, return_matrix){
                                          res <- exp(return_matrix %*% theta_vector)
                                          return(res)
                                        }
                                     )
                                     , private = list(
                                        lambda_opt = NULL
                                        , description = "Functions for fitting and recovering exponential-tilting SDFs. \nThey correspond to setting power = 0 in the Cressie-Read function family."
                                     ))