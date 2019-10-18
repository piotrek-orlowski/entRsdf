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

xlogx_functions <- R6::R6Class("xlogx_functions"
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
                                 description = "Functions for fitting and recovering minimum-KL-distance SDFs. \nThey correspond to setting power = 0 in the Cressie-Read function family."
                               ))