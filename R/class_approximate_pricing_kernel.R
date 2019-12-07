cv_pricing_kernel_constructor <- function(excess_returns = tibble::tibble(date = anytime::anydate(NA_real_))
                                          , type = c("kullback-leibler", "exponential-tilting", "cressie-read")
                                          , penalty_par
                                          , num_folds = 5L){
  
  # Call the super class initializer
  super$initialize(type = type
                   , excess_returns = excess_returns)
  
  # complementary sdf series and weights
  private$full_sdf_series <- private$sdf_series
  private$complementary_pfolio_wts <- private$pfolio_wts
  private$complementary_pfolio_wts_df <- private$pfolio_wts_df
  
  
  # set up functions
  private$entropy_foos <- switch(type
                                 , "kullback-leibler" = distance_et_functions$new()
                                 , "exponential-tilting" = distance_el_functions$new()
                                 , "cressie-read" = distance_cressie_read_functions$new())
  
  private$best_penalty_par <- NA_real_
  
  # Set up fields specific to the cross-validated kernel:
  # Number of folds and penalty parameter (lambda) vector
  private$num_folds <- num_folds
  private$penalty_par <- penalty_par
}

window_cv_pricing_kernel_constructor <- function(excess_returns = tibble::tibble(date = anytime::anydate(NA_real_))
                                                 , type = c("kullback-leibler", "exponential-tilting", "cressie-read")
                                                 , penalty_par
                                                 , num_folds = 5L
                                                 , sample_type = c("expanding", "rolling")
                                                 , sample_span = 180L){
  # Call the super class initializer
  super$initialize(excess_returns = excess_returns
                   , type = type
                   , penalty_par = penalty_par
                   , num_folds = num_folds
                   )
  
  # Fill in remaining elements Thu Oct 17 19:36:24 2019 ------------------------------
  private$sample_type <- sample_type[1L]
  private$sample_span <- sample_span
  
  # In this case the SDF data frame should only hold the series after the initial estimation period of length sample_span
  private$sdf_series <- tibble::tibble(date = excess_returns$date[-(1L:sample_span)]
                                       , sdf = NA_real_)
  
  # We have different weights for every date, and we have to pack them into a df that is long num_assets x (num_dates - sample_span) (the latter for eliminating the initial estimation period)
  private$pfolio_wts_df <- expand.grid(date = unique(excess_returns$date[-(1L:sample_span)])
                                       , portfolio = setdiff(colnames(excess_returns), "date")
                                       , weight = NA_real_) %>% 
    tibble::as_tibble()
  
  # matrix for portfolio weights, again num_dates - sample_span rows, num_assets columns
  private$pfolio_wts <- matrix(NA_real_
                               , nrow = nrow(excess_returns) - sample_span
                               , ncol = ncol(excess_returns) - 1L)
  
  # complementary sdf series and weights
  private$full_sdf_series <- private$sdf_series
  private$complementary_pfolio_wts <- private$pfolio_wts
  private$complementary_pfolio_wts_df <- private$pfolio_wts_df
  
  # vector of normalizing constants of length num_dates - sample_span
  private$normalizing_constant <- matrix(NA_real_
                                         , nrow = nrow(excess_returns) - sample_span
                                         , ncol = 1L)
  
  private$best_penalty_par <- matrix(NA_real_
                                 , nrow = nrow(excess_returns) - sample_span
                                 , ncol = 1L)
  
  # function for determining indexing vector in rolling and expanding window
  private$window_function <- ifelse(test = sample_type == "rolling"
                                    , yes = function(index_){
                                      (index_ - private$sample_span):(index_ - 1L)
                                    }
                                    , no = function(index_){
                                      1L:(index_ - 1L)
                                    })
  
}

cv_pricing_kernel <- R6::R6Class("cv_pricing_kernel"
                                    , inherit = pricing_kernel
                                    , public = list(
                                      initialize = cv_pricing_kernel_constructor
                                      , fit = function(solver_trace = FALSE, ...){
                                        
                                        # make cluster and first check environment variable for its size
                                        if(is.na(as.numeric(Sys.getenv("NUM_CORES")))){
                                          par_cluster <- parallel::makeCluster(detectCores(TRUE))  
                                        } else {
                                          par_cluster <- parallel::makeCluster(as.numeric(Sys.getenv("NUM_CORES")))
                                        }
                                        sink <- parallel::clusterEvalQ(par_cluster
                                                     , {
                                                       library(entRsdf)
                                                       if(require("RevoUtilsMath")){
                                                         setMKLthreads(1L) 
                                                       }
                                                       NULL
                                                     })
                                        
                                        # Fri Dec 06 14:04:05 2019 ------------------------------
                                        # set up unpenalized pricing kernel which will be a reference for the CV
                                        full_pricing_kernel <- pricing_kernel$new(type = super$get_type()
                                                                                  , excess_returns = super$get_excess_returns())
                                        full_pricing_kernel$fit()
                                        
                                        # Create folds in returns Thu Oct 17 23:41:01 2019 ------------------------------
                                        return_df <- self$get_excess_returns()
                                        set.seed(142L)
                                        return_df <- return_df %>% 
                                          # Folds are assigned by random draw from a uniform 
                                          # dplyr::mutate(foldid = floor(private$num_folds * runif(n())))
                                          # Or consecutively because random does not work for large number of folds
                                          dplyr::mutate(foldid = floor(1:n()/ceiling(n()/private$num_folds)))
                                        # Go across folds Thu Oct 17 23:53:14 2019 ------------------------------
                                        all_folds <- 0L:(private$num_folds-1L)
                                        
                                        # Set up options for optimizer
                                        def_opts <- nloptr::nl.opts()
                                        def_opts$algorithm <- "NLOPT_LD_LBFGS"
                                        
                                        # for each fold, and along the penalty path, fit SDFs and save:
                                        # 1) theta
                                        # 2) lambda for each theta
                                        coefficients_by_fold <- lapply(all_folds
                                         , function(fold){
                                           # for fitting you use data NOT in the fold
                                           return_matrix <- return_df %>% 
                                             dplyr::filter(foldid != fold) %>% 
                                             dplyr::select(-date, -foldid) %>% 
                                             as.matrix()
                                           
                                           theta_extended_init <- rep(0.0, 2L * ncol(return_matrix))
                                           
                                           # matrices for storing thetas and lambdas
                                           theta_store <- matrix(0, length(private$penalty_par), ncol(return_matrix))
                                           lambda_store <- matrix(0, length(private$penalty_par), ncol(return_matrix))
                                           
                                           # make a copy of the foos object for export to cluster
                                           foos_copy <- private$entropy_foos$clone(deep = TRUE)
                                           
                                           # export the necessary objects to the cluster:
                                           #  - foos_copy
                                           #  - theta_extended_init
                                           #  - def_opts
                                           #  - return_matrix
                                           parallel::clusterExport(par_cluster, c("foos_copy"
                                                                                  , "theta_extended_init"
                                                                                  , "def_opts"
                                                                                  , "return_matrix")
                                                                   , envir = environment())
                                           
                                           penalty_split <- split(private$penalty_par, ceiling(seq_along(private$penalty_par)/ceiling(length(private$penalty_par)/length(par_cluster))))
                                           
                                           optimisation_list <- parallel::parLapplyLB(
                                             cl = par_cluster
                                             , X = penalty_split
                                             , fun = function(mu_list){
                                                temp_sol <- rep(0, length(theta_extended_init))
                                                lapply(mu_list, function(mu) {
                                                  outer_sol <- tryCatch(
                                                      nloptr::nloptr(x0 = temp_sol
                                                                              , eval_f = foos_copy$objective
                                                                              , lb = rep(0, length(temp_sol))
                                                                              , opts = def_opts
                                                                              , return_matrix = return_matrix
                                                                              , penalty_value = mu)
                                                  , error = function(e) rep(NA_real_, length(temp_sol))
                                                  )
                                                  if(!any(is.na(outer_sol))){
                                                    temp_sol <<- outer_sol$solution  
                                                  }
                                                  sol_packed <- outer_sol$solution[1L:ncol(return_matrix)]
                                                  sol_packed <- sol_packed - outer_sol$solution[(ncol(return_matrix)+1L):(2*ncol(return_matrix))]
                                                  rbind(theta_packed = sol_packed
                                                        , lambda_full = foos_copy$get_lambda_stored())
                                                }
                                                )
                                               })
                                            
                                            optimisation_coeffs <- abind::abind(lapply(optimisation_list, abind::abind, along = 3L), along = 3L)
                                            
                                            theta_store[,] <- t(optimisation_coeffs[1L,,])
                                            lambda_store[,] <- t(optimisation_coeffs[2L,,])
                                             
                                           return(list(theta_compact_matrix = theta_store
                                                       , lambda_matrix = lambda_store))
                                                                         })
                                        # Based on coefficients estimated on each fold, construct SDFs from the data that was left out from the estimation, and use them to price the fold
                                        # Each element of this list contains a vector of sums of squared pricing errors (across assets) for all lambdas
                                        # These are fold-wise cv curves
                                        cv_criterion_by_fold <- lapply(all_folds
                                                                      , private$cv_criterion
                                                                      , return_df = return_df
                                                                      , coefficients_by_fold = coefficients_by_fold
                                                                      , cv_target = full_pricing_kernel$get_sdf_series())
                                        
                                        # put all fold-wise cv curves in a matrix
                                        # each row = cv crit values for a given lambda
                                        cv_criterion <- do.call(cbind, cv_criterion_by_fold)
                                        
                                        # average for each penalty (by row)
                                        cv_criterion <- apply(cv_criterion, 1L, mean, na.rm=TRUE)
                                        
                                        # pick penalty where the squared pricing errors are lowest
                                        best_penalty <- private$penalty_par[which.min(cv_criterion)]
                                        
                                        # estimate model for that penalty on all data
                                        # start from average theta for that penalty
                                        average_theta <- sapply(coefficients_by_fold, function(cf_list){
                                          cf_list$theta_compact_matrix[which.min(cv_criterion), ]
                                        })
                                        average_theta <- apply(average_theta, 1L, mean, na.rm=TRUE)
                                        approximate_sdf_theta <- nloptr::nloptr(x0 = private$theta_unpack(-full_pricing_kernel$get_pfolio_wts())
                                                                                , eval_f = private$entropy_foos$objective
                                                                                , lb = rep(0,2L*length(average_theta))
                                                                                , opts = def_opts
                                                                                , return_matrix = return_df %>% 
                                                                                  dplyr::select(-date, -foldid) %>% 
                                                                                  as.matrix()
                                                                                , penalty_value = best_penalty)$solution
                                        
                                        # generate penalised SDF
                                        approximate_sdf_series <- private$entropy_foos$sdf_recovery(theta_vector = private$theta_pack(approximate_sdf_theta)
                                                                                                    , return_matrix = return_df %>% 
                                                                                                      dplyr::select(-date, -foldid) %>% 
                                                                                                      as.matrix())
                                        # normalise penalised SDF and assign to slot
                                        private$sdf_series$sdf <- approximate_sdf_series[,1L] # / mean(approximate_sdf_series)
                                        
                                        # assign weights to slots
                                        private$pfolio_wts <- private$theta_pack(approximate_sdf_theta)
                                        private$pfolio_wts_df$weight <- private$theta_pack(approximate_sdf_theta)
                                        
                                        # assign best lambdas
                                        private$complementary_pfolio_wts <- private$entropy_foos$get_lambda_stored()
                                        private$complementary_pfolio_wts_df$weight <- private$entropy_foos$get_lambda_stored()
                                        
                                        # assign best penalty
                                        private$best_penalty_par <- best_penalty
                                        
                                        # evaluate full SDF and assign
                                        full_sdf_series <- approximate_sdf_series - 1.0 + 
                                          private$entropy_foos$sdf_recovery(theta_vector = private$entropy_foos$get_lambda_stored()
                                                                            , return_matrix = return_df %>% 
                                                                              dplyr::select(-date, -foldid) %>% 
                                                                              as.matrix())
                                        
                                        private$full_sdf_series <- tibble::tibble(date = private$sdf_series$date, sdf = as.numeric(full_sdf_series))
                                        
                                        # update fitted status
                                        private$fitted <- TRUE
                                        
                                        # stop cluster
                                        parallel::stopCluster(par_cluster)
                                        
                                        invisible(self)
                                      }
                                      , get_penalty_par = function(){
                                        private$penalty_par
                                      }
                                      , set_penalty_par = function(new_penalty){
                                        private$penalty_par <- new_penalty
                                      }
                                      , get_best_penalty = function(){
                                          private$best_penalty_par
                                      }
                                      , get_complementary_pf_wts = function(){
                                          private$complementary_pfolio_wts
                                      }
                                      , get_full_sdf_series = function(){
                                        private$full_sdf_series
                                      }
                                      , full_asset_pricing = function(new_excess_returns = NULL){
                                        if(is.null(new_excess_returns)){
                                          new_excess_returns <- private$excess_returns
                                        }
                                        return_sdf_matrix <- private$full_sdf_series %>% 
                                          dplyr::inner_join(new_excess_returns) %>% 
                                          dplyr::select(-date) %>% 
                                          as.matrix()
                                        pricing_error <- apply(X = return_sdf_matrix[,-1L]
                                                               , MARGIN = 2L
                                                               , FUN = function(x) 1.0 / length(x) * crossprod(x, return_sdf_matrix[,1L]))
                                        pricing_error <- tibble::tibble(portfolio = setdiff(colnames(new_excess_returns), "date")
                                                                        , pricing_error = pricing_error)
                                        return(pricing_error)
                                      }
                                    )
                                    , private = list(
                                      full_sdf_series = NULL
                                      , num_folds = NULL
                                      , penalty_par = NULL
                                      , best_penalty_par = NULL
                                      , fold_descr = NULL
                                      , entropy_foos = NULL
                                      , complementary_pfolio_wts = NULL
                                      , complementary_pfolio_wts_df = NULL
                                      , cv_criterion = function(fold, return_df, coefficients_by_fold, cv_target){
                                        # pick returns IN the fold for evaluating the fit
                                        return_fold <- return_df %>% 
                                          dplyr::filter(foldid == fold)
                                        
                                        return_matrix <- return_fold %>% 
                                          dplyr::select(-date, -foldid) %>% 
                                          as.matrix()
                                        
                                        theta_matrix <- coefficients_by_fold[[fold+1L]]$theta_compact_matrix
                                        # theta_matrix <- apply(theta_matrix, 1L, private$theta_unpack)
                                        # for each coefficient vector, create the sdf from return sample
                                        sdf_by_lambda <- apply(X = theta_matrix
                                                               , MARGIN = 1L
                                                               , FUN = private$entropy_foos$sdf_recovery
                                                               , return_matrix = return_matrix)
                                        # for each penalised SDF in the fold, check how far it is from full-sample unpenalised SDF
                                        cv_target <- cv_target %>% 
                                          dplyr::inner_join(return_fold %>% select(date))
                                        # for each penalised sdf in the fold, evaluate pricing errors
                                        squared_fit_error <- apply(X = sdf_by_lambda
                                                                       , MARGIN = 2L
                                                                       , FUN = function(sdf){
                                                                         sum(abs(sdf - cv_target$sdf))
                                                                       })
                                        squared_fit_error
                                      }
                                      , theta_unpack = function(theta_compact){
                                          num_par <- length(theta_compact)
                                          theta_extended <- numeric(2L * num_par)
                                          theta_extended[1L:num_par] <- pmax(theta_compact, 0.0)
                                          theta_extended[(num_par+1L):(2L*num_par)] <- - pmin(theta_compact, 0.0)
                                          theta_extended
                                      }
                                      , theta_pack = function(theta_extended){
                                          num_par <- length(theta_extended) / 2L
                                          theta_compact <- numeric(num_par)
                                          theta_compact[] <- theta_extended[1L:num_par] - theta_extended[(num_par+1L):(2L*num_par)]
                                          theta_compact
                                      }
                                    ))

lev_pricing_kernel <- R6::R6Class("lev_pricing_kernel"
                                  , inherit = cv_pricing_kernel
                                  , public = list(
                                    initialize = function(maximum_leverage, ...){
                                      super$initialize(...)
                                      private$maximum_leverage <- maximum_leverage
                                    }
                                  )
                                  , private = list(
                                    cv_criterion = function(fold, return_df, coefficients_by_fold){
                                      # recover portfolio coefficients (matrix size of num assets x penalty par)
                                      loc_coefs <- coefficients_by_fold[[fold + 1L]]$theta_compact_matrix
                                      # calculate excess leverage
                                      excess_leverage <- apply(X = loc_coefs
                                                               , MARGIN = 1L
                                                               , function(vec){
                                                                 (sum(abs(vec)) - private$maximum_leverage)^2
                                                               })}
                                    , maximum_leverage = NULL
                                  ))

window_cv_pricing_kernel <- R6::R6Class("window_cv_pricing_kernel"
                                        , inherit = cv_pricing_kernel
                                        , public = list(
                                          initialize = window_cv_pricing_kernel_constructor
                                          , fit = function(solver_trace = FALSE, ...){
                                              
                                              # make cluster and first check environment variable for its size
                                              if(is.na(as.numeric(Sys.getenv("NUM_CORES")))){
                                                par_cluster <- parallel::makeCluster(parallel::detectCores(TRUE))  
                                              } else {
                                                par_cluster <- parallel::makeCluster(as.numeric(Sys.getenv("NUM_CORES")))
                                              }
                                              parallel::clusterEvalQ(par_cluster
                                                           , {
                                                             library(entRsdf)
                                                             if(require("RevoUtilsMath")){
                                                               setMKLthreads(1L) 
                                                             }
                                                             NULL
                                                           })
                                            
                                              return_df <- self$get_excess_returns()
                                              # make period index
                                              fitting_index <- (private$sample_span + 1L):nrow(private$excess_returns)
                                              
                                              # apply over index:
                                              penalised_fits <- lapply(fitting_index, function(index_){
                                                # make window indexer
                                                fitting_window <- private$window_function(index_)
                                                # subset returns
                                                return_df <- return_df[fitting_window,]
                                                
                                                # Fri Dec 06 14:04:05 2019 ------------------------------
                                                # set up unpenalized pricing kernel which will be a reference for the CV
                                                full_pricing_kernel <- pricing_kernel$new(type = super$get_type()
                                                                                          , excess_returns = return_df)
                                                full_pricing_kernel$fit()
                                                
                                                # Set up options for optimizer
                                                def_opts <- nloptr::nl.opts()
                                                def_opts$algorithm <- "NLOPT_LD_LBFGS"
                                                # Create folds
                                                set.seed(142L)
                                                return_df <- return_df %>% 
                                                  # Folds are assigned by random draw from a uniform 
                                                  dplyr::mutate(foldid = floor(private$num_folds * runif(n())))
                                                # cv on limited sample
                                                # Go across folds Thu Oct 17 23:53:14 2019 ------------------------------
                                                all_folds <- 0L:(private$num_folds-1L)
                                                # Fit glmnet on every fold and save to list
                                                coefficients_by_fold <- lapply(all_folds
                                                   , function(fold){
                                                     # for fitting you use data NOT in the fold
                                                     return_matrix <- return_df %>% 
                                                       dplyr::filter(foldid != fold) %>% 
                                                       dplyr::select(-date, -foldid) %>% 
                                                       as.matrix()
                                                     
                                                     theta_extended_init <- rep(0.0, 2L * ncol(return_matrix))
                                                     
                                                     # matrices for storing thetas and lambdas
                                                     theta_store <- matrix(0, length(private$penalty_par), ncol(return_matrix))
                                                     lambda_store <- matrix(0, length(private$penalty_par), ncol(return_matrix))
                                                     
                                                     # make a copy of the foos object for export to cluster
                                                     foos_copy <- private$entropy_foos$clone(deep = TRUE)
                                                     
                                                     # export the necessary objects to the cluster:
                                                     #  - foos_copy
                                                     #  - theta_extended_init
                                                     #  - def_opts
                                                     #  - return_matrix
                                                     parallel::clusterExport(par_cluster, c("foos_copy"
                                                                                            , "theta_extended_init"
                                                                                            , "def_opts"
                                                                                            , "return_matrix")
                                                                             , envir = environment())
                                                     
                                                     penalty_split <- split(private$penalty_par, ceiling(seq_along(private$penalty_par)/ceiling(length(private$penalty_par)/length(par_cluster))))
                                                     
                                                     optimisation_list <- parallel::parLapplyLB(
                                                       cl = par_cluster
                                                       , X = penalty_split
                                                       , fun = function(mu_list){
                                                         temp_sol <- rep(0, length(theta_extended_init))
                                                         lapply(mu_list, function(mu) {
                                                           outer_sol <- tryCatch(
                                                             nloptr::nloptr(x0 = temp_sol
                                                                                         , eval_f = foos_copy$objective
                                                                                         , lb = rep(0, length(temp_sol))
                                                                                         , opts = def_opts
                                                                                         , return_matrix = return_matrix
                                                                                         , penalty_value = mu)
                                                             , error = function(e) list(solution=rep(NA_real_, length(temp_sol)))
                                                           )
                                                           if(!any(is.na(outer_sol$solution))){
                                                             temp_sol <<- outer_sol$solution  
                                                           }
                                                           sol_packed <- outer_sol$solution[1L:ncol(return_matrix)]
                                                           sol_packed <- sol_packed - outer_sol$solution[(ncol(return_matrix)+1L):(2*ncol(return_matrix))]
                                                           rbind(theta_packed = sol_packed
                                                                 , lambda_full = foos_copy$get_lambda_stored())
                                                         }
                                                         )
                                                       })
                                                     
                                                     optimisation_coeffs <- abind::abind(lapply(optimisation_list, abind::abind, along = 3L), along = 3L)
                                                     
                                                     theta_store[,] <- t(optimisation_coeffs[1L,,])
                                                     lambda_store[,] <- t(optimisation_coeffs[2L,,])
                                                     
                                                     return(list(theta_compact_matrix = theta_store
                                                                 , lambda_matrix = lambda_store))
                                                                               })
                                                # Based on coefficients estimated on each fold, construct SDFs from the data that was left out from the estimation, and use them to price the fold
                                                # Each element of this list contains a vector of sums of squared pricing errors (across assets) for all lambdas
                                                # These are fold-wise cv curves
                                                cv_criterion_by_fold <- lapply(all_folds
                                                                               , private$cv_criterion
                                                                               , return_df = return_df
                                                                               , coefficients_by_fold = coefficients_by_fold
                                                                               , cv_target = full_pricing_kernel$get_sdf_series())
                                                # put all fold-wise cv curves in a matrix
                                                # each row = cv crit values for a given lambda
                                                cv_criterion <- do.call(cbind, cv_criterion_by_fold)
                                                
                                                # average for each penalty (by row)
                                                cv_criterion <- apply(cv_criterion, 1L, mean, na.rm=TRUE)
                                                
                                                # pick penalty where the squared pricing errors are lowest
                                                best_penalty <- private$penalty_par[which.min(cv_criterion)]
                                                
                                                # estimate model for that penalty on all data
                                                # start from average theta for that penalty
                                                average_theta <- sapply(coefficients_by_fold, function(cf_list){
                                                  cf_list$theta_compact_matrix[which.min(cv_criterion), ]
                                                })
                                                average_theta <- apply(average_theta, 1L, mean, na.rm=TRUE)
                                                approximate_sdf_theta <- nloptr::nloptr(x0 = private$theta_unpack(average_theta)
                                                                                        , eval_f = private$entropy_foos$objective
                                                                                        , lb = rep(0,2L*length(average_theta))
                                                                                        , opts = def_opts
                                                                                        , return_matrix = return_df %>% 
                                                                                          dplyr::select(-date, -foldid) %>% 
                                                                                          as.matrix()
                                                                                        , penalty_value = best_penalty)$solution
                                                
                                                # generate penalised SDF
                                                approximate_sdf_series <- private$entropy_foos$sdf_recovery(theta_vector = private$theta_pack(approximate_sdf_theta)
                                                                                                            , return_matrix = self$get_excess_returns() %>% 
                                                                                                              dplyr::select(-date) %>% 
                                                                                                              as.matrix() %>% 
                                                                                                              .[c(fitting_window, index_),])
                                                
                                                # assign best lambdas
                                                complementary_pfolio_wts <- private$entropy_foos$get_lambda_stored()
                                                
                                                
                                                # evaluate full SDF and assign
                                                full_sdf_series <- approximate_sdf_series - 1.0 + 
                                                  private$entropy_foos$sdf_recovery(theta_vector = private$entropy_foos$get_lambda_stored()
                                                                                    , return_matrix = self$get_excess_returns() %>% 
                                                                                      dplyr::select(-date) %>% 
                                                                                      as.matrix() %>% 
                                                                                      .[c(fitting_window, index_),])
                                                
                                                
                                                # return:
                                                #   - norm const
                                                #   - wts
                                                #   - best lambda
                                                #   - skip sdf value because we will construct it later
                                                list(theta_vector = private$theta_pack(approximate_sdf_theta)
                                                     , lambda_vector = private$entropy_foos$get_lambda_stored()
                                                     , best_penalty = best_penalty
                                                     , sdf_rescale = mean(head(approximate_sdf_series,-1), na.rm=TRUE)
                                                     , sdf_value = tail(approximate_sdf_series, 1L)
                                                     , full_sdf_value = tail(full_sdf_series, 1L)
                                                     , full_sdf_rescale = mean(head(approximate_sdf_series,-1)), na.rm=TRUE)
                                              })
                                              
                                              # assign time series of sdf, norm consts, wts
                                              approximate_sdf_oos <- sapply(penalised_fits, function(x) x$sdf_value)
                                              
                                              # assigne time series of full sdf
                                              full_sdf_oos <- sapply(penalised_fits, function(x) x$full_sdf_value)
                                              private$full_sdf_series <- tibble::tibble(date = private$sdf_series$date, sdf = as.numeric(full_sdf_oos))
                                              
                                              pfolio_wts <- sapply(penalised_fits, function(x) as.numeric(x$theta_vector))
                                              
                                              private$sdf_series$sdf <- approximate_sdf_oos
                                              
                                              # be careful when assigning pfolio weights, dimension (1 + num_assets) x num_periods
                                              private$pfolio_wts[,] <- t(pfolio_wts)
                                              private$pfolio_wts_df$weight <- as.numeric(t(pfolio_wts))
                                              
                                              # be careful when assigning pfolio weights, dimension (1 + num_assets) x num_periods
                                              # private$complementary_pfolio_wts <- private$pfolio_wts
                                              complementary_pfolio_wts <- sapply(penalised_fits, function(x) as.numeric(x$lambda_vector))
                                              private$complementary_pfolio_wts[,] <- t(complementary_pfolio_wts)
                                              # private$complementary_pfolio_wts_df <- private$pfolio_wts_df
                                              private$complementary_pfolio_wts_df$weight <-  as.numeric(t(complementary_pfolio_wts))
                                              
                                              # penalties
                                              private$best_penalty_par[,] <- sapply(penalised_fits, function(x) x$best_penalty)
                                              
                                              # stop cluster
                                              parallel::stopCluster(par_cluster)
                                              
                                              invisible(self)
                                          })
                                        , private = list(
                                          sample_span = NULL
                                          , sample_type = NULL
                                          , window_function = NULL
                                        ))

fs_pr_cv_pricing_kernel <- R6::R6Class("fs_pr_cv_pricing_kernel"
                                       , inherit = cv_pricing_kernel
                                       , private = list(
                                         cv_criterion = function(fold, return_df, coefficients_by_fold, cv_target){
                                           # evaluate the fit on whole sample to have a more precise estimate
                                           return_matrix <- return_df %>% 
                                             # dplyr::filter(foldid == fold) %>% 
                                             dplyr::select(-date, -foldid) %>% 
                                             as.matrix()
                                           
                                           theta_matrix <- coefficients_by_fold[[fold+1L]]$theta_compact_matrix
                                           # theta_matrix <- apply(theta_matrix, 1L, private$theta_unpack)
                                           # for each coefficient vector, create the sdf from return sample
                                           sdf_by_lambda <- apply(X = theta_matrix
                                                                  , MARGIN = 1L
                                                                  , FUN = private$entropy_foos$sdf_recovery
                                                                  , return_matrix = return_matrix)
                                           # for each penalised sdf in the fold, evaluate pricing errors
                                           squared_pricing_error <- apply(X = sdf_by_lambda
                                                                          , MARGIN = 2L
                                                                          , FUN = function(sdf){
                                                                            res <- apply(X = return_matrix
                                                                                         , MARGIN = 2L
                                                                                         , FUN = function(ret){
                                                                                           1.0 / length(ret) * crossprod(ret, sdf) # crossprod is fastest for average of products
                                                                                         })
                                                                            # here we tried to put the discrepancy of sdf from 1 as criterion, but that does not work out too well
                                                                            # res <- c(res, mean(sdf - 1.0))
                                                                            sum(res^2)
                                                                          })
                                           squared_pricing_error
                                         }
                                       ))

fs_pr_cv_pricing_kernel <- R6::R6Class("fs_pr_cv_pricing_kernel"
                                       , inherit = cv_pricing_kernel
                                       , private = list(
                                         cv_criterion = function(fold, return_df, coefficients_by_fold, cv_target){
                                           # evaluate the fit on whole sample to have a more precise estimate
                                           return_matrix <- return_df %>% 
                                             # dplyr::filter(foldid == fold) %>% 
                                             dplyr::select(-date, -foldid) %>% 
                                             as.matrix()
                                           
                                           theta_matrix <- coefficients_by_fold[[fold+1L]]$theta_compact_matrix
                                           # theta_matrix <- apply(theta_matrix, 1L, private$theta_unpack)
                                           # for each coefficient vector, create the sdf from return sample
                                           sdf_by_lambda <- apply(X = theta_matrix
                                                                  , MARGIN = 1L
                                                                  , FUN = private$entropy_foos$sdf_recovery
                                                                  , return_matrix = return_matrix)
                                           # for each penalised sdf in the fold, evaluate pricing errors
                                           squared_pricing_error <- apply(X = sdf_by_lambda
                                                                          , MARGIN = 2L
                                                                          , FUN = function(sdf){
                                                                            res <- apply(X = return_matrix
                                                                                         , MARGIN = 2L
                                                                                         , FUN = function(ret){
                                                                                           1.0 / length(ret) * crossprod(ret, sdf) # crossprod is fastest for average of products
                                                                                         })
                                                                            # here we tried to put the discrepancy of sdf from 1 as criterion, but that does not work out too well
                                                                            # res <- c(res, mean(sdf - 1.0))
                                                                            sum(res^2)
                                                                          })
                                           squared_pricing_error
                                         }
                                       ))


window_lev_pricing_kernel = R6::R6Class("window_lev_pricing_kernel"
                                        , inherit = window_cv_pricing_kernel
                                        , public = list(
                                          initialize = function(maximum_leverage, ...){
                                            super$initialize(...)
                                            private$maximum_leverage <- maximum_leverage
                                          }
                                        )
                                        , private = list(
                                          cv_criterion = function(fold, return_df, coefficients_by_fold){
                                            # recover portfolio coefficients (matrix size of num assets x penalty par)
                                            loc_coefs <- coefficients_by_fold[[fold+1]]$theta_compact_matrix
                                            
                                            # calculate excess leverage
                                            excess_leverage <- apply(X = loc_coefs
                                                                     , MARGIN = 1L
                                                                     , function(vec){
                                                                       (sum(abs(vec)) - private$maximum_leverage)^2
                                                                     })}
                                          , maximum_leverage = NULL
                                        ))
