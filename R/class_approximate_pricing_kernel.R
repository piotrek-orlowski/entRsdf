cv_pricing_kernel_constructor <- function(excess_returns = tibble::tibble(date = anytime::anydate(NA_real_))
                                          , type = c("kullback-leibler", "exponential-tilting", "cressie-read")
                                          , penalty_par
                                          , num_folds = 5L
                                          , fit_full = FALSE){
  
  # Initialize target SDF and fit if if requested
  private$full_sdf <- pricing_kernel$new(type = type[1L]
                                         , excess_returns = excess_returns)
  if(fit_full){
    private$full_sdf$fit()
  }
  
  # Call the super class initializer
  super$initialize(type = type
                   , excess_returns = excess_returns)
  
  # Set up fields specific to the cross-validated kernel:
  # Number of folds and penalty parameter (lambda) vector
  private$num_folds <- num_folds
  private$penalty_par <- penalty_par
}

window_cv_pricing_kernel_constructor <- function(excess_returns = tibble::tibble(date = anytime::anydate(NA_real_))
                                                 , type = c("kullback-leibler", "exponential-tilting", "cressie-read")
                                                 , penalty_par
                                                 , num_folds = 5L
                                                 , fit_full = FALSE
                                                 , sample_type = c("expanding", "rolling")
                                                 , sample_span = 180L){
  # Call the super class initializer
  super$initialize(excess_returns = excess_returns
                   , type = type
                   , penalty_par = penalty_par
                   , fit_full = fit_full
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
  
  # vector of normalizing constants of length num_dates - sample_span
  private$normalizing_constant <- matrix(NA_real_
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
                                        # Check if full SDF was fitted Thu Oct 17 22:55:49 2019 ------------------------------
                                        if(!private$full_sdf$get_fitted_status()){
                                          private$full_sdf$fit(solver_trace, ...)
                                        }
                                        # Evaluate constant + portfolio return (arg of exponent) Thu Oct 17 23:49:55 2019 ------------------------------
                                        # (watch out, this is not the same as objective value!)
                                        target_portfolio <- private$full_sdf$get_normalizing_constant() + 
                                          as.matrix(private$full_sdf$get_excess_returns() %>% select(-date)) %*% private$full_sdf$get_pfolio_wts()
                                        # Create folds in returns Thu Oct 17 23:41:01 2019 ------------------------------
                                        return_df <- private$full_sdf$get_excess_returns()
                                        set.seed(142L)
                                        return_df <- return_df %>% 
                                          # Folds are assigned by random draw from a uniform 
                                          dplyr::mutate(foldid = floor(private$num_folds * runif(n())))
                                        # Go across folds Thu Oct 17 23:53:14 2019 ------------------------------
                                        all_folds <- 0L:(private$num_folds-1L)
                                        # Fit glmnet on every fold and save to list
                                        glmnet_by_fold <- lapply(all_folds
                                                                 , function(fold){
                                                                   # for fitting you use data NOT in the fold
                                                                     glmnet::glmnet(y = target_portfolio[return_df$foldid != fold]
                                                                                    , x = return_df %>% 
                                                                                      dplyr::filter(foldid != fold) %>% 
                                                                                      dplyr::select(-date, -foldid) %>% 
                                                                                      as.matrix()
                                                                                    , lambda = private$penalty_par)
                                                                 })
                                        # Based on coefficients estimated on each fold, construct SDFs from the data that was left out from the estimation, and use them to price the fold
                                        # Each element of this list contains a vector of sums of squared pricing errors (across assets) for all lambdas
                                        # These are fold-wise cv curves
                                        glmnet_criterion_by_fold <- lapply(all_folds
                                                                      , private$cv_criterion
                                                                      , return_df = return_df
                                                                      , glmnet_by_fold = glmnet_by_fold)
                                        # put all fold-wise cv curves in a matrix
                                        # each row = cv crit values for a given lambda
                                        glmnet_criterion <- do.call(cbind, glmnet_criterion_by_fold)
                                        # average for each lambda (by row)
                                        glmnet_criterion <- apply(glmnet_criterion, 1L, mean)
                                        # pick lambda where the squared pricing errors are lowest
                                        best_penalty <- private$penalty_par[which.min(glmnet_criterion)]
                                        # estimate model for that lambda on all data
                                        glmnet_approximation <- glmnet::glmnet(y = target_portfolio
                                                                               , x = return_df %>% 
                                                                                  dplyr::select(-date, -foldid) %>% 
                                                                                  as.matrix()
                                                                               , lambda = best_penalty)
                                        # generate penalised SDF
                                        approximate_sdf_series <- private$entropy_foos$sdf_recovery(theta_vector = as.matrix(glmnet:::coef.glmnet(glmnet_approximation))
                                                                                                    , return_matrix = return_df %>% 
                                                                                                      dplyr::select(-date, -foldid) %>% 
                                                                                                      as.matrix())
                                        # normalise penalised SDF and assign to slot
                                        private$sdf_series$sdf <- approximate_sdf_series[,1L] / mean(approximate_sdf_series)
                                        
                                        # assign weights to slots
                                        private$pfolio_wts <- as.numeric(glmnet:::coef.glmnet(glmnet_approximation))[-1L]
                                        private$pfolio_wts_df$weight <- as.numeric(glmnet:::coef.glmnet(glmnet_approximation))[-1L]
                                        
                                        # update fitted status
                                        private$fitted <- TRUE
                                        
                                        invisible(self)
                                      }
                                      , get_penalty_par = function(){
                                        private$penalty_par
                                      }
                                      , set_penalty_par = function(new_penalty){
                                        private$penalty_par <- new_penalty
                                      }
                                      , use_full_sdf = function(){
                                        invisible(private$full_sdf)
                                      }
                                    )
                                    , private = list(
                                      full_sdf = NULL
                                      , num_folds = NULL
                                      , penalty_par = NULL
                                      , fold_descr = NULL
                                      , entropy_foos = NULL
                                      , cv_criterion = function(fold, return_df, glmnet_by_fold){
                                        # pick returns IN the fold for evaluating the fit
                                        return_matrix <- return_df %>% 
                                          dplyr::filter(foldid == fold) %>% 
                                          dplyr::select(-date, -foldid) %>% 
                                          as.matrix()
                                        # recover intercepts (matrix size of penalty par)
                                        loc_coefs <- glmnet_by_fold[[fold + 1L]]$a0
                                        # recover coefficients (matrix size of num coefs x penalty par)
                                        loc_coefs <- rbind(loc_coefs, glmnet_by_fold[[fold + 1L]]$beta)
                                        # for each coefficient vector, create the sdf from return sample
                                        sdf_by_lambda <- apply(X = loc_coefs
                                                               , MARGIN = 2L
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

lev_pricing_kernel <- R6::R6Class("lev_pricing_kernel"
                                  , inherit = cv_pricing_kernel
                                  , public = list(
                                    initialize = function(maximum_leverage, ...){
                                      super$initialize(...)
                                      private$maximum_leverage <- maximum_leverage
                                    }
                                  )
                                  , private = list(
                                    cv_criterion = function(fold, return_df, glmnet_by_fold){
                                      # recover portfolio coefficients (matrix size of num assets x penalty par)
                                      loc_coefs <- glmnet_by_fold[[fold + 1L]]$beta
                                      # calculate excess leverage
                                      excess_leverage <- apply(X = loc_coefs
                                                               , MARGIN = 2L
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
                                              return_df <- private$excess_returns
                                              # make period index
                                              fitting_index <- (private$sample_span + 1L):nrow(private$excess_returns)
                                              
                                              # apply over index:
                                              penalised_fits <- lapply(fitting_index, function(index_){
                                                # make window indexer
                                                fitting_window <- private$window_function(index_)
                                                # subset returns
                                                return_df <- return_df[fitting_window,]
                                                # plug them into private$full_sdf$set_excess_returns()
                                                suppressWarnings(private$full_sdf$set_excess_returns(return_df))
                                                # fit the full sdf on limited sample
                                                private$full_sdf$fit()
                                                # Recover const + portfolio
                                                target_portfolio <- private$full_sdf$get_normalizing_constant() + 
                                                  as.matrix(return_df %>% select(-date)) %*% private$full_sdf$get_pfolio_wts()
                                                # Create folds
                                                set.seed(142L)
                                                return_df <- return_df %>% 
                                                  # Folds are assigned by random draw from a uniform 
                                                  dplyr::mutate(foldid = floor(private$num_folds * runif(n())))
                                                # cv on limited sample
                                                # Go across folds Thu Oct 17 23:53:14 2019 ------------------------------
                                                all_folds <- 0L:(private$num_folds-1L)
                                                # Fit glmnet on every fold and save to list
                                                glmnet_by_fold <- lapply(all_folds
                                                                         , function(fold){
                                                                           # for fitting you use data NOT in the fold
                                                                           glmnet::glmnet(y = target_portfolio[return_df$foldid != fold]
                                                                                          , x = return_df %>% 
                                                                                            dplyr::filter(foldid != fold) %>% 
                                                                                            dplyr::select(-date, -foldid) %>% 
                                                                                            as.matrix()
                                                                                          , lambda = private$penalty_par)
                                                                         })
                                                # Based on coefficients estimated on each fold, construct SDFs from the data that was left out from the estimation, and use them to price the fold
                                                # Each element of this list contains a vector of sums of squared pricing errors (across assets) for all lambdas
                                                # These are fold-wise cv curves
                                                glmnet_criterion_by_fold <- lapply(all_folds
                                                                                   , private$cv_criterion
                                                                                   , return_df = return_df
                                                                                   , glmnet_by_fold = glmnet_by_fold)
                                                # put all fold-wise cv curves in a matrix
                                                # each row = cv crit values for a given lambda
                                                glmnet_criterion <- do.call(cbind, glmnet_criterion_by_fold)
                                                # average for each lambda (by row)
                                                glmnet_criterion <- apply(glmnet_criterion, 1L, mean)
                                                # pick lambda where the squared pricing errors are lowest
                                                best_penalty <- private$penalty_par[which.min(glmnet_criterion)]
                                                # estimate model for that lambda on all data
                                                glmnet_approximation <- glmnet::glmnet(y = target_portfolio
                                                                                       , x = return_df %>% 
                                                                                         dplyr::select(-date, -foldid) %>% 
                                                                                         as.matrix()
                                                                                       , lambda = best_penalty)
                                                
                                                # generate penalised SDF in sample
                                                approximate_sdf_series <- private$entropy_foos$sdf_recovery(theta_vector = as.matrix(glmnet:::coef.glmnet(glmnet_approximation))
                                                                                                            , return_matrix = return_df %>% 
                                                                                                                dplyr::select(-date, -foldid) %>% 
                                                                                                                as.matrix())
                                                
                                                # generate one-step ahead penalised sdf value
                                                approximate_sdf_oos <- private$entropy_foos$sdf_recovery(theta_vector = as.matrix(glmnet:::coef.glmnet(glmnet_approximation))
                                                                                                         , private$excess_returns %>% 
                                                                                                            dplyr::select(-date) %>% 
                                                                                                            as.matrix() %>% 
                                                                                                            .[index_,,drop=FALSE])
                                                
                                                # return:
                                                #   - norm const
                                                #   - wts
                                                #   - best lambda
                                                #   - skip sdf value because we will construct it later
                                                list(glmnet_coefs = glmnet:::coef.glmnet(glmnet_approximation)
                                                     , best_penalty = best_penalty
                                                     , sdf_rescale = mean(approximate_sdf_series)
                                                     , sdf_value = approximate_sdf_oos)
                                              })
                                              # assign time series of sdf, norm consts, wts
                                              approximate_sdf_oos <- sapply(penalised_fits, function(x) x$sdf_value)
                                              pfolio_wts <- sapply(penalised_fits, function(x) as.numeric(x$glmnet_coefs))
                                              
                                              private$sdf_series$sdf <- approximate_sdf_oos
                                              # be careful when assigning pfolio weights, dimension (1 + num_assets) x num_periods
                                              private$pfolio_wts[,] <- t(pfolio_wts[-1L,])
                                              private$pfolio_wts_df$weight <- as.numeric(t(pfolio_wts[-1L,]))
                                              
                                              invisible(self)
                                          })
                                        , private = list(
                                          sample_span = NULL
                                          , sample_type = NULL
                                          , window_function = NULL
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
                                          cv_criterion = function(fold, return_df, glmnet_by_fold){
                                            # recover portfolio coefficients (matrix size of num assets x penalty par)
                                            loc_coefs <- glmnet_by_fold[[fold + 1L]]$beta
                                            # calculate excess leverage
                                            excess_leverage <- apply(X = loc_coefs
                                                                     , MARGIN = 2L
                                                                     , function(vec){
                                                                       (sum(abs(vec)) - private$maximum_leverage)^2
                                                                     })}
                                          , maximum_leverage = NULL
                                        ))
