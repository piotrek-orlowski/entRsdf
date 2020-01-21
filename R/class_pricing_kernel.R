pricing_kernel_constructor <- function(type = c("kullback-leibler", "exponential-tilting", "cressie-read")
                                       , excess_returns = tibble::tibble(date = anytime::anydate(NA_real_))
                                       , ...){
  
  # Thu Oct 17 18:44:41 2019 ------------------------------
  # make sure excess returns are a wide data frame with column "date" and numeric contents
  if(!("date" %in% colnames(excess_returns))){
    stop("The provided excess_returns data.frame does not contain a 'date' column")
  }
  if(nrow(excess_returns) != length(unique(excess_returns$date))){
    stop("The provided excess_returns data.frame is either not wide, or has non-unique dates")
  }
  if(!all(sapply(excess_returns %>% dplyr::select(-date) %>% as.list, is.numeric))){
    stop("The provided excess_returns data.frame contains non-numeric data")
  }
  private$excess_returns <- excess_returns
  
  private$type <- type[1L]
  private$entropy_foos <- switch(type
                              , "kullback-leibler" = et_functions$new()
                              , "exponential-tilting" = logx_functions$new()
                              , "cressie-read" = cressie_read_functions$new(...))
  # prepare data frame to hold sdf series
  private$sdf_series <- tibble::tibble(date = excess_returns$date
                                    , sdf = NA_real_)
  # prepare numeric to hold sdf portfolio weights
  private$pfolio_wts <- matrix(NA_real_
                           , nrow = ncol(excess_returns) - 1L
                           , ncol = 1L)
  # prepare df to hold sdf portfolio weights alongside portfolio names
  private$pfolio_wts_df <- tibble::tibble(portfolio = setdiff(colnames(excess_returns), "date")
                                      , weight = rep(NA_real_, ncol(excess_returns) - 1L))
  # prepare slot for holding the normalisation constant
  private$normalizing_constant <- matrix(NA_real_
                                     , ncol = 1L
                                     , nrow = 1L)
  private$fitted <- FALSE
}

window_pricing_kernel_constructor <- function(excess_returns = tibble::tibble(date = anytime::anydate(NA_real_))
                                              , type = c("kullback-leibler", "exponential-tilting", "cressie-read")
                                              , sample_type = c("expanding", "rolling")
                                              , sample_span = 180L){
  
  if(sample_span >= nrow(excess_returns)){
    stop("The provided excess_returns data is shorter than the minimum requested sample_span")
  }
  
  # Initialize super-class Thu Oct 17 19:30:08 2019 ------------------------------
  super$initialize(type[1L], excess_returns)
  
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

pricing_kernel <- R6::R6Class("pricing_kernel"
                              , public = list(
                                initialize = pricing_kernel_constructor # Constructor, will assign all fields
                                , get_excess_returns = function(){
                                  private$excess_returns
                                }
                                , get_excess_returns_tidy = function(){
                                  private$excess_returns %>%
                                    tidyr::gather(portfolio, return, -date)
                                }
                                , get_sdf_series = function(){
                                  private$sdf_series
                                }
                                , get_type = function(){
                                  private$type
                                }
                                , get_pfolio_wts = function(){
                                  private$pfolio_wts
                                }
                                , get_pfolio_wts_df = function(){
                                  private$pfolio_wts_df
                                }
                                , get_fitted_status = function(){
                                  private$fitted
                                }
                                , get_normalizing_constant = function(){
                                  private$normalizing_constant
                                }
                                , get_objective = function(){
                                  -private$objective
                                }
                                , set_excess_returns = function(excess_returns){
                                  # A method for updating the excess return set
                                  # Check if it's the same assets
                                  same_assets <- sort(colnames(private$excess_returns)) == sort(colnames(excess_returns))
                                  if(!same_assets){
                                    warning("The new return dataset contains a different set of portfolios \nfrom the one originally supplied. \nall weight parameters will be reset (no warm start when re-fitting).")
                                  }
                                  private$excess_returns <- excess_returns
                                  private$fitted <- FALSE
                                  private$sdf_series <- tibble::tibble(date = excess_returns$date
                                                                       , sdf = NA_real_)
                                  
                                  if(!same_assets){
                                    # reainitialize weight-holding slots if a different set of assets is supplied
                                    private$pfolio_wts <- matrix(NA_real_
                                                                 , nrow = ncol(excess_returns) - 1L
                                                                 , ncol = 1L)
                                    
                                    private$pfolio_wts_df <- tibble::tibble(portfolio = setdiff(colnames(excess_returns), "date")
                                                                            , weight = rep(NA_real_, ncol(excess_returns) - 1L))
                                  }
                                  
                                  private$normalizing_constant <- matrix(NA_real_
                                                                         , ncol = 1L
                                                                         , nrow = 1L)
                                  warning("The pricing kernel was reset, call the fit method.")
                                  invisible(self)
                                }
                                , fit = function(solver_trace = FALSE, ...){
                                  if(any(is.na(private$pfolio_wts)) | any(is.na(private$normalizing_constant))){
                                    const_and_wts <- rep(1.0, ncol(private$excess_returns) - 1L) / (ncol(private$excess_returns) - 1L)
                                  } else {
                                    const_and_wts <- c(private$pfolio_wts)
                                  }
                                  # Convert returns to matrix Thu Oct 17 18:54:04 2019 ------------------------------
                                  return_matrix <- private$excess_returns %>% 
                                    dplyr::select(-date) %>% 
                                    as.matrix()
                                  
                                  # Call the fitting routine Thu Oct 17 18:56:50 2019 ------------------------------
                                  
                                  const_and_wts <- solve_entropy_problem(entropy_foos = private$entropy_foos
                                                                         , excess_return_matrix = return_matrix
                                                                         , theta_vector_init = const_and_wts
                                                                         , solver_trace = solver_trace
                                                                         , ...
                                                                         )
                                  
                                  private$pfolio_wts <- const_and_wts
                                  private$pfolio_wts_df$weight <- const_and_wts
                                  private$objective <- private$entropy_foos$objective(const_and_wts, return_matrix)
                                  private$sdf_series$sdf <- as.numeric(private$entropy_foos$sdf_recovery(const_and_wts, return_matrix))
                                  private$normalizing_constant <- mean(private$sdf_series$sdf)
                                  
                                  private$fitted <- TRUE
                                  
                                  invisible(self)
                                }
                                , asset_pricing = function(new_excess_returns = NULL){
                                  if(is.null(new_excess_returns)){
                                    new_excess_returns <- private$excess_returns
                                  }
                                  return_sdf_matrix <- private$sdf_series %>% 
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
                                entropy_foos = NULL
                                , excess_returns = NULL
                                , type = NULL
                                , sdf_series = NULL
                                , pfolio_wts = NULL
                                , pfolio_wts_df = NULL
                                , normalizing_constant = NULL
                                , fitted = NULL
                                , objective = NULL
                              )
                              )

window_pricing_kernel <- R6::R6Class("window_pricing_kernel"
                                     , inherit = pricing_kernel
                                     , public = list(
                                       initialize = window_pricing_kernel_constructor
                                       , fit = function(solver_trace = FALSE, ...){
                                         return_matrix <- private$excess_returns %>% 
                                           dplyr::select(-date) %>% 
                                           as.matrix()
                                         # Determine the LAST sample indices of fitting iterations Thu Oct 17 19:31:14 2019 -----------------------------
                                         fitting_index <- (private$sample_span + 1L):nrow(private$excess_returns)
                                         # Initialize variable for previous solution Thu Oct 17 19:55:09 2019 ------------------------------
                                         previous_solution <- rep(1.0, ncol(private$excess_returns))
                                         # Iterate on the number of fits Thu Oct 17 19:32:08 2019 ------------------------------
                                         const_and_wts <- apply(X = matrix(fitting_index)
                                                                , MARGIN = 1L
                                                                , FUN = function(index_){
                                                                  fitting_window <- private$window_function(index_)
                                                                  return_matrix <- return_matrix[fitting_window,]
                                                                  # use previous solution as starting point
                                                                  # this speeds things up by 60% on 500 x 14 return sample
                                                                  local_solution <- solve_entropy_problem(entropy_foos = private$entropy_foos
                                                                                                          , excess_return_matrix = return_matrix
                                                                                                          , theta_vector_init = previous_solution
                                                                                                          , solver_trace = solver_trace
                                                                                                          , ...)
                                                                  # save most recent solution
                                                                  previous_solution <<- local_solution
                                                                  local_solution
                                                                })
                                         # Assign to slots Thu Oct 17 21:38:25 2019 ------------------------------
                                         private$normalizing_constant[,1L] <- const_and_wts[1L,]
                                         private$pfolio_wts[,] <- t(const_and_wts[-1L,])
                                         private$objective <- apply(X = matrix(fitting_index)
                                                                    , MARGIN = 1L
                                                                    , FUN = function(index_){
                                                                      fitting_window <- private$window_function(index_)
                                                                      return_matrix <- return_matrix[fitting_window,]
                                                                      private$entropy_foos$objective(const_and_wts[,index_-private$sample_span], return_matrix)
                                                                    })
                                         private$sdf_series$sdf <- apply(X = matrix(fitting_index)
                                                                         , MARGIN = 1L
                                                                         , FUN = function(index_){
                                                                           private$entropy_foos$sdf_recovery(const_and_wts[, index_-private$sample_span], private$excess_returns[index_,] %>% select(-date) %>% as.matrix())
                                                                         })
                                         invisible(self)
                                       }
                                     )
                                     , private = list(
                                       sample_span = NULL
                                       , sample_type = NULL
                                       , window_function = NULL
                                     ))
