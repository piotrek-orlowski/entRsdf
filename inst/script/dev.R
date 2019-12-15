library(dplyr)
library(entRsdf)
load("data/fx_portfolios.RData")

test_assets <- test_assets %>% 
  tidyr::drop_na()

# zz <- pricing_kernel$new(excess_returns = test_assets , type = "kullback-leibler")
# 
# zz$fit()
# 
# lhs <- zz$get_normalizing_constant() + as.matrix(test_assets[,-1L]) %*% zz$get_pfolio_wts()
# 
# glmtest <- glmnet::glmnet(x = as.matrix(test_assets[,-1L])
#                           , y = lhs)
# 
# 
# type <- "kullback-leibler"
# num_folds <- 3L
# excess_returns <- test_assets
# penalty_par <- exp(-seq(1,7,length.out = 100))
# 
# private$num_folds <- num_folds
# private$penalty_par <- penalty_par
# 
# zz <- cv_pricing_kernel$new(excess_returns = test_assets
#                                    , type = type
#                                    , penalty_par = penalty_par
#                                    , num_folds = 5L
#                                    , fit_full = FALSE
#                                    )
# 
# zz$fit()
# 
# ss <- lev_pricing_kernel$new(excess_returns = test_assets
#                             , type = type
#                             , penalty_par = penalty_par
#                             , num_folds = 5L
#                             , fit_full = FALSE
#                             , maximum_leverage = 50
# )
# 
# ss$fit()
# 
# 
# 
# mm <- window_cv_pricing_kernel$new(excess_returns = test_assets
#                                    , type = type
#                                    , penalty_par = penalty_par
#                                    , num_folds = 3L
#                                    , fit_full = FALSE
#                                    , sample_type = "expanding"
#                                    , sample_span = 180L)
# 
# system.time(
#   mm$fit()
# )
# 
# 
# ll <- window_lev_pricing_kernel$new(excess_returns = test_assets
#                                     , type = type
#                                     , penalty_par = penalty_par
#                                     , num_folds = 3L
#                                     , fit_full = FALSE
#                                     , sample_type = "expanding"
#                              , maximum_leverage = 50
# )
# 
# system.time(
#   ll$fit() 
# )

#### CHECK RCPP OBJECTIVE ####

et_distance_objective(lambda_exact = rep(1,2)
                      , theta_extended = rep(0,4)
                      , return_matrix = returns[,1:2]
                      , mu_penalty = 0.05)

library(numDeriv)
returns <- test_assets %>% select(-date) %>% tidyr::drop_na() %>% as.matrix()

lambda_opt <- cccp::getx(solve_inner_et_problem(theta = rep(0,4), returns = returns[,1:2], mu_penalty = 0.05))

cccp::getstate(solve_inner_et_problem(theta = rep(0,4), returns = returns[,1:2], mu_penalty = 0.05))["pobj"]

#### INNER OBJECTIVE SPEED: LBFGS vs CCCP ####

returns <- test_assets %>% select(-date) %>% tidyr::drop_na() %>% as.matrix()
theta_vector <- c(rep(1,14), rep(0,14))
return_matrix <- returns
penalty_value <- 0.005

lambda_opt <- rep(1,14)

system.time({
  for(i in 1:50){
    optimisation_result <- cccp::cccp(f0 = function(x) et_distance_objective(x
                                                                             , theta_vector
                                                                             , return_matrix
                                                                             , penalty_value)$objective
                                      , g0 = function(x) et_distance_lambda_gradient(x
                                                                                     , theta_vector
                                                                                     , return_matrix
                                                                                     , penalty_value)
                                      , h0 = function(x) et_distance_lambda_hessian(x
                                                                                    , theta_vector
                                                                                    , return_matrix
                                                                                    , penalty_value)
                                      , x0 = lambda_opt
                                      , optctrl = cccp::ctrl(trace = FALSE)
    )
}})

def_opts <- nloptr::nl.opts()
def_opts$algorithm <- "NLOPT_LD_LBFGS"

system.time({
  for(i in 1:100){
    optimisation_result <- nloptr::nloptr(x0 = lambda_opt
                                          , eval_f = et_distance_objective
                                          , opts = def_opts
                                          , theta_extended = theta_vector
                                          , mu_penalty = penalty_value
                                          , return_matrix = return_matrix
    )
  }})

#### CHECK SPEEDS ####

system.time({
  for(i in 1:10000){
    et_distance_objective(lambda_exact = lambda_opt
                          , theta_extended = rep(1,28)
                          , return_matrix = return_matrix
                          , mu_penalty = penalty_value)
  }
})

#### CHECK GRADIENT ####

foos <- distance_et_functions$new()
foos$objective(theta_vector = rep(1,4)
               , return_matrix = returns[,1:2]
               , penalty_value = 0.05)

et_distance_theta_gradient(foos$get_lambda_stored()
                           , rep(1,4)
                           , returns[,1:2]
                           , 0.05)

foos$gradient(theta_vector = rep(1,4)
              , return_matrix = returns[,1:2]
              , penalty_value = 0.05)

num_gradient <- function(theta_vector, return_matrix, penalty_value){
  foos <- distance_et_functions$new()
  foos$objective(theta_vector = theta_vector
                 , return_matrix = return_matrix
                 , penalty_value = penalty_value)$objective
  grad(func = function(...) foos$objective(...)$objective
       , x = theta_vector
       , return_matrix = return_matrix
       , penalty_value = penalty_value)
}

num_gradient(rep(1,4), returns[,1:2], penalty_value = 0.05)

grad(func = function(x) et_distance_objective(lambda_exact = foos$get_lambda_stored()
                                                 , theta_extended = x
                                                 , return_matrix = returns[,1:2]
                                                 , mu_penalty = 0.05)
        , x = rep(0,4)
)

#### TEST OF INNER GRAD_BASED OPTIMIZER ###

def_opts <- nloptr::nl.opts()
def_opts$algorithm <- "NLOPT_LD_LBFGS"


returns <- test_assets %>% select(-date) %>% tidyr::drop_na() %>% sample_frac(size=0.66) %>% as.matrix()
foos <- distance_et_functions$new()
penalty_path <- exp(seq(log(0.01), log(0.000001), length.out = 100))

par_store <- matrix(0, 100, 14)

system.time({
  outer_sol <- nloptr::nloptr(x0 = rep(0,28)
                              , eval_f = foos$objective
                              , lb = rep(0,28)
                              , opts = def_opts
                              , return_matrix = returns
                              , penalty_value = penalty_path[1]
                              )
  
  sol <- outer_sol$solution
  sol_compact <- numeric(length(sol)/2)
  sol_compact[sol[1:(length(sol)/2)]>=0] <- sol[1:(length(sol)/2)][sol[1:(length(sol)/2)]>=0]
  sol_compact[sol[(length(sol)/2+1):length(sol)]>0] <- - sol[(length(sol)/2+1):length(sol)][sol[(length(sol)/2+1):length(sol)]>0]
  
  counter <- 1
  par_store[counter,] <- sol_compact
  
  
  for(mu in penalty_path[-1]){
    outer_sol <- nloptr::nloptr(x0 = outer_sol$solution
                                , eval_f = foos$objective
                                , lb = rep(0,28)
                                , opts = def_opts
                                , return_matrix = returns
                                , penalty_value = mu
                                )
    sol <- outer_sol$solution
    sol_compact <- numeric(length(sol)/2)
    sol_compact[sol[1:(length(sol)/2)]>=0] <- sol[1:(length(sol)/2)][sol[1:(length(sol)/2)]>=0]
    sol_compact[sol[(length(sol)/2+1):length(sol)]>0] <- - sol[(length(sol)/2+1):length(sol)][sol[(length(sol)/2+1):length(sol)]>0]
    
    counter <- counter + 1
    par_store[counter, ] <- sol_compact
  }
})

system.time(
  outer_sol <- nloptr::nloptr(x0 = outer_sol$solution
                              , eval_f = foos$objective
                              , lb = rep(0,28)
                              , opts = def_opts
                              , return_matrix = returns
                              , penalty_value = 0.00004
                              , store_lambda = TRUE)
)

outer_sol <- nloptr::lbfgs(x0 = rep(1,4)
                           , fn = function(...) foos$objective(...)$objective
                           , gr = function(...) foos$objective(...)$gradient
                           , lower = rep(0,4)
                           , return_matrix = returns[,1:2]
                           , penalty_value = 0.0003
                           )


#### HESSIAN TSTS ####
# we dropped it because problem not conv in general, and cccp seems to prefer conv

foos <- distance_et_functions$new()
foos$objective(theta_vector = rep(0,4)
               , return_matrix = returns[,1:2]
               , penalty_value = 0.05)

num_hessian <- function(theta_vector, return_matrix, penalty_value){
  foos <- distance_et_functions$new()
  foos$objective(theta_vector = theta_vector
                 , return_matrix = return_matrix
                 , penalty_value = penalty_value)
  hessian(foos$objective
          , x = theta_vector
          , return_matrix = return_matrix
          , penalty_value = penalty_value)
}


foos$hessian(theta_vector = rep(0,4)
             , return_matrix = returns[,1:2]
             , penalty_value = 0.05)

num_hessian(rep(0,4), returns[,1:2], 0.05)

hessian(func = function(x) et_distance_objective(lambda_exact = foos$get_lambda_stored()
                                     , theta_extended = x
                                     , return_matrix = returns[,1:2]
                                     , mu_penalty = 0.05)
        , x = rep(0,4)
        )

# with sol
hessian(func = function(x) {
    foos$objective(x, returns[,1:2], 0.05)
    et_distance_objective(lambda_exact = foos$get_lambda_stored()
                                                 , theta_extended = x
                                                 , return_matrix = returns[,1:2]
                                                 , mu_penalty = 0.05)
  }
        , x = rep(0,4)
)

#### DEBUG CV_KERNEL ####

cv_debug <- cv_pricing_kernel$new(excess_returns = test_assets[,1:7]
                                  , type = "kullback-leibler"
                                  , penalty_par = exp(seq(log(0.01), log(0.000001), length.out = 15))
                                  , num_folds = 3L)
## debugging fit
debug(cv_debug$fit)
cv_debug$fit()
undebug(cv_debug$fit)

# tested: fit works up until evaluating criterion per fold

## debugging cv_crit
debug(cv_debug$fit)
cv_debug$fit()
debug(private$cv_criterion)

# tested: works well

## debugging assignment
debug(cv_debug$fit)
cv_debug$fit()

lev_debug <- lev_pricing_kernel$new(excess_returns = test_assets[,1:7]
                                    , type = "kullback-leibler"
                                    , penalty_par = exp(seq(log(0.01), log(0.000001), length.out = 100))
                                    , num_folds = 3L
                                    , maximum_leverage = 5)

lev_debug$fit()

#### PLAY AROUND ####

zz <- exp(seq(log(0.5), log(0.00000001), length.out = 300))

zz <- tail(zz[zz<=0.03], 200)

try(rm(cv_debug))

cv_target <- pricing_kernel$new(excess_returns = test_assets
                                  , type = "kullback-leibler")
cv_target$fit()

cv_debug <- cv_pricing_kernel$new(excess_returns = test_assets %>% head(-1)
                                  , type = "kullback-leibler"
                                  , penalty_par = zz
                                  , num_folds = 3L)

system.time(
  cv_debug$fit()  
)

#### DEV ROLLING ###

try(rm(cv_debug))

Sys.setenv(NUM_CORES = 8)

cv_target_rolling <- window_pricing_kernel$new(excess_returns = test_assets %>% head(185)
                                , type = "kullback-leibler"
                                , sample_type = "expanding"
                                , sample_span = 169)

cv_target_rolling$fit()

cv_debug <- window_cv_pricing_kernel$new(excess_returns = test_assets %>% head(185)
                                         , type = "kullback-leibler"
                                         # , penalty_par = exp(seq(log(0.01), log(0.000001), length.out = 100))
                                         , penalty_par = zz
                                         , num_folds = 3L
                                         , sample_type = "expanding"
                                         , sample_span = 169)
system.time(
  cv_debug$fit()
)


cv_debug <- window_lev_pricing_kernel$new(excess_returns = test_assets
                                         , type = "kullback-leibler"
                                         , penalty_par = exp(seq(log(0.01), log(0.000001), length.out = 10))
                                         , num_folds = 3L
                                         , sample_type = "expanding"
                                         , sample_span = 517
                                         , maximum_leverage = 20)
system.time(
  cv_debug$fit()
)

#### Multi-class ####

chi2_rolling <- window_chi2_cv_pricing_kernel$new(excess_returns = test_assets
                                                  , type = "kullback-leibler"
                                                  , penalty_par = zz
                                                  , num_folds = 3L
                                                  , sample_type = "expanding"
                                                  , sample_span = 180)

system.time(chi2_rolling$fit())


fs_rolling <- window_fs_pr_cv_pricing_kernel$new(excess_returns = test_assets %>% head(385)
                                                  , type = "kullback-leibler"
                                                  , penalty_par = zz
                                                  , num_folds = 3L
                                                  , sample_type = "expanding"
                                                  , sample_span = 373)

system.time(fs_rolling$fit())

cv_reg <- window_cv_pricing_kernel$new(excess_returns = test_assets %>% head(185)
                                                 , type = "kullback-leibler"
                                                 , penalty_par = zz
                                                 , num_folds = 3L
                                                 , sample_type = "expanding"
                                                 , sample_span = 182)

