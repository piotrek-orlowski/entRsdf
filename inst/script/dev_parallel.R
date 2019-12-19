library(dplyr)
library(entRsdf)
load("data/fx_portfolios.RData")

test_assets <- test_assets %>% 
  tidyr::drop_na()



#### TEST OF PARALLEL SPEEDUP ###

library(parallel)

cl <- makeCluster(4)
clusterEvalQ(cl, {
  library(dplyr)
  library(entRsdf)
  load("data/fx_portfolios.RData")
  
  test_assets <- test_assets %>% 
    tidyr::drop_na()
  
  def_opts <- nloptr::nl.opts()
  def_opts$algorithm <- "NLOPT_LD_LBFGS"
  
  
  returns <- test_assets %>% select(-date) %>% tidyr::drop_na() %>% sample_frac(size=0.66) %>% as.matrix()
  foos <- distance_et_functions$new()
  NULL
})


penalty_path <- exp(seq(log(0.01), log(0.000001), length.out = 100))

path_split <- split(penalty_path, ceiling(seq_along(penalty_path)/25))

### splitting penalty path into chunks to profit from warm starts seems the best from both wordls
#
# sequential : 33 sec
# sequential w/ warm start: 12 sec
# parallel chunk w/ warm start: 4 sec
system.time(
  res_list <- parLapplyLB(cl, X = path_split
              , fun = function(mu_list){
                temp_sol <- rep(0,28)
                lapply(mu_list, function(mu){
                  outer_sol <- nloptr::nloptr(x0 = temp_sol
                                              , eval_f = foos$objective
                                              , lb = rep(0,28)
                                              , opts = def_opts
                                              , return_matrix = returns
                                              , penalty_value = mu
                  )
                  temp_sol <<- outer_sol$solution
                  outer_sol$solution
                })
              })
)

system.time({
  temp_sol <- rep(0,28)
  res_list <- lapply(penalty_path, function(mu){
    outer_sol <- nloptr::nloptr(x0 = rep(0,28)
                                , eval_f = foos$objective
                                , lb = rep(0,28)
                                , opts = def_opts
                                , return_matrix = returns
                                , penalty_value = mu
                                )
    temp_sol <<- outer_sol$solution
    outer_sol$solution
    })
  })


def_opts <- nloptr::nl.opts()
def_opts$algorithm <- "NLOPT_LD_LBFGS"


returns <- test_assets %>% select(-date) %>% tidyr::drop_na() %>% sample_frac(size=0.66) %>% as.matrix()
foos <- distance_et_functions$new()


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
