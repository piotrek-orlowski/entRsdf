#!/usr/bin/env Rscript

library(entRsdf)
library(dplyr)

Sys.setenv("NUM_CORES" = 1)

# args

args = commandArgs(trailingOnly=TRUE)

# data
data("fx_portfolios")

print(head(test_assets))

test_assets <- test_assets %>% tidyr::drop_na()

# penalties
penalty_path <- exp(seq(log(0.5), log(0.000001), length.out = 300))
penalty_path <- head(penalty_path[penalty_path<=0.01], 100)

sdf_class <- args[1]
sdf_parameter <- args[2]

if(grepl("lev", sdf_class)){
  sdf_create_string <- sprintf("sdf <- %s$new(excess_returns = test_assets, type = 'kullback-leibler', penalty_par = penalty_path, num_folds = 3L, maximum_leverage = as.numeric(%s))", sdf_class, sdf_parameter)
  } else {
    sdf_create_string <- sprintf("sdf <- %s$new(excess_returns = test_assets, type = 'kullback-leibler', penalty_par = penalty_path, num_folds = as.numeric(%s))", sdf_class, sdf_parameter)
    }

eval(parse(text = sdf_create_string))

sdf$fit()

readr::write_csv(sdf$get_sdf_series(), path = sprintf("sdf-ts-%s-%s.csv", sdf_class, sdf_parameter))
save(sdf, file = sprintf("sdf-object-%s-%s.csv", sdf_class, sdf_parameter))