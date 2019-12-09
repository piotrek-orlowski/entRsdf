#!/usr/bin/env Rscript

library(entRsdf)
library(dplyr)

Sys.setenv("NUM_CORES" = 12)

# args

args = commandArgs(trailingOnly=TRUE)

# make R cluster from NODELIST envir variable
# (see docs.computecanada.ca/wiki/R)
#
# The cluster will be made inside the SDF

# data
data("fx_portfolios")

print(head(test_assets))

test_assets <- test_assets %>% tidyr::drop_na()

# penalties
penalty_path <- exp(seq(log(0.5), log(0.000001), length.out = 300))

sdf_class <- args[1]
sdf_parameter <- args[2]

print("SDf clas:")
print(sdf_class)

print("sdf par")
print(sdf_parameter)

if(grepl("lev", sdf_class)){
  sdf_create_string <- sprintf("sdf <- %s$new(excess_returns = test_assets, type = 'kullback-leibler', penalty_par = penalty_path, num_folds = 3L, maximum_leverage = as.numeric(%s), sample_type = 'expanding', sample_span = 180)", sdf_class, sdf_parameter)
  } else {
    sdf_create_string <- sprintf("sdf <- %s$new(excess_returns = test_assets, type = 'kullback-leibler', penalty_par = penalty_path, num_folds = as.numeric(%s), sample_type = 'expanding', sample_span = 180)", sdf_class, sdf_parameter)
    }

eval(parse(text = sdf_create_string))

sdf$fit()

readr::write_csv(sdf$get_sdf_series(), path = sprintf("sdf-ts-%s-%s.csv", sdf_class, sdf_parameter))
save(sdf, file = sprintf("sdf-object-%s-%s.RData", sdf_class, sdf_parameter))
