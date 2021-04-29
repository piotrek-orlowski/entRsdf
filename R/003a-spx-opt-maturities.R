# Fetch SPX option data in 5-minute intervals=

library(dplyr)
library(lubridate)
library(data.table)
library(hms)

Sys.setenv(TZ = "GMT")

library(foreach)
library(doParallel)
library(parallel)

dropbox_path <- Sys.getenv("DROPBOX_PATH")

setwd(dropbox_path)
setwd("aago-hf-tail-risk/spx-opt-returns/")

opt_file_list <- list.files("d:/datasets-preprocessed/cboe-2004-2020-option-intervals-5min/", pattern = "csv",  full.names = TRUE)

cl <- makeCluster(8)
clusterEvalQ(cl, {
  library(dplyr)
  library(lubridate)
  library(data.table)
  library(hms)
  library(RQuantLib)
  
  Sys.setenv(TZ = "GMT")
})

res <- parLapply(cl, opt_file_list, function(opt_file) {
  opt_date <- stringi::stri_extract_last_regex(str = opt_file, pattern = "\\d+-\\d+-\\d+")
  
  # Load data and select necessary columns
  opts <- data.table::fread(opt_file)
  opts <- opts %>% select(expiration) %>% distinct()
  
  # Option maturities to date with factor trick
  opts$expiration <- as.factor(opts$expiration)
  opts$expiration <- as.Date(levels(opts$expiration), format = "%Y-%m-%d")[as.integer(opts$expiration)]
  
  # Pick maturity closest to 1M
  opts <- opts[abs(as.numeric(expiration - as.Date(opt_date)) - 30) == min(abs(as.numeric(expiration - as.Date(opt_date)) - 30))]
  
  # Return only date and maturity
  tibble(date = opt_date, exdate = unique(opts$expiration))
})

res <- bind_rows(res)

data.table::fwrite(res, "res-spx-quote/option-maturities.csv")