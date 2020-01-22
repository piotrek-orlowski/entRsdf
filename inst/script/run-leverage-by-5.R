library(entRsdf)
library(dplyr)

data("fx_portfolios")

lev_span <- seq(50,100,by=5)

sdf_list <- lapply(lev_span
                   , FUN = function(lev){

                     sdf <- window_lev_pricing_kernel$new(excess_returns = test_assets
                                                          , type = "kullback-leibler"
                                                          , sample_type = "expanding"
                                                          , sample_span = 180
                                                          , maximum_leverage = lev)

                     sdf$fit()

                     readr::write_csv(sdf$get_sdf_series(), path = sprintf("%sErik_Valeri_Piotr/Analysis/regularized-sdf/006-leverage-expanding/sdf-ts-%s-%s.csv", Sys.getenv("DROPBOX_PATH"), "lev", lev))
                     save(sdf, file = sprintf("%sErik_Valeri_Piotr/Analysis/regularized-sdf/006-leverage-expanding/sdf-object-%s-%s.RData", Sys.getenv("DROPBOX_PATH"), "lev", lev))

                                          print(sprintf("Finished lev %d", lev))
                     
                     sdf
                   })
# 
# 
# library(entRsdf)
# library(dplyr)
# 
# data("fx_portfolios")
# 
# lev_span <- seq(5,100,by=5)
# 
# sdf_list <- lapply(lev_span
#                    , FUN = function(lev){
#                      
#                      sdf <- window_lev_pricing_kernel$new(excess_returns = test_assets %>% 
#                                                             select(date, contains("carry"))
#                                                           , type = "kullback-leibler"
#                                                           , sample_type = "expanding"
#                                                           , sample_span = 180
#                                                           , maximum_leverage = lev)
#                      
#                      sdf$fit()
#                      
#                      readr::write_csv(sdf$get_sdf_series(), path = sprintf("%sErik_Valeri_Piotr/Analysis/regularized-sdf/007-leverage-carry-expanding/sdf-ts-%s-%s.csv", Sys.getenv("DROPBOX_PATH"), "lev", lev))
#                      save(sdf, file = sprintf("%sErik_Valeri_Piotr/Analysis/regularized-sdf/007-leverage-carry-expanding/sdf-object-%s-%s.RData", Sys.getenv("DROPBOX_PATH"), "lev", lev))
#                      
#                      print(sprintf("Finished lev %d", lev))
#                      
#                      sdf
#                    })
