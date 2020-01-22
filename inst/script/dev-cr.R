library(entRsdf)
library(dplyr)
library(numDeriv)

data("fx_portfolios")

ret_mat <- test_assets %>% .[,2:4] %>% tidyr::drop_na() %>% as.matrix()
ret_mat <- ret_mat + 1.0

sdf_expect <- 1.0
cr_pow <- -0.5

obj <- function(theta_vector, return_matrix){
  return_matrix <- return_matrix - 1.0 / sdf_expect
  res <- cr_pow * return_matrix %*% theta_vector
  res <- res + sdf_expect^cr_pow
  res <- - 1.0 / (1.0 + cr_pow) * res^((cr_pow + 1.0)/cr_pow)
  res <- mean(res)
  res
}

theta <- c(0.2,-5,2)

obj(theta, ret_mat)

gr <- function(theta_vector, return_matrix){
  return_matrix <- return_matrix - 1.0 / sdf_expect
  res <- cr_pow * return_matrix %*% theta_vector
  res <- res + sdf_expect^cr_pow
  res <- -res^(1.0 / cr_pow)
  res <- apply(return_matrix, 2, function(x) res * x)
  res <- apply(res, 2, mean)
  res
}

gr(theta, ret_mat)

grad(func = obj
     , x = theta
     , return_matrix = ret_mat) # OK!

he <- function(theta_vector, return_matrix){
  return_matrix <- return_matrix - 1.0 / sdf_expect
  res <- cr_pow * return_matrix %*% theta_vector
  res <- res + sdf_expect^cr_pow
  res <- -res^((1.0 - cr_pow)/cr_pow)
  res <- apply(cbind(res, return_matrix), 1, function(x) x[1] * (x[-1] %*% t(x[-1])))
  res <- array(res, c(length(theta_vector), length(theta_vector), nrow(return_matrix)))
  res <- apply(res, c(1,2), mean)
  res
}

he(theta, ret_mat)

hessian(func = obj
     , x = theta
     , return_matrix = ret_mat) - he(theta, ret_mat) ## ok!


hell_sdf <- pricing_kernel$new(type = "cressie-read"
                               , excess_returns = test_assets[,c(1,2,3)]
                               , cressie_read_power = -0.5
                               , sdf_mean = 1.0)

hell_sdf$get_excess_returns()
hell_sdf$fit()

#### Test with Caio's data ####
dropbox_path <- Sys.getenv("DROPBOX_PATH")
setwd(dropbox_path)

sp500_hf <- readxl::read_xlsx(path = "aago-hf-tail-risk/Files Gustavo/SPX Data.xlsx"
                              , sheet = "SPX"
                              , n_max = 200)


sp500_hf <- sp500_hf %>% filter(as.Date(Date) == "2015-11-12") %>% 
  arrange(Date) %>% 
  mutate(SPX = (lead(Fechamento) - Fechamento)/Fechamento) %>%
  tidyr::drop_na() %>%
  select(Date, SPX) %>% 
  rename(date = Date)

hell_sdf <- pricing_kernel$new(type = "cressie-read"
                               , excess_returns = sp500_hf
                               , cressie_read_power = -0.5
                               , sdf_mean = 1.0)

hell_sdf$get_excess_returns()
hell_sdf$fit()
hell_sdf$get_pfolio_wts()
