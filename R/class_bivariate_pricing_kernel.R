bivariate_pricing_kernel_constructor <- function(home_returns,
                                                 foreign_returns,
                                                 sdf_powers) {
  private$home_returns <- home_returns
  private$foreign_returns <- foreign_returns
  private$type <- "bivariate_pricing_kernel"

  private$fitted <- FALSE

  private$entropy_foos <- bivariate_mgf_functions$new(sdf_powers = sdf_powers)
}

bivariate_pricing_kernel <- R6::R6Class("bivariate_pricing_kernel",
  public = list(
    initialize = bivariate_pricing_kernel_constructor # Constructor, will assign all fields
    , get_returns = function() {
      private$excess_returns
    },
    get_returns_tidy = function() {
      private$excess_returns %>%
        tidyr::gather(portfolio, return, -date)
    },
    get_sdf_series = function() {
      private$sdf_series
    },
    get_type = function() {
      private$type
    },
    get_pfolio_wts = function() {
      private$pfolio_wts
    },
    get_pfolio_wts_df = function() {
      private$pfolio_wts_df
    },
    get_fitted_status = function() {
      private$fitted
    },
    get_objective = function() {

    },
    set_returns = function(home_returns, foreign_returns) {

    },
    fit = function(solver_trace = FALSE, ...) {
      home_return_matrix <- private$home_returns %>%
        select(-date) %>%
        as.matrix()

      foreign_return_matrix <- private$foreign_returns %>%
        select(-date) %>%
        as.matrix()

      wts <- solve_bivariate_entropy_problem(
        entropy_foos = private$entropy_foos,
        home_return_matrix = home_return_matrix,
        foreign_return_matrix = foreign_return_matrix,
        solver_trace = solver_trace,
        ...
      )

      home_wts <- head(wts, ncol(home_return_matrix) - 1L)
      foreign_wts <- tail(wts, ncol(foreign_return_matrix) - 1L)

      home_wts <- c(1 - sum(home_wts), home_wts)
      foreign_wts <- c(1 - sum(foreign_wts), foreign_wts)

      names(home_wts) <- colnames(home_return_matrix)
      names(foreign_wts) <- colnames(foreign_return_matrix)

      private$pfolio_wts_home <- home_wts
      private$pfolio_wts_foreign <- foreign_wts
    },
    asset_pricing = function(new_excess_returns = NULL) {

    }
  ),
  private = list(
    entropy_foos = NULL,
    home_returns = NULL,
    foreign_returns = NULL,
    type = NULL,
    sdf_series = NULL,
    pfolio_wts_home = NULL,
    pfolio_wts_foreign = NULL,
    pfolio_wts_df = NULL,
    fitted = NULL,
    objective = NULL
  )
)