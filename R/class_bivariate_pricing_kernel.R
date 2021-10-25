bivariate_pricing_kernel_constructor <- function(home_returns,
                                                 foreign_returns,
                                                 sdf_powers) {
  # Make sure dates are aligned
  home_dates <- home_returns %>% dplyr::select(date)
  foreign_dates <- foreign_returns %>% dplyr::select(date)
  
  suppressMessages({
    common_dates <- home_dates %>% dplyr::inner_join(foreign_dates)

    home_returns <- home_returns %>% 
    dplyr::inner_join(common_dates)

  foreign_returns <- foreign_returns %>% 
    dplyr::inner_join(common_dates)
  })

  # Assign returns
  private$home_returns <- home_returns
  private$foreign_returns <- foreign_returns
  private$type <- "bivariate_pricing_kernel"

  private$fitted <- FALSE

  # Set powers and create functions
  private$entropy_foos <- bivariate_mgf_functions$new(sdf_powers = sdf_powers)
}

bivariate_pricing_kernel <- R6::R6Class("bivariate_pricing_kernel",
  public = list(
    initialize = bivariate_pricing_kernel_constructor # Constructor, will assign all fields
    , get_returns = function() {
      home_returns <- private$home_returns %>%
        mutate(country = "home")

      foreign_returns <- private$foreign_returns %>%
        mutate(country = "foreign")

      bind_rows(home_returns, foreign_returns)
    },
    get_returns_tidy = function() {
      returns_wide <- self$get_returns()
      returns_tidy <- returns_wide %>%
        tidyr::pivot_longer(
          names_to = "asset",
          values_to = "return",
          cols = -all_of("date", "country"))

      returns_tidy
    },
    get_sdf_series = function(wide=FALSE) {
      sdf_df <- bind_rows(private$home_sdf_series, private$foreign_sdf_series)
      if (wide) {
        sdf_df <- sdf_df %>%
          tidyr::pivot_wider(
            names_from = "country",
            values_from = "sdf")
      }
      sdf_df
    },
    get_type = function() {
      private$type
    },
    get_pfolio_wts = function() {
      list(
        home = private$home_pfolio_wts,
        foreign = private$foreign_pfolio_wts
      )
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

      ## Turn return dfs into matrices
      home_return_matrix <- private$home_returns %>%
        select(-date) %>%
        as.matrix()

      foreign_return_matrix <- private$foreign_returns %>%
        select(-date) %>%
        as.matrix()

      ## Solve optimisation problem: infimum of MGF
      optimisation_problem <- solve_bivariate_entropy_problem_nloptr(
        entropy_foos = private$entropy_foos,
        home_return_matrix = home_return_matrix,
        foreign_return_matrix = foreign_return_matrix,
        solver_trace = solver_trace,
        ...
      )

      ## Collect information about solution
      # wts <- cccp::getx(optimisation_problem)
      # mgf <- cccp::getstate(optimisation_problem)["pobj"]
      wts <- optimisation_problem$par
      mgf <- optimisation_problem$value
      names(mgf) <- NULL

      ## Create and store SDF
      sdf_matrix <- private$entropy_foos$sdf_recovery(
        wts,
        home_return_matrix,
        foreign_return_matrix
      )

      home_sdf_series <- tibble::tibble(
        date = private$home_returns$date,
        country = "home",
        sdf = sdf_matrix[, 1L]
        )

      foreign_sdf_series <- tibble::tibble(
        date = private$foreign_returns$date,
        country = "foreign",
        sdf = sdf_matrix[, 2L]
        )

      private$home_sdf_series <- home_sdf_series
      private$foreign_sdf_series <- foreign_sdf_series

      ## Create and store wts
      home_wts <- head(wts, ncol(home_return_matrix) - 1L)
      foreign_wts <- tail(wts, ncol(foreign_return_matrix) - 1L)

      home_wts <- c(1.0 - sum(home_wts), home_wts)
      foreign_wts <- c(1.0 - sum(foreign_wts), foreign_wts)

      names(home_wts) <- colnames(home_return_matrix)
      names(foreign_wts) <- colnames(foreign_return_matrix)

      private$home_pfolio_wts <- home_wts
      private$foreign_pfolio_wts <- foreign_wts

      ## Store mgf value, cgf value of both returns and SDF
      sdf_powers <- private$entropy_foos$get_sdf_powers()
      private$return_mgf <- mgf
      private$return_cgf <- log(mgf)
      private$sdf_cgf <- (1.0 - sum(sdf_powers)) * private$return_cgf
      private$sdf_mgf <- exp(private$sdf_cgf)

      # Calculate and store dispersion bound
      risk_free_returns <- cbind(home_return_matrix[, 1L], foreign_return_matrix[, 1L])
      risk_free_returns <- 1.0 / risk_free_returns
      average_risk_free_returns <- log(apply(risk_free_returns, 2L, mean))

      disp_bound <- crossprod(sdf_powers, average_risk_free_returns) # without minus because calculated as mean of SDF
      disp_bound <- disp_bound - private$sdf_cgf
      disp_bound <- disp_bound / sum(sdf_powers * (1.0 - sdf_powers))

      private$sdf_dispersion_bound <- disp_bound

      ## Set fitted flag
      private$fitted <- TRUE

      invisible(self)
    },
    asset_pricing = function(new_home_returns = NULL, new_foreign_returns = NULL) {

    },
    get_sdf_summary = function() {
      sdf_powers <- private$entropy_foos$get_sdf_powers()
      sdf_cgf <- private$sdf_cgf
      sdf_dispersion_bound <- private$sdf_dispersion_bound

      return(c(sdf_powers, cgf = sdf_cgf, dispersion = sdf_dispersion_bound))
    }
  ),
  private = list(
    entropy_foos = NULL,
    home_returns = NULL,
    foreign_returns = NULL,
    type = NULL,
    home_sdf_series = NULL,
    foreign_sdf_series = NULL,
    home_pfolio_wts = NULL,
    foreign_pfolio_wts = NULL,
    return_mgf = NULL,
    return_cgf = NULL,
    sdf_mgf = NULL,
    sdf_cgf = NULL,
    sdf_dispersion_bound = NULL,
    pfolio_wts_df = NULL,
    fitted = NULL,
    objective = NULL
  )
)