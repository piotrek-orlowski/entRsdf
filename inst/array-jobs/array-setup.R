library(dplyr)

job_array_df <- expand.grid(sdf_class = "window_cv_pricing_kernel", sdf_param = 3:10)

job_array_df <- bind_rows(job_array_df, expand.grid(sdf_class = "window_lev_pricing_kernel", sdf_param = seq(6,60,by=2)))

readr::write_delim(job_array_df, "array-args", col_names = FALSE, quote_escape = FALSE)
