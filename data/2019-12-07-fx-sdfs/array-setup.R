library(dplyr)

job_array_df <- expand.grid(sdf_class = "window_cv_pricing_kernel", sdf_param = 3:12)

# job_array_df <- bind_rows(job_array_df, expand.grid(sdf_class = "window_chi2_cv_pricing_kernel", sdf_param = 3:12))

job_array_df <- bind_rows(job_array_df, expand.grid(sdf_class = "window_fs_pr_cv_pricing_kernel", sdf_param = 3:12))


readr::write_delim(job_array_df, "array-args", col_names = FALSE, quote_escape = FALSE)
