

ltm_regularize_spidf <- function(spidf_obj, method = "linear", use_cloud_mask = TRUE, ...) {
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  supported_methods <- c("linear", "spline", "stine", "mean", "kalman", "locf", "nocb")
  if (!(method %in% supported_methods)) {
    stop(paste0("Interpolation method '", method, "' is not supported.\n",
                "Choose one of: ", paste(supported_methods, collapse = ", ")))
  }
  
  # Step 1: Prepare and clean the time series
  df <- spidf_obj %>%
    mutate(ts = as.Date(ti)) %>%
    arrange(ts)

  # Step 2: Create daily time sequence
  start_date <- get_range_start(spidf_obj)
  end_date   <- get_range_end(spidf_obj)
  regular_dates <- seq(from = start_date, to = end_date, by = "1 day")
  
  regular_df <- data.frame(ts = regular_dates)
  
  # Step 3: Join original data with daily dates
  df_regular <- full_join(regular_df, df, by = "ts") %>%
    arrange(ts)
  
  # Choose the correct column based on the flag
  spi_ <- if (use_cloud_mask) df_regular$masked_vals else df_regular$spi
  
  
  # Step 4: Interpolate using selected method
  interpolated_vals <- switch(method,
                              "linear" = na_interpolation(spi_, option = "linear", ...),
                              "spline" = na_interpolation(spi_, option = "spline", ...),
                              "stine"  = na_interpolation(spi_, option = "stine", ...),
                              "mean"   = na_mean(spi_, ...),
                              "kalman" = na_kalman(spi_, model = "StructTS", ...),
                              "locf"   = na_locf(spi_, ...),
                              "nocb"   = na_locf(spi_, option = "nocb", ...)
  )
  
  # Step 5: Construct the new spidf object
  out_df <- df_regular %>%
    mutate(original_spi = spi) %>% 
    mutate(spi = interpolated_vals) %>% 
    mutate(ti = ts) %>%
    select(ti, spi, original_spi, masked_vals, cloud_mask)
  
  # Copy and update metadata
  out_df <- ltm_copy_metadata(out_df, spidf_obj)
  attr(out_df, "regularized") <- TRUE
  attr(out_df, "regularize_method") <- method
  
  
  # Set class
  class(out_df) <- c("spidf", class(out_df))
  
  return(out_df)
}


ltm_apply_moving_quantile <- function(spidf_obj, quant = 0.95, win_size = 9) {
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  if (!isTRUE(is_regularized(spidf_obj))) {
    stop("Moving quantile can only be applied to a regularized spidf object.")
  }
  
  # Split the window evenly before and after
  before_days <- floor(win_size / 2)
  after_days  <- ceiling(win_size / 2)
  
  df <- spidf_obj %>%
    mutate(
      spi_mov_wind = slide_index_dbl(
        .x = spi,
        .i = ti,
        .f = ~ quantile(.x, probs = quant, na.rm = TRUE),
        .before = days(before_days),
        .after = days(after_days),
        .complete = FALSE
      )
    )
  
  # Copy and store new metadata
  df <- ltm_copy_metadata(df, spidf_obj)
  
  attr(df, "mov_window") <- TRUE
  attr(df, "mov_window_quantile") <- quant
  attr(df, "mov_window_size_days") <- win_size
  class(df) <- class(spidf_obj)
  
  return(df)
}



ltm_apply_whitaker <- function(spidf_obj, lambda = 5000, quantile_thresh = 0.35, 
                           use_weights = TRUE, min_weight = 0.01) {
  
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  # Check if series is regularized
  if (!isTRUE(is_regularized(spidf_obj))) {
    stop("Whittaker smoother requires a regularized (daily) spidf object.")
  }
  
  df <- spidf_obj
  
  # -----------------------------
  # Create weights (optional)
  # -----------------------------
  if (use_weights) {
    
    threshold <- quantile(df$spi, probs = quantile_thresh, na.rm = TRUE)
    
    df <- df %>%
      mutate(
        spi_weight = case_when(
          spi >= threshold ~ 1,
          !is.na(spi) ~ min_weight,
          TRUE ~ 0
        )
      )
    
    df$spi_smooth <- whit2(y = df$spi, w = df$spi_weight, lambda = lambda)
    
    # Smooth spi_mov_wind if it exists
    if ("spi_mov_wind" %in% names(df)) {
      
      # DO NOT USE WHEIGTS HERE!!!
      df$spi_mov_smooth <- whit2(y = df$spi_mov_wind, lambda = lambda)
      
    }
  } else {
    # No weights: standard smoothing
    df$spi_smooth <- whit2(y = df$spi, lambda = lambda)
    
    if ("spi_mov_wind" %in% names(df)) {
      df$spi_mov_smooth <- whit2(y = df$spi_mov_wind, lambda = lambda)
    }
  }
  
  # -----------------------------
  # Preserve & update metadata
  # -----------------------------

  df <- ltm_copy_metadata(df, spidf_obj)
  
  attr(df, "whit_smoothing") <- TRUE
  attr(df, "whit_lambda") <- lambda
  attr(df, "whit_quantile_threshold") <- if (use_weights) quantile_thresh else NA
  attr(df, "whit_weights_used") <- use_weights
  

  class(df) <- class(spidf_obj)
  return(df)
}



