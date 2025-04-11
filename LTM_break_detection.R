


# Print method for ts_breaks_run class
print.ts_breaks_run <- function(x, ...) {
  cat("Break detection run summary:\n")
  cat("-----------------------------\n")
  cat(sprintf("%-22s : %s\n", "Method", x$method))
  cat(sprintf("%-22s : %s\n", "Data type", x$data_type))
  cat(sprintf("%-22s : %s\n", "Season adjusted", ifelse(x$season_adj, "Yes", "No")))
  cat(sprintf("%-22s : %s\n", "Breaks detected", ifelse(x$has_breaks, "Yes", "No")))
  cat(sprintf("%-22s : %s\n", "Valid breaks detected", ifelse(x$has_valid_breaks, "Yes", "No")))
  
  cat("-----------------------------\n")
  if (x$has_breaks) {
    cat(sprintf("%-22s : %d\n", "Number of breaks", length(x$breaks_indices)))
    cat(sprintf("%-22s : %s\n", "Break indices", paste(x$breaks_indices, collapse = ", ")))
    cat(sprintf("%-22s : %s\n", "Break dates", paste(x$breaks_dates, collapse = ", ")))
    cat(sprintf("%-22s : %s\n", "Change magnitudes", paste(round(x$break_magn, 4), collapse = ", ")))
  } else {
    cat("No breaks were identified in this run.\n")
  }
  
  cat("-----------------------------\n")
  if (x$has_valid_breaks) {
    cat(sprintf("%-22s : %d\n", "Number of valid breaks", length(x$breaks_indices)))
    cat(sprintf("%-22s : %s\n", "Valid break indices", paste(x$breaks_indices, collapse = ", ")))
    cat(sprintf("%-22s : %s\n", "Valid break dates", paste(x$breaks_dates, collapse = ", ")))
    cat(sprintf("%-22s : %s\n", "Valid change magnitudes", paste(round(x$break_magn, 4), collapse = ", ")))
  } else {
    cat("No valid breaks were identified in this run.\n")
  }
  
  cat("\nFunction call:\n")
  print(x$call)
}

# Wrapper function for cpm Exponential detector

ltm_cpm_detect_breaks <- function(spidf,
                                  season_adj = TRUE, 
                                  ts_name = "spi",
                                  s_window = 30,
                                  cpm_method = "Exponential",
                                  ARL0 = 500, 
                                  thresh_date,
                                  thresh_change = -10,
                                  tresh_int       = NULL,
                                  thresh_fun = median,
                                  ...){
  
  stopifnot(inherits(spidf, "spidf"))
  
  #valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  
  if (!(ts_name %in% VALID_DATA_TYPES)) {
    stop("Invalid ts_name! Must be one of: ", paste(VALID_DATA_TYPES, collapse = ", "))
  }
  
  
  yts <- ltm_spidf_to_ts(spidf, ts_name)
  dts <- ltm_get_dates(spidf, remove_leap=TRUE)
  
  if(season_adj){
    decomp <- stl(yts, s.window = s_window)
    
    ysa <- seasadj(decomp)
  }else{
    ysa <- yts
  }
  
  startup <- ltm_days_between(get_range_start(spidf),
                          thresh_date)
  
  cp <- detectChangePoint(ysa,
                          cpmType = cpm_method,
                          ARL0    = ARL0, 
                          startup = startup,...)
  
  has_breaks <- cp$changeDetected

  
  if(has_breaks){
    
    brk <- cp$changePoint
    break_date <- dts[brk]
    
    if(is.null(tresh_int)){
      pre_break  <- thresh_fun(ysa[1:(brk-1)])
      post_break <- thresh_fun(ysa[(brk+1):length(ysa)])
      
    }else{
      pre_start <- max(1, brk - tresh_int)
      pre_end   <- max(1, brk - 1)
      post_start <- min(length(ysa), brk + 1)
      post_end   <- min(length(ysa), brk + tresh_int)
      
      pre_break  <- thresh_fun(ysa[pre_start:pre_end])
      post_break <- thresh_fun(ysa[post_start:post_end])
    }
    
    break_magn <- ((post_break - pre_break) / pre_break)*100
    
    if(
      (break_date >= thresh_date) && 
      (break_magn <= thresh_change) #&&
    ){
      
      out <-       
        list(method           = "cpm",
             data_type        = ts_name,
             has_breaks       = TRUE,
             has_valid_breaks = TRUE,
             break_magn       = break_magn,
             breaks_indices   = brk,
             breaks_dates     = break_date,
             output_object    = cp,
             #time_series_data = yts,
             season_adj       = season_adj,
             call             = match.call())
      class(out) <- "ts_breaks_run"
      return(out)
    }else{
      
      out <- list(method="cpm",
                  data_type        = ts_name,
                  has_breaks       = TRUE,
                  has_valid_breaks = FALSE,
                  break_magn       = break_magn,
                  breaks_indices   = brk,
                  breaks_dates     = break_date,
                  output_object    = cp,
                  #time_series_data = yts,
                  season_adj       = season_adj,
                  call             = match.call())
      
      class(out) <- "ts_breaks_run"
      return(out)
    }
  }
  else{
    
    out <- list(method           = "cpm",
                data_type        = ts_name,
                has_breaks       = FALSE,
                has_valid_breaks = FALSE,
                break_magn       = NA,
                breaks_indices   = NA,
                breaks_dates     = NA,
                output_object    = cp,
                #time_series_data = yts,
                season_adj       = season_adj,
                call             = match.call())
    
    class(out) <- "ts_breaks_run"
    return(out)
  }
}


# Wrapper function for running the bfast01 algo
ltm_bfast01_detect_breaks <- function(
    spidf,
    ts_name = "spi",
    formula = response ~ harmon + trend,
    s_window   = 30,
    test = "OLS-MOSUM",
    level = 0.05,
    aggregate = all,
    trim = NULL,
    bandwidth = 0.15,
    functional = "max",
    order = 3, 
    thresh_date,
    thresh_change = -10,
    tresh_int       = NULL,
    thresh_fun = median,
    ...){
  
  
  stopifnot(inherits(spidf, "spidf"))
  
  #valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  
  if (!(ts_name %in% VALID_DATA_TYPES)) {
    stop("Invalid ts_name! Must be one of: ", paste(VALID_DATA_TYPES, collapse = ", "))
  }
  
  
  yts <- ltm_spidf_to_ts(spidf, ts_name)
  dts <- ltm_get_dates(spidf, remove_leap=TRUE)
  
  bf01 <- bfast01(
    data = yts,
    formula = formula,
    test = test,
    level = level,
    aggregate = aggregate,
    trim = trim,
    bandwidth = bandwidth,
    functional = functional,
    order = order, ...)
  
  
  #break_date <- breakdates(bf01)
  brk <- bf01$breakpoints
  
  
  if(is.na(brk)){
    
    out <- list(method="bfast01",
                data_type = ts_name,
                has_breaks = FALSE,
                has_valid_breaks = FALSE,
                break_magn = NA,
                breaks_indices = NA,
                breaks_dates = NA,
                output_object = bf01,
                #time_series_data = yts,
                season_adj = FALSE,
                call = match.call())
    class(out) <- "ts_breaks_run"
    
    return(out)
    
  }else{
    
    
    brk <- bf01$breakpoints
    break_date <- dts[brk]
    fitted_ts <- ltm_copy_ts(yts, values=fitted(bf01))
    decomp <- stl(fitted_ts, s.window = s_window)
    fitted_ts <- seasadj(decomp)
    
    if(is.null(tresh_int)){
      pre_break  <- thresh_fun(fitted_ts[1:(brk-1)])
      post_break <- thresh_fun(fitted_ts[(brk+1):length(fitted_ts)])
      
    }else{
      pre_start <- max(1, brk - tresh_int)
      pre_end   <- max(1, brk - 1)
      post_start <- min(length(fitted_ts), brk + 1)
      post_end   <- min(length(fitted_ts), brk + tresh_int)
      
      pre_break  <- thresh_fun(fitted_ts[pre_start:pre_end])
      post_break <- thresh_fun(fitted_ts[post_start:post_end])
    }
    
    break_magn <- ((post_break - pre_break) / pre_break)*100
    
    if(
      (break_date >= thresh_date) && 
      (break_magn <= thresh_change) #&&
    ){
      
      out <- list(method="bfast01",
                  data_type = ts_name,
                  has_breaks = TRUE,
                  has_valid_breaks = TRUE,
                  break_magn = break_magn,
                  breaks_indices = brk,
                  breaks_dates = break_date,
                  output_object = bf01,
                  #time_series_data = yts,
                  season_adj = FALSE,
                  call = match.call())
      class(out) <- "ts_breaks_run"
      return(out)
      
    }else{
      out <- list(method="bfast01",
                  data_type = ts_name,
                  has_breaks = TRUE,
                  has_valid_breaks = FALSE,
                  break_magn = break_magn,
                  breaks_indices = brk,
                  breaks_dates = break_date,
                  output_object = bf01,
                  #time_series_data = yts,
                  season_adj = FALSE,
                  call = match.call())
      class(out) <- "ts_breaks_run"
      return(out)
    }
  }
}



ltm_ed_detect_breaks <- function(spidf,
                                 ts_name = "spi",
                                 sig_lvl   = 0.05, 
                                 R         = 1000, 
                                 k         = 1, 
                                 min_size  = 30, 
                                 alpha     = 1, 
                                 season_adj = TRUE,
                                 s_window   = 30,
                                 thresh_change  = -10, 
                                 thresh_date, 
                                 tresh_int       = NULL,
                                 thresh_fun = median,
                                 ...){
  
  stopifnot(inherits(spidf, "spidf"))
  
  #valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  
  if (!(ts_name %in% VALID_DATA_TYPES)) {
    stop("Invalid ts_name! Must be one of: ", paste(VALID_DATA_TYPES, collapse = ", "))
  }
  
  yts <- ltm_spidf_to_ts(spidf, ts_name)
  dts <- ltm_get_dates(spidf, remove_leap=TRUE)
  
  
  if(season_adj){
    decomp <- stl(yts, s.window = s_window)
    
    ysa <- seasadj(decomp)
  }else{
    ysa <- yts
  }
  
  mts <- matrix(ysa,nrow=length(ysa),1)
  
  ed <- e.divisive(mts, sig.lvl=sig_lvl, R=R, k=k, 
                   min.size=min_size, alpha=alpha)
  
  has_breaks <- ifelse(length(unique(ed$cluster))>1,TRUE,FALSE)

  
  if(has_breaks){
    
    brk <- ed$estimates[2]
    break_date <- dts[brk]
    
    if(is.null(tresh_int)){
      pre_break  <- thresh_fun(ysa[1:(brk-1)])
      post_break <- thresh_fun(ysa[(brk+1):length(ysa)])
      
    }else{
      pre_start <- max(1, brk - tresh_int)
      pre_end   <- max(1, brk - 1)
      post_start <- min(length(ysa), brk + 1)
      post_end   <- min(length(ysa), brk + tresh_int)
      
      pre_break  <- thresh_fun(ysa[pre_start:pre_end])
      post_break <- thresh_fun(ysa[post_start:post_end])
    }
    
    break_magn <- ((post_break - pre_break) / pre_break)*100
    
    if(
      (break_date >= thresh_date) && 
      (break_magn <= thresh_change) #&&
    ){
      
      out <-       
        list(method="ed",
             data_type = ts_name,
             has_breaks = TRUE,
             has_valid_breaks = TRUE,
             break_magn = break_magn,
             breaks_indices = brk,
             breaks_dates = break_date,
             output_object = ed,
             #time_series_data = yts,
             season_adj = season_adj,
             call = match.call())
      class(out) <- "ts_breaks_run"
      return(out)
    }else{
      
      out <- list(method="ed",
                  data_type = ts_name,
                  has_breaks = TRUE,
                  has_valid_breaks = FALSE,
                  break_magn = break_magn,
                  breaks_indices = brk,
                  breaks_dates = break_date,
                  output_object = ed,
                  #time_series_data = yts,
                  season_adj = season_adj,
                  call = match.call())
      
      class(out) <- "ts_breaks_run"
      return(out)
      
    }
  }
  else{
    
    out <- list(method="ed",
                data_type = ts_name,
                has_breaks = FALSE,
                has_valid_breaks = FALSE,
                break_magn = NA,
                breaks_indices = NA,
                breaks_dates = NA,
                output_object = ed,
                #time_series_data = yts,
                season_adj = season_adj,
                call = match.call())
    
    class(out) <- "ts_breaks_run"
    return(out)
  }
}


mcpSummary <- function (object, width = 0.95, digits = 2, 
                        prior = FALSE, ...) {
  fit = object
  mcp:::assert_mcpfit(fit)
  mcp:::assert_numeric(width, lower = 0, upper = 1)
  mcp:::assert_integer(digits, lower = 0)
  mcp:::assert_logical(prior)
  mcp:::assert_ellipsis(...)
  samples = mcp:::mcmclist_samples(fit, prior = prior, error = FALSE)
  
  if (!is.null(samples)) {
    
    result = mcp:::get_summary(fit, width, varying = FALSE, prior = prior)
    
    return(result)
  }else{
    return(NULL)
  }
}



diagnose_uncertainty <- function(df,
                                 w_width = 1,
                                 w_rhat  = 1,
                                 w_neff  = 1,
                                 cutoff_low = 0.01,
                                 cutoff_med = 0.05) {
  # Check columns exist
  req_cols <- c("mean", "lower", "upper", "Rhat", "n.eff")
  if (!all(req_cols %in% names(df))) {
    stop(paste("Dataframe must contain columns:",
               paste(req_cols, collapse = ", ")))
  }
  
  # Calculate relative credible interval width
  width_ratio <- abs(df$upper - df$lower) / abs(df$mean)
  
  # Calculate how far Rhat is from 1
  rhat_dev <- abs(df$Rhat - 1)
  
  # Inverse of effective sample size
  inv_neff <- 1 / df$n.eff
  
  # Combine them into a single uncertainty metric for each parameter
  param_uncert <- w_width * width_ratio + w_rhat * rhat_dev + w_neff * inv_neff
  
  # Average across all parameters to produce final uncertainty score
  uncertainty_score <- mean(param_uncert, na.rm = TRUE)
  
  # Qualitative classification
  classification <- dplyr::case_when(
    uncertainty_score < cutoff_low ~ "Low",
    uncertainty_score < cutoff_med ~ "Medium",
    TRUE                           ~ "High"
  )
  
  # Return a list with both numeric and qualitative summaries
  list(
    uncertainty_score = uncertainty_score,
    classification = classification
  )
}


ltm_mcp_detect_breaks <- function(spidf,
                                ts_name = "spi",
                                season_adj = TRUE,
                                s_window   = 30,
                                thresh_change  = -10, 
                                thresh_date,
                                
                                sample = "both",  
                                n_chains = 3,
                                n_cores = 3,
                                n_adapt = 500,
                                n_iter = 1000,
                                downsample=NULL,
                                ...){
  
  stopifnot(inherits(spidf, "spidf"))
  
  #valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  
  if (!(ts_name %in% VALID_DATA_TYPES)) {
    stop("Invalid ts_name! Must be one of: ", paste(VALID_DATA_TYPES, collapse = ", "))
  }
  
  # Prepare the time series data to fit the model
  yts <- ltm_spidf_to_ts(spidf, ts_name)
  dts <- ltm_get_dates(spidf, remove_leap=TRUE)
  
  if(!is.null(downsample)){
    ds_ids <- seq(1,length(yts),by=downsample)
    yts <- ltm_spidf_to_ts(spidf[ds_ids,], ts_name,freq = downsample)
    dts <- ltm_get_dates(spidf[ds_ids,], remove_leap=TRUE)
  }else{
    # Prepare the time series data to fit the model
    yts <- ltm_spidf_to_ts(spidf, ts_name)
    dts <- ltm_get_dates(spidf, remove_leap=TRUE)
    
  }
  
  if(season_adj){
    decomp <- stl(yts, s.window = s_window)
    
    ysa <- seasadj(decomp)
  }else{
    ysa <- yts
  }
  
  # Setup model data
  model_df <- data.frame(y = yts, ysa = ysa, 
                         jd = yday(dts), di=1:length(ysa))
  
  # Setup model
  model <- list(ysa~1, 
                ~1)
  
  # Define priors
  qts <- round(quantile(ysa, probs=c(0.75,0.25)), 3)
  std <- round(sd(ysa), 3)
  
  my_priors <- list(
    # First intercept: prior set at the the 75% percentile 
    # Truncate to [0,1] to reflect NDVI expectable domain.
    int_1 = paste0("dnorm(",qts[1],", ",std,") T(0, 1)"),
    
    # Second intercept: prior set at the the 25% percentile
    # but forced strictly below int_1 via T( , int_1)
    int_2 = paste0("dnorm(",qts[2],", ",std,") T(0, int_1)"),
    
    # Single change point with equal p at any point in time (uniform prior)
    cp_1  = paste0("dunif(1, ",length(ysa),")")
  )
  
  # 4) Fit the model
  fit_mcp <- mcp(
    model,
    data   = model_df,
    par_x  = "di",      # The variable representing x-axis
    prior  = my_priors,
    adapt  = n_adapt,  # Optional: controls tuning of step size
    sample = sample,  # This ensures both posterior and prior sampling
    chains = n_chains,
    cores  = n_cores,
    iter   = n_iter,
    ...
  )
  
  ## Get model parameters
  pars <- as.data.frame(mcpSummary(fit_mcp))
  brk <- as.integer(round(pars$mean[1]))
  break_date <- dts[brk]

  # Calculate the break magnitude by comparing the two intercepts
  pre_break  <- pars[2,2]
  post_break <- pars[3,2]
  break_magn <- ((post_break - pre_break) / pre_break)*100
  
  if(
    (break_date >= thresh_date) && 
    (break_magn <= thresh_change) #&&
  ){
    
    out <-       
      list(method="mcp",
           data_type = ts_name,
           has_breaks = TRUE,
           has_valid_breaks = TRUE,
           break_magn = break_magn,
           breaks_indices = brk,
           breaks_dates = break_date,
           output_object = fit_mcp,
           #time_series_data = yts,
           season_adj = season_adj,
           uncertainty_diag = diagnose_uncertainty(pars),
           pars = pars,
           call = match.call())
    class(out) <- "ts_breaks_run"
    return(out)
  }else{
    
    out <- list(method="mcp",
                data_type = ts_name,
                has_breaks = TRUE,
                has_valid_breaks = FALSE,
                break_magn = break_magn,
                breaks_indices = brk,
                breaks_dates = break_date,
                output_object = fit_mcp,
                #time_series_data = yts,
                season_adj = season_adj,
                uncertainty_diag = diagnose_uncertainty(pars),
                pars = pars,
                call = match.call())
    
    class(out) <- "ts_breaks_run"
    return(out)
    
  }
}


ltm_strucchange_detect_breaks <- function(
    spidf,
    ts_name = "spi",
    season_adj = TRUE,
    s_window   = 30,
    h          = 0.15,
    breaks     = 1,
    thresh_date,
    thresh_change = -10,
    tresh_int     = NULL,
    thresh_fun    = median,
    ...
) {
  # 1) Basic sanity checks
  
  stopifnot(inherits(spidf, "spidf"))
  
  #valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  
  if (!(ts_name %in% VALID_DATA_TYPES)) {
    stop("Invalid ts_name! Must be one of: ", paste(VALID_DATA_TYPES, collapse = ", "))
  }
  
  # 2) Convert spidf to a time-series object and get the date vector
  yts <- ltm_spidf_to_ts(spidf, ts_name)
  dts <- ltm_get_dates(spidf, remove_leap = TRUE)
  
  # 3) (Optional) Seasonal adjustment via STL
  if (season_adj) {
    decomp <- stl(yts, s.window = s_window)
    ysa <- seasadj(decomp)
  } else {
    ysa <- yts
  }
  
  # 4) Use breakpoints from 'strucchange' to detect structural breaks
  #    We'll fit a simple intercept-only model y ~ 1.
  bp <- strucchangeRcpp::breakpoints(ysa ~ 1,
                                 h = h,
                                 breaks = breaks,
                                 ...)
  
  # 5) Determine if a break was found
  #    If the best model is 'breaks = 0', we typically get NA in bp$breakpoints.
  if (all(is.na(bp$breakpoints))) {
    # No break found
    out <- list(
      method           = "stc",
      data_type        = ts_name,
      has_breaks       = FALSE,
      has_valid_breaks = FALSE,
      break_magn       = NA,
      breaks_indices   = NA,
      breaks_dates     = NA,
      output_object    = bp,
      season_adj       = season_adj,
      call             = match.call()
    )
    class(out) <- "ts_breaks_run"
    return(out)
  }
  
  # 6) For simplicity, we only consider the first break if multiple found
  #    If you allow multiple breaks, you can adapt the logic below.
  brk <- bp$breakpoints[1]
  
  # 7) Identify the break date
  break_date <- dts[brk]
  
  # 8) Evaluate pre/post intervals and compute the threshold function
  if (is.null(tresh_int)) {
    pre_break  <- thresh_fun(ysa[1:(brk - 1)])
    post_break <- thresh_fun(ysa[(brk + 1):length(ysa)])
  } else {
    pre_start  <- max(1, brk - tresh_int)
    pre_end    <- max(1, brk - 1)
    post_start <- min(length(ysa), brk + 1)
    post_end   <- min(length(ysa), brk + tresh_int)
    
    pre_break  <- thresh_fun(ysa[pre_start:pre_end])
    post_break <- thresh_fun(ysa[post_start:post_end])
  }
  
  # 9) Compute break magnitude
  break_magn <- ((post_break - pre_break) / pre_break) * 100
  
  # 10) Decide if the break is "valid" based on date + magnitude thresholds
  if ((break_date >= thresh_date) && (break_magn <= thresh_change)) {
    out <- list(
      method           = "stc",
      data_type        = ts_name,
      has_breaks       = TRUE,
      has_valid_breaks = TRUE,
      break_magn       = break_magn,
      breaks_indices   = brk,
      breaks_dates     = break_date,
      output_object    = bp,
      season_adj       = season_adj,
      call             = match.call()
    )
    class(out) <- "ts_breaks_run"
    return(out)
  } else {
    out <- list(
      method           = "stc",
      data_type        = ts_name,
      has_breaks       = TRUE,
      has_valid_breaks = FALSE,
      break_magn       = break_magn,
      breaks_indices   = brk,
      breaks_dates     = break_date,
      output_object    = bp,
      season_adj       = season_adj,
      call             = match.call()
    )
    class(out) <- "ts_breaks_run"
    return(out)
  }
}

ltm_wbs_detect_breaks <- function(spidf,
                                  ts_name = "spi",
                                  season_adj = TRUE,
                                  s_window = 30,
                                  thresh_date,
                                  thresh_change = -10,
                                  tresh_int = NULL,
                                  thresh_fun = median,
                                  num_intervals = 1000,
                                  ...)
{
  stopifnot(inherits(spidf, "spidf"))
  valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  if (!(ts_name %in% valid_data_types)) {
    stop("Invalid ts_name! Must be one of: ", paste(valid_data_types, collapse = ", "))
  }
  
  # Prepare time series
  yts <- ltm_spidf_to_ts(spidf, ts_name)
  dts <- ltm_get_dates(spidf, remove_leap = TRUE)
  
  # (Optional) Seasonal adjustment
  if (season_adj) {
    decomp <- stl(yts, s.window = s_window)
    ysa <- seasadj(decomp)
  } else {
    ysa <- yts
  }
  
  # Fit Wild Binary Segmentation
  wbs_fit <- wbs::wbs(ysa, numIntervals = num_intervals)
  cp_fit <- wbs::changepoints(wbs_fit, penalty = "ssic.penalty", ...)  # or "sic", "ssic.penalty", etc.
  
  # If no breakpoints
  if (length(cp_fit$cpt.th) == 0) {
    out <- list(
      method           = "wbs",
      data_type        = ts_name,
      has_breaks       = FALSE,
      has_valid_breaks = FALSE,
      break_magn       = NA,
      breaks_indices   = NA,
      breaks_dates     = NA,
      output_object    = cp_fit,
      season_adj       = season_adj,
      call             = match.call()
    )
    class(out) <- "ts_breaks_run"
    return(out)
  }
  
  # cp_fit$cpt.th might be a list of solutions or a numeric vector
  if (is.list(cp_fit$cpt.th)) {
    # Extract the first numeric solution
    brk_vec <- cp_fit$cpt.th[[1]]
  } else {
    # Already numeric
    brk_vec <- cp_fit$cpt.th
  }
  
  # If there's more than one break, we pick the first for single-break logic
  if (length(brk_vec) == 0) {
    # If there's truly no breaks in the first solution
    out <- list(
      method           = "wbs",
      data_type        = ts_name,
      has_breaks       = FALSE,
      has_valid_breaks = FALSE,
      break_magn       = NA,
      breaks_indices   = NA,
      breaks_dates     = NA,
      output_object    = cp_fit,
      season_adj       = season_adj,
      call             = match.call()
    )
    class(out) <- "ts_breaks_run"
    return(out)
  }
  
  # Take the first numeric index
  brk <- brk_vec[1]
  break_date <- dts[brk]
  
  # Evaluate pre- and post-break windows
  if (is.null(tresh_int)) {
    pre_break  <- thresh_fun(ysa[1:(brk - 1)])
    post_break <- thresh_fun(ysa[(brk + 1):length(ysa)])
  } else {
    pre_start  <- max(1, brk - tresh_int)
    pre_end    <- max(1, brk - 1)
    post_start <- min(length(ysa), brk + 1)
    post_end   <- min(length(ysa), brk + tresh_int)
    
    pre_break  <- thresh_fun(ysa[pre_start:pre_end])
    post_break <- thresh_fun(ysa[post_start:post_end])
  }
  
  break_magn <- ((post_break - pre_break) / pre_break) * 100
  valid_break <- (break_date >= thresh_date) && (break_magn <= thresh_change)
  
  out <- list(
    method           = "wbs",
    data_type        = ts_name,
    has_breaks       = TRUE,
    has_valid_breaks = valid_break,
    break_magn       = break_magn,
    breaks_indices   = brk,
    breaks_dates     = break_date,
    output_object    = cp_fit,
    season_adj       = season_adj,
    call             = match.call()
  )
  class(out) <- "ts_breaks_run"
  return(out)
}



