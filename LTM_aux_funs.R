

ltm_spidf_to_zoo <- function(spidf_obj, yname){
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  ts_obj <- zoo(spidf_obj[,yname], order.by = spidf_obj$ti)
  
  return(ts_obj)
  
}

ltm_spidf_to_ts <- function(spidf_obj, yname, freq = 365){
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  spidf_obj <- ltm_remove_leap_days(spidf_obj)
  
  ts_obj <- ts(spidf_obj[,yname],
                 start = c(year(min(spidf_obj$ti)), yday(min(spidf_obj$ti))),
                 frequency = freq)
  
  return(ts_obj)
  
}

ltm_remove_leap_days <- function(spidf_obj) {
  
  stopifnot(inherits(spidf_obj, "spidf"))
  stopifnot("ti" %in% names(spidf_obj))  # Ensure the column exists
  
  out_df <- spidf_obj[!(format(spidf_obj$ti, "%m-%d") == "02-29"), ]
  
  out_df <- ltm_copy_metadata(out_df,spidf_obj)
  class(out_df) <- c("spidf", "data.frame")
  
  return(out_df)
}


ltm_get_dates <- function(spidf_obj, remove_leap=TRUE) {
  stopifnot(inherits(spidf_obj, "spidf"))
  stopifnot("ti" %in% names(spidf_obj))  # Ensure the column exists
  
  if(remove_leap){
    out_ti <- as.Date(spidf_obj[!(format(spidf_obj$ti, "%m-%d") == "02-29"), "ti"])
  }else{
    out_ti <- as.Date(spidf_obj[, "ti"])
  }
  return(out_ti)
}


ltm_get_duplicate_date_indices <- function(spidf_obj) {
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  stopifnot("ti" %in% names(spidf_obj))
  
  # Ensure ti is Date class
  ti_dates <- as.Date(spidf_obj$ti)
  
  # Get duplicated indices: first and subsequent occurrences
  duplicated_indices <- which(duplicated(ti_dates) | duplicated(ti_dates, fromLast = TRUE))
  
  return(duplicated_indices)
}


ltm_remove_duplicates <- function(spidf_obj) {
  
  stopifnot(inherits(spidf_obj, "spidf"))
  
  stopifnot(all(c("ti", "spi") %in% names(spidf_obj)))
  
  out_df <- spidf_obj %>%
    mutate(ti = as.Date(ti)) %>%
    group_by(ti) %>%
    slice_max(order_by = spi, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(ti) %>% 
    as.data.frame()
  
  out_df <- ltm_copy_metadata(out_df,spidf_obj)
  class(out_df) <- c("spidf", "data.frame")
  
  return(out_df)
}

ltm_copy_metadata <- function(x,from){
  
  #stopifnot(inherits(x, "spidf"))
  #stopifnot(inherits(from, "spidf"))
  
  # Copy metadata
  attr(x, "lat") <- get_latitude(from)
  attr(x, "lon") <- get_longitude(from)
  attr(x, "start_date") <- get_start_date(from)
  attr(x, "end_date") <- get_end_date(from)
  attr(x, "range_start") <- get_range_start(from)
  attr(x, "range_end") <- get_range_end(from)
  
  attr(dt, "tree_id") <- get_tree_id(from)
  
  attr(x, "spi") <- get_spi(from)
  attr(x, "proc_level") <- get_proc_level(from)
  attr(x, "crs_code") <- get_crs_code(from)
  
  attr(x, "regularized") <- is_regularized(from)
  attr(x, "regularize_method") <- get_regularize_method(from)
  
  attr(x, "mov_window") <- has_mov_window(from)
  attr(x, "mov_window_quantile") <- get_mov_window_quantile(from)
  attr(x, "mov_window_size_days") <- get_mov_window_size_days(from)
  
  attr(x, "whit_smoothing") <- is_whit_smoothed(from)
  attr(x, "whit_lambda") <- get_whit_lambda(from)
  attr(x, "whit_quantile_threshold") <- get_whit_quantile_threshold(from)
  attr(x, "whit_weights_used") <- get_whit_weights_used(from)
  
  return(x)
}


ltm_days_between <- function(start_date, end_date) {
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  
  as.integer(difftime(end_date, start_date, units = "days"))
}

ltm_convert_breakdate_to_date <- function(breakDate) {
  year <- floor(breakDate)
  fractional_year <- breakDate - year
  day_of_year <- round(fractional_year * 365)
  
  # Generate the date
  formal_date <- as.Date(day_of_year - 1, origin = paste0(year, "-01-01"))
  
  return(formal_date)
}


ltm_copy_ts <- function(ts_to_copy, values){
  
  stopifnot(is.ts(ts_to_copy))
  
  if (length(values) != length(ts_to_copy)){
    stop("Length of values vector must match length of ts_to_copy")
  }
  
  new_ts <- ts(values,
               start = start(ts_to_copy),
               end = end(ts_to_copy),
               frequency = frequency(ts_to_copy))
  
  return(new_ts)
}



ltm_summarize_break_df <- function(break_df) {
  
  if (nrow(break_df) == 0 || !any(break_df$has_valid_breaks)) {
    return("No valid breaks were detected.")
  }
  
  df_valid <- subset(break_df, has_valid_breaks == TRUE & !is.na(break_date))
  
  summary_list <- df_valid |>
    dplyr::group_by(method, data_type) |>
    dplyr::summarise(
      num_breaks = dplyr::n(),
      break_date = min(break_date),
      #last_break = max(break_date),
      avg_magnitude = round(mean(break_magn, na.rm = TRUE), 2),
      .groups = "drop" 
    )
  
  # Format into a human-readable summary string
  summary_text <- paste0(
    "Break detection summary:\n",
    "-------------------------\n",
    paste(
      apply(summary_list, 1, function(row) {
        sprintf(". %s on %s: %d break(s), on %s (change = %.2f%%)",
                toupper(row["method"]),
                row["data_type"],
                as.integer(row["num_breaks"]),
                format(as.Date(row["break_date"]), "%Y-%m-%d"),
                #format(as.Date(row["last_break"]), "%Y-%m-%d"),
                as.numeric(row["avg_magnitude"]))
      }),
      collapse = "\n"
    )
  )
  
  return(summary_text)
}

Percentile90 <- function(x) quantile(x,probs=0.9,na.rm=TRUE)


ltm_loc_file <- file.path("ltm_cache", "current_location.txt")

ltm_read_location <- function() {
  # Check if the file exists
  if (!file.exists(ltm_loc_file)) {
    return(NULL)
  }
  # read the single line
  loc_str <- readLines(ltm_loc_file, warn = FALSE)
  if (length(loc_str) > 0) {
    return(loc_str[1])
  } else {
    return(NULL)
  }
}

ltm_store_location <- function(lat, lon) {
  # Ensure the directory exists
  if (!dir.exists("ltm_cache")) dir.create("ltm_cache", recursive = TRUE)
  
  # Round or format to avoid floating point artifacts
  loc_str <- paste(round(lat, 6), round(lon, 6), sep = ",")
  
  writeLines(loc_str, con = ltm_loc_file)
}

