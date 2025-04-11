

ltm_ts_breaks <- function(spidf_ts) {

  
  if (!inherits(spidf_ts, "spidf")) {
    stop("Error: 'spidf_ts' must be of class 'spidf'.")
  }
    
  structure(
    list(
      spidf_ts = spidf_ts,
      algorithms = list()
    ),
    class = "ts_breaks"
  )
}

ltm_add_runs <- function(ts_breaks_obj, ...) {
  
  if (!inherits(ts_breaks_obj, "ts_breaks")) {
    stop("Error: 'ts_breaks_obj' must be of class 'ts_breaks'.")
  }
  
  run_objects <- list(...)
  
  #valid_data_types <- c("spi", "spi_mov_wind", "spi_smooth", "spi_mov_smooth")
  
  for (run_object in run_objects) {
    if (!inherits(run_object, "ts_breaks_run")) {
      stop("Error: Each 'run_object' must be of class 'ts_breaks_run'.")
    }
    
    algorithm_name <- run_object$method
    if (is.null(algorithm_name) || !nzchar(algorithm_name)) {
      stop("Error: Each 'run_object' must have a non-empty 'method' field.")
    }
    
    if (!algorithm_name %in% names(ts_breaks_obj$algorithms)) {
      ts_breaks_obj$algorithms[[algorithm_name]] <- list()
    }
    
    if (!(run_object$data_type %in% VALID_DATA_TYPES)) {
      stop("Error: 'data_type' must be one of: ", paste(VALID_DATA_TYPES, collapse = ", "))
    }
    
    new_run <- list(
      method           = run_object$method,
      data_type        = run_object$data_type,
      has_breaks       = run_object$has_breaks,
      has_valid_breaks = run_object$has_valid_breaks,
      break_magn       = if (!is.null(run_object$break_magn)) run_object$break_magn else NA,
      breaks_indices   = run_object$breaks_indices,
      breaks_dates     = as.Date(run_object$breaks_dates),
      output_object    = run_object$output_object,
      spidf_ts         = ts_breaks_obj$spidf_ts,
      season_adj       = if (!is.null(run_object$season_adj)) run_object$season_adj else FALSE,
      call             = run_object$call
    )
    
    run_id <- paste0("run-", sprintf("%02d", length(ts_breaks_obj$algorithms[[algorithm_name]]) + 1))
    ts_breaks_obj$algorithms[[algorithm_name]][[run_id]] <- new_run
  }
  
  return(ts_breaks_obj)
}


ltm_get_algorithms <- function(ts_breaks_obj) {
  stopifnot(inherits(ts_breaks_obj, "ts_breaks"))
  names(ts_breaks_obj$algorithms)
}

ltm_get_runs <- function(ts_breaks_obj, algorithm_name) {
  stopifnot(inherits(ts_breaks_obj, "ts_breaks"))
  if (!algorithm_name %in% names(ts_breaks_obj$algorithms)) {
    stop(paste("Algorithm/method", algorithm_name, "does not exist in this object."))
  }
  names(ts_breaks_obj$algorithms[[algorithm_name]])
}

ltm_get_run_details <- function(ts_breaks_obj, algorithm_name, run_id) {
  stopifnot(inherits(ts_breaks_obj, "ts_breaks"))
  if (!algorithm_name %in% names(ts_breaks_obj$algorithms)) {
    stop(paste("Algorithm/method", algorithm_name, "does not exist in this object."))
  }
  if (!run_id %in% names(ts_breaks_obj$algorithms[[algorithm_name]])) {
    stop(paste("Run ID", run_id, "does not exist for algorithm", algorithm_name))
  }
  ts_breaks_obj$algorithms[[algorithm_name]][[run_id]]
}

print.ts_breaks <- function(x, ...) {
  cat("Time Series Break Detection Object (ts_breaks):\n\n")
  cat("Total algorithms used:", length(x$algorithms), "\n")
  
  if (length(x$algorithms) > 0) {
    cat("Algorithms and runs:\n")
    
    for (alg in names(x$algorithms)) {
      cat(" - ", alg, " (", length(x$algorithms[[alg]]), " run(s))\n",sep = "")
      
      for (run in names(x$algorithms[[alg]])) {
        run_obj <- x$algorithms[[alg]][[run]]
        
        cat("    .", run,
            "[method:", run_obj$method,
            "| data_type:", run_obj$data_type,
            "| breaks:", ifelse(run_obj$has_breaks, "Yes", "No"),
            "| valid_breaks:", ifelse(isTRUE(run_obj$has_valid_breaks), "Yes", "No"),
            "]\n")
        
        if (isTRUE(run_obj$has_valid_breaks)) {
          # Show first valid break index and date (or all if desired)
          cat("       > Valid break index:", paste(run_obj$breaks_indices, collapse = ", "), "\n")
          cat("       > Valid break date :", paste(run_obj$breaks_dates, collapse = ", "), "\n")
        }
      }
    }
  } else {
    cat(" No algorithms added yet.\n")
  }
  
  invisible(x)
}


as.data.frame.ts_breaks <- function(ts_breaks_obj, ...) {

  if (!inherits(ts_breaks_obj, "ts_breaks")) {
    stop("Error: 'ts_breaks_obj' must be of class 'ts_breaks'.")
  }
    
  # Initialize empty data frame
  breaks_df <- data.frame(
    algorithm        = character(0),
    run_id           = character(0),
    method           = character(0),
    data_type        = character(0),
    has_breaks       = logical(0),
    has_valid_breaks = logical(0),
    break_index      = integer(0),
    break_date       = as.Date(character(0)),
    break_magn       = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Loop over each algorithm (== method name)
  for (alg_name in names(ts_breaks_obj$algorithms)) {
    alg_runs <- ts_breaks_obj$algorithms[[alg_name]]
    
    for (run_name in names(alg_runs)) {
      run <- alg_runs[[run_name]]
      
      # If the run has breaks, expand to multiple rows
      if (run$has_breaks && length(run$breaks_indices) > 0) {
        if (length(run$breaks_indices) != length(run$breaks_dates)) {
          stop(paste("Mismatch between breaks_indices and breaks_dates in",
                     alg_name, run_name))
        }
        
        run_df <- data.frame(
          algorithm        = rep(alg_name, length(run$breaks_indices)),
          run_id           = rep(run_name, length(run$breaks_indices)),
          method           = rep(run$method, length(run$breaks_indices)),
          data_type        = rep(run$data_type, length(run$breaks_indices)),
          has_breaks       = rep(TRUE, length(run$breaks_indices)),
          has_valid_breaks = rep(run$has_valid_breaks, length(run$breaks_indices)),
          break_index      = run$breaks_indices,
          break_date       = as.Date(run$breaks_dates),
          break_magn       = as.numeric(run$break_magn),
          stringsAsFactors = FALSE
        )
      } else {
        # If no breaks or an empty list of breaks
        run_df <- data.frame(
          algorithm        = alg_name,
          run_id           = run_name,
          method           = run$method,
          data_type        = run$data_type,
          has_breaks       = FALSE,
          has_valid_breaks = FALSE,
          break_index      = NA_integer_,
          break_date       = as.Date(NA),
          break_magn       = NA,
          stringsAsFactors = FALSE
        )
      }
      
      breaks_df <- rbind(breaks_df, run_df)
    }
  }
  
  return(breaks_df)
}


