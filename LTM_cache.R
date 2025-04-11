

# Function to generate cache file name
ltm_cache_file_name <- function(lat, lon, start_date, end_date, spi, proc_level) {
  start_date <- format(as.Date(start_date), "%Y-%m-%d")
  end_date   <- format(as.Date(end_date),   "%Y-%m-%d")
  name_base  <- sprintf("ltm_%s_%s_%s_%s_%s_%s",
                        round(lat, 5), round(lon, 5),
                        start_date, end_date,
                        spi, proc_level)
  paste0("ltm_cache/", name_base, ".rds")
}

# Check cache function
ltm_check_cache <- function(lat, lon, start_date, end_date, spi, proc_level) {
  start_date <- format(as.Date(start_date), "%Y-%m-%d")
  end_date   <- format(as.Date(end_date),   "%Y-%m-%d")
  pattern    <- sprintf("ltm_%s_%s_%s_%s_%s_%s\\.rds$",
                        round(lat, 5), round(lon, 5),
                        start_date, end_date,
                        spi, proc_level)
  files <- list.files("ltm_cache", pattern = pattern, full.names = TRUE)
  if (length(files) > 0) return(files[1])
  return(NULL)
}

