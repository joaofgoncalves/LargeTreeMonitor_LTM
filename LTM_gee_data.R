

SPECTRAL_INDICES_LIST <<- c("NDVI","EVI")

PROC_LEVELS_LIST <<- c("L1","L1C","L2","L2A")

ltm_get_s2_imgcol <- function(proc_level = "L2A") {
  
  ## ----------------------------------------------------------------------- ##
  ## Sentinel-2a/b mission ----
  ## ----------------------------------------------------------------------- ##
  
  if (proc_level %in% c("L2A", "L2")) {
    
    # Harmonized Sentinel-2 MSI: MultiSpectral Instrument, Level-2A
    imCol <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED") %>%
      ee$ImageCollection$select(c("B2", "B3", "B4", "B5", "B6", "B7",
                                  "B8", "B8A", "B11", "B12", "QA60"),
                                # Renamed bands
                                c("Blue", "Green", "Red", "RE1", "RE2", "RE3",
                                  "NIR", "RE4", "SWIR1", "SWIR2", "QA60"))
    return(imCol)
    
  } else if (proc_level %in% c("L1C", "L1")) {
    
    # Harmonized Sentinel-2 MSI: MultiSpectral Instrument, Level-1C
    imCol <- ee$ImageCollection("COPERNICUS/S2_HARMONIZED") %>%
      ee$ImageCollection$select(c("B2","B3","B4","B5","B6","B7",
                                  "B8","B8A","B11","B12","QA60"),
                                # Renamed bands
                                c("Blue", "Green", "Red","RE1","RE2",
                                  "RE3","NIR","RE4","SWIR1", "SWIR2","QA60"))
    return(imCol)
    
  } else {
    stop("Sentinel-2 - Processing level value in proc_level is not supported")
  }
}


ltm_s2_mask_clouds <- function(img){
  
  # Select quality layer
  qa <- img$select("QA60")
  
  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloudBitMask <- ee$Number(1 %<<% 10)  # 1 << 10
  cirrusBitMask <- ee$Number(1 %<<% 11) # 1 << 11
  
  # Both flags should be set to zero, indicating clear conditions.
  mask <- ((qa$bitwiseAnd(cloudBitMask))$eq(0))$And(qa$bitwiseAnd(cirrusBitMask)$eq(0))
  
  return(ee$Image(img$updateMask(mask)))
}

ltm_s2_clouds <- function(img){
  
  # Select QA60 band
  qa <- img$select("QA60")
  
  # Define bit masks for clouds and cirrus
  cloudBitMask <- ee$Number(1 %<<% 10)  # 1 << 10
  cirrusBitMask <- ee$Number(1 %<<% 11) # 1 << 11
  
  # Identify cloud pixels (clouds or cirrus)
  cloud_mask <- qa$bitwiseAnd(cloudBitMask)$neq(0)$Or(qa$bitwiseAnd(cirrusBitMask)$neq(0))
  
  # Create binary mask: 1 for cloud pixels, 0 for clear pixels
  binary_cloud_mask <- cloud_mask$rename("cloud_mask")$uint8()
  
  # Preserve original image properties and timestamp
  return(binary_cloud_mask$copyProperties(img, img$propertyNames())$set('system:time_start', img$get('system:time_start')))
  
}


# Function to compute NDVI and extract value at the point
ltm_calc_ndvi <- function(img) {
  
  ndvi <- img$normalizedDifference(c('NIR', 'Red'))$rename('NDVI')
  ndvi <- ndvi$copyProperties(img, img$propertyNames())$set('system:time_start', img$get('system:time_start'))
  
  qa <- img$select("QA60")
  cloudBitMask <- ee$Number(1 %<<% 10)  # 1 << 10
  cirrusBitMask <- ee$Number(1 %<<% 11) # 1 << 11
  cloud_mask <- qa$bitwiseAnd(cloudBitMask)$neq(0)$Or(qa$bitwiseAnd(cirrusBitMask)$neq(0))
  binary_cloud_mask <- cloud_mask$rename("cloud_mask")$uint8()
  binary_cloud_mask <- binary_cloud_mask$copyProperties(img, img$propertyNames())$set('system:time_start', img$get('system:time_start'))
  
  return(ee$Image(ndvi)$addBands(ee$Image(binary_cloud_mask)))
}

ltm_calc_evi <- function(img){
  
  evi <- ee$Image(img$expression(
    '2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))',
    list(
      NIR = img$select('NIR'),
      Red = img$select('Red'),
      Blue = img$select('Blue')
    ))
  )$rename('EVI')
  evi <- evi$copyProperties(img, img$propertyNames())$set('system:time_start', img$get('system:time_start'))
  
  qa <- img$select("QA60")
  cloudBitMask <- ee$Number(1 %<<% 10)  # 1 << 10
  cirrusBitMask <- ee$Number(1 %<<% 11) # 1 << 11
  cloud_mask <- qa$bitwiseAnd(cloudBitMask)$neq(0)$Or(qa$bitwiseAnd(cirrusBitMask)$neq(0))
  binary_cloud_mask <- cloud_mask$rename("cloud_mask")$uint8()
  binary_cloud_mask <- binary_cloud_mask$copyProperties(img, img$propertyNames())$set('system:time_start', img$get('system:time_start'))
  
  return(ee$Image(evi)$addBands(ee$Image(binary_cloud_mask)))
}




get_timerange <- function(x) {
  #stopifnot(inherits(x, "spidf"))
  
  if (!"ti" %in% names(x)) {
    stop("The object does not contain a 'ti' column.")
  }
  
  # Convert to Date
  ti_dates <- as.Date(x$ti)
  
  # Handle empty or invalid cases
  if (length(ti_dates) == 0 || all(is.na(ti_dates))) {
    return(c(start_date = NA_character_, end_date = NA_character_))
  }
  
  range_dates <- range(ti_dates, na.rm = TRUE)
  
  return(c(
    start_date = as.Date(format(range_dates[1], "%Y-%m-%d")),
    end_date   = as.Date(format(range_dates[2], "%Y-%m-%d"))
  ))
}


ltm_s2_get_data_point <- function(lat, lon, start_date, end_date, spi = "NDVI",
                                  proc_level = "L2A", crs_code = "EPSG:4326", 
                                  rm_duplicates = TRUE, tree_id = NULL){
  
  
  if(!(spi %in% SPECTRAL_INDICES_LIST)){
    stop("The spi is not listed as a valid spectral index. Available indices are: ",
         paste(SPECTRAL_INDICES_LIST, collapse=", "))
  }
  
  if(!(proc_level %in% PROC_LEVELS_LIST)){
    stop("The proc_level is not listed as a valid spectral index. Available indices are: ",
         paste(PROC_LEVELS_LIST, collapse=", "))
  }
  
  if(ltm_check_gee_status() != "CONNECTED"){
    stop("Google Earth Engine (GEE) has not been initialized. Please run ltm_start_gee() 
         before calling this function.")
  }
  
  # Create the point geometry to extract data
  target_point <- ee$Geometry$Point(c(lon, lat))
  
  # Load the Sentinel-2 L2A image collection (surface reflectance product)
  s2_collection <- 
    ltm_get_s2_imgcol(proc_level) %>% 
    ee$ImageCollection$filterBounds(target_point) %>% 
    ee$ImageCollection$filterDate(start_date, end_date)
  
  
  dts <- try(ee_get_date_ic(s2_collection, time_end = FALSE), silent = TRUE)
  
  if(inherits(dts, "try-error")){
    stop("An error occurred while getting image dates from GEE", call. = FALSE)
  }
  
  if(spi == "NDVI"){
    # Map the NDVI extraction function over the image collection
    s2_collection <- s2_collection %>%
      ee$ImageCollection$map(ltm_calc_ndvi)
  }else if(spi == "EVI"){
    # Map the EVI extraction function over the image collection
    s2_collection <- s2_collection %>%
      ee$ImageCollection$map(ltm_calc_evi)
  }else{
    stop("Spectral index defined in spi parameter does not exist")
  }

  ## Spectral index
  
  # Map using reduceRegion instead of sample
  s2_spi_features <- s2_collection$map(function(img) {

      img_val <- img$reduceRegion(
      reducer = ee$Reducer$first(),
      geometry = target_point,
      scale = 10,
      crs = crs_code
    )
    
    # Set the result as property of a feature
    ee$Feature(NULL)$copyProperties(img, img$propertyNames())$set(img_val)
  })
  
  # Filter out any features where 'cloud_mask' is null
  spi_values <- ee$FeatureCollection(s2_spi_features$filter(ee$Filter$notNull(list("cloud_mask",spi))))
  

  # Convert safely to sf - get data from point
  spi_values_list <- try(ee_as_sf(spi_values) %>% st_drop_geometry(), silent = TRUE)
  
  #spi_mask_list <- try(ee_as_sf(s2_mask_values), silent = TRUE)
  
  
  if(inherits(spi_values_list, "try-error")){
    stop("An error occurred while getting point spectral data from GEE", call. = TRUE)
  }

  dt <- left_join(dts,spi_values_list[,c(spi,"cloud_mask","system.id")], 
                  by=c("id"="system.id"))
  
  colnames(dt) <- c("id","ti","spi","cloud_mask")
  
  dt <- dt %>% 
    mutate(masked_vals = ifelse(cloud_mask==1, NA, spi)) %>% 
    select(1,2,5,3,4)
  
  
  # Assign metadata as attributes
  attr(dt, "lat") <- lat
  attr(dt, "lon") <- lon
  attr(dt, "start_date") <- as.Date(start_date)
  attr(dt, "end_date") <- as.Date(end_date)

  attr(dt, "tree_id") <- tree_id
  
  attr(dt, "range_start") <- as.Date(get_timerange(dt)[1])
  attr(dt, "range_end") <- as.Date(get_timerange(dt)[2])
  
  attr(dt, "spi") <- spi
  attr(dt, "proc_level") <- proc_level
  attr(dt, "crs_code") <- crs_code
  
  attr(dt, "regularized") <- FALSE
  
  attr(dt, "mov_window") <- FALSE
  attr(dt, "mov_window_quantile") <- NA
  attr(dt, "mov_window_size_days") <- NA
  
  attr(dt, "whit_smoothing") <- FALSE
  attr(dt, "whit_lambda") <- NA
  attr(dt, "whit_quantile_threshold") <- NA
  attr(dt, "whit_weights_used") <- NA
  
  # Assign class
  class(dt) <- c("spidf", class(dt))
  
  if(rm_duplicates){
    dt <- ltm_remove_duplicates(dt)
  }
    
  return(dt)
}

print.spidf <- function(x, ...) {
  cat("Spectral Index Data Frame (spidf):\n")
  cat(" Location (lon/lat): [", get_longitude(x), ", ", get_latitude(x), "]\n", sep = "")
  cat(" Period asked: ", format(get_start_date(x)), " to ", format(get_end_date(x)), "\n", sep = "")
  cat(" Period available: ", format(get_range_start(x)), " to ", format(get_range_end(x)), "\n", sep = "")
  cat(" Spectral Index: ", get_spi(x), "\n", sep = "")
  cat(" Processing Level: ", get_proc_level(x), "\n", sep = "")
  cat(" Coordinate System: ", get_crs_code(x), "\n", sep = "")
  cat(" Is the time series regularized?: ", ifelse(isTRUE(is_regularized(x)), "Yes", "No"), "\n", sep = "")
  cat(" Moving window?: ", ifelse(isTRUE(has_mov_window(x)), "Yes", "No"), "\n", sep = "")
  cat(" Whittaker smoothing?: ", ifelse(isTRUE(is_whit_smoothed(x)), "Yes", "No"), "\n\n", sep = "")
  
  # Call default data.frame print method
  NextMethod("print", x, ...)
}


## -------------------------------------------------------------
## GET Functions 
## -------------------------------------------------------------

# Core metadata
get_latitude <- function(x) UseMethod("get_latitude")
get_latitude.spidf <- function(x) attr(x, "lat")

get_longitude <- function(x) UseMethod("get_longitude")
get_longitude.spidf <- function(x) attr(x, "lon")

get_start_date <- function(x) UseMethod("get_start_date")
get_start_date.spidf <- function(x) attr(x, "start_date")

get_end_date <- function(x) UseMethod("get_end_date")
get_end_date.spidf <- function(x) attr(x, "end_date")

get_spi <- function(x) UseMethod("get_spi")
get_spi.spidf <- function(x) attr(x, "spi")

get_proc_level <- function(x) UseMethod("get_proc_level")
get_proc_level.spidf <- function(x) attr(x, "proc_level")

get_crs_code <- function(x) UseMethod("get_crs_code")
get_crs_code.spidf <- function(x) attr(x, "crs_code")

# State checkers
is_regularized <- function(x) UseMethod("is_regularized")
is_regularized.spidf <- function(x) attr(x, "regularized")

get_regularize_method <- function(x) UseMethod("get_regularize_method")
get_regularize_method.spidf <- function(x) attr(x, "regularize_method")

has_mov_window <- function(x) UseMethod("has_mov_window")
has_mov_window.spidf <- function(x) attr(x, "mov_window")

is_whit_smoothed <- function(x) UseMethod("is_whit_smoothed")
is_whit_smoothed.spidf <- function(x) attr(x, "whit_smoothing")

# Moving window metadata
get_mov_window_quantile <- function(x) UseMethod("get_mov_window_quantile")
get_mov_window_quantile.spidf <- function(x) attr(x, "mov_window_quantile")

get_mov_window_size_days <- function(x) UseMethod("get_mov_window_size_days")
get_mov_window_size_days.spidf <- function(x) attr(x, "mov_window_size_days")

# Whittaker smoothing metadata
get_whit_lambda <- function(x) UseMethod("get_whit_lambda")
get_whit_lambda.spidf <- function(x) attr(x, "whit_lambda")

get_whit_quantile_threshold <- function(x) UseMethod("get_whit_quantile_threshold")
get_whit_quantile_threshold.spidf <- function(x) attr(x, "whit_quantile_threshold")

get_whit_weights_used <- function(x) UseMethod("get_whit_weights_used")
get_whit_weights_used.spidf <- function(x) attr(x, "whit_weights_used")

# Get start of actual time series range
get_range_start <- function(x) UseMethod("get_range_start")
get_range_start.spidf <- function(x) as.Date(attr(x, "range_start"))

# Get end of actual time series range
get_range_end <- function(x) UseMethod("get_range_end")
get_range_end.spidf <- function(x) as.Date(attr(x, "range_end"))

get_tree_id <- function(x) UseMethod("get_tree_id")
get_tree_id.spidf <- function(x) attr(x, "tree_id")


## -------------------------------------------------------------
## SET Functions 
## -------------------------------------------------------------

# Set latitude attribute
set_latitude <- function(x, value) UseMethod("set_latitude")
set_latitude.spidf <- function(x, value) {
  attr(x, "lat") <- value
  x
}

# Set longitude attribute
set_longitude <- function(x, value) UseMethod("set_longitude")
set_longitude.spidf <- function(x, value) {
  attr(x, "lon") <- value
  x
}

# Set start_date attribute
set_start_date <- function(x, value) UseMethod("set_start_date")
set_start_date.spidf <- function(x, value) {
  attr(x, "start_date") <- as.Date(value)
  x
}

# Set end_date attribute
set_end_date <- function(x, value) UseMethod("set_end_date")
set_end_date.spidf <- function(x, value) {
  attr(x, "end_date") <- as.Date(value)
  x
}

# Set spi attribute
set_spi <- function(x, value) UseMethod("set_spi")
set_spi.spidf <- function(x, value) {
  attr(x, "spi") <- value
  x
}

# Set proc_level attribute
set_proc_level <- function(x, value) UseMethod("set_proc_level")
set_proc_level.spidf <- function(x, value) {
  attr(x, "proc_level") <- value
  x
}

# Set crs_code attribute
set_crs_code <- function(x, value) UseMethod("set_crs_code")
set_crs_code.spidf <- function(x, value) {
  attr(x, "crs_code") <- value
  x
}

# Set regularized attribute
set_regularized <- function(x, value) UseMethod("set_regularized")
set_regularized.spidf <- function(x, value) {
  attr(x, "regularized") <- as.logical(value)
  x
}

# Set regularize_method attribute
set_regularize_method <- function(x, value) UseMethod("set_regularize_method")
set_regularize_method.spidf <- function(x, value) {
  attr(x, "regularize_method") <- value
  x
}

# Set mov_window attribute
set_mov_window <- function(x, value) UseMethod("set_mov_window")
set_mov_window.spidf <- function(x, value) {
  attr(x, "mov_window") <- as.logical(value)
  x
}

# Set whit_smoothing attribute
set_whit_smoothed <- function(x, value) UseMethod("set_whit_smoothed")
set_whit_smoothed.spidf <- function(x, value) {
  attr(x, "whit_smoothing") <- as.logical(value)
  x
}

# Set mov_window_quantile attribute
set_mov_window_quantile <- function(x, value) UseMethod("set_mov_window_quantile")
set_mov_window_quantile.spidf <- function(x, value) {
  attr(x, "mov_window_quantile") <- as.numeric(value)
  x
}

# Set mov_window_size_days attribute
set_mov_window_size_days <- function(x, value) UseMethod("set_mov_window_size_days")
set_mov_window_size_days.spidf <- function(x, value) {
  attr(x, "mov_window_size_days") <- as.integer(value)
  x
}

# Set whit_lambda attribute
set_whit_lambda <- function(x, value) UseMethod("set_whit_lambda")
set_whit_lambda.spidf <- function(x, value) {
  attr(x, "whit_lambda") <- as.numeric(value)
  x
}

# Set whit_quantile_threshold attribute
set_whit_quantile_threshold <- function(x, value) UseMethod("set_whit_quantile_threshold")
set_whit_quantile_threshold.spidf <- function(x, value) {
  attr(x, "whit_quantile_threshold") <- as.numeric(value)
  x
}

# Set whit_weights_used attribute
set_whit_weights_used <- function(x, value) UseMethod("set_whit_weights_used")
set_whit_weights_used.spidf <- function(x, value) {
  attr(x, "whit_weights_used") <- as.logical(value)
  x
}

# Set range_start attribute
set_range_start <- function(x, value) UseMethod("set_range_start")
set_range_start.spidf <- function(x, value) {
  attr(x, "range_start") <- as.Date(value)
  x
}

# Set range_end attribute
set_range_end <- function(x, value) UseMethod("set_range_end")
set_range_end.spidf <- function(x, value) {
  attr(x, "range_end") <- as.Date(value)
  x
}

# Set tree_id attribute
set_tree_id <- function(x, value) UseMethod("set_tree_id")
set_tree_id.spidf <- function(x, value) {
  attr(x, "tree_id") <- value
  x
}
