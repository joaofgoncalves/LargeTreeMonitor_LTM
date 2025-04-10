
## TreeTrakR / BigTreeTrakR
## TreeWatchR / BigTreeWatchR


library(sf)
library(terra)
library(rgee)
library(tidyverse)
library(rgeeExtra)
library(imputeTS)

ee_Initialize(user = "joaofgoncalves", gcs = TRUE, drive = TRUE)

#latitude <- 39.39
#longitude <- -8.22

latitude <- 41.720898
longitude <- -8.747039


# Create the point geometry for the centroid
centroid_point <- ee$Geometry$Point(c(longitude, latitude))


buffer_distance <- 20

# Create a buffer around the point (in meters)
buffer_region <- centroid_point$buffer(buffer_distance)

# Print the point object to check
print(centroid_point)

# Define the time range for the year 2023
start_date <- '2015-01-01'
end_date <- '2024-12-31'

# Load the Sentinel-2 L2A image collection (surface reflectance product)
sentinel_collection <- ee$ImageCollection('COPERNICUS/S2_SR_HARMONIZED')$
  filterBounds(centroid_point)$
  filterDate(start_date, end_date)
  #filter(ee$Filter$lt('SCL', 4)) # Filter out clouds and shadows, SCL = Scene Classification Layer



###################################################
## POINT SAMPLE
###################################################


# Function to compute NDVI and extract value at the point
compute_ndvi_at_point <- function(image) {
  
  ndvi <- image$normalizedDifference(c('B8', 'B4'))$rename('NDVI') # B8 = NIR, B4 = Red
  
  value <- ndvi$sample(region = centroid_point, scale = 10, projection = 'EPSG:32629')$
    first()
    #get('NDVI')
  return(value)
}

# Map the NDVI extraction function over the image collection
ndvi_values <- ee$FeatureCollection(sentinel_collection$map(compute_ndvi_at_point))

print(ndvi_values)

# Convert the result to an R data frame for easier handling
ndvi_values_list <- ee_as_sf(ndvi_values)



###################################################
## BUFFER REGION
###################################################


# Function to compute the mean NDVI for the buffer region
compute_mean_ndvi <- function(image) {
  # Compute NDVI: NDVI = (NIR - Red) / (NIR + Red)
  ndvi <- image$normalizedDifference(c('B8', 'B4'))$rename('NDVI') # B8 = NIR, B4 = Red
  
  # Calculate the mean NDVI over the buffer region
  mean_ndvi <- ndvi$reduceRegion(
    reducer = ee$Reducer$mean(),          # Use the mean reducer
    geometry = buffer_region,             # Apply it over the buffer region
    scale = 10,                           # Pixel scale
    maxPixels = 1e9                       # Max pixels to consider
  )
  
  # Return the mean NDVI value
  return(ee$Feature(NULL, mean_ndvi))
}

# Map the NDVI aggregation function over the image collection
ndvi_mean_features <- sentinel_collection$map(compute_mean_ndvi)

# Convert the results to a Feature Collection
ndvi_feature_collection <- ee$FeatureCollection(ndvi_mean_features)

# Print the NDVI feature collection
print(ndvi_feature_collection)

# Convert the result to an R data frame for easier handling
# (Using ee_as_sf to convert ee feature collection to sf object in R)
ndvi_mean_values_sf <- ee_as_sf(ndvi_feature_collection)


###################################################

cbind(ndvi_values_list,ndvi_mean_values_sf)



###################################################


# Print the time series of NDVI values
print(ndvi_values_list)


dts <- ee_get_date_ic(sentinel_collection, time_end = FALSE)


dt <- data.frame(ndvi = ndvi_values_list$NDVI,ts = dts$time_start)

g <- ggplot(dt,aes(x=ts, y=ndvi)) + 
  geom_point() + 
  geom_line()

plot(g)



###################################################



      
# Combine into a data frame
df <- data.frame(ts = dt$ts, ndvi = dt$ndvi)

df <- df %>% 
  mutate(dts = as.Date(ts)) %>% 
  arrange(ts)


start_date <- floor_date(min(df$ts), unit = "day")
end_date <- ceiling_date(max(df$ts), unit = "day")
regular_dates <- seq(from = start_date, to = end_date, by = "1 days")

# Create a data frame with regular timestamps
regular_df <- data.frame(ts = regular_dates)


df_regular <- full_join(regular_df, df, by = c("ts"="dts")) %>% 
  arrange(ts)

#df_regular$ndvi[df_regular$ndvi < 0.1] <- NA


df_regular$ndvi_interp <- na_interpolation(df_regular$ndvi, option = "linear")


# 
# df_regular$ndvi_interp[df_regular$ndvi_interp < 0.15] <- NA
# 
# 
# df_regular$ndvi_interp <- na_seadec(df_regular$ndvi_interp, find_frequency = TRUE)
# 



ggplot(df_regular, aes(x = `ts`)) +
  geom_line(aes(y = ndvi_interp), color = "blue", size = 1) +
  geom_point(aes(y = ndvi_interp), color = "red", size = 2) +
  geom_point(aes(y = ndvi), color = "yellow", size = 2) +
  
  labs(
    title = "Regularized NDVI Time Series (5-Day Intervals)",
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()



library(phenofit)


df_regular <- df_regular %>%
  rename(y = ndvi_interp) %>%
  select(ts, y)


# Define a threshold NDVI value
threshold <- quantile(df_regular$y,probs=0.35)

# Assign weights: 1 for NDVI >= threshold, lower weight for NDVI < threshold
df_regular <- df_regular %>%
  mutate(
    w = case_when(
      y >= threshold ~ 1,
      !is.na(y) ~ 0.01,  # Lower weight for potential cloud-contaminated values
      TRUE ~ 0          # Weight zero for missing values
    )
  )

# 
# df_regular <- df_regular %>%
#   mutate(
#     w = ifelse(!is.na(y), y / max(y, na.rm = TRUE), 0)
#   )
# 
# 
# 
# y_min <- min(df_regular$y, na.rm = TRUE)
# y_max <- max(df_regular$y, na.rm = TRUE)
# 
# # Rescale NDVI values to [0, 1] to use as weights
# df_regular <- df_regular %>%
#   mutate(
#     w = (y - y_min) / (y_max - y_min)
#   )

# Prepare the data list as required by phenofit
input_data <- df_regular %>%
  select(ts, y, w) %>%
  rename(t = ts, y = y, w = w)

# Convert date to numeric time (e.g., days since start)
input_data$t_num <- as.numeric(difftime(input_data$t, min(input_data$t), units = "days"))

smoothed_result <- whit2(y = df_regular$y, w = df_regular$w, lambda = 5000)

smoothed_result2 <- whit2(y = df_regular$y, lambda = 5000)


# Add the smoothed values to your data frame
df_regular$ndvi_smooth <- smoothed_result
df_regular$ndvi_smooth2 <- smoothed_result2



library(slider)

# Calculate the maximum NDVI in a forward 30-day window
df_regular <- df_regular %>%
  mutate(
    ndvi_mov_qt = slide_index_dbl(
      .x = y,
      .i = ts,
      .f = ~ quantile(.x, probs=0.95, na.rm = TRUE),
      .before = days(5),
      .after = days(4),
      .complete = FALSE
    )
  )


smoothed_result3<- whit2(y = df_regular$ndvi_mov_qt, lambda = 1E5)
df_regular$ndvi_smooth3 <- smoothed_result3


ggplot(df_regular, aes(x = ts)) +
  geom_line(aes(y = y), color = "grey70", size = 0.1, na.rm = TRUE) +
  geom_point(aes(y = y), color = "darkgreen", size = 0.5, na.rm = TRUE) +
  #geom_line(aes(y = ndvi_smooth), color = "blue", size = 1) +
  #geom_line(aes(y = ndvi_smooth2), color = "orange", size = 1) +
  geom_line(aes(y = ndvi_smooth3), color = "red", size = 1) +
  labs(
    title = "NDVI Time Series with Whittaker Smoother",
    x = "Date",
    y = "NDVI"
  ) +
  theme_minimal()



# Generic function to create a time series object from dates and values (adjusted for satellite data)
create_ts_from_dates <- function(dates, values) {
  # Ensure dates are of Date class
  dates <- as.Date(dates)
  
  # Order dates and corresponding values
  ord <- order(dates)
  dates <- dates[ord]
  values <- values[ord]
  
  # Calculate median difference between dates to determine frequency
  median_diff <- median(diff(dates))
  
  print(median_diff)
  
  # Calculate frequency based on actual median interval (use median_diff in days)
  freq <- round(365 / as.numeric(median_diff))
  
  # Determine start year and fractional period based on median interval
  start_year <- as.numeric(format(min(dates), "%Y"))
  days_into_year <- as.numeric(format(min(dates), "%j"))
  
  # Calculate fractional start period
  start_period <- 1 + (days_into_year - 1) / as.numeric(median_diff)
  
  ts_obj <- ts(values, start = c(start_year, start_period), frequency = freq)
  
  return(ts_obj)
}

# Subset data by code

create_ts_from_dates(dates = df_regular$ts, values = df_regular$ndvi_smooth3)


# Make the time series
yts <- ts(df_regular$ndvi_smooth3, start = c(2015,37), frequency = 73)



ecpDivisiveBreakDetect(yts,
                                   sig.lvl   = 0.05, 
                                   R         = 1000, 
                                   k         = 1, 
                                   min.size  = 30, 
                                   alpha     = 1, 
                                   seasonAdj = TRUE,
                                   threshPerc    = 15, 
                                   detectionTime = 185, 
                                   halfLen       = 20, ...)


