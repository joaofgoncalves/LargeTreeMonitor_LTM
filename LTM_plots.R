# 
# 
# ltm_plot_spidf_ts <- function(spidf_obj, tree_id="") {
#   stopifnot(inherits(spidf_obj, "spidf"))
#   
#   # Build a named vector of possible series
#   series_cols <- c(
#     "spi" = ifelse(is_regularized(spidf_obj), "Daily interpolated", "Original SPI"),
#     "spi_mov_wind" = "Moving-wind. quantile",
#     "spi_smooth" = "Whittaker-smooth.",
#     "spi_mov_smooth" = "Whittaker-smooth. mov. wind."
#   )
#   
#   # Always plot at least the spi column
#   selected_cols <- intersect(names(series_cols), names(spidf_obj))
#   
#   # Prepare long-format data
#   df_long <- spidf_obj %>%
#     select(ti, all_of(selected_cols)) %>%
#     pivot_longer(
#       cols = all_of(selected_cols),
#       names_to = "type",
#       values_to = "value",
#       values_drop_na = TRUE
#     ) %>%
#     mutate(
#       type = recode(type, !!!series_cols)
#     )
#   
#   # Define consistent color palette
#   color_map <- c(
#     "Original SPI" = "#7570b3",
#     "Daily interpolated" = "#2c7fb8",
#     "Moving-wind. quantile" = "red",
#     "Whittaker-smooth." = "darkgreen",
#     "Whittaker-smooth. mov. wind." = "orange"
#   )
#   
#   # Plot
#   ggplot(df_long, aes(x = ti, y = value, color = type)) +
#     geom_line(linewidth = 0.6) +
#     geom_point(size = 1, alpha = 0.7) +
#     scale_color_manual(
#       name = "Series Type",
#       values = color_map[names(color_map) %in% unique(df_long$type)]
#     ) +
#     labs(
#       title = paste("Time Series of", get_spi(spidf_obj)),
#       subtitle = paste0(
#         "Location: [", round(get_longitude(spidf_obj), 5), ", ", round(get_latitude(spidf_obj), 5), "] | ",
#         "Period: ", format(get_range_start(spidf_obj)), " to ", format(get_range_end(spidf_obj)), " | ",
#         "Proc. Level: ", get_proc_level(spidf_obj),
#         ifelse(tree_id!=""&&tree_id!="---",paste0("| Tree ID: ",tree_id),"")
#       ),
#       x = "Date",
#       y = get_spi(spidf_obj)
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(legend.position = "bottom")
# }
# 
# 
# 
# 

# In LTM_plots.R (or wherever your plotting code resides)
ltm_plot_spidf_ts <- function(spidf_obj, tree_id = "",
                              df_breaks = NULL) {
  stopifnot(inherits(spidf_obj, "spidf"))
  
  # 1) Define labels for your known time-series columns
  series_cols <- c(
    "spi"            = ifelse(is_regularized(spidf_obj), "Daily interpolated", "Original SPI"),
    "spi_mov_wind"   = "Moving-wind. quantile",
    "spi_smooth"     = "Whittaker-smooth.",
    "spi_mov_smooth" = "Whittaker-smooth. mov. wind."
  )
  
  # 2) Identify which series columns exist in spidf_obj
  selected_cols <- intersect(names(series_cols), names(spidf_obj))
  
  # 3) Reshape data to long format for ggplot
  df_long <- spidf_obj %>%
    dplyr::select(ti, dplyr::all_of(selected_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(selected_cols),
      names_to = "type",
      values_to = "value",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(
      type = dplyr::recode(type, !!!series_cols)
    )
  
  # ------------------------------------------------------
  # 4) Optional break data

  # ------------------------------------------------------
  df_plot_breaks <- NULL
  if (!is.null(df_breaks)) {
    df_plot_breaks <- df_breaks %>%
      dplyr::filter(
        .data[["has_breaks"]] == TRUE,
        .data[["has_valid_breaks"]] == TRUE
      )
    # Ensure Break Date is a true Date
    if ("break_date" %in% colnames(df_plot_breaks)) {
      df_plot_breaks[["break_date"]] <- as.Date(df_plot_breaks[["break_date"]])
    }
  }
  
  # ------------------------------------------------------
  # 5) Define color palettes for time series & algorithms
  #    Okabe-Ito palette: color-blind-safe
  # ------------------------------------------------------
  color_map_ts <- c(
    "Original SPI"                  = "#0072B2",  # Blue
    "Daily interpolated"            = "#56B4E9",  # Sky Blue
    "Moving-wind. quantile"         = "#D55E00",  # Vermilion
    "Whittaker-smooth."             = "#009E73",  # Bluish Green
    "Whittaker-smooth. mov. wind."  = "#E69F00"   # Orange
  )
  
  color_map_algs <- c(
    "cpm"     = "#CC79A7",  # Reddish Purple
    "ed"      = "#F0E442",  # Yellow
    "bfast01" = "#999999",  # Gray
    "mcp"     = "#000000"   # Black
  )
  
  used_ts_types <- unique(df_long$type)  
  used_algs     <- if (!is.null(df_plot_breaks)) unique(df_plot_breaks[["algorithm"]]) else character(0)
  
  color_map_ts   <- color_map_ts[names(color_map_ts) %in% used_ts_types]
  color_map_algs <- color_map_algs[names(color_map_algs) %in% used_algs]
  
  # Merge to single palette
  combined_color_map <- c(color_map_ts, color_map_algs)
  
  # ------------------------------------------------------
  # 6) Build the ggplot (time-series lines)
  # ------------------------------------------------------
  p <- ggplot(df_long, aes(x = ti, y = value, color = type)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1, alpha = 0.7)
  
  # ------------------------------------------------------
  # 7) Add vertical lines for valid breaks (if any)
  # ------------------------------------------------------
  if (!is.null(df_plot_breaks) && nrow(df_plot_breaks) > 0) {
    p <- p +
      geom_vline(
        data = df_plot_breaks,
        aes(xintercept = break_date, color = algorithm),
        #inherit.aes = FALSE,
        size = 0.8,
        alpha = 0.6
      )
  }
  
  # ------------------------------------------------------
  # 8) Finish the plot with combined scale + labels
  # ------------------------------------------------------
  p <- p +
    scale_color_manual(values = combined_color_map) +
    labs(
      color = "Legend",
      title = paste("Time Series of", get_spi(spidf_obj)),
      subtitle = paste0(
        "Location: [", round(get_longitude(spidf_obj), 5), ", ",
        round(get_latitude(spidf_obj), 5), "] | ",
        "Period: ", format(get_range_start(spidf_obj)), " to ",
        format(get_range_end(spidf_obj)), " | ",
        "Proc. Level: ", get_proc_level(spidf_obj),
        ifelse(tree_id != "" && tree_id != "---",
               paste0(" | Tree ID: ", tree_id),
               "")
      ),
      x = "Date",
      y = get_spi(spidf_obj)
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  return(p)
}
