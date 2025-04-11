

ltm_plot_spidf_ts <- function(spidf_obj, tree_id = "",
                              df_breaks = NULL) {
  stopifnot(inherits(spidf_obj, "spidf"))
  
  # Define labels for known time-series columns
  series_cols <- c(
    "spi"            = ifelse(is_regularized(spidf_obj), "Daily interpolated", "Original SPI"),
    "spi_mov_wind"   = "Moving-wind. quantile",
    "spi_smooth"     = "Whittaker-smooth.",
    "spi_mov_smooth" = "Whittaker-smooth. mov. wind."
  )
  
  # Identify existing series columns
  selected_cols <- intersect(names(series_cols), names(spidf_obj))
  
  # Reshape to long format
  df_long <- spidf_obj %>%
    dplyr::select(ti, dplyr::all_of(selected_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(selected_cols),
      names_to = "type",
      values_to = "value",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate(type = dplyr::recode(type, !!!series_cols))
  
  # Optional breaks data
  df_plot_breaks <- NULL
  if (!is.null(df_breaks)) {
    df_plot_breaks <- df_breaks %>%
      dplyr::filter(has_breaks & has_valid_breaks) %>%
      dplyr::mutate(break_date = as.Date(break_date))
  }
  
  # Color palettes
  color_map_ts <- c(
    "Original SPI"                  = "#0072B2",
    "Daily interpolated"            = "#56B4E9",
    "Moving-wind. quantile"         = "#D55E00",
    "Whittaker-smooth."             = "#009E73",
    "Whittaker-smooth. mov. wind."  = "#E69F00"
  )
  
  color_map_algs <- c(
    "cpm"     = "#CC79A7",
    "ed"      = "#F0E442",
    "bfast01" = "#999999",
    "mcp"     = "#000000",
    "stc"     = "#882255",
    "wbs"     = "#117733"
  )
  
  line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  
  # Filter colors
  color_map_ts   <- color_map_ts[names(color_map_ts) %in% unique(df_long$type)]
  color_map_algs <- color_map_algs[names(color_map_algs) %in% unique(df_plot_breaks$algorithm)]
  line_types <- line_types[1:length(unique(df_plot_breaks$algorithm))]
  
  # Plot
  p <- ggplot(df_long, aes(x = ti, y = value)) +
    geom_line(aes(color = type), linewidth = 0.6) +
    geom_point(aes(color = type), size = 0.8, alpha = 0.7)
  
  if (!is.null(df_plot_breaks)) {
    p <- p +
      geom_vline(
        data = df_plot_breaks,
        aes(xintercept = break_date, linetype = algorithm, color = algorithm),
        linewidth = 1,
        alpha = 0.6,
        show.legend = TRUE
      )
  }
  
  p <- p +
    scale_color_manual(
      name = "Time series",
      values = c(color_map_ts, color_map_algs),
      breaks = c(names(color_map_ts), names(color_map_algs)),
      # values = c(color_map_ts),
      #breaks = c(names(color_map_ts))
    ) +
    scale_linetype_manual(
      name = "Break detection algorithm",
      #values = rep("dotted", length(color_map_algs)),
      values = line_types,
      breaks = names(color_map_algs)
    ) +
    labs(
      title = paste("Time Series of", get_spi(spidf_obj)),
      subtitle = paste0(
        "Location: [", round(get_longitude(spidf_obj), 5), ", ",
        round(get_latitude(spidf_obj), 5), "] | ",
        "Period: ", format(get_range_start(spidf_obj)), " to ",
        format(get_range_end(spidf_obj)), " | ",
        "Proc. Level: ", get_proc_level(spidf_obj),
        ifelse(tree_id != "" && tree_id != "---", paste0(" | Tree ID: ", tree_id), "")
      ),
      x = "Date",
      y = get_spi(spidf_obj)
    ) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(0.2, 'cm'),
      legend.title = element_text(face = "bold")
    ) +
    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    )
  
  return(p)
}



