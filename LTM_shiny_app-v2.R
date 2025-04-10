

library(shiny)
library(shinyjs)
library(shinyalert)
library(readr)
library(ggplot2)

source("LTM_gee_auth.R")
source("LTM_gee_data.R")
source("LTM_pre_proc.R")
source("LTM_plots.R")
source("LTM_aux_funs.R")
source("LTM_break_detection.R")
source("LTM_breaks_class.R")
source("LTM_cache.R")



# Ensure cache directory exists
if (!dir.exists("ltm_cache")) dir.create("ltm_cache")

tree_ids_col <- "cid"
tree_list <- read_csv("GV_2025_list_v20250331-v1.csv")
tree_ids  <- c("---", sort(unique(tree_list[,tree_ids_col] %>% pull)))



ui <- fluidPage(
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("
      .resizable-sidebar {
        resize: horizontal;
        overflow: auto;
        min-width: 200px;
        max-width: 450px;
        padding-right: 10px;
        border-right: 1px solid #ccc;
        background-color: #4a5e68;
        color: #1b4040;
      }
      .main-panel {
        padding-left: 20px;
      }
      .form-group, .shiny-input-container {
        margin-bottom: 8px;
        font-size: 13px;
      }
      .shiny-input-label,
      .control-label,
      .checkbox label,
      .radio label {
        color: #1b4040 !important;
        font-weight: 500;
      }
      .form-control, .selectize-input {
        background-color: #e0ebeb !important;
        color: #2c3e50 !important;
        border: 1px solid #95a5a6;
        font-size: 13px;
        height: 32px;
      }
      .form-control:focus,
      .selectize-input:focus {
        border-color: #6ba292;
        box-shadow: none;
      }
      .selectize-dropdown-content {
        background-color: #f4f6f6;
        color: #2c3e50;
      }
      .btn {
        font-size: 13px;
        padding: 6px 10px;
        background-color: #6ba292;
        color: #ffffff;
        border: none;
        font-weight: bold;
        margin-top: 4px;
      }
      .btn:hover {
        background-color: #5a8d7c;
      }
      .ltm-banner {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        background-color: #354b5e;
        color: #ffffff;
        padding: 10px 20px;
        z-index: 1000;
        border-bottom: 1px solid #ccc;
      }
      .ltm-banner h2 {
        margin: 3px;
        padding: 3px;
        font-size: 25px;
        color: #ffffff;
      }
      body {
        padding-top: 70px;
      }
      
      .scrollable-sidebar {
        position: fixed;
        top: 70px; /* height of your banner */
        bottom: 0;
        left: 0;
        width: 350px;
        overflow-y: auto;
        overflow-x: hidden;
        z-index: 999;
      }
      
      .main-panel {
        margin-left: 330px; /* leave space for fixed sidebar */
        padding: 20px;
      }
      
    "))
  ),
  
  div(class = "ltm-banner", titlePanel("LTM | Large Tree Monitoring")),
  
  # column(
  #   width = 3,
  #   div(
  #     class = "resizable-sidebar",
  #     sidebarPanel(
  fluidRow(
    div(
      class = "resizable-sidebar scrollable-sidebar",
      sidebarPanel(
        width = 12,  # use full width inside the resizable container
        
        # Google Earth Engine input
        textInput("user_name", "GEE User Name", value = ""),
        
        # Location or Tree ID from pre-existing list
        selectInput("location_id", "Select Location ID", 
                    choices = tree_ids, selectize = TRUE, selected = "---"),
        
        numericInput("latitude", "Latitude", value = 41.720898),
        numericInput("longitude", "Longitude", value = -8.747039),
        
        dateInput("start_date", "Start Date", value = "2015-01-01"),
        dateInput("end_date", "End Date", value = "2024-12-31"),
        selectInput("spi", "Spectral Index", choices = c("NDVI", "EVI")),
        selectInput("proc_level", "Processing Level", choices = c("L1C","L2A")),
        actionButton("fetch_data", "Fetch Data"),
        
        # Data regularization
        hr(),
        selectInput("regularize_method", "Regularization Method", 
                    choices = c("linear", "spline", "stine", 
                                "mean", "kalman", "locf", "nocb")),
        checkboxInput("use_cloud_mask", "Use Cloud Mask", value = TRUE),
        actionButton("regularize", "Regularize Time Series"),
        
        # Moving window analysis
        hr(),
        
        selectInput("quant", "Moving Window Quantile", 
                    choices = c(0.95, 0.99, 0.9, 0.75)),
        numericInput("win_size", "Window Size (days)", 
                     value = 15, min = 5, max = 180),
        actionButton("apply_mov_quantile", "Apply Moving Window Quantile"),
        
        # Whittaker time series smoothing
        hr(),
        numericInput("lambda", "Whittaker Lambda", 
                     value = 20000, min = 1),
        numericInput("quantile_thresh", "Whittaker Quantile Threshold", 
                     value = 0.35, min = 0.01, max = 1),
        checkboxInput("use_weights", "Use Weights in Whittaker Smoothing", value = TRUE),
        actionButton("apply_whittaker", "Apply Whittaker Smoother"),
        
        
        # Break Detection component
        # -- NEW: Break Detection UI inputs --
        hr(),
        checkboxGroupInput(
          "ts_names", "Select Time Series Type(s)",
          choices = c(
            
            "Original time series (daily)" = "spi", 
            "Moving window series"         = "spi_mov_wind", 
            "Smoothed series"              = "spi_smooth", 
            "Moving window smoothed"       = "spi_mov_smooth"),
          selected = c("spi")  # Default selection
        ),
        
        numericInput("s_window", "Window size (STL decomposition)", 
                     value = 30, min = 5, max = 365),
        
        checkboxGroupInput(
          inputId = "break_methods", 
          label   = "Break-Detection Methods",
          choices = c(
            "Sequential Change Point Model" = "cpm",
            "Energy Divisive"               = "ed",
            "Bfast"                         = "bfast",
            "Bayesian Change Points Model"  = "mcp"
          ),
          selected = NULL
        ),
        
        dateInput("break_thresh_date", "Break Threshold Date", 
                  value = "2016-01-01"),
        numericInput("thresh_change", "% Change threshold", 
                     value = -10, max = 0),
        
        radioButtons(
          inputId = "thresh_fun", 
          label   = "Threshold change function",
          choices = c("Median","90% percentile"),
          selected = "Median"),
        
        actionButton("analyze_breaks", "Analyze Breaks")
      )
    )
  ),
  # column(
  #   width = 9,
  #   div(
  #     class = "main-panel",
  #     uiOutput("loading"),
  #     verbatimTextOutput("status"),
  #     plotOutput("plot_spidf", height = "500px"),
  #     hr(),
  #     #plotOutput("plot_breaks", height = "500px"),
  #     tableOutput("break_results_table"),
  #     uiOutput("break_summary_ui")
  #   )
  # )
  
  div(
    class = "main-panel",
    uiOutput("loading"),
    verbatimTextOutput("status"),
    plotOutput("plot_spidf", height = "500px"),
    hr(),
    #tableOutput("break_results_table"),
    
    div(
      style = "margin: 20px;",
      tableOutput("break_results_table")
    ),
    
    uiOutput("break_summary_ui")
  )
)

server <- function(input, output, session) {
  
  spidf_data <- reactiveVal(NULL)
  
  # ---- Data fetch logic ----
  observeEvent(input$fetch_data, {
    
    output$loading <- renderUI({
      tagList(
        tags$img(src = "loading.gif", height = "90px"),
        tags$span("Checking local cache or fetching GEE data server... Please wait")
      )
    })
    
    shinyjs::delay(1500, {
      if (ltm_check_gee_status() != "CONNECTED") {
        output$status <- renderText("Initializing GEE connection...")
        ltm_start_gee(input$user_name)
      }
      
      date_error <- FALSE
      # Validate dates for Sentinel L1C / L2A
      if ((input$proc_level %in% c("L1","L1C")) && (input$start_date < as.Date("2015-01-01"))) {
        shinyalert(
          title = "Error in start date",
          text = "The start date cannot be earlier than 2015-01-01 for L1/L1C Sentinel-2 data",
          type = "error"
        )
        output$loading <- renderUI({ NULL })
        date_error <- TRUE
      }
      if ((input$proc_level %in% c("L2","L2A")) && (input$start_date < as.Date("2017-01-01"))) {
        shinyalert(
          title = "Error in start date",
          text = "The start date cannot be earlier than 2017-01-01 for L2/L2A Sentinel-2 data",
          type = "error"
        )
        output$loading <- renderUI({ NULL })
        date_error <- TRUE
      }
      
      if(!date_error) {
        # Use either chosen location or manual lat/lon
        if (input$location_id != "---") {
          loc_row <- tree_list[tree_list$cid == input$location_id, ]
          if (nrow(loc_row) == 1) {
            lat <- loc_row$lat
            lon <- loc_row$lon
          } else {
            showNotification("Selected location ID not found in list.", type = "error")
            return(NULL)
          }
        } else {
          lat <- input$latitude
          lon <- input$longitude
        }
        
        cached_file <- ltm_check_cache(
          lat, lon,
          input$start_date, input$end_date,
          input$spi, input$proc_level
        )
        
        if (!is.null(cached_file)) {
          output$status <- renderText(paste("Data loaded from cache:", cached_file))
          spidf_data(readRDS(cached_file))
          output$loading <- renderUI({NULL})
          
        } else {
          shinyjs::delay(1500, {
            df <- ltm_s2_get_data_point(
              lat        = lat,
              lon        = lon,
              start_date = format(input$start_date, "%Y-%m-%d"),
              end_date   = format(input$end_date,   "%Y-%m-%d"),
              spi        = input$spi,
              proc_level = input$proc_level
            )
            
            file_name <- ltm_cache_file_name(
              lat, lon,
              input$start_date, input$end_date,
              input$spi, input$proc_level
            )
            
            saveRDS(df, file_name)
            output$status <- renderText(paste("Data fetched and cached at:", file_name))
            spidf_data(df)
            output$loading <- renderUI({NULL})
          })
        }
      }
    })
  })
  
  # ---- Regularize logic ----
  observeEvent(input$regularize, {
    req(spidf_data())
    spidf_data(
      ltm_regularize_spidf(
        spidf_data(),
        method         = input$regularize_method,
        use_cloud_mask = input$use_cloud_mask
      )
    )
    output$status <- renderText("Time series regularized.")
  })
  
  # ---- Moving window quantile ----
  observeEvent(input$apply_mov_quantile, {
    req(spidf_data())
    if (!isTRUE(is_regularized(spidf_data()))) {
      shinyalert(
        title = "Action Required",
        text  = "Please regularize the time series before applying a moving window quantile.",
        type  = "warning"
      )
      # still proceed or return invisible?
    }
    spidf_data(
      ltm_apply_moving_quantile(
        spidf_data(),
        quant    = as.numeric(input$quant),
        win_size = input$win_size
      )
    )
    output$status <- renderText("Moving window quantile applied.")
  })
  
  # ---- Whittaker smoother ----
  observeEvent(input$apply_whittaker, {
    req(spidf_data())
    spidf_data(
      ltm_apply_whitaker(
        spidf_data(),
        lambda          = input$lambda,
        quantile_thresh = input$quantile_thresh,
        use_weights     = input$use_weights
      )
    )
    output$status <- renderText("Whittaker smoother applied.")
  })
  
  # ---- Main plot (without break info) ----
  # output$plot_spidf <- renderPlot({
  #   req(spidf_data())
  #   ltm_plot_spidf_ts(spidf_data(), tree_id = input$location_id)
  # })
  
  # ---- Main plot ----
  output$plot_spidf <- renderPlot({
    req(spidf_data())
    
    # Retrieve the break results if available
    br_obj <- breakResults()   # This reactiveVal stores a list with df_breaks
    df_breaks <- NULL
    if (!is.null(br_obj)) {
      df_breaks <- br_obj$df_breaks
    }
    
    # Now call the updated plotting function with the new argument
    ltm_plot_spidf_ts(
      spidf_obj = spidf_data(),
      tree_id   = input$location_id,
      df_breaks = df_breaks    # <--- pass the breaks data if available
    )
  })
  
  
  
  # A reactive to hold the aggregated results (ts_breaks object + data frame)
  breakResults <- reactiveVal(NULL)
  
  # ---- Break Analysis Logic ----
  observeEvent(input$analyze_breaks, {
    
    
    output$loading <- renderUI({
      tagList(
        tags$img(src = "loading.gif", height = "90px"),
        tags$span(paste0("Performing break analysis. ",
                  "This may take a few seconds to some minutes... Please wait \n",
                  ifelse("mcp" %in% input$break_methods, 
                         "(Bayesian Change Points Model is very slow!)")

                  ))
      )
    })
    
    
    req(spidf_data())
    
    # If the user didn't select any methods or ts_name, just exit
    if (length(input$break_methods) == 0 || length(input$ts_names) == 0) {
      showNotification("No break methods or time-series types selected.", type = "warning")
      return(NULL)
    }
    
    ######
    ######
    
    shinyjs::delay(1000, {
      
      # Initialize a ts_breaks object
      tsb_obj <- ltm_ts_breaks(spidf_data())
      
      # Common parameters for all detection functions
      season_adj    <- TRUE
      s_window      <- input$s_window # 30
      thresh_change <- input$thresh_change #-10
      thresh_date   <- input$break_thresh_date
      tresh_int     <- NULL
      #thresh_fun    <- median
      
      if(input$thresh_fun == "Median"){
        thresh_fun <- median
      }else if(input$thresh_fun == "90% percentile"){
        thresh_fun <- Percentile90
      }
      
 
      # Run each selected method for each selected ts_name
      for (tsn in input$ts_names) {
        
        for (method in input$break_methods) {
          
          # We'll set a few method-specific defaults. Adjust if needed.
          if (method == "cpm") {
            # Example defaults
            run_obj <- ltm_cpm_detect_breaks(
              spidf       = spidf_data(),
              season_adj  = season_adj,
              ts_name     = tsn,
              s_window    = s_window,
              cpm_method  = "Exponential",
              ARL0        = 500,
              thresh_date = thresh_date,
              thresh_change = thresh_change,
              tresh_int     = tresh_int,
              thresh_fun    = thresh_fun
            )
            
          } else if (method == "ed") {
            # Example defaults
            run_obj <- ltm_ed_detect_breaks(
              spidf       = spidf_data(),
              ts_name     = tsn,
              season_adj  = season_adj,
              s_window    = s_window,
              sig.lvl     = 0.05,
              R           = 1000,
              k           = 1,
              min_size    = 30,
              alpha       = 1,
              thresh_date = thresh_date,
              thresh_change = thresh_change,
              tresh_int   = tresh_int,
              thresh_fun  = thresh_fun
            )
            
          } else if (method == "bfast") {
            # Example defaults
            run_obj <- ltm_bfast01_detect_breaks(
              spidf       = spidf_data(),
              ts_name     = tsn,
              formula     = response ~ harmon + trend,
              s_window    = s_window,
              test        = "OLS-MOSUM",
              level       = 0.05,
              aggregate   = all,
              trim        = NULL,
              bandwidth   = 0.15,
              functional  = "max",
              order       = 3,
              thresh_date = thresh_date,
              thresh_change = thresh_change,
              tresh_int   = tresh_int,
              thresh_fun  = thresh_fun
            )
            
          } else if (method == "mcp") {
            # Example defaults
            run_obj <- ltm_mcp_detect_breaks(
              spidf         = spidf_data(),
              ts_name       = tsn,
              season_adj    = season_adj,
              s_window      = s_window,
              thresh_change = thresh_change,
              thresh_date   = thresh_date,
              sample        = "both",
              n_chains      = 3,
              n_cores       = 3,
              n_adapt       = 1000,
              n_iter        = 2000,
              downsample    = 5
            )
          } else {
            next
          }
          
          # Add run to the ts_breaks object
          tsb_obj <- ltm_add_runs(tsb_obj, run_obj)
        }
      }
      
      # Convert to data frame and store
      df_breaks <- as.data.frame(tsb_obj)
      
      # Save them in a reactive
      breakResults(list(
        tsb_obj   = tsb_obj,
        df_breaks = df_breaks
      ))
      
      output$loading <- renderUI({ NULL })
      
      output$status <- renderText("Break analysis complete.")
      
    })
    
    ######
    ######
    
    
  })
  
  # ---- Table output with all break results ----
  output$break_results_table <- renderTable({
    br <- breakResults()
    req(br)
    
    df <- br$df_breaks
    
    # Format date columns to readable string format
    if ("break_date" %in% colnames(df)) {
      df$break_date <- format(df$break_date, "%Y-%m-%d")
    }
    
    # Remove columns
    df <- df %>% select(-method, -break_index)
    
    # Rename columns
    colnames(df) <- c(
      "Algorithm",
      "Run ID",
      "Data Type",
      "Breaks Found?",
      "Valid Breaks?",
      "Break Date",
      "Break Change"
    )
    
    df
  })
  
  
  output$break_summary_ui <- renderUI({
    
    br <- breakResults()
    req(br)
    
    txt <- ltm_summarize_break_df(br$df_breaks)
    
    div(
      style = "background-color: #f9f9f9; margin:20px; border: 1px solid #ccc; 
             padding: 15px; border-radius: 8px; font-size: 14px; 
             font-family: monospace; white-space: pre-wrap;",
      txt
    )
  })
  
}

shinyApp(ui = ui, server = server)
