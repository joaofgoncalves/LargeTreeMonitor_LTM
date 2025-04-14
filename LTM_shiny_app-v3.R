
library(bcp)
library(bfast)
library(bitops)
library(cpm)
library(dplyr)
library(ecp)
library(forecast)
library(ggplot2)
library(imputeTS)
library(leaflet)
library(lubridate)
library(mcp)
library(pastecs)
library(phenofit)
library(pracma)
library(readr)
library(rgee)
library(rgeeExtra)
library(sf)
library(shiny)
library(shinyalert)
library(shinyjs)
library(slider)
library(terra)
library(tidyverse)
library(jsonlite)
library(patchwork)  # for side-by-side layouts


source("LTM_gee_auth.R")
source("LTM_gee_data.R")
source("LTM_pre_proc.R")
source("LTM_plots.R")
source("LTM_aux_funs.R")
source("LTM_break_detection.R")
source("LTM_breaks_class.R")
source("LTM_cache.R")
source("LTM_CheckConsts.R")


params <<- fromJSON("params.json")

# Ensure cache directory exists
if (!dir.exists("ltm_cache")) dir.create("ltm_cache")


tree_list_path <- params$tree_list_file$tree_locs_file

if(tree_list_path!="" && !is.na(tree_list_path) && 
   (length(tree_list_path)!=0) && !is.null(tree_list_path)){
  
  tree_ids_col <- params$tree_list_file$tree_ids_col
  tree_lat <- params$tree_list_file$lat_col
  tree_lon <- params$tree_list_file$lon_col
  
  tree_list <- read_csv(tree_list_path)
  tree_ids  <- c("---", sort(unique(tree_list[,tree_ids_col] %>% pull)))
  
}else{
  tree_ids  <- c("---")
}



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
  
  #div(class = "ltm-banner", titlePanel("🌳 LTM | Large Tree Monitoring")),
  
  div(
    class = "ltm-banner",
    div(style = "display: flex; align-items: center;",
        tags$img(src = "LTM_logo-v1.png", height = "60px", style = "margin-right: 5px;"),
        h2("Large Tree Monitoring")
    )
  ),
  
  
  fluidRow(
    div(
      class = "resizable-sidebar scrollable-sidebar",
      sidebarPanel(
        width = 12,  # use full width inside the resizable container
        
        # Google Earth Engine input
        #textInput("user_name", "GEE User Name", value = ""),
        
        # Location or Tree ID from pre-existing list
        selectInput("location_id", "Select Location ID", 
                    choices = tree_ids, selectize = TRUE, selected = "---"),
        
        numericInput("latitude", "Latitude", value = 41.720898),
        numericInput("longitude", "Longitude", value = -8.747039),
        
        actionButton("choose_coords", "", 
                     icon = icon("map-location-dot"), 
                     style = "margin-bottom: 8px;"),
        
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
        
        actionButton("regularize", "Regularize/interpolate Time Series"),
        
        # Moving window analysis
        hr(),
        
        selectInput("quant", "Moving Window Quantile", 
                    choices = c(0.95, 0.9, 0.75)),
        numericInput("win_size", "Window Size (days)", 
                     value = 15, min = 5, max = 180),
        
        actionButton("apply_mov_quantile", "Apply Moving Window"),
        
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
          # choices = c(
          #   "Sequential Change Point Model" = "cpm",
          #   "Energy Divisive"               = "ed",
          #   "Bfast"                         = "bfast",
          #   "Bayesian Change Points Model"  = "mcp",
          #   "Structural Changes"            = "stc",
          #   "Wild Binary Segmentation"      = "wbs"
          # ),
          choiceNames = list(
            HTML("Change Point Model (cpm) <span title='Detects sequential changes using the cpm package'> ℹ️</span>"),
            HTML("Energy Divisive (ed) <span title='ecp package e.divisive method: clustering-based method'> ℹ️</span>"),
            HTML("BFAST <span title='Breaks in season-trend using bfast'> ℹ️</span>"),
            HTML("Bayesian Change Points (mcp) <span title='Bayesian model-based detection using mcp'> ℹ️</span>"),
            HTML("Structural Changes (stc) <span title='strucchange package for structural breaks'> ℹ️</span>"),
            HTML("Wild Binary Segmentation (wbs) <span title='Fast segmentation of time series'> ℹ️</span>")
          ),
          choiceValues = c("cpm", "ed", "bfast", "mcp", "stc", "wbs"),
          
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
        
        actionButton("analyze_breaks", "Analyze Breaks"),
        
        ## Parameter editing
        hr(),
        p("⚙️ Configure parameters"),
        actionButton("edit_params", "Edit parameters"),
        verbatimTextOutput("param_save_status")
        
        
      )
    )
  ),
  
  ## Right-side with analysis results
  
  div(
    class = "main-panel",
    
    uiOutput("loading"),
    
    verbatimTextOutput("status"),
    
    # Time series plot + small thumbnail map side by side
    fluidRow(
      column(
        width = 9,
        plotOutput("plot_spidf", height = "500px")
      ),
      column(
        width = 3,
        leafletOutput("location_map", height = "500px")  # <-- Leaflet map output
      )
    ),
    
    uiOutput("download_data_ui"),  # <- dynamic button only when data is ready
    
    
    hr(),
    
    verbatimTextOutput("break_analysis_status"),
    
    div(
      style = "margin: 20px;",
      tableOutput("break_results_table"),
      
      # Add this download button next to the table
      #downloadButton("download_breaks_csv", "Download table to csv file")
      uiOutput("download_breaks_ui")  # <- dynamic button only when data is ready
      
    ),
    
    # Two new plots side by side (patchwork) - only if valid breaks exist
    fluidRow(
      column(
        width = 12,
        plotOutput("valid_breaks_plots", height = "400px")
      )
    ),
    
    uiOutput("break_summary_ui")
  )
)



server <- function(input, output, session) {
  
  spidf_data <- reactiveVal(NULL)
  
  # A new reactiveVal to store the actual lat/lon used in the last fetch
  selectedCoord <- reactiveVal(NULL)
  
  
  ## ------------------------------------------------------------------------ ##
  #  ---- Coordinate picker with reactive val to update main coordinates ----
  ## ------------------------------------------------------------------------ ##
  
  # Reactive value to store the coordinate chosen in the modal
  pickedCoords <- reactiveVal(list(lat = NULL, lon = NULL))
  
  observeEvent(input$choose_coords, {
    # Reset the picked coordinates
    pickedCoords(list(lat = NULL, lon = NULL))
    
    # Re-render the leaflet output for the modal to ensure no marker is present
    output$coord_picker_map <- renderLeaflet({
      lat <- input$latitude
      lon <- input$longitude
      leaflet() %>%
        addTiles(group = "OpenStreetMap") %>%
        addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
        setView(lng = lon, lat = lat, zoom = 11) %>%
        addLayersControl(
          baseGroups = c("OpenStreetMap", "Satellite"),
          options = layersControlOptions(collapsed = FALSE)
        )
    })
    
    showModal(modalDialog(
      title = "Select location",
      size = "xl",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_coords", "Confirm selection")
      ),
      leafletOutput("coord_picker_map", height = "600px")
    ))
  })
  
  output$coord_picker_map <- renderLeaflet({
    lat <- input$latitude
    lon <- input$longitude
    
    leaflet() %>%
      addTiles(group = "OpenStreetMap") %>%
      addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
      setView(lng = lon, lat = lat, zoom = 13) %>%
      addLayersControl(
        baseGroups = c("OpenStreetMap", "Satellite"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  observeEvent(input$coord_picker_map_click, {
    
    click <- input$coord_picker_map_click
    
    if (!is.null(click)) {
      
      pickedCoords(list(lat = click$lat, lon = click$lng))
      
      leafletProxy("coord_picker_map") %>%
        clearMarkers() %>%
        addMarkers(lng = click$lng, lat = click$lat, popup = "Chosen location")
    }
  })
  
  observeEvent(input$confirm_coords, {
    coords <- pickedCoords()
    
    if (!is.null(coords$lat) && !is.null(coords$lon)) {
      
      updateNumericInput(session, "latitude", value = coords$lat)
      updateNumericInput(session, "longitude", value = coords$lon)
      
      selectedCoord(list(lat = coords$lat, lon = coords$lon))
    }
    
    removeModal()
  })
  
  
  ## -------------------------------------------------------------- ##
  # ---- Data fetch logic ----
  ## -------------------------------------------------------------- ##
  
  observeEvent(input$fetch_data, {
    
    output$loading <- renderUI({
      tagList(
        tags$img(src = "loading.gif", height = "90px"),
        tags$span("🔄 Checking local cache or fetching GEE data server... Please wait")
      )
    })
    
    shinyjs::delay(1000, {
      if (ltm_check_gee_status() != "CONNECTED") {
        output$status <- renderText("🔄 Initializing GEE connection...")
        #ltm_start_gee(input$user_name)
        ltm_start_gee(params$gee$gee_username)
        
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
          
          #loc_row <- tree_list[tree_list$cid == input$location_id, ]
          loc_row <- tree_list[tree_list %>% pull(tree_ids_col) == input$location_id, ]
          
          if (nrow(loc_row) == 1) {
            # lat <- loc_row$lat
            # lon <- loc_row$lon
            lat <- loc_row %>% pull(tree_lat)
            lon <- loc_row %>% pull(tree_lon)
            
          } else {
            #showNotification("Selected location ID not found in list.", type = "error")
            
            shinyalert(
              title = "Error in start date",
              text = "Selected location ID not found in list",
              type = "error"
            )
            
            return(NULL)
          }
        } else {
          lat <- input$latitude
          lon <- input$longitude
        }
        
        # Store the lat/lon used (for the map, etc.)
        selectedCoord(list(lat = lat, lon = lon))
        
        
        ## Check if it is necessary to clean break analysis outputs if a 
        # new set of coordinates have appeared including the nulling the 
        ## breaks in the plot and the tables in the bottom part. 
        ###########################################################################
        
        last_loc <- ltm_read_location()
        current_loc <- paste(round(lat, 6), round(lon, 6), sep = ",")
        
        # If location changed, reset break results
        if (is.null(last_loc) || (!identical(current_loc, last_loc))) {
          # Clear break analysis objects & outputs
          breakResults(NULL)
          output$break_analysis_status <- renderText("")
        }
        
        # 3) Update the stored location
        ltm_store_location(lat, lon)
        
        ###########################################################################

        cached_file <- ltm_check_cache(
          lat, lon,
          input$start_date, input$end_date,
          input$spi, input$proc_level
        )
        
        if (!is.null(cached_file)) {
          output$status <- renderText(paste("✅ Data loaded from cache:", cached_file))
          spidf_data(readRDS(cached_file))
          output$loading <- renderUI({NULL})
          
        } else {
          
          shinyjs::delay(1000, {
            df <- ltm_s2_get_data_point(
              lat        = lat,
              lon        = lon,
              start_date = format(input$start_date, "%Y-%m-%d"),
              end_date   = format(input$end_date,   "%Y-%m-%d"),
              spi        = input$spi,
              proc_level = input$proc_level,
              tree_id = ifelse(input$location_id != "" && input$location_id != "---",
                               input$location_id, NULL)
            )
            
            file_name <- ltm_cache_file_name(
              lat, lon,
              input$start_date, input$end_date,
              input$spi, input$proc_level
            )
            
            saveRDS(df, file_name)
            
            output$status <- renderText(paste("✅ Data fetched and cached at:", file_name,"💾 "))
            
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
    output$status <- renderText("✅ Time series regularized.")
  })
  
  # ---- Moving window quantile ----
  observeEvent(input$apply_mov_quantile, {
    
    req(spidf_data())
    
    if (!isTRUE(is_regularized(spidf_data()))) {
      shinyalert(
        title = "Action Required",
        text  = "Please regularize/interpolate the time series before 
        applying the moving window analysis.",
        type  = "error"
      )
      return(NULL)
    }
    
    # Perform the moving window analysis with quantiles
    spidf_data(
      ltm_apply_moving_quantile(
        spidf_data(),
        quant    = as.numeric(input$quant),
        win_size = input$win_size
      )
    )
    
    output$status <- renderText("✅ Moving window quantile applied.")
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
    output$status <- renderText("✅ Whittaker smoother applied.")
  })
  
  # ---- Main plot ----
  output$plot_spidf <- renderPlot({
    req(spidf_data())
    
    # Retrieve the break results if available
    br_obj <- breakResults()   # This reactiveVal stores a list with df_breaks
    df_breaks <- NULL
    if (!is.null(br_obj)) {
      df_breaks <- br_obj$df_breaks
    }
    
    ltm_plot_spidf_ts(
      spidf_obj = spidf_data(),
      tree_id   = input$location_id,
      df_breaks = df_breaks
    )
  })
  
  
  ## Download button for time series data
  
  output$download_data_ui <- renderUI({
    
    df <- spidf_data()#breakResults()
    req(df)
    
    downloadButton("download_data_csv", "Download time series to CSV")
  })
  
  output$download_data_csv <- downloadHandler(
    filename = function() {
      paste0("TimeSeriesData_",format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      # Safely retrieve the breaks data frame
      df <- spidf_data()#breakResults()
      req(df)  # Ensure breakResults() is not NULL
      #df <- br$df_breaks
      
      # Write the CSV
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  

  # ---- Small Leaflet map showing the selected lat/lon ----
  output$location_map <- renderLeaflet({
    coords <- selectedCoord()
    req(coords)
    
    leaflet() %>%
      addTiles(group = "OpenStreetMap") %>%
      addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
      setView(lng = coords$lon, lat = coords$lat, zoom = 13) %>%
      addMarkers(lng = coords$lon, lat = coords$lat, popup = "Tree Location") %>%
      addLayersControl(
        baseGroups = c("OpenStreetMap", "Satellite"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  # A reactive to hold the aggregated results (ts_breaks object + data frame)
  breakResults <- reactiveVal(NULL)
  
  
  #################################################
  
  observeEvent(input$edit_params, {
    showModal(modalDialog(
      title = "Edit parameters (JSON)",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_params", "💾 Save Changes")
      ),
      textAreaInput(
        inputId = "json_editor",
        label = NULL,
        value = "",  # initially blank, will update below
        width = "100%",
        height = "500px"
      )
    ))
    
    shinyjs::delay(100, {
      updateTextAreaInput(
        session,
        inputId = "json_editor",
        value = as.character(toJSON(params, pretty = TRUE, auto_unbox = TRUE))
      )
    })
  })
  
  
  observeEvent(input$save_params, {
    tryCatch({
      new_params <- fromJSON(input$json_editor, simplifyVector = TRUE)
      
      # Save to file
      write(toJSON(new_params, pretty = TRUE, auto_unbox = TRUE),
            file = "params.json")
      
      # Update the reactive global params (if you want real-time usage)
      params <<- new_params
      
      removeModal()
      output$param_save_status <- renderText("✅ Parameters saved")
    }, error = function(e) {
      output$param_save_status <- renderText(
        paste("❌ Error parsing JSON:", e$message)
      )
    })
  })
  
  #################################################
  
  
  # ---- Break Analysis Logic ----
  observeEvent(input$analyze_breaks, {
    
    message("[DEBUG] Analyze_breaks triggered.")
    message("[DEBUG] Break_methods: ", paste(input$break_methods, collapse=", "))
    message("[DEBUG] ts_names: ", paste(input$ts_names, collapse=", "))
    
    req(spidf_data())
    message("[DEBUG] spidf_data rows: ", nrow(spidf_data()))  # or however your spidf is structured
    
    if (length(input$break_methods) == 0 || length(input$ts_names) == 0) {
      #showNotification("No break methods or time-series types selected.", type = "error")
      #message("[DEBUG] => no break methods or no timeseries, returning NULL.")
      
      shinyalert(
        title = "Action Required",
        text  = "No break methods or time-series types selected.",
        type  = "error"
      )
      
      return(NULL)
    }
    
    if(!all(c("spi_smooth", "spi_mov_smooth") %in% colnames(spidf_data()))){
      
      #showNotification("Time series smoothing has to be performed first", type = "error")
      
      shinyalert(
        title = "Action Required",
        text  = "Time series smoothing (and maybe moving window analysis) have to be performed first",
        type  = "error"
      )
      
      return(NULL)
      
    }
    
    output$loading <- renderUI({
      tagList(
        tags$img(src = "loading.gif", height = "90px"),
        tags$span(paste0("Performing break analysis. ",
                         "This may take a few seconds to some minutes... Please wait \n",
                         ifelse("mcp" %in% input$break_methods, 
                                "(Bayesian Change Points Model is very slow!)", "")
        ))
      )
    })
    
    req(spidf_data())
    
    # If the user didn't select any methods or ts_name, just exit
    if (length(input$break_methods) == 0 || length(input$ts_names) == 0) {
      #showNotification("No break methods or time-series types selected.", type = "warning")
      
      shinyalert(
        title = "Action Required",
        text  = "No break methods or time-series types were selected",
        type  = "error"
      )
      
      return(NULL)
    }
    
    shinyjs::delay(1000, {
      
      tsb_obj <- ltm_ts_breaks(spidf_data())
      
      # Common parameters for all detection functions
      season_adj    <- TRUE
      s_window      <- input$s_window
      thresh_change <- input$thresh_change
      thresh_date   <- input$break_thresh_date
      tresh_int     <- NULL
      
      if(input$thresh_fun == "Median"){
        thresh_fun <- median
      } else if(input$thresh_fun == "90% percentile"){
        thresh_fun <- Percentile90
      }
      
      ##
      ## Loop through all time series and algorithms for break detection selected
      ##
      ##
      for (tsn in input$ts_names) {
        for (method in input$break_methods) {
          
          if (method == "cpm") {
            run_obj <- ltm_cpm_detect_breaks(
              spidf         = spidf_data(),
              season_adj    = season_adj,
              ts_name       = tsn,
              s_window      = s_window,
              cpm_method    = params$cpm$cpm_method,
              ARL0          = params$cpm$ARL0,
              thresh_date   = thresh_date,
              thresh_change = thresh_change,
              tresh_int     = tresh_int,
              thresh_fun    = thresh_fun
            )
            
          } else if (method == "ed") {
            run_obj <- ltm_ed_detect_breaks(
              spidf         = spidf_data(),
              ts_name       = tsn,
              season_adj    = season_adj,
              s_window      = s_window,
              sig.lvl       = params$ed$sig.lvl,
              R             = params$ed$R,
              k             = params$ed$k,
              min_size      = params$ed$min_size,
              alpha         = params$ed$alpha,
              thresh_date   = thresh_date,
              thresh_change = thresh_change,
              tresh_int     = tresh_int,
              thresh_fun    = thresh_fun
            )
            
          } else if (method == "bfast") {
            run_obj <- ltm_bfast01_detect_breaks(
              spidf         = spidf_data(),
              ts_name       = tsn,
              formula       = as.formula(params$bfast$formula),
              s_window      = s_window,
              test          = params$bfast$test,
              level         = params$bfast$level,
              aggregate     = all,
              trim          = NULL,
              bandwidth     = params$bfast$bandwidth,
              functional    = params$bfast$functional,
              order         = params$bfast$order,
              thresh_date   = thresh_date,
              thresh_change = thresh_change,
              tresh_int     = tresh_int,
              thresh_fun    = thresh_fun
            )
 
            
          } else if (method == "mcp") {
            run_obj <- ltm_mcp_detect_breaks(
              spidf         = spidf_data(),
              ts_name       = tsn,
              season_adj    = season_adj,
              s_window      = s_window,
              thresh_change = thresh_change,
              thresh_date   = thresh_date,
              sample        = params$mcp$sample,
              n_chains      = params$mcp$n_chains,
              n_cores       = params$mcp$n_cores,
              n_adapt       = params$mcp$n_adapt,
              n_iter        = params$mcp$n_iter,
              downsample    = params$mcp$downsample
            )
            
          } else if (method == "stc") {
            run_obj <- ltm_strucchange_detect_breaks(
              spidf         = spidf_data(),
              ts_name       = tsn,
              season_adj    = season_adj,
              s_window      = s_window,
              h             = params$stc$h,
              breaks        = params$stc$breaks,
              thresh_date   = thresh_date,
              thresh_change = thresh_change,
              tresh_int     = tresh_int,
              thresh_fun    = thresh_fun
            )
            
          } else if (method == "wbs") {
            run_obj <- ltm_wbs_detect_breaks(
              spidf         = spidf_data(),
              ts_name       = tsn,
              season_adj    = season_adj,
              s_window      = s_window,
              thresh_date   = thresh_date,
              thresh_change = thresh_change,
              tresh_int     = tresh_int,
              thresh_fun    = thresh_fun,
              num_intervals = params$wbs$num_intervals
            )
          }

          else {
            next
          }
          
          tsb_obj <- ltm_add_runs(tsb_obj, run_obj)
        }
      }
      
      # Convert the ts_breaks objects into a more readable data frame
      df_breaks <- as.data.frame(tsb_obj)
      message("[DEBUG] Analysis done. Setting breakResults.")
      
      
      # Set the break analysis reactive value
      breakResults(list(
        tsb_obj   = tsb_obj,
        df_breaks = df_breaks
      ))
      
      # Remove the ongoing analysis panel
      output$loading <- renderUI({ NULL })
      
      # Make an output appear once finishing the analysis
      output$break_analysis_status <- renderText("✅ Break analysis complete.")
    })
    
  })
  
  # --------------------------------------------- #
  # ---- Table output with all break results ----
  # --------------------------------------------- #
  
  output$break_results_table <- renderTable({
    br <- breakResults()
    req(br)
    
    df <- br$df_breaks
    
    # Format as dates to avoid conversion to integers
    if ("break_date" %in% colnames(df)) {
      df$break_date <- format(df$break_date, "%Y-%m-%d")
    }
    
    # Remove columns
    df <- df %>% dplyr::select(-method, -break_index)
    
    # Replace values for easier reading
    # Replace time series data names
    replacements <- 
      c("spi" = "Original time series (daily)",
        "spi_mov_wind" = "Moving window series",     
        "spi_smooth" = "Smoothed series",            
        "spi_mov_smooth" = "Moving window smoothed" )
    df$data_type <- recode(df$data_type, !!!replacements)
    
    # Replace model names
    da_replacements <- 
      c(
        "cpm"   = "Change Point Model (cpm)",
        "ed"    = "Energy Divisive (ed)",
        "bfast" = "Break Detection Seasonal & Trend (bfast)",
        "mcp"   = "Bayesian Multiple Change Points Model (mcp)",
        "stc"   = "Structural Changes (stc)",
        "wbs"   = "Wild Binary Segmentation (wbs)"
      )
    df$algorithm <- recode(df$algorithm, !!!da_replacements)
    
    # Replace boolean values
    df$has_breaks <- ifelse(df$has_breaks, "Yes", "No")
    df$has_valid_breaks <- ifelse(df$has_valid_breaks, "Yes", "No")
    
    # Rename columns for easier reading
    colnames(df) <- c(
      "Algorithm",
      "Run ID",
      "Data type",
      "Breaks found?",
      "Valid breaks?",
      "Break date",
      "Break change"
    )
    
    df
  })
  
  ## Download button
  
  output$download_breaks_ui <- renderUI({
    br <- breakResults()
    req(br)
    
    downloadButton("download_breaks_csv", "Download breaks to CSV")
  })
  
  output$download_breaks_csv <- downloadHandler(
    filename = function() {
      paste0("BreakResults_",format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      # Safely retrieve the breaks data frame
      br <- breakResults()
      req(br)  # Ensure breakResults() is not NULL
      df <- br$df_breaks
      
      # Write the CSV
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  
  output$valid_breaks_plots <- renderPlot({
    # 1) Access the results
    br <- breakResults()
    req(br)  # makes sure breakResults() is not NULL
    df <- br$df_breaks
    
    # 2) Call the function
    
    if(any(df$has_valid_breaks == TRUE)){
      plot_valid_breaks(df)
    }else{
      return(NULL)
    }
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
