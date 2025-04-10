

library(shiny)
library(shinyjs)
library(shinyalert)

# Ensure cache directory exists
if (!dir.exists("ltm_cache")) dir.create("ltm_cache")

source("LTM_cache.R")

gv_list <- read_csv("GV_2025_list_v20250331-v1.csv")

gv_ids <- c("---",sort(unique(gv_list$cid)))

# UI
ui <- fluidPage(
  
  useShinyjs(),  # <--- activate JS
  
  #useShinyalert(),
  
  tags$head(
    tags$style(HTML("
.resizable-sidebar {
      resize: horizontal;
      overflow: auto;
      min-width: 200px;
      max-width: 350px;
      padding-right: 10px;
      border-right: 1px solid #ccc;
      background-color: #4a5e68; /* Muted slate blue */
      color: #1b4040;
    }

    .main-panel {
      padding-left: 20px;
    }

    .form-group,
    .shiny-input-container {
      margin-bottom: 8px;
      font-size: 13px;
    }

    .shiny-input-label,
    .control-label,
    .checkbox label,
    .radio label {
      color: #1b4040 !important;  /* Soft white */
      font-weight: 500;
    }

    .form-control,
    .selectize-input {
      background-color: #e0ebeb !important; /* Pale blue-gray */
      color: #2c3e50 !important;            /* Deep charcoal-blue text */
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
      background-color: #f4f6f6; /* Light dropdown */
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
        background-color: #354b5e;  /* dark blue aligned with theme */
        color: #ffffff;             /* white text */
        padding: 10px 20px;
        z-index: 1000;
        border-bottom: 1px solid #ccc;
      }
      
      .ltm-banner h2 {
        margin: 3px;
        padding: 3px;
        font-size: 25px;
        color: #ffffff;  /* ensure title text is white */
      }
      
      body {
        padding-top: 70px;  /* spacing for fixed header */
      }
    "))
  ),
  
  #titlePanel("LTM | Large Tree Monitoring"),
  div(class = "ltm-banner", titlePanel("LTM | Large Tree Monitoring")),
  
  
  column(
    width = 3,
    div(
      class = "resizable-sidebar",
      sidebarPanel(
        width = 12,  # take full width within resizable container
        textInput("user_name", "GEE User Name", value = ""),
        
        selectInput("location_id", "Select Location ID", choices = gv_ids, selectize = TRUE, selected = "---"#,
                    #options = list(maxOptions = 1000, placeholder = 'Start typing an ID...')
        ),
        
        numericInput("latitude", "Latitude", value = 41.720898),
        numericInput("longitude", "Longitude", value = -8.747039),
        
        dateInput("start_date", "Start Date", value = "2015-01-01"),
        dateInput("end_date", "End Date", value = "2024-12-31"),
        selectInput("spi", "Spectral Index", choices = c("NDVI", "EVI")),
        selectInput("proc_level", "Processing Level", choices = c("L1C","L2A" )),
        actionButton("fetch_data", "Fetch Data"),
        hr(),
        selectInput("regularize_method", "Regularization Method", 
                    choices = c("linear", "spline", "stine", "mean", "kalman", "locf", "nocb")),
        checkboxInput("use_cloud_mask", "Use Cloud Mask", value = TRUE),
        actionButton("regularize", "Regularize Time Series"),
        hr(),
        selectInput("quant", "Moving Window Quantile", choices = c(0.95, 0.99, 0.9, 0.75)),
        numericInput("win_size", "Window Size (days)", value = 15, min = 5, max = 180),
        actionButton("apply_mov_quantile", "Apply Moving Window Quantile"),
        hr(),
        numericInput("lambda", "Whittaker Lambda", value = 20000, min = 1),
        numericInput("quantile_thresh", "Whittaker Quantile Threshold", value = 0.35, min = 0.01, max = 1),
        checkboxInput("use_weights", "Use Weights in Whittaker Smoothing", value = TRUE),
        actionButton("apply_whittaker", "Apply Whittaker Smoother")
      )
    )
  ),
  column(
    width = 9,
    div(
      class = "main-panel",
      uiOutput("loading"),
      verbatimTextOutput("status"),
      plotOutput("plot_spidf", height = "600px")
    )
  )
  
  
)

# Server
server <- function(input, output, session) {
  
  spidf_data <- reactiveVal(NULL)
  
  
  observeEvent(input$fetch_data, {
    
    #output$status <- renderText("Checking local cache. Please wait...")
    
    output$loading <- renderUI({
      tagList(
        tags$img(src = "loading.gif", height = "90px"),
        tags$span("Checking local cache or fetching GEE data... Please wait")
      )
    })
    
    # Let UI render before heavy computation
    shinyjs::delay(1500, {
      if (ltm_check_gee_status() != "CONNECTED") {
        output$status <- renderText("Initializing GEE connection...")
        ltm_start_gee(input$user_name)
      }
      
      ################################################
      date_error <- FALSE
      
      # Issue error message for dates
      
      if((input$proc_level %in% c("L1","L1C")) && (input$start_date < as.Date("2015-01-01"))){
        
        shinyalert(
          title = "Error in start date",
          text = "The start date of the analysis cannot be earlier than 2015-01-01 for L1/L1C Sentinel-2 data",
          type = "error"
        )
        output$loading <- renderUI({ NULL })
        
        date_error <- TRUE
        #return(invisible(NULL)) 
      }
      if((input$proc_level %in% c("L2","L2A")) && (input$start_date < as.Date("2017-01-01"))){
        
        shinyalert(
          title = "Error in start date",
          text = "The start date of the analysis cannot be earlier than 2017-01-01 for L2/L2A Sentinel-2 data",
          type = "error"
        )
        output$loading <- renderUI({ NULL })
        
        #return(invisible(NULL))
        date_error <- TRUE
      }
      
      ################################################
      
      
      if(!date_error){
        
        
        # Get coordinates from selected location or fallback to manual inputs
        if (input$location_id != "---") {
          
          # Get specific row and its coordinates
          loc_row <- gv_list[gv_list$cid == input$location_id, ]
          
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
          #input$latitude, input$longitude,
          lat, lon,
          input$start_date, input$end_date,
          input$spi, input$proc_level
        )
        
        if (!is.null(cached_file)) {
          
          
          output$status <- renderText(paste("Data loaded from cache:", cached_file))
          spidf_data(readRDS(cached_file))
          output$loading <- renderUI({NULL})
          
        } else {
          #output$status <- renderText("No cache found. Fetching data from GEE...")
          shinyjs::delay(1500, {
            df <- ltm_s2_get_data_point(
              lat        = lat, # input$latitude,
              lon        = lon, # input$longitude,
              start_date = format(input$start_date, "%Y-%m-%d"),
              end_date   = format(input$end_date,   "%Y-%m-%d"),
              spi        = input$spi,
              proc_level = input$proc_level
            )
            
            file_name <- ltm_cache_file_name(
              #input$latitude, input$longitude,
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
  
  # Regularize
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
  
  # Moving window quantile
  observeEvent(input$apply_mov_quantile, {
    
    req(spidf_data())
    
    
    if (!isTRUE(is_regularized(spidf_data()))) {
      shinyalert(
        title = "Action Required",
        text = "Please regularize the time series before applying a moving window quantile.",
        type = "warning"
      )
      #return(invisible(NULL)) 
    }
    
    spidf_data(
      ltm_apply_moving_quantile(
        spidf_data(),
        quant   = as.numeric(input$quant),
        win_size = input$win_size
      )
    )
    output$status <- renderText("Moving window quantile applied.")
  })
  
  # Whittaker smoother
  observeEvent(input$apply_whittaker, {
    req(spidf_data())
    spidf_data(
      ltm_apply_whitaker(
        spidf_data(),
        lambda         = input$lambda,
        quantile_thresh = input$quantile_thresh,
        use_weights     = input$use_weights
      )
    )
    output$status <- renderText("Whittaker smoother applied.")
  })
  
  # Plot
  output$plot_spidf <- renderPlot({
    req(spidf_data())
    ltm_plot_spidf_ts(spidf_data(),tree_id = input$location_id)
  })
}

shinyApp(ui = ui, server = server)
