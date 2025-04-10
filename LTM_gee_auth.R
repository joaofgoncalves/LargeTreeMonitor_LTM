
ltm_start_gee <- function(user_name = NULL, ...){
  
  out <- try(ee_Initialize(user = user_name, ...), silent = TRUE)
  
  if(inherits(out, "try-error")){
    options(LTM_GEE_STATUS = "NOT_CONNECTED")
    return(FALSE)
  }
  else{
    options(LTM_GEE_STATUS = "CONNECTED")
    return(TRUE)
  }
}

ltm_check_gee_status <- function() {
  
  status <- getOption("LTM_GEE_STATUS")
  
  if (is.null(status)) {
    warning("GEE has not been initialized yet.")
    return(FALSE)
  }
  return(status)
}

