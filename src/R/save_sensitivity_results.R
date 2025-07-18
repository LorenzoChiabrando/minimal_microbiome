
process_sensitivity_analysis <- function(wd, model_name = "BacPN", save_type = "both") {
  # Validate save_type parameter
  valid_save_types <- c("both", "trace", "flux", "none")
  if (!save_type %in% valid_save_types) {
    stop("save_type must be one of: 'both', 'trace', 'flux', or 'none'")
  }
  
  load(paste0(wd, "/BacPN_analysis/BacPN-analysis.RData"))
  
  # Helper function to extract number from filename
  extract_file_number <- function(filename, pattern) {
    as.numeric(gsub(pattern = pattern, replacement = "", x = basename(filename)))
  }
  
  # Helper function to load and validate config
  validate_config <- function(n, config) {
    if(!exists("config") || is.null(config[[2]][[n]][[3]][1])) {
      stop(paste("Invalid config structure for n =", n))
    }
    return(config[[2]][[n]][[3]][1])
  }
  
  # Helper function to save results based on type
  save_results <- function(results, type) {
    if (save_type %in% c("both", type)) {
      saveRDS(results, 
              paste0(wd, "/", model_name, "_analysis/sensitivity_", type, ".rds"))
      message(paste("Saved", type, "results successfully"))
    }
  }
  
  # Process trace files
  process_trace_files <- function(wd, config) {
    trace_files <- sample(list.files(pattern = "\\.trace$", 
                                     path = paste0(wd, "/", model_name, "_analysis"), 
                                     full.names = TRUE))
    
    lapply(trace_files, function(file) {
      n <- extract_file_number(file, 
                               pattern = paste0("(", model_name, "-analysis-)|(.trace)"))
      
      trace <- read.table(file, header = TRUE)
      
      # Process trace data
      sensitivity_trace <- trace %>% 
        dplyr::select(Time, E_coli, E_coli_biomass_e) %>%
        tidyr::gather(key = "Places", value = "Marking", -Time)
      
      # Add configuration information
      cbind(sensitivity_trace,
            config = rep(as.character(file), nrow(sensitivity_trace)),
            param_sen = rep(validate_config(n, config), nrow(sensitivity_trace)))
    })
  }
  
  # Process flux files
  process_flux_files <- function(wd, config) {
    flux_files <- sample(list.files(pattern = "\\.flux$", 
                                    path = paste0(wd, "/", model_name, "_analysis"), 
                                    full.names = TRUE))
    
    lapply(flux_files, function(file) {
      tryCatch({
        n <- extract_file_number(file, 
                                 pattern = paste("(", model_name, "-analysis-)|(-1-Escherichia_coli_str_K_12_substr_MG1655.flux)", sep = ""))
        
        # Read and process flux data
        flux <- read.table(file, header = FALSE)
        reactions_names <- gsub("\\(|\\)", "", gsub("\\(", "_", flux[1, -c(1, 2)]))
        colnames(flux) <- c("Time", "Obj_0", reactions_names)
        
        # Transform flux data
        temp_flux <- flux[, c("Time", reactions_names)] %>%
          tidyr::gather(key = "Reaction", value = "Flux", -Time)
        
        # Validate configuration and create result
        param_sen_value <- validate_config(n, config)
        
        cbind(temp_flux,
              config = rep(as.character(file), nrow(temp_flux)),
              param_sen = rep(param_sen_value, nrow(temp_flux)))
        
      }, error = function(e) {
        warning(paste("Error processing file", file, ":", e$message))
        return(NULL)
      })
    })
  }
  
  # Main execution
  tryCatch({
    results <- list()
    
    # Process trace files if needed
    if (save_type %in% c("both", "trace")) {
      sensitivity_trace <- process_trace_files(wd, config)
      sensitivity_trace <- Filter(Negate(is.null), sensitivity_trace)
      results$trace <- sensitivity_trace
      save_results(sensitivity_trace, "trace")
    }
    
    # Process flux files if needed
    if (save_type %in% c("both", "flux")) {
      sensitivity_flux <- process_flux_files(wd, config)
      sensitivity_flux <- Filter(Negate(is.null), sensitivity_flux)
      results$flux <- sensitivity_flux
      save_results(sensitivity_flux, "flux")
    }
    
    return(results)
    
  }, error = function(e) {
    stop(paste("Error in sensitivity analysis:", e$message))
  })
}