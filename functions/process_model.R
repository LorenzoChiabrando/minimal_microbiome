
# Define function to process one model
process_model <- function(model) {
  # Extract model information
  organism <- model$organism
  abbr <- model$abbreviation[2]
  FBAmodel <- ifelse(is.null(model$FBAmodel), organism, model$FBAmodel)
  
  # Log processing step
  cat(sprintf("\nProcessing %s (%s)\n", organism, abbr))
  cat(paste(rep("-", nchar(organism) + nchar(abbr) + 13), collapse = ""), "\n")
  
  # Construct file paths
  input_dir <- file.path(wd, "input", FBAmodel)
  mat_file <- file.path(input_dir, paste0(FBAmodel, ".mat"))
  
  # Validate input file
  if (!file.exists(mat_file)) {
    cat(sprintf("ERROR: Input file not found: %s\n", mat_file))
    return(list(
      status = "error",
      message = sprintf("Input file not found: %s", mat_file),
      organism = organism,
      abbr = abbr
    ))
  }
  
  # Generate FBA model
  cat(sprintf("Generating FBA model for %s...\n", organism))
  model_obj <- tryCatch({
    FBA4Greatmod.generation(fba_mat = mat_file, input_dir = input_dir)
  }, error = function(e) {
    cat(sprintf("ERROR: Failed to generate model: %s\n", e$message))
    return(NULL)
  })
  
  # If model generation failed, return error
  if ( is.null(model_obj) ) {
    return(list(
      status = "error",
      message = "Failed to generate model",
      organism = organism,
      abbr = abbr
    ))
  }
  
  # Set biomass parameters
  cat("Setting biomass parameters...\n")
  model_obj <- setBiomassParameters(
    model_obj, 
    bioMax = model$biomass$max, 
    bioMean = model$biomass$mean, 
    bioMin = model$biomass$min
  )
  
  # Save model in compiled_models directory
  output_dir <- file.path(wd, "compiled_models")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_file <- file.path(output_dir, paste0(abbr, "_model.txt"))
  capture_output <- capture.output(model_obj)
  writeLines(capture_output, output_file)
  
  cat(sprintf("Model saved to %s\n", output_file))
  
  return(list(
    status = "success",
    message = sprintf("Successfully processed %s (%s)", organism, abbr),
    organism = organism,
    abbr = abbr,
    model_file = output_file
  ))
}
