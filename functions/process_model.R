
# Define function to process one model
process_model <- function(model) {
  # Extract model information
  organism <- model$organism
  abbr <- model$abbreviation
  FBAmodel <- ifelse(is.null(model$FBAmodel), organism, model$FBAmodel)
  
  # Log processing step
  cat(sprintf("\nProcessing %s (%s)\n", organism, abbr))
  cat(paste(rep("-", nchar(organism) + nchar(abbr) + 13), collapse = ""), "\n")
  
  # Construct file paths
  input_dir <- file.path(wd, "input", organism)
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
    FBA4Greatmod.generation(fba_mat = mat_file)
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
  
  # Create output directories if they don't exist
  reactions_dir <- file.path(input_dir)
  if (!dir.exists(reactions_dir)) dir.create(reactions_dir, recursive = TRUE)
  
  # Write metadata files
  cat("Writing metadata files...\n")
  reactions_file <- file.path(reactions_dir, "reactions_metadata.csv")
  write.csv(reactions_df, file = reactions_file, row.names = FALSE)
  
  metabolites_file <- file.path(reactions_dir, "metabolites_metadata.csv")
  write.csv(met_metadata, file = metabolites_file, row.names = FALSE)
  
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
    reactions_file = reactions_file,
    metabolites_file = metabolites_file,
    model_file = output_file
  ))
}
