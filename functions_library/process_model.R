# Define function to process one model
process_model <- function(model, hypernode_name) {
  # Extract model information
  organism <- model$organism
  abbr <- model$abbreviation[2]
  FBAmodel <- ifelse(is.null(model$FBAmodel), organism, model$FBAmodel)
  
  # Log processing step
  cat(sprintf("\nProcessing %s (%s)\n", organism, abbr))
  cat(paste(rep("-", nchar(organism) + nchar(abbr) + 13), collapse = ""), "\n")
  
  # Correct paths
  mat_file <- file.path(wd, "metabolic_networks_library", paste0(FBAmodel, ".mat"))
  input_dir <- file.path(wd, "hypernodes", hypernode_name, paste0("metabolic_networks_", hypernode_name), FBAmodel)
  output_dir <- file.path(input_dir, paste0("compiled_", FBAmodel))
  
  # Validate .mat file
  if (!file.exists(mat_file)) {
    cat(sprintf("ERROR: Input file not found: %s\n", mat_file))
    return(list(
      status = "error",
      message = sprintf("Input file not found: %s", mat_file),
      organism = organism,
      abbr = abbr
    ))
  }
  
  # Create input and output directories if needed
  if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Copy the .mat file into the model-specific input_dir (if needed)
  file.copy(mat_file, file.path(input_dir, paste0(FBAmodel, ".mat")))
  
  # Generate FBA model
  cat(sprintf("Generating FBA model for %s...\n", organism))
  model_obj <- tryCatch({
    FBA4Greatmod.generation(
      fba_mat = file.path(input_dir, paste0(FBAmodel, ".mat")),
      input_dir = input_dir
    )
  }, error = function(e) {
    cat(sprintf("ERROR: Failed to generate model: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(model_obj)) {
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
  
  # Save model
  output_file <- file.path(output_dir, paste0(abbr, "_model.txt"))
  capture_output <- capture.output(model_obj)
  writeLines(capture_output, output_file)
  
  cat(sprintf("Model saved to %s\n", output_file))
  
  process_results <- list(
    status = "success",
    message = sprintf("Successfully processed %s (%s)", organism, abbr),
    organism = organism,
    abbr = abbr,
    model_file = output_file
  )
  
  # Proper printing
  status_icon <- if (process_results$status == "success") "✓" else "✗"
  msg         <- if (process_results$status == "success")
    "Successfully processed"
  else
    paste("ERROR -", process_results$message)
  
  cat(sprintf("%s  %s (%s): %s\n",
              status_icon,
              process_results$organism,
              process_results$abbr,
              msg))
  
  return(process_results)

}

