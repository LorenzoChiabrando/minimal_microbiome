
#' Initialize Working Environment
#'
#' Sets up the working environment by sourcing required functions and creating directories
#'
#' @param wd Character string specifying the working directory path
#' @return Character string containing the destination directory path
#' @export
initialize_environment <- function(wd) {
  # Source required functions
  required_functions <- c(
    "FBAgreatmodeClass.R",
    "class_generation.R",
    "readMat.R"
  )
  
  for (func in required_functions) {
    source_path <- file.path(wd, "epimod_FBAfunctions/R", func)
    if (!file.exists(source_path)) {
      stop(paste("Required function file not found:", source_path))
    }
    source(source_path)
  }
  
  # Create results directory
  results_dir <- file.path(wd, "results/")
  dir.create(results_dir, showWarnings = FALSE)
  
  return(results_dir)
}

#' Process FBA Model
#'
#' Processes a single FBA model by generating the model, setting parameters,
#' and writing output files
#'
#' @param model_tags List containing model information (name and type)
#' @param wd Character string specifying the working directory path
#' @param dest_dir Character string specifying the destination directory path
#' @param bio_params List containing biomass parameters (max, mean, min)
#' @return Numeric value indicating processing time in seconds
#' @export
process_model <- function(model_tags, 
                         wd, 
                         dest_dir, 
                         bio_params = FBA_PARAMETERS$biomass) {
  
  # Extract model information
  model_name <- model_tags$name
  model_type <- model_tags$type
  
  # Build and verify input path
  matfile <- file.path(wd, 'input', model_type, paste0(model_name, ".mat"))
  if (!file.exists(matfile)) {
    stop(paste("Input file not found:", matfile))
  }
  
  # Generate and configure model
  message("Generating model for: ", model_name)
  model <- FBA4Greatmod.generation(fba_mat = matfile)
  model <- setBiomassParameters(
    model, 
    bio_params$max, 
    bio_params$mean, 
    bio_params$min
  )
  
  # Handle output files
  output_rds_path <- file.path(wd, "input/models", model_type)
  move_rds_files(output_rds_path)
  
  # Write FBA file and measure time
  start_time <- Sys.time()
  writeFBAfile(model, model_name, dest_dir)
  end_time <- Sys.time()
  
  return(as.numeric(difftime(end_time, start_time, units = "secs")))
}

#' Move RDS Files to Destination
#'
#' Helper function to move generated RDS files to their final location
#'
#' @param output_path Character string specifying the destination path for RDS files
#' @return Invisible NULL
#' @export
move_rds_files <- function(output_path) {
  files_to_move <- list.files(pattern = "\\.rds$", full.names = TRUE)
  if (length(files_to_move) > 0) {
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    sapply(files_to_move, function(f) {
      file.rename(from = f, to = file.path(output_path, basename(f)))
    })
  }
  invisible(NULL)
}

#' Process FBA Models
#'
#' Main function to process all FBA models in batch
#'
#' @param working_dir Optional character string specifying the working directory path
#' @param models List of model specifications to process
#' @return Numeric value indicating total processing time in seconds
#' @export
process_fba_models <- function(
    working_dir = NULL, 
    models = list(FBA_PARAMETERS$model_tags)) {
  # Initialize working directory
  wd <- working_dir %||% getwd()
  setwd(wd)
  
  # Initialize environment
  dest_dir <- initialize_environment(wd)
  
  # Process all models
  total_time <- 0
  for (model_info in models) {
    elapsed_time <- process_model(model_info, wd, dest_dir)
    total_time <- total_time + elapsed_time
    message(
      sprintf(
        "Time to write %s.txt: %.2f seconds", 
        model_info$name, 
        elapsed_time
      )
    )
  }
  
  message(sprintf("Total time elapsed: %.2f seconds", total_time))
  return(total_time)
}
