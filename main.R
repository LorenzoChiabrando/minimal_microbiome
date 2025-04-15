
library(devtools)
library(dplyr)
library(R.matlab)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(scales)
library(rlang)
library(parallel)
library(foreach)
library(doParallel)
# remove.packages("epimod")
# install_github("https://github.com/qBioTurin/epimod", ref="epimod_pFBA")
library(epimod)
# downloadContainers()

wd = getwd()

# Load custom R scripts for Flux Balance Analysis (FBA) functions and utilities
source(paste0(wd, "/epimod_FBAfunctions/R/FBAgreatmodeClass.R"))    # Core FBA class functions
source(paste0(wd, "/epimod_FBAfunctions/R/class_generation.R"))     # Functions for generating FBA models
source(paste0(wd, "/epimod_FBAfunctions/R/readMat.R"))              # MATLAB file reading utilities
source(paste0(wd, "/epimod_FBAfunctions/R/ex_bounds_module.R"))
source(paste0(wd, "/epimod_FBAfunctions/inst/diets/Script4Diets.R"))
source(paste0(wd, "/functions/compiling_mn.R"))

# Define the output directory for results
dest_dir <- paste0(wd, "/results/")

# Define model configuration
model_name <- "Minimal_EcCb"

# Define bacterial models with consistent parameters
bacterial_models <- list(
  list(
    FBAmodel = "Escherichia_coli_str_K_12_substr_MG1655",
    organism = "Escherichia_coli",
    abbreviation = "Ec",
    biomass = list(max = 1.172, mean = 0.489, min = 0.083)
  ),
  list(
    FBAmodel = "Clostridium_butyricum_DSM_10702",
    organism = "Clostridium_butyricum",
    abbreviation = "Cb",
    biomass = list(max = 1.5, mean = 0.6, min = 0.1)  # Customize these values
  )
)

# Start timing
start_time <- Sys.time()

# Process each bacterial model sequentially 
process_results <- list()

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
  
  # Get model components for classification
  S <- model_obj@S
  react_id <- model_obj@react_id
  met_id <- model_obj@met_id
  lb <- model_obj@lowbnd
  ub <- model_obj@uppbnd
  
  # Generate reaction classification
  cat("Classifying reactions...\n")
  n_rxns <- length(react_id)
  reaction_types <- character(n_rxns)
  reaction_subtypes <- character(n_rxns)
  
  # Your reaction classification code here
  # This is a placeholder - add your actual classification logic
  for (i in seq_len(n_rxns)) {
    if (grepl("^EX_", react_id[i])) {
      reaction_types[i] <- "boundary"
      reaction_subtypes[i] <- "exchange"
    } else if (grepl("^DM_", react_id[i])) {
      reaction_types[i] <- "boundary"
      reaction_subtypes[i] <- "demand"
    } else if (grepl("^SINK_", react_id[i])) {
      reaction_types[i] <- "boundary"
      reaction_subtypes[i] <- "sink"
    } else {
      reaction_types[i] <- "core"
      reaction_subtypes[i] <- "internal"
    }
  }
  
  # Create reactions dataframe
  cat("Creating reaction metadata...\n")
  reactions_df <- data.frame(
    abbreviation = react_id,
    type = reaction_types,
    subtype = reaction_subtypes,
    lowbnd = lb,
    uppbnd = ub,
    stringsAsFactors = FALSE
  )
  
  # Generate metabolite classification
  cat("Classifying metabolites...\n")
  n_mets <- length(met_id)
  is_core <- logical(n_mets)
  is_boundary <- logical(n_mets)
  
  # Your metabolite classification code here
  # This is a placeholder - add your actual classification logic
  for (i in seq_len(n_mets)) {
    # Simple placeholder logic
    is_core[i] <- TRUE
    is_boundary[i] <- FALSE
  }
  
  # Create metabolites dataframe
  cat("Creating metabolite metadata...\n")
  met_metadata <- data.frame(
    id = met_id,
    is_core = is_core,
    is_boundary = is_boundary,
    stringsAsFactors = FALSE
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

# Process each model one at a time
for (i in seq_along(bacterial_models)) {
  process_results[[i]] <- process_model(bacterial_models[[i]])
}

# Display summary
cat("\nProcessing Summary:\n")
cat("-----------------\n")

for (result in process_results) {
  if (result$status == "success") {
    cat(sprintf("✓ %s (%s): Successfully processed\n", result$organism, result$abbr))
  } else {
    cat(sprintf("✗ %s (%s): ERROR - %s\n", result$organism, result$abbr, result$message))
  }
}

# Calculate and display elapsed time
end_time <- Sys.time()
elapsed <- end_time - start_time
cat(sprintf("\nAll models processed in %.2f seconds\n", as.numeric(elapsed)))

met_places = c("ac_e", "ppa_e", "but_e", "glc__D_e", "ltcs_e")
FBA_possible_reactions  = ...

# Define the FBA models and associated counts
bacteria_txt <- c(paste0(wd, "/compiled_mdoels/", bacterial_models[[1]][["abbreviation"]], ".txt"), 
                     paste0(wd, "/compiled_mdoels/", bacterial_models[[1]][["abbreviation"]], ".txt"))

biomass_0 = 1 # pgWD
biomass_0 = 1 # pgWD

bacteria_counts <- 1000000 # (cell)
bacteria_counts <- 1000000 # (cell)

run_full_ex_bounds(
  extraction_output  = "extracted_ex_reactions.txt",
  bacteria_files     = bacteria_txt,
  output_dir         = "results_ex_reactions",
  bacteria_counts    = c(bacteria_counts*biomass, bacteria_counts2*biomass),
  fba_upper_bound    = c(1000, 1000),
  fba_reactions      = c("EX_glc_D_e", "EX_lcts_e", "EX_ppa_e", "EX_ac_e"),
  not_shared_base_bound = C,
  reaction_version   = "r"
)

bounds_file_path = paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")
irr_exchange_bounds = read.csv(bounds_file_path, head = T)

molar = 10 # mmol/mL (1000 mM)
V = 0.001 # mL (1 microL)
C = molar*V # mmol

diet_medium = process_medium_data(media_wd = paste0(wd, "/epimod_FBAfunctions/inst/diets/vmh"),
                                  medium_name = "EU_average",
                                  bacteria_counts = bacteria_counts,
                                  biomass = biomass)

# Create a copy of the original dataframe
result_df <- irr_exchange_bounds

# Remove the "_r" suffix from base_upper_bounds for matching
cleaned_bounds <- sub("_r$", "", irr_exchange_bounds$base_upper_bounds)

# For each row in irr_exchange_bounds
for( i in 1:nrow(irr_exchange_bounds) ) {
  # Find matching exchange reaction in diet_medium
  match_idx <- which(diet_medium$exchange_reaction == cleaned_bounds[i])
  
  # If there's a match, substitute the value
  if(length(match_idx) > 0) {
    result_df[[2]][i] <- diet_medium$Flux_in_mmol_human_day[match_idx]
  }
}

write.csv(result_df, bounds_file_path, row.names = F, quote = F)

irr_exchange_bounds = readLines(bounds_file_path)
irr_exchange_bounds[1] = paste0("base_upper_bounds,", C)
writeLines(irr_exchange_bounds, bounds_file_path)

Bacteria_Parameters <- read.csv(paste0(wd, "/input/Bacteria_Parameters.csv"), head = F)
Bacteria_Parameters[1] <- 0.21 # starv
Bacteria_Parameters[2] <- 1 # dup
Bacteria_Parameters[3] <- 0.018 # death

write.table(Bacteria_Parameters,
            paste0(wd, "/input/Bacteria_Parameters.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = F)

# Adjust General_functions.cpp
delta = 1e+10 # density

# Read the C++ file
cpp_file_path <- paste0(wd, "/input/General_functions.cpp")
cpp_content <- readLines(cpp_file_path)

# Define the new variables to insert
new_variables <- c(
  paste0("double V = ", V, "; // (mL) (1 microL)"),
  paste0("long long int delta = ", delta, "; // (cell/mL)")
)

# Find position after includes and namespace (looking for first empty line)
insert_position <- which(cpp_content == "")[1]

# Insert the new variables after the includes
cpp_content <- c(
  cpp_content[1:17],
  new_variables,
  cpp_content[20:length(cpp_content)]
)

# Write the modified content back to the file
writeLines(cpp_content, cpp_file_path)

R_funct_file_path <- paste0(wd, "/functions/functions.R")
R_funct_content <- readLines(R_funct_file_path)

# yini.names <- c("E_coli", "E_coli_biomass_e", "lcts_e", "glc_D_e", "Clost_buty", "Clost_buty_biomass_e")
# Insert the new variables after the includes
R_funct_content <- c(
  R_funct_content[1:6],
  paste0("  y_ini <- c(", bacteria_counts, ", ", 
         bacteria_counts2, ", ", 10, ", " , 10, ", ", biomass, ", ", biomass, ")"),
  R_funct_content[8:length(R_funct_content)]
)

# Write the modified content back to the file
writeLines(R_funct_content, R_funct_file_path)

start_time <- Sys.time()

model.generation(net_fname = paste0(wd, "/net/", model_name, ".PNPRO"),
                 transitions_fname = paste0(wd, "/input/General_functions.cpp"),
                 fba_fname = paste0(wd, "/results/", mn_name, ".txt"))

end_time <- Sys.time()

elapsed_time <- end_time - start_time
units(elapsed_time) <- "mins"
cat("Time elapsed:", elapsed_time, "minutes\n")

system(paste0("mv ", 
              model_name, ".def ", 
              model_name, ".fbainfo ", 
              model_name, ".net ", 
              model_name, ".PlaceTransition ", 
              model_name, ".solver ", wd, "/net/"))

start_time <- Sys.time()

model.analysis(
  solver_fname    = paste0(wd, "/net/", model_name, ".solver"),
  parameters_fname= paste0(wd, "/input/initData.csv"),
  functions_fname = paste0(wd, "/functions/functions.R"),
  debug           = TRUE,
  f_time          = 10,
  s_time          = 1,
  i_time          = 0,
  rchn            = 1e-06,
  event_function  = NULL,
  fba_fname       = paste0(wd, "/results/", mn_name, ".txt"),
  user_files      = c(
    paste0(wd, "/input/Bacteria_Parameters.csv"),
    paste0(wd, "/net/", model_name, ".fbainfo"),
    paste0(wd, "/results_ex_reactions/EX_upper_bounds_FBA.csv"),
    paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")
  )
)


end_time <- Sys.time()

elapsed_time <- end_time - start_time
units(elapsed_time) <- "mins"
cat("Time elapsed:", elapsed_time, "minutes\n")

p = plot_analysis(reactions_of_interest = c("EX_biomass_e_f", "EX_biomass_e_r",
                                            "EX_glc_D_e_r", "EX_glc_D_e_f",
                                            "EX_lcts_e_r", "EX_lcts_e_f",
                                            "EX_but_e_r", "EX_but_e_f",
                                            "EX_ppa_e_r", "EX_ppa_e_f",
                                            "EX_ac_e_r", "EX_ac_e_f"), ncol_reactions = 3,
                  place2plot = c(abbr, paste0(abbr, "_biomass_e")), ncol_places = 1)

#Error in plot_analysis(reactions_of_interest = c("EX_biomass_e_f", "EX_biomass_e_r",  : 
#could not find function "plot_analysis"

ggsave(paste0(model_name, "_",  "Analysis_results.pdf"), 
       p[[1]] + p[[2]] + (p[[3]] / p[[4]] / p[[5]]), 
       width = 22, height = 9)

file.remove(list.files(path = wd, pattern = "\\.log$", full.names = TRUE))
file.remove(list.files(path = wd, pattern = "\\ID$", full.names = TRUE))
file.remove(list.files(path = wd, pattern = "\\StatusFile$", full.names = TRUE))
gc()

Bacteria_Parameters <- read.csv(paste0(wd, "/input/Bacteria_Parameters.csv"), head = F)

starv_int = c(0.15, 0.5)
dup = Bacteria_Parameters[2]
death = Bacteria_Parameters[3]

init_data <- read.csv(paste0(wd, "/input/initData.csv"), head = F)
init_data[2, ] <- paste0("g; Bacteria_Parameters.csv; psensitivty; n = 1; min = ", 
                         starv_int[1], "; max = ", starv_int[2], 
                         "; dup = ", dup, "; death = ", death)
write.table(init_data,
            paste0(wd, "/input/initDataSensitivity.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = F)

start_time <- Sys.time()

model.sensitivity(solver_fname = paste0(wd, "/net/", model_name, ".solver"),
                  fba_fname = paste0(wd, "/results/", mn_name, ".txt"),
                  i_time = 0, 
                  f_time = 72, 
                  s_time = 0.25,
                  n_config = parallel::detectCores()*2,
                  debug = T,
                  rchn = 1e-06,
                  event_function = NULL,
                  target_value = c(abbr, paste0(abbr, "_biomass_e")),
                  parameters_fname = paste0(wd, "/input/initDataSensitivity.csv"),
                  parallel_processors = parallel::detectCores(),
                  functions_fname = paste0(wd, "/functions/functions.R"),
                  user_files = c(paste0(wd, "/net/", model_name, ".fbainfo"),
                                 paste0(wd, "/results_ex_reactions/EX_upper_bounds_FBA.csv"),
                                 paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")))

end_time <- Sys.time()

elapsed_time <- end_time - start_time
units(elapsed_time) <- "mins"
cat("Time elapsed:", elapsed_time, "minutes\n")

# results <- process_sensitivity_analysis(wd = wd, save_type = "both")
# Save only trace files
trace_results <- process_sensitivity_analysis(wd = wd, save_type = "trace")

# Get the plots
plots_places <- analyze_trace_sensitivity(wd = wd, model_name = model_name)

ggsave(paste0(model_name, "_",  "Sensitivity_places.pdf"), 
       plots_places[[1]] + plots_places[[2]], 
       width = 8, height = 3)

# Run the analysis
results <- analyze_flux_sensitivity(
  wd = wd,
  model_name = model_name,
  mn_name = mn_name,
  reactions_of_interest = c("EX_biomass_e_f", "EX_biomass_e_r",
                            "EX_glc_D_e_r", "EX_glc_D_e_f",
                            "EX_lcts_e_r", "EX_lcts_e_f",
                            "EX_but_e_r", "EX_but_e_f",
                            "EX_ppa_e_r", "EX_ppa_e_f",
                            "EX_ac_e_r", "EX_ac_e_f"))

# Access results
# sensitivity_data <- results$sensitivity_data
# combined_data <- results$combined_data
plots <- results$plots

ggsave(paste0(model_name, "_",  "Sensitivity_fluxes.pdf"), 
       ((plots[[1]] | plots[[2]]) /
          (plots[[3]] | plots[[4]]) /
          (plots[[5]] | plots[[6]])) |
         ((plots[[7]] | plots[[8]]) /
            (plots[[9]] | plots[[10]]) /
            (plots[[11]] | plots[[12]])),
       width = 12, height = 7)

file.remove(list.files(path = wd, pattern = "\\.log$", full.names = TRUE))
file.remove(list.files(path = wd, pattern = "\\ID$", full.names = TRUE))
file.remove(list.files(path = wd, pattern = "\\StatusFile$", full.names = TRUE))
gc()
