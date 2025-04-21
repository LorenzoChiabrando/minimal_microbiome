
wd = getwd()

source(paste0(wd, "/install_and_setup.R"))

# Load custom R scripts for Flux Balance Analysis (FBA) functions and utilities
source(paste0(wd, "/epimod_FBAfunctions/R/FBAgreatmodeClass.R"))    # Core FBA class functions
source(paste0(wd, "/epimod_FBAfunctions/R/class_generation.R"))     # Functions for generating FBA models
source(paste0(wd, "/epimod_FBAfunctions/R/readMat.R"))              # MATLAB file reading utilities
source(paste0(wd, "/epimod_FBAfunctions/R/ex_bounds_module.R"))
source(paste0(wd, "/epimod_FBAfunctions/inst/diets/Script4Diets.R"))
source(paste0(wd, "/functions/compiling_mn.R"))
source(paste0(wd, "/functions/process_model.R"))

# Define the output directory for results
dest_dir <- paste0(wd, "/results/")

# Define model configuration
model_name <- "Minimal_EcCb"

# Define bacterial models with consistent parameters
bacterial_models <- list(
  list(
    FBAmodel = "Escherichia_coli_SE11",
    organism = "Escherichia_coli",
    abbreviation = "Ec",
    biomass = list(max = 1.172, mean = 0.489, min = 0.083)
  ),
  list(
    FBAmodel = "Clostridium_butyricum_DSM_10702",
    organism = "Clostridium_butyricum",
    abbreviation = "Cb",
    biomass = list(max = 1.5, mean = 0.6, min = 0.1)
  )
)

# Start timing
start_time <- Sys.time()

# Process each bacterial model sequentially 
process_results <- list()

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

molar = 10 # mmol/mL (1000 mM)
V = 0.001 # mL (1 microL)
C = molar*V # mmol


# Define metabolite places to project across models
metabolite_places <- c("ac_e", "ppa_e", "but_e", "glc__D_e", "lcts_e")

output_dir <- file.path(wd, "results_ex_reactions")

# Initialize output lists
bacteria_files <- character()
bacteria_counts <- numeric()
fba_upper_bounds <- numeric()
shared_reactions <- list()

# Log file for tracking projectable metabolites
log_file <- file.path(output_dir, "projection_log.txt")
cat("Metabolite Projection Analysis\n", file = log_file)
cat("============================\n\n", file = log_file, append = TRUE)

# Process each bacterial model
for ( i in seq_along(bacterial_models) ) {
  model <- bacterial_models[[i]]
  organism <- model$organism
  abbr <- model$abbreviation
  
  cat(sprintf("Analyzing model: %s (%s)\n", organism, abbr), file = log_file, append = TRUE)
  
  # Set file path for compiled model
  model_file <- file.path("compiled_models", paste0(abbr, "_model.txt"))
  
  # Verify model file exists
  if (!file.exists(model_file)) {
    cat(sprintf("Error: Model file not found: %s\n", model_file), file = log_file, append = TRUE)
    next
  }
  
  # Load metabolite metadata
  meta_file <- file.path("input", organism, "metabolites_metadata.csv")
  if (!file.exists(meta_file)) {
    cat(sprintf("Error: Metabolite metadata not found: %s\n", meta_file), file = log_file, append = TRUE)
    next
  }
  
  # Load reaction metadata
  rxn_file <- file.path("input", organism, "reactions_metadata.csv")
  if (!file.exists(rxn_file)) {
    cat(sprintf("Error: Reaction metadata not found: %s\n", rxn_file), file = log_file, append = TRUE)
    next
  }
  
  # Read metadata
  metabolite_data <- read.csv(meta_file)
  reaction_data <- read.csv(rxn_file, stringsAsFactors = FALSE)
  
  # Check which metabolite places can be projected to this model
  cat("Checking metabolite projections:\n", file = log_file, append = TRUE)
  projectable_metabolites <- character()
  
  for (met in metabolite_places) {
    # Check if metabolite exists in this model
    if (any(grepl(met, metabolite_data$id))) {
      projectable_metabolites <- c(projectable_metabolites, met)
      cat(sprintf("  ✓ %s: Found in model\n", met), file = log_file, append = TRUE)
    } else {
      cat(sprintf("  ✗ %s: Not found in model\n", met), file = log_file, append = TRUE)
    }
  }
  
  # Find boundary reactions associated with projectable metabolites
  boundary_reactions <- reaction_data$abbreviation[reaction_data$type == "boundary"]
  
  for (met in projectable_metabolites) {
    
    related_rxns <- character()
    
    # 1. Check reaction equations if available in metadata
    if ("equation" %in% names(reaction_data)) {
      # Get boundary reactions where the metabolite appears in the equation
      eq_matches <- character()
      for (i in which(reaction_data$type == "boundary")) {
        # Check if metabolite is in the equation as a standalone term
        # This uses word boundaries to avoid partial matches
        if (grepl(paste0("\\b", met, "\\b"), reaction_data$equation[i])) {
          eq_matches <- c(eq_matches, reaction_data$abbreviation[i])
        }
      }
      related_rxns <- c(related_rxns, eq_matches)
    }
    
    # Remove duplicates
    related_rxns <- unique(related_rxns)
    
    if (length(related_rxns) > 0) {
      model_shared_reactions <- c(unique(model_shared_reactions), related_rxns)
      cat(sprintf("  Found %d reactions for %s: %s\n", 
                  length(related_rxns), met, paste(related_rxns, collapse=", ")), 
          file = log_file, append = TRUE)
    } else {
      # If no reactions found through analysis, log this information
      cat(sprintf("  No boundary reactions found for %s through direct analysis\n", met), 
          file = log_file, append = TRUE)
      
    }
  }
  
  # Add to output lists if projectable metabolites were found
  if (length(projectable_metabolites) > 0) {
    bacteria_files <- c(bacteria_files, model_file)
    
    # Use default values or those from the model if available
    biomass_value <- ifelse(!is.null(model$biomass$mean), 
                            model$biomass$mean, 
                            default_biomass)
    
    bacteria_count <- ifelse(!is.null(model$initial_count), 
                             model$initial_count, 
                             default_bacteria_count)
    
    bacteria_counts <- c(bacteria_counts, bacteria_count * biomass_value)
    fba_upper_bounds <- c(fba_upper_bounds, default_fba_upper_bound)
    shared_reactions[[abbr]] <- model_shared_reactions
  }
  
  cat("\n", file = log_file, append = TRUE)
}

# Remove duplicate entries for each abbreviation
shared_reactions <- lapply(shared_reactions, unique)
# Compute the joint (union) of all reactions
joint_reactions <- Reduce(union, shared_reactions)

# Only run if we found compatible models
if (length(bacteria_files) > 0) {
  
  # Log summary
  cat("\nSummary:\n", file = log_file, append = TRUE)
  cat(sprintf("Total bacterial models analyzed: %d\n", length(bacterial_models)), file = log_file, append = TRUE)
  cat(sprintf("Compatible models found: %d\n", length(bacteria_files)), file = log_file, append = TRUE)
  cat(sprintf("Total boundary reactions: %d\n", length(joint_reactions)), file = log_file, append = TRUE)
  cat(sprintf("Output files saved to: %s\n", output_dir), file = log_file, append = TRUE)
  
} else {
  cat("No compatible models found for the specified metabolites.\n", file = log_file, append = TRUE)
}

run_full_ex_bounds(
  extraction_output  = "extracted_ex_reactions.txt",
  bacteria_files     = bacteria_files,
  output_dir         = "results_ex_reactions",
  bacteria_counts    = bacteria_counts,
  fba_upper_bound    = rep(1000, length(joint_reactions)),
  fba_reactions      = joint_reactions,
  not_shared_base_bound = C,
  reaction_version   = "r"
)

bounds_file_path = paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")
irr_exchange_bounds = read.csv(bounds_file_path, head = T)

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
