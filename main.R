
library(devtools)
library(dplyr)
library(R.matlab)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(scales)
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

model_name = "BacPN_FBAdouble"
mn_name = "Escherichia_coli_str_K_12_substr_MG1655"
org_name = "Escherichia_coli"
abbr = "E_coli"

mn_name2 = "Clostridium_butyricum_DSM_10702"
org_name2 = "Clostridium_butyricum"
abbr2 = "Clost_buty"

# Define the models to be processed
models <- list(list(name = mn_name, type = org_name), list(name = mn_name2, type = org_name2))

# Define the output directory for results
dest_dir <- paste0(wd, "/results/")

# Start timing
start_time <- Sys.time()

# Loop through each model for FBA processing
for (model_info in models) {
  fba_name <- model_info$name    # Extract the model name
  model_type <- model_info$type    # Extract the model type
  
  # Construct file paths for FBA input
  fba_fname <- paste0(fba_name, ".txt")
  matfile <- paste0(wd, '/input/', model_type, "/", fba_name, ".mat")
  cat("METABOLIC NETWORK: ", fba_name, " ( ", model_type, " )\n")
  
  # Check if the .mat file exists
  if (!file.exists(matfile)) {
    stop(paste("File not found:", matfile))  # Stop execution if file is missing
  }
  
  # Generate the model for the specified .mat file
  cat("Generating model for:", fba_name, "\n")
  model <- FBA4Greatmod.generation(fba_mat = matfile)
  
  # Set biomass parameters for the FBA model
  model <- setBiomassParameters(model, bioMax = 1.172, bioMean = 0.489, bioMin = 0.083)
  
  # Prepare the output path for .rds files
  output_rds_path <- paste0(wd, "/input/", model_type)
  files_to_move <- list.files(pattern = "\\.rds$", full.names = TRUE)
  
  # Move .rds files to the appropriate directory
  if (length(files_to_move) > 0) {
    sapply(files_to_move, function(f) {
      file.rename(from = f, to = file.path(output_rds_path, basename(f)))
    })
  } else {
    cat("No .rds files found for", fba_name, "\n")
  }
}

# End timing
end_time <- Sys.time()

# Force to minutes
elapsed_time <- end_time - start_time
units(elapsed_time) <- "mins"
cat("Time elapsed:", elapsed_time, "minutes\n")

FBA_PARAMETERS <- list(biomass = list(max = 1.172, mean = 0.489, min = 0.083),
                       model_tags = list(name = mn_name, type = org_name))

start_time <- Sys.time()

process_fba_models(working_dir = NULL, list(FBA_PARAMETERS$model_tags))

end_time <- Sys.time()

elapsed_time <- end_time - start_time
units(elapsed_time) <- "mins"
cat("Time elapsed:", elapsed_time, "minutes\n")

# This function extracts all the specified EX_ reactions from the provided models (bacteria_files)
# and calculates their lower bound as (bacteria_counts*biomass)/MW.

# NotFBA exchange reactions LOWER BOUNDS are denoting metabolite concentration units. 
# Essential unit conversions include:
#
# - Biomass in grams of dry weight (gDW)
# - Environmental volume in milliliters (mL)
# - Metabolite concentration in millimolar (mM = mmol/L = 0.001 mmol/mL)
# 
# The environment is defined by volume (V) bacterial cell counts, and medium composition. 
# Conversion between molar concentrations and absolute quantities follows:

# Define the FBA models and associated counts
bacteria_models <- c(paste0(wd, "/results/", mn_name, ".txt"), paste0(wd, "/results/", mn_name2, ".txt"))
bacteria_counts <- 1000 # (cells)
bacteria_counts2 <- 500 # (cells)

molar = 10 # mmol/mL (1000 mM)
V = 0.001 # mL (1 microL)
C = molar*V # mmol
biomass = 1 # pgWD

run_full_ex_bounds(
  extraction_output  = "extracted_ex_reactions.txt",
  bacteria_files     = bacteria_models,
  output_dir         = "results_ex_reactions",
  bacteria_counts    = c(bacteria_counts*biomass, bacteria_counts2*biomass),
  fba_upper_bound    = c(100,50),
  fba_reactions      = c("EX_glc_D_e", "EX_lcts_e"),
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

# Insert the new variables after the includes
R_funct_content <- c(
  R_funct_content[1:6],
  paste0("  y_ini <- c(", bacteria_counts, ", ", biomass, ")"),
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
  f_time          = 72,
  s_time          = 0.25,
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
