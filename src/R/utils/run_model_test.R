# UGM_Microbiota/src/R/utils/run_model_test.R

# --- 0. Setup Environment and Load Libraries ---
# Ensure you are in the UGM_Microbiota/ root directory when running this script
# setwd("path/to/your/UnifiedGreatMod_Microbiota_Community_Models/")

# Install/load epimod_FBAfunctions_enhanced package
if (!requireNamespace("epimod_FBAfunctions_enhanced", quietly = TRUE)) {
  message("Installing epimod_FBAfunctions_enhanced from local source...")
  install.packages("src/R/epimod_FBAfunctions_enhanced", repos = NULL, type = "source")
  library(epimod_FBAfunctions_enhanced)
} else {
  library(epimod_FBAfunctions_enhanced)
}

if (!requireNamespace("epimod", quietly = TRUE)) {
  message("epimod package not found. Please install it (e.g., install.packages('epimod')).")
  stop("epimod package is required for this test.")
} else {
  library(epimod)
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml")
  library(yaml)
} else {
  library(yaml)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
  library(jsonlite)
} else {
  library(jsonlite)
}

# --- Function to run a generic model test ---
run_model_test <- function(model_label,
                           base_dir = "models",
                           pnpro_sub_dir = "pnpro",
                           data_sub_dir = "data",
                           config_sub_dir = "config",
                           output_sub_dir = "output",
                           solver_type = "LSODA",
                           sim_f_time = 24,
                           sim_s_time = 1,
                           sim_timeout = "5m",
                           epsilon_fba = 1e-6) { # Add epsilon for FBA automation
  
  message(paste0("\n--- Starting test for model: ", model_label, " ---"))
  
  # --- 1. Define Paths Dynamically ---
  model_path <- file.path(base_dir, model_label)
  
  # Check if model directory exists
  if (!dir.exists(model_path)) {
    stop(paste("Model directory not found:", model_path))
  }
  
  pnpro_file <- file.path(model_path, pnpro_sub_dir, paste0(model_label, ".PNPRO"))
  config_yaml_file <- file.path(model_path, config_sub_dir, paste0("config_", model_label, ".yaml"))
  boundary_json_file <- file.path(model_path, config_sub_dir, "boundary_conditions.json")
  
  raw_mats_dir <- file.path(model_path, data_sub_dir, "raw_mats")
  biounits_dir <- file.path(model_path, data_sub_dir, "biounits")
  
  solver_output_dir <- file.path(model_path, output_sub_dir, "solver")
  simulation_results_dir <- file.path(model_path, output_sub_dir, "simulation_results")
  
  # Create output directories if they don't exist
  if (!dir.exists(biounits_dir)) dir.create(biounits_dir, recursive = TRUE)
  if (!dir.exists(solver_output_dir)) dir.create(solver_output_dir, recursive = TRUE)
  if (!dir.exists(simulation_results_dir)) dir.create(simulation_results_dir, recursive = TRUE)
  
  # --- 2. Read Configuration Files ---
  message("\n--- Reading configuration files ---")
  
  config_yaml <- tryCatch(yaml::read_yaml(config_yaml_file), 
                          error = function(e) stop(paste("Error reading YAML config:", e$message)))
  boundary_json <- tryCatch(jsonlite::read_json(boundary_json_file), 
                            error = function(e) stop(paste("Error reading JSON boundary conditions:", e$message)))
  
  # --- 3. Test .mat to .txt Conversion (for all cellular units in config_yaml) ---
  message("\n--- Step 3: Testing .mat to .txt conversion for all cellular units ---")
  
  processed_txt_files <- list()
  initial_counts_list <- list()
  
  for (unit in config_yaml$cellular_units) {
    model_name <- unit$model_name
    mat_file_path <- file.path(raw_mats_dir, paste0(model_name, ".mat"))
    txt_file_path <- file.path(biounits_dir, paste0(gsub(" ", "_", tolower(model_name)), "_model.txt"))
    
    if (!file.exists(mat_file_path)) {
      stop(paste("MAT file not found for", model_name, ":", mat_file_path))
    }
    
    message(paste("Converting", mat_file_path, "to", txt_file_path))
    fba_model <- FBA4Greatmod.generation(mat_file_path)
    
    # Apply biomass parameters from config_yaml
    if (!is.null(unit$biomass)) {
      fba_model <- setBiomassParameters(fba_model, 
                                        bioMax = unit$biomass$max, 
                                        bioMean = unit$biomass$mean, 
                                        bioMin = unit$biomass$min)
    }
    writeFBAfile(fba_model, file = txt_file_path)
    
    processed_txt_files[[model_name]] <- txt_file_path
    initial_counts_list[[unit$label]] <- unit$initial_count # Store initial population count
    initial_counts_list[[paste0(unit$label, "Biom")]] <- unit$biomass$mean # Store initial biomass mean
  }
  message("--- .mat to .txt conversion complete. Inspect files in data/biounits/ ---")
  
  # --- 4. Prepare initial conditions for simulation ---
  # This part needs to be customized based on how your PNPRO places are named
  # and what initial values they expect beyond population counts and biomass.
  # For the minimal doublet, it seems to expect 'lac_L_e', 'ac_e', 'but_e', 'ppa_e' from boundary_conditions.json
  
  initial_conditions_sim <- c()
  
  # Add initial counts for cellular units
  for (label in names(initial_counts_list)) {
    initial_conditions_sim[label] <- initial_counts_list[[label]]
  }
  
  # Add initial concentrations for boundary metabolites, assuming they map directly to PN places
  # You need to manually specify what metabolites from boundary_conditions.json
  # map to places in your PNPRO (e.g., 'lac_L_e', 'ac_e', etc.) and their initial concentrations.
  # For a simple test, we can hardcode for the known metabolites in minimal_doublet.PNPRO
  # In a real scenario, these would come from 'initial_data.csv' or be derived.
  
  # Example for minimal_doublet:
  # Initial concentrations from config_minimal_doublet.yaml "boundary_metabolites"
  # These are conceptual in the yaml, but their initial values might come from initial_data.csv or be fixed.
  # For this test, let's assume some non-zero starting values for simple testing.
  initial_conditions_sim["lac_L_e"] <- 10 # Example starting concentration
  initial_conditions_sim["ac_e"] <- 10    # Example starting concentration
  initial_conditions_sim["but_e"] <- 10   # Example starting concentration
  initial_conditions_sim["ppa_e"] <- 10   # Example starting concentration
  
  # Sort initial_conditions_sim by name to ensure consistent ordering if epimod expects it (good practice)
  initial_conditions_sim <- initial_conditions_sim[order(names(initial_conditions_sim))]
  
  message("\n--- Prepared initial conditions for simulation: ---")
  print(initial_conditions_sim)
  
  # --- 5. Generate Solver from PNPRO (triggers C++ FBA automation) ---
  message("\n--- Step 5: Generating solver from PNPRO (including FBA automation) ---")
  
  solver_name_path <- file.path(solver_output_dir, model_label)
  
  # The FBA[] command in the PNPRO will use the .txt files generated in Step 3.
  # The epsilon_fba parameter can be passed to the model_generation, which
  # will then be used by the C++ FBA automation logic.
  model_generation(
    net_fname = pnpro_file, 
    out_fname = solver_name_path,
    epsilon = epsilon_fba # Passing epsilon to the model_generation
  )
  
  message(paste0("--- Solver generation complete. Check for '", model_label, ".solver' in ", solver_output_dir, " ---"))
  if (file.exists(paste0(solver_name_path, ".solver"))) {
    message(paste0("Solver file '", model_label, ".solver' successfully created."))
  } else {
    stop("ERROR: Solver file was not created. Check model_generation logs for errors.")
  }
  
  # --- 6. Run Basic Simulation to Test FBA Logic ---
  message("\n--- Step 6: Running a basic simulation to test FBA logic ---")
  
  tryCatch({
    sim_results <- model_analysis(
      solver_fname = paste0(solver_name_path, ".solver"),
      f_time = sim_f_time, 
      s_time = sim_s_time,
      ini_v = initial_conditions_sim,
      timeout = sim_timeout, 
      out_fname = file.path(simulation_results_dir, paste0(model_label, "_test_simulation_output")),
      solver_type = solver_type
    )
    message(paste0("--- Basic simulation complete for ", model_label, ". Check output in ", simulation_results_dir, " ---"))
    message("First few rows of simulation results:")
    print(head(sim_results))
    
  }, error = function(e) {
    message(paste0("ERROR during simulation for ", model_label, ":"))
    message(e$message)
    stop("Simulation failed. Check the solver compilation and FBA logic.")
  })
  
  message(paste0("\n--- Test for model: ", model_label, " Concluded ---"))
  message("Please manually verify:")
  message(paste0("1. The sizes of .txt files in `", biounits_dir, "` (should be small, indicating sparse format)."))
  message("2. The presence of biomass parameters (BioMax, BioMean, BioMin) in the generated .txt files.")
  message(paste0("3. The successful creation of the .solver file in `", solver_output_dir, "`."))
  message(paste0("4. The simulation output in `", simulation_results_dir, "` for reasonable initial trends."))
}

# --- Example Usage (call this from your main R script or directly) ---
# To run the EcCb test:
# source("src/R/utils/run_model_test.R") # Load the function first
# run_model_test(model_label = "minimal_doublet")
#
# You would ensure your directory structure is like:
# models/minimal_doublet/
# ├── pnpro/minimal_doublet.PNPRO
# ├── data/raw_mats/Escherichia_coli_SE11.mat
# ├── data/raw_mats/Clostridium_butyricum_DSM_10702.mat
# ├── config/config_minimal_doublet.yaml
# └── config/boundary_conditions.json

# # From your UGM_Microbiota/ root directory:
# source("src/R/utils/run_model_test.R") 
# 
# # To run the test for your minimal_doublet model:
# run_model_test(model_label = "minimal_doublet", 
#                sim_f_time = 24, # Example: 24 hours simulation
#                sim_s_time = 0.1, # Example: record every 0.1 hours
#                epsilon_fba = 1e-6) # Example: strict FBA re-calculation
