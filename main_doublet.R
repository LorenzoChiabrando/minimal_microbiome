
wd <- getwd()

source(paste0(wd, "/functions/install_and_setup.R"))

# Load necessary scripts
sapply(c("/epimod_FBAfunctions/R/FBAgreatmodeClass.R", 
         "/epimod_FBAfunctions/R/class_generation.R", 
         "/epimod_FBAfunctions/R/readMat.R", 
       "/epimod_FBAfunctions/R/ex_bounds_module.R"),
       function(f) source(paste0(wd, f)))

hypernode = "minimal_doublet.PNPRO"

# Create the corresponding directory inside "hypernodes/"
hypernode_dirname <- tools::file_path_sans_ext(hypernode)
full_dir_path <- file.path(wd, "hypernodes", hypernode_dirname)

if (!dir.exists(full_dir_path)) {
  dir.create(full_dir_path, recursive = TRUE)
  cat(sprintf("Created directory: %s\n", full_dir_path))
} else {
  cat(sprintf("Directory already exists: %s\n", full_dir_path))
}

# List of subdirectory suffixes
suffixes <- c("configs", "projections", "validations", "functions", "metabolic_networks")

# Generate full paths
subdirs <- file.path(full_dir_path, paste0(suffixes, "_", hypernode_dirname))

# Create each subdirectory if not existing
for (subdir in subdirs) {
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", subdir))
  } else {
    cat(sprintf("Directory already exists: %s\n", subdir))
  }
}

# Correct path to the YAML config under new structure
cfg <- yaml::read_yaml(file.path(
  wd, 
  "hypernodes", 
  "minimal_doublet", 
  "configs_minimal_doublet", 
  "config_minimal_doublet.yaml"
))

# 1) Pull out the pieces
model_names    <- vapply(cfg$organisms, `[[`, character(1), "model_name")
biomass_params <- lapply(cfg$organisms, `[[`, "biomass")
pop_params     <- lapply(cfg$organisms, `[[`, "population")
initial_counts <- as.numeric(
  vapply(cfg$organisms, `[[`, character(1), "initial_count")
)

# Load setup_models
source(file.path(wd, "functions_library", "setup_models.R"))

# 2) Build the list (derive_abbrs() runs internally)
bacterial_models <- make_bacterial_models(
  model_names    = model_names,
  biomass_params = biomass_params,
  pop_params     = pop_params,
  initial_counts = initial_counts
)

# 3) Write out organisms_parameters_minimal_doublet.csv
write_bac_params(
  bacterial_models,
  file.path(
    wd,
    "hypernodes",
    "minimal_doublet",
    "configs_minimal_doublet",
    "organisms_parameters_minimal_doublet.csv"
  )
)

# 4) Now you have
pnpro_path       <- paste0(full_dir_path, "/", hypernode)
metabolite_places<- cfg$metabolite_places

source(file.path(wd, "functions_library", "process_model.R"))

process_results <- lapply(bacterial_models, function(model) process_model(model, hypernode_name = hypernode_dirname))

invisible(
  lapply(process_results, function(r) {
    status_icon <- if (r$status == "success") "✓" else "✗"
    msg         <- if (r$status == "success")
      "Successfully processed"
    else
      paste("ERROR -", r$message)
    cat(sprintf("%s  %s (%s): %s\n",
                status_icon,
                r$organism, r$abbr,
                msg))
  })
)

source(paste0(wd, "/functions/project_boundary_reactions.R"))

results_projection <- project_boundary_reactions(
  bacterial_models   = bacterial_models,
  metabolite_places  = metabolite_places,
  output_dir_projections = paste0(wd, "/net/config/projection_info_hypernode_minimal_doublet")
)

cat( capture.output(results_projection$bounds), sep = "\n" )

source(paste0(wd, "/functions/validate_pnpro.R"))

validation = validate_pnpro(pnpro_path,
                            bacterial_models,
                            metabolite_places,
                            validation_dir = paste0(wd, "/net/config/validation_data_hypernode_minimal_doublet"),
                            results_projection)



# build the PNPRO from your repaired arcs
generate_pnpro(arc_df <- readr::read_csv(paste0(wd, "/net/config/validation_data_hypernode_minimal_doublet/hypernode_minimal_doublet_arc_df_repaired.csv")), 
               pnpro_out = file.path(wd, "net/hypernode_minimal_doublet.PNPRO"))

source(paste0(wd, "/functions/layout_generator.R"))

layout_pnpro(input_pnpro = paste0(wd, "/net/hypernode_minimal_doublet.PNPRO"), 
             output_pnpro = paste0(wd, "/net/hypernode_minimal_doublet_layout.PNPRO"))
  
##################
## continuing main
##################

# Set parameters
not_projected_met_molar <- 1000 # mmol/mL 
V <- 0.001       # mL (1 microL)
C <- not_projected_met_molar * V   # mmol
delta <- 1e+10   # density

# Run extraction with bounds
run_full_ex_bounds(extraction_output = "extracted_ex_reactions.txt", bacteria_files = bacteria_files,
                   output_dir_projections = "results_ex_reactions", bacteria_counts = bacteria_counts,
                   fba_upper_bound = rep(1000, length(joint_reactions)), fba_reactions = joint_reactions,
                   not_shared_base_bound = C, reaction_version = "r")

# Process diet medium data
bounds_file_path <- paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv")
irr_exchange_bounds <- read.csv(bounds_file_path, head = T)

diet_medium <- process_medium_data(media_wd = paste0(wd, "/epimod_FBAfunctions/inst/diets/vmh"),
                                   medium_name = "EU_average", bacteria_counts = bacteria_counts,
                                   biomass = biomass)

# Update bounds with diet data
result_df <- irr_exchange_bounds
cleaned_bounds <- sub("_r$", "", irr_exchange_bounds$base_upper_bounds)
for(i in 1:nrow(irr_exchange_bounds)) {
  match_idx <- which(diet_medium$exchange_reaction == cleaned_bounds[i])
  if(length(match_idx) > 0) result_df[[2]][i] <- diet_medium$Flux_in_mmol_human_day[match_idx]
}
write.csv(result_df, bounds_file_path, row.names = F, quote = F)

# Adjust the bounds file
irr_exchange_bounds <- readLines(bounds_file_path)
irr_exchange_bounds[1] <- paste0("base_upper_bounds,", C)
writeLines(irr_exchange_bounds, bounds_file_path)

# Update C++ file with parameters
cpp_file_path <- paste0(wd, "/input/General_functions.cpp")
cpp_content <- readLines(cpp_file_path)
new_vars <- c(paste0("double V = ", V, "; // (mL)"), paste0("long long int delta = ", delta, "; // (cell/mL)"))
cpp_content <- c(cpp_content[1:17], new_vars, cpp_content[20:length(cpp_content)])
writeLines(cpp_content, cpp_file_path)

# Update R functions file
R_funct_file_path <- paste0(wd, "/functions/functions.R")
R_funct_content <- readLines(R_funct_file_path)
R_funct_content[7] <- paste0("  y_ini <- c(", bacteria_counts, ", ", bacteria_counts2, ", ", 
                             10, ", ", 10, ", ", biomass, ", ", biomass, ")")
writeLines(R_funct_content, R_funct_file_path)

# Generate model
start_time <- Sys.time()
model.generation(net_fname = paste0(wd, "/net/", model_name, ".PNPRO"),
                 transitions_fname = paste0(wd, "/input/General_functions.cpp"),
                 fba_fname = paste0(wd, "/results/", mn_name, ".txt"))
system(paste0("mv ", paste(c(model_name, model_name, model_name, model_name, model_name), 
                           c(".def", ".fbainfo", ".net", ".PlaceTransition", ".solver"), 
                           sep="", collapse=" "), " ", wd, "/net/"))
cat("Model generation time:", difftime(Sys.time(), start_time, units="mins"), "minutes\n")

# Analyze model
start_time <- Sys.time()
model.analysis(
  solver_fname = paste0(wd, "/net/", model_name, ".solver"),
  parameters_fname = paste0(wd, "/input/initData.csv"),
  functions_fname = paste0(wd, "/functions/functions.R"),
  debug = TRUE, f_time = 10, s_time = 1, i_time = 0, rchn = 1e-06,
  event_function = NULL, fba_fname = paste0(wd, "/results/", mn_name, ".txt"),
  user_files = c(paste0(wd, "/input/Bacteria_Parameters.csv"), 
                 paste0(wd, "/net/", model_name, ".fbainfo"),
                 paste0(wd, "/results_ex_reactions/EX_upper_bounds_FBA.csv"),
                 paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv"))
)
cat("Model analysis time:", difftime(Sys.time(), start_time, units="mins"), "minutes\n")
