
# 0) Setup
wd <- getwd()
source(file.path(wd, "functions_library", "install_and_setup.R"))

# Load core FBA functions
fba_scripts <- c("FBAgreatmodeClass.R", "class_generation.R", "readMat.R")
invisible(lapply(fba_scripts, function(f) {
  source(file.path(wd, "epimod_FBAfunctions", "R", f))
}))

# 1) Define hypernode
hypernode       <- "minimal_doublet.PNPRO"
hypernode_name  <- tools::file_path_sans_ext(hypernode)
hypernode_root  <- file.path(wd, "hypernodes", hypernode_name)

# 2) Create the three canonical subfolders
dirs <- c("config", "src", "output")
for (d in dirs) {
  dir.create(file.path(hypernode_root, d), recursive = TRUE, showWarnings = FALSE)
  cat("Ensured directory:", file.path(hypernode_root, d), "\n")
}

config_dir <- file.path(hypernode_root, "config")
src_dir    <- file.path(hypernode_root, "src")
out_dir    <- file.path(hypernode_root, "output")

# 3) Read hypernode-specific YAML
cfg <- yaml::read_yaml(file.path(
  config_dir,
  paste0("config_", hypernode_name, ".yaml")
))

# 4) Extract parameters
model_names    <- vapply(cfg$organisms, `[[`, character(1), "model_name")
biomass_params <- lapply(cfg$organisms, `[[`, "biomass")
pop_params     <- lapply(cfg$organisms, `[[`, "population")
initial_counts <- as.numeric(vapply(cfg$organisms, `[[`, character(1), "initial_count"))
metabolite_places <- cfg$metabolite_places

# 5) Build organism models and dump parameters
source(file.path(wd, "functions_library", "setup_models.R"))
bacterial_models <- make_bacterial_models(
  model_names, biomass_params, pop_params, initial_counts
)

write_bac_params(
  bacterial_models,
  file.path(config_dir, "organisms_parameters.csv")
)

# 6) Process each model
source(file.path(wd, "functions_library", "process_model.R"))
process_results <- lapply(
  bacterial_models,
  function(m) process_model(m, hypernode_name = hypernode_name)
)

# 7) Project boundary reactions → single output folder
source(file.path(wd, "functions_library", "project_boundary_reactions.R"))
proj_res <- project_boundary_reactions(
  bacterial_models      = bacterial_models,
  metabolite_places     = metabolite_places,
  output_dir_projections = out_dir
)
cat(capture.output(proj_res$bounds), sep = "\n")

# 8) Validate PNPRO → same output folder
source(file.path(wd, "functions_library", "validate_pnpro.R"))
validation <- validate_pnpro(
  pnpro2validate    = file.path(wd, "petri_nets_library", "blank.PNPRO"),
  metadata_path     = hypernode_root,
  bacterial_models  = bacterial_models,
  metabolite_places = metabolite_places,
  validation_dir    = out_dir,
  pnpro_name        = hypernode
)

# 9) Generate the final PNPRO → output folder
source(file.path(wd, "functions_library", "generate_pnpro.R"))
generate_pnpro(
  arc_df   = readr::read_csv(
    file.path(out_dir, paste0(hypernode, "_arc_df_repaired.csv"))
  ),
  pnpro_out = file.path(out_dir, hypernode)
)

# 10) Render layout
system(paste(
  "python3",
  shQuote(file.path(wd, "functions_library", "render_pnpro_layout.py")),
  shQuote(file.path(wd, "petri_nets_library", "blank.PNPRO")),
  shQuote(file.path(out_dir, hypernode))
))

# 1) Read the complete JSON file into a list
cfg <- fromJSON(file.path(config_dir, "boundary_conditions.json"), simplifyVector = TRUE)

# 2) Extract global parameters
background_met <- cfg$background_met
volume                   <- cfg$volume
cell_density             <- cfg$cell_density

C = background_met * volume

# 3) Build a data.frame of exchange bounds
exchange_bounds <- as.data.frame(cfg$exchange_bounds, stringsAsFactors = FALSE)

# 4) (Optional) Convert columns to appropriate types
exchange_bounds$value <- as.numeric(exchange_bounds$value)

##################
## continuing main
##################

# Run extraction with bounds
source(file.path(wd, "epimod_FBAfunctions", "R", "ex_bounds_module.R"))
run_full_ex_bounds(extraction_output = "extracted_ex_reactions.txt", 
                   bacteria_files = bacteria_files,
                   output_dir_projections = "results_ex_reactions", 
                   bacteria_counts = bacteria_counts,
                   fba_upper_bound = rep(1000, length(joint_reactions)), 
                   fba_reactions = joint_reactions,
                   not_shared_base_bound = C, reaction_version = "r")

irr_exchange_bounds <- read.csv(paste0(wd, "/results_ex_reactions/EX_upper_bounds_nonFBA.csv"), head = T)

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
model.generation(net_fname = paste0(wd, "/net/", model_name, ".PNPRO"),
                 transitions_fname = paste0(wd, "/input/General_functions.cpp"),
                 fba_fname = paste0(wd, "/results/", mn_name, ".txt"))
system(paste0("mv ", paste(c(model_name, model_name, model_name, model_name, model_name), 
                           c(".def", ".fbainfo", ".net", ".PlaceTransition", ".solver"), 
                           sep="", collapse=" "), " ", wd, "/net/"))

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
