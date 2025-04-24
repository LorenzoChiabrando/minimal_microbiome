
wd <- getwd()
dest_dir <- paste0(wd, "/results/")
output_dir_projections <- file.path(wd, "results_ex_reactions")

# Load necessary scripts
sapply(c("/functions/install_and_setup.R", "/functions/process_model.R", 
         "/epimod_FBAfunctions/R/FBAgreatmodeClass.R", "/epimod_FBAfunctions/R/class_generation.R", 
         "/epimod_FBAfunctions/R/readMat.R"), function(f) source(paste0(wd, f)))

# Model configuration
model_name <- "Minimal_EcCb"
metabolite_places <- c("glc__D_e", "lcts_e")
model_names <- c("Escherichia_coli_SE11", "Clostridium_butyricum_DSM_10702")
model_abbr_candidates <- list()

for (file_name in model_names) {
  abbr1 <- str_match(file_name, "^([^_]+)_model$")
  model_abbr_candidates[[file_name]] <- if(!is.na(abbr1[1,1])) abbr1[1,2] else character(0)
  
  name_parts <- strsplit(tools::file_path_sans_ext(file_name), "_")[[1]]
  if (length(name_parts) > 1) {
    abbr2 <- tolower(paste0(substr(name_parts, 1, 1), collapse=""))
    model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], abbr2)
    
    if (length(name_parts) == 2) {
      model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], 
                                              tolower(paste0(substr(name_parts[1], 1, 1), substr(name_parts[2], 1, 1))))
    }
  }
}

# Define bacterial models
bacterial_models <- list(
  list(FBAmodel = "Escherichia_coli_SE11", organism = "Escherichia_coli",
       abbreviation = model_abbr_candidates$Escherichia_coli_SE11, txt_file = "Escherichia_coli_SE11",
       biomass = list(max = 1.172, mean = 0.489, min = 0.083),
       bac_pop_p = list(starv = 0.21, dup = 1, death = 0.018), initial_count = 1e+06),
  list(FBAmodel = "Clostridium_butyricum_DSM_10702", organism = "Clostridium_butyricum",
       abbreviation = model_abbr_candidates$Clostridium_butyricum_DSM_10702, txt_file = "Clostridium_butyricum_DSM_10702",
       biomass = list(max = 1.5, mean = 0.6, min = 0.1),
       bac_pop_p = list(starv = 0.21, dup = 1, death = 0.018), initial_count = 1e+06)
)

# Update parameters file
Bacteria_Parameters <- data.frame(do.call(rbind, lapply(bacterial_models, function(m) as.vector(unlist(m$bac_pop_p)))))
write.table(Bacteria_Parameters, paste0(wd, "/input/Bacteria_Parameters.csv"), 
            sep = ",", row.names = FALSE, col.names = FALSE, quote = F)

# Process models and display summary
process_results <- lapply(bacterial_models, process_model)
cat("\nProcessing Summary:\n-----------------\n")
lapply(process_results, function(r) {
  cat(sprintf("%s %s (%s): %s\n", 
              ifelse(r$status == "success", "✓", "✗"), 
              r$organism, r$abbr, 
              ifelse(r$status == "success", "Successfully processed", paste("ERROR -", r$message))))
})

# Load additional scripts and configure environment
source(paste0(wd, "/functions/validating_projection.R"))
source(paste0(wd, "/functions/pnpro_validation.R"))
source(paste0(wd, "/functions/repairing_pnpro.R"))

# Set parameters
molar <- 10      # mmol/mL 
V <- 0.001       # mL (1 microL)
C <- molar * V   # mmol
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