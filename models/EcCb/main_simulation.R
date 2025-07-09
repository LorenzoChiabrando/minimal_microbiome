
# setwd("~/UGM_Microbiota")

wd = getwd()

source(file.path("src/R/utils/install_and_setup.R"))

install.packages("src/R/epimod_FBAfunctions", repos = NULL, type = "source")

# Load core FBA functions
fba_scripts <- c("FBAgreatmodeClass.R", "class_generation.R", "readMat.R")
invisible(lapply(fba_scripts, function(f) {
  source(file.path("src/R/epimod_FBAfunctions/R/", f))
}))

hypernode       <- "minimal_doublet.PNPRO"
hypernode_name  <- tools::file_path_sans_ext(hypernode)
hypernode_root  <- file.path(wd, "hypernodes", hypernode_name)

dirs <- c("config", "src", "output")
for (d in dirs) {
  dir.create(file.path(hypernode_root, d), recursive = TRUE, showWarnings = FALSE)
  cat("Ensured directory:", file.path(hypernode_root, d), "\n")
}

config_dir <- file.path(hypernode_root, "config")
src_dir    <- file.path(hypernode_root, "src")
out_dir    <- file.path(hypernode_root, "output")

cfg <- yaml::read_yaml(file.path(config_dir, paste0("config_", hypernode_name, ".yaml")))

model_names    <- vapply(cfg$cellular_units, `[[`, character(1), "model_name")
biomass_params <- lapply(cfg$cellular_units, `[[`, "biomass")
pop_params     <- lapply(cfg$cellular_units, `[[`, "population")
initial_counts <- as.numeric(vapply(cfg$cellular_units, `[[`, character(1), "initial_count"))
boundary_metabolites <- cfg$boundary_metabolites

source(file.path(wd, "functions_library", "setup_models.R"))
biounit_models <- make_biounit_models(model_names, biomass_params, pop_params, initial_counts)
write_population_params(biounit_models, file.path(config_dir, "population_parameters.csv"))

source(file.path(wd, "functions_library", "process_model.R"))
process_results <- lapply(
  biounit_models,
  function(m) process_model(m, hypernode_name = hypernode_name)
)

source(file.path(wd, "functions_library", "project_boundary_reactions.R"))
proj_res <- project_boundary_reactions(
  biounit_models = biounit_models,
  boundary_metabolites = boundary_metabolites,
  out_dir = out_dir,
  hypernode_name = hypernode_name
)
cat(capture.output(proj_res$bounds), sep = "\n")

source(file.path(wd, "functions_library", "validate_pnpro.R"))
validation <- validate_pnpro(
  pnpro2validate    = file.path(wd, "petri_nets_library", "blank.PNPRO"),
  hypernode_root    = hypernode_root,
  biounit_models    = biounit_models,
  boundary_metabolites = boundary_metabolites,
  out_dir           = out_dir,
  hypernode_name    = hypernode_name
)

source(file.path(wd, "functions_library", "generate_pnpro.R"))
generate_pnpro(
  arc_df   = readr::read_csv(file.path(out_dir, "repaired_arcs.csv")),
  pnpro_out = file.path(wd, "petri_nets_library", hypernode)
)

cfg <- fromJSON(file.path(config_dir, "boundary_conditions.json"), simplifyVector = TRUE)

background_met <- cfg$background_met
volume <- cfg$volume
cell_density <- cfg$cell_density
projected_base_ub <- cfg$fba_upper_bound
projected_base_lb <- cfg$fba_lower_bound
not_projected_base_lb = background_met * volume
not_projected_base_ub = 1000

source(file.path(wd, "epimod_FBAfunctions", "R", "ex_bounds_module.R"))
run_full_ex_bounds(hypernode_name = hypernode_name,
                   biounit_models = biounit_models,
                   projected_base_ub = projected_base_ub,
                   projected_base_lb = projected_base_lb,
                   not_projected_base_lb = not_projected_base_lb,
                   not_projected_base_ub = not_projected_base_ub)

# Load not-projected bounds CSV
not_proj_df <- read.csv(file.path(wd, "hypernodes", hypernode_name, "output", "ub_bounds_not_projected.csv"),
                        header = FALSE, col.names = c("reaction", "FBAmodel", "upper_bound"),
                        stringsAsFactors = FALSE)

proj_df <- read.csv(file.path(wd, "hypernodes", hypernode_name, "output", "ub_bounds_projected.csv"),
                    header = FALSE, col.names = c("reaction", "FBAmodel", "upper_bound"),
                    stringsAsFactors = FALSE)

# Load diet-based reaction bounds
reaction_bounds <- cfg$exchange_bounds
initial_counts <- sapply(biounit_models, function(m) m$initial_count)
names(initial_counts) <- sapply(biounit_models, function(m) m$FBAmodel)

# Convert reaction_bounds to "_r" and compute new upper bounds
reaction_bounds$reaction_r <- paste0(reaction_bounds$reaction, "_r")
reaction_bounds$background_met <- as.numeric(reaction_bounds$value)

# Create a lookup of updated upper bounds for each organism
updated_rows <- list()
for (i in seq_len(nrow(reaction_bounds))) {
  rxn <- paste0(reaction_bounds$reaction[i], "_r")
  met_val <- reaction_bounds$value[i]
  
  # Find organisms in not_proj_df using this reaction
  subset_df <- not_proj_df[not_proj_df$reaction == rxn, ]
  orgs <- unique(subset_df$FBAmodel)
  
  if (length(orgs) == 0) next  # no match → skip
  total_count <- sum(initial_counts[orgs])
  
  ub <- abs((met_val * volume) / total_count)
  
  for (org in orgs) {
    updated_rows[[length(updated_rows) + 1]] <- data.frame(
      reaction = rxn,
      FBAmodel = org,
      upper_bound = ub,
      stringsAsFactors = FALSE
    )
  }
}
diet_df <- do.call(rbind, updated_rows)

# Merge with not-projected bounds, overwrite where applicable
not_proj_df_updated <- merge(not_proj_df, diet_df,
                             by = c("reaction", "FBAmodel"), all.x = TRUE, suffixes = c("", "_diet"))

not_proj_df_updated$upper_bound <- ifelse(!is.na(not_proj_df_updated$upper_bound_diet),
                                          not_proj_df_updated$upper_bound_diet,
                                          not_proj_df_updated$upper_bound)

# Drop the temp column
not_proj_df_updated$upper_bound_diet <- NULL

# Apply replacement to _r reactions
for (i in seq_len(nrow(proj_df))) {
  rxn <- proj_df$reaction[i]
  if (grepl("_r$", rxn)) {
    base_rxn <- sub("_r$", "", rxn)
    matched_row <- reaction_bounds[reaction_bounds$reaction == base_rxn, ]
    if (nrow(matched_row) == 1) {
      proj_df$upper_bound[i] <- matched_row$background_met
    }
  }
}

# Write result
write.csv(not_proj_df_updated, file.path(wd, "hypernodes", hypernode_name, "output", "ub_bounds_not_projected.csv"),
          row.names = FALSE, quote = FALSE)

# Save the updated version if needed
write.csv(proj_df,
          file.path(wd, "hypernodes", hypernode_name, "output", "ub_bounds_projected.csv"),
          row.names = FALSE, quote = FALSE)

cat("✔ Diet-adjusted bounds saved to ub_bounds_not_projected_diet.csv\n")

source(file.path(wd, "functions_library/generate_cpp_from_arcs.R"))
generate_cpp_from_arcs(
  arcs_csv     = file.path(wd, "hypernodes", hypernode_name, "output", "repaired_arcs.csv"),
  cpp_template = file.path(wd, "functions_library/general_functions_template.cpp"),
  output_cpp = file.path(wd, "hypernodes", hypernode_name, "src", paste0("general_functions_", hypernode_name, ".cpp")),
  volume       = volume,
  cell_density = cell_density
)

source(file.path(wd, "functions_library/generate_R_from_pnpro.R"))
generate_R_from_pnpro(
  pnpro_file       = file.path(wd, "petri_nets_library", hypernode),
  r_template       = file.path(wd, "functions_library/functions_hypernode_template.R"),
  biounit_models = biounit_models,
  output_r         = file.path(wd, "hypernodes", hypernode_name, "src", paste0("functions_", hypernode_name, ".R"))
)

output_csv = file.path(wd, "hypernodes", hypernode_name, "config", "initial_data.csv")
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
writeLines(c("i; init; init.gen;"),con = output_csv)

fba_files <- sapply(biounit_models, function(m) {
  file.path(wd, "hypernodes", hypernode_name, "biounits", m$FBAmodel, m$txt_file)
})

# Then call model.generation with that vector:
model.generation(net_fname = file.path(wd, "petri_nets_library", hypernode),
                 transitions_fname = file.path(
                   wd, "hypernodes", hypernode_name, "src",
                   paste0("general_functions_", hypernode_name, ".cpp")),
                 fba_fname = fba_files)

# 1) create the 'gen' directory under this hypernode
gen_dir <- file.path(wd, "hypernodes", hypernode_name, "gen")
dir.create(gen_dir, recursive = TRUE, showWarnings = FALSE)

# 2) collect the four output filenames
out_files <- c(
  paste0(hypernode_name, ".solver"),
  paste0(hypernode_name, ".def"),
  paste0(hypernode_name, ".fbainfo"),
  paste0(hypernode_name, ".net"),
  paste0(hypernode_name, ".PlaceTransition")
)

# 3) move each file into gen/
for(f in out_files) {
  src <- file.path(wd, f)
  dst <- file.path(gen_dir, f)
  if (file.exists(src)) {
    file.rename(src, dst)
  } else {
    warning("Expected output not found: ", src)
  }
}

message("Generation outputs moved to: ", gen_dir)

model.analysis(
  solver_fname     = file.path(wd, "hypernodes", hypernode_name, "gen", paste0(hypernode_name, ".solver")),
  parameters_fname = file.path(wd, "hypernodes", hypernode_name, "config", "initial_data.csv"),
  functions_fname  = file.path(wd, "hypernodes", hypernode_name, "src", paste0("functions_", hypernode_name, ".R")),
  debug            = TRUE,
  i_time           = 0,
  f_time           = 10,
  s_time           = 1,
  atol             = 1e-06,
  rtol             = 1e-06,
  event_function   = NULL,
  fba_fname        = fba_files,
  user_files       = c(
    file.path(wd, "hypernodes", hypernode_name, "config", "population_parameters.csv"),
    file.path(wd, "hypernodes", hypernode_name, "gen", paste0(hypernode_name, ".fbainfo")),
    file.path(wd, "hypernodes", hypernode_name, "output", "ub_bounds_projected.csv"),
    file.path(wd, "hypernodes", hypernode_name, "output", "ub_bounds_not_projected.csv")
  )
)

if (COL) {
  # === Perspective on Coloring in the Framework ===
  #
  # The "coloring" step reflects the long-term goal of transitioning from
  # fully unrolled, organism-specific Petri Nets (one subnet per biounit)
  # to a **Colored Petri Net (CPN)** formalism.
  #
  # In a Colored Petri Net, a single place or transition can represent many
  # instances by assigning a "color" to each token or arc. In our context:
  #    - Colors correspond to biounits
  #    - Places like "n" and "biomass_e" would become typed with color sets
  #    - Transitions like Duplication, Death, FBA, etc., become parameterized
  #      by color, avoiding explicit duplication per biounit.
  #
  # The benefit: this abstraction compresses the visual and semantic complexity
  # of the Petri Net, making hybrid models of large microbial communities
  # much more scalable and understandable — especially during visual inspection.
  #
  # However, GreatSPN does not permit manual crafting of Colored PNPRO files easily.
  # Instead, this section assumes that:
  #    - A validated, editor-created template (CPN-style) exists
  #    - A Python layout script can intelligently fold the uncolored structure
  #      into a colored one by pattern matching and transformation
  #
  # This logic below calls such a script, passing in:
  #    - The original uncolored PNPRO
  #    - The CPN-style blank template
  #    - An output file path for the rendered Colored PNPRO
  
  # === File paths for layouting ===
  layout_script <- file.path(wd, "functions_library", "render_pnpro_colours.py")
  blank_file    <- file.path(wd, "petri_nets_library", "blank_coloured.PNPRO")  # must be valid CPN base
  input_file    <- file.path(wd, "petri_nets_library", hypernode)
  output_file   <- file.path(wd, "petri_nets_library", paste0(hypernode_name, "_coloured.PNPRO"))
  
  # === Launch the color layout generation ===
  system(paste(
    "python3",
    shQuote(layout_script),
    shQuote(blank_file),
    shQuote(input_file),
    shQuote(output_file)
  ))
  
  cat("✓ Colored layout rendered to", output_file, "\n")
}
