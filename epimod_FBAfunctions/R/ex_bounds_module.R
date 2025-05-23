# ============================================================
#  File: ex_bounds_module.R
#
#  Description:
#    This R module supports the configuration of boundary constraints
#    for hybrid metabolic models by:
#      1. Extracting boundary reactions from model-specific metadata files
#      2. Generating CSV files with upper bounds for each reaction variant,
#         using cell-population-normalized scaling logic.
#
#    Convention:
#      - "_r": treated as import; assigned an upper bound.
#      - "_f": treated as export; assigned an upper bound.
#
#    For projected reactions:
#      - "_r" → ub = projected_base_lb / total_cell_count
#      - "_f" → ub = projected_base_ub
#
#    For not-projected reactions:
#      - "_r" → ub = not_projected_base_lb / total_cell_count
#      - "_f" → ub = not_projected_base_ub
# ============================================================

# ------------------------------------------------------------------
# 1) EXTRACT BOUNDARY REACTIONS
# ------------------------------------------------------------------

extract_boundary_reactions <- function(files, output_file) {
  all_reactions <- character()
  for (f in files) {
    model_dir <- dirname(f)
    metadata_file <- file.path(model_dir, "reactions_metadata.csv")
    if (!file.exists(metadata_file)) {
      warning("Metadata not found for: ", f)
      next
    }
    meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
    if (!all(c("type", "abbreviation") %in% names(meta))) {
      warning("Missing required columns in metadata for: ", f)
      next
    }
    boundary_filtered <- meta$abbreviation[meta$type == "boundary"]
    all_reactions <- c(all_reactions, boundary_filtered)
  }
  all_reactions <- sort(unique(all_reactions))
  con <- file(output_file, open = "w")
  for (rxn in all_reactions) {
    writeLines(paste0(rxn, "_r"), con = con)
    writeLines(paste0(rxn, "_f"), con = con)
  }
  close(con)
  message("Boundary reactions written to: ", output_file)
}

# ------------------------------------------------------------------
# 2) PROCESS BOUNDARY REACTIONS TO GENERATE FLUX BOUNDS
# ------------------------------------------------------------------

process_boundary_reactions <- function(
    output_file,
    bacterial_models,
    output_dir,
    bacteria_counts,
    projected_base_lb,
    not_projected_base_lb,
    projected_base_ub,
    not_projected_base_ub,
    reaction_versions = c("both", "r", "f")
) {
  reaction_versions <- match.arg(reaction_versions)
  n_bact <- length(bacterial_models)
  if (!file.exists(output_file)) stop("Reactions list not found: ", output_file)
  if (length(bacteria_counts) != n_bact) stop("Mismatch in bacteria_counts length.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  reaction_metadata_all <- do.call(rbind, lapply(bacterial_models, function(model) {
    meta_path <- file.path("hypernodes", hypernode_name, "biounits", model$FBAmodel, "reactions_metadata.csv")
    if (!file.exists(meta_path)) stop("Missing metadata for model: ", model$FBAmodel)
    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    meta$FBAmodel <- model$FBAmodel
    meta
  }))
  
  cat("✅ Total reactions in combined metadata:", nrow(reaction_metadata_all), "\n")
  cat("✅ Unique abbreviations in metadata:", length(unique(reaction_metadata_all$abbreviation)), "\n")
  print(head(unique(reaction_metadata_all$abbreviation), 10))
  
  model_names <- sapply(bacterial_models, function(m) m$FBAmodel)
  names(bacteria_counts) <- model_names
  
  reactions <- readLines(output_file, warn = FALSE)
  
  cat("✅ Total reactions read from file:", length(reactions), "\n")
  print(head(reactions, 10))
  
  get_base_name <- function(rxn) sub("(_r|_f)$", "", rxn)
  
  get_base_name <- function(rxn) sub("(_r|_f)$", "", rxn)
  base_rxns <- unique(sapply(reactions, get_base_name))
  
  matches <- base_rxns %in% reaction_metadata_all$abbreviation
  cat("✅ Reactions with metadata match:", sum(matches), "/", length(base_rxns), "\n")
  
  projected <- list()
  not_projected <- list()
  
  for (reaction in reactions) {
    if (grepl("EX_biomass_e", reaction)) next
    base_rxn <- get_base_name(reaction)
    rxn_type <- ifelse(grepl("_f$", reaction), "f", "r")
    
    if ((reaction_versions == "r" && rxn_type == "f") ||
        (reaction_versions == "f" && rxn_type == "r")) next
    
    used_models <- unique(reaction_metadata_all$FBAmodel[reaction_metadata_all$abbreviation == base_rxn])
    relevant_counts <- bacteria_counts[used_models]
    total_count <- sum(relevant_counts)
    if (total_count == 0) next
    
    if (base_rxn %in% reaction_metadata_all$abbreviation) {
      value <- if (rxn_type == "r") projected_base_lb / total_count else projected_base_ub
      projected[[length(projected) + 1]] <- data.frame(reaction = reaction, bound = value, stringsAsFactors = FALSE)
    } else {
      value <- if (rxn_type == "r") not_projected_base_lb / total_count else not_projected_base_ub
      not_projected[[length(not_projected) + 1]] <- data.frame(reaction = reaction, bound = value, stringsAsFactors = FALSE)
    }
  }
  
  df_proj <- do.call(rbind, projected)
  df_nonproj <- do.call(rbind, not_projected)
  
  # ------------------------------------------------------------------
  # Diagnostics before writing files
  # ------------------------------------------------------------------
  
  cat("📊 Diagnostics: Projected Reactions\n")
  cat("  Total:", nrow(df_proj), "\n")
  cat("  Unique:", length(unique(df_proj$reaction)), "\n")
  cat("  _r type:", sum(grepl("_r$", df_proj$reaction)), "\n")
  cat("  _f type:", sum(grepl("_f$", df_proj$reaction)), "\n")
  print(head(df_proj, 5))
  
  cat("\n📊 Diagnostics: Not Projected Reactions\n")
  cat("  Total:", nrow(df_nonproj), "\n")
  cat("  Unique:", length(unique(df_nonproj$reaction)), "\n")
  cat("  _r type:", sum(grepl("_r$", df_nonproj$reaction)), "\n")
  cat("  _f type:", sum(grepl("_f$", df_nonproj$reaction)), "\n")
  print(head(df_nonproj, 5))
  
  write.table(df_proj, file = file.path(output_dir, "ub_bounds_projected.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(df_nonproj, file = file.path(output_dir, "ub_bounds_not_projected.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  message("✔ Written:\n - ub_bounds_projected.csv\n - ub_bounds_not_projected.csv")
}

# ------------------------------------------------------------------
# 3) MASTER RUNNER FUNCTION
# ------------------------------------------------------------------

run_full_ex_bounds <- function(
    hypernode_name,
    bacterial_models,
    projected_base_lb,
    projected_base_ub,
    not_projected_base_lb,
    not_projected_base_ub,
    reaction_versions
) {
  bacteria_counts <- sapply(bacterial_models, function(x) x$initial_count)
  names(bacteria_counts) <- sapply(bacterial_models, function(x) x$FBAmodel)
  
  output_dir_projections <- file.path(getwd(), "hypernodes", hypernode_name, "output")
  output_file <- file.path(output_dir_projections, "extracted_boundary_reactions.txt")
  
  txt <- sapply(bacterial_models, function(x) x$txt_file)
  names(txt) <- sapply(bacterial_models, function(x) x$FBAmodel)
  files <- mapply(function(model, file) {
    file.path("hypernodes", hypernode_name, "biounits", model, file)
  }, names(txt), txt, USE.NAMES = TRUE)
  
  extract_boundary_reactions(files = files, output_file = output_file)
  
  process_boundary_reactions(
    output_file = output_file,
    bacterial_models = bacterial_models,
    output_dir = output_dir_projections,
    bacteria_counts = bacteria_counts,
    projected_base_lb = projected_base_lb,
    not_projected_base_lb = not_projected_base_lb,
    projected_base_ub = projected_base_ub,
    not_projected_base_ub = not_projected_base_ub,
    reaction_versions = reaction_versions
  )
}
