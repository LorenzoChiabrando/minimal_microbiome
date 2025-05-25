# ============================================================
#  File: ex_bounds_module.R
#
#  Description:
#    Configures flux constraints for boundary reactions in hybrid
#    Petri net models by reading annotated metadata and population
#    dynamics. Defines:
#      - Projected bounds (modeled as places)
#      - Not-projected bounds (modeled as constraints)
# ============================================================

process_boundary_reactions <- function(
    hypernode_name,
    biounit_models,
    output_dir,
    bacteria_counts,
    projected_base_lb,
    projected_base_ub,
    not_projected_base_lb,
    not_projected_base_ub,
    projected_reactions_df
) {
  # Load and expand boundary reactions
  metadata_all <- do.call(rbind, lapply(biounit_models, function(model) {
    meta_path <- file.path("hypernodes", hypernode_name, "biounits", model$FBAmodel, "reactions_metadata.csv")
    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    meta$FBAmodel <- model$FBAmodel
    meta$txt_file <- model$txt_file
    meta
  }))
  
  metadata_all <- subset(metadata_all, type == "boundary")
  
  expanded <- do.call(rbind, lapply(seq_len(nrow(metadata_all)), function(i) {
    row <- metadata_all[i, , drop = FALSE]
    base_rxn <- row$abbreviation
    data.frame(
      reaction     = rep(base_rxn, 2),
      abbreviation = c(paste0(base_rxn, "_r"), paste0(base_rxn, "_f")),
      name         = c(paste0(row$name, " (reverse)"), paste0(row$name, " (forward)")),
      equation     = c(paste("IMPORT:", row$equation), paste("EXPORT:", row$equation)),
      subtype      = rep(row$subtype, 2),
      FBAmodel     = rep(row$FBAmodel, 2),
      txt_file     = rep(row$txt_file, 2),
      upper_bound  = c(abs(row$lowbnd), row$uppbnd),
      stringsAsFactors = FALSE
    )
  }))
  
  expanded$projected <- mapply(function(rxn, org) {
    isTRUE(any(projected_reactions_df$reaction == rxn & projected_reactions_df$FBAmodel == org))
  }, expanded$abbreviation, expanded$FBAmodel)
  
  cat("Total reactions in expanded:", nrow(expanded), "\n")
  cat("Marked as projected:", sum(expanded$projected), "\n")
  print(head(expanded[expanded$projected, ], 5))
  
  # Compute bounds
  projected <- list()
  non_projected <- list()
  
  for (i in seq_len(nrow(expanded))) {
    row <- expanded[i, ]
    if (grepl("EX_biomass_e", row$abbreviation)) next
    
    rxn_type <- if (grepl("_f$", row$abbreviation)) "f" else "r"
    
    if (row$projected) {
      upper_bound <- if (rxn_type == "r") abs(projected_base_lb) else projected_base_ub
      projected[[length(projected) + 1]] <- data.frame(
        reaction     = row$abbreviation,
        FBAmodel     = row$FBAmodel,
        upper_bound  = upper_bound,
        stringsAsFactors = FALSE
      )
    } else {
      if (rxn_type == "r") {
        if (row$subtype == "exchange") {
          # Shared resource: divide by total cell count of all organisms using it
          match_rows <- expanded[expanded$abbreviation == row$abbreviation & !expanded$projected, ]
          orgs <- unique(match_rows$FBAmodel)
          total_count <- sum(bacteria_counts[orgs])
          upper_bound <- abs(not_projected_base_lb / total_count)
        } else {
          # Organism-specific: divide only by this organism's count
          cell_count <- bacteria_counts[row$FBAmodel]
          upper_bound <- abs(not_projected_base_lb / cell_count)
        }
      } else {
        upper_bound <- not_projected_base_ub
      }
      
      non_projected[[length(non_projected) + 1]] <- data.frame(
        reaction     = row$abbreviation,
        FBAmodel     = row$FBAmodel,
        upper_bound  = upper_bound,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df_proj <- if (length(projected)) do.call(rbind, projected) else data.frame()
  df_nonproj <- if (length(non_projected)) do.call(rbind, non_projected) else data.frame()
  
  write.table(df_proj, file.path(output_dir, "ub_bounds_projected.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(df_nonproj, file.path(output_dir, "ub_bounds_not_projected.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  message("✔ Bounds generated:\n - ub_bounds_projected.csv\n - ub_bounds_not_projected.csv")
}

# ------------------------------------------------------------------
# Master function to set bounds using metadata and projection
# ------------------------------------------------------------------

run_full_ex_bounds <- function(
    hypernode_name,
    biounit_models,
    projected_base_lb,
    projected_base_ub,
    not_projected_base_lb,
    not_projected_base_ub
) {
  output_dir <- file.path(getwd(), "hypernodes", hypernode_name, "output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  bacteria_counts <- sapply(biounit_models, function(x) x$initial_count)
  names(bacteria_counts) <- sapply(biounit_models, function(x) x$FBAmodel)
  
  # Load projection info
  reaction_bounds_path <- file.path(output_dir, "reaction_bounds.csv")
  if (!file.exists(reaction_bounds_path)) stop("Missing: ", reaction_bounds_path)
  reaction_bounds_df <- read.csv(reaction_bounds_path, stringsAsFactors = FALSE)
  
  # Map abbreviations to full FBAmodel names
  abbr_map <- setNames(
    sapply(biounit_models, function(x) x$FBAmodel),
    sapply(biounit_models, function(x) x$abbreviation[2])
  )
  
  # Build _r/_f projection table with resolved FBAmodel
  projected_reactions_df <- do.call(rbind, lapply(seq_len(nrow(reaction_bounds_df)), function(i) {
    row <- reaction_bounds_df[i, ]
    full_model <- abbr_map[[row$organism]]
    if (is.null(full_model)) {
      warning("No FBAmodel found for abbreviation: ", row$organism)
      return(NULL)
    }
    data.frame(
      reaction = c(paste0(row$reaction, "_r"), paste0(row$reaction, "_f")),
      type     = c("r", "f"),
      FBAmodel = rep(full_model, 2),
      bound    = c(row$lower_bound, row$upper_bound),
      stringsAsFactors = FALSE
    )
  }))
  
  # Launch constraint assignment
  process_boundary_reactions(
    hypernode_name          = hypernode_name,
    biounit_models        = biounit_models,
    output_dir              = output_dir,
    bacteria_counts         = bacteria_counts,
    projected_base_lb       = projected_base_lb,
    projected_base_ub       = projected_base_ub,
    not_projected_base_lb   = not_projected_base_lb,
    not_projected_base_ub   = not_projected_base_ub,
    projected_reactions_df  = projected_reactions_df
  )
}
