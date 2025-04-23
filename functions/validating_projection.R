# Initialize output lists
bacteria_files <- character()
bacteria_counts <- numeric()
projected_lower_bounds <- numeric()
projected_upper_bounds <- numeric()
shared_reactions <- list()
population_biomass <- numeric()  # Added this which was referenced but not initialized

# Enhanced tracking for reaction exclusivity
reaction_organism_map <- list()  # Maps reactions to organisms that have them
exclusive_reactions <- list()    # Will store organism-exclusive reactions
reaction_bounds_info <- list()   # Store bounds information for each reaction by organism

# Log file for tracking projectable metabolites
log_file <- file.path(output_dir_projections, "projection_log.txt")
cat("Metabolite Projection Analysis\n", file = log_file)
cat("============================\n\n", file = log_file, append = TRUE)

# Process each bacterial model
for ( i in seq_along(bacterial_models) ) {
  model <- bacterial_models[[i]]
  organism <- model$organism
  abbr <- model$abbreviation
  
  cat(sprintf("Analyzing model: %s (%s)\n", organism, abbr), file = log_file, append = TRUE)
  
  # Set file path for compiled model - updated to use FBAmodel property directly
  model_file <- file.path("compiled_models", paste0(model$txt_file, ".txt"))
  
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
  model_shared_reactions <- character()  # Initialize outside the loop
  
  for (met in projectable_metabolites) {
    related_rxns <- character()
    
    # 1. Check reaction equations if available in metadata
    if ("equation" %in% names(reaction_data)) {
      # Get boundary reactions where the metabolite appears in the equation
      eq_matches <- character()
      for (j in which(reaction_data$type == "boundary")) {
        # Check if metabolite is in the equation as a standalone term
        # This uses word boundaries to avoid partial matches
        if (grepl(paste0("\\b", met, "\\b"), reaction_data$equation[j])) {
          eq_matches <- c(eq_matches, reaction_data$abbreviation[j])
        }
      }
      related_rxns <- c(related_rxns, eq_matches)
    }
    
    # Remove duplicates
    related_rxns <- unique(related_rxns)
    
    if (length(related_rxns) > 0) {
      model_shared_reactions <- c(model_shared_reactions, related_rxns)
      cat(sprintf("  Found %d reactions for %s: %s\n", 
                  length(related_rxns), met, paste(related_rxns, collapse=", ")), 
          file = log_file, append = TRUE)
    } else {
      # If no reactions found through analysis, log this information
      cat(sprintf("  No boundary reactions found for %s through direct analysis\n", met), 
          file = log_file, append = TRUE)
    }
  }
  
  # Remove duplicates from model_shared_reactions
  model_shared_reactions <- unique(model_shared_reactions)
  
  # Update reaction-organism mapping
  for (rxn in model_shared_reactions) {
    if (is.null(reaction_organism_map[[rxn]])) {
      reaction_organism_map[[rxn]] <- character()
    }
    reaction_organism_map[[rxn]] <- c(reaction_organism_map[[rxn]], abbr)
    
    # Store bounds information
    rxn_idx <- match(rxn, reaction_data$abbreviation)
    if (!is.na(rxn_idx)) {
      lb <- reaction_data$lowbnd[rxn_idx]
      ub <- reaction_data$uppbnd[rxn_idx]
      
      if (is.null(reaction_bounds_info[[rxn]])) {
        reaction_bounds_info[[rxn]] <- list()
      }
      reaction_bounds_info[[rxn]][[abbr]] <- list(lower = lb, upper = ub)
    }
  }
  
  # Log bounds information for this model's reactions
  cat("\nReaction Bounds for Model:", abbr, "\n", file = log_file, append = TRUE)
  cat("----------------------------------\n", file = log_file, append = TRUE)
  for (rxn in model_shared_reactions) {
    rxn_idx <- match(rxn, reaction_data$abbreviation)
    if (!is.na(rxn_idx)) {
      lb <- reaction_data$lowbnd[rxn_idx]
      ub <- reaction_data$uppbnd[rxn_idx]
      cat(sprintf("  %s: Lower = %g, Upper = %g\n", rxn, lb, ub), file = log_file, append = TRUE)
    }
  }
  
  # Add to output lists if projectable metabolites were found
  if (length(projectable_metabolites) > 0) {
    bacteria_files <- c(bacteria_files, model_file)
    
    # Use default values or those from the model if available
    biomass_value <- ifelse(!is.null(model$biomass$mean), 
                            model$biomass$mean, 
                            0)
    
    bacteria_count <- ifelse(!is.null(model$initial_count), 
                             model$initial_count, 
                             0)
    
    population_biomass <- c(population_biomass, bacteria_count * biomass_value)
    shared_reactions[[abbr]] <- model_shared_reactions
    
    projected_upper_bounds <- c(projected_upper_bounds, 
                                reaction_data$uppbnd[match(model_shared_reactions, reaction_data$abbreviation)])
    projected_lower_bounds <- c(projected_lower_bounds, 
                                reaction_data$lowbnd[match(model_shared_reactions, reaction_data$abbreviation)])
  }
  
  cat("\n", file = log_file, append = TRUE)
}

# Remove duplicate entries for each abbreviation
shared_reactions <- lapply(shared_reactions, unique)

# Make reaction-organism mappings unique
reaction_organism_map <- lapply(reaction_organism_map, unique)

# Compute the joint (union) of all reactions
joint_reactions <- Reduce(union, shared_reactions)

# Identify organism-exclusive reactions
for (abbr in names(shared_reactions)) {
  # Get reactions that appear only for this organism
  exclusive <- character()
  
  for (rxn in shared_reactions[[abbr]]) {
    if (length(reaction_organism_map[[rxn]]) == 1 && 
        reaction_organism_map[[rxn]] == abbr) {
      exclusive <- c(exclusive, rxn)
    }
  }
  
  exclusive_reactions[[abbr]] <- exclusive
}

# Identify reactions shared across all organisms
all_organisms <- names(shared_reactions)
common_reactions <- character()

for (rxn in joint_reactions) {
  orgs <- reaction_organism_map[[rxn]]
  if (length(orgs) == length(all_organisms) && 
      all(sort(orgs) == sort(all_organisms))) {
    common_reactions <- c(common_reactions, rxn)
  }
}

# Only run if we found compatible models
if (length(bacteria_files) > 0) {
  # Log summary
  cat("\nSummary:\n", file = log_file, append = TRUE)
  cat("============================\n", file = log_file, append = TRUE)
  cat(sprintf("Total bacterial models analyzed: %d\n", length(bacterial_models)), file = log_file, append = TRUE)
  cat(sprintf("Compatible models found: %d\n", length(bacteria_files)), file = log_file, append = TRUE)
  cat(sprintf("Total shared boundary reactions: %d\n", length(joint_reactions)), file = log_file, append = TRUE)
  
  # Log detailed information about each reaction
  cat("\nDetailed Reaction Information:\n", file = log_file, append = TRUE)
  cat("============================\n\n", file = log_file, append = TRUE)
  
  # Log all joint reactions
  cat("All Shared Boundary Reactions:\n", file = log_file, append = TRUE)
  for (rxn in sort(joint_reactions)) {
    organisms_with_rxn <- paste(reaction_organism_map[[rxn]], collapse=", ")
    cat(sprintf("  %s: Found in organisms [%s]\n", rxn, organisms_with_rxn), file = log_file, append = TRUE)
  }
  
  # Log organism-exclusive reactions
  cat("\nOrganism-Exclusive Reactions:\n", file = log_file, append = TRUE)
  for (abbr in names(exclusive_reactions)) {
    if (length(exclusive_reactions[[abbr]]) > 0) {
      cat(sprintf("  %s exclusive reactions: %s\n", 
                  abbr, paste(exclusive_reactions[[abbr]], collapse=", ")), 
          file = log_file, append = TRUE)
    } else {
      cat(sprintf("  %s has no exclusive reactions\n", abbr), file = log_file, append = TRUE)
    }
  }
  
  # Log reactions shared across all organisms
  cat("\nReactions Shared Across All Organisms:\n", file = log_file, append = TRUE)
  if (length(common_reactions) > 0) {
    cat(sprintf("  %s\n", paste(common_reactions, collapse=", ")), file = log_file, append = TRUE)
    
    # Add bounds information for common reactions
    cat("\nBounds for Common Reactions:\n", file = log_file, append = TRUE)
    for (rxn in common_reactions) {
      cat(sprintf("  %s:\n", rxn), file = log_file, append = TRUE)
      for (abbr in names(reaction_bounds_info[[rxn]])) {
        bounds <- reaction_bounds_info[[rxn]][[abbr]]
        cat(sprintf("    %s: Lower = %g, Upper = %g\n", 
                    abbr, bounds$lower, bounds$upper), 
            file = log_file, append = TRUE)
      }
    }
  } else {
    cat("  No reactions are shared across all organisms\n", file = log_file, append = TRUE)
  }
  
  cat(sprintf("\nOutput files saved to: %s\n", output_dir_projections), file = log_file, append = TRUE)
  
  # Create output files for reaction projections
  reactions_file <- file.path(output_dir_projections, "extracted_ex_reactions.txt")
  writeLines(joint_reactions, reactions_file)
  
  # Create files for organism-exclusive reactions if any exist
  for (abbr in names(exclusive_reactions)) {
    if (length(exclusive_reactions[[abbr]]) > 0) {
      exclusive_file <- file.path(output_dir_projections, paste0("exclusive_reactions_", abbr, ".txt"))
      writeLines(exclusive_reactions[[abbr]], exclusive_file)
    }
  }
  
  # Create file for common reactions if any exist
  if (length(common_reactions) > 0) {
    common_file <- file.path(output_dir_projections, "common_reactions.txt")
    writeLines(common_reactions, common_file)
  }
  
  # Create CSV files with reaction bounds
  bounds_data <- data.frame(
    reaction = character(),
    organism = character(),
    lower_bound = numeric(),
    upper_bound = numeric(),
    stringsAsFactors = FALSE
  )
  
  row_index <- 1
  for (rxn in names(reaction_bounds_info)) {
    for (abbr in names(reaction_bounds_info[[rxn]])) {
      bounds <- reaction_bounds_info[[rxn]][[abbr]]
      bounds_data[row_index, "reaction"] <- rxn
      bounds_data[row_index, "organism"] <- abbr
      bounds_data[row_index, "lower_bound"] <- bounds$lower
      bounds_data[row_index, "upper_bound"] <- bounds$upper
      row_index <- row_index + 1
    }
  }
  
  # Write bounds data to CSV
  write.csv(bounds_data, file.path(output_dir_projections, "reaction_bounds.csv"), row.names = FALSE)
  
  # Create CSV files with FBA bound values
  for (abbr in names(shared_reactions)) {
    rxns <- shared_reactions[[abbr]]
    if (length(rxns) > 0) {
      # Create data for this organism's bounds
      org_bounds <- data.frame(
        reaction = character(),
        lower_bound = numeric(),
        upper_bound = numeric(),
        stringsAsFactors = FALSE
      )
      
      row_index <- 1
      for (rxn in rxns) {
        if (!is.null(reaction_bounds_info[[rxn]][[abbr]])) {
          bounds <- reaction_bounds_info[[rxn]][[abbr]]
          org_bounds[row_index, "reaction"] <- rxn
          org_bounds[row_index, "lower_bound"] <- bounds$lower
          org_bounds[row_index, "upper_bound"] <- bounds$upper
          row_index <- row_index + 1
        }
      }
      
      # Write organism-specific bounds to CSV
      write.csv(org_bounds, 
                file.path(output_dir_projections, paste0("EX_upper_bounds_", abbr, ".csv")), 
                row.names = FALSE)
    }
  }
  
} else {
  cat("No compatible models found for the specified metabolites.\n", file = log_file, append = TRUE)
}

# Return vectors for use in validation
get_projection_vectors <- function() {
  return(list(
    all_boundary_reactions = shared_reactions,         # All reactions by organism
    exclusive_reactions = exclusive_reactions,         # Organism-exclusive reactions
    common_reactions = common_reactions,               # Reactions shared across all organisms
    joint_reactions = joint_reactions,                 # Union of all reactions
    reaction_organism_map = reaction_organism_map,     # Which organisms have each reaction
    reaction_bounds_info = reaction_bounds_info        # Bounds information for each reaction
  ))
}