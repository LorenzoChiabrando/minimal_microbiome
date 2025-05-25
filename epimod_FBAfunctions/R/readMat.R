#' @title Read matlab model
#' @description The generic function \code{readMATmod} imports matlab cobra (ONLY) models into ....
#' @param fba_mat Full path to matlab model file
#' @details Returns ...
#' @import R.matlab stringr
#' @export

FBAmat.read = function(fba_mat, input_dir){
  
  data = R.matlab::readMat(fba_mat)
  dat.mat = data[[1]]
  
  mod.var = dimnames(dat.mat)[[1]]
  
  # 1) model name
  if( "modelID" %in% mod.var ){
    mod.id = as.character(dat.mat[[which(mod.var == "modelID")]])
  } else {
    mod.id = names(data)[1]
  }
  
  if( "modelName" %in% mod.var ){
    mod.name = as.character(dat.mat[[which(mod.var == "modelName")]])
  } else {
    mod.name = mod.id
  }
  
  if( "description" %in% mod.var ){
    mod.desc = as.character(dat.mat[[which(mod.var == "description")]])
  } else {
    mod.desc = mod.id
  }
  
  if( "rev" %in% mod.var ){
    mod.react_rev = as.vector(dat.mat[[which(mod.var=="rev")]]) == TRUE
  } else {
    mod.react_rev = as.vector(dat.mat[[which(mod.var=="lb")]]) < 0
  }
  
  #' Generate human-readable reaction equations based on stoichiometry and flux bounds
  #'
  #' @param S Stoichiometric matrix (metabolites × reactions)
  #' @param met_ids Character vector of metabolite IDs
  #' @param rxn_ids Character vector of reaction IDs
  #' @param lb Numeric vector of lower bounds for each reaction
  #' @param ub Numeric vector of upper bounds for each reaction
  #'
  #' @return A character vector of formatted equations, one per reaction
  generate_reaction_equations <- function(S, met_ids, rxn_ids, lb, ub) {
    n_rxns <- length(rxn_ids)
    equations <- character(n_rxns)
    
    for (j in seq_len(n_rxns)) {
      coefs <- S[, j]
      participating <- which(coefs != 0)
      
      # Skip if no non-zero stoichiometry
      if (length(participating) == 0) {
        equations[j] <- ""
        next
      }
      
      # ------------------------
      # Structure: Reactants and Products
      # ------------------------
      
      reactants <- participating[coefs[participating] < 0]
      products  <- participating[coefs[participating] > 0]
      
      # Format each side as a string with coefficients
      format_side <- function(indices, coefs) {
        if (length(indices) == 0) return("")
        paste(
          sapply(indices, function(i) {
            coef <- abs(coefs[i])
            if (coef == 1) met_ids[i] else paste(coef, met_ids[i])
          }),
          collapse = " + "
        )
      }
      
      lhs <- format_side(reactants, coefs)
      rhs <- format_side(products,  coefs)
      
      # ------------------------
      # Flux Constraints → Arrow Direction
      # ------------------------
      
      # The arrow used in the equation reflects the *permitted flux direction*
      # based on bounds — even if only products are present (as in EX reactions)
      
      arrow <- if (lb[j] < 0 && ub[j] > 0) {
        "<=>"
      } else if (lb[j] >= 0 && ub[j] > 0) {
        "-->"
      } else if (ub[j] <= 0 && lb[j] < 0) {
        "<--"
      } else {
        # degenerate or undefined case
        " ?? "
      }
      
      # ------------------------
      # Handle edge cases like EX reactions with only one metabolite
      # ------------------------
      
      # Examples:
      # - if lhs is empty and rhs is "ala__L_e", then write: --> ala__L_e
      # - if lhs is "ala__L_e" and rhs is empty, then: ala__L_e -->
      # - if both sides empty: "??"
      
      if (lhs == "" && rhs != "") {
        equations[j] <- paste0(arrow, " ", rhs)
      } else if (lhs != "" && rhs == "") {
        equations[j] <- paste0(lhs, " ", arrow)
      } else if (lhs == "" && rhs == "") {
        equations[j] <- paste0(rxn_ids[j], " ", arrow)
      } else {
        equations[j] <- paste(lhs, arrow, rhs)
      }
    }
    
    return(equations)
  }
  
  # Extract model components
  S <- Matrix(dat.mat[[which(mod.var == "S")]], sparse = TRUE)
  
  # Process reaction data
  react_id <- unlist(dat.mat[[which(mod.var == "rxns")]])
  react_name <- unlist(dat.mat[[which(mod.var == "rxnNames")]])
  
  # Clean reaction IDs
  react_id <- gsub("\\(", "_", react_id)
  react_id <- gsub("\\[", "_", react_id)
  react_id <- gsub("\\)", "", react_id)
  react_id <- gsub("\\]", "", react_id)
  react_id <- gsub("_c_", "_c", react_id)
  
  # Extract bounds
  lb <- as.vector(dat.mat[[which(mod.var == "lb")]])
  ub <- as.vector(dat.mat[[which(mod.var == "ub")]])
  
  # Extract objective coefficients
  obj_coef <- as.vector(dat.mat[[which(mod.var == "c")]])
  
  # Extract subsystems
  subsystems <- sapply(dat.mat[[which(mod.var == "subSystems")]], unlist)
  
  # Extract metabolite data
  met_id <- unlist(dat.mat[[which(mod.var == "mets")]])
  met_name <- unlist(dat.mat[[which(mod.var == "metNames")]])
  
  met_comp <- stringr::str_extract_all(met_id, "(?<=\\[)[a-z](?=\\])")
  
  if( all(sapply(met_comp, length) == 0 ) ){
    met_comp <- stringr::str_extract_all(met_id, "(?<=_)[a-z][0-9]?(?=$)")
  }
  
  mod.mod_compart <- unique(unlist(met_comp))
  mod.met_comp    <- match(met_comp, mod.mod_compart)
  
  original_met_ids = met_id
  
  met_id <- gsub("_", "__", met_id)
  met_id <- gsub("\\(", "_", met_id)
  met_id <- gsub("\\[", "_", met_id)
  met_id <- gsub("\\)", "", met_id)
  met_id <- gsub("\\]", "", met_id)
  
  # Generate reaction equations
  equations <- generate_reaction_equations(S, met_id, react_id, lb, ub)
  
  # Extract gene-reaction associations
  if ("rxnGeneMat" %in% mod.var) {
    gene_mat <- dat.mat[[which(mod.var == "rxnGeneMat")]]
    gene_id <- unname(unlist(dat.mat[[which(mod.var == "genes")]]))
    
    # Process GPR rules
    gpr_rules <- sapply(sapply(dat.mat[[which(mod.var == "grRules")]], unlist),
                        function(entry) {
                          if (length(entry) == 0) "" else unname(entry)
                        })
    
    # Calculate gene association status
    is_gene_associated <- apply(gene_mat != 0, 1, any)
  } else {
    # Default values if gene data not available
    gpr_rules <- rep("", length(react_id))
    is_gene_associated <- rep(FALSE, length(react_id))
  }
  
  extract_kegg_ids <- function(kegg_data) {
    result <- vector("character", length(kegg_data))
    
    for (i in seq_along(kegg_data)) {
      if (length(kegg_data[[i]][[1]]) == 0) {
        result[i] <- NA_character_
      } else {
        result[i] <- kegg_data[[i]][[1]]
      }
    }
  }
  
  # Extract KEGG IDs if available
  if ("rxnKEGGID" %in% mod.var) {
    rxn_kegg <- extract_kegg_ids(dat.mat[[which(mod.var == "rxnKEGGID")]])
  } else {
    rxn_kegg <- rep(NA, length(react_id))
  }
  
  # Direct reaction classification without helper functions
  n_rxns <- length(react_id)
  reaction_types <- character(n_rxns)
  reaction_subtypes <- character(n_rxns)
  
  for (i in seq_len(n_rxns)) {
    # Extract stoichiometric coefficients for this reaction
    stoich_col <- S[, i]
    participating_mets_indices <- which(stoich_col != 0)
    num_participants <- length(participating_mets_indices)
    
    # Calculate mass balance (sum of stoichiometric coefficients)
    mass_balance <- sum(stoich_col)
    
    # Check if reaction is mass-balanced (within numerical tolerance)
    is_balanced <- abs(mass_balance) < 1e-10
    
    # Get compartment information for participating metabolites
    if (num_participants > 0) {
      participating_met_ids <- met_id[participating_mets_indices]
      
      # Extract compartment information
      compartments <- stringr::str_extract(participating_met_ids, "\\[[a-z]\\]$")
      if (all(is.na(compartments))) {
        compartments <- stringr::str_extract(participating_met_ids, "_[a-z][0-9]?$")
        compartments <- gsub("^_", "", compartments)
      }
      if (all(is.na(compartments))) {
        compartments <- rep(NA_character_, length(participating_met_ids))
      }
      
      # Count unique non-NA compartments
      unique_compartments <- unique(compartments[!is.na(compartments)])
      spans_compartments <- length(unique_compartments) > 1
    } else {
      spans_compartments <- FALSE
    }
    
    # Check reversibility
    is_reversible <- lb[i] < 0 && ub[i] > 0
    
    # Determine if reaction is a boundary reaction
    # Boundary reactions are usually not mass balanced and have few participants
    if (!is_balanced && (num_participants <= 2)) {
      reaction_types[i] <- "boundary"
      
      # Get metabolite information for classification
      participating_met_ids <- met_id[participating_mets_indices]
      stoich_values <- stoich_col[participating_mets_indices]
      
      # Check if any metabolites are in extracellular compartment
      is_extracellular <- any(grepl("_e$|_e[0-9]$|\\[e\\]$", participating_met_ids))
      
      # Exchange reactions: involve extracellular metabolites and metabolic network boundary
      if (is_extracellular && num_participants == 1) {
        # Classic exchange reaction with extracellular metabolite
        reaction_subtypes[i] <- "exchange"
      } 
      else if (is_extracellular && num_participants > 1) {
        # Some models represent exchanges with multiple metabolites including extracellular ones
        reaction_subtypes[i] <- "exchange"
      }
      # Demand reactions: irreversible consumption with intracellular metabolites
      else if (!is_extracellular && !is_reversible && lb[i] >= 0) {
        # Check if it's a consumption process (negative coefficient for the metabolite)
        if (all(stoich_values < 0)) {
          reaction_subtypes[i] <- "demand"
        } else {
          # Production-only processes are sometimes classified as special demand reactions
          reaction_subtypes[i] <- "demand"
        }
      }
      # Sink reactions: allow both production and consumption of intracellular metabolites
      else if (!is_extracellular && is_reversible) {
        reaction_subtypes[i] <- "sink"
      }
      else {
        # Default classification for boundary reactions that don't fit the specific criteria
        # This might include nonstandard representations or special cases
        if (is_extracellular) {
          # Prioritize extracellular connection if present
          reaction_subtypes[i] <- "exchange"
        } else if (!is_reversible) {
          reaction_subtypes[i] <- "demand"
        } else {
          reaction_subtypes[i] <- "sink"
        }
      }
    } 
    else {
      # Core reactions are mass balanced
      reaction_types[i] <- "core"
      
      if (spans_compartments) {
        # Transport reactions involve multiple compartments
        reaction_subtypes[i] <- "transport"
      } 
      else if (num_participants > 1) {
        # Internal metabolic reactions
        reaction_subtypes[i] <- "internal"
      } else {
        # Default classification
        reaction_subtypes[i] <- "core"
      }
    }
  
  }
  # Create a list with the classification results
  reaction_classification <- list(
    type = reaction_types,
    subtype = reaction_subtypes
  )
  
  # Replace empty elements in the subsystems vector with NA
  subsystems <- sapply(subsystems, function(x) {
    if (length(x) == 0) NA else x
  }, simplify = TRUE, USE.NAMES = FALSE)
  
  # Create a lookup table for replacements (order by length to avoid partial matches)
  met_lookup <- data.frame(
    original = original_met_ids,
    processed = met_id,
    stringsAsFactors = FALSE
  )
  
  # Sort by length (descending) to ensure longer patterns are matched first
  met_lookup <- met_lookup[order(-nchar(met_lookup$original)), ]
  
  # Load required package
  library(stringr)
  
  # Precompute and clean metabolite patterns once
  clean_original <- gsub("\\(|\\[", "_", met_lookup$original)
  clean_original <- gsub("\\)|\\]", "", clean_original)
  pattern_keys <- paste0("\\b", clean_original, "\\b")
  replacements <- met_lookup$processed
  names(replacements) <- pattern_keys
  
  # Define all prefixes and compile a single regex
  prefixes <- c("EX_", "DM_", "SINK_", "ex_", "dm_", "sink_")
  prefix_pattern <- paste0("^(?:", paste(prefixes, collapse = "|"), ")")
  
  # Vectorized function to update all reaction IDs
  update_reaction_ids <- function(rx_ids) {
    # Extract any prefix (or NA if none)
    found_prefix <- str_extract(rx_ids, prefix_pattern)
    
    # Remove prefix for processing
    no_prefix <- str_replace(rx_ids, prefix_pattern, "")
    
    # Perform all metabolite replacements at once
    replaced <- str_replace_all(no_prefix, replacements)
    
    # Replace NA prefixes with empty string
    found_prefix[is.na(found_prefix)] <- ""
    
    # Reattach prefixes
    paste0(found_prefix, replaced)
  }
  
  # Apply to your vector of IDs
  react_id <- update_reaction_ids(react_id)
  
  # Create reactions dataframe
  reactions_df <- data.frame(
    abbreviation = react_id,
    name = react_name,
    equation = equations,
    lowbnd = lb,
    uppbnd = ub,
    obj_coef = obj_coef,
    subsystem = subsystems,
    type = reaction_classification$type,
    subtype = reaction_classification$subtype,
    gpr_rule = gpr_rules,
    is_gene_associated = is_gene_associated,
    stringsAsFactors = FALSE
  )
  
  if ("metKEGGID" %in% mod.var) {
    met_kegg <- extract_kegg_ids(dat.mat[[which(mod.var == "metKEGGID")]])
  } else {
    met_kegg <- rep(NA, length(met_id))
  }
  
  # Create a data frame to store metabolite metadata
  n_mets <- length(met_id)
  met_metadata <- data.frame(
    id = met_id,
    name = met_name,
    compartment = character(n_mets),
    is_core = logical(n_mets),
    is_boundary = logical(n_mets),
    is_external_boundary = logical(n_mets),
    is_internal_boundary = logical(n_mets),
    stringsAsFactors = FALSE
  )
  
  # Extract compartment information
  for (i in seq_len(n_mets)) {
    # Try different compartment notation patterns
    if (grepl("\\[[a-z]\\]$", met_id[i])) {
      # [c], [e], etc. notation
      met_metadata$compartment[i] <- sub(".*\\[([a-z])\\]$", "\\1", met_id[i])
    } else if (grepl("_[a-z][0-9]?$", met_id[i])) {
      # _c, _e, _c1, etc. notation
      met_metadata$compartment[i] <- sub(".*_([a-z][0-9]?)$", "\\1", met_id[i])
    } else {
      # No recognized compartment notation
      met_metadata$compartment[i] <- NA_character_
    }
  }
  
  # Initialize sets for metabolite classification
  R_core_indices <- which(reaction_types == "core")
  R_boundary_indices <- which(reaction_types == "boundary")
  R_exchange_indices <- which(reaction_subtypes == "exchange")
  R_dm_sink_indices <- which(reaction_subtypes %in% c("demand", "sink"))
  
  # Generate metabolite participation in reaction types
  C_core <- rep(FALSE, n_mets)
  C_boundary <- rep(FALSE, n_mets)
  C_boundary_ex <- rep(FALSE, n_mets)
  C_boundary_int <- rep(FALSE, n_mets)
  
  # Iterate through all metabolites
  for (i in seq_len(n_mets)) {
    # Find which reactions this metabolite participates in
    participating_reactions <- which(S[i,] != 0)
    
    # Check if metabolite participates in core reactions
    # C_core = { m ∈ C | ∃ r ∈ R_core : m ∈ supp(r) }
    C_core[i] <- any(participating_reactions %in% R_core_indices)
    
    # Check if metabolite participates in boundary reactions
    # C_b = { m ∈ C | ∃ r ∈ R_b : m ∈ supp(r) }
    C_boundary[i] <- any(participating_reactions %in% R_boundary_indices)
    
    # Check if metabolite participates in exchange reactions
    # C_b_ex = { m ∈ C_b | ∃ r ∈ R_ex : m ∈ supp(r) }
    C_boundary_ex[i] <- any(participating_reactions %in% R_exchange_indices)
    
    # Check if metabolite participates in demand/sink reactions
    # C_b_int = { m ∈ C_b | ∃ r ∈ (R_dm ∪ R_sink) : m ∈ supp(r) }
    C_boundary_int[i] <- any(participating_reactions %in% R_dm_sink_indices)
  }
  
  # Update metabolite metadata
  met_metadata$is_core <- C_core
  met_metadata$is_boundary <- C_boundary
  met_metadata$is_external_boundary <- C_boundary_ex
  met_metadata$is_internal_boundary <- C_boundary_int
  
  # Add categorical classification column
  met_metadata$primary_class <- "unclassified"
  met_metadata$primary_class[C_core & !C_boundary] <- "core_only"
  met_metadata$primary_class[!C_core & C_boundary] <- "boundary_only"
  met_metadata$primary_class[C_core & C_boundary] <- "core_and_boundary"
  
  # Add boundary subtype classification
  met_metadata$boundary_subtype <- NA_character_
  met_metadata$boundary_subtype[C_boundary_ex & !C_boundary_int] <- "external_only"
  met_metadata$boundary_subtype[!C_boundary_ex & C_boundary_int] <- "internal_only"
  met_metadata$boundary_subtype[C_boundary_ex & C_boundary_int] <- "external_and_internal"
  
  # Create a quantitative view of compartment distribution
  compartment_table <- table(met_metadata$compartment, useNA = "always")
  compartment_percent <- 100 * prop.table(compartment_table)
  
  # Display summary statistics
  cat("Metabolite Classification Summary:\n")
  cat("----------------------------------\n")
  cat("Total metabolites:", n_mets, "\n")
  cat("Core metabolites:", sum(C_core), sprintf("(%.1f%%)", 100*sum(C_core)/n_mets), "\n")
  cat("Boundary metabolites:", sum(C_boundary), sprintf("(%.1f%%)", 100*sum(C_boundary)/n_mets), "\n")
  cat("  - External boundary:", sum(C_boundary_ex), "\n")
  cat("  - Internal boundary:", sum(C_boundary_int), "\n")
  cat("Both core and boundary:", sum(C_core & C_boundary), "\n")
  cat("\nCompartment Distribution:\n")
  
  for (i in seq_along(compartment_table)) {
    comp_name <- names(compartment_table)[i]
    if (is.na(comp_name)) comp_name <- "Unknown"
    cat(sprintf("  - %s: %d (%.1f%%)\n", 
                comp_name, 
                compartment_table[i], 
                compartment_percent[i]))
  }
  
  # Calculate counts for each classification
  total_reactions <- length(reaction_types)
  core_reactions <- sum(reaction_types == "core")
  boundary_reactions <- sum(reaction_types == "boundary")
  
  # Core reaction subtypes
  transport_reactions <- sum(reaction_subtypes == "transport")
  internal_reactions <- sum(reaction_subtypes == "internal")
  other_core_reactions <- sum(reaction_types == "core" & !reaction_subtypes %in% c("transport", "internal"))
  
  # Boundary reaction subtypes
  exchange_reactions <- sum(reaction_subtypes == "exchange")
  demand_reactions <- sum(reaction_subtypes == "demand")
  sink_reactions <- sum(reaction_subtypes == "sink")
  other_boundary_reactions <- sum(reaction_types == "boundary" & 
                                    !reaction_subtypes %in% c("exchange", "demand", "sink"))
  
  # Create a table of reaction subtypes
  subtype_table <- table(reaction_subtypes, useNA = "always")
  subtype_percent <- 100 * prop.table(subtype_table)
  
  # Identify reactions spanning multiple compartments
  compartment_spanning_counts <- numeric(0)
  unique_compartment_counts <- numeric(0)
  
  # Helper function to extract compartments from metabolite IDs
  extract_compartments_from_ids <- function(met_ids) {
    compartments <- stringr::str_extract(met_ids, "\\[[a-z]\\]$")
    if (all(is.na(compartments))) {
      compartments <- stringr::str_extract(met_ids, "_[a-z][0-9]?$")
      compartments <- gsub("^_", "", compartments)
    }
    return(compartments)
  }
  
  # Loop through reactions to count compartment spans
  for (i in seq_len(total_reactions)) {
    if (reaction_types[i] == "core") {
      # Only examine core reactions
      stoich_col <- S[, i]
      participating_mets_indices <- which(stoich_col != 0)
      
      if (length(participating_mets_indices) > 0) {
        # Get compartments for participating metabolites
        participating_met_ids <- met_id[participating_mets_indices]
        compartments <- extract_compartments_from_ids(participating_met_ids)
        
        # Count unique compartments
        unique_comps <- unique(compartments[!is.na(compartments)])
        n_comps <- length(unique_comps)
        
        # Update counts
        unique_compartment_counts[n_comps] <- 
          if(is.na(unique_compartment_counts[n_comps])) 1 
        else unique_compartment_counts[n_comps] + 1
      }
    }
  }
  
  # Display summary statistics
  cat("\nReaction Classification Summary:\n")
  cat("----------------------------------\n")
  cat("Total reactions:", total_reactions, "\n\n")
  
  cat("Core reactions:", core_reactions, 
      sprintf("(%.1f%%)", 100*core_reactions/total_reactions), "\n")
  cat("  - Transport reactions:", transport_reactions, 
      sprintf("(%.1f%% of core)", 100*transport_reactions/core_reactions), "\n")
  cat("  - Internal reactions:", internal_reactions, 
      sprintf("(%.1f%% of core)", 100*internal_reactions/core_reactions), "\n")
  if (other_core_reactions > 0) {
    cat("  - Other core reactions:", other_core_reactions, 
        sprintf("(%.1f%% of core)", 100*other_core_reactions/core_reactions), "\n")
  }
  
  cat("\nBoundary reactions:", boundary_reactions, 
      sprintf("(%.1f%%)", 100*boundary_reactions/total_reactions), "\n")
  cat("  - Exchange reactions:", exchange_reactions, 
      sprintf("(%.1f%% of boundary)", 100*exchange_reactions/boundary_reactions), "\n")
  cat("  - Demand reactions:", demand_reactions, 
      sprintf("(%.1f%% of boundary)", 100*demand_reactions/boundary_reactions), "\n")
  cat("  - Sink reactions:", sink_reactions, 
      sprintf("(%.1f%% of boundary)", 100*sink_reactions/boundary_reactions), "\n")
  if (other_boundary_reactions > 0) {
    cat("  - Other boundary reactions:", other_boundary_reactions, 
        sprintf("(%.1f%% of boundary)", 100*other_boundary_reactions/boundary_reactions), "\n")
  }
  
  cat("\nReaction Subtype Distribution:\n")
  for (i in seq_along(subtype_table)) {
    subtype_name <- names(subtype_table)[i]
    if (is.na(subtype_name)) subtype_name <- "Unclassified"
    if (subtype_table[i] > 0) {
      cat(sprintf("  - %s: %d (%.1f%%)\n", 
                  subtype_name, 
                  subtype_table[i], 
                  subtype_percent[i]))
    }
  }
  
  if (length(unique_compartment_counts) > 0) {
    cat("\nCompartment Span Distribution (Core Reactions):\n")
    for (i in seq_along(unique_compartment_counts)) {
      if (!is.na(unique_compartment_counts[i]) && unique_compartment_counts[i] > 0) {
        if (i == 1) {
          cat(sprintf("  - Single compartment: %d (%.1f%%)\n", 
                      unique_compartment_counts[i], 
                      100*unique_compartment_counts[i]/core_reactions))
        } else {
          cat(sprintf("  - Spanning %d compartments: %d (%.1f%%)\n", 
                      i,
                      unique_compartment_counts[i], 
                      100*unique_compartment_counts[i]/core_reactions))
        }
      }
    }
  }
  
  # Create output directories if they don't exist
  reactions_dir <- file.path(input_dir)
  if (!dir.exists(reactions_dir)) dir.create(reactions_dir, recursive = TRUE)
  
  # Write metadata files
  cat("Writing metadata files...\n")
  reactions_file <- file.path(reactions_dir, "reactions_metadata.csv")
  write.csv(reactions_df, file = reactions_file, row.names = FALSE)
  
  metabolites_file <- file.path(reactions_dir, "metabolites_metadata.csv")
  write.csv(met_metadata, file = metabolites_file, row.names = FALSE)
  
  mod.S = S
  return(
    list(
      S = mod.S,
      lowbnd = lb,
      uppbnd = ub,
      react_id = reactions_df$abbreviation,
      met_id = met_id,
      obj_coef = dat.mat[[which(mod.var == "c")]],
      gene_assoc = reactions_df$is_gene_associated
    )
  )
}
