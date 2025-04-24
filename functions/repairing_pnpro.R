# ─────────────────────────────────────────────────────────────────────────────
# Auto-Repair PNPRO from Validation Log
# ─────────────────────────────────────────────────────────────────────────────

library(xml2)
library(stringr)
library(purrr)
library(jsonlite)
library(xml2)
library(stringr)
library(dplyr)
library(readr)

file_path <- paste0(wd, "/net/Minimal_EcCb.PNPRO")
log_file <- paste0(wd, "/net/validation_logs/", "Minimal_EcCb.PNPRO_validation.log")
output_path = paste0(wd, "/net/", "Minimal_EcCb.PNPRO")

# Set default output path if not provided
if (is.null(output_path)) {
  dir_name <- dirname(file_path)
  base_name <- tools::file_path_sans_ext(basename(file_path))
  output_path <- file.path(dir_name, paste0(base_name, "_repaired.PNPRO"))
}

# Read XML content
cat("Reading XML file:", file_path, "\n")
xml_content <- tryCatch(read_xml(file_path), error = function(e) {
  stop("Failed to parse XML: ", e$message)
})

# Read log file
cat("Reading validation log:", log_file, "\n")
log_lines <- readLines(log_file)

# Import arc data frame
arcs_csv_path <- file.path(dirname(log_file), paste0(tools::file_path_sans_ext(basename(file_path)), "_filtered_arcs.csv"))
cat("Reading arc data from:", arcs_csv_path, "\n")
arc_df <- read_csv(arcs_csv_path, show_col_types = FALSE)

# Extract all fixes from the validation log
fixes <- list()

# 1. Detect organism abbreviations and model files from the arc_df
# Analyze commands in arc_df to identify bacterial models and their abbreviations
command_patterns <- list(
  # Extract model file paths and organism info from FBA commands
  fba_pattern = "FBA\\[ *\"([^\"]+)\" *, *\"([^\"]+)\" *, *\\d+ *, *\"([^\"]+)\" *, *\"([^\"]+)\"",
  # Extract organism indexes from Call commands
  call_pattern = "Call\\[\"[^\"]+\", [^,]+, (\\d+)\\]"
)

# Extract all unique model files and their associated transitions
model_data <- arc_df %>%
  filter(!is.na(command) & str_detect(command, "FBA\\[")) %>%
  mutate(
    model_file = str_match(command, command_patterns$fba_pattern)[,2],
    reaction = str_match(command, command_patterns$fba_pattern)[,3],
    bacteria_count = str_match(command, command_patterns$fba_pattern)[,4],
    biomass_place = str_match(command, command_patterns$fba_pattern)[,5]
  ) %>%
  filter(!is.na(model_file)) %>%
  dplyr::select(transition, model_file, bacteria_count, biomass_place) %>%
  distinct()

organism_info <- list()
if (nrow(model_data) > 0) {
  # First, try to extract standard abbreviations from validation log
  # Look for patterns like "Missing place 'n_XX'" or "should be 'n_YY'"
  abbr_patterns <- c(
    "Missing place 'n_([a-z]+)'",
    "Missing place 'biomass_e_([a-z]+)'",
    "should be 'n_([a-z]+)'",
    "should be 'biomass_e_([a-z]+)'"
  )
  
  known_abbrs <- c()
  for (pattern in abbr_patterns) {
    matches <- str_match_all(log_lines, pattern)[[1]]
    if (nrow(matches) > 0) {
      for (i in 1:nrow(matches)) {
        if (!is.na(matches[i, 2])) {
          known_abbrs <- c(known_abbrs, matches[i, 2])
        }
      }
    }
  }
  known_abbrs <- unique(known_abbrs)
  
  # Create an association between model files and possible organism abbreviations
  model_abbr_candidates <- list()
  
  # Group model files to identify patterns
  unique_models <- unique(model_data$model_file)
  
  for (file_name in unique_models) {
    # Extract possible abbreviations from model file name
    # Try different patterns:
    # 1. Standard pattern: name_model.txt -> name
    # 2. Organism name with underscore: Organism_name -> on
    # 3. Extract first letter of each word
    
    # Pattern 1: Get text before _model.txt
    abbr1 <- str_match(file_name, "^([^_]+)_model\\.txt$")
    if (!is.na(abbr1[1,1])) {
      model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], abbr1[1,2])
    }
    
    # Pattern 2: Extract letters from organism names with underscores
    name_parts <- strsplit(tools::file_path_sans_ext(file_name), "_")[[1]]
    if (length(name_parts) > 1) {
      # Extract first letter of each part
      abbr2 <- tolower(paste0(substr(name_parts, 1, 1), collapse=""))
      model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], abbr2)
      
      # If there are only two parts, they might be genus/species - try first 2 letters of each
      if (length(name_parts) == 2) {
        abbr3 <- tolower(paste0(substr(name_parts[1], 1, 1), substr(name_parts[2], 1, 1)))
        model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], abbr3)
      }
    }
    
    # Check transitions associated with this model file
    transitions <- model_data$transition[model_data$model_file == file_name]
    
    # Try to extract abbreviations from transition names
    for (trans in transitions) {
      # Extract suffix if it follows the pattern name_suffix (e.g., Dup_ec -> ec)
      suffix_match <- str_match(trans, "_([a-z]{2,3})$")
      if (!is.na(suffix_match[1,1])) {
        model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], suffix_match[1,2])
      }
      
      # Extract prefix if it looks like a short code (e.g., ec_biomass -> ec)
      prefix_match <- str_match(trans, "^([a-z]{2,3})_")
      if (!is.na(prefix_match[1,1])) {
        model_abbr_candidates[[file_name]] <- c(model_abbr_candidates[[file_name]], prefix_match[1,2])
      }
    }
    
    # Deduplicate candidates
    model_abbr_candidates[[file_name]] <- unique(model_abbr_candidates[[file_name]])
    
    # Prioritize abbreviations that match known ones from log
    if (length(known_abbrs) > 0 && length(model_abbr_candidates[[file_name]]) > 0) {
      for (abbr in known_abbrs) {
        if (abbr %in% model_abbr_candidates[[file_name]]) {
          # This is our best guess - move it to the front
          model_abbr_candidates[[file_name]] <- c(
            abbr, 
            model_abbr_candidates[[file_name]][model_abbr_candidates[[file_name]] != abbr]
          )
        }
      }
    }
  }
  
  # Now process each model file and create organism info
  for (i in 1:nrow(model_data)) {
    file_name <- model_data$model_file[i]
    bacteria_count <- model_data$bacteria_count[i]
    biomass_place <- model_data$biomass_place[i]
    
    # Try to determine the abbreviation for this model
    abbr <- NA
    
    # First try abbreviations from the model file name
    if (!is.null(model_abbr_candidates[[file_name]]) && length(model_abbr_candidates[[file_name]]) > 0) {
      # Take the first one (highest priority)
      abbr <- model_abbr_candidates[[file_name]][1]
    } else {
      # Fallback: try to generate an abbreviation from the filename
      base_name <- tools::file_path_sans_ext(basename(file_name))
      if (grepl("_", base_name)) {
        # Take first letter of each part
        name_parts <- strsplit(base_name, "_")[[1]]
        abbr <- tolower(paste0(substr(name_parts, 1, 1), collapse=""))
      } else {
        # Take first two letters
        abbr <- tolower(substr(base_name, 1, 2))
      }
    }
    
    if (!is.na(abbr) && nchar(abbr) > 0) {
      # Create organism info entry if it doesn't exist
      if (is.null(organism_info[[abbr]])) {
        # Determine new model file name based on abbreviation
        new_model_file <- paste0(abbr, "_model.txt")
        
        # Check if current bacteria_count is already in the correct format
        new_count_place <- ifelse(
          grepl(paste0("^n_", abbr, "$"), bacteria_count), 
          bacteria_count, 
          paste0("n_", abbr)
        )
        
        # Check if current biomass_place is already in the correct format
        new_biomass_place <- ifelse(
          grepl(paste0("^biomass_e_", abbr, "$"), biomass_place), 
          biomass_place, 
          paste0("biomass_e_", abbr)
        )
        
        organism_info[[abbr]] <- list(
          old_model_file = file_name,
          new_model_file = new_model_file,
          old_count_place = bacteria_count,
          new_count_place = new_count_place,
          old_biomass_place = biomass_place,
          new_biomass_place = new_biomass_place
        )
      }
    }
  }
  
  # Fallback: If we couldn't determine organism info, use a simple counter
  if (length(organism_info) == 0) {
    for (i in seq_along(unique_models)) {
      file_name <- unique_models[i]
      # Generate a numeric abbreviation
      abbr <- sprintf("o%d", i)
      
      # Find an example transition for this model
      example_idx <- which(model_data$model_file == file_name)[1]
      bacteria_count <- model_data$bacteria_count[example_idx]
      biomass_place <- model_data$biomass_place[example_idx]
      
      organism_info[[abbr]] <- list(
        old_model_file = file_name,
        new_model_file = paste0(abbr, "_model.txt"),
        old_count_place = bacteria_count,
        new_count_place = paste0("n_", abbr),
        old_biomass_place = biomass_place,
        new_biomass_place = paste0("biomass_e_", abbr)
      )
    }
  }
}

# Debug info
cat("Detected organism abbreviations:", paste(names(organism_info), collapse=", "), "\n")
for (abbr in names(organism_info)) {
  cat("Organism", abbr, ":\n")
  cat("  Model file:", organism_info[[abbr]]$old_model_file, "->", organism_info[[abbr]]$new_model_file, "\n")
  cat("  Count place:", organism_info[[abbr]]$old_count_place, "->", organism_info[[abbr]]$new_count_place, "\n")
  cat("  Biomass place:", organism_info[[abbr]]$old_biomass_place, "->", organism_info[[abbr]]$new_biomass_place, "\n")
}

# 2. Build a mapping of all transitions and places
# This helps us standardize naming conventions

# Get all transitions from arc_df
all_transitions <- unique(c(
  arc_df$transition[arc_df$direction == "INPUT"],
  arc_df$transition[arc_df$direction == "OUTPUT"]
))

# Get all places from arc_df
all_places <- unique(arc_df$place)

# Create mapping for transitions
transition_mapping <- list()

# First, extract transition patterns from validation logs
missing_trans_pattern <- "\\[ERROR\\]\\s+Required Transitions: Missing transition: ([^ ]+)"
missing_transitions <- str_match_all(log_lines, missing_trans_pattern)[[1]][,2]

# Extract metabolite names from logs
metabolite_pattern <- "\\[ERROR|INFO\\]\\s+Shared Metabolites: (?:Missing|Found) metabolite place: ([^ ]+)"
metabolites <- unique(str_match_all(log_lines, metabolite_pattern)[[1]][,2])

# Categorize transitions from arc_df
trans_categories <- list()

# Identify Call function transitions (population dynamics)
call_transitions <- arc_df %>%
  filter(!is.na(command) & str_detect(command, "Call\\[")) %>%
  distinct(transition, command)

for (i in 1:nrow(call_transitions)) {
  trans <- call_transitions$transition[i]
  cmd <- call_transitions$command[i]
  
  # Extract function type (Starvation, Duplication, Death)
  func_match <- str_match(cmd, "Call\\[\"([^\"]+)\"")
  if (!is.na(func_match[1,2])) {
    func_type <- func_match[1,2]
    # Store in categories
    if (is.null(trans_categories[[func_type]])) {
      trans_categories[[func_type]] <- c()
    }
    trans_categories[[func_type]] <- c(trans_categories[[func_type]], trans)
  }
}

# Identify FBA transitions
fba_transitions <- arc_df %>%
  filter(!is.na(command) & str_detect(command, "FBA\\[")) %>%
  distinct(transition, command)

for (i in 1:nrow(fba_transitions)) {
  trans <- fba_transitions$transition[i]
  cmd <- fba_transitions$command[i]
  
  # Extract reaction type
  reaction_match <- str_match(cmd, "FBA\\[ *\"[^\"]+\" *, *\"([^\"]+)\"")
  if (!is.na(reaction_match[1,2])) {
    reaction <- reaction_match[1,2]
    
    # Categorize by reaction type
    reaction_type <- NA
    
    # Check if it's a biomass reaction
    if (reaction == "EX_biomass_e") {
      # Check if it's input or output (generally has "in" or "out" in the name)
      if (grepl("in", trans, ignore.case = TRUE)) {
        reaction_type <- "biomass_in"
      } else if (grepl("out", trans, ignore.case = TRUE)) {
        reaction_type <- "biomass_out"
      } else {
        reaction_type <- "biomass_other"
      }
    } else {
      # For metabolite reactions, extract the metabolite name
      # Common pattern: EX_[metabolite]_e
      metab_match <- str_match(reaction, "EX_([^_]+)_e")
      if (!is.na(metab_match[1,2])) {
        metabolite <- metab_match[1,2]
        
        # Check if it's input or output
        if (grepl("in", trans, ignore.case = TRUE)) {
          reaction_type <- paste0(metabolite, "_in")
        } else if (grepl("out", trans, ignore.case = TRUE)) {
          reaction_type <- paste0(metabolite, "_out")
        } else {
          reaction_type <- paste0(metabolite, "_other")
        }
      } else {
        # Fallback for other reaction types
        reaction_type <- reaction
      }
    }
    
    # Store in categories
    if (!is.na(reaction_type)) {
      if (is.null(trans_categories[[reaction_type]])) {
        trans_categories[[reaction_type]] <- c()
      }
      trans_categories[[reaction_type]] <- c(trans_categories[[reaction_type]], trans)
    }
  }
}

# Now create standardized naming patterns for each category
standardized_patterns <- list()

# For population dynamics
if (length(trans_categories[["Starvation"]]) > 0) standardized_patterns[["Starvation"]] <- "Starv_%s"
if (length(trans_categories[["Duplication"]]) > 0) standardized_patterns[["Duplication"]] <- "Dup_%s"
if (length(trans_categories[["Death"]]) > 0) standardized_patterns[["Death"]] <- "Death_%s"

# For FBA reactions
if (length(trans_categories[["biomass_in"]]) > 0) standardized_patterns[["biomass_in"]] <- "EX_biomass_e_in_%s"
if (length(trans_categories[["biomass_out"]]) > 0) standardized_patterns[["biomass_out"]] <- "EX_biomass_e_out_%s"

# For metabolite reactions, use standard pattern
for (metab in metabolites) {
  # Handle double underscore in glucose
  metab_clean <- gsub("__", "_", metab)
  base_name <- gsub("_e$", "", metab_clean)
  
  # Create patterns for this metabolite's reactions
  in_key <- paste0(base_name, "_in")
  out_key <- paste0(base_name, "_out")
  
  # For metabolites like "glc__D_e" -> "EX_glc__D_e_in_%s"
  if (grepl("__", metab)) {
    standardized_patterns[[in_key]] <- paste0("EX_", metab, "_in_%s")
    standardized_patterns[[out_key]] <- paste0("EX_", metab, "_out_%s")
  } else {
    # For regular metabolites
    standardized_patterns[[in_key]] <- paste0("EX_", metab, "_in_%s")
    standardized_patterns[[out_key]] <- paste0("EX_", metab, "_out_%s")
  }
}

# Check for missing transitions mentioned in validation log
for (missing_trans in missing_transitions) {
  # Extract pattern and organism
  pattern_match <- str_match(missing_trans, "^([^_]+)_(.*)$|^(.*)_([^_]+)$")
  if (!is.na(pattern_match[1,1])) {
    # Try to extract pattern and organism portion
    if (!is.na(pattern_match[1,2]) && !is.na(pattern_match[1,3])) {
      # Format is prefix_suffix
      pattern <- pattern_match[1,2]
      org <- pattern_match[1,3]
      
      standardized_patterns[[pattern]] <- paste0(pattern, "_%s")
    } else if (!is.na(pattern_match[1,4]) && !is.na(pattern_match[1,5])) {
      # Format is prefix_suffix where suffix is org
      pattern <- pattern_match[1,4]
      org <- pattern_match[1,5]
      
      standardized_patterns[[pattern]] <- paste0(pattern, "_%s")
    }
  }
}

# Now create the actual mapping for each transition
for (trans in all_transitions) {
  new_name <- trans
  
  # Try to determine which organism this transition belongs to
  best_abbr <- NA
  abbr_confidence <- 0
  
  for (abbr in names(organism_info)) {
    # Check if abbr appears at the end of the transition name
    if (grepl(paste0("_", abbr, "$"), trans)) {
      best_abbr <- abbr
      abbr_confidence <- 5
      break
    }
    
    # Check if this organism's count or biomass place is connected to this transition
    org_info <- organism_info[[abbr]]
    
    # Check if this transition appears in arc_df connected to this organism's places
    count_arcs <- sum(arc_df$transition == trans & 
                        (arc_df$place == org_info$old_count_place | 
                           arc_df$place == org_info$new_count_place))
    
    biomass_arcs <- sum(arc_df$transition == trans & 
                          (arc_df$place == org_info$old_biomass_place | 
                             arc_df$place == org_info$new_biomass_place))
    
    curr_confidence <- count_arcs + biomass_arcs
    
    if (curr_confidence > abbr_confidence) {
      best_abbr <- abbr
      abbr_confidence <- curr_confidence
    }
  }
  
  # If we couldn't determine the organism, try to extract it from the name
  if (is.na(best_abbr) || abbr_confidence == 0) {
    # Try to extract org portion from transition name
    for (abbr in names(organism_info)) {
      # Common patterns: org_reaction, reaction_org, prefix_reaction_org
      if (grepl(paste0("^", abbr, "_"), trans, ignore.case = TRUE) || 
          grepl(paste0("_", abbr, "$"), trans, ignore.case = TRUE) ||
          grepl(paste0("_", abbr, "_"), trans, ignore.case = TRUE)) {
        best_abbr <- abbr
        abbr_confidence <- 2
        break
      }
      
      # Fuzzy match: does this transition contain a string related to the organism?
      full_name <- tools::file_path_sans_ext(basename(organism_info[[abbr]]$old_model_file))
      if (grepl(full_name, trans, ignore.case = TRUE)) {
        best_abbr <- abbr
        abbr_confidence <- 1
        break
      }
    }
  }
  
  # If we found an organism, try to standardize the transition name
  if (!is.na(best_abbr)) {
    # Identify which category this transition belongs to
    matched_category <- NA
    
    for (cat in names(trans_categories)) {
      if (trans %in% trans_categories[[cat]]) {
        matched_category <- cat
        break
      }
    }
    
    # Fallback: try to guess category by looking at the name
    if (is.na(matched_category)) {
      # Try to match with standard transition types
      if (grepl("starv", trans, ignore.case = TRUE)) {
        matched_category <- "Starvation"
      } else if (grepl("dup", trans, ignore.case = TRUE)) {
        matched_category <- "Duplication"
      } else if (grepl("death", trans, ignore.case = TRUE)) {
        matched_category <- "Death"
      } else if (grepl("biomass.*in", trans, ignore.case = TRUE)) {
        matched_category <- "biomass_in"
      } else if (grepl("biomass.*out", trans, ignore.case = TRUE)) {
        matched_category <- "biomass_out"
      } else {
        # Check for metabolite reactions
        for (metab in metabolites) {
          metab_clean <- gsub("__", "_", metab)
          if (grepl(metab_clean, trans, ignore.case = TRUE)) {
            base_name <- gsub("_e$", "", metab_clean)
            if (grepl("in", trans, ignore.case = TRUE)) {
              matched_category <- paste0(base_name, "_in")
            } else if (grepl("out", trans, ignore.case = TRUE)) {
              matched_category <- paste0(base_name, "_out")
            }
            break
          }
        }
      }
    }
    
    # If we identified a category and it has a standardized pattern, apply it
    if (!is.na(matched_category) && !is.null(standardized_patterns[[matched_category]])) {
      new_name <- sprintf(standardized_patterns[[matched_category]], best_abbr)
    } else {
      # Fallback: keep current name but ensure proper organism suffix
      if (!grepl(paste0("_", best_abbr, "$"), trans)) {
        if (grepl("_[a-z]+$", trans)) {
          # Replace existing suffix
          new_name <- gsub("_[a-z]+$", paste0("_", best_abbr), trans)
        } else {
          # Add suffix
          new_name <- paste0(trans, "_", best_abbr)
        }
      }
    }
  }
  
  # If there's a name change, add to mapping
  if (new_name != trans) {
    transition_mapping[[trans]] <- new_name
  }
}

# Debug: show transition mapping
cat("Transition mappings:\n")
for (old_name in names(transition_mapping)) {
  cat(sprintf("  %s -> %s\n", old_name, transition_mapping[[old_name]]))
}

# Create mapping for places
place_mapping <- list()

# Extract metabolite standardizations from logs
metabolite_renames <- str_match_all(log_lines, "\\[ERROR\\]\\s+Shared Metabolites: Missing metabolite place: ([^ ]+)")[[1]]
existing_metabolites <- str_match_all(log_lines, "\\[INFO\\]\\s+Shared Metabolites: Found metabolite place: ([^ ]+)")[[1]]

# Combine all metabolite names for standardization checks
all_metabolites <- c()
if (nrow(metabolite_renames) > 0) {
  all_metabolites <- c(all_metabolites, metabolite_renames[,2])
}
if (nrow(existing_metabolites) > 0) {
  all_metabolites <- c(all_metabolites, existing_metabolites[,2])
}
all_metabolites <- unique(all_metabolites)

# Extract explicit place renaming suggestions from log
explicit_renames <- str_match_all(log_lines, "\\[INFO\\]\\s+Fixes: Rename place '([^']+)' to '([^']+)'")[[1]]
if (nrow(explicit_renames) > 0) {
  for (i in 1:nrow(explicit_renames)) {
    place_mapping[[explicit_renames[i,2]]] <- explicit_renames[i,3]
  }
}

# Create a map of place types based on arc connections
place_types <- list()

# Check each place connected to FBA transitions to determine its type
for (i in 1:nrow(arc_df)) {
  if (!is.na(arc_df$command[i]) && grepl("FBA\\[", arc_df$command[i])) {
    place <- arc_df$place[i]
    cmd <- arc_df$command[i]
    
    # Extract parameters from FBA command
    params <- str_match(cmd, "FBA\\[ *\"([^\"]+)\" *, *\"([^\"]+)\" *, *\\d+ *, *\"([^\"]+)\" *, *\"([^\"]+)\"")
    
    if (!is.na(params[1,1])) {
      reaction <- params[1,2]
      count_place <- params[1,4]
      biomass_place <- params[1,5]
      
      # Determine place type
      if (place == count_place) {
        place_types[[place]] <- "count"
      } else if (place == biomass_place) {
        place_types[[place]] <- "biomass"
      } else {
        # Check if this is a boundary metabolite
        for (metab in all_metabolites) {
          # Compare with and without double underscore
          metab_single <- gsub("__", "_", metab)
          
          if (place == metab || place == metab_single) {
            place_types[[place]] <- "metabolite"
            # Store the standard form (with double underscore if needed)
            place_types[[paste0(place, "_standard")]] <- metab
            break
          }
        }
      }
    }
  }
}

# Create standard form for each place
for (place in all_places) {
  new_name <- place
  
  # Check if this is a place we should standardize
  if (place %in% names(place_mapping)) {
    # Already in mapping from explicit renames
    new_name <- place_mapping[[place]]
  } else {
    # Determine place type and apply standardization
    place_type <- place_types[[place]]
    
    if (!is.null(place_type)) {
      if (place_type == "count" || place_type == "biomass") {
        # Find which organism this belongs to
        for (abbr in names(organism_info)) {
          org_info <- organism_info[[abbr]]
          
          if (place_type == "count" && place == org_info$old_count_place) {
            new_name <- org_info$new_count_place
            break
          } else if (place_type == "biomass" && place == org_info$old_biomass_place) {
            new_name <- org_info$new_biomass_place
            break
          }
        }
      } else if (place_type == "metabolite") {
        # Use standard form with proper underscore format
        standard_form <- place_types[[paste0(place, "_standard")]]
        if (!is.null(standard_form)) {
          new_name <- standard_form
        }
      }
    } else {
      # Try to infer type from place name
      
      # Check if this might be a count place
      if (grepl("^[Nn]_", place) || place == "E_coli" || place == "Clost_buty" || 
          grepl("count", place, ignore.case = TRUE)) {
        
        # Try to extract organism abbreviation
        abbr <- NA
        
        # Try to extract from standard format n_XX
        abbr_match <- str_match(place, "^[Nn]_([a-z]+)$")
        if (!is.na(abbr_match[1,1])) {
          abbr <- abbr_match[1,2]
        } else {
          # Try common patterns
          for (org_abbr in names(organism_info)) {
            if (grepl(org_abbr, place, ignore.case = TRUE)) {
              abbr <- org_abbr
              break
            }
          }
        }
        
        if (!is.na(abbr)) {
          new_name <- paste0("n_", abbr)
        }
      }
      # Check if this might be a biomass place 
      else if (grepl("biomass", place, ignore.case = TRUE)) {
        # Try to extract organism abbreviation
        abbr <- NA
        
        # Try to extract from standard format biomass_e_XX
        abbr_match <- str_match(place, "^biomass_e_([a-z]+)$")
        if (!is.na(abbr_match[1,1])) {
          abbr <- abbr_match[1,2]
        } else {
          # Try to extract from non-standard format
          for (org_abbr in names(organism_info)) {
            if (grepl(org_abbr, place, ignore.case = TRUE)) {
              abbr <- org_abbr
              break
            }
          }
        }
        
        if (!is.na(abbr)) {
          new_name <- paste0("biomass_e_", abbr)
        }
      }
      # Check if this is a metabolite that needs standardization
      else {
        # Look for underscore issues and other metabolite patterns
        for (standard_metab in all_metabolites) {
          # Check various transformations of the metabolite name
          metab_variants <- c(
            standard_metab,
            gsub("__", "_", standard_metab),  # Replace double with single
            gsub("_", "__", standard_metab)   # Replace single with double
          )
          
          if (place %in% metab_variants) {
            new_name <- standard_metab
            break
          }
        }
      }
    }
  }
  
  # Add to mapping if name changed
  if (new_name != place) {
    place_mapping[[place]] <- new_name
  }
}

# Debug: show place mapping
cat("Place mappings:\n")
for (old_name in names(place_mapping)) {
  cat(sprintf("  %s -> %s\n", old_name, place_mapping[[old_name]]))
}

# Create a mapping for command replacements
command_updates <- list()

# Extract command patterns for parsing
command_patterns <- list(
  # Extract components from FBA commands
  fba_pattern = "FBA\\[ *\"([^\"]+)\" *, *\"([^\"]+)\" *, *(\\d+) *, *\"([^\"]+)\" *, *\"([^\"]+)\"(?:, *\"([^\"]+)\")?\\]",
  # Extract components from Call commands
  call_pattern = "Call\\[\"([^\"]+)\", *([^,]+), *(\\d+)\\]"
)

# Process all commands in arc_df
for (i in 1:nrow(arc_df)) {
  if (is.na(arc_df$command[i])) next
  
  old_trans <- arc_df$transition[i]
  cmd <- arc_df$command[i]
  
  # Get new transition name if it exists
  new_trans <- ifelse(!is.null(transition_mapping[[old_trans]]), 
                      transition_mapping[[old_trans]], 
                      old_trans)
  
  # Determine which organism this transition belongs to
  trans_abbr <- NA
  for (abbr in names(organism_info)) {
    if (grepl(paste0("_", abbr, "$"), new_trans)) {
      trans_abbr <- abbr
      break
    }
  }
  
  # If we couldn't determine the organism from the new name, try the old name
  if (is.na(trans_abbr)) {
    for (abbr in names(organism_info)) {
      if (grepl(paste0("_", abbr, "$"), old_trans) || 
          grepl(abbr, old_trans, ignore.case = TRUE)) {
        trans_abbr <- abbr
        break
      }
    }
  }
  
  # Update the command based on its type
  if (grepl("^FBA\\[", cmd)) {
    # This is an FBA command
    params <- str_match(cmd, command_patterns$fba_pattern)
    
    if (!is.na(params[1,1])) {
      old_model <- params[1,2]
      reaction <- params[1,3]
      scaling <- params[1,4]
      old_count <- params[1,5]
      old_biomass <- params[1,6]
      biomass_flag <- params[1,7] # Optional parameter
      
      # Get new values
      new_model <- old_model
      new_count <- old_count
      new_biomass <- old_biomass
      
      # If we know the organism, use its standard names
      if (!is.na(trans_abbr) && !is.null(organism_info[[trans_abbr]])) {
        new_model <- organism_info[[trans_abbr]]$new_model_file
        new_count <- organism_info[[trans_abbr]]$new_count_place
        new_biomass <- organism_info[[trans_abbr]]$new_biomass_place
      } else {
        # Otherwise use mappings if available
        if (!is.null(place_mapping[[old_count]])) {
          new_count <- place_mapping[[old_count]]
        }
        if (!is.null(place_mapping[[old_biomass]])) {
          new_biomass <- place_mapping[[old_biomass]]
        }
        
        # Try to infer model file from transition
        for (abbr in names(organism_info)) {
          if (grepl(abbr, new_trans, ignore.case = TRUE)) {
            new_model <- organism_info[[abbr]]$new_model_file
            break
          }
        }
      }
      
      # Reconstruct the command
      new_cmd <- paste0('FBA[ "', new_model, '", "', reaction, '", ', scaling, ', "', 
                        new_count, '", "', new_biomass, '"')
      
      # Add biomass flag if it exists
      if (!is.na(biomass_flag)) {
        new_cmd <- paste0(new_cmd, ', "', biomass_flag, '"')
      }
      
      # Close the command
      new_cmd <- paste0(new_cmd, ']')
      
      # Add to updates if changed
      if (cmd != new_cmd) {
        command_updates[[paste(old_trans, i, sep="_")]] <- list(
          transition = old_trans,
          old_command = cmd,
          new_command = new_cmd,
          row = i
        )
      }
    }
  } else if (grepl("^Call\\[", cmd)) {
    # This is a Call command
    params <- str_match(cmd, command_patterns$call_pattern)
    
    if (!is.na(params[1,1])) {
      function_name <- params[1,2]
      parameter_expr <- params[1,3]
      org_index <- as.integer(params[1,4])
      
      # Determine if organism index needs to be updated
      new_org_index <- org_index
      
      # If we know the organism, check if its index has changed
      if (!is.na(trans_abbr)) {
        # Find the position of this organism in the organism_info list
        abbr_pos <- match(trans_abbr, names(organism_info)) - 1  # 0-based index
        
        if (!is.na(abbr_pos) && abbr_pos != org_index) {
          new_org_index <- abbr_pos
        }
      }
      
      # Update any place references in the parameter expression if needed
      # (This is a more complex case that would require parsing the parameter expression)
      
      # Reconstruct the command
      new_cmd <- paste0('Call["', function_name, '", ', parameter_expr, ', ', new_org_index, ']')
      
      # Add to updates if changed
      if (cmd != new_cmd) {
        command_updates[[paste(old_trans, i, sep="_")]] <- list(
          transition = old_trans,
          old_command = cmd,
          new_command = new_cmd,
          row = i
        )
      }
    }
  }
}
