# Setup paths and defaults
file_path      <- paste0(wd, "/net/", model_name, ".PNPRO")
bacterial_models
metabolite_places
model_dir      <- "compiled_models"
log_file       <- NULL

# Default log file
if (is.null(log_file)) {
  log_dir  <- file.path(dirname(file_path), "validation_logs")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  log_file <- file.path(log_dir, paste0(basename(file_path), "_validation.log"))
}

# Initialize log
cat("PNPRO Validation Report
", file = log_file)
cat("=====================

", file = log_file, append = TRUE)
cat("File: ", file_path, "
", file = log_file, append = TRUE)
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "

", file = log_file, append = TRUE)

# Validation summary
validation_summary <- list(errors = 0, warnings = 0, info = 0, status = "PASS")

# Logging helpers
log_issue <- function(level, message, section = NULL) {
  prefix <- switch(level,
                   "ERROR"   = "[ERROR] ",
                   "WARNING" = "[WARNING] ",
                   "INFO"    = "[INFO] ",
                   "")
  if (!is.null(section)) message <- paste0(section, ": ", message)
  cat(prefix, message, "\n", file = log_file, append = TRUE)
  validation_summary[[tolower(level)]] <- validation_summary[[tolower(level)]] + 1
  if (level == "ERROR" && validation_summary$status != "FAIL") validation_summary$status <- "FAIL"
  if (level == "WARNING" && validation_summary$status == "PASS") validation_summary$status <- "WARN"
}
log_section <- function(title) {
  cat("\n", title, "\n", file = log_file, append = TRUE)
  cat(paste(rep("-", nchar(title)), collapse = ""), "\n", file = log_file, append = TRUE)
}

# Normalize names for case-insensitive matching
# Note: trimws to remove any leading/trailing whitespace
normalize <- function(x) tolower(trimws(x))

# Fetch and parse XML
projection_data <- tryCatch(get_projection_vectors(), error = function(e) {
  log_issue("WARNING", paste("Could not get projection vectors:", e$message), "Initialization")
  NULL
})
xml_content <- tryCatch(xml2::read_xml(file_path), error = function(e) {
  log_issue("ERROR", paste("Failed to parse XML file:", e$message), "XML Parsing")
  return(validation_summary)
})

# Extract components
places      <- xml2::xml_find_all(xml_content, "//place")
transitions <- xml2::xml_find_all(xml_content, "//transition")
place_names <- xml2::xml_attr(places, "name")
transition_names <- xml2::xml_attr(transitions, "name")
transition_delays <- xml2::xml_attr(transitions, "delay")
# Clean vectors for insensitive compare
place_names_lc      <- normalize(place_names)
transition_names_lc <- normalize(transition_names)

# Step counters
step_num <- 0

# === Validation Steps ===

# Step 0: Metabolic Network Validation
step_num <- step_num + 1
log_section(sprintf("%d. Metabolic Network Validation", step_num))
# Check each FBA model file directly
for (model in bacterial_models) {
  fba_file <- sprintf("%s.txt", model$txt_file)
  fba_path <- file.path(model_dir, fba_file)
  if (!file.exists(fba_path)) {
    log_issue("ERROR", sprintf("Missing FBA model file: %s", fba_path), "Metabolic Network")
  } else {
    log_issue("INFO", sprintf("Found FBA model file: %s", fba_path), "Metabolic Network")
  }
}
# Also verify every FBA transition references an existing .txt
fba_indices <- grepl("^FBA\\[", transition_delays)
for (i in which(fba_indices)) {
  trans_name <- transition_names[i]
  delay_attr <- transition_delays[i]
  m <- regmatches(delay_attr, regexec('FBA\\[\\s*"([^".]+\\.txt)"', delay_attr))[[1]]
  if (length(m) >= 2) {
    fba_path <- file.path(model_dir, m[2])
    if (!file.exists(fba_path)) {
      log_issue("ERROR",
                sprintf("Transition '%s': referenced FBA file not found: %s", trans_name, fba_path),
                "Metabolic Network")
    }
  } else {
    log_issue("WARNING",
              sprintf("Transition '%s': cannot parse FBA filename: %s", trans_name, delay_attr),
              "Metabolic Network")
  }
}

# Step 1: Place Validation
step_num <- step_num + 1
log_section(sprintf("%d. Place Validation", step_num))
for (model in bacterial_models) {
  abbr <- normalize(model$abbreviation)
  log_issue("INFO", sprintf("Validating places for %s", model$organism), "Places")
  expected <- c(sprintf("n_%s", abbr), sprintf("biomass_e_%s", abbr))
  for (place in expected) {
    if (!(place %in% place_names_lc)) {
      log_issue("ERROR", sprintf("Missing place '%s'", place), "Places")
    }
  }
  # Warn about near matches
  matches <- grep(abbr, place_names_lc, value = TRUE)
  extras  <- setdiff(matches, expected)
  if (length(extras) > 0) {
    log_issue("WARNING",
              sprintf("Unexpected places related to '%s': %s", abbr, paste(extras, collapse=", ")),
              "Places")
  }
}

# Step 2: Required Transitions Validation
step_num <- step_num + 1
log_section(sprintf("%d. Required Transitions Validation", step_num))
for (model in bacterial_models) {
  abbr <- normalize(model$abbreviation)
  log_issue("INFO", sprintf("Validating transitions for %s", model$organism), "Transitions")
  required <- c(
    sprintf("starv_%s", abbr),
    sprintf("death_%s", abbr),
    sprintf("dup_%s", abbr),
    sprintf("ex_biomass_e_out_%s", abbr),
    sprintf("ex_biomass_e_in_%s", abbr)
  )
  # Vectorized check
  exists <- required %in% transition_names_lc
  for (j in seq_along(required)) {
    level <- if (exists[j]) "INFO" else "ERROR"
    log_issue(level,
              sprintf("%s transition: %s",
                      if (exists[j]) "Found" else "Missing",
                      required[j]),
              "Required Transitions")
  }
  # Unexpected near matches
  near <- grep(abbr, transition_names_lc, value = TRUE)
  extras <- setdiff(near, required)
  if (length(extras) > 0) {
    log_issue("WARNING",
              sprintf("Unexpected transitions related to '%s': %s", abbr, paste(extras, collapse=", ")),
              "Required Transitions")
  }
}
# Step 3: Validate FBA transition patterns
step_num <- step_num + 1
log_section(sprintf("%d. FBA Transition Pattern Validation", step_num))

# Match any transition whose name starts with EX_, DM_ or SINK_ and mentions in/out
fba_pattern <- "(EX|DM|SINK)_[A-Za-z0-9_]+_(in|out)"
candidates   <- transitions[
  grepl(fba_pattern, transition_names, ignore.case = TRUE)
]
names_lc   <- normalize(transition_names)

if (length(candidates) > 0) {
  log_issue("INFO", sprintf("Found %d boundary reaction transitions", length(candidates)),
            "FBA Patterns")
  
  for (trans in candidates) {
    trans_name <- xml2::xml_attr(trans, "name")
    delay_attr <- xml2::xml_attr(trans, "delay")
    section    <- "FBA Patterns"
    
    # 1) Must invoke FBA[...]
    if (!grepl("^FBA\\[", delay_attr)) {
      log_issue("ERROR",
                sprintf("Transition '%s' does not use an FBA command", trans_name),
                section)
      next
    }
    
    # 2) Biomass‐flag for biomass reactions
    if (grepl("biomass_e", trans_name, ignore.case = TRUE) &&
        !grepl("\"true\"", delay_attr)) {
      log_issue("WARNING",
                sprintf("Biomass transition '%s' should include a biomass-flag = 'true'", trans_name),
                section)
    }
    
    # 3) Model file reference (*.txt) must exist
    m <- regmatches(delay_attr, regexec('FBA\\[\\s*"([^"]+\\.txt)"', delay_attr))[[1]]
    if (length(m) >= 2) {
      fba_file <- m[2]
      fba_path <- file.path(model_dir, fba_file)
      if (!file.exists(fba_path)) {
        log_issue("ERROR",
                  sprintf("Transition '%s': referenced model file not found: %s",
                          trans_name, fba_path),
                  section)
      }
    } else {
      log_issue("ERROR",
                sprintf("Transition '%s': could not parse model filename in delay: %s",
                        trans_name, delay_attr),
                section)
    }
  }
} else {
  log_issue("WARNING", "No boundary reaction transitions matching EX_/DM_/SINK_..._(in|out) found",
            "FBA Patterns")
}

# Step 4: Shared Metabolite Validation
step_num <- step_num + 1
log_section(sprintf("%d. Shared Metabolite Validation", step_num))

# Normalize metabolite place names for case‐insensitive compare
metabolite_places_lc <- normalize(metabolite_places)

for (i in seq_along(metabolite_places)) {
  orig_met <- metabolite_places[i]
  met_lc   <- metabolite_places_lc[i]
  
  if (!(met_lc %in% place_names_lc)) {
    log_issue(
      "ERROR",
      sprintf("Missing metabolite place: %s", orig_met),
      "Shared Metabolites"
    )
  } else {
    log_issue(
      "INFO",
      sprintf("Found metabolite place: %s", orig_met),
      "Shared Metabolites"
    )
  }
  
  # Warn about near‐matches (typos, case differences, extra underscores)
  near <- grep(met_lc, place_names_lc, value = TRUE)
  extras <- setdiff(near, met_lc)
  if (length(extras) > 0) {
    log_issue(
      "WARNING",
      sprintf("Unexpected variant(s) for '%s': %s", orig_met, paste(extras, collapse = ", ")),
      "Shared Metabolites"
    )
  }
}

# Step 9: FBA Command Syntax Validation
step_num <- step_num + 1
log_section(sprintf("%d. FBA Command Syntax Validation", step_num))
fba_idx <- grepl("^FBA\\[", transition_delays)
if (any(fba_idx)) {
  log_issue("INFO", sprintf("Found %d FBA command transitions", sum(fba_idx)), "FBA Syntax")
  for (i in which(fba_idx)) {
    trans_name <- transition_names[i]
    delay_attr <- transition_delays[i]
    section    <- "FBA Syntax"
    
    fba_match <- regexec('FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9.]+) *, *"([^"]+)" *, *"([^"]+)" *\\]', delay_attr)
    parts     <- regmatches(delay_attr, fba_match)[[1]]
    
    # Validate match length
    if (length(parts) < 6) {
      log_issue("ERROR", sprintf("Invalid FBA syntax in '%s': %s", trans_name, delay_attr), section)
      next
    }
    
    # Extract parameters
    bacterial_file        <- parts[2]
    reaction_to_associate <- parts[3]
    scaling_constant      <- as.numeric(parts[4])
    bacterial_count_place <- parts[5]
    average_biomass_place <- parts[6]
    biomass_flag          <- if (length(parts) >= 7) parts[7] else "false"
    
    # 1) Validate model file exists
    model_found <- any(sapply(bacterial_models, function(m) bacterial_file == paste0(m$FBAmodel, ".txt")))
    if (!model_found) {
      log_issue("ERROR", sprintf("Unknown model file '%s' in %s", bacterial_file, trans_name), section)
    }
    
    # 2) Reaction naming
    if (!grepl("^(EX|DM|SINK)_", reaction_to_associate, ignore.case = TRUE)) {
      log_issue("WARNING", sprintf("Unusual reaction name '%s' in %s", reaction_to_associate, trans_name), section)
    }
    
    # 3) Scaling constant
    if (is.na(scaling_constant) || scaling_constant <= 0) {
      log_issue("ERROR", sprintf("Non-positive scaling constant '%s' in %s", parts[4], trans_name), section)
    }
    
    # 4) Place existence
    for (place in c(bacterial_count_place, average_biomass_place)) {
      if (!(place %in% place_names)) {
        log_issue("ERROR", sprintf("Non-existent place '%s' in %s", place, trans_name), section)
      }
    }
    
    # 5) Biomass flag consistency
    if (grepl("biomass", reaction_to_associate, ignore.case = TRUE) && biomass_flag != "true") {
      log_issue("WARNING", sprintf("Missing biomassFlag='true' for biomass reaction in %s", trans_name), section)
    } else if (biomass_flag == "true" && !grepl("biomass", reaction_to_associate, ignore.case = TRUE)) {
      log_issue("WARNING", sprintf("biomassFlag=true on non-biomass reaction %s", trans_name), section)
    }
  }
} else {
  log_issue("WARNING", "No FBA command transitions found", "FBA Syntax")
}

# Step 10: Call Function Parameter Validation
log_section("11. Call Function Parameter Validation")

call_transitions <- transitions[grepl("Call\\[", xml2::xml_attr(transitions, "delay"))]

if (length(call_transitions) > 0) {
  log_issue("INFO", sprintf("Found %d Call function transitions", length(call_transitions)), 
            "Call Functions")
  
  # Check if Bacteria_Parameters.csv exists
  bacteria_params_file <- file.path(str_remove(dirname(file_path), "/net"), "input", "Bacteria_Parameters.csv")
  if ( !file.exists(bacteria_params_file) ) {
    log_issue("ERROR", "Bacteria_Parameters.csv file not found", "Call Functions")
  } else {
    # Load bacteria parameters
    bacteria_params <- tryCatch({
      read.csv(bacteria_params_file, header = FALSE)
      
      # Log parameters dimensions
      log_issue("INFO", sprintf("Bacteria_Parameters.csv: %d rows, %d columns", 
                                nrow(bacteria_params), ncol(bacteria_params)), "Call Functions")
      
      # Check dimensions match bacterial models
      if (nrow(bacteria_params) < length(bacterial_models)) {
        log_issue("ERROR", sprintf("Bacteria_Parameters.csv has fewer rows (%d) than bacterial models (%d)", 
                                   nrow(bacteria_params), length(bacterial_models)), "Call Functions")
      }
      
      # Process each Call function
      for (trans in call_transitions) {
        trans_name <- xml2::xml_attr(trans, "name")
        delay_attr <- xml2::xml_attr(trans, "delay")
        
        # Extract function call parameters using regex
        call_match <- regexec('Call\\["([^"]+)",\\s*FromTable\\["Bacteria_Parameters.csv",\\s*([0-9]+),\\s*([0-9]+)\\],\\s*([0-9]+)\\]', delay_attr)
        call_parts <- regmatches(delay_attr, call_match)
        
        if (length(call_parts) == 0 || length(call_parts[[1]]) < 5) {
          log_issue("ERROR", sprintf("Invalid Call function syntax: %s", delay_attr), "Call Functions")
          next
        }
        
        # Extract the parameters
        parts <- call_parts[[1]]
        function_name <- parts[2]
        row_index <- as.numeric(parts[3]) + 1  # Convert to 1-indexed for R
        col_index <- as.numeric(parts[4]) + 1  # Convert to 1-indexed for R
        organism_index <- as.numeric(parts[5])
        
        # Validate each parameter
        valid_functions <- c("Death", "Duplication", "Starvation")
        if (!(function_name %in% valid_functions)) {
          log_issue("WARNING", sprintf("Unknown function name: %s", function_name), "Call Functions")
        }
        
        if (row_index <= 0 || row_index > nrow(bacteria_params)) {
          log_issue("ERROR", sprintf("Invalid row index %d: Bacteria_Parameters.csv has %d rows", 
                                     row_index-1, nrow(bacteria_params)), "Call Functions")
        }
        
        if (col_index <= 0 || col_index > ncol(bacteria_params)) {
          log_issue("ERROR", sprintf("Invalid column index %d: Bacteria_Parameters.csv has %d columns", 
                                     col_index-1, ncol(bacteria_params)), "Call Functions")
        }
        
        if (organism_index < 0 || organism_index >= length(bacterial_models)) {
          log_issue("ERROR", sprintf("Invalid organism index %d: There are %d bacterial models", 
                                     organism_index, length(bacterial_models)), "Call Functions")
        } else {
          # Check organism-transition consistency if valid index
          expected_abbr <- bacterial_models[[organism_index + 1]]$abbreviation
          
          if (grepl("_", trans_name)) {
            parts <- strsplit(trans_name, "_")[[1]]
            if (length(parts) > 1) {
              abbr_in_name <- parts[length(parts)]
              
              if (abbr_in_name != expected_abbr) {
                log_issue("WARNING", sprintf("Organism index mismatch: %s uses index %d (%s) but name suggests %s", 
                                             trans_name, organism_index, expected_abbr, abbr_in_name), 
                          "Call Functions")
              }
            }
          }
        }
        
        # Check parameter values if indices are valid
        if (row_index <= nrow(bacteria_params) && col_index <= ncol(bacteria_params)) {
          param_value <- bacteria_params[row_index, col_index]
          
          if (function_name == "Death" && (is.na(param_value) || param_value <= 0)) {
            log_issue("WARNING", sprintf("Death rate should be positive: %g", param_value), 
                      "Call Functions")
          } else if (function_name == "Duplication" && (is.na(param_value) || param_value <= 0)) {
            log_issue("WARNING", sprintf("Duplication rate should be positive: %g", param_value), 
                      "Call Functions")
          } else if (function_name == "Starvation" && (is.na(param_value) || param_value <= 0)) {
            log_issue("WARNING", sprintf("Starvation rate should be positive: %g", param_value), 
                      "Call Functions")
          }
        }
      }
    }, error = function(e) {
      log_issue("ERROR", sprintf("Error reading Bacteria_Parameters.csv: %s", e$message), 
                "Call Functions")
      return(NULL)
    })
  }
} else {
  log_issue("WARNING", "No Call function transitions found", "Call Functions")
}

# Step X: Validate shared metabolite places
log_section("12. Shared Metabolite Validation") 
for (met in metabolite_places) {
  if (!any(place_names == met)) {
    log_issue("ERROR", sprintf("Missing metabolite place: %s", met), "Metabolites")
  } else {
    log_issue("INFO", sprintf("Found metabolite place: %s", met), "Metabolites")
  }
}

# Step 12: Validate place-transition connectivity
log_section("12. Place-Transition Connectivity Validation")

# Extract edges (arcs) from the PNPRO file
edges <- xml2::xml_find_all(xml_content, "//arc")
edge_count <- length(edges)

log_issue("INFO", sprintf("Found %d arcs in the Petri net", edge_count), "Connectivity")

# Create lookup tables for place-transition connections
input_arcs <- list()  # Transitions with arcs coming from specific places
output_arcs <- list() # Transitions with arcs going to specific places

# Populate the lookup tables
for (edge in edges) {
  head <- xml2::xml_attr(edge, "head")  # Destination (transition)
  tail <- xml2::xml_attr(edge, "tail")  # Source (place)
  kind <- xml2::xml_attr(edge, "kind")  # Type of arc (INPUT or OUTPUT)
  
  if (kind == "INPUT") {
    # This is an arc from a place to a transition
    if (is.null(input_arcs[[head]])) {
      input_arcs[[head]] <- character()
    }
    input_arcs[[head]] <- c(input_arcs[[head]], tail)
  } else if (kind == "OUTPUT") {
    # This is an arc from a transition to a place
    if (is.null(output_arcs[[head]])) {
      output_arcs[[head]] <- character()
    }
    output_arcs[[head]] <- c(output_arcs[[head]], tail)
  }
}

# Identify boundary transitions based on naming conventions
boundary_trans <- transition_names[grepl("(EX|DM|SINK)_.*_(in|out)_[A-Za-z0-9]+$|[A-Za-z0-9]+_(EX|DM|SINK)_.*_(in|out)$", transition_names, ignore.case = TRUE)]

if (length(boundary_trans) > 0) {
  log_issue("INFO", sprintf("Found %d boundary reaction transitions to validate", length(boundary_trans)), "Connectivity")
  
  for (trans_name in boundary_trans) {
    # Parse the transition name to understand its components
    # Try to identify: reaction, direction (in/out), and organism
    
    # Pattern 1: PREFIX_metabolite_direction_organism
    pattern1 <- regexec("((EX|DM|SINK)_[A-Za-z0-9_]+)_(in|out)_([A-Za-z0-9]+)$", trans_name, ignore.case = TRUE)
    components1 <- regmatches(trans_name, pattern1)
    
    # Pattern 2: organism_PREFIX_metabolite_direction
    pattern2 <- regexec("([A-Za-z0-9]+)_((EX|DM|SINK)_[A-Za-z0-9_]+)_(in|out)$", trans_name, ignore.case = TRUE)
    components2 <- regmatches(trans_name, pattern2)
    
    # Use whichever pattern matched
    if (length(components1[[1]]) >= 5) {
      reaction <- components1[[1]][2]  # PREFIX_metabolite
      direction <- components1[[1]][4]  # in or out
      organism <- components1[[1]][5]  # organism abbreviation
    } else if (length(components2[[1]]) >= 5) {
      organism <- components2[[1]][2]  # organism abbreviation
      reaction <- components2[[1]][3]  # PREFIX_metabolite
      direction <- components2[[1]][4]  # in or out
    } else {
      log_issue("WARNING", sprintf("Boundary transition naming doesn't match expected patterns: %s", trans_name), "Connectivity")
      next
    }
    
    # Instead of assuming the metabolite name is just the reaction name minus prefix,
    # we check what metabolite places are actually connected to this transition
    
    # Get connected metabolite places
    connected_metabolite_places <- character()
    
    if (direction == "in") {
      # For "in" transitions, check input connections (metabolites consumed)
      if (!is.null(input_arcs[[trans_name]])) {
        connected_metabolite_places <- input_arcs[[trans_name]][input_arcs[[trans_name]] %in% metabolite_places]
      }
    } else if (direction == "out") {
      # For "out" transitions, check output connections (metabolites produced)
      if (!is.null(output_arcs[[trans_name]])) {
        connected_metabolite_places <- output_arcs[[trans_name]][output_arcs[[trans_name]] %in% metabolite_places]
      }
    }
    
    # Check if we found any connected metabolite places
    if (length(connected_metabolite_places) == 0) {
      log_issue("ERROR", sprintf("Boundary transition %s (%s) is not properly connected to any metabolite place", 
                                 trans_name, direction), "Connectivity")
    } else if (length(connected_metabolite_places) > 1) {
      log_issue("WARNING", sprintf("Boundary transition %s (%s) is connected to multiple metabolite places: %s", 
                                   trans_name, direction, paste(connected_metabolite_places, collapse=", ")), 
                "Connectivity")
    } else {
      log_issue("INFO", sprintf("Boundary transition %s (%s) is correctly connected to metabolite place: %s", 
                                trans_name, direction, connected_metabolite_places[1]), "Connectivity")
    }
    
    # Check biomass place connectivity
    biomass_place <- sprintf("biomass_e_%s", organism)
    if (direction == "in") {
      # For "in" transitions, check for output to biomass
      if (is.null(output_arcs[[trans_name]]) || !any(output_arcs[[trans_name]] == biomass_place)) {
        log_issue("ERROR", sprintf("Transition %s should have an output arc TO place %s", 
                                   trans_name, biomass_place), "Connectivity")
      } else {
        log_issue("INFO", sprintf("Transition %s correctly outputs to place %s", 
                                  trans_name, biomass_place), "Connectivity")
      }
    } else if (direction == "out") {
      # For "out" transitions, check for input from biomass
      if (is.null(input_arcs[[trans_name]]) || !any(input_arcs[[trans_name]] == biomass_place)) {
        log_issue("ERROR", sprintf("Transition %s should have an input arc FROM place %s", 
                                   trans_name, biomass_place), "Connectivity")
      } else {
        log_issue("INFO", sprintf("Transition %s correctly has input from place %s", 
                                  trans_name, biomass_place), "Connectivity")
      }
    }
    
    # Check for bacterial count place connectivity
    bacterial_place <- sprintf("N_%s", organism)
    if ((is.null(input_arcs[[trans_name]]) || !any(input_arcs[[trans_name]] == bacterial_place)) &&
        (is.null(output_arcs[[trans_name]]) || !any(output_arcs[[trans_name]] == bacterial_place))) {
      log_issue("WARNING", sprintf("Transition %s might be missing connectivity with bacterial count place %s", 
                                   trans_name, bacterial_place), "Connectivity")
    }
  }
} else {
  log_issue("WARNING", "No boundary transitions found with standard naming patterns", "Connectivity")
}

# Rest of the validation remains the same...

# Return summary object
validation_test = list(
  status = validation_summary$status,
  error_count = validation_summary$errors,
  warning_count = validation_summary$warnings,
  info_count = validation_summary$info,
  log_file = log_file
)
