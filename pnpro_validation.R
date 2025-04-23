
validate_pnpro <- function(file_path, 
                           bacterial_models, 
                           metabolite_places, 
                           model_dir = "compiled_models", 
                           log_file = NULL) {
  
  # Set default log file if not provided
  if (is.null(log_file)) {
    log_dir <- file.path(dirname(file_path), "validation_logs")
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    log_file <- file.path(log_dir, paste0(basename(file_path), "_validation.log"))
  }
  
  # Initialize log file
  cat("PNPRO Validation Report\n", file = log_file)
  cat("=====================\n\n", file = log_file, append = TRUE)
  cat("File: ", file_path, "\n", file = log_file, append = TRUE)
  cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = log_file, append = TRUE)
  
  # Initialize validation counters
  validation_summary <- list(
    errors = 0,
    warnings = 0,
    info = 0,
    status = "PASS"
  )
  
  # Helper function to log issues
  log_issue <- function(level, message, section = NULL) {
    prefix <- switch(level,
                     "ERROR" = "[ERROR] ",
                     "WARNING" = "[WARNING] ",
                     "INFO" = "[INFO] ",
                     "")
    
    if (!is.null(section)) {
      message <- paste0(section, ": ", message)
    }
    
    cat(prefix, message, "\n", file = log_file, append = TRUE)
    
    # Update counters
    validation_summary[[tolower(level)]] <- validation_summary[[tolower(level)]] + 1
    
    # Update status
    if (level == "ERROR" && validation_summary$status != "FAIL") {
      validation_summary$status <- "FAIL"
    } else if (level == "WARNING" && validation_summary$status == "PASS") {
      validation_summary$status <- "WARN"
    }
  }
  
  # Helper function to log section headers
  log_section <- function(title) {
    cat("\n", title, "\n", file = log_file, append = TRUE)
    cat(paste(rep("-", nchar(title)), collapse = ""), "\n", file = log_file, append = TRUE)
  }
  
  # Get projection vectors for validating exchange reactions
  projection_data <- tryCatch({
    get_projection_vectors()
  }, error = function(e) {
    log_issue("WARNING", paste("Could not get projection vectors:", e$message), "Initialization")
    return(NULL)
  })
  
  # Extract components from PNPRO file
  xml_content <- tryCatch({
    xml2::read_xml(file_path)
  }, error = function(e) {
    log_issue("ERROR", paste("Failed to parse XML file:", e$message), "XML Parsing")
    cat("\nValidation Summary:\n", file = log_file, append = TRUE)
    cat("  Errors:", validation_summary$errors, "\n", file = log_file, append = TRUE)
    cat("  Warnings:", validation_summary$warnings, "\n", file = log_file, append = TRUE)
    cat("  Status:", validation_summary$status, "\n", file = log_file, append = TRUE)
    return(validation_summary)
  })
  
  # Extract components
  places <- xml2::xml_find_all(xml_content, "//place")
  transitions <- xml2::xml_find_all(xml_content, "//transition")
  place_names <- xml2::xml_attr(places, "name")
  transition_names <- xml2::xml_attr(transitions, "name")
  transition_delays <- xml2::xml_attr(transitions, "delay")
  
  # Log basic model stats
  log_section("Model Statistics")
  cat("Places:", length(places), "\n", file = log_file, append = TRUE)
  cat("Transitions:", length(transitions), "\n", file = log_file, append = TRUE)
  cat("Bacterial Models:", length(bacterial_models), "\n\n", file = log_file, append = TRUE)
  
  # === VALIDATION STEPS ===
  
  # Step 0: Validate metabolic networks
  log_section("1. Metabolic Network Validation")
  for (model in bacterial_models) {
    fba_model <- model$FBAmodel
    model_file <- file.path(model_dir, paste0(fba_model, ".txt"))
    
    # Check if file exists
    if (!file.exists(model_file)) {
      log_issue("ERROR", sprintf("File not found: %s", model_file), "Metabolic Network")
      next
    }
    
    # Check file format
    tryCatch({
      lines <- readLines(model_file, n = 5)
      
      # Basic content validation
      if (length(lines) < 3 || !grepl("^\\s*[A-Za-z0-9_\\-]+(\\s+[A-Za-z0-9_\\-]+)*\\s*$", lines[1])) {
        log_issue("ERROR", sprintf("Invalid reactions line format: %s", model_file), "Metabolic Network")
      }
      
      if (!grepl("^\\s*(min|max)\\s+[0-9]+\\s+[0-9]+\\s*$", lines[2])) {
        log_issue("ERROR", sprintf("Invalid optimization line format: %s", model_file), "Metabolic Network")
      }
      
      if (!grepl("^\\s*[0-9\\.eE\\+\\-]+(\\s+[0-9\\.eE\\+\\-]+)*\\s*$", lines[3])) {
        log_issue("ERROR", sprintf("Invalid objective coefficients format: %s", model_file), "Metabolic Network")
      }
      
    }, error = function(e) {
      log_issue("ERROR", sprintf("Error reading file %s: %s", model_file, e$message), "Metabolic Network")
    })
  }
  
  # Step 1: Validate place names
  log_section("2. Place Validation")
  for (model in bacterial_models) {
    abbr <- model$abbreviation
    
    # Check bacterial place
    bacterial_place <- sprintf("N_%s", abbr)
    if (!any(place_names == bacterial_place)) {
      log_issue("ERROR", sprintf("Missing bacterial place: %s", bacterial_place), "Places")
    } else {
      log_issue("INFO", sprintf("Found bacterial place: %s", bacterial_place), "Places")
    }
    
    # Check biomass place
    biomass_place <- sprintf("biomass_e_%s", abbr)
    if (!any(place_names == biomass_place)) {
      log_issue("ERROR", sprintf("Missing biomass place: %s", biomass_place), "Places")
    } else {
      log_issue("INFO", sprintf("Found biomass place: %s", biomass_place), "Places")
    }
  }
  
  # Step 2: Validate required transitions
  log_section("3. Required Transitions Validation")
  for (model in bacterial_models) {
    abbr <- model$abbreviation
    
    # Check necessary transitions
    required_transitions <- c(
      sprintf("Starv_%s", abbr),
      sprintf("Death_%s", abbr),
      sprintf("Dup_%s", abbr),
      sprintf("EX_biomass_e_out_%s", abbr),
      sprintf("EX_biomass_e_in_%s", abbr)
    )
    
    cat(sprintf("Checking required transitions for %s (%s):\n", model$organism, abbr), 
        file = log_file, append = TRUE)
    
    for (req_trans in required_transitions) {
      if (!any(transition_names == req_trans)) {
        log_issue("ERROR", sprintf("Missing transition: %s", req_trans), "Required Transitions")
      } else {
        log_issue("INFO", sprintf("Found transition: %s", req_trans), "Required Transitions")
      }
    }
  }
  
  # Step 3: Validate biomass exchange transitions regardless of naming convention
  log_section("4. Biomass Exchange Validation")
  has_biomass_in <- any(grepl("biomass_e.*in", transition_names, ignore.case = TRUE))
  has_biomass_out <- any(grepl("biomass_e.*out", transition_names, ignore.case = TRUE))
  
  if (!has_biomass_in) {
    log_issue("ERROR", "Missing biomass uptake transition (containing 'biomass_e' and 'in')", "Biomass Exchange")
  } else {
    log_issue("INFO", "Found biomass uptake transitions", "Biomass Exchange")
  }
  
  if (!has_biomass_out) {
    log_issue("ERROR", "Missing biomass secretion transition (containing 'biomass_e' and 'out')", "Biomass Exchange")
  } else {
    log_issue("INFO", "Found biomass secretion transitions", "Biomass Exchange")
  }
  
  # Step 4: For each bacterial model, validate their biomass transitions exist
  log_section("5. Organism-Specific Biomass Transitions")
  for (model in bacterial_models) {
    # Extract possible identifiers to look for in transition names
    identifiers <- c(model$organism, model$abbreviation)
    model_has_biomass <- FALSE
    
    for (id in identifiers) {
      model_biomass_pattern <- paste0(".*", id, ".*biomass_e.*|.*biomass_e.*", id, ".*")
      if (any(grepl(model_biomass_pattern, transition_names, ignore.case = TRUE))) {
        model_has_biomass <- TRUE
        break
      }
    }
    
    if (!model_has_biomass) {
      log_issue("ERROR", sprintf("Missing biomass exchange transitions for %s", model$organism), 
                "Organism Biomass")
    } else {
      log_issue("INFO", sprintf("Found biomass exchange transitions for %s", model$organism), 
                "Organism Biomass")
    }
  }
  
  # Step 5: Validate FBA transition patterns
  log_section("6. FBA Transition Pattern Validation")
  fba_pattern <- ".*EX_.*_(in|out).*|.*(in|out)_EX_.*"
  potential_fba_transitions <- transitions[grepl(fba_pattern, transition_names, ignore.case = TRUE)]
  
  if (length(potential_fba_transitions) > 0) {
    log_issue("INFO", sprintf("Found %d potential FBA transitions", length(potential_fba_transitions)), 
              "FBA Patterns")
    
    for (trans in potential_fba_transitions) {
      trans_name <- xml2::xml_attr(trans, "name")
      delay_attr <- xml2::xml_attr(trans, "delay")
      
      # Check if it uses FBA command
      if (!grepl("FBA\\[", delay_attr)) {
        log_issue("ERROR", sprintf("Transition %s does not use FBA command", trans_name), "FBA Patterns")
        next
      }
      
      # Check biomass flag for biomass transitions
      if (grepl("biomass_e", trans_name, ignore.case = TRUE)) {
        if (!grepl("\"true\"", delay_attr)) {
          log_issue("WARNING", sprintf("Biomass transition %s should have biomass flag set to 'true'", 
                                       trans_name), "FBA Patterns")
        }
      }
      
      # Check for model file reference
      model_found <- FALSE
      for (model in bacterial_models) {
        if (grepl(sprintf('"%s\\.txt"', model$FBAmodel), delay_attr)) {
          model_found <- TRUE
          break
        }
      }
      
      if (!model_found) {
        log_issue("ERROR", sprintf("Transition %s does not reference any valid model file", 
                                   trans_name), "FBA Patterns")
      }
    }
  } else {
    log_issue("WARNING", "No potential FBA transitions found", "FBA Patterns")
  }
  
  # Step 6: Validate shared metabolite places
  log_section("7. Shared Metabolite Validation")
  for (met in metabolite_places) {
    if (!any(place_names == met)) {
      log_issue("ERROR", sprintf("Missing metabolite place: %s", met), "Shared Metabolites")
    } else {
      log_issue("INFO", sprintf("Found metabolite place: %s", met), "Shared Metabolites")
    }
  }
  
  # Step 7: Validate that appropriate exchange reactions for shared metabolites are present
  log_section("8. Exchange Reaction Validation")
  
  if (!is.null(projection_data) && !is.null(projection_data$joint_reactions)) {
    for (rxn in projection_data$joint_reactions) {
      # Skip biomass reactions as they're already validated
      if (grepl("biomass", rxn, ignore.case = TRUE)) {
        next
      }
      
      # Look for transitions handling this reaction
      has_reaction_transition <- FALSE
      
      # Check both standard naming and any variations
      reaction_patterns <- c(
        paste0("^", rxn, "_(in|out)_[A-Za-z0-9]+$"),
        paste0("^[A-Za-z0-9]+_", rxn, "_(in|out)$"),
        paste0("^.*_", rxn, "_.*$")
      )
      
      for (pattern in reaction_patterns) {
        if (any(grepl(pattern, transition_names, ignore.case = TRUE))) {
          has_reaction_transition <- TRUE
          break
        }
      }
      
      # Alternative check - look in the FBA command text
      if (!has_reaction_transition) {
        for (i in seq_along(transition_delays)) {
          if (grepl(paste0('"', rxn, '"'), transition_delays[i])) {
            has_reaction_transition <- TRUE
            break
          }
        }
      }
      
      # Flag missing transitions
      if (!has_reaction_transition) {
        orgs_with_rxn <- projection_data$reaction_organism_map[[rxn]]
        
        if (length(orgs_with_rxn) > 0) {
          log_issue("ERROR", sprintf("Missing transitions for exchange reaction '%s' needed by organisms: %s", 
                                     rxn, paste(orgs_with_rxn, collapse=", ")), "Exchange Reactions")
        }
      } else {
        log_issue("INFO", sprintf("Found transitions for exchange reaction: %s", rxn), 
                  "Exchange Reactions")
      }
    }
  } else {
    log_issue("WARNING", "Projection data not available, skipping exchange reaction validation", 
              "Exchange Reactions")
  }
  
  # Step 8: Validate organism-specific reactions
  log_section("9. Organism-Specific Reaction Validation")
  
  if (!is.null(projection_data) && !is.null(projection_data$exclusive_reactions)) {
    for (abbr in names(projection_data$exclusive_reactions)) {
      org_specific_rxns <- projection_data$exclusive_reactions[[abbr]]
      
      if (length(org_specific_rxns) > 0) {
        log_issue("INFO", sprintf("Found %d exclusive reactions for %s", 
                                  length(org_specific_rxns), abbr), "Exclusive Reactions")
        
        for (rxn in org_specific_rxns) {

          # Update the code to skip reactions with multiple boundary reaction patterns
          if (!grepl("^(EX|DM|SINK)_", rxn, ignore.case = TRUE)) {
            next
          }
          
          # Look for organism-specific transitions
          has_org_specific <- FALSE
          
          # Check naming patterns
          reaction_patterns <- c(
            paste0("^", rxn, "_(in|out)_", abbr, "$"),
            paste0("^", abbr, "_", rxn, "_(in|out)$"),
            paste0("^.*", rxn, ".*", abbr, ".*$")
          )
          
          for (pattern in reaction_patterns) {
            if (any(grepl(pattern, transition_names, ignore.case = TRUE))) {
              has_org_specific <- TRUE
              break
            }
          }
          
          # Alternative check - look in the FBA command text
          if (!has_org_specific) {
            for (i in seq_along(transition_delays)) {
              if (grepl(paste0('"', rxn, '"'), transition_delays[i]) && 
                  grepl(paste0('"', abbr, '"'), transition_delays[i])) {
                has_org_specific <- TRUE
                break
              }
            }
          }
          
          # Flag missing organism-specific transitions
          if (!has_org_specific) {
            log_issue("WARNING", sprintf("Missing transitions for organism-specific exchange reaction '%s' for %s", 
                                         rxn, abbr), "Exclusive Reactions")
          } else {
            log_issue("INFO", sprintf("Found transitions for organism-specific exchange reaction: %s for %s", 
                                      rxn, abbr), "Exclusive Reactions")
          }
        }
      } else {
        log_issue("INFO", sprintf("No exclusive reactions found for %s", abbr), "Exclusive Reactions")
      }
    }
  } else {
    log_issue("WARNING", "Projection data not available, skipping organism-specific reaction validation", 
              "Exclusive Reactions")
  }
  
  # Step 9: FBA Command Syntax Validation
  log_section("10. FBA Command Syntax Validation")
  
  fba_transitions <- transitions[grepl("FBA\\[", xml2::xml_attr(transitions, "delay"))]
  
  if (length(fba_transitions) > 0) {
    log_issue("INFO", sprintf("Found %d FBA command transitions", length(fba_transitions)), 
              "FBA Syntax")
    
    for (trans in fba_transitions) {
      trans_name <- xml2::xml_attr(trans, "name")
      delay_attr <- xml2::xml_attr(trans, "delay")
      
      # Extract FBA command parameters using regex
      fba_match <- regexec('FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9.]+) *, *"([^"]+)" *, *"([^"]+)"(?: *, *"(true|false)")? *\\]', delay_attr)
      fba_parts <- regmatches(delay_attr, fba_match)
      
      if (length(fba_parts) == 0 || length(fba_parts[[1]]) < 6) {
        # Try alternative pattern with fewer parameters
        fba_match <- regexec('FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9.]+) *, *"([^"]+)" *, *"([^"]+)" *\\]', delay_attr)
        fba_parts <- regmatches(delay_attr, fba_match)
      }
      
      # Check if we could extract the FBA command parts
      if (length(fba_parts) == 0 || length(fba_parts[[1]]) < 6) {
        log_issue("ERROR", sprintf("Invalid FBA command syntax: %s", delay_attr), "FBA Syntax")
        next
      }
      
      # Extract the parameters
      parts <- fba_parts[[1]]
      bacterial_file <- parts[2]
      reaction_to_associate <- parts[3]
      scaling_constant <- as.numeric(parts[4])
      bacterial_count_place <- parts[5]
      average_biomass_place <- parts[6]
      biomass_flag <- if(length(parts) >= 7) parts[7] else "false"
      
      # Validate each parameter
      model_found <- FALSE
      for (model in bacterial_models) {
        if (bacterial_file == paste0(model$FBAmodel, ".txt")) {
          model_found <- TRUE
          break
        }
      }
      
      if (!model_found) {
        log_issue("ERROR", sprintf("Unknown model file: %s", bacterial_file), "FBA Syntax")
      }
      
      if (!grepl("^(EX|DM|SINK)_", reaction_to_associate, ignore.case = TRUE)) {
        log_issue("WARNING", sprintf("Unusual reaction naming: %s", reaction_to_associate), "FBA Syntax")
      }
      
      if (is.na(scaling_constant) || scaling_constant <= 0) {
        log_issue("ERROR", sprintf("Invalid scaling constant: %s", parts[4]), "FBA Syntax")
      }
      
      if (!(bacterial_count_place %in% place_names)) {
        log_issue("ERROR", sprintf("Non-existent bacterial count place: %s", bacterial_count_place), 
                  "FBA Syntax")
      }
      
      if (!(average_biomass_place %in% place_names)) {
        log_issue("ERROR", sprintf("Non-existent average biomass place: %s", average_biomass_place), 
                  "FBA Syntax")
      }
      
      if (grepl("biomass", reaction_to_associate, ignore.case = TRUE) && biomass_flag != "true") {
        log_issue("WARNING", sprintf("Biomass reaction should have biomass flag set to 'true': %s", 
                                     trans_name), "FBA Syntax")
      } else if (biomass_flag == "true" && !grepl("biomass", reaction_to_associate, ignore.case = TRUE)) {
        log_issue("WARNING", sprintf("Non-biomass reaction has biomass flag set to 'true': %s", 
                                     trans_name), "FBA Syntax")
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
  return(list(
    status = validation_summary$status,
    error_count = validation_summary$errors,
    warning_count = validation_summary$warnings,
    info_count = validation_summary$info,
    log_file = log_file
  ))
}