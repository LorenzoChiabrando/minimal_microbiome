# ─────────────────────────────────────────────────────────────────────────────
# Setup paths and defaults
# ─────────────────────────────────────────────────────────────────────────────

file_path <- file.path(wd, "net", paste0(model_name, ".PNPRO"))
model_dir <- "compiled_models"

# Default log file
log_dir  <- file.path(dirname(file_path), "validation_logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, paste0(basename(file_path), "_validation.log"))

# Helper: centralized logging
log_msg <- function(txt)    cat(txt, file = log_file, append = TRUE)
log_section <- function(title){
  log_msg("\n"    )
  log_msg(title )
  log_msg("\n"    )
  log_msg(str_dup("-", nchar(title)))
  log_msg("\n\n")
}

# Function to create log subsections
log_subsection <- function(title) {
  log_msg("\n  ")
  log_msg(title)
  log_msg("\n  ")
  log_msg(str_dup("-", nchar(title)))
  log_msg("\n")
}

log_issue <- function(level, message, section = NULL) {
  prefix <- switch(level,
                   "ERROR"   = "[ERROR] ",
                   "WARNING" = "[WARNING] ",
                   "INFO"    = "[INFO] ",
                   "")
  if (!is.null(section)) message <- paste0(section, ": ", message)
  cat(prefix, message, "\n", file = log_file, append = TRUE)
  invisible(TRUE)
}

# Name‐normalization
normalize <- function(x) tolower(trimws(x))

# ─────────────────────────────────────────────────────────────────────────────
# Import Bacteria_Parameters.csv
# ─────────────────────────────────────────────────────────────────────────────

cat("Importing Bacteria_Parameters.csv...\n")
bacteria_parameters <- read.csv(file.path(wd, "input", "Bacteria_Parameters.csv"),
                                header = FALSE,
                                col.names = c("Starvation","Duplication","Death"),
                                stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Import metadata for each organism
# ─────────────────────────────────────────────────────────────────────────────

all_metabolites <- list()
all_reactions   <- list()
boundary_mets   <- list()

for (m in bacterial_models) {
  abbr   <- m$abbreviation
  org    <- m$organism
  dir    <- file.path(wd, "input", org)
  cat("Loading metadata for", org, "...\n")
  all_metabolites[[abbr]] <- read.csv(file.path(dir,"metabolites_metadata.csv"),
                                      stringsAsFactors=FALSE)
  all_reactions[[abbr]]   <- read.csv(file.path(dir,"reactions_metadata.csv"),
                                      stringsAsFactors=FALSE)
  boundary_mets[[abbr]]   <- filter(all_metabolites[[abbr]], grepl("_e$", id))
}

# ─────────────────────────────────────────────────────────────────────────────
# Initialize validation report
# ─────────────────────────────────────────────────────────────────────────────

cat("Writing validation report to:", log_file, "\n")
cat("PNPRO Validation Report\n=======================\n",
    file = log_file)
cat(sprintf("File: %s\nDate: %s\n\n",
            file_path, format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    file = log_file, append = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Fetch & parse XML
# ─────────────────────────────────────────────────────────────────────────────

xml_content <- tryCatch(read_xml(file_path), error = function(e){
  log_issue("ERROR", paste("Failed to parse XML:", e$message), "XML Parsing")
  stop("Aborting due to XML error")
})

# Extract place & transition names
places            <- xml_find_all(xml_content, "//place")
transitions       <- xml_find_all(xml_content, "//transition")
place_names       <- xml_attr(places,      "name")
transition_names  <- xml_attr(transitions, "name")
transition_delays <- xml_attr(transitions, "delay")
place_names_lc    <- normalize(place_names)
transition_names_lc <- normalize(transition_names)

# ─────────────────────────────────────────────────────────────────────────────
# 0. FBA Model & Transition Validation
# ─────────────────────────────────────────────────────────────────────────────

# Step counters
step_num <- 0

# Step 0: FBA Model & Transition Validation
log_section(sprintf("%d. FBA Model & Transition Validation", step_num))

# Pre–check all declared model files
for (model in bacterial_models) {
  fba_path <- file.path(model_dir, paste0(model$txt_file, ".txt"))
  if (!file.exists(fba_path)) {
    log_issue("ERROR", sprintf("Missing model file: %s", fba_path), "FBA Models")
  } else {
    log_issue("INFO", sprintf("Found model file: %s", fba_path), "FBA Models")
  }
}

# Then validate every FBA transition (syntax + semantic tests)
fba_idx <- grepl("^FBA\\[", transition_delays)
for  (i in which(fba_idx) ) {
  trans_name <- transition_names[i]
  delay_attr <- transition_delays[i]
  section    <- "FBA Validation"
  
  # Validate full syntax and extract parameters
  pattern <- 'FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9.]+) *, *"([^"]+)" *, *"([^"]+)"(?: *, *"(true|false)")? *\\]'
  match   <- regexec(pattern, delay_attr)
  parts   <- regmatches(delay_attr, match)[[1]]
  
  if ( (length(parts) - 1) < 6 ) {
    log_issue("ERROR", sprintf("Invalid FBA syntax in '%s': %s", trans_name, delay_attr), section)
    next
  }
  
  bacterial_file        <- parts[2]
  reaction_to_associate <- parts[3]
  scaling_constant      <- as.numeric(parts[4])
  bacterial_count_place <- parts[5]
  average_biomass_place <- parts[6]
  biomass_flag          <- if (length(parts) >= 7) parts[7] else "false"
  
  # === Validation Rules ===
  
  # 1. Model file must be known
  model_found <- any(sapply(bacterial_models, function(m) bacterial_file == paste0(m$FBAmodel, ".txt")))
  if (!model_found) {
    log_issue("ERROR", sprintf("Unknown model file '%s' in %s", bacterial_file, trans_name), section)
  }
  
  # 2. Reaction naming convention
  if (!grepl("^(EX|DM|SINK)_", reaction_to_associate, ignore.case = TRUE)) {
    log_issue("WARNING", sprintf("Unusual reaction name '%s' in %s", reaction_to_associate, trans_name), section)
  }
  
  # 3. Scaling constant should be numeric and > 0
  if ( is.na(scaling_constant) || scaling_constant <= 0 ) {
    log_issue("ERROR", sprintf("Non-positive scaling constant '%s' in %s", parts[4], trans_name), section)
  }
  
  # 4) Validate standard place naming for count and biomass
  # Determine organism abbreviation from context
  abbr <- NULL
  for (m in bacterial_models) {
    # Match on model file base or abbreviation in transition name
    if (grepl(m$abbreviation, trans_name, ignore.case = TRUE) ||
        grepl(m$FBAmodel, bacterial_file, ignore.case = TRUE)) {
      abbr <- normalize(m$abbreviation)
      break
    }
  }
  if (!is.null(abbr)) {
    expected_count   <- sprintf("n_%s", abbr)
    expected_biomass <- sprintf("biomass_e_%s", abbr)
    if (normalize(bacterial_count_place) != expected_count) {
      log_issue("ERROR", sprintf("Transition '%s': bacterialCountPlace '%s' should be '%s'", 
                                 trans_name, bacterial_count_place, expected_count),
                section)
    }
    if (normalize(average_biomass_place) != expected_biomass) {
      log_issue("ERROR", sprintf("Transition '%s': averageBiomassPlace '%s' should be '%s'", 
                                 trans_name, average_biomass_place, expected_biomass),
                section)
    }
  } else {
    log_issue("WARNING", sprintf("Could not infer abbreviation for transition '%s' for place validation", 
                                 trans_name), section)
  }
  
  # 5. Biomass flag logic
  if (grepl("biomass", reaction_to_associate, ignore.case = TRUE) && biomass_flag != "true") {
    log_issue("WARNING", sprintf("Missing biomassFlag='true' for biomass reaction in %s", trans_name), section)
  } else if (biomass_flag == "true" && !grepl("biomass", reaction_to_associate, ignore.case = TRUE)) {
    log_issue("WARNING", sprintf("biomassFlag=true on non-biomass reaction %s", trans_name), section)
  }
}

# If no FBA transitions found:
if (!any(fba_idx)) {
  log_issue("WARNING", "No FBA command transitions found in the net", "FBA Models")
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Place Validation
# ─────────────────────────────────────────────────────────────────────────────

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

# ─────────────────────────────────────────────────────────────────────────────
# 2. Required Transitions Validation
# ─────────────────────────────────────────────────────────────────────────────

# Step 2: Required Transitions Validation
step_num <- step_num + 1
log_section(sprintf("%d. Required Transitions Validation", step_num))
for (model in bacterial_models) {
  abbr <- normalize(model$abbreviation)
  log_issue("INFO", sprintf("Validating transitions for %s", model$organism), "Transitions")
  required <- c(
    sprintf("Starv_%s", abbr),
    sprintf("Death_%s", abbr),
    sprintf("Dup_%s", abbr),
    sprintf("EX_biomass_e_out_%s", abbr),
    sprintf("EX_biomass_e_in_%s", abbr)
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

# ─────────────────────────────────────────────────────────────────────────────
# 3. Shared Metabolite Validation
# ─────────────────────────────────────────────────────────────────────────────

# Step 3: Shared Metabolite Validation
step_num <- step_num + 1
log_section(sprintf("%d. Shared Metabolite Validation", step_num))

metabolite_places_lc <- metabolite_places

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

# ─────────────────────────────────────────────────────────────────────────────
# 4. Call Function Validation
# ─────────────────────────────────────────────────────────────────────────────

# Step 4: Call Function Validation
step_num <- step_num + 1
log_section(sprintf("%d. Call Function Validation", step_num))

# Identify Call[...] transitions
call_idx <- grepl("^Call\\[", transition_delays)
call_list <- which(call_idx)

if (length(call_list) > 0) {
  log_issue("INFO", sprintf("Found %d Call function transitions", length(call_list)), "Call Functions")
  
  # Pre-load Bacteria_Parameters.csv if needed
  uses_table <- any(grepl("FromTable\\[", transition_delays))
  if (uses_table) {
    params_path <- file.path(dirname(dirname(file_path)), "input", "Bacteria_Parameters.csv")
    if (!file.exists(params_path)) {
      log_issue("ERROR", sprintf("Missing parameters file: %s", params_path), "Call Functions")
      params_tbl <- NULL
    } else {
      params_tbl <- tryCatch(
        read.csv(params_path, header = FALSE, stringsAsFactors = FALSE),
        error = function(e) {
          log_issue("ERROR", paste("Error reading parameters:", e$message), "Call Functions")
          NULL
        }
      )
      if (!is.null(params_tbl)) {
        log_issue("INFO", sprintf("Parameters loaded: %d rows, %d cols", nrow(params_tbl), ncol(params_tbl)), "Call Functions")
        if (nrow(params_tbl) < length(bacterial_models)) {
          log_issue("ERROR", sprintf("Parameter rows (%d) < models (%d)", nrow(params_tbl), length(bacterial_models)), "Call Functions")
        }
      }
    }
  }
  
  # Get registered C++ functions
  registered <- tryCatch(c("Starvation", "Duplication", "Death"), error = function(e) {
    log_issue("WARNING", paste("Could not fetch registered functions:", e$message), "Call Functions")
    character()
  })
  
  # Validate each Call transition
  for (i in call_list) {
    trans_name <- transition_names[i]
    delay_attr <- transition_delays[i]
    section    <- "Call Functions"
    
    # Generic regex: function name + args string
    call_re <- '^Call\\[\\s*"([^"\\\\]+)"\\s*,\\s*(.*)\\s*\\]$'
    parts  <- regmatches(delay_attr, regexec(call_re, delay_attr, perl = TRUE))[[1]]
    if (length(parts) < 3) {
      log_issue("ERROR", sprintf("Invalid Call syntax: %s", delay_attr), section)
      next
    }
    fn_name <- parts[2]
    args_str <- parts[3]
    
    # Check function registration
    if (!(fn_name %in% registered)) {
      log_issue("WARNING", sprintf("Unregistered function '%s' in %s", fn_name, trans_name), section)
    }
    
    # Split args
    args <- trimws(strsplit(args_str, ",")[[1]])
    
    # Handle FromTable calls
    if (grepl("FromTable\\[", args_str)) {
      tbl_re <- 'FromTable\\[\\"Bacteria_Parameters.csv\\"\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]'
      tbl_parts <- regmatches(args_str, regexec(tbl_re, args_str, perl = TRUE))[[1]]
      if (length(tbl_parts) < 3) {
        log_issue("ERROR", sprintf("Malformed FromTable args: %s", args_str), section)
      } else if (!is.null(params_tbl)) {
        r <- as.integer(tbl_parts[2]) + 1
        c <- as.integer(tbl_parts[3]) + 1
        if (r < 1 || r > nrow(params_tbl)) {
          log_issue("ERROR", sprintf("Row index %d out of bounds", r-1), section)
        }
        if (c < 1 || c > ncol(params_tbl)) {
          log_issue("ERROR", sprintf("Col index %d out of bounds", c-1), section)
        }
      }
    }
  }
} else {
  log_issue("WARNING", "No Call function transitions found", "Call Functions")
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. Optimized Arc Extraction & FBA/Call Connectivity Logging
# ─────────────────────────────────────────────────────────────────────────────

step_num <- step_num + 1
log_section(sprintf("%d. Special Transitions & Arc Connectivity", step_num))

# Identify special transitions
is_call <- str_detect(transition_delays, "^Call\\[")
is_fba  <- str_detect(transition_delays, "^FBA\\[")
special <- transition_names[is_call | is_fba]

# Build arc‐dataframe in one go
arc_df <- xml_find_all(xml_content, "//arc") %>%
  map_df(~{
    a <- xml_attrs(.x)
    tibble(
      head         = a["head"],
      tail         = a["tail"],
      kind         = a["kind"],
      multiplicity = as.numeric(
        if_else(is.na(a["mult"]), "1", a["mult"])
      )
    )
  })%>%
  filter((kind=="INPUT"  & head %in% special) |
           (kind=="OUTPUT" & tail %in% special)) %>%
  transmute(
    transition   = if_else(kind=="INPUT", head, tail),
    direction    = kind,
    place        = if_else(kind=="INPUT", tail, head),
    multiplicity,
    command      = transition_delays[match(transition, transition_names)]
  )

log_section("Full arc_df contents")
# Convert to data.frame so print() won’t truncate
capture.output(
  print(as.data.frame(arc_df)),
  file   = log_file,
  append = TRUE
)

# Full per‐transition connectivity log (as before)
arc_df %>%
  arrange(transition, direction, place) %>%
  group_by(transition) %>%
  summarize(
    cmd   = first(command),
    conns = paste0(direction, " ", place, " (", multiplicity, ")", collapse = ", "),
    .groups = "drop"
  )

# Export CSV for downstream analysis
write_csv(arc_df,
          file.path(dirname(log_file),
                    paste0(tools::file_path_sans_ext(basename(file_path)), "_filtered_arcs.csv"))
)


# ─────────────────────────────────────────────────────────────────────────────
# 6. Arc Connectivity Pattern Validation
# ─────────────────────────────────────────────────────────────────────────────

# Import (IN) metabolite transition validation
validate_import_arcs <- function(arcs_df, trans_name, met, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    # This might be OK - not all metabolites need import transitions
    log_issue("INFO", sprintf("No import transition found for %s (%s)", met, trans_name), "Arc Connectivity")
    return(TRUE)
  }
  
  # Expected connections
  count_place <- paste0("n_", abbr)
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check input from metabolite
  has_input_from_met <- any(trans_arcs$direction == "INPUT" & 
                              tolower(trans_arcs$place) == met)
  
  # Check input from count place
  has_input_from_count <- any(trans_arcs$direction == "INPUT" & 
                                tolower(trans_arcs$place) == count_place)
  
  # Check output to count place
  has_output_to_count <- any(trans_arcs$direction == "OUTPUT" & 
                               tolower(trans_arcs$place) == count_place)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Check output to biomass
  has_output_to_biomass <- any(trans_arcs$direction == "OUTPUT" & 
                                 tolower(trans_arcs$place) == biomass_place)
  
  # Log any missing connections
  if (!has_input_from_met) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, met), 
              "Arc Connectivity")
  }
  
  if (!has_input_from_count) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, count_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_count) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, count_place), 
              "Arc Connectivity")
  }
  
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_biomass) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  return(has_input_from_met && has_input_from_count && has_output_to_count && 
           has_input_from_biomass && has_output_to_biomass)
}

# Export (OUT) metabolite transition validation
validate_export_arcs <- function(arcs_df, trans_name, met, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    # This might be OK - not all metabolites need export transitions
    log_issue("INFO", sprintf("No export transition found for %s (%s)", met, trans_name), "Arc Connectivity")
    return(TRUE)
  }
  
  # Expected connections
  count_place <- paste0("n_", abbr)
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check output to metabolite
  has_output_to_met <- any(trans_arcs$direction == "OUTPUT" & 
                             tolower(trans_arcs$place) == met)
  
  # Check input from count place
  has_input_from_count <- any(trans_arcs$direction == "INPUT" & 
                                tolower(trans_arcs$place) == count_place)
  
  # Check output to count place
  has_output_to_count <- any(trans_arcs$direction == "OUTPUT" & 
                               tolower(trans_arcs$place) == count_place)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Check output to biomass
  has_output_to_biomass <- any(trans_arcs$direction == "OUTPUT" & 
                                 tolower(trans_arcs$place) == biomass_place)
  
  # Log any missing connections
  if (!has_output_to_met) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, met), 
              "Arc Connectivity")
  }
  
  if (!has_input_from_count) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, count_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_count) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, count_place), 
              "Arc Connectivity")
  }
  
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_biomass) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  return(has_output_to_met && has_input_from_count && has_output_to_count && 
           has_input_from_biomass && has_output_to_biomass)
}

# Duplication transition validation
validate_dup_arcs <- function(arcs_df, trans_name, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    log_issue("ERROR", sprintf("No arcs found for %s", trans_name), "Arc Connectivity")
    return(FALSE)
  }
  
  # Expected connections
  count_place <- paste0("n_", abbr)
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check input from count
  has_input_from_count <- any(trans_arcs$direction == "INPUT" & 
                                tolower(trans_arcs$place) == count_place)
  
  # Check output to count with multiplicity 2
  has_output_to_count <- any(trans_arcs$direction == "OUTPUT" & 
                               tolower(trans_arcs$place) == count_place & 
                               trans_arcs$multiplicity == 2)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Check output to biomass
  has_output_to_biomass <- any(trans_arcs$direction == "OUTPUT" & 
                                 tolower(trans_arcs$place) == biomass_place)
  
  # Log any missing connections
  if (!has_input_from_count) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, count_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_count) {
    log_issue("ERROR", sprintf("%s missing output arc to %s with multiplicity 2", 
                               trans_name, count_place), "Arc Connectivity")
  }
  
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_biomass) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  return(has_input_from_count && has_output_to_count && 
           has_input_from_biomass && has_output_to_biomass)
}

# Death transition validation
validate_death_arcs <- function(arcs_df, trans_name, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    log_issue("ERROR", sprintf("No arcs found for %s", trans_name), "Arc Connectivity")
    return(FALSE)
  }
  
  # Expected connections
  count_place <- paste0("n_", abbr)
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check input from count
  has_input_from_count <- any(trans_arcs$direction == "INPUT" & 
                                tolower(trans_arcs$place) == count_place)
  
  # Check no output to count (cell is removed)
  has_no_output_to_count <- !any(trans_arcs$direction == "OUTPUT" & 
                                   tolower(trans_arcs$place) == count_place)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Check output to biomass
  has_output_to_biomass <- any(trans_arcs$direction == "OUTPUT" & 
                                 tolower(trans_arcs$place) == biomass_place)
  
  # Log any missing connections
  if (!has_input_from_count) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, count_place), 
              "Arc Connectivity")
  }
  
  if (!has_no_output_to_count) {
    log_issue("WARNING", sprintf("%s should not have output arc to %s", 
                                 trans_name, count_place), "Arc Connectivity")
  }
  
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_biomass) {
    log_issue("ERROR", sprintf("%s missing output arc to %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  return(has_input_from_count && has_no_output_to_count && 
           has_input_from_biomass && has_output_to_biomass)
}

# Starvation transition validation
validate_starv_arcs <- function(arcs_df, trans_name, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    log_issue("ERROR", sprintf("No arcs found for %s", trans_name), "Arc Connectivity")
    return(FALSE)
  }
  
  # Expected connections
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Check no output to biomass (biomass is consumed)
  has_no_output_to_biomass <- !any(trans_arcs$direction == "OUTPUT" & 
                                     tolower(trans_arcs$place) == biomass_place)
  
  # Log any missing connections
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_no_output_to_biomass) {
    log_issue("WARNING", sprintf("%s should not have output arc to %s", 
                                 trans_name, biomass_place), "Arc Connectivity")
  }
  
  return(has_input_from_biomass && has_no_output_to_biomass)
}

# Biomass exchange - IN transition validation
validate_biomass_in_arcs <- function(arcs_df, trans_name, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    log_issue("ERROR", sprintf("No arcs found for %s", trans_name), "Arc Connectivity")
    return(FALSE)
  }
  
  # Expected connections
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Biomass IN should not have output back to biomass (or multiplicity of 0)
  has_no_output <- !any(trans_arcs$direction == "OUTPUT" & 
                          tolower(trans_arcs$place) == biomass_place & 
                          trans_arcs$multiplicity > 0)
  
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_no_output) {
    log_issue("WARNING", sprintf("%s should not have output arc to %s", 
                                 trans_name, biomass_place), "Arc Connectivity")
  }
  
  return(has_input_from_biomass && has_no_output)
}

# Biomass exchange - OUT transition validation
validate_biomass_out_arcs <- function(arcs_df, trans_name, abbr) {
  # Filter arcs for this transition
  trans_arcs <- filter(arcs_df, tolower(transition) == trans_name)
  
  if (nrow(trans_arcs) == 0) {
    log_issue("ERROR", sprintf("No arcs found for %s", trans_name), "Arc Connectivity")
    return(FALSE)
  }
  
  # Expected connections
  biomass_place <- paste0("biomass_e_", abbr)
  
  # Check input from biomass
  has_input_from_biomass <- any(trans_arcs$direction == "INPUT" & 
                                  tolower(trans_arcs$place) == biomass_place)
  
  # Check output to biomass with multiplicity 2
  has_output_to_biomass <- any(trans_arcs$direction == "OUTPUT" & 
                                 tolower(trans_arcs$place) == biomass_place & 
                                 trans_arcs$multiplicity == 2)
  
  if (!has_input_from_biomass) {
    log_issue("ERROR", sprintf("%s missing input arc from %s", trans_name, biomass_place), 
              "Arc Connectivity")
  }
  
  if (!has_output_to_biomass) {
    log_issue("ERROR", sprintf("%s missing output arc to %s with multiplicity 2", 
                               trans_name, biomass_place), "Arc Connectivity")
  }
  
  return(has_input_from_biomass && has_output_to_biomass)
}

step_num <- step_num + 1
log_section(sprintf("%d. Arc Connectivity Pattern Validation", step_num))

# Process each transition type and validate connectivity patterns
for (model in bacterial_models) {
  abbr <- normalize(model$abbreviation)
  log_issue("INFO", sprintf("Validating arc patterns for %s", model$organism), "Arc Connectivity")
  
  # 1. Biomass Exchange Transitions
  log_subsection("Biomass Exchange Transitions")
  # For EX_biomass_e_out
  biomass_out_name <- paste0("EX_biomass_e_out_", abbr)
  validate_biomass_out_arcs(arc_df, biomass_out_name, abbr)
  
  # For EX_biomass_e_in
  biomass_in_name <- paste0("EX_biomass_e_in_", abbr)
  validate_biomass_in_arcs(arc_df, biomass_in_name, abbr)
  
  # 2. Metabolite Boundary Transitions
  log_subsection("Metabolite Boundary Transitions")
  for (met in metabolite_places_lc) {
    # For import transitions
    import_name <- paste0("EX_", met, "_in_", abbr)
    validate_import_arcs(arc_df, import_name, met, abbr)
    
    # For export transitions
    export_name <- paste0("EX_", met, "_out_", abbr)
    validate_export_arcs(arc_df, export_name, met, abbr)
  }
  
  # 3. Population Transitions
  log_subsection("Population Transitions")
  # Duplication transition
  dup_name <- paste0("Dup_", abbr)
  validate_dup_arcs(arc_df, dup_name, abbr)
  
  # Death transition
  death_name <- paste0("Death_", abbr)
  validate_death_arcs(arc_df, death_name, abbr)
  
  # Starvation transition
  starv_name <- paste0("Starv_", abbr)
  validate_starv_arcs(arc_df, starv_name, abbr)
}

# ─────────────────────────────────────────────────────────────────────────────
# 8. Model Completeness Summary
# ─────────────────────────────────────────────────────────────────────────────

step_num <- step_num + 1
log_section(sprintf("%d. Model Completeness Summary", step_num))

# Read the log file to count errors and warnings
log_content <- readLines(log_file)
error_count <- sum(grepl("\\[ERROR\\]", log_content))
warning_count <- sum(grepl("\\[WARNING\\]", log_content))

# Calculate completeness score
completeness <- 100 * (1 - error_count/(error_count + 50))  # 50 is a scaling factor
completeness <- max(0, min(100, completeness))  # Bound between 0-100%

log_issue("INFO", sprintf("Model has %d errors and %d warnings", 
                          error_count, warning_count), "Completeness")
log_issue("INFO", sprintf("Estimated model completeness: %.1f%%", 
                          completeness), "Completeness")

# Overall model status
if (error_count == 0) {
  if (warning_count == 0) {
    log_issue("INFO", "Model status: FULLY VALIDATED ✓", "Completeness")
  } else {
    log_issue("INFO", "Model status: VALID WITH WARNINGS ⚠", "Completeness")
  }
} else if (completeness >= 75) {
  log_issue("INFO", "Model status: NEARLY COMPLETE (address errors) ⚠", "Completeness")
} else if (completeness >= 40) {
  log_issue("INFO", "Model status: PARTIALLY COMPLETE (significant errors) ⚠", "Completeness")
} else {
  log_issue("INFO", "Model status: MAJOR ISSUES (consider rebuilding) ✗", "Completeness")
}

# ─────────────────────────────────────────────────────────────────────────────
# 9. Fix Suggestions
# ─────────────────────────────────────────────────────────────────────────────
step_num <- step_num + 1
log_section(sprintf("%d. Fix Suggestions", step_num))

if (error_count > 0 || warning_count > 0) {
  log_issue("INFO", "The following fixes are suggested:", "Fixes")
  
  # Place naming suggestions
  for (model in bacterial_models) {
    abbr <- normalize(model$abbreviation)
    expected_places <- c(sprintf("n_%s", abbr), sprintf("biomass_e_%s", abbr))
    
    for (place in expected_places) {
      if (!(place %in% place_names_lc)) {
        # Look for similar names that might need renaming
        similar <- grep(substr(place, 1, 3), place_names_lc, value = TRUE)
        if (length(similar) > 0) {
          log_issue("INFO", sprintf("Rename place '%s' to '%s'", similar[1], place), "Fixes")
        } else {
          log_issue("INFO", sprintf("Add missing place '%s'", place), "Fixes")
        }
      }
    }
  }
  
  # FBA command fixes
  for (i in which(fba_idx)) {
    trans_name <- transition_names[i]
    delay_attr <- transition_delays[i]
    
    # Find model file issues
    pattern <- 'FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9.]+) *, *"([^"]+)" *, *"([^"]+)"(?: *, *"(true|false)")? *\\]'
    match <- regexec(pattern, delay_attr)
    parts <- regmatches(delay_attr, match)[[1]]
    
    if (length(parts) > 1) {
      bacterial_file <- parts[2]
      reaction_to_associate <- parts[3]
      bacterial_count_place <- parts[5]
      average_biomass_place <- parts[6]
      biomass_flag <- if (length(parts) >= 7) parts[7] else "false"
      
      # Check model file path
      model_found <- FALSE
      correct_model_file <- ""
      
      for (model in bacterial_models) {
        if (bacterial_file == paste0(model$FBAmodel, ".txt")) {
          model_found <- TRUE
          break
        } else if (grepl(model$abbreviation, trans_name, ignore.case = TRUE)) {
          correct_model_file <- paste0(model$txt_file, ".txt")
        }
      }
      
      if (!model_found && correct_model_file != "") {
        # Suggest correct model file
        correct_command <- gsub(
          pattern = sprintf('"%s"', bacterial_file),
          replacement = sprintf('"%s"', correct_model_file),
          x = delay_attr
        )
        log_issue("INFO", sprintf("In transition '%s', replace command with:\n%s", 
                                  trans_name, correct_command), "Fixes")
      }
      
      # Check if biomass flag is needed but missing
      if (grepl("biomass", reaction_to_associate, ignore.case = TRUE) && 
          (length(parts) < 7 || biomass_flag != "true")) {
        correct_command <- gsub(
          pattern = '\\]$',
          replacement = ', "true"]',
          x = delay_attr
        )
        log_issue("INFO", sprintf("In transition '%s', add biomass flag:\n%s", 
                                  trans_name, correct_command), "Fixes")
      }
      
      # Check place naming in command
      abbr <- NULL
      for (model in bacterial_models) {
        if (grepl(model$abbreviation, trans_name, ignore.case = TRUE)) {
          abbr <- normalize(model$abbreviation)
          break
        }
      }
      
      if (!is.null(abbr)) {
        expected_count <- sprintf("n_%s", abbr)
        expected_biomass <- sprintf("biomass_e_%s", abbr)
        
        if (normalize(bacterial_count_place) != expected_count) {
          correct_command <- gsub(
            pattern = sprintf('"%s"', bacterial_count_place),
            replacement = sprintf('"%s"', expected_count),
            x = delay_attr
          )
          log_issue("INFO", sprintf("In transition '%s', fix bacterial count place:\n%s", 
                                    trans_name, correct_command), "Fixes")
        }
        
        if (normalize(average_biomass_place) != expected_biomass) {
          correct_command <- gsub(
            pattern = sprintf('"%s"', average_biomass_place),
            replacement = sprintf('"%s"', expected_biomass),
            x = delay_attr
          )
          log_issue("INFO", sprintf("In transition '%s', fix biomass place:\n%s", 
                                    trans_name, correct_command), "Fixes")
        }
      }
    }
  }
  
  # Call command syntax fixes
  for (i in which(call_idx)) {
    trans_name <- transition_names[i]
    delay_attr <- transition_delays[i]
    
    # Generic regex: function name + args string
    call_re <- '^Call\\[\\s*"([^"\\\\]+)"\\s*,\\s*(.*)\\s*\\]$'
    parts <- regmatches(delay_attr, regexec(call_re, delay_attr, perl = TRUE))[[1]]
    if (length(parts) >= 3) {
      fn_name <- parts[2]
      args_str <- parts[3]
      
      # Check if the FromTable syntax is correct
      if (grepl("FromTable\\[", args_str)) {
        tbl_re <- 'FromTable\\[\\"([^"]+)\\"\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]'
        tbl_parts <- regmatches(args_str, regexec(tbl_re, args_str, perl = TRUE))[[1]]
        
        if (length(tbl_parts) < 4) {
          log_issue("INFO", sprintf("In transition '%s', fix FromTable syntax:\nCall[\"%s\", FromTable[\"Bacteria_Parameters.csv\", 0, 1], 0]", 
                                    trans_name, fn_name), "Fixes")
        } else {
          # Get organism index from transition name
          org_idx <- NA
          for (j in seq_along(bacterial_models)) {
            if (grepl(bacterial_models[[j]]$abbreviation, trans_name, ignore.case = TRUE)) {
              org_idx <- j - 1  # 0-based indexing
              break
            }
          }
          
          if (!is.na(org_idx)) {
            # Check if row index is correct
            row_idx <- as.integer(tbl_parts[3])
            if (row_idx != org_idx) {
              correct_args <- gsub(
                pattern = sprintf(', %d,', row_idx),
                replacement = sprintf(', %d,', org_idx),
                x = args_str
              )
              log_issue("INFO", sprintf("In transition '%s', fix organism index:\nCall[\"%s\", %s]", 
                                        trans_name, fn_name, correct_args), "Fixes")
            }
          }
        }
      }
      
      # Check if the last parameter (organism index) is present
      if (!grepl(',\\s*[0-9]+\\s*$', args_str)) {
        # Get organism index from transition name
        org_idx <- NA
        for (j in seq_along(bacterial_models)) {
          if (grepl(bacterial_models[[j]]$abbreviation, trans_name, ignore.case = TRUE)) {
            org_idx <- j - 1  # 0-based indexing
            break
          }
        }
        
        if (!is.na(org_idx)) {
          correct_command <- gsub(
            pattern = '\\]$',
            replacement = sprintf(', %d]', org_idx),
            x = delay_attr
          )
          log_issue("INFO", sprintf("In transition '%s', add organism index:\n%s", 
                                    trans_name, correct_command), "Fixes")
        }
      }
    }
  }
  
  # Arc connectivity and multiplicity fixes
  # Process the arc dataframe to find connectivity issues
  arc_issues <- data.frame(
    transition = character(),
    issue_type = character(),
    place = character(),
    fix = character(),
    stringsAsFactors = FALSE
  )
  
  # Check for biomass out transitions (should have multiplicity 2)
  for (model in bacterial_models) {
    abbr <- normalize(model$abbreviation)
    biomass_place <- paste0("biomass_e_", abbr)
    biomass_out_trans <- paste0("EX_biomass_e_out_", abbr)
    
    # Find relevant arcs
    biomass_out_arcs <- arc_df %>%
      filter(tolower(transition) == biomass_out_trans,
             tolower(place) == biomass_place,
             direction == "OUTPUT")
    
    if (nrow(biomass_out_arcs) > 0) {
      if (any(biomass_out_arcs$multiplicity != 2)) {
        arc_issues <- rbind(arc_issues, data.frame(
          transition = biomass_out_trans,
          issue_type = "multiplicity",
          place = biomass_place,
          fix = "Set output multiplicity to 2",
          stringsAsFactors = FALSE
        ))
      }
    } else {
      arc_issues <- rbind(arc_issues, data.frame(
        transition = biomass_out_trans,
        issue_type = "missing_arc",
        place = biomass_place,
        fix = "Add output arc to biomass_e with multiplicity 2",
        stringsAsFactors = FALSE
      ))
    }
    
    # Check for duplication transitions (should have multiplicity 2 to cell count)
    dup_trans <- paste0("dup_", abbr)
    count_place <- paste0("n_", abbr)
    
    dup_arcs <- arc_df %>%
      filter(tolower(transition) == dup_trans,
             tolower(place) == count_place,
             direction == "OUTPUT")
    
    if (nrow(dup_arcs) > 0) {
      if (any(dup_arcs$multiplicity != 2)) {
        arc_issues <- rbind(arc_issues, data.frame(
          transition = dup_trans,
          issue_type = "multiplicity",
          place = count_place,
          fix = "Set output multiplicity to 2",
          stringsAsFactors = FALSE
        ))
      }
    } else {
      arc_issues <- rbind(arc_issues, data.frame(
        transition = dup_trans,
        issue_type = "missing_arc",
        place = count_place,
        fix = "Add output arc to cell count with multiplicity 2",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Log all arc issues
  if (nrow(arc_issues) > 0) {
    log_issue("INFO", "Arc connectivity and multiplicity issues:", "Fixes")
    for (i in 1:nrow(arc_issues)) {
      log_issue("INFO", sprintf("Transition '%s' - %s - Place '%s': %s", 
                                arc_issues$transition[i],
                                arc_issues$issue_type[i],
                                arc_issues$place[i],
                                arc_issues$fix[i]), "Fixes")
    }
  }
  
  # Check for metabolite-related arcs
  for (model in bacterial_models) {
    
    
    related_rxns <- character()
    
    # 1. Check reaction equations if available in metadata
    if ("equation" %in% names(all_reactions[[model$abbreviation]])) {
      # Get boundary reactions where the metabolite appears in the equation
      eq_matches <- character()
      for (j in which(all_reactions[[model$abbreviation]]$type == "boundary")) {
        # Check if metabolite is in the equation as a standalone term
        # This uses word boundaries to avoid partial matches
        if (grepl(paste0("\\b", met, "\\b"), all_reactions[[model$abbreviation]]$equation[j])) {
          eq_matches <- c(eq_matches, all_reactions[[model$abbreviation]]$abbreviation[j])
        }
      }
      related_rxns <- c(related_rxns, eq_matches)
    }
    
    abbr <- normalize(model$abbreviation)
    
    for (met in metabolite_places_lc) {
      # Check import transitions
      import_trans <- paste0(related_rxns, "_in_", abbr)
      
      # Should have input from metabolite
      met_input_arcs <- arc_df %>%
        filter(tolower(transition) == import_trans,
               tolower(place) == met,
               direction == "INPUT")
      
      if (nrow(met_input_arcs) == 0) {
        arc_issues <- rbind(arc_issues, data.frame(
          transition = import_trans,
          issue_type = "missing_arc",
          place = met,
          fix = "Add input arc from metabolite",
          stringsAsFactors = FALSE
        ))
      }
      
      # Check export transitions
      export_trans <-  paste0(related_rxns, "_out_", abbr)
      
      # Should have output to metabolite
      met_output_arcs <- arc_df %>%
        filter(tolower(transition) == export_trans,
               tolower(place) == met,
               direction == "OUTPUT")
      
      if (nrow(met_output_arcs) == 0) {
        arc_issues <- rbind(arc_issues, data.frame(
          transition = export_trans,
          issue_type = "missing_arc",
          place = met,
          fix = "Add output arc to metabolite",
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Additional arc issue fixes
  if (nrow(arc_issues) > 0) {
    log_issue("INFO", "Additional arc fixes:", "Fixes")
    unique_transitions <- unique(arc_issues$transition)
    
    for (trans in unique_transitions) {
      trans_issues <- arc_issues[arc_issues$transition == trans, ]
      log_issue("INFO", sprintf("For transition '%s', make these arc changes:", trans), "Fixes")
      
      for (i in 1:nrow(trans_issues)) {
        log_issue("INFO", sprintf("  • %s with place '%s'", 
                                  trans_issues$fix[i], 
                                  trans_issues$place[i]), "Fixes")
      }
    }
  }
  
} else {
  log_issue("INFO", "No fixes needed - model is fully compliant", "Fixes")
}