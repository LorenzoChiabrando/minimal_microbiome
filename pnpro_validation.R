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

# Step 0: FBA Model & Transition Validation
step_num <- step_num + 1
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

# Step 3: Shared Metabolite Validation
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

# Identify “registered” Call[...] and FBA[...] transitions by their delay attribute
call_transitions <- transition_names[grepl("^Call\\[", transition_delays)]
fba_transitions  <- transition_names[grepl("^FBA\\[",  transition_delays)]
interesting_trans <- c(call_transitions, fba_transitions)

arcs = edges
# Build a data.frame of all arcs
arc_df <- data.frame(
  head = xml2::xml_attr(arcs, "head"),
  tail = xml2::xml_attr(arcs, "tail"),
  kind = xml2::xml_attr(arcs, "kind"),
  stringsAsFactors = FALSE
)

# Filter for arcs involving our transitions and reshape
filtered_arcs <- arc_df %>%
  dplyr::filter(
    (kind == "INPUT"  & head %in% interesting_trans) |
      (kind == "OUTPUT" & tail %in% interesting_trans)
  ) %>%
  dplyr::mutate(
    transition = ifelse(kind == "INPUT", head, tail),
    direction  = kind,
    place      = ifelse(kind == "INPUT", tail, head)
  ) %>%
  dplyr::select(transition, direction, place)