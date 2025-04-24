#' Validate a PNPRO file for FBA-integrated bacterial models
#'
#' Applies the standard rules for FBA-driven Petri nets:
#'  • FBA[...] syntax & model-file existence
#'  • Place naming (n_<abbr>, biomass_e_<abbr>)
#'  • Required transitions (Starv_, Death_, Dup_, EX_biomass_e_in_, EX_biomass_e_out_)
#'  • Shared metabolite places
#'  • Arc connectivity for every FBA or Call transition
#'
#' @param pnpro_path     Path to your .PNPRO file
#' @param bacterial_models  A list of lists, each with
#'        $abbr     — character abbreviation (e.g. "ec")
#'        $FBAmodel — model base name (no `.txt`)
#' @param metabolite_places  Character vector of expected shared metabolite place IDs
#' @param model_dir      Directory where you keep compiled_models/*.txt  
#' @param log_dir        Directory to write validation logs (CSV outputs)
#' @return Invisibly returns a list with:
#'         - issues: tibble of validation issues
#'         - arc_df: tibble of arcs for FBA/Call transitions
validate_pnpro <- function(pnpro_path,
                           bacterial_models,
                           metabolite_places = character(),
                           model_dir,
                           log_dir) {
  library(xml2)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(readr)
  
  issues <- list()
  log_issue <- function(level, section, message, object = NA_character_) {
    issues <<- append(issues, list(tibble(
      level   = level,
      section = section,
      message = message,
      object  = object
    )))
  }
  
  # Parse XML
  xml <- tryCatch(read_xml(pnpro_path),
                  error = function(e) stop("Failed to parse XML: ", e$message))
  places      <- xml_find_all(xml, "//place")
  transitions <- xml_find_all(xml, "//transition")
  place_names       <- xml_attr(places,      "name") %>% tolower()
  transition_names  <- xml_attr(transitions, "name")
  transition_delays <- xml_attr(transitions, "delay")
  
  # Build arc_df for FBA/Call transitions
  is_call <- str_detect(transition_delays, "^Call\\[")
  is_fba  <- str_detect(transition_delays, "^FBA\\[")
  special <- transition_names[is_call | is_fba]
  arc_df <- xml_find_all(xml, "//arc") %>%
    map_df(~{
      a <- xml_attrs(.x)
      tibble(
        head         = a["head"],
        tail         = a["tail"],
        kind         = a["kind"],
        multiplicity = as.integer(dplyr::coalesce(a["mult"], "1"))
      )
    }) %>%
    filter((kind=="INPUT"  & head %in% special) |
             (kind=="OUTPUT" & tail %in% special)) %>%
    transmute(
      transition   = if_else(kind=="INPUT", head, tail),
      direction    = kind,
      place        = if_else(kind=="INPUT", tail, head),
      multiplicity,
      command      = transition_delays[match(transition, transition_names)]
    )
  
  # 0. FBA models exist on-disk
  for (m in bacterial_models) {
    fpath <- file.path(model_dir, paste0(m$FBAmodel, ".txt"))
    if (!file.exists(fpath)) {
      log_issue("ERROR", "FBA Models",
                sprintf("Missing model file: %s", fpath))
    }
  }
  
  # 0b. Validate each FBA[...] transition
  fba_idx <- which(str_detect(transition_delays, "^FBA\\["))
  pat <- 'FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9.]+) *, *"([^"]+)" *, *"([^"]+)"(?: *, *"(true|false)")? *\\]'
  for (i in fba_idx) {
    nm    <- transition_names[i]
    delay <- transition_delays[i]
    parts <- regmatches(delay, regexec(pat, delay))[[1]]
    if (length(parts) < 6) {
      log_issue("ERROR", "FBA Validation",
                sprintf("Bad syntax in %s: %s", nm, delay), nm)
      next
    }
    # model-file check
    model_file <- parts[2]
    if (!any(sapply(bacterial_models,
                    function(m) paste0(m$FBAmodel, ".txt")==model_file))) {
      log_issue("ERROR", "FBA Validation",
                sprintf("Unknown model file '%s' in %s", model_file, nm), nm)
    }
    # scaling-constant > 0
    sc_str <- parts[3]
    if (!grepl("^[0-9]+(?:\\.[0-9]+)?$", sc_str)) {
      log_issue("ERROR", "FBA Validation",
                sprintf("Invalid scaling constant '%s' in %s", sc_str, nm), nm)
    } else {
      sc <- as.numeric(sc_str)
      if (sc <= 0) {
        log_issue("ERROR", "FBA Validation",
                  sprintf("Non-positive scaling constant '%s' in %s", sc_str, nm), nm)
      }
    }
    # place-names
    count_p   <- parts[4]
    biomass_p <- parts[5]
    abbr <- NULL
    for (m in bacterial_models) {
      if (grepl(m$abbr[2], nm, ignore.case=TRUE) ||
          grepl(m$FBAmodel, model_file, ignore.case=TRUE)) {
        abbr <- tolower(m$abbr[2]); break
      }
    }
    if (!is.null(abbr)) {
      exp_count <- sprintf("n_%s", abbr)
      exp_biom  <- sprintf("biomass_e_%s", abbr)
      if (tolower(count_p)   != exp_count)
        log_issue("ERROR","FBA Validation",
                  sprintf("Count place '%s' ≠ '%s'", count_p, exp_count), nm)
      if (tolower(biomass_p) != exp_biom)
        log_issue("ERROR","FBA Validation",
                  sprintf("Biomass place '%s' ≠ '%s'", biomass_p, exp_biom), nm)
    } else {
      log_issue("WARNING","FBA Validation",
                sprintf("Could not infer abbreviation for %s", nm), nm)
    }
  }
  if (length(fba_idx) == 0) {
    log_issue("WARNING", "FBA Models", "No FBA[...] transitions found")
  }
  
  # 1. Place naming
  for (m in bacterial_models) {
    abbr <- tolower(m$abbr[2])
    cplace <- sprintf("n_%s", abbr)
    bplace <- sprintf("biomass_e_%s", abbr)
    if (!cplace %in% place_names)
      log_issue("ERROR","Places",sprintf("Missing place '%s'", cplace))
    if (!bplace %in% place_names)
      log_issue("ERROR","Places",sprintf("Missing place '%s'", bplace))
  }
  
  # 2. Required transitions
  for (m in bacterial_models) {
    abbr <- tolower(m$abbr[2])
    req  <- c(sprintf("Starv_%s",abbr),
              sprintf("Death_%s",abbr),
              sprintf("Dup_%s",abbr),
              sprintf("EX_biomass_e_in_%s",abbr),
              sprintf("EX_biomass_e_out_%s",abbr))
    for (r in req) {
      if (!tolower(r) %in% tolower(transition_names))
        log_issue("ERROR","Required Transitions",
                  sprintf("Missing %s", r))
    }
  }
  
  # 3. Shared metabolites
  for (met in metabolite_places) {
    if (!tolower(met) %in% place_names)
      log_issue("ERROR","Shared Metabolites",
                sprintf("Missing metabolite place: %s", met))
  }
  
  # 4. Arc connectivity
  for (nm in special) {
    sub <- arc_df %>% filter(transition==nm)
    if (nrow(sub)==0) {
      log_issue("ERROR","Arc Connectivity",
                sprintf("No arcs attached to %s", nm), nm)
    } else {
      if (!any(sub$direction=="INPUT"))
        log_issue("ERROR","Arc Connectivity",
                  sprintf("%s has no INPUT arc", nm), nm)
      if (!any(sub$direction=="OUTPUT"))
        log_issue("ERROR","Arc Connectivity",
                  sprintf("%s has no OUTPUT arc", nm), nm)
    }
  }
  
  # compile results
  issues_df <- bind_rows(issues)
  
  # write CSV logs
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive=TRUE)
  base <- tools::file_path_sans_ext(basename(pnpro_path))
  write_csv(issues_df, file.path(log_dir, paste0(base, "_issues.csv")))
  write_csv(arc_df,    file.path(log_dir, paste0(base, "_arc_df.csv")))
  
  invisible(list(issues=issues_df, arc_df=arc_df))
}
