# ─────────────────────────────────────────────────────────────────────────────
# Auto-Repair PNPRO from Validation Log
# ─────────────────────────────────────────────────────────────────────────────

library(xml2)
library(stringr)
library(purrr)
library(jsonlite)

# ─────────────────────────────────────────────────────────────────────────────
# Inputs (must exist in your environment)
# ─────────────────────────────────────────────────────────────────────────────
# log_file         — path to your validation log
# arc_df           — data.frame of arcs extracted earlier
# transition_names — character vector of transition @name
# transition_delays— character vector of transition @delay
# bacterial_models — list(...) as defined up-front
# xml_content      — xml2 document read from your .PNPRO
# file_path        — original PNPRO file path
# model_name       — base name for output file
# ─────────────────────────────────────────────────────────────────────────────

xml_content <- tryCatch(read_xml(file_path), error = function(e){
  log_issue("ERROR", paste("Failed to parse XML:", e$message), "XML Parsing")
  stop("Aborting due to XML error")
})

log_lines <- readLines(log_file)

# Export CSV for downstream analysis
arc_df = read_csv(file.path(dirname(log_file), paste0(tools::file_path_sans_ext(basename(file_path)), "_filtered_arcs.csv")))

# ─────────────────────────────────────────────────────────────────────────────
# 1) PLACE RENAMES
#    Lines like: Rename place 'old' to 'new'
# ─────────────────────────────────────────────────────────────────────────────
ren_lns <- grep("Rename place '", log_lines, value = TRUE)
ren_mat <- str_match(ren_lns, "Rename place '([^']+)' to '([^']+)'")
fixes_rename <- map(seq_len(nrow(ren_mat)), ~list(
  type     = "rename_place",
  old_name = ren_mat[.x,2],
  new_name = ren_mat[.x,3]
))

# ─────────────────────────────────────────────────────────────────────────────
# 2) FBA COMMAND FIXES
#    Driven by: [ERROR] FBA Validation: Unknown model file 'BAD.txt' in TRANS
# ─────────────────────────────────────────────────────────────────────────────
fba_err <- str_match(log_lines, "Unknown model file '([^']+)' in ([^ ]+)")
fba_err <- fba_err[!is.na(fba_err[,1]), , drop = FALSE]
fixes_fba <- map(seq_len(nrow(fba_err)), function(i) {
  bad_file <- fba_err[i,2]
  trans    <- fba_err[i,3]
  # pick model by abbreviation in transition name
  idx      <- which(map_lgl(bacterial_models,
                            ~ str_detect(trans, fixed(.x$abbreviation, ignore_case = TRUE))))
  model    <- bacterial_models[[ idx[1] ]]
  correct  <- paste0(model$FBAmodel, ".txt")
  old_delay<- transition_delays[transition_names == trans]
  new_delay<- str_replace(old_delay, fixed(bad_file), correct)
  list(
    type        = "fix_fba_command",
    transition  = trans,
    new_command = new_delay
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# 3) MISSING ARC ADDITIONS
#    Driven by: Transition 'T' - missing_arc - Place 'P': Add (in|out) arc … multiplicity N
# ─────────────────────────────────────────────────────────────────────────────
ma_lns <- grep("missing_arc", log_lines, value = TRUE)
ma_mat <- str_match(
  ma_lns,
  "Transition '([^']+)' - missing_arc - Place '([^']+)': Add (input|output) arc .* multiplicity (\\d+)"
)
ma_mat <- ma_mat[!is.na(ma_mat[,1]), , drop = FALSE]
fixes_add_arc <- map(seq_len(nrow(ma_mat)), ~{
  trans <- ma_mat[.x,2]
  place <- ma_mat[.x,3]
  dir   <- toupper(ma_mat[.x,4])
  mult  <- as.integer(ma_mat[.x,5])
  list(
    type     = "fix_arc_mult",
    head     = if(dir=="INPUT")  place else trans,
    tail     = if(dir=="INPUT")  trans else place,
    kind     = dir,
    new_mult = mult
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# 4) EXISTING ARC MULTIPLICITY CORRECTIONS
#    Driven by lines like: Set output multiplicity to N
# ─────────────────────────────────────────────────────────────────────────────
mx_lns <- grep("multiplicity to \\d+", log_lines, value = TRUE)
mx_mat <- str_match(
  mx_lns,
  "Transition '([^']+)' - multiplicity - Place '([^']+)': Set (input|output) multiplicity to (\\d+)"
)
mx_mat <- mx_mat[!is.na(mx_mat[,1]), , drop = FALSE]
fixes_mult <- map(seq_len(nrow(mx_mat)), ~{
  trans <- mx_mat[.x,2]
  place <- mx_mat[.x,3]
  dir   <- toupper(mx_mat[.x,4])
  mult  <- as.integer(mx_mat[.x,5])
  list(
    type     = "fix_arc_mult",
    head     = if(dir=="INPUT")  place else trans,
    tail     = if(dir=="INPUT")  trans else place,
    kind     = dir,
    new_mult = mult
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# 5) COMBINE & DEDUPE FIXES
# ─────────────────────────────────────────────────────────────────────────────
all_fixes <- c(fixes_rename, fixes_fba, fixes_add_arc, fixes_mult)
# remove duplicates by JSON string
uniq_json <- unique(map_chr(all_fixes, ~ toJSON(.x, auto_unbox=TRUE)))
fixes <- map(uniq_json, ~ fromJSON(.x))

# ─────────────────────────────────────────────────────────────────────────────
# 6) REPAIR FUNCTION
# ─────────────────────────────────────────────────────────────────────────────
repair_xml <- function(xml_content, fixes, fixed_path) {
  # deep copy via tempfile
  tmp <- tempfile(fileext = ".xml"); on.exit(unlink(tmp), add=TRUE)
  write_xml(xml_content, file = tmp)
  modified <- read_xml(tmp)
  for (fix in fixes) {
    switch(fix$type,
           rename_place = {
             n <- xml_find_first(modified, sprintf("//place[@name='%s']", fix$old_name))
             if (!inherits(n, "xml_missing")) xml_set_attr(n, "name", fix$new_name)
           },
           fix_fba_command = {
             n <- xml_find_first(modified, sprintf("//transition[@name='%s']", fix$transition))
             if (!inherits(n, "xml_missing")) xml_set_attr(n, "delay", fix$new_command)
           },
           fix_call_command = {
             n <- xml_find_first(modified, sprintf("//transition[@name='%s']", fix$transition))
             if (!inherits(n, "xml_missing")) xml_set_attr(n, "delay", fix$new_command)
           },
           fix_arc_mult = {
             xp <- sprintf("//arc[@head='%s' and @tail='%s' and @kind='%s']",
                           fix$head, fix$tail, fix$kind)
             arc <- xml_find_first(modified, xp)
             if (!inherits(arc, "xml_missing")) xml_set_attr(arc, "mult", as.character(fix$new_mult))
           },
           {
             warning("Unknown fix type: ", fix$type)
           }
    )
  }
  write_xml(modified, file = fixed_path)
  invisible(fixed_path)
}

# ─────────────────────────────────────────────────────────────────────────────
# 7) APPLY ALL FIXES
# ─────────────────────────────────────────────────────────────────────────────
fixed_path <- file.path(dirname(file_path), paste0(model_name, "_auto.PNPRO"))
repair_xml(xml_content, fixes, fixed_path)
cat("Repaired PNPRO written to:", fixed_path, "\n")
