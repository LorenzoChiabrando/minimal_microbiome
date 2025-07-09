# generate_cpp_from_arcs.R

#' Generate a hypernode‐specific C++ file by filling in:
#'   • V                (culture volume in mL)
#'   • delta            (cell density in cell/mL)
#'   • bacteria_names   (vector of "n_<id>")
#'   • bacteriaBiomass_names (vector of "biomass_e_<id>")
#'
#' @param arcs_csv     Path to repaired_arcs.csv
#' @param cpp_template Path to general_functions_template.cpp
#' @param output_cpp   Path to write the instantiated C++ (defaults to overwrite cpp_template)
#' @param volume       Numeric: culture volume in mL
#' @param cell_density Numeric: cell density in cell/mL
#' @return (invisibly) the character vector of lines written
generate_cpp_from_arcs <- function(arcs_csv,
                                   cpp_template,
                                   output_cpp    = cpp_template,
                                   volume,
                                   cell_density) {
  # 1) read arcs and discover organism codes
  arcs  <- read.csv(arcs_csv, stringsAsFactors = FALSE)
  places <- unique(arcs$place)
  orgs   <- unique(sub("^n_", "", places[grepl("^n_", places)]))
  
  # 2) load the C++ template
  lines <- readLines(cpp_template)
  
  # 3) inject volume and cell density
  lines <- sub(
    "^\\s*double\\s+V\\s*=.*?;",
    sprintf("double V = %s;", volume),
    lines
  )
  lines <- sub(
    "^\\s*long\\s+long\\s+int\\s+delta\\s*=.*?;",
    sprintf("long long int delta = %s;", cell_density),
    lines
  )
  
  # 4) helper to rewrite a vector<string> block, matching std::vector or vector
  rewrite_block <- function(lines, field, values) {
    # match either "vector<string>" or "std::vector<string>"
    start_pat <- sprintf("^\\s*(?:std::)?vector<[^>]*string>\\s+%s\\s*=\\s*\\{", field)
    start_i   <- grep(start_pat, lines)
    if (length(start_i) != 1) stop("cannot find start of ", field)
    end_i     <- grep("^\\s*\\};", lines)
    end_i     <- end_i[end_i > start_i][1]
    if (is.na(end_i))  stop("cannot find end of ", field)
    
    # build new block
    new_block <- c(
      lines[start_i],
      paste0("    ", paste(sprintf('"%s"', values), collapse = ", ")),
      "};"
    )
    
    # splice in
    c(
      lines[1:(start_i-1)],
      new_block,
      lines[(end_i+1):length(lines)]
    )
  }
  
  # 5) inject bacteria_names
  lines <- rewrite_block(
    lines,
    "bacteria_names",
    paste0("n_", orgs)
  )
  
  # 6) inject bacteriaBiomass_names
  lines <- rewrite_block(
    lines,
    "bacteriaBiomass_names",
    paste0("biomass_e_", orgs)
  )
  
  # 7) write out
  writeLines(lines, output_cpp)
  invisible(lines)
}
