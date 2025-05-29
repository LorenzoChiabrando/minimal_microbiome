# generate_R_from_pnpro.R

#' Generate a hypernode‐specific R script from:
#'   • a .PNPRO file (to get place ordering)
#'   • the functions_hypernode_template.R (with two stubs)
#'   • your biounit_models list (to get cell & biomass initials)
#'
#' @param pnpro_file       Path to the .PNPRO XML file
#' @param r_template       Path to functions_hypernode_template.R
#' @param biounit_models List of model objects (each with $initial_count and $biomass$mean)
#' @param output_r         Path for the generated R script (defaults to overwrite r_template)
#' @return Invisibly, the character vector of lines written
generate_R_from_pnpro <- function(pnpro_file,
                                  r_template,
                                  biounit_models,
                                  output_r = r_template) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Please install the 'xml2' package (install.packages('xml2'))")
  }
  
  # 1) Parse the PNPRO and extract place names in order
  doc         <- xml2::read_xml(pnpro_file)
  places_node <- xml2::xml_find_all(doc, ".//nodes/place")
  place_names <- xml2::xml_attr(places_node, "name")
  
  # 2) Build the initial marking vector:
  #    - for "n_<id>": biounit_models[[…]]$initial_count
  #    - for "biomass_e_<id>": biounit_models[[…]]$biomass$mean
  #    - otherwise: 1
  init_vals <- vapply(place_names, function(pl) {
    if (grepl("^n_", pl)) {
      id  <- sub("^n_", "", pl)
      idx <- which(sapply(biounit_models, function(m) id %in% sub("^n_", "", m$abbreviation)))
      return(biounit_models[[idx]]$initial_count)
    }
    if (grepl("^biomass_e_", pl)) {
      id  <- sub("^biomass_e_", "", pl)
      idx <- which(sapply(biounit_models, function(m) id %in% sub("^n_", "", m$abbreviation)))
      return(biounit_models[[idx]]$biomass$mean)
    }
    1
  }, numeric(1))
  
  # 3) Read the template
  lines <- readLines(r_template)
  
  # 4) Replace the stub for yini.names
  stub_names <- sprintf(
    "  yini.names <- c(%s)",
    paste(shQuote(place_names), collapse = ", ")
  )
  lines <- sub("^\\s*yini\\.names\\s*<-.*$", stub_names, lines)
  
  # 5) Replace the stub for y_ini
  stub_vals <- sprintf(
    "  y_ini <- c(%s)",
    paste(init_vals, collapse = ", ")
  )
  lines <- sub("^\\s*y_ini\\s*<-.*$", stub_vals, lines)
  
  # 6) Write out the generated file
  writeLines(lines, output_r)
  invisible(lines)
}
