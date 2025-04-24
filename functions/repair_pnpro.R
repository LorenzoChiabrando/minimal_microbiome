#' Repair a PNPRO Petri net by regenerating from scratch
#'
#' Uses:
#'  - `<base>_issues.csv`, `<base>_arc_df.csv` in `log_dir`
#'  - original `.PNPRO` metadata (project & measures blocks)
#'  - `bacterial_models` list with `$abbr` and `$FBAmodel`
#'  - `metabolite_places` vector of shared metabolites
#'
#' @param pnpro_path        Path to the original `.PNPRO` file
#' @param bacterial_models  List of lists, each containing:
#'                           $abbr     — two‐letter abbreviation string
#'                           $FBAmodel — base name of model file (no `.txt`)
#' @param metabolite_places Character vector of boundary metabolite place IDs
#' @param log_dir           Directory containing `<base>_issues.csv` and `<base>_arc_df.csv`
#' @param out_path          Path to write the repaired `.PNPRO` (defaults to `_repaired.PNPRO`)
#' @return Invisibly TRUE on success
repair_pnpro <- function(pnpro_path,
                         bacterial_models,
                         metabolite_places = character(),
                         log_dir = dirname(pnpro_path),
                         out_path = NULL) {
  library(xml2); library(dplyr); library(readr); library(stringr)
  
  # derive base filename
  base <- tools::file_path_sans_ext(basename(pnpro_path))
  
  # read original arc & issues
  issues <- read_csv(file.path(log_dir, paste0(base, '_issues.csv')), show_col_types = FALSE)
  arc_df  <- read_csv(file.path(log_dir, paste0(base, '_arc_df.csv')), show_col_types = FALSE)
  
  # generate arc_df_repaired by adding back missing arcs
  missing <- issues %>%
    filter(section == 'Arc Connectivity' & str_detect(message, 'missing')) %>%
    transmute(
      transition   = object,
      direction    = if_else(str_detect(message, 'INPUT'), 'INPUT', 'OUTPUT'),
      place        = str_extract(message, "(?<=arc (to|from) )[A-Za-z0-9_]+"),
      multiplicity = 1L,
      command      = NA_character_
    )
  
  arc_df_repaired <- bind_rows(arc_df, missing) %>%
    distinct(transition, direction, place, multiplicity, command)
  
  # save for audit
  write_csv(arc_df_repaired,
            file.path(log_dir, paste0(base, '_arc_df_repaired.csv')))
  
  # parse original for metadata blocks
  doc_orig  <- read_xml(pnpro_path)
  proj_node <- xml_find_first(doc_orig, '/project')
  proj_attrs<- xml_attrs(proj_node)
  gspn_node <- xml_find_first(proj_node, './gspn')
  gspn_attrs<- xml_attrs(gspn_node)
  measures  <- xml_find_first(proj_node, './measures')
  
  # default output path
  if (is.null(out_path)) out_path <- file.path(dirname(pnpro_path), paste0(base,'_repaired.PNPRO'))
  
  # build new XML doc
  new_doc <- xml_new_root('project')
  xml_set_attrs(new_doc, proj_attrs)
  new_gspn <- xml_add_child(new_doc, 'gspn')
  xml_set_attrs(new_gspn, gspn_attrs)
  
  # 1) nodes
  nodes_parent <- xml_add_child(new_gspn, 'nodes')
  # shared metabolites
  for (met in metabolite_places) {
    plc <- xml_add_child(nodes_parent, 'place')
    xml_set_attrs(plc, c(name=met, x='0.0', y='0.0', 'label-x'='0.0', 'label-y'='0.0'))
  }
  # bacterial count/biomass places
  for (m in bacterial_models) {
    abbr <- tolower(m$abbr)
    for (pfx in c(paste0('n_',abbr), paste0('biomass_e_',abbr))) {
      plc <- xml_add_child(nodes_parent, 'place')
      xml_set_attrs(plc, c(name=pfx, x='0.0', y='0.0', 'label-x'='0.0', 'label-y'='0.0'))
    }
  }
  
  # 2) transitions
  # dynamic population transitions
  dyn_tr <- unlist(lapply(bacterial_models, function(m) {
    ab <- tolower(m$abbr)
    c(paste0('Starv_',ab), paste0('Death_',ab), paste0('Dup_',ab),
      paste0('EX_biomass_e_in_',ab), paste0('EX_biomass_e_out_',ab))
  }))
  
  for (tr in unique(dyn_tr)) {
    tnode <- xml_add_child(nodes_parent, 'transition')
    xml_set_attrs(tnode, c(name=tr, type='EXP', rotation='0.0', x='0.0', y='0.0', delay='0'))
  }
  # FBA/Call-driven transitions
  for (tr in unique(arc_df_repaired$transition)) {
    if (!tr %in% dyn_tr) {
      cmd <- arc_df_repaired$command[arc_df_repaired$transition==tr][1]
      tnode <- xml_add_child(nodes_parent, 'transition')
      xml_set_attrs(tnode, c(name=tr, type='EXP', rotation='0.0', x='0.0', y='0.0', delay=cmd))
    }
  }
  
  # 3) measures block
  xml_add_child(new_doc, measures)
  
  # 4) edges
  edges_parent <- xml_add_child(new_gspn, 'edges')
  apply(arc_df_repaired, 1, function(r) {
    r <- as.list(r)
    head <- if (r$direction=='INPUT') r$transition else r$place
    tail <- if (r$direction=='INPUT') r$place      else r$transition
    arc <- xml_add_child(edges_parent, 'arc')
    xml_set_attrs(arc, c(head=head, tail=tail, kind=r$direction, mult=as.character(r$multiplicity)))
    invisible(NULL)
  })
  
  write_xml(new_doc, out_path, options=c('format','no_declaration'))
  invisible(TRUE)
}
