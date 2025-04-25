repair_pnpro <- function(arc_df_repaired,
                         project_name,
                         gspn_name,
                         output_file) {
  # 1) collect unique places & transitions
  places        <- unique(arc_df_repaired$place)
  trans_cmd_tbl <- unique(arc_df_repaired[, c("transition","command")])
  
  # 2) start XML
  doc  <- xml2::xml_new_root("project", name = project_name, version = "121")
  gspn <- xml2::xml_add_child(doc, "gspn", name = gspn_name)
  
  # 3) places + transitions under <nodes>
  nodes <- xml2::xml_add_child(gspn, "nodes")
  for (p in places) {
    xml2::xml_add_child(nodes, "place",
                        name     = p,
                        x        = "0.0", y = "0.0",
                        `label-x` = "0.0", `label-y` = "0.0"
    )
  }
  for (i in seq_len(nrow(trans_cmd_tbl))) {
    row <- trans_cmd_tbl[i, ]
    xml2::xml_add_child(nodes, "transition",
                        name     = row$transition,
                        type     = "EXP",
                        rotation = "0.0",
                        x        = "0.0", y = "0.0",
                        delay    = row$command
    )
  }
  
  # 4) arcs under <edges> — now with type="normal"
  edges <- xml2::xml_add_child(gspn, "edges")
  for (i in seq_len(nrow(arc_df_repaired))) {
    r <- arc_df_repaired[i, ]
    if (r$direction == "INPUT") {
      head <- r$transition; tail <- r$place
    } else {
      head <- r$place;      tail <- r$transition
    }
    xml2::xml_add_child(edges, "arc",
                        head = head, tail = tail,
                        kind = r$direction,
                        type = "normal",
                        mult = as.character(r$multiplicity)
    )
  }
  
  # 5) mandatory <measures> block
  m <- xml2::xml_add_child(doc, "measures",
                           `gspn-name`     = gspn_name,
                           name            = "Measures",
                           `simplified-UI` = "false"
  )
  xml2::xml_add_child(m, "assignments")
  xml2::xml_add_child(m, "greatspn")
  f <- xml2::xml_add_child(m, "formulas")
  xml2::xml_add_child(f, "formula",
                      comment  = "Basic statistics of the toolchain execution.",
                      language = "STAT"
  )
  xml2::xml_add_child(f, "formula",
                      comment  = "All the basic Petri net measures",
                      language = "ALL"
  )
  
  # 6) write file
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  xml2::write_xml(doc, output_file, options = "format")
  message("Repaired PNPRO written to: ", output_file)
}
