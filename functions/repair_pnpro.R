repair_pnpro <- function(pnpro_path, bacterial_models, metabolite_places,
                          arc_df_repaired, rx_meta) {

  # write it out
  write_xml(new_doc, sub("\\.PNPRO$", "_final.PNPRO", pnpro_path),
            options=c("format","no_declaration"))
}