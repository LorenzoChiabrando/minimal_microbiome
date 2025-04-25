reorganize_pnpro <- function(input_file,
                             arc_df_repaired,
                             output_file) {
  # load xml2 + stringr
  if (!requireNamespace("xml2", quietly=TRUE) ||
      !requireNamespace("stringr", quietly=TRUE)) {
    stop("Please install.packages(c('xml2','stringr'))")
  }
  
  # 1) parse the PNPRO
  doc   <- xml2::read_xml(input_file)
  nodes <- xml2::xml_find_first(doc, "//gspn/nodes")
  
  # 2) derive places & organisms from your repaired arc‐df
  counts      <- unique(arc_df_repaired$place[grepl("^n_",          arc_df_repaired$place)])
  biomasses   <- unique(arc_df_repaired$place[grepl("^biomass_e_",  arc_df_repaired$place)])
  metabolites <- setdiff(unique(arc_df_repaired$place), c(counts, biomasses))
  abbrs       <- sub("^n_","", counts)
  
  N <- length(abbrs)
  M <- length(metabolites)
  
  # 3) compute canvas size dynamically
  metabolite_spacing <- 15
  module_spacing     <- 20
  margin             <- 5
  
  width  <- max((M+1)*metabolite_spacing, (N+1)*module_spacing) + margin*2
  height <- margin +            # bottom margin
    module_spacing +    # space for modules
    (module_spacing*0.8) + # for Call[] offsets
    module_spacing +    # gap between layers
    module_spacing +    # for metabolites
    margin              # top margin
  
  # 4) compute Y positions
  top_y    <- height - margin - module_spacing   # metabolites
  bottom_y <- margin + module_spacing            # modules
  
  # 5) X positions
  xs_met <- seq(margin + metabolite_spacing,
                width  - margin - metabolite_spacing,
                length.out = M)
  xs_mod <- seq(margin + module_spacing,
                width  - margin - module_spacing,
                length.out = N)
  
  module_radius <- module_spacing * 0.6
  call_offsets  <- c(
    Starv =  module_radius*0.8,
    Dup   =  0,
    Death = -module_radius*0.8
  )
  
  # 6) reposition <place> nodes
  place_nodes <- xml2::xml_find_all(nodes, "place")
  place_names <- xml2::xml_attr(place_nodes, "name")
  
  # 6a) metabolites
  for (j in seq_along(metabolites)) {
    p  <- metabolites[j]
    nd <- place_nodes[place_names == p]
    xml2::xml_set_attr(nd, "x",        as.character(xs_met[j]))
    xml2::xml_set_attr(nd, "y",        as.character(top_y))
    xml2::xml_set_attr(nd, "label-x",  as.character(xs_met[j]))
    xml2::xml_set_attr(nd, "label-y",  as.character(top_y + 20))
  }
  
  # 6b) modules: count & biomass
  for (i in seq_along(abbrs)) {
    cx <- xs_mod[i]; cy <- bottom_y
    # count place
    cp <- place_nodes[place_names == counts[i]]
    xml2::xml_set_attr(cp, "x", as.character(cx - module_radius/2))
    xml2::xml_set_attr(cp, "y", as.character(cy))
    # biomass place
    bp <- place_nodes[place_names == biomasses[i]]
    xml2::xml_set_attr(bp, "x", as.character(cx + module_radius/2))
    xml2::xml_set_attr(bp, "y", as.character(cy))
  }
  
  # 7) reposition <transition> nodes
  trans_nodes <- xml2::xml_find_all(nodes, "transition")
  trans_names <- xml2::xml_attr(trans_nodes, "name")
  
  for (k in seq_along(trans_names)) {
    tn    <- trans_names[k]
    node  <- trans_nodes[k]
    delay <- xml2::xml_attr(node, "delay")
    
    # A) Call[…] transitions
    for (fn in names(call_offsets)) {
      if (startsWith(tn, paste0(fn, "_"))) {
        ab  <- sub(paste0(fn, "_"), "", tn)
        idx <- match(ab, abbrs)
        xml2::xml_set_attr(node, "x", as.character(xs_mod[idx]))
        xml2::xml_set_attr(node, "y", as.character(bottom_y + call_offsets[[fn]]))
      }
    }
    
    # B) FBA[…] transitions
    if (startsWith(delay, "FBA[")) {
      parts       <- stringr::str_match(
        delay,
        'FBA\\s*\\["[^"]+","([^"]+)",[^,]+,"([^"]+)"'
      )
      reaction    <- parts[2]
      count_place <- parts[3]
      ab          <- sub("^n_","", count_place)
      idx_o       <- match(ab, abbrs)
      
      if (reaction == "EX_biomass_e") {
        bp_nd <- place_nodes[place_names == paste0("biomass_e_",ab)]
        x0    <- as.numeric(xml2::xml_attr(bp_nd,"x"))
        y0    <- as.numeric(xml2::xml_attr(bp_nd,"y"))
        vof   <- if (grepl("_in_", tn)) +module_radius*0.8 else -module_radius*0.8
        xml2::xml_set_attr(node, "x", as.character(x0))
        xml2::xml_set_attr(node, "y", as.character(y0 + vof))
      } else {
        mets <- unique(
          arc_df_repaired$place[
            arc_df_repaired$transition == tn &
              arc_df_repaired$place %in% metabolites
          ]
        )
        if (length(mets) != 1)
          stop("Cannot find unique metabolite for ", tn)
        met <- mets[[1]]
        xm  <- xs_met[ match(met,       metabolites) ]
        xo  <- xs_mod[idx_o]
        mx  <- (xm + xo)/2
        my  <- (top_y + bottom_y)/2
        vo  <- if (grepl("_in_", tn)) +2 else -2
        xml2::xml_set_attr(node, "x", as.character(mx))
        xml2::xml_set_attr(node, "y", as.character(my + vo))
      }
    }
  }
  
  # 8) write out the new PNPRO
  dir.create(dirname(output_file), recursive=TRUE, showWarnings=FALSE)
  xml2::write_xml(doc, output_file, options="format")
  message("Laid out PNPRO: ", output_file,
          sprintf("  [%d organisms × %d metabolites, canvas %dx%d]",
                  N, M, width, height))
}
