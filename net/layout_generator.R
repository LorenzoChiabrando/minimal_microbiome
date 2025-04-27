# layout_generator.R
# Generates a hierarchical, hypernode-based layout for PNPRO files
# - Exposes layout_pnpro(input, output)
# - Adds community, organism, FBA wrapper, module, and boundary boxes
# - Arranges places and transitions inside each module in a grid layout

library(dplyr)
library(xml2)

#' Layout a PNPRO with hypernode/community boxes and auto-arrange nodes
#'
#' @param input_pnpro Path to original .PNPRO file
#' @param output_pnpro Path where the laid-out .PNPRO will be written
layout_pnpro <- function(input_pnpro, output_pnpro) {
  # Parse existing document
  doc <- read_xml(input_pnpro)
  nds <- xml_find_first(doc, "//nodes")

  # Identify organisms and boundary counts
  fba_in <- xml_find_all(doc, "//transition[starts-with(@delay, 'FBA[') and contains(@name,'_in_')]")
  df <- data.frame(
    abbr = sub('^.*_in_', '', xml_attr(fba_in,'name')),
    stringsAsFactors = FALSE
  )
  bnds <- df %>% count(abbr, name='n_boundary')

  # Dimensions
  module_w <- 9; module_h <- 15
  react_w <- 4; react_h_unit <- 3; pad <- 1
  bnds <- bnds %>% mutate(
    wrapper_w = react_w * n_boundary + pad*2,
    wrapper_h = max(module_h, react_h_unit * n_boundary) + pad*2
  )
  n_org <- nrow(bnds)

  # Community box
  comm_w <- sum(bnds$wrapper_w) + (n_org+1)*pad
  comm_h <- max(bnds$wrapper_h) + pad*2
  xml_add_child(nds, 'text-box',
    name='community_box', shape='ROUND_RECTANGLE',
    `fill-color`='#f3f6ff', `border-color`='#000000', `text-color`='#000000',
    width=as.character(comm_w), height=as.character(comm_h), x='0', y='0',
    `label-x`='0', `label-y`='0'
  )

  # Boundary metabolites pool on left
  all_mets <- xml_find_all(doc, "//transition[starts-with(@delay,'FBA[')]") %>%
    xml_attr('name') %>% grep('_in_', ., value=TRUE) %>%
    sub('_in_.*$','',.) %>% unique()
  met_h <- react_h_unit * length(all_mets) + pad*2
  xml_add_child(nds,'text-box',
    name='boundary_metabolites_box', shape='ROUND_RECTANGLE',
    `fill-color`='#fff3fc', `border-color`='#000000', `text-color`='#000000',
    width=as.character(react_w+pad), height=as.character(met_h),
    x=as.character(pad), y=as.character(comm_h - met_h - pad),
    `label-x`='0', `label-y`='0'
  )

  # Place wrappers and nested module/boundary boxes
  x_ptr <- pad + react_w + pad
  layouts <- list()
  for(i in seq_len(n_org)){
    org <- bnds$abbr[i]; w <- bnds$wrapper_w[i]; h <- bnds$wrapper_h[i]
    y_base <- comm_h - h - pad
    # wrapper
    xml_add_child(nds,'text-box', name=paste0(org,'_fba_box'), shape='ROUND_RECTANGLE',
      `fill-color`='#ffffff', `border-color`='#000000', `text-color`='#000000',
      width=as.character(w), height=as.character(h), x=as.character(x_ptr), y=as.character(y_base),
      `label-x`='0', `label-y`='0'
    )
    # core module
    xml_add_child(nds,'text-box', name=paste0(org,'_box'), shape='ROUND_RECTANGLE',
      `fill-color`='#fff3f3', `border-color`='#000000', `text-color`='#000000',
      width=as.character(module_w), height=as.character(module_h),
      x=as.character(x_ptr + pad), y=as.character(y_base + h - module_h - pad),
      `label-x`='0', `label-y`='0'
    )
    # boundary reactions
    br_h <- h - module_h - pad*2
    xml_add_child(nds,'text-box', name=paste0('boundary_reactions_',org,'_box'), shape='ROUND_RECTANGLE',
      `fill-color`='#f4fff3', `border-color`='#000000', `text-color`='#000000',
      width=as.character(react_w), height=as.character(br_h),
      x=as.character(x_ptr + w - react_w - pad),
      y=as.character(y_base + pad),`label-x`='0', `label-y`='0'
    )
    layouts[[org]] <- list(x=x_ptr+pad, y=y_base+pad, w=w-pad*2, h=h-pad*2)
    x_ptr <- x_ptr + w + pad
  }

  # Arrange places & transitions inside each organism wrapper in grid
  places <- xml_find_all(doc,'//place')
  trans  <- xml_find_all(doc,'//transition')
  for(org in names(layouts)){
    L <- layouts[[org]]
    pl <- places[grepl(paste0('_',org,'($|_)'), xml_attr(places,'name'))]
    tr <- trans[grepl(paste0('_',org,'($|_)'), xml_attr(trans,'name'))]
    # grid dims for places
    if(length(pl)>0){
      nc <- ceiling(sqrt(length(pl)))
      xgap <- (module_w - pad*2) / (nc+1)
      ygap <- (module_h - pad*2) / (nc+1)
      for(j in seq_along(pl)){
        cx <- L$x + ((j-1) %% nc +1) * xgap
        cy <- L$y + L$h - ((floor((j-1)/nc)+1) * ygap)
        xml_set_attr(pl[[j]],'x',as.character(cx)); xml_set_attr(pl[[j]],'y',as.character(cy))
      }
    }
    # grid dims for transitions
    if(length(tr)>0){
      nt <- ceiling(sqrt(length(tr)))
      xgap2 <- (module_w - pad*2)/(nt+1)
      ygap2 <- (module_h - pad*2)/(nt+1)
      for(j in seq_along(tr)){
        cx <- L$x + ((j-1) %% nt +1) * xgap2
        cy <- L$y + ((floor((j-1)/nt)+1) * ygap2)
        xml_set_attr(tr[[j]],'x',as.character(cx)); xml_set_attr(tr[[j]],'y',as.character(cy))
      }
    }
  }

  # Write out the updated PNPRO
  write_xml(doc, output_pnpro)
}

# CLI
if(sys.nframe()==0){
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)!=2)stop('Usage: Rscript layout_generator.R in.PNPRO out.PNPRO')
  layout_pnpro(args[1],args[2])
}

