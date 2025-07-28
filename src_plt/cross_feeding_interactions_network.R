
# This study presents a constraint-based approach to model minimal microbiota 
# with specific functions. Using SIHUMIx as an example, we demonstrate how microbial 
# communities process metabolic resources and produce various compounds.
# 
# Our workflow quantifies community growth rates and metabolite production while preserving 
# individual species functionality. The model identifies keystone species and reveals 
# an emergent metabolic interaction network where fermentation products are exchanged 
# between populations, with secretion and uptake determined by overall metabolic flux.

## example 4 debug 
# wd = "~/Documents/FBA_project"
# cross_fed_met = c("ac_e","ppa_e","but_e", "glu_L_e", "lac_L_e"),
# t_frame = 1:48
# col_bac = "#f06b368c"
# col_met = "lightblue"
# shape_bac = "square"
# shape_met = "circle"

plotting_int_net = function(wd,
                            cross_fed_met,
                            t_frame,
                            col_bac,
                            col_met,
                            shape_bac,
                            shape_met
                            ) {
  
  library(dplyr)
  library(igraph)
  library(colourvalues)
  
  bac_model <- list.files(paste0(wd, "/code/SIHUMI_test"), pattern = "\\.mat$")
  
  files <- list.files(paste0(wd, "/results_SIHUMI"), pattern = "\\.flux$")
  # Initialize an empty list to store dataframes
  list_of_dataframes <- list()
  
  # Loop through each file
  for (f in files) {
    subflux <- utils::read.table(paste0(wd, "/results_SIHUMI/", f), header = TRUE) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      mutate(SourceFile = f)
    
    list_of_dataframes[[f]] <- subflux
  }
  
  combined_df <- bind_rows(list_of_dataframes)
  
  mets = c()
  react2plot = c()
  
  met_interaction_net = vector(mode = "list", length = length(bac_model))
  names(met_interaction_net) = stringr::str_remove(bac_model, ".mat")
  
  for(bac in bac_model) {
    
    all_react <- readRDS(
      paste0(wd, "/code/SIHUMI_test/additional_data_metabolic_models/", 
             stringr::str_remove(bac, ".mat"), "/all_react.rds"))
    
    reaction_data = all_react[[1]]
    
    equation <- readRDS(
      paste0(wd, "/code/SIHUMI_test/additional_data_metabolic_models/", 
             stringr::str_remove(bac, ".mat"), "/equation.rds"))
    
    reaction_data$Eq = equation[reaction_data$ReactionPos]
    
    for(met in cross_fed_met){
      pattern <- paste0("\\b", met, "\\b")
      r = reaction_data$React_ID[grepl(pattern, reaction_data$Eq)]
      react2plot = c(react2plot, r)
      mets = c(mets, rep(met, length(r)))
    }
  }
  
  df_react_met <- data.frame(mets = mets, react2plot = react2plot)
  cross_fed_reaction <- unique(df_react_met$react2plot)
  
  exchanges <- cross_fed_reaction[which(grepl("EX_", cross_fed_reaction))]
  
  fluxes_ex_all <- combined_df %>%
    filter(Reaction %in% exchanges) %>%
    mutate(SourceFile = stringr::str_remove(SourceFile, ".flux")) %>%
    mutate(SourceFile = stringr::str_remove(SourceFile, "result-")) %>%
    mutate(Reaction = stringr::str_remove(Reaction, "EX_")) %>%
    rename("Metabolite" = "Reaction") %>%
    mutate(Time = as.character(round(Time, 2)))
  
  for (time_point in t_frame) {
    
    fluxes_ex <- combined_df %>%
      filter(Reaction %in% exchanges) %>%
      mutate(SourceFile = stringr::str_remove(SourceFile, ".flux")) %>%
      mutate(SourceFile = stringr::str_remove(SourceFile, "result-")) %>%
      mutate(Reaction = stringr::str_remove(Reaction, "EX_")) %>%
      rename("Metabolite" = "Reaction") %>%
      mutate(Time = as.character(round(Time, 2))) %>%
      filter(Time == time_point) %>%
      mutate(Time = as.double(Time)) %>%
      select(-Time) %>%
      mutate(Type = ifelse(Flux > 0, "Producer", "Consumer")) %>%
      mutate(SourceFile = stringr::str_remove_all(SourceFile, ".flux")) %>%
      mutate(SourceFile = stringr::str_remove_all(SourceFile, "result-")) %>%
      mutate(SourceFile = stringr::str_extract(SourceFile, "^[^_]+_[^_]+")) %>%
      mutate(SourceFile = stringr::str_replace_all(SourceFile, "_", " ")) %>%
      rename("Microorganism" = "SourceFile") %>%
      filter(abs(Flux) >= 1e-06)
    
    # Correcting the edges and nodes
    edges <- fluxes_ex %>%
      dplyr::mutate(Type = ifelse(Flux > 0, "produces", "consumes")) %>%
      dplyr::select(Microorganism, Metabolite, Type) %>%
      dplyr::mutate(Microorganism_type = ifelse(
        Type == "produces",
        paste(Microorganism, "(+)", sep = " "),
        paste(Microorganism, "(-)", sep = " "))) %>% 
      dplyr::select(c(Metabolite, Microorganism_type))
    
    nodes <- data.frame(name = unique(c(as.character(edges$Microorganism_type), as.character(edges$Metabolite))),
                        type = c(rep("Microorganism", length(unique(edges$Microorganism_type))),
                                 rep("Metabolite", length(unique(edges$Metabolite)))))
    
    df_edges<-data.frame(x = NA, y = NA)
    
    for (i in 1:nrow(edges)) {
      if (grepl("\\(\\+\\)", edges$Microorganism_type[i])){df_edges[i,]<-c(edges$Microorganism_type[i],edges$Metabolite[i])}
      else{df_edges[i,]<-edges[i,]}
      }
    
    # Create the graph
    g <- graph_from_data_frame(d = df_edges, vertices = nodes, directed = T)
    
    # E(g)$direction <- ifelse(grepl("\\(\\+\\)", edges$Microorganism_type), "produces", "consumes")
    # 
    # for ( i in seq_along(E(g)) ) {
    #   if ( E(g)$direction[i] == "consumes" ) {
    #     # Set direction from Metabolite to Microorganism
    #     # No action needed since we start with directed edges; simply ensure the from-to is correct
    #   } else {
    #     # If "produces", we reverse the direction from Microorganism to Metabolite
    #     edge <- ends(g, E(g)[i])  # Get the nodes involved in the edge
    #     g <- delete_edges(g, E(g)[i])  # Remove the existing edge
    #     # Add the edge in reverse order
    #     g <- add_edges(g, c(edge[2], edge[1]))
    #   }
    # }
    
    # Set colors for each node type
    V(g)$color <- ifelse(V(g)$type == "Microorganism", col_bac, col_met)
    # Set the shape for each layer of nodes (use consistent shapes for simplicity)
    V(g)$shape <- ifelse(V(g)$type == "Microorganism", shape_bac, shape_met)
    # Set sizes for the nodes using igraph's degree function
    V(g)$size <- (igraph::degree(g, mode = "all") * 2) + 1
    # Set the font size for vertex labels to make it more readable
    V(g)$label.cex <- 1
    V(g)$label.font <- 1
    
    # Set sizes for the edges using the "fluxes_ex" dataframe
    # E(g)$width <- log2(abs(fluxes_ex$Flux) + 1) * 3
    # Set a uniform edge width
    E(g)$width <- 2
    
    # Normalize the Flux values for transparency (between 0 and 1)
    # max_flux <- max(abs(fluxes_ex_all$Flux)) * 0.05
    # 
    # Adjust edge color and width for according to flux sign
    # E(g)$color <- ifelse(grepl("\\(\\+\\)", edges$Microorganism_type), linkcol_producer, linkcol_consumer)
    #
    # transparency <- abs(fluxes_ex$Flux) / max_flux  # scale between 0 and 1
    # 
    # df_edges$flux<-fluxes_ex$Flux 
    # df_edges$flux_d<-fluxes_ex$Flux* 0.05
    # df_edges$tr<-transparency
    # 
    # # Set edge colors with transparency based on Flux values
    # E(g)$color <- mapply(function(col, alpha) adjustcolor(col, alpha.f = alpha), 
    #                      E(g)$color, transparency)
    
    flux_interval = data.frame(x = c(-unique(abs(fluxes_ex_all$Flux)), unique(abs(fluxes_ex_all$Flux))))
    
    # Apply the palette to your data
    pal <- colour_values(flux_interval$x,
                         palette = "RdBu",
                         na_colour = "#808080FF",
                         n_summaries = 1000)
    
    flux_interval$col = pal$colours
    
    flux_col = data.frame(fluxes = fluxes_ex$Flux) %>%
      left_join(flux_interval, by=c("fluxes" = "x"))
    
    E(g)$color<- flux_col$col
    
    # Define a layout to align nodes in three columns
    Metabolite_nodes <- unique(fluxes_ex$Metabolite)
    Producer_nodes <- unique(edges$Microorganism_type[grepl("\\(\\+\\)", edges$Microorganism_type)])
    Consumer_nodes <- unique(edges$Microorganism_type[grepl("\\(\\-\\)", edges$Microorganism_type)])
    
    # Define the layout
    layout <- matrix(NA, nrow = vcount(g), ncol = 2)
    # Assign positions for "Producer"
    layout[V(g)$name %in% Producer_nodes, 1] <- 1
    layout[V(g)$name %in% Producer_nodes, 2] <- seq(1, sum(V(g)$name %in% Producer_nodes))
    # Assign positions for "Metabolite"
    layout[V(g)$name %in% Metabolite_nodes, 1] <- 2
    layout[V(g)$name %in% Metabolite_nodes, 2] = 
      seq(1, sum(V(g)$name %in% Metabolite_nodes)) - round((length(Metabolite_nodes) - length(Producer_nodes))/2, 0)
    # Assign positions for "Consumer"
    layout[V(g)$name %in% Consumer_nodes, 1] <- 3
    layout[V(g)$name %in% Consumer_nodes, 2] <- seq(1, sum(V(g)$name %in% Consumer_nodes))
    
    pdf(file = paste0(wd, "/code/pictures/crossfeeding_plotting/int_net", "_", time_point, "_h", ".pdf"), 
        width = 8, height = 6)
    
    # Create a plot with a cleaner layout
    plot(g, 
         layout = layout, 
         frame = FALSE,
         vertex.color = V(g)$color, 
         vertex.size = V(g)$size,
         edge.color = E(g)$color, 
         edge.width = E(g)$width,
         vertex.label.cex = 0.6,
         edge.arrow.size = 0.6, 
         vertex.frame.color = "lightgrey",
         vertex.label.color = "black",
         vertex.label.family = "Helvetica",
         main = paste0("UnifiedGreatMod-predicted cross-feeding interactions; t = ", time_point, " [h]"))
    
    fields::image.plot(zlim = c(min(flux_interval$x), max(flux_interval$x)), 
                       col = pal$summary_colours, 
                       legend.only = TRUE,
                       horizontal = F,  # Vertical legend
                       legend.shrink = 0.25,  # Adjust size as needed
                       legend.width = 0.35,   # Adjust width as needed
                       axis.args = list(at = c(pal$summary_values[1], "0", pal$summary_values[1000]),
                                        labels = c(pal$summary_values[1], "0", pal$summary_values[1000]),
                                        cex.axis = 0.75, font.axis = 1),  # Adjust label font size
                       legend.line = 10,  # Adjust line position
                       # Adjust title font size and position
                       legend.args=list(text = "Flux Intensity [mmol/gDW*h]",
                                        col="black", cex = 0.8, side = 2, line=0.7))
    
    
    graphics.off()
    
  }
  
  pdf_files <- list.files(paste0(wd, "/code/pictures/crossfeeding_plotting"), 
                          pattern = "int_net_.*_h.pdf", 
                          full.names = TRUE)
  
  # Extract the frame numbers and order the files accordingly
  frame_order <- as.numeric(gsub(".*_(\\d+)_h\\.pdf$", "\\1", basename(pdf_files)))
  ordered_files <- pdf_files[order(frame_order)]
  
  # Initialize a list to store the images
  images <- list()
  
  # Convert each PDF to image and store it in the list
  for (file in ordered_files) {
    # Convert the first page of the PDF to an image
    img <- magick::image_read_pdf(file, density = 300) # Set the density for better resolution
    images <- c(images, img)
  }
  
  # Combine images into a GIF
  gif <- magick::image_animate(magick::image_join(images), fps = 100) # Set fps (frames per second)
  
  gganimate::anim_save(paste0(wd, "/code/pictures/crossfeeding_plotting.gif"), animation = gif)
  
}

# example usage 
plotting_int_net(wd = "~/Documents/FBA_project",
                 cross_fed_met = c(
                   "ac_e","ppa_e","but_e",
                   "M03134_e", "caproic_e", "isobut_e",
                   "isoval_e", "isocapr_e", "isobut_e",
                   "glu_L_e", "lac_L_e", "for_e", "ade_e",
                   "gua_e", "nac_e", "thymd_e", "ptrc_e"),
                 t_frame = 1:48,
                 col_bac = "#f06b368c",
                 col_met = "lightblue",
                 shape_bac = "square",
                 shape_met = "circle")
