# Ensure necessary packages are installed and loaded
# install.packages(c("dplyr", "igraph", "colourvalues", "stringr", "readr", "magick", "fields", "yaml"))
library(dplyr)
library(igraph)
library(colourvalues)
library(stringr)
library(readr)
library(magick)
library(fields) # For image.plot legend
library(yaml) # For reading YAML config

# Helper function to make a dataframe of reactions irreversible
# This function is crucial for aligning reaction data with _f/_r fluxes.
make_reactions_irreversible <- function(df_reactions_meta) {
  df_irreversible <- tibble() 
  for (i in 1:nrow(df_reactions_meta)) {
    r <- df_reactions_meta[i, ] 
    is_reversible <- (r$lowbnd < 0 && r$uppbnd > 0) || str_detect(r$equation, "<=>")
    
    if (is_reversible) {
      forward_abbr <- paste0(r$abbreviation, "_f")
      forward_eq <- str_replace(r$equation, "<=>", "->")
      forward_row <- r %>% mutate(abbreviation = forward_abbr, equation = forward_eq, lowbnd = 0)
      
      reverse_abbr <- paste0(r$abbreviation, "_r")
      eq_parts <- str_split(r$equation, "<=>")[[1]]
      if (length(eq_parts) == 2) {
        reverse_eq <- paste(str_trim(eq_parts[2]), "->", str_trim(eq_parts[1]))
      } else {
        reverse_eq <- r$equation 
        warning(paste("Complex reversible equation for", r$abbreviation, ". Cannot reliably reverse:", r$equation))
      }
      reverse_row <- r %>% mutate(abbreviation = reverse_abbr, equation = reverse_eq, lowbnd = 0, uppbnd = r$uppbnd) 
      
      df_irreversible <- bind_rows(df_irreversible, forward_row, reverse_row)
    } else { 
      if (r$lowbnd < 0 && r$uppbnd <= 0) {
        df_irreversible <- bind_rows(df_irreversible, r %>% mutate(abbreviation = paste0(r$abbreviation, "_r")))
      } else { 
        df_irreversible <- bind_rows(df_irreversible, r %>% mutate(abbreviation = paste0(r$abbreviation, "_f")))
      }
    }
  }
  return(df_irreversible)
}


# Define the main plotting function for cross-feeding interaction networks
plotting_int_net <- function(
    case_name,           # String: Name of the case study folder (e.g., "EcCb_SCFA")
    t_frame,             # Numeric vector: Time points (in hours) at which to generate network snapshots
    col_bac,             # String: Single color (hex code or name) for ALL bacteria nodes (e.g., "blue")
    col_met,             # String: Single color (hex code or name) for ALL metabolite nodes (e.g., "lightblue")
    shape_bac = "square", # String: Shape for bacteria nodes (set to square as per image analysis)
    shape_met = "circle", # String: Shape for metabolite nodes (set to circle as per image analysis)
    cross_fed_met,       # String vector: List of metabolite IDs to consider in the cross-feeding network
    min_flux_threshold = 1e-6, # Numeric: Minimum absolute flux to draw an edge
    fixed_edge_width = 2, # Numeric: Fixed width for all edges (as per user's latest correction)
    output_filename_prefix = "int_net", # String: Prefix for output PDF/GIF files
    plot_filename_suffix = ".pdf" # String: Suffix for individual plot files (e.g., ".pdf", ".png")
) {
  
  # --- Internal Paths & Config Loading ---
  base_dir <- file.path(".", case_name)
  analysis_dir <- file.path(base_dir, paste0(case_name, "_analysis"))
  config_path <- file.path(base_dir, "config", paste0(case_name, ".yaml"))
  
  output_plot_dir <- file.path(base_dir, "plots") 
  frames_output_dir <- file.path(output_plot_dir, "crossfeeding_plotting_frames")
  if (!dir.exists(frames_output_dir)) {
    dir.create(frames_output_dir, recursive = TRUE)
  }
  
  config <- yaml::read_yaml(config_path)
  
  label <- sapply(config$cellular_units, function(x) x$label)
  model_name <- sapply(config$cellular_units, function(x) x$model_name)
  
  # --- Load and Prepare Flux Data (All fluxes) ---
  flux_files <- list.files(analysis_dir, pattern = "\\.flux$", full.names = TRUE)
  
  list_of_raw_flux_dataframes <- list()
  for (f in flux_files) {
    current_label <- sub(paste0(".*-(", paste(label, collapse = "|"), ")_model\\.flux$"), "\\1", basename(f))
    current_full_name <- model_name[match(current_label, label)] 
    
    subflux_raw <- utils::read.table(f, header = TRUE) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      mutate(Microorganism_Label = current_label,      
             Microorganism_FullName = current_full_name) 
    list_of_raw_flux_dataframes[[f]] <- subflux_raw
  }
  df_flux_raw <- bind_rows(list_of_raw_flux_dataframes)
  
  # --- Load Trace Data (Internal) for Node Sizing ---
  trace_file <- file.path(analysis_dir, paste0(case_name, "-analysis-1.trace"))
  subtrace <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  # --- Load Reaction Metadata (Internal) ---
  all_reaction_metadata_list <- list()
  for (m_name_full in model_name) {
    reactions_meta_file <- file.path(base_dir, "biounits", m_name_full, "reactions_metadata.csv")
    if (file.exists(reactions_meta_file)) {
      all_reaction_metadata_list[[m_name_full]] <- read_csv(reactions_meta_file, show_col_types = FALSE)
    } else {
      warning(paste("Reactions metadata not found for", m_name_full, "at", reactions_meta_file))
    }
  }
  reaction_data_full <- bind_rows(all_reaction_metadata_list) %>%
    distinct(abbreviation, .keep_all = TRUE) %>%
    filter(type == "boundary") # Filter to only boundary reactions as per user's earlier snippet
  
  reaction_data_full_irreversible <- make_reactions_irreversible(reaction_data_full)
  
  # NEW: Add 'Metabolite' and 'Equation' columns to df_flux_raw by joining with the irreversible reaction data.
  df_flux_raw <- df_flux_raw %>%
    left_join(reaction_data_full_irreversible %>% select(Reaction = abbreviation, Equation = equation, Subtype = subtype), 
              by = "Reaction") %>%
    mutate(Metabolite = ifelse(Subtype == "exchange", 
                               str_extract(Equation, "\\b\\d*\\s*([a-zA-Z0-9_]+_e)\\b"), 
                               NA_character_), 
           Metabolite = str_remove(Metabolite, "^\\d+\\s*"),
           Metabolite = str_trim(Metabolite)
    ) %>%
    select(-Equation, -Subtype) 
  
  # --- Flux Transformation for effective_flux_for_color ---
  df_flux_raw <- df_flux_raw %>%
    mutate(
      cleaned_flux = ifelse(Flux < 0, 0, Flux) # Set numerical errors (negative fluxes) to 0
    ) %>%
    mutate(
      effective_flux_for_color = case_when(
        str_detect(Reaction, "_f$") ~ cleaned_flux,  # Forward/Outflow: Use the cleaned value as positive
        str_detect(Reaction, "_r$") ~ -cleaned_flux, # Reverse/Inflow: Make the cleaned value negative
        TRUE ~ cleaned_flux # Fallback
      )
    ) %>%
    select(-cleaned_flux) 
  
  # --- Determine Global Flux Min/Max for Consistent Color Scale ---
  relevant_exchanges_for_coloring <- c()
  for(met_id in cross_fed_met) {
    pattern <- paste0("\\b", met_id, "\\b")
    matched_reactions <- reaction_data_full_irreversible %>%
      filter(grepl(pattern, equation, ignore.case = TRUE) & subtype == "exchange") %>%
      pull(abbreviation)
    relevant_exchanges_for_coloring <- c(relevant_exchanges_for_coloring, matched_reactions)
  }
  relevant_exchanges_for_coloring <- unique(relevant_exchanges_for_coloring)
  
  relevant_fluxes_for_range <- df_flux_raw %>%
    filter(Reaction %in% relevant_exchanges_for_coloring) %>%
    pull(effective_flux_for_color)
  
  if(length(relevant_fluxes_for_range) > 0) {
    global_min_flux <- min(relevant_fluxes_for_range)
    global_max_flux <- max(relevant_fluxes_for_range)
  } else {
    global_min_flux <- -1 
    global_max_flux <- 1
  }
  abs_max_val <- max(abs(global_min_flux), abs(global_max_flux))
  z_lim_sym <- c(-abs_max_val, abs_max_val)
  
  
  # --- Loop through time points to generate individual network plots ---
  all_graph_pdf_paths <- list() 
  
  for (time_point in t_frame) {
    actual_time_points_in_data <- unique(df_flux_raw$Time)
    if (length(actual_time_points_in_data) == 0) {
      message(paste("df_flux_raw has no time points. Skipping plot for time =", time_point, "h."))
      next
    }
    closest_data_time_point_index <- which.min(abs(actual_time_points_in_data - time_point))
    closest_data_time_point <- actual_time_points_in_data[closest_data_time_point_index]
    
    fluxes_ex <- df_flux_raw %>%
      filter(dplyr::near(Time, closest_data_time_point, tol = 1e-9)) %>% 
      filter(Reaction %in% relevant_exchanges_for_coloring) %>%          
      filter(abs(effective_flux_for_color) >= min_flux_threshold) 
    
    if (nrow(fluxes_ex) == 0) {
      message(paste("No significant cross-feeding fluxes at time =", time_point, "h (closest actual time in data was", closest_data_time_point, "h). Skipping plot for this frame."))
      next
    }
    
    # --- Prepare Edges and Nodes for igraph based on the reference picture's layout ---
    # Edges: define 'from' and 'to' based on '_f' (outflow) or '_r' (inflow) suffix in Reaction ID.
    edges_df <- fluxes_ex %>%
      mutate(
        from_node = case_when(
          str_detect(Reaction, "_f$") ~ Microorganism_FullName, # Outflow: Microorganism is source
          str_detect(Reaction, "_r$") ~ Metabolite,             # Inflow: Metabolite is source
          TRUE ~ NA_character_ 
        ),
        to_node = case_when(
          str_detect(Reaction, "_f$") ~ Metabolite,             # Outflow: Metabolite is target
          str_detect(Reaction, "_r$") ~ Microorganism_FullName, # Inflow: Microorganism is target
          TRUE ~ NA_character_ 
        )
      ) %>%
      select(from_node, to_node, Flux = effective_flux_for_color, Metabolite, Reaction, Microorganism_FullName, Microorganism_Label) 
    
    # Nodes: Create unique nodes and assign types, and add roles for labels
    nodes_df <- data.frame(name = unique(c(edges_df$from_node, edges_df$to_node))) %>%
      mutate(type = ifelse(name %in% model_name, "Microorganism", "Metabolite")) %>%
      # Add role-based suffix for microorganism labels: (+) for producer, (-) for consumer
      mutate(role_suffix = case_when(
        name %in% unique(edges_df$from_node[str_detect(edges_df$Reaction, "_f$")]) & type == "Microorganism" ~ "(+)", 
        name %in% unique(edges_df$to_node[str_detect(edges_df$Reaction, "_r$")]) & type == "Microorganism" ~ "(-)", 
        TRUE ~ "" 
      ))
    
    # Create the graph object
    g <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = T)
    
    # --- Set Node Attributes ---
    V(g)$color <- ifelse(V(g)$type == "Microorganism", col_bac, col_met) # Use single color arguments
    V(g)$shape <- ifelse(V(g)$type == "Microorganism", shape_bac, shape_met) 
    V(g)$label.cex <- 0.8
    V(g)$label.font <- 1
    V(g)$label.color <- "black"
    V(g)$frame.color = "lightgrey" 
    # Set node labels as "Short Label (+/-)" for microorganisms, and ID for metabolites
    V(g)$label <- sapply(V(g)$name, function(node_name) {
      node_info <- nodes_df %>% filter(name == node_name)
      if (node_info$type == "Microorganism") {
        paste0(label[match(node_info$name, model_name)], node_info$role_suffix) 
      } else {
        node_name # Use metabolite ID
      }
    })
    
    # NEW: --- Determine Node Sizes Based on Marking Values (from subtrace) ---
    # 1. Filter subtrace for the current `closest_data_time_point`.
    current_markings <- subtrace %>% 
      filter(dplyr::near(Time, closest_data_time_point, tol = 1e-9))
    
    # 2. Extract marking for each node in the graph 'g'
    node_sizes_raw <- sapply(V(g)$name, function(node_name) {
      node_type <- nodes_df %>% filter(name == node_name) %>% pull(type)
      if (node_type == "Microorganism") {
        mic_label_short <- label[match(node_name, model_name)] 
        # Prioritize population size ('n_') for scaling, then biomass ('biomass_e_').
        pop_marking <- current_markings %>% filter(Places == paste0("n_", mic_label_short)) %>% pull(Marking)
        if (length(pop_marking) > 0 && pop_marking[1] > 0) return(pop_marking[1])
        
        biomass_marking <- current_markings %>% filter(Places == paste0("biomass_e_", mic_label_short)) %>% pull(Marking)
        if (length(biomass_marking) > 0 && biomass_marking[1] > 0) return(biomass_marking[1])
        
        return(NA_real_) 
      } else if (node_type == "Metabolite") {
        # For metabolites, use their concentration. Node name is already the metabolite ID.
        met_marking <- current_markings %>% filter(Places == node_name) %>% pull(Marking)
        if (length(met_marking) > 0) return(met_marking[1])
        return(NA_real_) 
      }
      return(NA_real_) 
    })
    
    # Scale raw marking values to appropriate visual sizes.
    # Use log1p(x) for log(1+x) which handles zero values gracefully.
    # Scale all sizes relative to the maximum observed log marking to fit a visual range.
    base_node_visual_size <- 5 # A minimum size for visibility
    max_log_marking <- max(log1p(node_sizes_raw[!is.na(node_sizes_raw) & node_sizes_raw > 0]), na.rm = TRUE)
    if (is.infinite(max_log_marking) || max_log_marking == 0) max_log_marking <- 1 # Avoid division by zero/infinity
    
    V(g)$size <- ifelse(is.na(node_sizes_raw) | node_sizes_raw <= 0, 
                        base_node_visual_size, # Default size for inactive/zero marking
                        base_node_visual_size + (log1p(node_sizes_raw) / max_log_marking) * 15 # Scale to a range (e.g., base + up to 15 more)
    )
    
    # Set Edge Attributes (color based on effective_flux_for_color, width configurable)
    # Assuming 'edges_df$Flux' contains the 'effective_flux_for_color' calculated previously.
    
    # Apply manual scaling of Flux values to the symmetrical range z_lim_sym for consistent coloring.
    # This approach ensures the color mapping is comparable across all time frames and centered at zero.
    scaled_flux_for_palette <- (edges_df$Flux - z_lim_sym[1]) / (z_lim_sym[2] - z_lim_sym[1])
    
    # Clamp values to [0, 1] range to handle any numerical edge cases or values slightly outside z_lim_sym.
    # This is important before passing to a palette function that expects values in a normalized range.
    scaled_flux_for_palette <- pmax(0, pmin(1, scaled_flux_for_palette))
    
    # Assign edge colors using colourvalues::colour_values.
    # FIX: The 'limits' argument is removed as it caused an error. Manual scaling before this step handles the range.
    # And explicitly select the '$colours' component from the returned list.
    E(g)$color <- colourvalues::colour_values(scaled_flux_for_palette, palette = "RdBu", n_summaries = 1000)$colours
    E(g)$width <- 1 # Fixed width for edges, as per user's latest correction
    E(g)$arrow.size <- 0.6 
    
    # --- Define a Layout (Strict 3-Column Layout: Producers | Metabolites | Consumers) ---
    layout_matrix <- matrix(NA, nrow = vcount(g), ncol = 2)
    rownames(layout_matrix) <- V(g)$name 
    
    active_producers_in_frame <- nodes_df %>% filter(role_suffix == "(+)") %>% pull(name)
    active_consumers_in_frame <- nodes_df %>% filter(role_suffix == "(-)") %>% pull(name)
    active_metabolite_nodes_in_graph <- nodes_df %>% filter(type == "Metabolite") %>% pull(name)
    
    layout_matrix[V(g)$name %in% model_name, 1] <- 1 
    layout_matrix[active_metabolite_nodes_in_graph, 1] <- 2
    layout_matrix[active_consumers_in_frame, 1] <- 3
    
    temp_layout_df <- data.frame(
      name = V(g)$name,
      x_pos = layout_matrix[V(g)$name, 1], 
      stringsAsFactors = FALSE
    ) %>%
      group_by(x_pos) %>%
      arrange(name) %>% 
      mutate(y_pos = seq(0, 1, length.out = n())) %>% 
      ungroup()
    
    layout_matrix[temp_layout_df$name, 2] <- temp_layout_df$y_pos # Corrected assignment
    
    current_layout <- layout_matrix[V(g)$name, ] 
    
    # --- Plotting to PDF for the current frame ---
    frame_pdf_path <- file.path(frames_output_dir, paste0(output_filename_prefix, "_", time_point, "_h", plot_filename_suffix))
    pdf(file = frame_pdf_path, width = 8, height = 6)
    
    plot(g, 
         layout = current_layout, 
         frame = FALSE,
         vertex.color = V(g)$color, 
         vertex.size = V(g)$size,
         edge.color = E(g)$color, 
         edge.width = E(g)$width,
         vertex.label = V(g)$label, 
         vertex.label.cex = V(g)$label.cex,
         edge.arrow.size = E(g)$arrow.size, 
         vertex.frame.color = V(g)$frame.color,
         vertex.label.color = V(g)$label.color,
         vertex.label.family = V(g)$label.family,
         main = paste0("UnifiedGreatMod-predicted cross-feeding interactions; t = ", time_point, " [h]"))
    
    fields::image.plot(zlim = z_lim_sym, 
                       col = colourvalues::get_palette("RdBu", n = 1000), 
                       legend.only = TRUE, horizontal = FALSE, 
                       legend.shrink = 0.5, legend.width = 0.35, 
                       axis.args = list(at = c(round(z_lim_sym[1], 2), 0, round(z_lim_sym[2], 2)), cex.axis = 0.75, font.axis = 1),
                       legend.args=list(text = "Flux [mmol/gDW*h]", col="black", cex = 0.8, side = 2, line=0.7))
    
    graphics.off()
    all_graph_pdf_paths[[as.character(time_point)]] <- frame_pdf_path 
  }
  
  # --- GIF Generation ---
  pdf_files_for_gif <- unlist(all_graph_pdf_paths)
  
  if (length(pdf_files_for_gif) > 0 && requireNamespace("magick", quietly = TRUE)) {
    message("Attempting to generate GIF animation. Requires Ghostscript installed (for PDF conversion).")
    images <- list()
    for (file_path in pdf_files_for_gif) {
      img <- tryCatch({
        magick::image_read_pdf(file_path, density = 300)
      }, error = function(e) {
        warning(paste("Could not read PDF for GIF (", file_path, "):", e$message))
        return(NULL)
      })
      if (!is.null(img)) {
        images <- c(images, list(img)) 
      }
    }
    
    if (length(images) > 0) {
      gif <- magick::image_animate(magick::image_join(images), fps = 2) 
      gif_path <- file.path(output_plot_dir, paste0(output_filename_prefix, "_crossfeeding.gif"))
      magick::image_write(gif, path = gif_path)
      message(paste("GIF animation saved to:", gif_path))
    } else {
      message("No images successfully processed for GIF. Check warnings above.")
    }
  } else {
    message("GIF animation skipped. No PDF frames generated or 'magick' package not available/Ghostscript not installed.")
  }
}