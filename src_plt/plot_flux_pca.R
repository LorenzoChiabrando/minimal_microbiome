plot_flux_pca <- function(case_name, flux_sampling_times, flux_th_l, flux_th_h, 
                          col_bac, scfa_list, entities_list, # New arguments
                          plot_filename = "pca_fluxomics_temporal_evolution.pdf") {
  
  # Ensure necessary packages are loaded.
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr) 
  library(readr)   
  library(yaml)    
  
  case_name <- case_name # This argument is now passed in
  base_dir <- file.path(".", case_name) # Base directory for the current case study
  analysis_dir <- file.path(base_dir, paste0(case_name, "_analysis"))
  config_path <- file.path(base_dir, "config", paste0(case_name, ".yaml"))
  
  config <- yaml::read_yaml(config_path)
  label <- sapply(config$cellular_units, function(x) x$label)     # Short labels (e.g., "ecs", "cbd1")
  
  model_name <- sapply(config$cellular_units, function(x) x$model_name) # Full model names (e.g., "Escherichia_coli_SE11")

  col_bac <- col_bac # Now passed as argument
  scfa_list <- scfa_list # Now passed as argument
  entities_list <- entities_list # Now passed as argument
  
  # --- Load Reaction Metadata (Internal) ---
  # Load reactions_metadata.csv for each microorganism and combine them
  all_reaction_metadata_list <- list() 
  for (m_name_full in model_name) { # Use full model names for constructing path
    reactions_meta_file <- file.path(base_dir, "biounits", m_name_full, "reactions_metadata.csv")
    if (file.exists(reactions_meta_file)) {
      all_reaction_metadata_list[[m_name_full]] <- read_csv(reactions_meta_file, show_col_types = FALSE)
    } else {
      warning(paste("Reactions metadata not found for", m_name_full, "at", reactions_meta_file))
    }
  }
  reaction_data_full <- bind_rows(all_reaction_metadata_list) %>%
    distinct(abbreviation, .keep_all = TRUE) # Ensure unique reactions across all models
  
  # --- Dynamically Determine react2plot from 'entities_list' argument ---
  # This section generates the 'react2plot_final' list based on the 'entities_list' passed as argument.
  reacts_for_plot_internal <- c() 
  for(entity_met_id in entities_list) {
    pattern <- paste0("\\b", entity_met_id, "\\b") 
    r_matched <- reaction_data_full$abbreviation[grepl(pattern, reaction_data_full$equation, ignore.case = TRUE)]
    reacts_for_plot_internal <- c(reacts_for_plot_internal, r_matched)
  }
  react2plot_final <- unique(reacts_for_plot_internal) 
  
  # --- Create species_metadata (Internal, uses 'col_bac' argument) ---
  species_metadata <- data.frame(
    label = label, # Short labels (e.g., "ecs", "cbd1")
    name_species = model_name, # Full names (e.g., "Escherichia_coli_SE11")
    col = col_bac, # Assigned colors from argument
    stringsAsFactors = FALSE
  )
  
  # --- Construct df_flux (combined flux data from all .flux files) (Internal) ---
  flux_files <- list.files(analysis_dir, pattern = "\\.flux$", full.names = TRUE)
  
  list_of_dataframes <- list()
  for (f in flux_files) {
    current_label <- sub(paste0(".*-(", paste(label, collapse = "|"), ")_model\\.flux$"), "\\1", basename(f)) 
    subflux <- utils::read.table(f, header = TRUE) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      mutate(Microorganism = current_label) # Add microorganism label for coloring
    list_of_dataframes[[f]] <- subflux
  }
  df_flux <- bind_rows(list_of_dataframes)
  
  # Filter df_flux to include only the reactions identified in 'react2plot_final'
  df_flux <- df_flux %>%
    filter(Reaction %in% react2plot_final) 
  
  # Ensure factor levels for 'Reaction' are consistent for facet_wrap plotting
  df_flux$Reaction <- factor(df_flux$Reaction, levels = unique(react2plot_final[react2plot_final %in% unique(df_flux$Reaction)]))
  
  final_time <- max(df_flux$Time)
  
  # Prepare data for PCA
  p_pca_data <- df_flux %>%
    mutate(Time_rounded = round(Time, 2)) %>%
    filter(Time_rounded %in% flux_sampling_times) %>% 
    filter(Flux >= flux_th_l & Flux <= flux_th_h) %>% 
    tidyr::unite("var", c(Time_rounded, Microorganism), sep = "&&", remove = FALSE) %>%
    select(-Time, -Microorganism) %>% 
    tidyr::spread(Reaction, Flux, fill = 0)
  
  # Set row names for PCA function
  p_pca_data_rownames <- p_pca_data %>% tibble::column_to_rownames("var")
  
  # Filter out columns (reactions) with zero standard deviation
  sd_check <- apply(p_pca_data_rownames, 2, sd)
  valid_cols <- names(sd_check[sd_check != 0])
  
  # Perform PCA
  if (length(valid_cols) > 1 && nrow(p_pca_data_rownames) > 1) { 
    q <- prcomp(p_pca_data_rownames[, valid_cols], scale. = TRUE, center = TRUE)
    
    # Extract PCA scores for plotting.
    pca_plot_df <- q$x[, 1:2] %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("var") %>% 
      tidyr::separate(var, into = c("Time", "Microorganism"), sep = "&&") %>% 
      mutate(Time = as.numeric(Time)) %>% 
      arrange(Microorganism, Time) %>% # Order for geom_path trajectories
      left_join(species_metadata, by = c("Microorganism" = "label")) # Join with species_metadata
    
    # Calculate explained variance for plot labels.
    explained_variance <- summary(q)$importance[2,] * 100
    
    # Identify start and end time points for highlighting.
    start_end_points <- pca_plot_df %>%
      group_by(Microorganism) %>%
      filter(Time == min(Time) | Time == max(Time)) %>%
      mutate(PointType = ifelse(Time == min(Time), "Start", "End")) %>%
      ungroup()
    
    pca_plot <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Microorganism)) +
      # Border for paths (white border layer)
      geom_path(aes(group = Microorganism, alpha = Time), 
                linewidth = 1.0, lineend = "round", color = "black") +
      # Main paths
      geom_path(aes(group = Microorganism, alpha = Time), 
                linewidth = 0.8, lineend = "round") +
      
      # Border for regular points (white border)
      geom_point(size = 2.4, color = "black") +
      # Main regular points
      geom_point(size = 2) +
      
      # Border for start/end points (white border)
      geom_point(data = start_end_points, aes(x = PC1, y = PC2, shape = PointType), 
                 size = 4.4, stroke = 1.5, color = "black") +
      # Main start/end points
      geom_point(data = start_end_points, aes(x = PC1, y = PC2, shape = PointType), 
                 size = 4, stroke = 1) +
      
      scale_color_manual(labels = unique(species_metadata$label), 
                         name = "Microorganism",
                         values = unique(species_metadata$col)) +
      scale_shape_manual(name = "Time Point Type", values = c(Start = 15, End = 17)) + 
      scale_alpha_continuous(range = c(0.4, 1), guide = "none") + 
      
      labs(title = "PCA of Fluxomics",
           subtitle = paste("Trajectories and Temporal Evolution: Features (reactions) used:", length(valid_cols)),
           x = paste0("PC1 (", round(explained_variance[1], 2), "%)"),
           y = paste0("PC2 (", round(explained_variance[2], 2), "%)")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right", 
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'dotted', color = "gray"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
    # Save the PCA plot
    output_dir <- file.path("plots")
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    ggsave(pca_plot, file = file.path(output_dir, plot_filename), width = 5, height = 3.75)
    
    message(paste("PCA plot saved to:", file.path(output_dir, plot_filename)))
  } else {
    message("PCA plot could not be generated due to insufficient valid data after filtering (e.g., too few non-constant reactions or samples).")
  }
}
