
# Solution 4: Custom spacing with clean aesthetics
create_alluvial_plot <- function(data, col_bac) {
  alluvial_data <- data %>%
    mutate(
      Reaction_Type = factor(Reaction_Type, levels = c("core", "boundary")),
      Species_Full = factor(Species_Full)
    )
  
  ggplot(alluvial_data, aes(axis1 = Species_Full, axis2 = Reaction_Type, 
                            axis3 = Reaction_Subtype, y = Count)) +
    geom_alluvium(aes(fill = Species_Full), 
                  alpha = 0.75, 
                  width = 1/12,
                  knot.pos = 0.25) +
    geom_stratum(alpha = 0.55, 
                 color = "black", 
                 width = 1/5,
                 linewidth = 0.2) +
    geom_text(stat = "stratum", 
              aes(label = after_stat(stratum)), 
              size = 3.5,
              # fontface = "bold",
              color = "black") +
    scale_x_discrete(limits = c("Taxa", "Type", "Subtype"), 
                     expand = c(0.4, 0.4)) +
    scale_fill_manual(labels = levels(alluvial_data$Species_Full), 
                      name = "Species", 
                      values = col_bac) +
    labs(title = NULL,
         y = "Number of Reactions") +
    theme_void() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold", 
                                 margin = margin(t = 15)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                margin = margin(b = 20)),
      legend.position = "right",
      plot.margin = margin(20, 60, 20, 60)
    ) +
    # Add y-axis label manually
    annotate("text", x = 0.5, y = max(alluvial_data$Count) * 2.25, 
             label = "Number of Reactions", angle = 90, 
             size = 4, fontface = "bold", color = "grey30")
}

create_sunburst_plot <- function(data) {
  # Prepare data for sunburst - create hierarchical paths
  sunburst_data <- data %>%
    mutate(
      # Create hierarchical path strings
      path = paste(Species_Full, Reaction_Type, Reaction_Subtype, Subsystem_Name, sep = "-"),
      # Clean subsystem names for better display
      Subsystem_Clean = gsub(" metabolism| degradation| biosynthesis", "", Subsystem_Name)
    ) %>%
    select(path, Count)
  
  # Create sunburst
  sunburst(sunburst_data, 
           count = TRUE,
           legend = list(w = 150, h = 25, s = 5, t = 25),
           colors = list(range = brewer.pal(8, "Set3")),
           legendOrder = NULL)
}

create_hierarchical_bars <- function(data) {
  
  plot_reactions_treemap_data <- data %>%
    group_by(Species_Full, Reaction_Type, Subsystem_Name) %>%
    summarise(Count = sum(Count), .groups = 'drop') %>%
    # Keep only top subsystems for readability
    group_by(Species_Full, Reaction_Type) %>%
    slice_max(Count, n = 18) %>%
    ungroup() %>%
    mutate(
      Subsystem_Short = case_when(
        nchar(Subsystem_Name) > 25 ~ paste0(substr(Subsystem_Name, 1, 22), "..."),
        TRUE ~ Subsystem_Name
      )
    )
  
  ggplot(plot_reactions_treemap_data, aes(x = reorder(Subsystem_Short, Count), y = Count, fill = Reaction_Type)) +
    geom_col(alpha = 0.8, color = "white", size = 0.3) +
    coord_flip() +
    facet_wrap(~Species_Full, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("core" = "#8B8B8B", "boundary" = "#F18F01"),
                      name = "Reaction Type") +
    labs(title = "Metabolic Subsystems",
         x = "Metabolic Subsystem",
         y = "Number of Reactions") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(b = 20)),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

create_metabolite_classification_plot <- function(data) {
  # Prepare classification data
  plot_data <- data %>%
    mutate(
      # Create meaningful classification labels
      Classification = case_when(
        is_core & !is_boundary ~ "Core Only",
        is_boundary & is_external_boundary ~ "External Boundary",
        is_boundary & is_internal_boundary ~ "Internal Boundary", 
        is_boundary & !is_external_boundary & !is_internal_boundary ~ "Boundary (Other)",
        TRUE ~ "Other"
      )
    ) %>%
    count(Species_Full, Classification, compartment, name = "Count") %>%
    # Create compartment labels
    mutate(
      Compartment_Label = case_when(
        compartment == "c" ~ "Cytoplasm",
        compartment == "e" ~ "External",
        compartment == "p" ~ "Periplasm",
        TRUE ~ compartment
      )
    )
  
  ggplot(plot_data, aes(x = Species_Full, y = Count, fill = Classification)) +
    geom_col(position = "stack", alpha = 0.8, color = "white", linewidth = 0.2) +
    facet_wrap(~Compartment_Label, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Core Only" = "#2E86AB",
        "External Boundary" = "#A23B72", 
        "Internal Boundary" = "#F18F01",
        "Boundary (Other)" = "#C73E1D",
        "Other" = "#8B8B8B"
      ),
      name = "Metabolite Type"
    ) +
    labs(
      title = "Metabolite Distribution",
      x = NULL,
      y = "Metabolites",
      subtitle = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(size = 11, margin = margin(b = 15)),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

create_summary_plot <- function(data, col_bac) {
  # Create summary statistics
  summary_data <- data %>%
    group_by(Species_Full) %>%
    summarise(
      Total_Metabolites = n(),
      Core_Only = sum(is_core & !is_boundary),
      External_Boundary = sum(is_boundary & is_external_boundary),
      Internal_Boundary = sum(is_boundary & is_internal_boundary),
      Cytoplasm = sum(compartment == "c"),
      External = sum(compartment == "e"),
      Periplasm = sum(compartment == "p", na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Convert to long format for plotting
    pivot_longer(cols = -c(Species_Full, Total_Metabolites), 
                 names_to = "Category", values_to = "Count") %>%
    mutate(
      Percentage = Count / Total_Metabolites * 100,
      Category_Type = ifelse(Category %in% c("Cytoplasm", "External", "Periplasm"), 
                             "Compartment", "Classification")
    )
  
  ggplot(summary_data, aes(x = Category, y = Count, fill = Species_Full)) +
    geom_col(position = "dodge", alpha = 0.8, color = "white", linewidth = 0.3) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3, fontface = "bold") +
    facet_wrap(~Category_Type, scales = "free_x") +
    scale_fill_manual(labels = levels(summary_data$Species_Full), name = "Species", values = col_bac) +
    labs(
      title = "Metabolites Summary",
      x = "Category",
      y = "Number of Metabolites",
      subtitle = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(size = 11, margin = margin(b = 15)),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    )
}

generate_metadata_plots <- function(hypernode_name, col_bac) {
  
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(sunburstR)
  library(ggalluvial)
  library(RColorBrewer)
  library(ggplot2)
  library(viridis)
  library(scales)

  config_yaml_path <- file.path(hypernode_name, "config", paste0(hypernode_name, ".yaml"))
  base_biounits_dir <- file.path(hypernode_name, "biounits")
  hypernode_config <- read_yaml(config_yaml_path) # Reads the YAML file
  
  # --- Initialize empty lists to store dataframes for all GEMs ---
  all_metabolites_list <- list()
  all_reactions_list <- list()
  
  # --- Dynamically Load Data for Each Cellular Unit (GEM) defined in the YAML ---
  # This loop iterates through each cellular unit defined in the YAML configuration.
  for (unit in hypernode_config$cellular_units) {
    model_name_dir <- unit$model_name # This is the specific folder name for the GEM (e.g., "Escherichia_coli_SE11")
    species_label <- unit$label       # This is the full species name for plots (e.g., "Escherichia coli SE11")
    metabolites_path <- file.path(base_biounits_dir, model_name_dir, "metabolites_metadata.csv")
    reactions_path <- file.path(base_biounits_dir, model_name_dir, "reactions_metadata.csv")
    
    # Check if files exist before attempting to read. If missing, a warning is issued, and the unit is skipped.
    if (!file.exists(metabolites_path) || !file.exists(reactions_path)) {
      warning(paste0("Metadata files not found for ", species_label, " (", species_label, "). Skipping this unit. Please ensure files are in: ", file.path(base_biounits_dir, model_name_dir)))
      next # Skip to the next unit in the loop
    }
    
    # Load metabolites_metadata.csv and add 'Species_Full' and 'species_label' columns.
    current_metabolites <- read_csv(metabolites_path) %>% 
      mutate(Species_Full = species_label, species_label = species_label) 
    
    # Load reactions_metadata.csv and add 'Species_Full' and 'species_label' columns.
    current_reactions <- read_csv(reactions_path) %>% 
      mutate(Species_Full = species_label, species_label = species_label) 
    
    # Store the processed dataframes in their respective lists, using the abbreviation as the list key.
    all_metabolites_list[[species_label]] <- current_metabolites
    all_reactions_list[[species_label]] <- current_reactions
  }
  
  # If no data could be loaded for any GEM, stop the function execution and inform the user.
  if (length(all_metabolites_list) == 0 || length(all_reactions_list) == 0) {
    stop("No GEM metadata could be loaded. Please check file paths and YAML configuration.")
  }
  
  all_metabolites <- bind_rows(all_metabolites_list)
  all_reactions <- bind_rows(all_reactions_list)
  
  ordered_species_labels <- sapply(hypernode_config$cellular_units, `[[`, "label")
  all_metabolites$Species_Full <- factor(all_metabolites$Species_Full, levels = ordered_species_labels)
  all_reactions$Species_Full <- factor(all_reactions$Species_Full, levels = ordered_species_labels)
  
  # USAGE EXAMPLES:
  # Most comprehensive overview:
  p1 <- create_metabolite_classification_plot(all_metabolites)
  ggsave(file.path("plots", "met_classification.pdf"), plot = p1, width = 5, height = 3)
  
  p2 <- create_summary_plot(data = all_metabolites, col_bac)
  ggsave(file.path("plots", "met_summary.pdf"), plot = p2, width = 4, height = 6)
  
  reactions_treemap_data <- all_reactions %>%
    # Ensure factors are correctly ordered for consistent treemap layout if needed.
    # The 'treemap' package typically uses its internal sorting if not explicitly factored.
    mutate(
      # 'type' column from all_reactions is used directly for Reaction_Type (Core/Boundary)
      Reaction_Type = type, 
      # 'subtype' column from all_reactions is used directly for Reaction_Subtype
      Reaction_Subtype = subtype,
      # 'subsystem' column from all_reactions is used directly for Subsystem_Name
      Subsystem_Name = subsystem
    ) %>%
    # Count reactions for each unique combination of these hierarchical levels
    count(Species_Full, Reaction_Type, Reaction_Subtype, Subsystem_Name, name = "Count") %>%
    # Filter out any zero counts, although 'count()' typically prevents this
    filter(Count > 0)
  
  p3 <- create_hierarchical_bars(reactions_treemap_data)
  ggsave(file.path("plots", "react_classification.pdf"), plot = p3, width = 6, height = 6.5)
  
  p4 <- create_alluvial_plot(reactions_treemap_data, col_bac)
  ggsave(file.path("plots", "react_alluv.pdf"), plot = p4, width = 7, height = 4)
  
  px = create_sunburst_plot(reactions_treemap_data)

}
