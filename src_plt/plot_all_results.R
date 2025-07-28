
plot_marking_and_flux_trends <- function(case_name, 
                                         col_bac,            # Vector of colors for bacteria
                                         col_met_places,     # Vector of colors for metabolites for bar plots
                                         met_to_plot,        # Vector of metabolite IDs for line/bar plots
                                         react2plot,         # Vector of reaction IDs to plot for flux trends
                                         num_sampling_points_rel_abun = 8,  # Number of time points for relative abundance bar plot
                                         num_sampling_points_met_plots = 8, # Number of time points for metabolite concentration bar plots
                                         plot_filename_suffix = ".pdf"      # Suffix for all output plot filenames
) {
  
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
  library(yaml)
  library(stringr)
  library(readr)
  
  # --- Internal Configuration Variables & Dynamic Paths ---
  # These are derived based on the 'case_name' argument
  base_dir <- file.path(".", case_name) 
  analysis_dir <- file.path(base_dir, paste0(case_name, "_analysis"))
  config_path <- file.path(base_dir, "config", paste0(case_name, ".yaml"))

  output_dir <- file.path("plots")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read configuration from YAML file
  config <- yaml::read_yaml(config_path)
  
  # Dynamically extract labels and model names from configuration
  label <- sapply(config$cellular_units, function(x) x$label)     # Short labels (e.g., "ecs", "cbd1")
  model_name <- sapply(config$cellular_units, function(x) x$model_name) # Full model names (e.g., "Escherichia_coli_SE11")
  met_to_plot <- config$boundary_metabolites 
  
  # --- Create species_metadata (Internal, uses 'col_bac' argument) ---
  species_metadata <- data.frame(
    label = label, 
    name_species = model_name, 
    col = col_bac, 
    stringsAsFactors = FALSE
  )
  
  # --- Load and Prepare Trace Data (Internal) ---
  trace_file <- file.path(analysis_dir, paste0(case_name, "-analysis-1.trace"))
  subtrace <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  # --- Load and Prepare Flux Data (Internal) ---
  flux_files <- list.files(analysis_dir, pattern = "\\.flux$", full.names = TRUE)
  
  list_of_dataframes <- list()
  for (f in flux_files) {
    current_label <- sub(paste0(".*-(", paste(label, collapse = "|"), ")_model\\.flux$"), "\\1", basename(f)) 
    subflux <- utils::read.table(f, header = TRUE) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      mutate(Microorganism = current_label) 
    list_of_dataframes[[f]] <- subflux
  }
  df_flux <- bind_rows(list_of_dataframes)
  
  # Filter df_flux to include only the reactions specified in 'react2plot' argument
  df_flux <- df_flux %>%
    filter(Reaction %in% react2plot) 
  
  # Ensure factor levels for 'Reaction' are consistent for facet_wrap plotting
  df_flux$Reaction <- factor(df_flux$Reaction, levels = unique(react2plot[react2plot %in% unique(df_flux$Reaction)]))
  
  # --- PLOTTING SECTION ---
  p_bac <- ggplot(subtrace %>% filter(Places %in% paste0("n_", label)),
                  aes(x = Time, y = Marking, color = Places, group = Places)) +
    # Border layer (slightly thicker, white/black)
    geom_line(linewidth = 1, color = "black") +  # White border
    geom_line(linewidth = 0.85) +
    scale_color_manual(labels = label, name = "Microorganism", values = col_bac) +
    theme_minimal(base_family = "Helvetica", base_size = 14) +
    labs(title=NULL, 
         x="Time (h)", y="Abundance (cell)", color="Taxa") +
    theme(plot.title = element_text(size=12, face="bold", margin=margin(b=10)),
          plot.subtitle = element_text(size=12, margin=margin(b=10)),
          axis.title = element_text(size=12), 
          axis.text = element_text(size=12),
          legend.title = element_text(size=10), 
          legend.text = element_text(size=10),
          legend.position = "right",
          panel.grid = element_blank(), axis.line = element_line(color = "black"))
  
  ggsave(p_bac, 
         file = file.path(output_dir, paste0("population_dynamics", plot_filename_suffix)), 
         width = 4.5, height = 3)
  
  ## 2. Plot Biomass Dynamics (Average per-cell biomass)
  p_biom <- ggplot(subtrace %>% filter(Places %in% paste0("biomass_e_", label)),
                   aes(x = Time, y = Marking, color = Places, group = Places)) +
    geom_line(linewidth = 1, color = "black") +  # White border
    geom_line(linewidth = 0.85) +
    scale_color_manual(labels=label, name = "Microorganism", values = col_bac) +
    theme_classic(base_size = 12) +
    labs(title=NULL, 
         x="Time (h)", y="Biomass (pg)", color="Taxa") +
    theme(plot.title = element_text(size=12, face="bold", margin=margin(b=10)),
          plot.subtitle = element_text(size=12, margin=margin(b=10)),
          axis.title = element_text(size=12), axis.text = element_text(size=14),
          legend.title = element_text(size=11), legend.text = element_text(size=11),
          legend.position = "right", axis.line = element_line(color = "black"))
  
  ggsave(p_biom, file = file.path(output_dir, paste0("biomass_dynamics", plot_filename_suffix)), width = 4.5, height = 3)
  
  ## 3. Plot Metabolite Concentrations (Extracellular pools)
  p_met <- ggplot(subtrace %>% filter(Places %in% met_to_plot),
                  aes(x = Time, y = Marking, color = Places, group = Places)) +
    geom_line(linewidth = 1, color = "black") +  # White border
    geom_line(linewidth=0.85) +
    scale_color_manual(name="Metabolite", values = col_met_places) + 
    theme_classic(base_size = 12) +
    labs(title=NULL, x="Time (h)", 
         y="Concentration (mmol/mL)", color="Metabolite") + 
    theme(plot.title = element_text(size=12, face="bold", margin=margin(b=10)),
          plot.subtitle = element_text(size=12, margin=margin(b=10)),
          axis.title = element_text(size=12), axis.text = element_text(size=14),
          legend.title = element_text(size=11), legend.text = element_text(size=11),
          legend.position = "right", axis.line = element_line(color = "black"))
  
  ggsave(p_met, file = file.path(output_dir, paste0("metabolite_concentrations", plot_filename_suffix)), width = 4.5, height = 3)
  
  all_unique_times <- unique(subtrace$Time)
  if (length(all_unique_times) > num_sampling_points_rel_abun) {
    indices_to_select <- round(seq(1, length(all_unique_times), length.out = num_sampling_points_rel_abun))
    sampling_times_rel_abun <- all_unique_times[indices_to_select]
  } else {
    sampling_times_rel_abun <- all_unique_times
  }
  relative_abundance <- subtrace %>% 
    filter(Places %in% paste0("n_", label)) %>%
    filter(Time %in% sampling_times_rel_abun) %>%
    group_by(Time) %>%
    mutate(total_marking = sum(Marking), rel_abundance = (Marking / total_marking) * 100) %>%
    mutate(Time = round(Time, 0)) %>%
    ungroup()
  
  p_rel_abun <- ggplot(relative_abundance, aes(x = as.factor(Time), y = rel_abundance, fill = Places)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(labels = paste0("n_", label), name = "Microorganism", values = col_bac) +
    labs(title = NULL, x = "Time (hours)", y = "Relative Abundance (%)") +
    theme_classic(base_size = 14, base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = 12, margin = margin(b = 10)),
          plot.subtitle = element_text(size = 12, margin=margin(b = 15)),
          axis.title.x = element_text(face = "bold", size = 9, margin = margin(t = 10)),
          axis.title.y = element_text(face = "bold", size = 9, margin = margin(r = 10)),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(face = "bold", size = 9),
          legend.text = element_text(size = 9),
          legend.position = "bottom", 
          legend.key.size = unit(0.4, "cm"),
          panel.grid.major = element_line(color = "grey80", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.8, color = "black"),
          axis.ticks = element_line(linewidth = 0.8),
          plot.margin = margin(20, 20, 20, 20))
  
  ggsave(p_rel_abun, file = file.path(output_dir, paste0("relative_abundance", plot_filename_suffix)), 
         width = length(sampling_times_rel_abun) * 0.25, height = 3)
  
  ## 5. Plot Metabolites Over Time (Bar Plots for selected time points)
  all_unique_times_met_plots <- unique(subtrace$Time)
  if (length(all_unique_times_met_plots) > num_sampling_points_met_plots) {
    indices_to_select <- round(seq(1, length(all_unique_times_met_plots), length.out = num_sampling_points_met_plots))
    sampling_times_met_plots <- all_unique_times_met_plots[indices_to_select]
  } else {
    sampling_times_met_plots <- all_unique_times_met_plots
  }
  barplots <- list()
  for (tp in sampling_times_met_plots) {
    tp_data <- subtrace %>% mutate(Time_rounded = round(Time, 0)) %>%
      filter(Time_rounded == round(tp, 0), Places %in% met_to_plot) # Use 'met_to_plot' argument
    if(nrow(tp_data) > 0) {
      p <- ggplot(tp_data, aes(x = Places, y = Marking, fill = Places)) +
        geom_bar(stat = "identity", position = "dodge", color = "darkgrey", width = 0.7) +
        scale_fill_manual(values = col_met_places) +
        theme_minimal(base_size = 14) +
        labs(title = paste("Time =", round(tp, 0), "h"), x = "Metabolite", y = "Concentration (mmol/mL)", fill = "Met") +
        theme(plot.title = element_text(size=14, face="bold", margin=margin(b=10)),
              axis.title = element_text(size=14), axis.text = element_text(size=12),
              legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=1),
              panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey80"),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + coord_flip()
      barplots[[as.character(tp)]] <- p
    }
  }
  if (length(barplots) > 0) {
    plot_cols_wrap <- ceiling(sqrt(length(barplots)))
    combined_plot <- wrap_plots(barplots, ncol = plot_cols_wrap) +
      plot_annotation(title = "Metabolites over time",
                      theme = theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 20), family = "Helvetica"),
                                    plot.margin = margin(10, 10, 10, 10)))
    ggsave(combined_plot, file = file.path(output_dir, paste0("metabolites_over_time_barplots", plot_filename_suffix)), width = plot_cols_wrap * 2.3, height = ceiling(length(barplots) / plot_cols_wrap) * 2.5)
  } else {
    message("No data points found for generating metabolite barplots. Check sampling_times_met_plots and filtering criteria.")
  }
  
  plot_cols_flux <- 4 # As per your snippet's variable
  num_facets_flux <- length(unique(df_flux$Reaction)) # Uses filtered df_flux
  plot_rows_flux <- ceiling(num_facets_flux / plot_cols_flux)
  
  p_ex <- ggplot(df_flux, aes(x = Time, y = Flux, color = Microorganism)) + # Color by Microorganism
    geom_line(linewidth = 0.85) +
    scale_color_manual(labels = label, name = "Microorganism", values = col_bac) +
    theme_minimal(base_size = 14) +
    labs(title="Dynamic FBA simulations", x="Time (h)", y="Metabolic flux (mmol/gDW*h)") +
    theme(legend.position = "right", legend.title = element_text(size = 12),
          legend.text = element_text(size = 12), plot.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 12), plot.subtitle = element_text(size = 14, margin=margin(b = 10)),
          axis.title = element_text(size = 16), panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          plot.margin = margin(15, 15, 15, 15)) +
    facet_wrap(~ Reaction, scales = "free_y", ncol = plot_cols_flux)
  
  ggsave(p_ex, file = file.path(output_dir, paste0("flux_simulations", plot_filename_suffix)), 
         width = plot_cols_flux * 2.5, height = plot_rows_flux * 1.5) 
}
