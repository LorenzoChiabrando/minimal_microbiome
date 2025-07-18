
plt_ana <- function(fba_name, 
                    reactions_of_interest,
                    place2plot) {
  
  library(ggplot2)
  
  # Create a custom publication theme
  custom_theme <- function(base_size) {
    theme_minimal(base_size = base_size) %+replace%
      theme(
        # Text elements
        text = element_text(color = "black"),
        plot.title = element_text(
          size = rel(1.3),
          face = "bold",
          margin = margin(b = 15),
          hjust = 0
        ),
        plot.subtitle = element_text(
          size = rel(1.1),
          margin = margin(b = 10),
          hjust = 0
        ),
        axis.title = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(0.9)),
        
        # Legend formatting
        legend.position = "bottom",
        legend.title = element_text(size = rel(0.9)),
        legend.text = element_text(size = rel(0.8)),
        legend.box.spacing = unit(0.5, "lines"),
        
        # Grid lines and panel
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA, linewidth = 0.5),
        
        # Facet formatting
        strip.text = element_text(
          size = rel(0.9),
          face = "bold",
          margin = margin(b = 5, t = 5)
        ),
        strip.background = element_rect(fill = "grey95", color = NA),
        
        # Plot margins
        plot.margin = margin(15, 15, 15, 15),
        
        complete = TRUE
      )
  }
  
  base_size = 12
  
  path = paste0(hypernode_name, "_analysis")

  trace_file <- list.files(path, 
                           pattern = "\\.trace$", 
                           full.names = TRUE)
  
  flux_file <- list.files(path, 
                          pattern = "\\.flux$", 
                          full.names = TRUE)
  
  if (length(trace_file) == 0 || length(flux_file) == 0) {
    stop("No .trace or .flux files found and none provided")
  }
  
  trace <- utils::read.table(trace_file, header = TRUE) %>%
    tidyr::gather(key = "Places", value = "Marking", -Time) %>%
    # Convert Places to factor to maintain order
    dplyr::mutate(Places = factor(Places, levels = unique(Places)))
  
  critical_times <- trace %>%
    dplyr::filter(Places == paste0("n_", fba_name) & Marking < 1) %>%
    dplyr::pull(Time)
  
  flux <- utils::read.table(flux_file, header = TRUE) %>%
    tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
    dplyr::mutate(Organism = fba_name) %>% 
    dplyr::mutate(Flux = dplyr::if_else(Time %in% critical_times, 0, Flux))
  
  # Assuming place2plot is your vector of place names to plot
  variable_names <- setNames(
    list(
      paste0(levels(trace$Places)[1], " (cell)"),
      paste0(levels(trace$Places)[2], " (pgDW)")
    ),
    c(paste0("n_", fba_name), paste0(fba_name, "_biomass_e"))
  )
  
  p_pl <- ggplot2::ggplot(trace, ggplot2::aes(Time, Marking)) +
    # Using a subtle, yet distinct line color for the primary trace
    ggplot2::geom_line(linewidth = 0.8, color = "#2c7bb6") + # A clean, scientific blue
    # Using slightly larger points with a softer outline and increased transparency for better visual flow
    ggplot2::geom_point(size = 1.5, shape = 21, fill = "#ab9", color = "white", alpha = 0.7) + # Light blue fill, white outline
    # Adding a shaded ribbon for visual emphasis of the trend, using a complementary light blue
    geom_ribbon(aes(ymin = min(Marking), ymax = Marking), alpha = 0.15, fill = "#a6bddb") + # A slightly darker, muted blue for the ribbon
    
    ggplot2::labs(title = NULL, x = "Time (h)", y = "Marking") +
    ggplot2::facet_grid(Places ~ ., scales = "free_y") +
  
    # Adopting a more refined minimal theme with adjusted base font size
    theme_minimal(base_size = 14) + # Increased base font size for better readability in papers
    theme(
      text = element_text(color = "black", family = "Arial"), # Specify a common, professional font like Arial
      plot.title = element_text(size = rel(1.3), face = "bold", hjust = 0.5, color = "#333333"), # Darker grey for title
      axis.title = element_text(size = rel(1.15), face = "bold", color = "#333333"), # Darker grey for axis titles
      axis.text = element_text(size = rel(1), color = "#555555"), # Softer grey for axis text
      
      # Fine-tuning grid lines for a cleaner look
      panel.grid.major = element_line(color = "gray90", linewidth = 0.1), # Slightly thicker, still subtle
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.25), # Slightly thicker, very subtle
      
      # Stripping text and background for facet labels
      strip.text = element_text(face = "bold", size = rel(1.05), color = "white"), # White text for better contrast
      strip.background = element_rect(fill = "#5e81a6", color = NA), # A darker, professional blue for facet backgrounds
      
      panel.spacing = unit(1.5, "lines"), # Reduced spacing slightly for a more compact look
      legend.position = "bottom",
      legend.title = element_text(face = "bold", color = "#333333"),
      legend.text = element_text(color = "#555555"),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15), # Slightly reduced margins
      
      # Adding a subtle light grey background to the entire plot area
      plot.background = element_rect(fill = "#f9f9f9", color = NA) 
    )
  
  # Process the data as before
  flux_filtered <- flux %>%
    dplyr::filter(Reaction %in% reactions_of_interest) %>%
    dplyr::mutate(
      Type = ifelse(grepl("_f$", Reaction), "forward", "reverse"),
      Base_Reaction = sub("_[fr]$", "", Reaction)
    )
  
  # Create list to store plots
  plot_list <- list()
  
  # Create plots for each base reaction
  # Ensure plot_list is initialized outside the loop if it's not already
  # plot_list <- list() 
  
  for(reaction in unique(sub("_[fr]$", "", reactions_of_interest))) {
    
    message(paste("Attempting to plot for Base_Reaction:", reaction)) # Debugging message
    
    f_data <- flux_filtered %>% dplyr::filter(Base_Reaction == reaction)
    
    if (nrow(f_data) == 0) {
      warning(paste("No data for Base_Reaction '", reaction, "'. Skipping plot.", sep = ""))
      next # Skip to the next iteration of the loop
    }
    
    # --- NEW CHECK FOR FACETING VARIABLE ---
    unique_facet_vars <- unique(f_data$Reaction)
    if (length(unique_facet_vars) == 0 || all(is.na(unique_facet_vars))) {
      warning(paste("Faceting variable 'Reaction' has no valid values for Base_Reaction '", reaction, "'. Skipping plot.", sep = ""))
      next # Skip if no valid facet levels
    }
    # --- END NEW CHECK ---
    
    # Get maximum y value for this reaction pair
    # Use f_data directly here as it's already filtered
    y_max <- f_data %>%
      dplyr::pull(Flux) %>%
      max(na.rm = TRUE)
    
    # Ensure y_max is not -Inf if all values were negative or missing
    if (is.infinite(y_max) && y_max < 0) {
      y_max <- 0 # Or some appropriate small positive value if fluxes must be non-negative
    } else if (is.na(y_max)) {
      y_max <- 1 # Default to 1 if max is NA (e.g., if Flux column only had NAs)
    }
    
    # Determine minimum y value for coord_cartesian to correctly display negative fluxes if any
    y_min <- f_data %>%
      dplyr::pull(Flux) %>%
      min(na.rm = TRUE)
    
    # Handle cases where y_min might be infinite or NA
    if (is.infinite(y_min) && y_min > 0) {
      y_min <- 0 # If min is +Inf (e.g., from empty data with min=TRUE), default to 0
    } else if (is.na(y_min)) {
      y_min <- 0 # Default if min is NA
    }
    
    # Ensure the lower limit for coord_cartesian is 0, or below 0 if min flux is genuinely negative
    lower_ylim <- min(0, y_min) # This correctly sets 0 as the lower bound unless actual data goes below 0
    
    plot_list[[reaction]] = ggplot(f_data, aes(x = Time, y = Flux, color = Organism)) +
      geom_line(linewidth = 1) +
      coord_cartesian(ylim = c(lower_ylim, y_max)) + # Use the dynamically determined lower_ylim
      labs(x = "Time (h)", y = "Flux (mmol/gDW*h)") +
      facet_wrap(~Reaction, ncol = 1) +
      custom_theme(base_size = base_size) +
      theme(
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "bottom"
      )
  }
  
  # Combine plots using patchwork
  final_plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = NULL,
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
    ) & theme(legend.position = "bottom")
  
  # Return enhanced plots
  return(list(
    places = p_pl,
    exchange = final_plot
  ))
}
