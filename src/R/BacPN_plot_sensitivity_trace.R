
analyze_trace_sensitivity <- function(wd, model_name) {
  
  # Helper function to create individual trace plots
  create_trace_plot <- function(data, place_name, parameter_name = "param_sen") {
    ggplot(data[data$Places == place_name,], 
           aes(x = as.numeric(Time), 
               y = as.numeric(Marking), 
               group = config, 
               color = param_sen)) +
      geom_line(alpha = 0.7, linewidth = 0.5) +
      scale_color_gradient(low = "#178", 
                           high = "#907", 
                           name = "Starv") +
      labs(
        title = place_name,
        x = "Time",
        y = "Marking",
        color = "Parameter Value"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90"),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  # Main execution
  tryCatch({
    # Read the saved trace data
    sensitivity_trace <- readRDS(file = paste0(wd, "/", model_name, "_analysis/sensitivity_trace.rds"))
    
    # Combine all list elements into a single data frame
    combined_traces <- do.call(rbind, sensitivity_trace)
    
    # Get unique places
    unique_places <- unique(combined_traces$Places)
    
    # Create individual plots
    plots <- lapply(unique_places, function(place) {
      create_trace_plot(combined_traces, place)
    })
    
    # Return plots list
    return(plots)
    
  }, error = function(e) {
    stop(paste("Error in trace sensitivity analysis:", e$message))
  })
}