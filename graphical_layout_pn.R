# ─────────────────────────────────────────────────────────────────────────────
# 7. Graphical Layout Analysis
# ─────────────────────────────────────────────────────────────────────────────
step_num <- step_num + 1
log_section(sprintf("%d. Graphical Layout Analysis", step_num))

# Extract position information
place_positions <- xml_find_all(xml_content, "//place") %>%
  map_df(~{
    attrs <- xml_attrs(.x)
    tibble(
      name = attrs["name"],
      x = as.numeric(attrs["x"]),
      y = as.numeric(attrs["y"])
    )
  })

transition_positions <- xml_find_all(xml_content, "//transition") %>%
  map_df(~{
    attrs <- xml_attrs(.x)
    tibble(
      name = attrs["name"],
      x = as.numeric(attrs["x"]),
      y = as.numeric(attrs["y"])
    )
  })

# Log position data
log_issue("INFO", sprintf("Extracted positions for %d places and %d transitions", 
                          nrow(place_positions), nrow(transition_positions)), 
          "Layout Analysis")

# Save position data to CSV files for external analysis
write_csv(place_positions,
          file.path(dirname(log_file), 
                    paste0(tools::file_path_sans_ext(basename(file_path)), 
                           "_place_positions.csv")))
write_csv(transition_positions,
          file.path(dirname(log_file), 
                    paste0(tools::file_path_sans_ext(basename(file_path)), 
                           "_transition_positions.csv")))

# Add basic layout validation
# Validate organism grouping 
for (model in bacterial_models) {
  abbr <- normalize(model$abbreviation)
  
  # Get places specific to this organism
  organism_places <- place_positions %>%
    filter(grepl(abbr, tolower(name)))
  
  if (nrow(organism_places) >= 2) {
    # Calculate spatial boundaries
    x_min <- min(organism_places$x)
    x_max <- max(organism_places$x)
    y_min <- min(organism_places$y)
    y_max <- max(organism_places$y)
    
    x_range <- x_max - x_min
    y_range <- y_max - y_min
    
    # Log spatial distribution
    log_issue("INFO", sprintf("%s places spatial distribution: x=[%.1f, %.1f], y=[%.1f, %.1f]", 
                              model$organism, x_min, x_max, y_min, y_max), "Layout Analysis")
    
    # Check if places are reasonably grouped (not too spread out)
    if (x_range > 20 || y_range > 20) {
      log_issue("WARNING", sprintf("%s places are spatially dispersed (x_range=%.1f, y_range=%.1f)", 
                                   model$organism, x_range, y_range), "Layout Analysis")
    } else {
      log_issue("INFO", sprintf("%s places are properly grouped together", model$organism), 
                "Layout Analysis")
    }
    
    # Check proximity between organism places and their transitions
    organism_transitions <- transition_positions %>%
      filter(grepl(abbr, tolower(name)))
    
    if (nrow(organism_transitions) > 0) {
      # Calculate "center of gravity" for organism places and transitions
      place_center_x <- mean(organism_places$x)
      place_center_y <- mean(organism_places$y)
      trans_center_x <- mean(organism_transitions$x)
      trans_center_y <- mean(organism_transitions$y)
      
      # Calculate distance between centers
      center_distance <- sqrt((place_center_x - trans_center_x)^2 + 
                                (place_center_y - trans_center_y)^2)
      
      log_issue("INFO", sprintf("%s transitions are %.1f units from places center", 
                                model$organism, center_distance), "Layout Analysis")
      
      if (center_distance > 15) {
        log_issue("WARNING", sprintf("%s transitions are distant from places (%.1f units)", 
                                     model$organism, center_distance), "Layout Analysis")
      }
    }
  } else {
    log_issue("WARNING", sprintf("Too few places found for %s (expected at least 2)", 
                                 model$organism), "Layout Analysis")
  }
}

# Check if shared metabolite places are positioned appropriately (centrally)
if (length(metabolite_places) > 0) {
  met_places <- place_positions %>%
    filter(tolower(name) %in% tolower(metabolite_places))
  
  if (nrow(met_places) > 0) {
    # Calculate center of all organism-specific places
    org_places <- place_positions %>%
      filter(!tolower(name) %in% tolower(metabolite_places))
    
    if (nrow(org_places) > 0) {
      org_center_x <- mean(org_places$x)
      org_center_y <- mean(org_places$y)
      
      # For each metabolite place, check if it's positioned relatively centrally
      for (i in 1:nrow(met_places)) {
        met_name <- met_places$name[i]
        met_x <- met_places$x[i]
        met_y <- met_places$y[i]
        
        # Calculate distance to organism center
        distance <- sqrt((met_x - org_center_x)^2 + (met_y - org_center_y)^2)
        
        log_issue("INFO", sprintf("Metabolite place '%s' is %.1f units from organism center", 
                                  met_name, distance), "Layout Analysis")
        
        if (distance > 25) {
          log_issue("WARNING", sprintf("Metabolite place '%s' is not centrally positioned", 
                                       met_name), "Layout Analysis")
        }
      }
    }
  }
}

# Calculate overall model dimensions and density
model_width <- max(place_positions$x, transition_positions$x) - 
  min(place_positions$x, transition_positions$x)
model_height <- max(place_positions$y, transition_positions$y) - 
  min(place_positions$y, transition_positions$y)
model_area <- model_width * model_height
element_density <- (nrow(place_positions) + nrow(transition_positions)) / model_area

log_issue("INFO", sprintf("Model dimensions: %.1f × %.1f units (area: %.1f sq. units)", 
                          model_width, model_height, model_area), "Layout Analysis")
log_issue("INFO", sprintf("Element density: %.4f elements per sq. unit", 
                          element_density), "Layout Analysis")

if (element_density < 0.005) {
  log_issue("WARNING", "Model elements are very sparse (consider more compact layout)", 
            "Layout Analysis")
} else if (element_density > 0.05) {
  log_issue("WARNING", "Model elements are very dense (may be difficult to read)", 
            "Layout Analysis")
} else {
  log_issue("INFO", "Model element density is appropriate", "Layout Analysis")
}