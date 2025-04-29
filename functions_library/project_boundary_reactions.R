project_boundary_reactions <- function(bacterial_models,
                                       metabolite_places,
                                       output_dir_projections) {
  # Build dataframe describing each model
  models_df <- tibble(model = bacterial_models) %>%
    mutate(
      abbr       = map_chr(model, ~ .x$abbreviation[2]),
      FBAmodel   = map_chr(model, ~ .x$FBAmodel),
      model_file = file.path("metabolic_networks_library", paste0(FBAmodel, ".mat")), # corrected
      meta_dir   = file.path("hypernodes", "minimal_doublet", paste0("metabolic_networks_", "minimal_doublet"), gsub(" ", "_", FBAmodel)) # generalized
    ) %>% filter(file.exists(model_file)) %>%
    mutate(
      metabolites = map(meta_dir, ~ read_csv(file.path(.x, "metabolites_metadata.csv"), show_col_types = FALSE)),
      reactions   = map(meta_dir, ~ read_csv(file.path(.x, "reactions_metadata.csv"),   show_col_types = FALSE))
    )
    
    boundary_df <- models_df %>%
    	dplyr::select(abbr, reactions) %>%
    	tidyr::unnest(reactions) %>%
    	dplyr::filter(type == "boundary") %>%
    	dplyr::select(abbr, reaction = abbreviation, lowbnd, uppbnd, equation)
  
  # Match projectable metabolites to boundary reactions
  shared_rxns_df <- models_df %>%
    unnest(metabolites) %>%
    filter(id %in% metabolite_places) %>%
    distinct(abbr, id) %>%
    left_join(boundary_df, by = "abbr", relationship = "many-to-many") %>%
    filter(str_detect(equation, str_c("\\b", id, "\\b"))) %>%
    distinct(abbr, reaction, lowbnd, uppbnd)
  
  # Organize projections
  shared_list <- shared_rxns_df %>%
    group_by(abbr) %>%
    summarize(rxns = list(reaction), .groups = "drop")
  
  all_models <- shared_list$abbr
  joint_rxns <- unique(shared_rxns_df$reaction)
  
  reaction_orgs <- shared_rxns_df %>%
    distinct(abbr, reaction) %>%
    group_by(reaction) %>%
    summarize(orgs = list(abbr), .groups = "drop")
  
  common_rxns <- reaction_orgs %>%
    filter(map_lgl(orgs, ~ setequal(.x, all_models))) %>%
    pull(reaction)
  
  exclusive_list <- shared_list %>%
    mutate(excl = map(rxns, ~ setdiff(.x, common_rxns)))
  
  # Output
  dir.create(output_dir_projections, showWarnings = FALSE, recursive = TRUE)
  
  write_lines(joint_rxns, file.path(output_dir_projections, "extracted_boundary_reactions.txt"))
  
  if (length(common_rxns) > 0)
    write_lines(common_rxns, file.path(output_dir_projections, "common_reactions.txt"))
  
  walk2(
    exclusive_list$excl, exclusive_list$abbr,
    ~ if (length(.x) > 0)
      write_lines(.x, file.path(output_dir_projections, sprintf("exclusive_reactions_%s.txt", .y)))
  )
  
bounds_tbl <- shared_rxns_df %>%
    dplyr::select(reaction, organism = abbr, lower_bound = lowbnd, upper_bound = uppbnd) %>%
    dplyr::distinct()
  
  write_csv(bounds_tbl, file.path(output_dir_projections, "reaction_bounds.csv"))
  
  bounds_tbl %>%
    group_by(organism) %>%
    group_walk(~ write_csv(.x, file.path(output_dir_projections, sprintf("per_bounds_%s.csv", .y$organism))))
  
# =============================
# NEW: Generate irreversible network version
# =============================
cat("\nGenerating irreversible network references...\n")

models_df %>%
  dplyr::select(abbr, reactions) %>%
  unnest(reactions) %>%
  group_by(abbr) %>%
  group_walk(~ {
    # Create irreversible reaction table
    irr_net <- .x %>%
      mutate(
        is_reversible = lowbnd < 0 & uppbnd > 0,
        forward_reaction = if_else(is_reversible, paste0(abbreviation, "_f"), abbreviation),
        reverse_reaction = if_else(is_reversible, paste0(abbreviation, "_r"), NA_character_),
        lb_f = if_else(is_reversible, 0, pmax(0, lowbnd)),
        ub_f = if_else(is_reversible, uppbnd, uppbnd),
        lb_r = if_else(is_reversible, 0, NA_real_),
        ub_r = if_else(is_reversible, -lowbnd, NA_real_)
      ) %>%
      dplyr::select(
        original_reaction = abbreviation,
        forward_reaction, lb_f, ub_f,
        reverse_reaction, lb_r, ub_r
      )
    
    # Reshape into long irreversible reaction table
    irr_long <- irr_net %>%
      pivot_longer(
        cols = c(forward_reaction, reverse_reaction),
        names_to = "direction",
        values_to = "reaction"
      ) %>%
      filter(!is.na(reaction)) %>%
      mutate(
        lb = if_else(direction == "forward_reaction", lb_f, lb_r),
        ub = if_else(direction == "forward_reaction", ub_f, ub_r)
      ) %>%
      dplyr::select(reaction, lb, ub)
    
    # Save to CSV
    irr_output_path <- file.path(output_dir_projections, sprintf("irreversible_reactions_%s.csv", unique(.y$abbr)))
    write_csv(irr_long, irr_output_path)
    
    cat(sprintf("  ✓ Irreversible network saved for %s: %s\n", unique(.y$abbr), irr_output_path))
  })

  # Return structured output
  invisible(list(
    shared_reactions    = set_names(shared_list$rxns, shared_list$abbr),
    common_reactions    = common_rxns,
    exclusive_reactions = set_names(exclusive_list$excl, exclusive_list$abbr),
    bounds              = bounds_tbl
  ))
}

