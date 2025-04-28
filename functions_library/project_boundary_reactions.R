
project_boundary_reactions <- function(bacterial_models,
                                       metabolite_places,
                                       output_dir_projections) {
  
  models_df <- tibble(model = bacterial_models) %>%
    mutate(
      abbr       = map_chr(model, ~ .x$abbreviation[2]),
      FBAmodel   = map_chr(model, ~ .x$FBAmodel),
      model_file = file.path("compiled_models", paste0(abbr, "_model.txt")),
      meta_dir   = file.path("input", FBAmodel)
    ) %>%
    filter(file.exists(model_file)) %>%
    mutate(
      metabolites = map(meta_dir, ~ read_csv(file.path(.x, "metabolites_metadata.csv"))),
      reactions   = map(meta_dir, ~ read_csv(file.path(.x, "reactions_metadata.csv")))
    )
  
  # 2. Find which metabolites are projectable
  projectable_df <- models_df %>%
    tidyr::unnest(metabolites) %>%
    dplyr::filter(id %in% metabolite_places) %>%
    dplyr::select(abbr, id) %>%
    dplyr::distinct()
  
  # 3. Pull out all boundary reactions
  boundary_df <- models_df %>%
    tidyr::unnest(reactions) %>%
    dplyr::filter(type == "boundary") %>%
    dplyr::select(abbr, reaction = abbreviation, lowbnd, uppbnd, equation)
  
  # 4. Match projectable metabolites to boundary reactions
  shared_rxns_df <- projectable_df %>%
    left_join(boundary_df,
              by           = "abbr",
              relationship = "many-to-many") %>%
    filter(str_detect(equation, str_c("\\b", id, "\\b"))) %>%
    distinct(abbr, reaction, lowbnd, uppbnd)
  
  # 5. Summarize: shared, joint, common, exclusive
  shared_list <- shared_rxns_df %>%
    group_by(abbr) %>%
    summarize(rxns = list(reaction), .groups="drop")
  
  all_models  <- shared_list$abbr
  joint_rxns  <- unique(shared_rxns_df$reaction)
  
  reaction_orgs <- shared_rxns_df %>%
    distinct(abbr, reaction) %>%
    group_by(reaction) %>%
    summarize(orgs = list(abbr), .groups="drop")
  
  common_rxns <- reaction_orgs %>%
    filter(map_lgl(orgs, ~ setequal(.x, all_models))) %>%
    pull(reaction)
  
  exclusive_list <- shared_list %>%
    mutate(excl = map(rxns, ~ setdiff(.x, common_rxns)))
  
  # 6. Write outputs into output_dir_projections
  dir.create(output_dir_projections, showWarnings = FALSE, recursive = TRUE)
  
  # 6a. joint & common & exclusive reaction files
  write_lines(joint_rxns, file.path(output_dir_projections, "extracted_boundary_reactions.txt"))
  if (length(common_rxns) > 0)
    write_lines(common_rxns, file.path(output_dir_projections, "common_reactions.txt"))
  
  walk2(
    exclusive_list$excl, exclusive_list$abbr,
    ~ if (length(.x) > 0)
      write_lines(.x, file.path(output_dir_projections,
                                sprintf("exclusive_reactions_%s.txt", .y)))
  )
  
  # 6b. bounds table & per‐model CSVs
  bounds_tbl <- shared_rxns_df %>%
    dplyr::select(reaction, organism = abbr,
           lower_bound = lowbnd, upper_bound = uppbnd) %>%
    dplyr::distinct()
  
  write_csv(bounds_tbl, file.path(output_dir_projections, "reaction_bounds.csv"))
  
  bounds_tbl %>%
    group_by(organism) %>%
    group_walk(~ write_csv(.x,
                           file.path(output_dir_projections,
                                     sprintf("pper_bounds_%s.csv", .y$organism))))
  
  # Return a summary list invisibly
  invisible(list(
    shared_reactions   = set_names(shared_list$rxns,   shared_list$abbr),
    common_reactions   = common_rxns,
    exclusive_reactions= set_names(exclusive_list$excl, exclusive_list$abbr),
    bounds             = bounds_tbl
  ))
}
