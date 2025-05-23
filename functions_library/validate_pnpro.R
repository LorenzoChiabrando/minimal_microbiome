
validate_pnpro <- function(pnpro2validate,
                           hypernode_root,
                           bacterial_models,
                           metabolite_places,
                           out_dir,
                           hypernode_name) {
  
  # Extract the 2nd abbreviation for each organism in order
  abbrs <- map_chr(bacterial_models, ~ .x$abbreviation[2])
  
  # ————————————————————————————————————————————
  # 0) Read & parse the PNPRO
  # ————————————————————————————————————————————
  xml <- read_xml(pnpro2validate)
  
  # grab all <transition> names/delays
  tnodes <- xml_find_all(xml, "//transition")
  t_names  <- xml_attr(tnodes, "name")
  t_delays <- xml_attr(tnodes, "delay")
  
  # build arc_df: only those arcs attached to either Call[…] or FBA[…]
  is_call <- str_detect(t_delays, "^Call\\[")
  is_fba  <- str_detect(t_delays, "^FBA\\[")
  special <- t_names[is_call | is_fba]
  
  arc_df <- xml_find_all(xml, "//arc") %>%
    map_df(~{
      a <- xml_attrs(.x)
      tibble(
        head         = a["head"],
        tail         = a["tail"],
        kind         = a["kind"],
        multiplicity = as.integer(coalesce(a["mult"], "1"))
      )
    }) %>%
    filter((kind=="INPUT"  & head %in% special) |
             (kind=="OUTPUT" & tail %in% special)) %>%
    transmute(
      transition   = if_else(kind=="INPUT", head, tail),
      direction    = kind,
      place        = if_else(kind=="INPUT", tail, head),
      multiplicity,
      command      = t_delays[match(transition, t_names)]
    )
  
  # ————————————————————————————————————————————
  # 1) Re-derive exactly which boundary reactions each species should have
  #    (same logic as project_boundary_reactions())
  # ————————————————————————————————————————————
  models_df <- tibble(model = bacterial_models) %>%
    mutate(
      abbr     = map_chr(model, ~ .x$abbreviation[2]),
      meta_dir = file.path(hypernode_root, "biounits", map_chr(model, ~ .x$FBAmodel))
    )
  
  # read metabolites & reactions metadata
  models_df <- models_df %>%
    mutate(
      metabolites = map(meta_dir, ~ read_csv(file.path(.x, "metabolites_metadata.csv"), show_col_types = FALSE)),
      reactions   = map(meta_dir, ~ read_csv(file.path(.x, "reactions_metadata.csv"), show_col_types = FALSE))
    )

  projectable_df <- models_df %>%
    tidyr::unnest(metabolites) %>%
    dplyr::filter(id %in% metabolite_places) %>%
    dplyr::select(abbr, met_id = id)
  
  # pull out boundary reactions
  boundary_df <- models_df %>%
    tidyr::unnest(reactions) %>%
    dplyr::filter(type=="boundary") %>%
    dplyr::select(abbr, reaction = abbreviation, equation)
  
  # match metabolite ↔ boundary reaction
shared_rxns_df <- projectable_df %>%
  left_join(boundary_df, by = "abbr", relationship = "many-to-many") %>%
  filter(str_detect(equation, paste0("\\b", met_id, "\\b"))) %>%
  distinct(abbr, met_id, reaction)
  
  # ————————————————————————————————————————————
  # 2) Parse all your FBA[…] commands out of the PN, keeping exact delay
  # ————————————————————————————————————————————
  fba_cmds <- tibble(
    transition = t_names,
    delay      = t_delays
  ) %>%
    dplyr::filter(str_detect(delay, "^FBA\\[")) %>%
    dplyr::mutate(
      # parse model, reaction, scale, countPlace, biomassPlace
      parts = str_match(
        delay,
        'FBA\\[ *"([^"]+)" *, *"([^"]+)" *, *([0-9\\.]+) *, *"([^"]+)" *, *"([^"]+)"'
      ),
      model    = parts[,2],
      reaction = parts[,3],
      abbr     = str_remove(parts[,5], "^n_")  # get the species abbr
    ) %>%
    dplyr::select(transition, delay, reaction, abbr)
  
  call_cmds <- tibble(
    transition = t_names,
    delay      = t_delays
  ) %>%
    dplyr::filter(str_detect(delay, "^Call\\[")) %>%
    dplyr::mutate(
      # parse: function name, full parameter expr, organism‐index
      parts = str_match(
        delay,
        'Call\\[\\s*"([^"]+)"\\s*,\\s*(.+)\\s*,\\s*([0-9]+)\\s*\\]'
      ),
      fun_name   = parts[,2],
      param_expr = parts[,3],
      org_index  = as.integer(parts[,4])
    ) %>%
    dplyr::select(
      transition,
      command    = delay,
      fun_name,
      param_expr,
      org_index
    )
  
  # append biomass reaction for every organism
  biomass_rxns <- tibble(
    abbr    = abbrs,
    met_id  = "biomass_e",       # no environmental metabolite
    reaction= "EX_biomass_e"
  )
  
  shared_rxns_df <- bind_rows(shared_rxns_df, biomass_rxns)
  
  # now this will include EX_biomass_e_in_<abbr> & EX_biomass_e_out_<abbr>
  repaired_fba_cmds <- shared_rxns_df %>%
    dplyr::mutate(
      model_file    = paste0(abbr, "_model.txt"),
      count_place   = paste0("n_",         abbr),
      biomass_place = paste0("biomass_e_", abbr),
      scaling       = 1L,
      is_biomass    = reaction == "EX_biomass_e"
    ) %>%
    crossing(direction = c("INPUT","OUTPUT")) %>%
    dplyr::mutate(
      transition = paste0(
        reaction,
        if_else(direction=="INPUT","_in_","_out_"),
        abbr
      ),
      place = if_else(
        reaction=="EX_biomass_e", 
        biomass_place,      # hook biomass reactions to biomass place
        met_id
      ),
      command = sprintf(
        'FBA["%s","%s",%d,"%s","%s"%s]',
        model_file,
        reaction,
        scaling,
        count_place,
        biomass_place,
        if_else(is_biomass, ', "true"', '')
      )
    ) %>%
    dplyr::select(transition, command, reaction, abbr, direction, place)

  # 2. Define column‐indices for each function
  func_cols <- c(Starvation = 0L, Duplication = 1L, Death = 2L)

  # 3. Map function names → transition‐prefixes
  prefix_map <- c(Starvation="Starv", Duplication="Dup", Death="Death")
  
  repaired_call_cmds <- purrr::imap_dfr(abbrs, function(abbr, idx) {
  org_index <- idx - 1L
  
  tibble(
    fun_name  = names(func_cols),
    col_index = unname(func_cols),
    org_index = org_index
  ) %>%
    dplyr::mutate(
      # always append the abbr, even for the first organism
      transition = paste0(prefix_map[fun_name], "_", abbr),
      command    = sprintf(
        'Call["%s", FromTable["organisms_parameters.csv", %d, %d], %d]',
        fun_name, org_index, col_index, org_index
      )
    ) %>%
    dplyr::select(transition, command, fun_name, org_index, col_index)
})

# ————————————————————————————————————————————
# 3) Build call_arcs from repaired_call_cmds
#    (hook each population‐process transition up to its count & biomass places)
# ————————————————————————————————————————————

call_arcs <- purrr::pmap_dfr(
  list(
    transition = repaired_call_cmds$transition,
    command    = repaired_call_cmds$command,
    fun_name   = repaired_call_cmds$fun_name,
    org_index  = repaired_call_cmds$org_index
  ),
  function(transition, command, fun_name, org_index) {
    abbr          <- abbrs[org_index + 1]
    count_place   <- paste0("n_",         abbr)
    biomass_place <- paste0("biomass_e_", abbr)
    
    switch(fun_name,
           # Starvation: INPUT from biomass only
           "Starvation" = tibble(
             transition, direction = "INPUT", place = biomass_place,
             multiplicity = 1L, command
           ),
           
           # Duplication: INPUT from count & biomass; OUTPUT back to both (count×2)
           "Duplication" = bind_rows(
             tibble(transition, direction = "INPUT",  place = count_place,   multiplicity = 1L, command),
             tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
             tibble(transition, direction = "OUTPUT", place = count_place,   multiplicity = 2L, command),
             tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 1L, command)
           ),
           
           # Death: INPUT from count & biomass; OUTPUT to biomass only
           "Death" = bind_rows(
             tibble(transition, direction = "INPUT",  place = count_place,   multiplicity = 1L, command),
             tibble(transition, direction = "INPUT",  place = biomass_place, multiplicity = 1L, command),
             tibble(transition, direction = "OUTPUT", place = biomass_place, multiplicity = 1L, command)
           )
    )
  }
)

# ————————————————————————————————————————————
# 4) Build the FBA arcs according to the connectivity patterns
# ————————————————————————————————————————————
fba_arcs <- purrr::pmap_dfr(
  list(
    transition = repaired_fba_cmds$transition,
    reaction   = repaired_fba_cmds$reaction,
    abbr       = repaired_fba_cmds$abbr,
    met_place  = repaired_fba_cmds$place,
    command    = repaired_fba_cmds$command
  ),
  function(transition, reaction, abbr, met_place, command) {
    count_place   <- paste0("n_",         abbr)
    biomass_place <- paste0("biomass_e_", abbr)
    
    if (reaction != "EX_biomass_e") {
      # Import vs Export
      if (endsWith(transition, paste0("_in_",  abbr))) {
        # IMPORT: consume metabolite + population/biomass, then re‐output population/biomass
        bind_rows(
          # INPUT arcs
          tibble(transition, direction="INPUT",  place=met_place,        multiplicity=1L, command),
          tibble(transition, direction="INPUT",  place=count_place,      multiplicity=1L, command),
          tibble(transition, direction="INPUT",  place=biomass_place,    multiplicity=1L, command),
          # OUTPUT arcs
          tibble(transition, direction="OUTPUT", place=count_place,      multiplicity=1L, command),
          tibble(transition, direction="OUTPUT", place=biomass_place,    multiplicity=1L, command)
        )
      } else {
        # EXPORT: consume population/biomass, then output metabolite + population/biomass
        bind_rows(
          # INPUT arcs
          tibble(transition, direction="INPUT",  place=count_place,      multiplicity=1L, command),
          tibble(transition, direction="INPUT",  place=biomass_place,    multiplicity=1L, command),
          # OUTPUT arcs
          tibble(transition, direction="OUTPUT", place=met_place,        multiplicity=1L, command),
          tibble(transition, direction="OUTPUT", place=count_place,      multiplicity=1L, command),
          tibble(transition, direction="OUTPUT", place=biomass_place,    multiplicity=1L, command)
        )
      }
    } else {
      # Biomass boundary
      if (endsWith(transition, paste0("_in_",  abbr))) {
        # BIOMASS IN: consume biomass
        tibble(transition, direction="INPUT",  place=biomass_place, multiplicity=1L, command)
      } else {
        # BIOMASS OUT: consume biomass, produce biomass×2
        bind_rows(
          tibble(transition, direction="INPUT",  place=biomass_place, multiplicity=1L, command),
          tibble(transition, direction="OUTPUT", place=biomass_place, multiplicity=2L, command)
        )
      }
    }
  }
)

# ————————————————————————————————————————————
# 5) Combine everything into arc_df_repaired
# ————————————————————————————————————————————
arc_df_repaired <- bind_rows(
  call_arcs,   # repaired Call[…] arcs
  fba_arcs     # repaired FBA[…] arcs
) %>%
  distinct(transition, direction, place, multiplicity, command)

# ————————————————————————————————————————————
# 6) Write out both raw and repaired arc tables
# ————————————————————————————————————————————

# raw (partial) arcs
write_csv(
  arc_df,
  file.path(out_dir, "raw_arcs.csv")
)

# full repaired arcs
write_csv(
  arc_df_repaired,
  file.path(out_dir, "repaired_arcs.csv")
)

# ————————————————————————————————————————————
# 7) Return both data frames
# ————————————————————————————————————————————
invisible(list(
  raw_arcs      = arc_df,
  repaired_arcs = arc_df_repaired
))

}
