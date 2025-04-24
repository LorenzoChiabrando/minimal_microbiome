
# Generate a shortlist of plausible abbreviations for one model name
derive_abbrs <- function(model_name) {
  # drop any trailing “_model”
  core <- str_remove(model_name, "_model$")
  parts <- strsplit(core, "_")[[1]]
  abbrs <- character()
  
  # candidate 1: everything before the first underscore
  abbrs <- c(abbrs, parts[1])
  
  # candidate 2: first letters of each part
  abbrs <- c(abbrs, paste0(substr(parts, 1,1), collapse=""))
  
  # candidate 3: if there are exactly two parts, 1st letter+2nd letter
  if(length(parts)==2)
    abbrs <- c(abbrs, paste0(substr(parts[1],1,1), substr(parts[2],1,1)))
  
  unique(tolower(abbrs))
}

# Build your bacterial_models list
make_bacterial_models <- function(model_names,
                                  biomass_params,
                                  pop_params,
                                  initial_counts) {
  stopifnot(length(model_names) == length(biomass_params),
            length(model_names) == length(pop_params),
            length(model_names) == length(initial_counts))
  
  lapply(seq_along(model_names), function(i) {
    mn <- model_names[i]
    abbrs <- derive_abbrs(mn)
    
    list(
      FBAmodel     = mn,
      organism     = gsub("_", " ", mn),
      abbreviation = abbrs,            # we pass the vector of candidates
      txt_file     = mn,
      biomass      = biomass_params[[i]],
      bac_pop_p    = pop_params[[i]],
      initial_count= initial_counts[i]
    )
  })
}

# Write out your Bacteria_Parameters.csv from the new list
write_bac_params <- function(bacterial_models, path) {
  df <- do.call(rbind, lapply(bacterial_models, function(x) {
    as.data.frame(unlist(x$bac_pop_p), stringsAsFactors=FALSE)
  }))
  write.table(df, path, row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
}
