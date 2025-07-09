process_model <- function(m, hypernode_name) {
  cat("\n========== PROCESSING BIOUNIT ==========\n")
  
  unit <- m$unit
  abbr <- m$abbreviation[2]
  FBAmodel <- ifelse(is.null(m$FBAmodel), unit, m$FBAmodel)
  
  mat_file <- file.path(wd, "metabolic_networks_library", paste0(FBAmodel, ".mat"))
  input_dir <- file.path(wd, "hypernodes", hypernode_name, "biounits", FBAmodel)
  output_file <- file.path(input_dir, paste0(abbr, "_model.txt"))
  parent_dir <- dirname(input_dir)
  
  if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE)
  file.copy(mat_file, file.path(input_dir, paste0(FBAmodel, ".mat")), overwrite = TRUE)
  
  model_obj <- FBA4Greatmod.generation(fba_mat = file.path(input_dir, paste0(FBAmodel, ".mat")),
                                       input_dir = input_dir)
  model_obj <- setBiomassParameters(model_obj,
                                    bioMax  = m$biomass$max,
                                    bioMean = m$biomass$mean,
                                    bioMin  = m$biomass$min)
  
  cat("📝 Writing model with writeFBAfile...\n")
  before_files <- list.files(parent_dir, pattern = "\\.txt$", full.names = TRUE)
  
  writeFBAfile(model_obj,
               fba_fname = paste0(abbr, "_model.txt"),
               dest_dir  = input_dir)
  
  after_files <- list.files(parent_dir, pattern = "\\.txt$", full.names = TRUE)
  new_txt <- setdiff(after_files, before_files)
  
  if (length(new_txt) == 1) {
    ghost_path <- new_txt[1]
    cat("🔁 Moving ghost file to correct location:\n  FROM:", ghost_path, "\n  TO  :", output_file, "\n")
    file.rename(ghost_path, output_file)
  } else {
    cat("❌ Could not locate new .txt model output.\n")
    return(list(status = "error", message = "Model file not created", unit = unit, abbr = abbr))
  }
  
  cat("✅ Model saved to:", output_file, "\n")
  return(list(
    status = "success",
    message = sprintf("Processed %s (%s)", unit, abbr),
    unit = unit,
    abbr = abbr,
    model_file = output_file
  ))
}
