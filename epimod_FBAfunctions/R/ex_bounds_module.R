# @Author: Mandringo02

# ============================================================
#  File: ex_bounds_module.R
#  Description: Un modulo in R che estrae le reazioni EX_ da
#               file di testo e le elabora per generare 
#               file CSV con vincoli calcolati dinamicamente.
# ============================================================

# ------------------------------------------------------------------
# 1) EXTRACT EX_ REACTIONS
#    Legge la prima riga di ciascun file in 'files' (un vettore 
#    di percorsi), divide per spazio, filtra le reazioni che 
#    iniziano con "EX_" e le scrive in 'output_file' con "_r" e "_f"
# ------------------------------------------------------------------
extract_ex_reactions <- function(
  files,
  output_file = file.path(getwd(), "all_ex_reactions.txt")
) {
  # Inizializza un vettore vuoto per memorizzare tutte le reazioni EX_
  all_ex_reactions <- character()
  
  # Itera su ciascun file in 'files'
  for (f in files) {
    if (!file.exists(f)) {
      warning(paste("File not found:", f, "- skipping."))
      next
    }
    # Legge solo la prima riga
    first_line <- readLines(f, n = 1, warn = FALSE)
    if (length(first_line) == 0) {
      # Se il file è vuoto o non ha una prima riga, salta
      next
    }
    
    # Divide la riga in base agli spazi
    reactions <- unlist(strsplit(first_line, " "))
    
    # Filtra solo le reazioni che iniziano con "EX_"
    ex_filtered <- reactions[grepl("^EX_", reactions)]
    
    # Accumula nel vettore principale
    all_ex_reactions <- c(all_ex_reactions, ex_filtered)
  }
  
  # Rimuove eventuali duplicati e ordina
  all_ex_reactions <- unique(all_ex_reactions)
  all_ex_reactions <- sort(all_ex_reactions)
  
  # Crea (o sovrascrive) il file di output
  con <- file(output_file, open = "w")
  
  # Per ciascuna reazione EX_, scrivi le versioni "_r" e "_f"
  for (rxn in all_ex_reactions) {
    writeLines(paste0(rxn, "_r"), con = con)
    writeLines(paste0(rxn, "_f"), con = con)
  }
  close(con)
  
  message("Le reazioni EX_ uniche sono state scritte in: ", output_file)
}


# ------------------------------------------------------------------
# 2) PROCESS EX_ REACTIONS TO CREATE CSV FILES
#    Questa funzione legge un file ('reaction_file') contenente righe 
#    del tipo "EX_something_r" o "EX_something_f". 
#    Vengono generati due CSV:
#      - FBA CSV:   contiene solo le reazioni le cui *base* 
#                   sono in 'fba_reactions'. Per queste, si imposta 
#                   SEMPRE la versione "_f" con l'upper bound 'fba_upper_bound'.
#      - nonFBA CSV: contiene le altre reazioni (non presenti in 'fba_reactions'),
#                    a cui vengono applicati i filtri (in base al parametro 
#                    reaction_versions) e il vincolo = not_shared_base_bound / bacteria_count.
# ------------------------------------------------------------------
process_ex_reactions <- function(
  reaction_file,
  bacteria_files,
  output_dir           = getwd(),
  fba_reactions        = character(),
  bacteria_counts      = c(1),
  not_shared_base_bound   = 1000,
  fba_upper_bound      = 0.015,
  reaction_versions    = c("both", "r", "f")  # parametro che controlla le direzioni per le non-FBA
) {
  # Forza uno dei tre valori consentiti per reaction_versions
  reaction_versions <- match.arg(reaction_versions)
  
  n_bact <- length(bacteria_files)
  if (n_bact < 1) {
    stop("Devi fornire almeno un modello batterico in 'bacteria_files'.")
  }
  if (!file.exists(reaction_file)) {
    stop("Il file delle reazioni non esiste: ", reaction_file)
  }
  if (length(bacteria_counts) != n_bact) {
    stop("La lunghezza di 'bacteria_counts' deve essere uguale a quella di 'bacteria_files'.")
  }
  if (length(not_shared_base_bound) != 1 && length(not_shared_base_bound) != n_bact) {
    stop("'not_shared_base_bound' deve essere un singolo numero o un vettore della stessa lunghezza dei batteri.")
  }
  if (length(fba_upper_bound) != 1 && length(fba_upper_bound) != n_bact) {
    stop("'fba_upper_bound' deve essere un singolo numero o un vettore della stessa lunghezza dei batteri.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_ub_fba_path    <- file.path(output_dir, "EX_upper_bounds_FBA.csv")
  output_ub_nonfba_path <- file.path(output_dir, "EX_upper_bounds_nonFBA.csv")
  
  f_fba    <- file(output_ub_fba_path, "w")
  f_nonfba <- file(output_ub_nonfba_path, "w")
  
  # Per le non-FBA, "replica" il valore base se necessario
  if (length(not_shared_base_bound) == 1) {
    base_values <- rep(not_shared_base_bound, n_bact)
  } else {
    base_values <- not_shared_base_bound
  }
  
  first_line <- paste(c("base_upper_bounds", base_values), collapse = ",")
  writeLines(first_line, con = f_nonfba)
  
  # Funzione per ottenere il nome base (senza _r o _f)
  get_base_name <- function(rxn) {
    sub("(_r|_f)$", "", rxn)
  }
  
  # Legge tutte le reazioni dal file
  reactions <- readLines(reaction_file, warn = FALSE)
  
  # Itera sulle reazioni
  for (reaction in reactions) {
    # Salta eventuali reazioni non desiderate (es. EX_biomass_e)
    if (grepl("EX_biomass_e", reaction)) next
    
    # Ottieni il nome base della reazione
    base_rxn <- get_base_name(reaction)
    
    # Se la reazione è una FBA (cioè, la sua base è in fba_reactions)
    if (base_rxn %in% fba_reactions) {
      # Per le FBA si vuole SEMPRE processare la versione _f
      if (!grepl("_f$", reaction)) next
      # Imposta l'upper bound per la FBA
      ub_values <- numeric(n_bact)
      for (i in seq_len(n_bact)) {
        ub_values[i] <- if (length(fba_upper_bound) == 1) {
          fba_upper_bound
        } else {
          fba_upper_bound[i]
        }
      }
      line_to_write <- paste(c(reaction, ub_values), collapse = ",")
      writeLines(line_to_write, con = f_fba)
    } else {
      # Per le non-FBA applica il filtro reaction_versions
      if (reaction_versions == "r" && grepl("_f$", reaction)) next
      if (reaction_versions == "f" && grepl("_r$", reaction)) next
      
      # Imposta l'upper bound per le reazioni non-FBA
      ub_values <- numeric(n_bact)
      for (i in seq_len(n_bact)) {
        ub_values[i] <- if (length(not_shared_base_bound) == 1) {
          not_shared_base_bound / bacteria_counts[i]
        } else {
          not_shared_base_bound[i] / bacteria_counts[i]
        }
      }
      line_to_write <- paste(c(reaction, ub_values), collapse = ",")
      writeLines(line_to_write, con = f_nonfba)
    }
  }
  
  # Chiude i file di output
  close(f_fba)
  close(f_nonfba)
  
  message("Elaborazione completata.\n  FBA file:    ", output_ub_fba_path,
          "\n  nonFBA file: ", output_ub_nonfba_path)
}


# ------------------------------------------------------------------
# 3) COMBINED FUNCTION
#    Una chiamata unica che esegue entrambe le operazioni:
#      - Estrazione delle reazioni EX_ dai file di testo
#      - Elaborazione delle reazioni per creare i file CSV
# ------------------------------------------------------------------
run_full_ex_bounds <- function(
  extraction_output  = "all_ex_reactions.txt",
  bacteria_files,
  output_dir         = getwd(),
  fba_reactions      = character(),
  bacteria_counts    = c(1),
  not_shared_base_bound = 1000,
  fba_upper_bound    = 0.015,
  reaction_versions  = c("both", "r", "f")
) {
  extract_ex_reactions(
    files       = bacteria_files,
    output_file = file.path(output_dir, extraction_output)
  )
  
  process_ex_reactions(
    reaction_file      = file.path(output_dir, extraction_output),
    bacteria_files     = bacteria_files,
    output_dir         = output_dir,
    fba_reactions      = fba_reactions,
    bacteria_counts    = bacteria_counts,
    not_shared_base_bound = not_shared_base_bound,
    fba_upper_bound    = fba_upper_bound,
    reaction_versions  = reaction_versions 
  )
}

