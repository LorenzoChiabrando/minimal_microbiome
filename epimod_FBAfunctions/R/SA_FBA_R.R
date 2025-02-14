
SA_FBA = function(result_dir,
                  model_file,
                  all_react_file,
                  param_values_file,
                  bounds_set_file,
                  fbasol_t_file,
                  fbasol_file,
                  Y_file,
                  Y2_file,
                  indices_file,
                  Y2_reaction,
                  # Define settings
                  # N = 2^13
                  N,
                  cores,
                  b1,
                  b2) {
 
  # Global Sensitivity Analysis using R
  #
  # We will investigate the sensitivty of the reaction objective ("biomass205") 
  # flux solution by varying the D parameters within the bounds (-10.0, 0.0) 
  # defined by the problem shown below.
 
  # FBA model reactions annotation
  all_react = readRDS(all_react_file)
  ex_react = all_react[[2]]
  all_react = all_react[[1]]
  
  # Define settings
  chunks <<- cores*500
  D = length(ex_react$react.id)
  params = paste("P", 1:D, sep = ".")
 
  # Create sample matrix using Sobol' Quasi Random Numbers.
  # The parameter values are sampled using a uniform distribution (its quantile function)
  mat = as.data.frame(apply(sobol_matrices(N = N, 
                                           params = params, 
                                           order = "first"), 2, 
                            function(x) qunif(x, b1, b2)))
  
  mat_list = split(mat, rep(1:(nrow(mat) %/% chunks + 1), 
                            each = chunks, length.out = nrow(mat)))
  
  for( i in 1:length(mat_list) ) {
    
    print(paste0('chunk: ', i, ' out of ', length(mat_list)))
    saveRDS(mat_list[[i]], param_values_file)
    
    ## FBA runcode
    ParallelFBA_sens(model_file = model_file,
                     param_values_file = param_values_file,
                     model_metadata_file = all_react_file,
                     result_dir = result_dir,
                     input_dir = paste0(wd, "/inst/input"),
                     reactIndx4secondSol = 
                       all_react$ReactionPos[which(
                         all_react$React_ID == Y2_reaction)],
                     cores = cores,
                     wd = wd)
    
    file.rename(bounds_set_file, paste0(result_dir, "/bounds_set", "_", i, ".rds"))
    file.rename(fbasol_file, paste0(result_dir, "/fbasol", "_", i, ".rds"))
    
    Y = readRDS(Y_file)
    Y2 = readRDS(Y2_file)
    
    if(i == 1) {
      Y_t = Y
      Y2_t = Y2
    } else {
      Y_t = c(Y_t, Y)
      Y2_t = c(Y2_t, Y2)
    }
    
    gc()
    
  }
  
  # saving data
  saveRDS(mat, param_values_file)
  saveRDS(Y_t, Y_file)
  saveRDS(Y2_t, Y2_file)
  
  param_values = readRDS(param_values_file)
  Y = readRDS(Y_file)
  Y2 = readRDS(Y2_file)
  
  gc()
  
  # Bootstraping is a statistical technique used to estimate the sampling distribution 
  # of an estimator by resampling the original data with replacement. 
  # It is used to construct bootstrap confidence intervals on sensitivity indices 
  # computed by polynomial chaos expansion
  
  # When parallel = "multicore" is used, each worker process inherits
  # the environment of the current session, including the workspace and the loaded namespaces and
  # attached packages (but not the random number seed: see below).
  
  # Compute and bootstrap Sobol' indices
  sensitivity_Y = sobol_indices(Y = Y, N = N, params = params, type = "norm",
                                boot = TRUE, R = 1000, parallel = "multicore")
  
  sensitivity_Y2 = sobol_indices(Y = Y2, N = N, params = params, 
                                 boot = TRUE, R = 1000, parallel = "multicore")
  
  sens.ind = rbind(
    "Y" = rbind(first = data.frame(parameter = params,
                                   value = dplyr::filter(
                                     sensitivity_Y$results, sensitivity == "Si")$original,
                                   low_ci = dplyr::filter(
                                     sensitivity_Y$results, sensitivity == "Si")$low.ci,
                                   high_ci = dplyr::filter(
                                     sensitivity_Y$results, sensitivity == "Si")$high.ci),
                total = data.frame(parameter = params,
                                   value = dplyr::filter(
                                     sensitivity_Y$results, sensitivity == "Ti")$original,
                                   low_ci = dplyr::filter(
                                     sensitivity_Y$results, sensitivity == "Ti")$low.ci,
                                   high_ci = dplyr::filter(
                                     sensitivity_Y$results, sensitivity == "Ti")$high.ci)), 
    "Y2" = rbind(first = data.frame(parameter = params,
                                    value = dplyr::filter(
                                      sensitivity_Y2$results, sensitivity == "Si")$original,
                                    low_ci = dplyr::filter(
                                      sensitivity_Y2$results, sensitivity == "Si")$low.ci,
                                    high_ci = dplyr::filter(
                                      sensitivity_Y2$results, sensitivity == "Si")$high.ci),
                 total = data.frame(parameter = params,
                                    value = dplyr::filter(
                                      sensitivity_Y2$results, sensitivity == "Ti")$original,
                                    low_ci = dplyr::filter(
                                      sensitivity_Y2$results, sensitivity == "Ti")$low.ci,
                                    high_ci = dplyr::filter(
                                      sensitivity_Y2$results, sensitivity == "Ti")$high.ci)))
  
  df_ind = sens.ind %>%
    rownames_to_column("sol") %>%
    separate(sol, into = c("sol", "order", "param"), sep="\\.")
  
  saveRDS(df_ind, indices_file)
  
  gc()
  
}

##### other functions

GenerateReactionEquations = function(model,
                                      makeClosedNetwork,
                                      entrydelim = ", ",
                                      extMetFlag = "b") {
  
  # required data structures
  equat  <- vector(mode = "character", length = dim(model@S)[2])
  compat <- vector(mode = "character", length = dim(model@S)[2])
  
  react_rev = c()
  
  b = cbind(lb = model@lowbnd, ub = model@uppbnd)
  
  for(i in 1:nrow(b)) {
    react_rev[i] = !any(b[i, ] == 0)
  }
  
  revers = ifelse(react_rev, "Reversible", "Irreversible")
  arrow = ifelse(react_rev, "<==>", "-->")
  
  # remove compartment flag if existing
  metab = sub("\\[\\w+\\]$", "", model@met_id)
  
  met.id = model@met_id
  
  if( sum(unique(stringr::str_sub(met.id, -1)) == "]") >= 1 ) {
    met.id = gsub("\\[",replacement = "_", met.id)
    met.id = gsub("\\]",replacement = "", met.id)
  }
  
  metcp = stringr::str_sub(met.id, -1)
  mod_compart = unique(stringr::str_sub(met.id, -1))
  
  met_comp = c()
  
  for (i in 1:length(met.id)) {
    comp = stringr::str_sub(met.id[i], -1)
    met_comp[i] = which(comp == mod_compart)
  }
  
  for (j in 1:dim(model@S)[2]) {
    
    column = model@S[, j]
    # row indices
    constr_ind = which(column != 0)
    # stoichiometric coefficients
    stcoef = column[constr_ind]
    
    # check if reaction is empty
    if (length(constr_ind) > 0) {
      
      comp = unique(mod_compart[met_comp[constr_ind]])
      
      # reaction involves more than one compartment
      if (length(comp) > 1) {
        if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
          metr = paste(model@met_id[constr_ind], 
                       "_", mod_compart[met_comp[constr_ind]], sep = "")
        } 
        else {
          metr = met.id[constr_ind]
        }
        compat[j] = paste(comp, collapse = entrydelim)
        compflag  = ""
      } 
      else {
        
        # Check if the current reaction is an exchange reaction.
        # In order to build a closed network, we need to add a 'boundary'
        # metabolite [b].
        
        if ((isTRUE(makeClosedNetwork)) && (length(constr_ind) == 1)) {
          if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
            metIN = paste(met.id[constr_ind], "_",
                          mod_compart[met_comp[constr_ind[1]]], sep = "")
          }
          else {
            metIN = met.id[constr_ind]
          }
          metr = c(metIN, paste(metab[constr_ind], "_", extMetFlag, sep = ""))
          constr_ind = c(constr_ind, constr_ind)
          stcoef =  c(stcoef, (stcoef * -1))
          compat[j] = comp
          compflag = ""
        }
        else {
          metr = metab[constr_ind]
          compat[j] <- comp
          # if yes, the metabolite id does not contain the metabolite compartment
          if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
            compflag  <- paste0("_", mod_compart[met_comp[constr_ind[1]]])
          }
          else {
            compflag  <- metcp[constr_ind[1]]
          }
          metr = paste0(metab[constr_ind], "_", comp)
          metr = strex::str_singleize(metr, paste0( "_", comp))
        }
      }
      
      educt <- vector(mode = "list")
      product <- vector(mode = "list")
      
      for (i in seq(along = constr_ind)) {
        if (stcoef[i] > 0) {
          stoich <- ifelse(stcoef[i] != 1,
                           paste("(", stcoef[i], ") ", sep = ""),
                           "")
          product[metr[i]] <- paste(stoich, metr[i], sep = "")
        }
        else {
          stoich <- ifelse(stcoef[i] != -1,
                           paste("(", (stcoef[i] * -1), ") ", sep = ""),
                           "")
          educt[metr[i]] <- paste(stoich, metr[i], sep = "")
        }
      }
      
      equattmp = paste(paste(educt, collapse = " + "),
                       arrow[j],
                       paste(product, collapse = " + "))
      
      equat[j] = equattmp
      
    }
  }
  
  return(list(equat  = equat,
              compat = compat,
              revers = revers,
              metab  = metab))
  
}

######

ExtractEx <- function(model, bigg.path) {
  
  StoichM = model@S
  
  # identifying columns with only one entry
  oneEntry = Matrix::colSums(StoichM != 0) == 1
  
  if (sum(oneEntry) > 0) {
    
    # exchange reactions can be with a -1 or 1
    ExReact = (Matrix::colSums(StoichM[, oneEntry, drop = FALSE] == 1) == 1 | 
                 Matrix::colSums(StoichM[, oneEntry, drop = FALSE] == -1) == 1)
    
    # finding Ex_reaction's IDs
    ex = c(1:dim(StoichM)[2])[oneEntry[ExReact]]
    
    # extracting metabolites involved in exchange reactions
    Met = which(Matrix::rowSums(
      abs(StoichM[, ex, drop = FALSE])) > 0)[StoichM[which(
        Matrix::rowSums(abs(StoichM[, ex, drop = FALSE])) > 0), 
        ex, drop = FALSE]+1]
    
    BIGGdata = read.delim2(bigg.path)
    
    reactions = data.frame(index = ex,
                           react.id = model@react_id[ex],
                           lb = model@lowbnd[ex], 
                           ub = model@uppbnd[ex],
                           react.eq = BIGGdata[ex, ]$equation
    )
    
  } else {
    
    warning("model without exchage ractions")
    reactions <- NULL
    
  }
  
  reactions$react.id = gsub("\\(", replacement = "_", reactions$react.id)
  reactions$react.id = gsub("\\[", replacement = "_", reactions$react.id)
  reactions$react.id = gsub("\\)", replacement = "", reactions$react.id)
  reactions$react.id = gsub("\\]", replacement = "", reactions$react.id)
  
  reactions$react.id = gsub("_c_", replacement = "_c", reactions$react.id)
  
  return(reactions)
  
}

######

FBAmodel.metadata <- function(model, model.name, wd, bigg.path) {
  
  br = ExtractEx(model, bigg.path = bigg.path)
  
  data = rbind(data.frame(value = br$lb,
                          dir = rep("lowbnd", length(br$lb)),
                          react.id = br$react.id,
                          react.index = br$index),
               data.frame(value = br$ub,
                          dir = rep("uppbnd", length(br$ub)),
                          react.id = br$react.id,
                          react.index = br$index))
  
  all.react = data.frame(React_ID = model@react_id,
                         ReactionType = rep(".", dim(model@S)[2]),
                         React_Lb = model@lowbnd,
                         React_Ub = model@uppbnd,
                         ReactionDistrict = rep(".", dim(model@S)[2]),
                         ReactionPos = 1:dim(model@S)[2])
  
  all.react$ReactionDistrict[data$react.index] = "boundary"
  all.react$ReactionDistrict[which(all.react$ReactionDistrict != "boundary")] = "core"
  
  a = dplyr::filter(all.react, grepl("EX_", React_ID))
  a$ReactionType = rep("Exchange", length(a$ReactionType))
  
  b = dplyr::filter(all.react, grepl("DM_", React_ID))
  b$ReactionType = rep("Demand/Sink", length(b$ReactionType))
  
  c = dplyr::filter(all.react, grepl("sink_", React_ID))
  c$ReactionType = rep("Demand/Sink", length(c$ReactionType))
  
  d = dplyr::filter(all.react, grepl("trans", React_ID))
  d$ReactionType = rep("Transcription", length(d$ReactionType))
  
  e = dplyr::filter(all.react, grepl("repl", React_ID))
  e$ReactionType = rep("Replication", length(e$ReactionType))
  
  f = dplyr::filter(all.react, grepl("pbios", React_ID))
  f$ReactionType = rep("Biosynthesis", length(f$ReactionType))
  
  g = dplyr::filter(all.react, !grepl("BIOMASS|biomass|rep|tra|EX_|DM_|sink_|pbios", React_ID))
  g$ReactionType = rep("Internals/Transporters", length(g$ReactionType))
  
  h = all.react[which(model@obj_coef == 1), ]
  h$ReactionType = rep("Objective", length(h$ReactionType))
  
  all.react = rbind(a, b, c, d, e, f, g, h)
  
  all.react$React_ID = gsub("\\(", replacement = "_", all.react$React_ID)
  all.react$React_ID = gsub("\\[", replacement = "_", all.react$React_ID)
  all.react$React_ID = gsub("\\)", replacement = "", all.react$React_ID)
  all.react$React_ID = gsub("\\]", replacement = "", all.react$React_ID)
  
  all.react$React_ID = gsub("_c_", replacement = "_c", all.react$React_ID)
  
  return(list(all.react, br))
  
}

######

findReactEq <- function(model, reagent, metadata) {
  
  # Returns a list of equations (formulas) of reactions in which the 
  # metabolite (reagent) is involved
  #
  # USAGE:
  #
  #   findReactEq(model, reagent)
  #
  # INPUTS:
  #    model:             Model structure
  #    reagent:           Metabolite
  #    metadata:          Reaction/Metabolite and equations data
  #
  # OUTPUTS:
  #     List of reactions + Reaction formulas corresponding and indexing
  
  # reactions nomenclature in accordance with VMH
  ReactionsNames = unlist(model@react_id)
  # reactant's nomenclature in accordance with VMH
  ReagentsNames = unlist(model@met_id)
  
  S = as.matrix(model@S)
  ReagentsIndex = which(ReagentsNames == reagent)
  ReactionIndex = which(S[ReagentsIndex, ] != 0)
  reactList = ReactionsNames[ReactionIndex]
  
  ReactEq = metadata[which(metadata$abbreviation %in% reactList), c(1, 2, 3, 4)]
  # sybil::printReaction(model, react = reactList, printOut = TRUE)
  
  return(ReactEq)
}

######

metadataEQ = function(model, model.name, 
                      prefix, suffix, extMetFlag = "b",
                      fielddelim = "\t", entrydelim = ", ",
                      makeClosedNetwork = FALSE,
                      onlyReactionList,
                      minimalSet,
                      fpath = SYBIL_SETTINGS("PATH_TO_MODEL")) {
  
  # validate model structure before writing
  validObject(model)
  
  # filenames
  if (missing(prefix)) {
    prefix = gsub("\\s+", "_", model.name)
  }
  
  if (missing(suffix)) {
    suffix = switch(fielddelim,
                    "\t" = { "tsv" },
                    ";"  = { "csv" },
                    ","  = { "csv" },
                    "|"  = { "dsv" },
                    { "dsv" }
    )
  }
  
  fnameR = paste(paste(prefix, "react", sep = "_"), suffix, sep = ".")
  fnameM = paste(paste(prefix, "met",   sep = "_"), suffix, sep = ".")
  fnameD = paste(paste(prefix, "desc",  sep = "_"), suffix, sep = ".")
  
  # path to output file
  tsvfileR = file.path(fpath, fnameR)
  tsvfileM = file.path(fpath, fnameM)
  tsvfileD = file.path(fpath, fnameD)
  
  #--------------------------------------------------------------------------#
  # reactions list
  #--------------------------------------------------------------------------#
  
  # create reaction strings
  rstr = GenerateReactionEquations(model,
                                   makeClosedNetwork,
                                   entrydelim,
                                   extMetFlag)
  
  met.id = model@met_id
  
  if( sum(unique(stringr::str_sub(met.id, -1)) == "]") > 1) {
    met.id = gsub("\\[",replacement = "_", met.id)
    met.id = gsub("\\]",replacement = "", met.id)
  }
  
  mod_compart = unique(stringr::str_sub(met.id, -1))
  
  met_comp = c()
  
  for (i in 1:length(met.id)) {
    comp = stringr::str_sub(met.id[i], -1)
    met_comp[i] = which(comp == mod_compart)
  }
  
  if (isTRUE(onlyReactionList)) {
    write.table(x = data.frame(equation = rstr$equat),
                row.names = FALSE, 
                file = tsvfileR, 
                sep = fielddelim)
  } else if (isTRUE(minimalSet)) {
    write.table(x = data.frame(abbreviation = model@react_id,
                               equation = rstr$equat,
                               lowbnd = model@lowbnd,
                               uppbnd =  model@uppbnd,
                               obj_coef = model@obj_coef),
                row.names = FALSE, 
                file = tsvfileR, 
                sep = fielddelim)
  } 
  else {
    
    write.table(x = data.frame(
      abbreviation = model@react_id,
      equation     = rstr$equat,
      reversible   = rstr$revers,
      compartment  = rstr$compat,
      lowbnd       = model@lowbnd,
      uppbnd       = model@uppbnd,
      obj_coef     = model@obj_coef),
      row.names = FALSE, 
      file = tsvfileR, 
      sep = fielddelim)
  }
  
  #--------------------------------------------------------------------------#
  # metabolites list
  #--------------------------------------------------------------------------#
  
  if ( (!isTRUE(onlyReactionList)) && (!isTRUE(minimalSet)) ) {
    
    metunq <- sort(unique(rstr$metab))
    mpos <- lapply(metunq, function(x) rstr$metab %in% x)
    
    # metabolite names
    metNames <- lapply(mpos, function(x) unique(model@met_id[x]))
    metNames <- unlist(lapply(metNames, paste, collapse = entrydelim))
    
    # metabolite compartments
    metCompart = lapply(mpos, function(x) mod_compart[met_comp[x]])
    metCompart = unlist(lapply(metCompart, paste, collapse = entrydelim))
    
    write.table(x = data.frame(
      abbreviation = metunq,
      name         = metNames,
      compartment  = metCompart),
      row.names = FALSE, file = tsvfileM, sep = fielddelim)
    
  }
  
  #--------------------------------------------------------------------------#
  # model description
  #--------------------------------------------------------------------------#
  
  if ((!isTRUE(onlyReactionList)) && (!isTRUE(minimalSet))) {
    
    # get id's of metabolites in different compartments
    # (one per compartment)
    metDiffComp = match(mod_compart, mod_compart[met_comp])
    metAbbrevComp = character(length(metDiffComp))
    
    # get the compartment abbreviations
    metALl = grepl("^.+(\\[\\w+\\])$", model@met_id[metDiffComp])
    metAbbrevComp[metALl] = sub("^.+(\\[\\w+\\])$", "\\1", model@met_id[metDiffComp[metALl]])
    metAbbrevComp[!metALl] = mod_compart[!metALl]
    
    # generate output format
    ma = paste(metAbbrevComp, collapse = entrydelim)
    mc = paste(mod_compart, collapse = entrydelim)
    
    write.table(x = data.frame(
      name         = model.name,
      id           = "FBAmodel",
      description  = model.name,
      compartment  = mc,
      abbreviation = ma,
      Nmetabolites = dim(model@S)[1],
      Nreactions   = dim(model@S)[2]),
      row.names = FALSE, file = tsvfileD, sep = fielddelim)
    
  }
  
  #--------------------------------------------------------------------------#
  # end
  #--------------------------------------------------------------------------#
  
  return(TRUE)
  
}

######

# loading gene expression data:
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

######

distr_stats = function(result_dir,
                       fbasol_file,
                       file_sd,
                       thrsd,
                       linethrsd) {
  
  fbasol = readRDS(fbasol_file)
  
  # Most variable fluxes
  fbasol = as.data.frame(fbasol)
  
  # Define function to calculate statistics for each Reaction
  get_stats = function(r) {
    
    d = dplyr::filter(fbasol, Reaction == r)
    Flux = tidyr::gather(d, key, value, -c(Reaction))$value
    Flux = as.numeric(Flux, digits = 12)
    sd = sd(Flux)
    E = mean(Flux)
    median = median(Flux)
    variance = var(Flux)
    c(r, as.double(sd), E, median, variance)
    
  }
  
  # use multicore, set to the number of our cores
  registerDoParallel((detectCores()))
  
  data = foreach(r = unique(fbasol$Reaction), 
                 .combine = rbind) %dopar% {
                   get_stats(r)
                 }
  
  colnames(data) = c("react.id", "sd", "E", "Median", "Variance")
  row.names(data) = NULL
  data = as.data.frame(data)
  data$sd = as.double(data$sd); data$E = as.double(data$E)
  data$Median = as.double(data$Median); data$Variance = as.double(data$Variance)
  
  names(all_react)[1] = names(data)[1]
  
  subdata = merge(x = data, y = all_react, by = "react.id", all.x = TRUE)
  
  mode = unique(as.double(subdata$E))[which.max(
    tabulate(match(as.double(subdata$E), unique(as.double(subdata$E)))))]
  
  # Plotting fluxes variance among configurations and differences compare to reference
  title.size = 8
  subtitle.size = 6
  text.size = 8 
  axistext.size = 6
  axistitle.size = 6
  legendsize = 0.25
  point.size = 0.5
  
  subdata.sd = dplyr::filter(subdata, as.double(subdata$sd) > thrsd)
  
  psdAll = ggplot(subdata, aes(x = as.double(sd), fill = "")) + 
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste0("SD Distr. | ", "number of reaction = ",
                        length(subdata$react.id)),
         subtitle = paste0("dotted line : sd = ", linethrsd, "\n", 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) + 
    geom_vline(xintercept = linethrsd, linetype = "dotted", color = "darkred", lwd = 0.5) +
    geom_histogram(aes(y = ..count..), alpha = 0.5, 
                   position = "identity", binwidth = 9,
                   col = "#651e3e", fill = "#651e3e")
  
  psd = ggplot(subdata.sd, aes(x = as.double(sd), fill = "")) + 
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste0("SD Distr. | ", "number of reaction = ",
                        length(subdata.sd$react.id)),
         subtitle = paste0("sd < ", thrsd, " was filtered out", "\n", 
                           "dotted line : sd = ", linethrsd, "\n" , 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) + 
    geom_vline(xintercept = linethrsd, linetype = "dotted", color = "darkred", lwd = 0.5) +
    geom_histogram(aes(y = ..count..), alpha = 0.5, 
                   position = "identity", binwidth = 9,
                   col = "#ceb5b7", fill = "#ceb5b7")
  
  psd_wrapAll = ggplot(subdata, aes(as.double(sd))) +
    geom_histogram(aes(y = (..count..), fill = ReactionDistrict)) +
    facet_wrap("ReactionDistrict", scales = "free") +
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste("SD Distr. | ", "number of reaction = ", 
                       length(subdata$react.id), sep = ""),
         subtitle = paste0("Number of ractions (core) = ", 
                           sum(subdata$ReactionDistrict == "core"), " \n",
                           "Number of ractions (boundary) = ", 
                           sum(subdata$ReactionDistrict == "boundary"), " \n", 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = text.size, color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          legend.key.size = unit(legendsize, 'cm'),
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold"))
  
  psd_wrap = ggplot(subdata.sd, aes(as.double(sd))) +
    geom_histogram(aes(y = (..count..), fill = ReactionDistrict)) +
    facet_wrap("ReactionDistrict", scales = "free") +
    labs(x = "sd (mmol/gDW*h)", y = "Counts",
         title = paste0("SD Distr. | ", "number of reaction = ", 
                        length(subdata.sd$react.id)),
         subtitle = paste0("Number of ractions (core) = ", 
                           sum(subdata.sd$ReactionDistrict == "core"), " \n",
                           "Number of ractions (boundary) = ", 
                           sum(subdata.sd$ReactionDistrict == "boundary"), " \n", 
                           paste0("config = ", (dim(param_values)[1] + 1)))) +
    theme(text = element_text(size = text.size, color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          legend.key.size = unit(legendsize, 'cm'),
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold"))
  
  ggsave( (psdAll|psd_wrapAll)/(psd|psd_wrap) , 
         file = paste0(result_dir, "/sd.pdf"), 
         width = 7, height = 4)
  
  subdata.sd_boundary = dplyr::filter(subdata.sd, subdata.sd$ReactionDistrict == "boundary")
  
  if( identical(which(all_react$React_ID == "sink_pheme_c"), integer(0)) ) {
    
    targetrnx = c("EX_trp_L_e", "EX_pro_L_e", 
                  "EX_val_L_e", "EX_ile_L_e", 
                  "EX_cys_L_e", "EX_leu_L_e",
                  "EX_lys_L_e",  "EX_biomass_e")
    
  } else {
    
    targetrnx = c("EX_trp_L_e", "EX_pro_L_e", 
                  "EX_val_L_e", "EX_ile_L_e", 
                  "EX_cys_L_e", "EX_leu_L_e", 
                  "EX_biomass_e", "sink_pheme_c")
    
  }
  
  fbasol = readRDS(paste0(result_dir, "/fbasol.rds"))
  fbasol_t = readRDS(paste0(result_dir, "/fbasol_template.rds"))
  
  subflux = data.frame(Reaction = rep(fbasol[, 1], (dim(fbasol)[2] - 1)),
                       Flux = gather(as.data.frame(fbasol), key, 
                                     value, -c(Reaction))$value,
                       config = str_remove(gather(as.data.frame(fbasol), key, 
                                                  value, -c(Reaction))$key, "config."))
  
  subflux = cbind(subflux, 
                  data.frame(lapply(param_values[, which(br$react.id %in% targetrnx)], 
                                    rep, each = dim(fbasol)[1])))
  
  subflux$Flux = as.double(subflux$Flux)
  
  df_plot = data.frame(ass.param = paste0("P.", which(br$react.id %in% targetrnx)),
                       react.index = which(br$react.id %in% targetrnx),
                       Reaction = br$react.id[which(br$react.id %in% targetrnx)],
                       col.data = seq(colnames(subflux[4:dim(subflux)[2]])) + 3,
                       col = brewer.pal(length(targetrnx),"Spectral"))
  
  subflux_t = cbind(data.frame(fbasol_t$abbreviation,
                               fbasol_t$flux,
                               rep(0, dim(fbasol)[1])),
                    data.frame(lapply(fbasol_t$lowbnd[df_plot$react.index],
                                      rep, each = dim(fbasol)[1])))
  
  colnames(subflux_t) = colnames(subflux)
  
  subflux = rbind(subflux, subflux_t)
  subflux = dplyr::filter(subflux, config != "0")[, c(1, 2, 3)]
  
  subflux.boun = 
    dplyr::filter(subflux, subflux$Reaction %in% dplyr::filter(all_react, ReactionDistrict == "boundary")$react.id)
  
  ssd = dplyr::filter(subflux.boun, subflux.boun$Reaction %in% subdata.sd_boundary$react.id)
  
  
  l = "#cc002d"; m = "#fffdbf"; h = "#112764"; c = "#653592"; pn = "darkred"
  
  title.size = 5; subtitle.size = 5; text.size = 5; axistext.size = 4; 
  axistitle.size = 5; legendsize = 0.25; point.size = 0.5
  
  react.subset = ggplot(ssd, aes(x = Flux, fill = ordered(Reaction, levels = unique(Reaction)))) +
    geom_density(aes(y = ..count..), alpha = 0.5, position = "identity", col = c, fill = c) +
    ggplot2::labs(x = "Flux (mmol/gDW*h)", y = "Count",
                  title = paste0("Reactions with sd > ", thrsd),
                  subtitle = paste0("Number of reactions: ", length(unique(ssd$Reaction)), " \n", 
                                    "Number of configurations: ", length(unique(ssd$config)))) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Condition block")) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = 14, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = 12, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) +
    facet_wrap(~ Reaction, scales = "free")
  
  ggsave(react.subset, file = paste0(result_dir, "/react_subsetsd.pdf"),
         width = sqrt(length(unique(ssd$Reaction)))*1.75, 
         height = sqrt(length(unique(ssd$Reaction))*1.75))
  
}

######

plotting_indeces = function(p, result_dir) {
  
  p_f = ggplot(dplyr::filter(p, order == "first"),
               aes(x = reorder(react.id, -value), y = value, fill = react.id)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1, position = position_dodge(0.5)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(p$react.id)))) +
    labs(x = "", y = "First-order sensitivity index") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, family = "Helvetica"),
          axis.text.y = element_text(face = "plain", size = 10, family = "Helvetica"),
          axis.title.x = element_text(face = "plain", size = 10, family = "Helvetica"),
          axis.title.y = element_text(face = "bold", size = 9, family = "Helvetica"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(size = 0.35)) 
  
  p_t = ggplot(dplyr::filter(p, order == "total"), 
               aes(x = reorder(react.id, -value), y = value, fill = react.id)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.1, position = position_dodge(0.5)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(p$react.id)))) +
    labs(x = "", y = "Total-order sensitivity index") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, family = "Helvetica"),
          axis.text.y = element_text(face = "plain", size = 10, family = "Helvetica"),
          axis.title.x = element_text(face = "plain", size = 10, family = "Helvetica"),
          axis.title.y = element_text(face = "bold", size = 9, family = "Helvetica"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          axis.line = element_line(size = 0.35))
  
  ggsave( (p_f / p_t), file = paste0(result_dir, "/SA_C_diff.pdf"), width = 7, height = 6)
  
  png(paste0(result_dir, "/key_r_indices.png"), height = 25*length(p$value), width = 650)
  gridExtra::grid.arrange(gridExtra::tableGrob(p))
  dev.off()
  
  # New facet label names for supp variable
  supp.labs = c("First-order indices", "Total-order indices")
  names(supp.labs) = c("first", "total")
  
  a = ggplot(dplyr::filter(df_ind, sol == "Y"), 
             aes(x = "", y = value, fill = order, color = order)) + 
    geom_boxplot(alpha = 0.2, size = 0.25) +
    coord_flip() +
    labs(title = paste0("Sensitivity indices evaluation | ", "config = " , nrow(param_values)), 
         subtitle = paste0("number of parameters : ", length(unique(df_ind$param))),
         x = "Output Y (biomass)", 
         y = "") +
    theme(legend.position="none",
          axis.ticks = element_line(size = 0.25),
          strip.background = element_blank(),
          line =  element_line(size = 0.75, color = "black"),
          axis.text = element_text(size = 5, color = "black"),
          axis.title = element_text(size = 6, face = "bold"),
          strip.text.x = element_text(size = 6, face = "bold", color = "#2A475E", hjust = 0.5),
          text = element_text(size = 8,  color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          plot.title = element_text(size = 8, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = 7, face = "bold", color = "#1b2838"),
          plot.title.position = "plot") +
    facet_wrap(~order, labeller = labeller(order = supp.labs)) +
    scale_fill_manual(values=c('blue', 'purple')) +
    scale_color_manual(values=c('black', 'black'))
  
  b = ggplot(dplyr::filter(df_ind, sol == "Y"), 
             aes(x = value, fill = order, color = order)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.2, size = 0.25) +
    labs(x = "Index Value", 
         y = "Count") +
    theme(legend.position="none",
          axis.ticks = element_line(size = 0.25),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          line =  element_line(size = 0.75, color = "black"),
          axis.text = element_text(size = 5, color = "black"),
          axis.title = element_text(size = 6, face = "bold"),
          text = element_text(size = 8,  color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          plot.title = element_text(size = 8, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = 6, face = "bold", color = "#1b2838"),
          plot.title.position = "plot") +
    facet_wrap(~order) + 
    scale_fill_manual(values=c('blue', 'purple')) +
    scale_color_manual(values=c('#1b2838', '#1b2838'))
  
  ggsave( (a/b) , file = paste0(result_dir, "/indices.pdf"), width = 3, height = 2.5)
  
  if( identical(which(all_react$React_ID == "sink_pheme_c"), integer(0)) ) {
    
    targetrnx = c("EX_trp_L_e", "EX_pro_L_e", 
                  "EX_val_L_e", "EX_ile_L_e", 
                  "EX_cys_L_e", "EX_leu_L_e",
                  "EX_lys_L_e",  "EX_biomass_e")
    
  } else {
    
    targetrnx = c("EX_trp_L_e", "EX_pro_L_e", 
                  "EX_val_L_e", "EX_ile_L_e", 
                  "EX_cys_L_e", "EX_leu_L_e", 
                  "EX_biomass_e", "sink_pheme_c")
    
  }
  
  ### Uncertainty of flux boundaries on a FBA model affects its objective
  
  fbasol = readRDS(paste0(result_dir, "/fbasol.rds"))
  fbasol_t = readRDS(paste0(result_dir, "/fbasol_template.rds"))
  
  subflux = data.frame(Reaction = rep(fbasol[, 1], (dim(fbasol)[2] - 1)),
                       Flux = gather(as.data.frame(fbasol), key, 
                                     value, -c(Reaction))$value,
                       config = str_remove(gather(as.data.frame(fbasol), key, 
                                                  value, -c(Reaction))$key, "config."))
  
  subflux = cbind(subflux, 
                  data.frame(lapply(param_values[, which(br$react.id %in% targetrnx)], 
                                    rep, each = dim(fbasol)[1])))
  
  subflux$Flux = as.double(subflux$Flux)
  
  df_plot = data.frame(ass.param = paste0("P.", which(br$react.id %in% targetrnx)),
                       react.index = which(br$react.id %in% targetrnx),
                       Reaction = br$react.id[which(br$react.id %in% targetrnx)],
                       col.data = seq(colnames(subflux[4:dim(subflux)[2]])) + 3,
                       col = brewer.pal(length(targetrnx),"Spectral"))
  
  subflux_t = cbind(data.frame(fbasol_t$abbreviation,
                               fbasol_t$flux,
                               rep(0, dim(fbasol)[1])),
                    data.frame(lapply(fbasol_t$lowbnd[df_plot$react.index],
                                      rep, each = dim(fbasol)[1])))
  
  colnames(subflux_t) = colnames(subflux)
  
  subflux = rbind(subflux, subflux_t)
  
  l = "#cc002d"; m = "#fffdbf"; h = "#112764"; c = "#653592"; pn = "darkred"
  
  title.size = 5; subtitle.size = 5; text.size = 5; axistext.size = 4; 
  axistitle.size = 5; legendsize = 0.25; point.size = 0.5
  
  PN.plots.violin = list()
  PN.plots.count = list()
  
  names(df_ind)[3] = "react.index"
  df_ind$react.index = as.integer(df_ind$react.index)
  
  # PN-related reactions indices
  PN_r = left_join(
    dplyr::filter(df_ind, sol == "Y" & 
                    parameter %in% df_plot$ass.param)[, c(2, 3, 5)],
    df_plot, by = "react.index")[, c(1, 3, 4, 5)]
  
  PN_r = PN_r[order(PN_r$value, decreasing = TRUE), ]
  
  for (i in 1:length(df_plot$Reaction)) {
    
    tmp = subflux %>% dplyr::filter(Reaction == df_plot$Reaction[i])
    tmp = dplyr::filter(tmp, config != "0")
    
    PN.plots.violin[[i]] = ggplot(dplyr::filter(tmp, config != "0"), 
                                  aes(x = Reaction, y = Flux)) +
      labs(x = df_plot$Reaction[i], y = "Flux (mmol/gDW*h)",
           title = paste0(df_plot$Reaction[i]),
           subtitle = paste0("number of reactions varied: ", 
                             length(param_values)),
           colour = names(tmp)[i + 3]) +
      theme(text = element_text(size = text.size, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
            plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
            plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
            plot.title.position = "plot",
            legend.key.size = unit(legendsize, 'cm'),
            axis.text = element_text(size = axistext.size, color = "black"),
            axis.title = element_text(size = axistitle.size, face = "bold")) +
      scale_color_gradient2(
        midpoint = mean(dplyr::filter(tmp, config != "0")[, names(tmp)[i + 3]]),
        low = l, 
        mid = m, 
        high = h) +
      geom_boxplot(alpha = 0.5) +
      geom_violin(scale = "width", alpha = 0.3, width = 1) +
      scale_x_discrete(name = "", 
                       label = rep( paste0(paste0(tmp$Reaction, "\n"),
                                           "(n = ", (dim(param_values)[1]), ")"),
                                    length(tmp$Reaction)))
    
    PN.plots.count[[i]] = ggplot(tmp, 
                                 aes(x = Flux, fill = ordered(Reaction, levels = unique(Reaction)))) + 
      geom_histogram(aes(y = ..count..), alpha = 0.5, bins = 25, fill = df_plot$col[i]) +
      labs(x = "Flux (mmol/gDW*h)", y = "Count",
           subtitle = paste0("number of reactions varied: ", length(param_values)),
           title = paste0("Fluxes distribution | ", df_plot$Reaction[i])) +
      guides(fill = guide_legend(title = "Condition block")) +
      theme(text = element_text(size = 8, color = "black"),
            panel.background = element_blank(),
            legend.position = "none",
            panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
            plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
            plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
            plot.title.position = "plot", 
            axis.text = element_text(size = axistext.size, color = "black"),
            axis.title = element_text(size = axistitle.size, face = "bold"))
    
  }
  
  png(paste0(result_dir, "/PN_r_indices.png"), height = 25*length(p$value), width = 350)
  gridExtra::grid.arrange(gridExtra::tableGrob(PN_r))
  dev.off()
  
  names(PN.plots.violin) = df_plot$Reaction
  names(PN.plots.count) = df_plot$Reaction
  
  ggsave(((PN.plots.violin[[1]] | PN.plots.violin[[2]] | 
             PN.plots.violin[[3]] | PN.plots.violin[[4]]) / 
            (PN.plots.violin[[5]] | PN.plots.violin[[6]] | 
               PN.plots.violin[[7]] | PN.plots.violin[[8]])),
         file = paste0(result_dir, "/FBA_Sens1.pdf"), 
         width = 5, height = 3)
  
  ggsave((PN.plots.count[[1]] | PN.plots.count[[2]] | 
            PN.plots.count[[3]] | PN.plots.count[[4]]) / 
           (PN.plots.count[[5]] | PN.plots.count[[6]] | 
              PN.plots.count[[7]] | PN.plots.count[[8]]),
         file = paste0(result_dir, "/FBA_Sens2.pdf"), 
         width = 6, height = 3.5)
  
  subflux = dplyr::filter(subflux, config != "0")[, c(1, 2, 3)]
  
  subflux.boun = 
    dplyr::filter(subflux, subflux$Reaction %in% 
                    dplyr::filter(all_react, ReactionDistrict == "boundary")$React_ID)
  
  dall_wrap = ggplot(subflux.boun, 
                     aes(x = Flux, fill = ordered(Reaction, levels = unique(Reaction)))) + 
    geom_histogram(aes(y = ..count..), alpha = 0.25, 
                   position = "identity", bin = 30,
                   binwidth = 0.05, col = c, fill = c) +
    labs(x = "Flux (mmol/gDW*h)", y = "Count") +
    guides(fill = guide_legend(title = "Condition block")) +
    theme(text = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
          plot.title = element_text(size = title.size, face = "bold", color = "#2a475e"),
          plot.subtitle = element_text(size = subtitle.size, face = "bold", color = "#1b2838"),
          plot.title.position = "plot", 
          axis.text = element_text(size = axistext.size, color = "black"),
          axis.title = element_text(size = axistitle.size, face = "bold")) +
    facet_wrap(~ Reaction, scales = "free")
  
  ggsave(dall_wrap, 
         file = paste0(result_dir, "/FBA_Sens_allreact.pdf"),
         width = 15, height = 15, limitsize = FALSE)
  
}
