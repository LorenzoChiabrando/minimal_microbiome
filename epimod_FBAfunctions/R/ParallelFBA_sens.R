
ParallelFBA_sens = function(model_file,
                            param_values_file,
                            model_metadata_file,
                            result_dir,
                            input_dir,
                            reactIndx4secondSol,
                            cores,
                            wd) {
  
  # loading FBA model
  load(model_file)
  
  # loading FBA model boundary reactions
  all_react = readRDS(model_metadata_file)
  br = all_react[[2]]
  
  equation = readRDS(paste0(input_dir, "/equation.rds"))
  officialName = readRDS(paste0(input_dir, "/officialName.rds"))
  geneAssociation = readRDS(paste0(input_dir, "/geneAssociation.rds"))
  
  subsystem = readRDS(paste0(input_dir, "/subsystem.rds"))
  subsystem[which(sapply(1:length(subsystem), function(i) length(subsystem[[i]])) == 0)] = ""
  
  # template solution
  fbasol_template_file = paste0(result_dir, "/fbasol_template.rds")
  
  modelFBA.t = data.frame(abbreviation = model@react_id,
                          lowbnd = model@lowbnd,
                          uppbnd = model@uppbnd,
                          obj_coef = model@obj_coef,
                          equation = equation,
                          officialName = officialName,
                          geneAssociation = geneAssociation,
                          subsystem = unlist(subsystem))
  
  saveRDS(find_fluxes_df(modelFBA.t), fbasol_template_file)
  
  # loading the parameter set generated with saltelli sampling
  param_values = readRDS(param_values_file)
  
  bounds_set = list()
  outputs = list()
  results = list()
  second_results = list()
  
  # Define the function to be applied to each element of the iterator
  ParFBA = function(i) {
    
    # loading functions
    source(paste0(wd, "/R/class_generation.R"))
    source(paste0(wd, "/R/FBAgreatmodeClass.R"))
    source(paste0(wd, "/R/readMat.R"))
    
    sample = param_values[i, ]
    
    for(r in 1:length(br$react.id)) {
      model = setConstraints(
        model,
        reaction.name = br$react.id[r],
        newConstraints = c(as.double(sample[r]),
                           br$ub[which(br$react.id == br$react.id[r])]))
    }
    
    modelFBA = data.frame(abbreviation = model@react_id,
                          lowbnd = model@lowbnd,
                          uppbnd = model@uppbnd,
                          obj_coef = model@obj_coef,
                          equation = equation,
                          officialName = officialName,
                          geneAssociation = geneAssociation,
                          subsystem = unlist(subsystem))
    
    fbasol = find_fluxes_df(modelFBA)
    
    outputs = fbasol
    results = fbasol$flux[which(fbasol$obj_coef == 1)]
    second_results = fbasol$flux[reactIndx4secondSol]
    
    return(list(outputs, results, second_results))
    
  }
  
  # initialising the list to pass to the apply functions
  i.list = sapply(seq_along(vector("list", nrow(param_values))), list)
  
  # detect the number of cores
  n.cores = cores
  
  system.time({
    clust = makeCluster(n.cores, type = "FORK")
    clusterExport(cl = clust, 
                  varlist = c("wd", "result_dir", "model", "all_react",
                              "param_values", "br", "reactIndx4secondSol",
                              "equation", "officialName", "geneAssociation", "subsystem",
                              "outputs", "results", "second_results", 
                              "find_fluxes_df"),
                  envir = environment())
    sol = parLapply(clust, i.list, ParFBA)})
  
  stopCluster(clust)
  
  bounds_set = data.frame(lb = sapply(sol, function(x) x[[1]][["lowbnd"]]),
                          ub = sapply(sol, function(x) x[[1]][["uppbnd"]]))
  
  # saving flux distribution
  fluxes = sapply(sol, function(x) x[[1]][["flux"]])
  all_react = all_react[[1]]
  fluxes = cbind(all_react$React_ID, fluxes)
  colnames(fluxes) = c("Reaction", paste("config", 1:(dim(fluxes)[2]-1), sep = "."))
  
  Y = sapply(sol, function(x) x[[2]])
  Y2 = sapply(sol, function(x) x[[3]])
  
  saveRDS(bounds_set, paste0(result_dir, "/bounds_set.rds")) 
  saveRDS(fluxes, paste0(result_dir, "/fbasol.rds"))
  saveRDS(Y, paste0(result_dir, "/Y.rds"))
  saveRDS(Y2, paste0(result_dir, "/Y2.rds"))
  
}
