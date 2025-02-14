#' @title Read matlab model
#' @description The generic function \code{readMATmod} imports matlab cobra (ONLY) models into ....
#' @param fba_mat Full path to matlab model file
#' @details Returns ...
#' @author Aucello Riccardo, Beccuti Marco, Cordero Francesca, Pernice Simone
#' @import R.matlab stringr
#' @export

FBAmat.read = function(fba_mat){
  
  data = R.matlab::readMat(fba_mat)
  dat.mat = data[[1]]
  
  mod.var = dimnames(dat.mat)[[1]]
  
  # 1) model name
  if( "modelID" %in% mod.var ){
    mod.id = as.character(dat.mat[[which(mod.var == "modelID")]])
  } else {
    mod.id = names(data)[1]
  }
  if( "modelName" %in% mod.var ){
    mod.name = as.character(dat.mat[[which(mod.var == "modelName")]])
  } else {
    mod.name = mod.id
  }
  if( "description" %in% mod.var ){
    mod.desc = as.character(dat.mat[[which(mod.var == "description")]])
  } else {
    mod.desc = mod.id
  }
  
  # 2) stoich matrix
  mod.S = Matrix(dat.mat[[which(mod.var == "S")]], sparse = T)
  
  # 3) rxn
  mod.react_id = unlist(dat.mat[[which(mod.var == "rxns")]])
  mod.react_name = unlist(dat.mat[[which(mod.var == "rxnNames")]])
  
  saveRDS(mod.react_name, file = "officialName.rds")
  
  if( "rev" %in% mod.var ){
    mod.react_rev = as.vector(dat.mat[[which(mod.var=="rev")]]) == TRUE
  } else {
    mod.react_rev = as.vector(dat.mat[[which(mod.var=="lb")]]) < 0
  }
  
  # 4) met
  mod.met_id     <- unlist(dat.mat[[which(mod.var=="mets")]])
  mod.met_name   <- unlist(dat.mat[[which(mod.var=="metNames")]])
  
  # 5) genes
	mod.gene_id <- unname(unlist(dat.mat[[which(mod.var=="genes")]]))

	if ("rxnGeneMat" %in% mod.var) {
		mod.GeneMat <- dat.mat[[which(mod.var == "rxnGeneMat")]]
		
		mod.genes <- apply(mod.GeneMat, 1, function(row) {
		  x <- unname(mod.gene_id[which(row != 0)])
		  if (length(x) == 0) "" else x
		})
		mod.gpr <- sapply(sapply(dat.mat[[which(mod.var == "grRules")]], unlist),
		                  function(entry) {
		                    if (length(entry) == 0) "" else unname(entry)
		                  })
		
		# Calcolo isGeneAssociated
		isGeneAssociated <- apply(mod.GeneMat != 0, 1, any)   # vettore TRUE/FALSE di lunghezza nReactions
		isGeneAssociated_num <- as.numeric(isGeneAssociated)  # convertito in 0/1

	} else {
		# if gene matrix not present create it
		mod.GeneMat <- Matrix(0, nrow = length(mod.react_id), ncol = length(mod.gene_id))
		mod.genes <- list()
		rules.list <- unlist(dat.mat[[which(mod.var == "rules")]], recursive = FALSE)
		
		if (length(rules.list) != length(mod.react_id))
		  stop("Length of rules not same as length of reactions.")
		
		for (i in seq_along(mod.react_id)) {
		  rule.tmp <- rules.list[[i]]
		  if (length(rule.tmp) == 0) {
		    mod.genes[[i]] <- ""
		    next
		  }
		  j <- as.numeric(
		    unlist(stringr::str_extract_all(rule.tmp, "(?<=x\\()[0-9]+?(?=\\))"))
		  )
		  mod.GeneMat[i, j] <- 1
		  mod.genes[[i]] <- mod.gene_id[j]
		}
		
		# Calcolo isGeneAssociated
		isGeneAssociated <- apply(mod.GeneMat != 0, 1, any)
		isGeneAssociated_num <- as.numeric(isGeneAssociated)
	}

	saveRDS(mod.gene_id, file = "geni.rds")

  
  # gpr rules needs to be converted (brackets + numbering)
	if( "grRules" %in% mod.var ){
		mod.gprRules <- sapply(sapply(dat.mat[[which(mod.var=="grRules")]], unlist), function(entry){
		  if (length(entry) == 0) {
		    ""
		  } else {
		    numbers <- as.numeric(unlist(stringr::str_extract_all(entry, "[0-9]+")))
		    
		    if (length(numbers) == 0) {
		      ""
		    } else {
		      dict <- as.character(numbers - min(numbers) + 1)
		      names(dict) <- as.character(numbers)
		      gsub("\\(([0-9]+)\\)","\\[\\1\\]", stringr::str_replace_all(entry, dict))
		    }
		  }
		})
	}else{ # if 'rules' is not present construct own rules from gpr+genes
    mod.gprRules <- sapply(seq_along(mod.gpr), function(i){
      if( mod.gpr[i] == "") ""
      else{
        genes <- unlist(mod.genes[i])
        dict <- paste0("x[",seq_along(genes),"]")
        names(dict) <- genes
        gsub("and","&", gsub("or","|",stringr::str_replace_all(mod.gpr[i], dict)))
      }
    })
    
  }
  
  saveRDS(mod.gprRules, file = "geneAssociation.rds")
  
  genes = c()
  
  for (i in 1:length(mod.gprRules)) {
    
    txt = mod.gprRules[i]
    txt1 = gsub("\\(",replacement = "", txt )
    txt1 = gsub("\\)",replacement = "", txt1 )
    
    genes = c(genes, unlist(strsplit(txt1, "( and | or )")))
    
  }
  
  genes = unique(genes)
  genes = gsub(" ", replacement = "", genes)
  genes = unlist(strsplit(genes, "(and|or)"))
  genes = unique(genes)
  
  saveRDS(genes, file = "genesFromGeneAss.rds")
  
  # 6) bounds
  mod.lb <- as.vector(dat.mat[[which(mod.var=="lb")]])
  mod.ub <- as.vector(dat.mat[[which(mod.var=="ub")]])
  
  # 7) compartments
  met_comp <- stringr::str_extract_all(mod.met_id, "(?<=\\[)[a-z](?=\\])")
  if( all(sapply(met_comp, length) == 0 ) ){
    met_comp <- stringr::str_extract_all(mod.met_id, "(?<=_)[a-z][0-9]?(?=$)")
  }
  mod.mod_compart <- unique(unlist(met_comp))
  mod.met_comp    <- match(met_comp, mod.mod_compart)
  
  # 8) subsystems
  
  sub = sapply(dat.mat[[which(mod.var == "subSystems")]], unlist)
  sub.unique = unique(sub)
  
  mod.subSys = Matrix(FALSE, nrow = length(sub), ncol = length(sub.unique), sparse=T)
  
  for(i in 1:length(sub)){
    
    j <- match(sub[i], sub.unique)
    mod.subSys[i,j] <- TRUE
    
  }
  
  colnames(mod.subSys) <- sub.unique
  
  saveRDS(sub, file = "subsystem.rds")
  
  # 9) metKEGGID
  
  if( "metKEGGID" %in% mod.var ) {
    
    met.KEGGID = dat.mat[[which(mod.var=="metKEGGID")]]
    met_KEGGID = c()
    
    for(i in 1:length(met.KEGGID)) {
      if(length(met.KEGGID[[i]][[1]]) == 0) {
        met_KEGGID[i] = NA
      } else {
        met_KEGGID[i] = met.KEGGID[[i]][[1]]
      }
    }
    
    met.KEGGID = data.frame(met_id = mod.met_id, met_KEGGID)
    saveRDS(met.KEGGID, file = "met_KEGGID.rds")
    
  }
  
  # 10) rxnKEGGID
  
  if( "rxnKEGGID" %in% mod.var ) {
    
    rxn.KEGGID = dat.mat[[which(mod.var=="rxnKEGGID")]]
    rxn_KEGGID = c()
    
    for(i in 1:length(rxn.KEGGID)) {
      if(length(rxn.KEGGID[[i]][[1]]) == 0) {
        rxn_KEGGID[i] = NA
      } else {
        rxn_KEGGID[i] = rxn.KEGGID[[i]][[1]]
      }
    }
    
    rxn.KEGGID = data.frame(react_id = mod.react_id, rxn_KEGGID)
    saveRDS(rxn.KEGGID, file = "rxn_KEGGID.rds")
    
  }
  
  # 11) rxnKeggOrthology
  
  if( "rxnKeggOrthology" %in% mod.var ) {
    
    rxn.KeggOrthology = dat.mat[[which(mod.var=="rxnKEGGID")]]
    rxn_KeggOrthology = c()
    
    for(i in 1:length(rxn.KeggOrthology)) {
      if(length(rxn.KeggOrthology[[i]][[1]]) == 0) {
        rxn_KeggOrthology[i] = NA
      } else {
        rxn_KeggOrthology[i] = rxn.KeggOrthology[[i]][[1]]
      }
    }
    
    rxn.KeggOrthology = data.frame(react_id = mod.react_id, rxn_KeggOrthology)
    saveRDS(rxn.KeggOrthology, file = "rxn_KeggOrthology.rds")
    
    # see -> https://webdav-r3lab.uni.lu/public/msp/Mol-Files/atomMappedExplicit/png/
    
  }
  
  # create new model
  # model = FBAgreatmod.classGeneration()
  #
  # slot(model,"S") <- as.matrix(mod.S)
  # slot(model,"lowbnd") <- mod.lb
  # slot(model,"uppbnd") <- mod.ub
  # slot(model,"react_id")  <- mod.react_id
  # slot(model,"obj_coef") <- c(dat.mat[[which(mod.var == "c")]])
  # model$mod_desc <- mod.desc
  # model$met_id <- mod.met_id
  # model$met_name <- mod.met_name
  # model$met_num <- length(mod.met_id)
  # model$react_name  <- mod.react_name
  # model$react_num <- length(mod.react_id)
  # model$react_rev <- mod.react_rev
  # model$genes <- mod.genes
  # model$gprRules <- mod.gprRules
  # model$gpr <- mod.gpr
  # model$allGenes <- mod.gene_id
  # model$mod_compart <- mod.mod_compart
  # model$met_comp <- mod.met_comp
  # model$subSys <- mod.subSys
  
  mod.react_id = gsub("\\(", replacement = "_", mod.react_id)
  mod.react_id = gsub("\\[", replacement = "_", mod.react_id)
  mod.react_id = gsub("\\)", replacement = "", mod.react_id)
  mod.react_id = gsub("\\]", replacement = "", mod.react_id)
  mod.react_id = gsub("_c_", replacement = "_c", mod.react_id)
  
  return(
    list(
      S = mod.S,
      lowbnd = mod.lb,
      uppbnd = mod.ub,
      react_id = mod.react_id,
      met_id = mod.met_id,
      obj_coef = dat.mat[[which(mod.var == "c")]],
      gene_assoc = isGeneAssociated_num
    )
  )
}
