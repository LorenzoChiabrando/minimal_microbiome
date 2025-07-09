
# libraries required
library(dplyr)
library(tidyr)
library(epiR)
library(tibble)
library(ggplot2)
library(readr)
library(sensitivity)
library(pander)
library(RColorBrewer)
library(patchwork)
library(fbar)
library(sensobol)
library(tibble)
library(foreach)
library(doParallel)
library(fbar)
library(stringr)

# setting working directory
# "./epimod_globalSAfunctions
wd = getwd()
setwd(wd)

# setting results directory
result_dir = paste0(wd, "/inst/res")
input_dir = paste0(wd, "/inst/input")

# loading epimod_globalSA functions
source(paste0(wd, "/R/readMat.R"))
source(paste0(wd, "/R/ParallelFBA_sens.R"))
source(paste0(wd, "/R/FBAgreatmodeClass.R"))
source(paste0(wd, "/R/class_generation.R"))
source(paste0(wd, "/R/SA_FBA_R.R"))

# setting GEMM tag
model.name = "CD196HemeSink"

# setting analysis FBA files names
mat_file = paste0(wd, "/inst/data/MATmodels/", model.name, ".mat")
model_file = paste0(input_dir, "/FBAmodel.RData")
# generate S3 R object FBA_model
model = FBA4Greatmod.generation(fba_mat = mat_file)
# saving RData FBA model
save(model, file = model_file)

files = c("geneAssociation.rds", 
          "geni.rds", "officialName.rds", 
          "subsystem.rds", "genesFromGeneAss.rds",
          "met_KEGGID.rds", "rxn_KEGGID.rds")

for(f in files) {
  system(paste0("cd ", wd, " && ", "mv ", f, " ", input_dir, "/", f))
}

metadataEQ(model = model, 
           model.name = model.name,
           prefix = "BIGGdata_FBAmodel",
           suffix = "tsv",
           extMetFlag = "b", 
           fielddelim = "\t", 
           entrydelim = ", ", 
           makeClosedNetwork = FALSE, 
           onlyReactionList = FALSE, 
           minimalSet = TRUE, 
           fpath = input_dir)

bigg.path = paste0(input_dir, "/BIGGdata_FBAmodel_react.tsv")
all_react = FBAmodel.metadata(model = model,
                              wd = wd, 
                              bigg.path = bigg.path)

all_react_file = paste0(input_dir, "/all_react.rds")
saveRDS(all_react, all_react_file)

# saving biochemical reaction equations
BIGGdata = read.delim2(bigg.path)

saveRDS(BIGGdata$equation, file = paste0(input_dir, "/equation.rds"))

# We will investigate the sensitivty of the reaction objective ("biomass205") 
# flux solution by varying the D parameters within the bounds (-10.0, 0.0) 
# defined by the problem shown below.
# 
# I’m going to use saltelli to generate a collection of parameter values to use for FBA. 
# Then sobol to analyze the results and generate the sensitivity analysis. 
# To use sample.saltelli, you have to provide the problem definition in a specific format.

SA_FBA(result_dir = result_dir,
       model_file = paste0(input_dir, "/FBAmodel.RData"),
       all_react_file = paste0(input_dir, "/all_react.rds"),
       param_values_file = paste0(result_dir, "/param_values.rds"),
       bounds_set_file = paste0(result_dir, "/bounds_set.rds"),
       fbasol_t_file = paste0(result_dir, "/fbasol_template.rds"),
       fbasol_file = paste0(result_dir, "/fbasol.rds"),
       Y_file = paste0(result_dir, "/Y.rds"),
       Y2_file = paste0(result_dir, "/Y2.rds"),
       indices_file = paste0(result_dir, "/df_ind.rds"),
       Y2_reaction = "EX_h2o_e",
       # Define settings
       # N = 2^13
       N = 2^13, 
       cores = 5,
       b1 = -10.0,
       b2 = 0.0)

param_values = readRDS(paste0(result_dir, "/param_values.rds"))
df_ind = readRDS(paste0(result_dir, "/df_ind.rds"))

df_ind_filtered = dplyr::filter(df_ind, sol == "Y" & value > 0.01)
df_ind_filtered$param = as.integer(df_ind_filtered$param)

# model annotation
all_react = readRDS(all_react_file)
br = all_react[[2]]
all_react = all_react[[1]]

br$param = seq(br$react.id)

p = left_join(df_ind_filtered, br[, c(2, 6)], by = "param")
p = p[order(p$value, decreasing = TRUE), ]

mat_list = split(
  param_values,
  rep(1:(nrow(param_values) %/% chunks + 1),
      each = chunks, length.out = nrow(param_values)))

for(i in 1:length(mat_list)) {
  
  bounds_set = readRDS(paste0(result_dir, "/bounds_set", "_", i, ".rds"))
  
  bounds_set_lb = bounds_set[, 1:(ncol(bounds_set)/2)]
  bounds_set_ub = bounds_set[, ((ncol(bounds_set)/2) + 1):ncol(bounds_set)]
  colnames(bounds_set_lb) = NULL
  colnames(bounds_set_ub) = NULL
  
  fbasol = readRDS(paste0(result_dir, "/fbasol", "_", i, ".rds"))
  colnames(fbasol) = NULL
  
  if(i == 1) {
    
    bounds_set_lb_t = bounds_set_lb
    bounds_set_ub_t = bounds_set_ub
    fbasol_t = fbasol
    
  } else {
    
    bounds_set_lb_t = cbind(bounds_set_lb_t, bounds_set_lb)
    bounds_set_ub_t = cbind(bounds_set_ub_t, bounds_set_ub)
    fbasol_t = cbind(fbasol_t, fbasol[, -1])
    
  }
}

colnames(bounds_set_lb_t) = paste("lb", 1:ncol(bounds_set_lb_t), sep = ".")
colnames(bounds_set_ub_t) = paste("ub", 1:ncol(bounds_set_ub_t), sep = ".")
saveRDS(cbind(bounds_set_lb_t, bounds_set_ub_t), paste0(result_dir, "/bounds_set.rds"))

colnames(fbasol_t) = c("Reaction", paste("config", 2:ncol(fbasol_t), sep = "."))
saveRDS(fbasol_t, paste0(result_dir, "/fbasol.rds"))

# list_files = list.files(result_dir)
# 
# file.remove(paste0(result_dir, "/", list_files[grepl("bounds_set", list_files)]))
# v_fbasol = paste0(result_dir, "/", list_files[grepl("fbasol", list_files)])
# file.remove(v_fbasol[-length(v_fbasol)])

plotting_indeces(p, result_dir)
  
distr_stats(result_dir = result_dir,
            fbasol_file = paste0(result_dir, "/fbasol.rds"),
            file_sd = paste0(result_dir, "/data.RData"),
            thrsd = 10,
            linethrsd = 100)
