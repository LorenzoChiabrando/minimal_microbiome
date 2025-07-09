library(epimodFBAfunctions)
library(readr)
library(dplyr)

### Data files
fba_mat_file <- system.file("data/MATmodels", "Ec_core.mat", package = "epimodFBAfunctions")
diet_file <- system.file("data/diets/vmh", "EU_average.tsv", package = "epimodFBAfunctions")

### Model
model = FBA4Greatmod.generation(fba_mat = fba_mat_file)

### Diet integration
diet = read_delim(diet_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)

diet$Reaction = gsub("\\[", replacement = "\\(", diet$Reaction)
diet$Reaction = gsub("\\]", replacement = "\\)", diet$Reaction)

diet = diet %>%
  rename(lwbnd = `Flux Value`) %>%
  mutate(lwbnd = -(lwbnd)/24, uppbwnd = 10)

model = setDiet(model, dietf = diet)
model = setDiet.name(model,diet_name = "EU_average")

for(r in model@react_id) model = setConstraints(model, reaction.name = r, newConstraints = c(-0.1, 10))

saveRDS(model, "./Ec_core.RDs")
writeFBAfile(model, "./Ec_core")
