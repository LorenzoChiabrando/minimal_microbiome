# -------------------------------------------------------------------------
# run.R – one-click demo for the “minimal-doublet” community
# -------------------------------------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github(
  "qBioTurin/epimod_FBAfunctions",
  ref      = "unified-epimod_FBAfunctions",
  upgrade  = "never"
)

library(epimodFBAfunctions)
# remove.packages("epimod")
# devtools::install_github("https://github.com/qBioTurin/epimod", ref="epimod_pFBA")
library(epimod)
# downloadContainers()

# -------------------------------------------------------------------------
# user-editable knobs
# -------------------------------------------------------------------------
hypernode_name <- "EcCb_SCFA"
base_dir       <- getwd()                       # where you launch run.R
mat_dir        <- "models"                      # user *.mat* live here
cfg_yaml       <- paste0("config/config_", hypernode_name, ".yaml")
bc_json        <- paste0("config/boundary_conditions_", hypernode_name, ".json")
init_csv       <- paste0("config/initial_data_", hypernode_name, ".csv")
overwrite_run  <- TRUE
debug_solver   <- FALSE                         # epimod::model.analysis()

# -------------------------------------------------------------------------
# 1. build the hyper-node (returns all useful sub-paths)
# -------------------------------------------------------------------------
paths <- build_hypernode(
  hypernode_name           = hypernode_name,
  config_yaml              = cfg_yaml,
  boundary_conditions_file = bc_json,
  initial_data             = init_csv,
  mat_dir                  = mat_dir,
  base_dir                 = base_dir,
  overwrite                = overwrite_run,
  debug                    = debug_solver          # passed to epimod later
)

cat("\n✔  Hyper-node ready in", fs::path_rel(fs::path_dir(paths$config), base_dir), "\n")

# -------------------------------------------------------------------------
# 2. run epimod solver
# -------------------------------------------------------------------------
net_file   <- fs::path(base_dir, "hypernodes/", hypernode_name, "petri_net/", paste0(hypernode_name, ".PNPRO"))
trans_file <- fs::path(paths$src,   paste0("general_functions_", hypernode_name, ".cpp"))

fba_files <- fs::dir_ls(paths$biounit,
                        recurse      = TRUE,
                        regexp       = "_model\\.txt$",
                        type         = "file")

epimod::model.generation(
  net_fname        = net_file,
  transitions_fname= trans_file,
  fba_fname        = fba_files
)

# -------------------------------------------------------------------------
# 3. move the solver artefacts into hypernodes/<name>/gen/
# -------------------------------------------------------------------------
solver_suffixes <- c(".solver", ".def", ".fbainfo", ".net", ".PlaceTransition")

fs::dir_create(paths$gen)        # already created inside build_hypernode(), but harmless

purrr::walk(solver_suffixes, function(ext) {
  src <- fs::path(base_dir, paste0(hypernode_name, ext))
  if (fs::file_exists(src)) {
    fs::file_move(src, paths$gen)
    message("✓ moved ", fs::path_file(src))
  } else {
    warning("solver file not found: ", src)
  }
})

# -------------------------------------------------------------------------
# 4. analyse the simulation
# -------------------------------------------------------------------------
epimod::model.analysis(
  solver_fname     = fs::path(paths$gen, paste0(hypernode_name, ".solver")),
  parameters_fname = fs::path(paths$config, paste0("initial_data_", hypernode_name, ".csv")),
  functions_fname  = fs::path(paths$src,   paste0("functions_", hypernode_name, ".R")),
  debug            = debug_solver,
  i_time = 0, f_time = 24, s_time = 0.25,
  atol  = 1e-6, rtol = 1e-6,
  fba_fname = fba_files,
  user_files = c(
    fs::path(paths$config, "population_parameters.csv"),
    fs::path(paths$gen,    paste0(hypernode_name, ".fbainfo")),
    fs::path(paths$output, "ub_bounds_projected.csv"),
    fs::path(paths$output, "ub_bounds_not_projected.csv")
  )
)


cat("\n🎉 All results gathered in", fs::path_rel(paths$gen, base_dir), "\n")
