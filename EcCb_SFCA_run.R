# -------------------------------------------------------------------------
# run.R – one-click demo for the “minimal-doublet” community
# -------------------------------------------------------------------------

# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")

# remove.packages("epimodFBAfunctionsGUI")
# remove.packages("epimod")
# remove.packages("epimodFBAfunctions")
remotes::install_github(
  "qBioTurin/epimod_FBAfunctions",
  ref      = "unified-epimod_FBAfunctions",
  upgrade  = "never"
)

remotes::install_github(
  "https://github.com/LorenzoChiabrando/epimodFBAfunctions_GUI",
  ref = "main")

remotes::install_github("https://github.com/LorenzoChiabrando/epimod_gui", ref="main")

library(epimodFBAfunctionsGUI)
library(epimodFBAfunctions)
# install.packages("devtools")
library(devtools)
library(epimod)
# downloadContainers()

epimodFBAfunctionsGUI::run_app()

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

source(file.path("plot_flux_pca.R"))
source(file.path("plot_all_results.R"))
source(file.path("plotting_int_net.R"))
source(file.path("plot_metadata.R"))

f_time = 24
s_time = 0.01

col_bac <- c("#341539", "#FFE60C")

generate_metadata_plots(hypernode_name, col_bac) 
  
my_col_met_places <- c("#282728", "#eaa380", "#ad4233", "#8887cd", "#086")

my_react2plot <- c("EX_biomass_e_f", "EX_biomass_e_r",
                   "EX_glc__D_e_f", "EX_glc__D_e_r", "EX_lac__D_e_r", "EX_lac__D_e_f",
                   "EX_ppa_e_f", "EX_ac_e_f", "EX_but_e_f", "EX_for_e_f",
                   "EX_ppa_e_r", "EX_ac_e_r", "EX_but_e_r", "EX_for_e_r")

my_scfa_list <- c("ac_e", "ac_c", "ppa_e", "ppa_c", "but_e", "but_c", "for_e", "for_c",
                  "M03134_e", "M03134_c", "caproic_e", "caproic_c",
                  "isobut_e", "isobut_c", "isoval_e", "isoval_c",
                  "isocapr_e", "isocapr_c", "isobut_e", "isobut_c")

my_entities_list <- c(my_scfa_list, "glu__L_e", "glu__L_c", 
                      "lac__L_e", "lac__L_c", "lac__D_e", "lac__D_c",
                      "ade_e", "ade_c", "gua_e", "gua_c", "nac_e", "nac_c",
                      "thymd_e", "thymd_c", "ptrc_e", "ptrc_c")

plot_marking_and_flux_trends(case_name = hypernode_name, 
                             col_bac = col_bac,
                             col_met_places = my_col_met_places,
                             react2plot = my_react2plot,
                             num_sampling_points_rel_abun = 12,
                             num_sampling_points_met_plots = 12,
                             plot_filename_suffix = ".pdf"
)

my_flux_sampling_times <- unique(round(seq(0, f_time, length.out = f_time*10), round(abs(log10(s_time)), 1)))

plot_flux_pca(
  case_name = hypernode_name,
  flux_sampling_times = my_flux_sampling_times,
  flux_th_l = -1000,
  flux_th_h = 1000,
  col_bac = col_bac,
  scfa_list = my_scfa_list,
  entities_list = my_entities_list,
  plot_filename = paste0(hypernode_name, "_pca_output.pdf")
)

my_cross_fed_met_list <- c("lcts_e", "glc__D_e", "ac_e","ppa_e", "but_e",
                           "M03134_e", "caproic_e", "isobut_e",
                           "isoval_e", "isocapr_e", "isobut_e",
                           "glu_L_e", "lac__L_e", "lac__D_e", "for_e", "ade_e",
                           "gua_e", "nac_e", "thymd_e", "ptrc_e")

plotting_int_net(
  case_name = hypernode_name,
  t_frame = my_flux_sampling_times,
  col_bac = "#f06b368c",
  col_met = "lightblue",
  shape_bac = "square",
  shape_met = "circle",
  cross_fed_met = my_cross_fed_met_list,
  output_filename_prefix = paste0(hypernode_name, "_CrossFeeding"),
  plot_filename_suffix = ".pdf"
)
