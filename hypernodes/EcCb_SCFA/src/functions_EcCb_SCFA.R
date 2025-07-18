# -----------------------------------------------------------------------------
# functions_hypernode_template.R
#
# Template for hypernode R functions.
# Run generate_R_from_pnpro() to fill in the two stub vectors below.
# -----------------------------------------------------------------------------

# init.gen: generate the initial marking for your PNPRO places
init.gen <- function() {
  # place names, in PNPRO order
  yini.names <- c('biomass_e_ecs', 'n_ecs', 'biomass_e_cbd1', 'n_cbd1', 'ac_e', 'but_e', 'for_e', 'glc__D_e', 'lcts_e', 'ppa_e')
  
  # initial marking vector (same length as yini.names)
  y_ini <- c(0.489, 1e+06, 0.6, 1e+06, 1, 1, 1, 1, 1, 1)
  
  # assign names & enforce PNPRO order
  names(y_ini) <- yini.names
  y_ini <- c(0.489, 1e+06, 0.6, 1e+06, 1, 1, 1, 1, 1, 1)
  
  return(y_ini)
}

# --------------------------------------------------------------
# (You can add other functions below that will use init.gen(), 
#  or additional functions for sensitivity, ect.
# --------------------------------------------------------------
