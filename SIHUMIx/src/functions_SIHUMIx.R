# -----------------------------------------------------------------------------
# functions_hypernode_template.R
#
# Template for hypernode R functions.
# Run generate_R_from_pnpro() to fill in the two stub vectors below.
# -----------------------------------------------------------------------------

# init.gen: generate the initial marking for your PNPRO places
init.gen <- function() {
  # place names, in PNPRO order
  yini.names <- c('biomass_e_acd1', 'n_acd1', 'biomass_e_btv5', 'n_btv5', 'biomass_e_bln', 'n_bln', 'biomass_e_bpd2', 'n_bpd2', 'biomass_e_cbd1', 'n_cbd1', 'biomass_e_crv0d1', 'n_crv0d1', 'biomass_e_ecs', 'n_ecs', 'biomass_e_lpspa1', 'n_lpspa1', 'ac_e', 'but_e', 'lac__L_e', 'ppa_e')
  
  # initial marking vector (same length as yini.names)
  y_ini <- c(1, 1e+06, 1, 1e+06, 1, 1e+06, 1, 1e+06, 1, 1e+06, 1, 1e+06, 1, 1e+06, 1, 1e+06, 1, 1, 1, 1)
  
  # assign names & enforce PNPRO order
  names(y_ini) <- yini.names
  #y_ini <- y_ini[yini.names]
  
  return(y_ini)
}

# --------------------------------------------------------------
# (You can add other functions below that will use init.gen(), 
#  or additional functions for sensitivity, ect.
# --------------------------------------------------------------
