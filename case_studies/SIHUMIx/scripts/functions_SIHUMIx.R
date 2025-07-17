
# init.gen: generate the initial marking for your PNPRO places
init.gen_SIHUMIx <- function() {
  # place names, in PNPRO order
  yini.names <- c('N_c1', 'N_c2', 'N_c3', 'N_c4', 'N_c5', 'N_c6', 'N_c7', 'N_c8',
                  'biomass_e_c1', 'biomass_e_c2', 'biomass_e_c3', 'biomass_e_c4', 'biomass_e_c5', 'biomass_e_c6', 'biomass_e_c7', 'biomass_e_c8',
                  'lac__L_e','ac_e','but_e','ppa_e')
  
  # initial marking vector (same length as yini.names)
  y_ini <- c(1e+06, 1e+06, 1e+06, 1e+06, 1e+06, 1e+06, 1e+06, 1e+06,
             1, 1, 1, 1, 1, 1, 1, 1, 
             5, 5, 5, 5)
  
  names(y_ini) <- yini.names
  
  return(y_ini)
}
