
# init.gen: generate the initial marking for your PNPRO places
init.gen <- function() {
  # place names, in PNPRO order
  yini.names <- c('N_c1', 'N_c2', 'biomass_e_c1', 'biomass_e_c2')
  
  # initial marking vector (same length as yini.names)
  y_ini <- c(1e+06, 1e+06, 1, 1)
  
  names(y_ini) <- yini.names
  
  return(y_ini)
}

