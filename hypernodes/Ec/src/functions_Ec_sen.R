# -----------------------------------------------------------------------------
# functions_hypernode_template.R
#
# Template for hypernode R functions.
# Run generate_R_from_pnpro() to fill in the two stub vectors below.
# -----------------------------------------------------------------------------

# init.gen: generate the initial marking for your PNPRO places
init.gen <- function() {
  # place names, in PNPRO order
  yini.names <- c('biomass_e_ecsk1sm', 'n_ecsk1sm')
  
  # initial marking vector (same length as yini.names)
  y_ini <- c(0.489, 1000)
  
  # assign names & enforce PNPRO order
  names(y_ini) <- yini.names
  y_ini <- c(0.489, 1000)
  
  return(y_ini)
}

psensitivty = function(n, min, max, dup, death) {
  starv = runif(n, min, max)
  return(matrix(c(starv, dup, death), ncol = 3))
}

# --------------------------------------------------------------
# (You can add other functions below that will use init.gen(), 
#  or additional functions for sensitivity, ect.
# --------------------------------------------------------------
