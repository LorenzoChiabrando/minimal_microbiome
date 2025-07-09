# -----------------------------------------------------------------------------
# functions_hypernode_template.R
#
# Template for hypernode R functions.
# Run generate_R_from_pnpro() to fill in the two stub vectors below.
# -----------------------------------------------------------------------------

# init.gen: generate the initial marking for your PNPRO places
init.gen <- function() {
  # place names, in PNPRO order
  yini.names <- c()  # ← stub: replaced by generate_R_from_pnpro()
  
  # initial marking vector (same length as yini.names)
  y_ini <- c()       # ← stub: replaced by generate_R_from_pnpro()
  
  # assign names & enforce PNPRO order
  names(y_ini) <- yini.names
  y_ini <- y_ini[yini.names]
  
  return(y_ini)
}

# --------------------------------------------------------------
# (You can add other functions below that will use init.gen(), 
#  or additional functions for sensitivity, ect.
# --------------------------------------------------------------
