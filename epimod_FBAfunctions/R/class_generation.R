#' @title Flux Balance Analysis model generation for GreatMod
#' @description ....
#'
#' @param fba_mat Matlab file .....
#' @param fba_model List representing the FBA model generated from *FBAfile.generation*.
#' @param S  Stoic matrix
#' @param ub Vector of the upper constraint of the fluxes
#' @param lb Vector of the lower constraint of the fluxes
#' @param obj_fun Objective function to maximize. Vector of zeros with length equal to the number of reactions,
#'  with one 1 in the position of the flux to maximize.
#' @param react_name The reactions name
#'
#' @details
#'
#' Problem: S*x = b -. size x = n and size of S = m X n
#' With this function the file to pass to the GLPKsolver is generated:
#' 1st row) n_row n_col GLP_MAX (or GLP_MIN if the obj has to max or min)
#' 2nd row) the coeff for defining the obj function (the number of coeff has to be == length of x)
#' m rows)   definition of the S row bounds (i.e. b)
#' n rows)   definition of the x bounds
#' m*n rows) definition of the S coeffs: row_index col_index value
#'
#' @author Aucello Riccardo, Beccuti Marco, Cordero Francesca, Pernice Simone
#'
#' @export

FBA4Greatmod.generation = function(fba_mat=NULL,
                                   fba_model = NULL,
                                   S = NULL, 
                                   ub = NULL,
                                   lb = NULL, 
                                   obj_fun = NULL,
                                   react_name = NULL
                                   ) 
                                   {
  #### Checking the existence of all elements in fba_mat file
  if( !is.null(fba_mat) && file.exists(fba_mat) ) 
    {
    
    modelMAT = FBAmat.read(fba_mat)
    model = FBA_greatmod(as.matrix(modelMAT$S),
                         modelMAT$uppbnd,
                         modelMAT$lowbnd,
                         c(modelMAT$obj_coef),
                         modelMAT$react_id,
                         modelMAT$met_id,
                         gene_assoc = modelMAT$gene_assoc 
												)
  }
  else if(!is.null(S) && !is.null(lb) && !is.null(ub) && !is.null(obj_fun) )
  {
    model = FBA_greatmod(S, ub, lb, obj_fun, react_name)
  }
  else if(!is.null(fba_model))
  {
    model <- fba_model
  }
  else
  {
    stop("The input paramenters are incorrect or missing.")
  }
  validityGreatModClass(model)
  return(model)
}

