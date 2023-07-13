#' Create unique alphanumeric atom labels
#' 
#' This function takes a vector of atom labels and returns a vector with a
#' unique alphanumeric label for each atom.
#' 
#' This function takes a vector of atom labels and returns a vector with unique
#' alphanumeric labels. If the input vector contains numeric labels, the
#' function converts them to their corresponding symbol using the
#' \code{\link[PeriodicTable]{symb}} function.
#' 
#' @param label Vector containing the atom labels
#' 
#' @examples
#'   examlab <- c("C", "H", "H", "H", "H")
#'   unilab <- uniqueAtomLabels(examlab)
#' 
#' @export
#' @importFrom PeriodicTable symb
uniqueAtomLabels <- function(label) {
  
  if(is.numeric(label)) {label <- PeriodicTable::symb(label)}
  
  dummy <- label
  unilab <- unique(label)
  
  for(i in unilab) {
    temp <- which(label == i)
    
    dummy[temp] <- paste(i, 1:length(temp), sep = "")
  }
  
  return(dummy)
}