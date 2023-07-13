#' Read xyz file
#'
#' This function reads data from an .xyz file and returns a data frame with the
#' parsed information.
#' 
#' This function reads data from an .xyz file and returns a data frame with the
#' parsed information. .xyz files are commonly used to represent molecular
#' structures, where each line contains the atomic symbol followed by its
#' Cartesian coordinates.
#'
#' @param path The name (with path) of the file to be read.
#' @param skip The number of lines to skip during the read process. Default is
#'   2, assuming the file contains a header and the number of atoms.
#' @param ... Additional arguments passed to the underlying read.table function.
#' 
#' @return \code{readxyz} constructs a \code{data.frame} with seven columns for
#'   the number of each atom (\code{Num}), atomic numbers (\code{AtomNum}), atom
#'   symbols (\code{Sym}), unique labels (\code{Lab}), and the Cartesian
#'   coordinates (\code{X}, \code{Y}, \code{Z}) from the read data frame. The
#'   column names are set accordingly.
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   XYZout <- readxyz("Path/To/File.xyz")
#' }
#'
#' @importFrom utils read.table
#' @export
readxyz <- function(path, skip = 2, ...) {
  periodicTable <- get("periodicTable")
  
  temp <- read.table(path, skip = skip, ...)
  
  num <- 1:nrow(temp)
  unilab <- uniqueAtomLabels(temp[, 1])
  atomnum <- match(temp[, 1], periodicTable$symb) - 1
  datdf <- data.frame(num, atomnum, temp[, 1], unilab, temp[, 2:4])
  colnames(datdf) <- c("Num", "AtomNum", "Sym", "Lab", "X", "Y", "Z")
  
  return(datdf)
}