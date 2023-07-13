#' Read vibrational data from ORCA 5.0 output
#'
#' This function reads vibrational data from a ORCA 5.0 (and potentially
#' later) output file.
#'
#' \code{readorcaout} reads vibrational data from a ORCA 5.0 (and potentially
#' later) output file. It extracts the frequencies, IR intensities, and Raman
#' activities from the file and returns them as a matrix.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' 
#' @return A matrix with three rows containing the vibrational frequency/band
#'   position in cm-1 (\code{BPosi}), the IR intensity (\code{IRint}) as well as
#'   the Raman activity (\code{RActi}) for each vibration.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readorcaout("Path/To/File.out")
#' }
#' 
#' @export
#' @importFrom gsubfn read.pattern
#' @importFrom stringr str_match str_split
readorcaout <- function(Path) {
  tmp <- readLines(Path)

  IRstart <- which(!is.na(as.vector(stringr::str_match(tmp, "IR SPECTRUM")))) + 6
  RAstart <- which(!is.na(as.vector(stringr::str_match(tmp, "RAMAN SPECTRUM  ")))) + 5
  IRend <- which(!is.na(as.vector(stringr::str_match(tmp, "\\* The epsilon \\(eps\\) is given for a Dirac delta lineshape\\.")))) - 2
  RAend <- which(!is.na(as.vector(stringr::str_match(tmp, "The first frequency considered to be a vibration is")))) - 2

  IRtmp <- do.call("rbind", lapply(stringr::str_split(trimws(tmp[IRstart:IRend]), "\\s+"), function(x) {as.numeric(x[c(2, 4)])}))
  RAtmp <- do.call("rbind", lapply(stringr::str_split(trimws(tmp[RAstart:RAend]), "\\s+"), function(x) {as.numeric(x[2:3])}))

  result <- t(cbind(IRtmp, RAtmp[, 2]))
  rownames(result) <- c("BPosi", "IRint", "RActi")
  return(result)
}

#' Read energies from ORCA 5.0 output
#' 
#' This function reads energy data from a ORCA 5.0 (and potentially later)
#' output file.
#' 
#' \code{readorcaener} is designed to read energy data from a ORCA 5.0 output file.
#' It can potentially work with earlier and later versions of ORCA as well. The
#' function returns the Electronic Energy, the Total Thermal Energy, the Total
#' Enthalpy as well as the Final Gibbs Free Energies in kJ/mol.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' 
#' @return A numeric vector containing the energy values in kilojoules per mole
#'   (kJ/mol).
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readorcaener("Path/To/File.out")
#' }
#' 
#' @export
readorcaener <- function(Path) {
  tmp <- paste(readLines(Path), collapse = "")
  
  HatreeEner <- c(stringr::str_match(tmp, "Electronic energy \\s+ ... \\s* +(-?\\d+\\.\\d+)")[, 2],
                  stringr::str_match(tmp, "Total thermal energy \\s+ +(-?\\d+\\.\\d+)")[, 2],
                  stringr::str_match(tmp, "Total Enthalpy \\s+ ... \\s* +(-?\\d+\\.\\d+)")[, 2],
                  stringr::str_match(tmp, "Final Gibbs free energy \\s+ ... \\s* +(-?\\d+\\.\\d+)")[, 2]
  )
  
  kJmolEner <- as.numeric(HatreeEner) * 2625.438
  
  return(kJmolEner)
}

#' Read calculation time from ORCA 5.0 output
#' 
#' This function reads the calculation time from a ORCA 5.0 (and potentially
#' earlier) output file.
#' 
#' \code{readorcacalctime} is designed to read the calculation time from a ORCA 5.0
#' output file. It can potentially work with earlier and later versions of ORCA
#' as well. The function returns the Total Run Time in s.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' 
#' @return A numeric vector containing the calculation time in seconds (s).
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readorcacalctime("Path/To/File.out")
#' }
#' 
#' @export
readorcacalctime <- function (Path) {
  tmp <- paste(readLines(Path), collapse = "")
  dhms <- stringr::str_match_all(tmp, "TOTAL RUN TIME: +(\\d+) days +(\\d+) hours +(\\d+) minutes +(\\d+) seconds +(\\d+) msec")[[1]][, 2:6]
  TimeInS <- sum(as.numeric(t(dhms)) * c(24*60*60, 60*60, 60, 1, 0.001))
  return(TimeInS)
}

