#' Read vibrational data from Gaussian 09 output
#'
#' This function reads vibrational data from a Gaussian 09 (and potentially
#' later) output file.
#'
#' \code{readgaussout} reads vibrational data from a Gaussian 09 (and potentially
#' later) output file. It extracts the frequencies, IR intensities, and Raman
#' activities from the file and returns them as a matrix.
#'
#' @param Path Name (with path) of the file from which the data are to be read.
#' @param logcont The content of a .log file as read in by readLines(). This
#'   parameter can be used to combine multiple readgaussout function calls by
#'   avoiding re-reading the file in each function. If logcont is provided, the
#'   function uses it as the content of the file instead of reading from Path.
#'
#' @return A matrix with three rows containing the vibrational frequency/band
#'   position in cm-1 (\code{BPosi}), the IR intensity (\code{IRint}) as well as
#'   the Raman activity (\code{RActi}) for each vibration.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgaussout("Path/To/File.out")
#' }
#'
#' @export
#' @importFrom gsubfn read.pattern
readgaussout <- function(Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- readLines(Path)
  } else {
    tmp <- logcont
  }

  BPosi <- stringr::str_match_all(tmp, "Frequencies -- +(-?\\d+.\\d+) *(-?\\d+.\\d+) *(-?\\d+.\\d+)")
  IRInt <- stringr::str_match_all(tmp, "IR Inten    -- +(-?\\d+.\\d+) *(-?\\d+.\\d+) *(-?\\d+.\\d+)")
  RActi <- stringr::str_match_all(tmp, "Raman Activ -- +(-?\\d+.\\d+) *(-?\\d+.\\d+) *(-?\\d+.\\d+)")
  
  BPosi <- as.numeric(unlist(lapply(BPosi, "[", 2:4))[!is.na(unlist(lapply(BPosi, "[", 2:4)))])
  IRInt <- as.numeric(unlist(lapply(IRInt, "[", 2:4))[!is.na(unlist(lapply(IRInt, "[", 2:4)))])
  RActi <- as.numeric(unlist(lapply(RActi, "[", 2:4))[!is.na(unlist(lapply(RActi, "[", 2:4)))])
  
  result <- matrix(c(BPosi, IRInt, RActi), 3, length(BPosi), byrow = TRUE)
  rownames(result) <- c("BPosi", "IRint", "RActi")
  
  return(result)
}

#' Read geometry and Mullikan charges from Gaussian 16 output
#' 
#' This function reads geometric and Mullikan charges data from a Gaussian 16
#' (and potentially earlier) output file.
#' 
#' \code{readgaussgeom} is designed to extract geometric and Mullikan charges data from
#' a Gaussian 16 output file, and potentially earlier versions. Gaussian 16 is a
#' widely used quantum chemistry software package for electronic structure
#' calculations.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' @param logcont The content of a .log file as read in by readLines(). May be
#'   used to combine multiple readgauss functions by avoiding to re-read the
#'   file in each function.
#'   
#' @return A data frame containing the extracted information from the Gaussian
#'   output file. The data frame includes the following columns:
#'   \describe{
#'   \item{Num}{The atom number.}
#'   \item{AtomNum}{The atomic number of the atom.}
#'   \item{Sym}{The atomic symbol of the atom.}
#'   \item{Lab}{Unique atom label assigned to each atom.}
#'   \item{Type}{The atom type.}
#'   \item{X}{The X-coordinate of the atom in the Cartesian coordinate system.}
#'   \item{Y}{The Y-coordinate of the atom in the Cartesian coordinate system.}
#'   \item{Z}{The Z-coordinate of the atom in the Cartesian coordinate system.}
#'   \item{MullikanCharge}{The Mullikan charge associated with each atom.}
#'   }
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgaussgeom("Path/To/File.out")
#' }
#' 
#' @export
#' @importFrom PeriodicTable symb
readgaussgeom <- function(Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- readLines(Path)
  } else {
    tmp <- logcont
  }
  
  lastgeom <- max(which(tmp == "                         Standard orientation:                         " | tmp == "                          Input orientation:                          ")) ## Standard orientation and input orientation are identical, the first is a rotation of the latter
  breaks <- which(tmp == " ---------------------------------------------------------------------")
  lastgeomend <- breaks[which(breaks == lastgeom+4)+1]-1
  optgeom <- paste(tmp[(lastgeom+5):(lastgeomend)], collapse = "")
  
  dat <- as.numeric(strsplit(optgeom, "\\s+")[[1]][-1])
  datmat <- matrix(dat, length(dat)/6, 6, byrow = TRUE)
  
  mulli <- max(which(tmp == " Mulliken charges:"))
  mullivec <- paste(tmp[(mulli+2):(mulli+1+nrow(datmat))], collapse = "")
  mullicharge <- as.numeric(grep("-?\\d+\\.\\d*", strsplit(mullivec, "\\s+")[[1]][-1],
                                 perl = TRUE, value = TRUE))
  
  unilab <- uniqueAtomLabels(datmat[, 2])
  
  datdf <- data.frame(datmat[, 1:2], PeriodicTable::symb(datmat[, 2]), unilab,
                      datmat[, 3:6], mullicharge)
  colnames(datdf) <- c("Num", "AtomNum", "Sym", "Lab", "Type",
                       "X", "Y", "Z","MullikanCharge")
  
  return(datdf)
}

#' Read energies from Gaussian 16 output
#' 
#' This function reads energy data from a Gaussian 16 (and potentially earlier)
#' output file.
#' 
#' \code{readgaussener} is designed to read energy data from a Gaussian 16 output file.
#' It can potentially work with earlier versions of Gaussian as well. The
#' function returns the Zero-point Energy, the Energy, the Enthalpies as well as
#' the Free Energies in kJ/mol.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' @param logcont The content of a .log file as read in by readLines(). May be
#'   used to combine multiple readgauss functions by avoiding to re-read the
#'   file in each function.
#'   
#' @return A numeric vector containing the energy values in kilojoules per mole
#'   (kJ/mol).
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgaussener("Path/To/File.out")
#' }
#' 
#' @export
readgaussener <- function(Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- paste(readLines(Path), collapse = "")
  } else {
    tmp <- paste(logcont, collapse = "")
  }

  HatreeEner <- c(stringr::str_match(tmp, "Sum of electronic and zero-point Energies= +(-?\\d+\\.\\d+)")[, 2],
                  stringr::str_match(tmp, "Sum of electronic and thermal Energies= +(-?\\d+\\.\\d+)")[, 2],
                  stringr::str_match(tmp, "Sum of electronic and thermal Enthalpies= +(-?\\d+\\.\\d+)")[, 2],
                  stringr::str_match(tmp, "Sum of electronic and thermal Free Energies= +(-?\\d+\\.\\d+)")[, 2]
  )
  
  kJmolEner <- as.numeric(HatreeEner) * 2625.438
  
  return(kJmolEner)
}

#' Read calculation time from Gaussian 16 output
#' 
#' This function reads the calculation time from a Gaussian 16 (and potentially
#' earlier) output file.
#' 
#' readgausscalctime is designed to read the calculation time from a Gaussian 16
#' output file. It can potentially work with earlier versions of Gaussian as
#' well. The function returns the Job CPU Time in s.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' @param logcont The content of a .log file as read in by readLines(). May be
#'   used to combine multiple readgauss functions by avoiding to re-read the
#'   file in each function.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgausscalctime("Path/To/File.out")
#' }
#'   
#' @return A numeric vector containing the calculation time in seconds (s).
#' 
#' @export
readgausscalctime <- function (Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- paste(readLines(Path), collapse = "")
  } else {
    tmp <- paste(logcont, collapse = "")
  }

  tmp <- paste(readLines(Path), collapse = "")
  dhms <- stringr::str_match_all(tmp, "Job cpu time: +(\\d+) days +(\\d+) hours +(\\d+) minutes +(\\d+\\.\\d+) seconds.")[[1]][, 2:5]
  TimeInS <- sum(as.numeric(t(dhms)) * c(24*60*60, 60*60, 60, 1))
  return(TimeInS)
}

#' Read atom motion during vibration from Gaussian 16 output
#' 
#' This function reads the atom motions from a Gaussian 16 (and potentially
#' earlier) output file.
#' 
#' \code{readgaussvibmotion} is designed to read the calculated motions of individual
#' atoms during vibrations from a Gaussian 16 output file. The function can
#' potentially work with earlier versions of Gaussian as well. The function
#' assumes 3N-6 vibrations, where N is the number of atoms. It might give
#' strange results for linear molecules with 3N-5 vibrations.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' @param logcont The content of a .log file as read in by readLines(). May be
#'   used to combine multiple readgauss functions by avoiding to re-read the
#'   file in each function.
#'   
#' @return A data frame containing the atom motion during vibrations.
#'   \describe{
#'     \item{Num}{The atom number.}
#'     \item{AtomNum}{The atomic number of the atom.}
#'     \item{Sym}{The atomic symbol of the atom.}
#'     \item{Lab}{Unique atom label assigned to each atom.}
#'     \item{Vib}{Vibration index.}
#'     \item{dX}{Displacement in the X direction.}
#'     \item{dY}{Displacement in the Y direction.}
#'     \item{dZ}{Displacement in the Z direction.}
#'   }
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgaussvibmotion("Path/To/File.out")
#' }
#' 
#' @export
#' @importFrom PeriodicTable symb
#' @importFrom dplyr arrange
#' @importFrom stringr str_match str_detect str_split
#' @importFrom PeriodicTable symb
#' @importFrom rlang .data
#' @importFrom stats complete.cases
readgaussvibmotion <- function(Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- readLines(Path)
  } else {
    tmp <- logcont
  }
  
  heads <- which(str_detect(tmp, "  Atom  AN      X      Y      Z.*"))
  natoms <- as.numeric(str_match(tmp[which(
    !is.na(str_match(tmp, " NAtoms= *(\\d+) NQM")))],
    " NAtoms= *(\\d+) NQM")[1, 2])
  
  temp <- tmp[as.vector(t(outer(heads, (1:natoms), "+")))]
  temp2 <- lapply(str_split(temp, "\\s+"), function(x) as.numeric(x[-1]))
  
  datmat <- do.call(rbind, lapply(temp2,
                                  function(x) {
                                    length(x) <- 11
                                    rbind(x[1:5], x[c(1:2, 6:8)], x[c(1:2, 9:11)])
                                  }
  ))
  
  unilab <- rep(rep(uniqueAtomLabels(datmat[seq(1, 3*natoms, 3), 2]), each = 3), (3*natoms - 6) / 3)
  
  datmat <- datmat[complete.cases(datmat), ]
  
  datdf <- data.frame(Num = datmat[, 1],
                      AtomNum = datmat[, 2],
                      Sym = symb(datmat[, 2]),
                      Lab = unilab,
                      Vib = rep(1:3, nrow(datmat) / 3) + rep(seq(0, (3 * natoms - 6) - 1, 3), each = 3 * natoms),
                      dX = datmat[, 3],
                      dY = datmat[, 4],
                      dZ = datmat[, 5])
  
  datdf <- arrange(datdf, .data$Vib, .data$Num)
  rownames(datdf) <- NULL
  
  return(datdf)
}

#' Read dipole moment from Gaussian 16 output
#' 
#' This function reads the total dipole moment as well as the X, Y and Z
#' components of the dipole moment from a Gaussian 16 (and potentially earlier)
#' output file.
#' 
#' \code{readgaussdipole} is designed to read the total dipole moment as well as the
#' X, Y and Z components of the dipole moment from a Gaussian 16
#' output file. It can potentially work with earlier versions of Gaussian as
#' well. The function returns the dipole moments in Debye.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' @param logcont The content of a .log file as read in by readLines(). May be
#'   used to combine multiple readgauss functions by avoiding to re-read the
#'   file in each function.
#'   
#' @return A numeric vector containing the dipole moments in Debye (D).
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgaussdipole("Path/To/File.out")
#' }
#' 
#' @export
#' @importFrom stringr str_match
readgaussdipole <- function(Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- readLines(Path)
  } else {
    tmp <- logcont
  }

    lastdipole <- max(which(tmp == " Dipole moment (field-independent basis, Debye):"))
    dipoleline <- tmp[lastdipole + 1]
    dipoles <- as.numeric(str_match(dipoleline, "\\s*X=\\s*(-?\\d+\\.\\d+)\\s*Y=\\s*(-?\\d+\\.\\d+)\\s*Z=\\s*(-?\\d+\\.\\d+)\\s*Tot=\\s*(-?\\d+\\.\\d+)")[, 2:5])
    return(dipoles)
}

#' Read frontier orbital energies from Gaussian 16 output
#' 
#' This function reads the frontier orbital energies of the Highest Occupied
#' Molecular Orbital (HOMO) and the Lowest Unoccupied Molecular Orbital (LUMO)
#' from a Gaussian 16 (and potentially earlier) output file.
#' 
#' \code{readgaussfroorb} is designed to read the frontier orbital energies of the
#' Highest Occupied Molecular Orbital (HOMO) and the Lowest Unoccupied Molecular
#' Orbital (LUMO) from a Gaussian 16. It can potentially work with earlier
#' versions of Gaussian as well. It also calculates the HOMO-LUMO energy
#' difference. The function returns the dipole moments in in kJ/mol and assumes
#' at max 5 orbitals per line in the output file.
#' 
#' @param Path Name (with path) of the file which the data are to be read from.
#' @param logcont The content of a .log file as read in by readLines(). May be
#'   used to combine multiple readgauss functions by avoiding to re-read the
#'   file in each function.
#'   
#' @return A numeric vector containing the energy values in kilojoules per mole
#'   (kJ/mol).
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   QMout <- readgaussfroorb("Path/To/File.out")
#' }
#' 
#' @export
#' @importFrom stringr str_match
readgaussfroorb <- function(Path, logcont = NULL) {

  if(is.null(logcont)) {
    tmp <- readLines(Path)
  } else {
    tmp <- logcont
  }

    homopos <- max(grep(" Alpha  occ. eigenvalues --", tmp))
    homoline <- tmp[homopos]
    homo <- as.numeric(str_match(homoline, " Alpha  occ\\. eigenvalues -- +(-?\\d+\\.\\d+) *(-?\\d+\\.\\d+)? *(-?\\d+\\.\\d+)? *(-?\\d+\\.\\d+)? *(-?\\d+\\.\\d+)?")[, 2:6])
    homo <- max(homo[!is.na(homo)])
    lumoline <- tmp[homopos+1]
    lumo <- as.numeric(str_match(lumoline, " Alpha virt\\. eigenvalues -- +(-?\\d+\\.\\d+) *(-?\\d+\\.\\d+)? *(-?\\d+\\.\\d+)? *(-?\\d+\\.\\d+)? *(-?\\d+\\.\\d+)?")[, 2])
    froorbener <- c(homo, lumo, lumo-homo) * 2625.438
    return(froorbener)
}