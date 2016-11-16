# Copyright (C) 2016 Lyron Winderbaum
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



#' @describeIn combine_peaklists Reads Peak-list Produced by
#'   \code{combine_peaklists}.
#' @export
load_peaklist <- function(o.name , o.path = "."){
  utils::read.table(file.path(o.path, paste(o.name, "_peaklist.txt",sep="")),
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

#' @describeIn combine_peaklists Reads Spectrum-list Produced by
#' \code{combine_peaklists}.
#' @export
load_speclist <- function(o.name , o.path = "."){
  utils::read.table(file.path(o.path, paste(o.name, "_speclist.txt",sep="")),
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

extract_capture <- function(x, rex, cap){
  as.numeric(substring(x, attr(rex,"capture.start")[,cap], {
    attr(rex,"capture.start")[,cap] + attr(rex,"capture.length")[,cap] - 1}))
}

#' Read and Combine Bruker MSI peaklist files.
#'
#' Peakpicking using Bruker mass spectrometry imaging software generally
#' produces a folder filled with many peaklists --- one for each spectrum.
#' \code{combine_peaklists} reads all such peaklist files in a given folder,
#' summarises and writes the relevant information to two tables.
#' The \code{load_*} functions are for reading the tables created by
#' \code{combine_peaklists}.
#'
#' Typically, Bruker raw data will be stored all in a single folder (often with
#' the ".d" extension), with a folder name that identifies the relevant
#' information about the run. The only compulsory argument to
#' \code{combine_peaklists} is a filepath to that folder, \code{i.path} under
#' arguments. Typically the peaklist files themselves will be in a subfolder of
#' this folder, often called "peaklists" (can be specified via \code{i.name}).
#' This is the only information with which \code{combine_peaklists} is concerned
#' --- any other files in said folder are ignored. The output files will by
#' default be named using the folder name provided in \code{i.path}, but this
#' can be overridden by supplying \code{o.name}. These output files are by
#' default created in the current working directory but this can be changed via
#' \code{o.path}.
#'
#' \code{combine_peaklists} creates and writes relevant information to two
#' tables: a peak-list, and a spectrum-list:
#' \itemize{
#'   \item Peak-list: A table in which each row corresponds to a peak. Consists
#'   of all original peaklists concatenated together with one additional
#'   column, \code{Acq}, identifying the spectrum from which a peak
#'   originated.
#'   \item Spectrum-list: A table in which each row corresponds to a spectrum.
#'   Columns contain information relevant at a spectrum level, specifically
#'   this table contains five columns total:
#'   \itemize{
#'     \item \code{fname}: The original peaklist filename.
#'     \item \code{R}:     The region numbers of the corresponding spectra.
#'     \item \code{X}:     The X-coordinates of the corresponding spectra.
#'     \item \code{Y}:     The Y-coordinates of the corresponding spectra.
#'     \item \code{Acq}:   An integer identifier used to cross-reference to
#'     the peak-list table.
#'   }
#' }
#'
#' Note that \code{Acq} will correspond to the order of acquisition of the
#' spectra if the spectra where acquired in increasing order of first region
#' number, second Y-coordinate and third X-coordinate, in that order of
#' priority. This order or acquisition is a common default on Bruker Flex
#' instruments.
#'
#' @section Assumptions:
#' \code{combine_peaklists} assumes that the data is in a particular format --
#' that produced by Bruker software peakpicking at the time of writing this
#' vignette, specifically that:
#' \itemize{
#'   \item Each spectrum is represented by a single peaklist text-file
#'   containing any peaks detected in that spectrum.
#'   \item All peaklist files (for any given run) are stored in the same folder.
#'   \item The filename of each peaklist file contains a match to the perl
#'   regular expression \code{R(?P<r>\\d{2,3})X(?P<x>\\d{3,4})Y(?P<y>\\d{3,4})}
#'   in which the named capture groups `r`, `x`, and `y`, correspond to the
#'   region number, x-coordinate, and y-coordinate of the corresponding
#'   spectrum.
#'   \item The contents of each peaklist file is formatted such that it can be
#'   correctly read with the use \code{utils::read.table(. , header = TRUE)}.
#'   TODO: Update this to be more specific and instead use \code{base::scan},
#'   and remove the \code{load_*} functions instead replacing them with
#'   suggestions for direct use of a function such as \code{utils::read.csv}.
#' }
#'
#' Note that \code{combine_peaklists} checks for duplicate spectra (peaklist
#' files with the same region number, X-coordinate, and Y-coordinate) and will
#' throw an error if it finds any.
#'
#' @param i.path Path to dataset of interest. This will usually be a folder with
#'   a subfolder (\code{i.name}) that contains the peaklist files.
#' @param i.name Name of the subfolder in \code{i.path} that contains the
#'   peaklist files.
#' @param o.path Path to where output files should be written.
#' @param o.name Name identifying the dataset of interest. If left as
#'   \code{NULL} for \code{combine_peaklists}, this will default to
#'   \code{basename(i.path)}, which is then what should be passed to the
#'   \code{load_*} functions.
#'
#' @return On successful completion returns the number of empty spectra found
#'   --- peaklist files with a header but no peaks. If no peaklist files are
#'   found at all returns -1 and a warning.
#'
#' @seealso \code{\link{basename}}
#'
#' @examples
#' combine_peaklists("A1")
#'
#' @export
combine_peaklists <- function(i.path,
                           i.name = "peaklists",
                           o.path = ".",
                           o.name = NULL){

  if (is.null(o.name)) o.name = basename(i.path)

  ###################################################
  # Find the peaklist files and extract their names #
  ###################################################
  pl.path = file.path(i.path, i.name)
  # Matches files to the regular expression for peaklist files.
  pl.fnames = list.files(path = pl.path, "R\\d{2,3}X\\d{3,4}Y\\d{3,4}.txt")
  if (length(pl.fnames) == 0){
    warning("dipps::combine_peaklists: No peaklist files found, aborting.")
    return(-1)
  }

  # Extracts the Region numbers and X,Y coordinates for all the spectra.
  rex = regexpr("R(?P<r>\\d{2,3})X(?P<x>\\d{3,4})Y(?P<y>\\d{3,4})",
                pl.fnames, perl = TRUE)
  df.spec <- data.frame(fname = pl.fnames,
                        R     = extract_capture(pl.fnames, rex, "r"),
                        X     = extract_capture(pl.fnames, rex, "x"),
                        Y     = extract_capture(pl.fnames, rex, "y"))
  df.spec = df.spec[order(df.spec$R, df.spec$Y, df.spec$X),]
  df.spec$Acq = 1:nrow(df.spec)
  utils::write.table(df.spec,
                     file = file.path(o.path,
                                      paste(o.name, "_speclist.txt", sep="")),
                     sep = "\t", row.names = FALSE, col.names = TRUE)

  pl_not_empty = logical(length=length(pl.fnames))
  header_not_written = TRUE

  # File to write output
  o.file = file(file.path(o.path, paste(o.name, "_peaklist.txt", sep="")), "w")

  # For each peaklist file
  for (spec_idx in 1:length(pl.fnames)){
    fname <- pl.fnames[spec_idx]

    # Read the peaklist file
    pl.cur <- utils::read.table(file.path(pl.path, fname), header=TRUE)
    # Check that it is not empty (no peaks)
    if (nrow(pl.cur) == 0) next
    # It's not empty!
    pl_not_empty[spec_idx] <- TRUE
    # Annotate peaks with the Acq of their parent peaklist.
    pl.cur <- transform(pl.cur, Acq = df.spec[df.spec$fname == fname,]$Acq)
    # Write peaklist
    if (header_not_written) {
      utils::write.table(pl.cur, o.file, sep="\t", row.names=FALSE)
      header_not_written = FALSE
    } else {
      utils::write.table(pl.cur, o.file, sep="\t",
                         row.names=FALSE, col.names=FALSE)
    }
  }
  close(o.file)
  # could save peaklist_not_empty here, if you needed it for anything later on.
  return(sum(!pl_not_empty))
}
