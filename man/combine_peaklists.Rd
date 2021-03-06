% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load_peaklist}
\alias{combine_peaklists}
\alias{load_peaklist}
\alias{load_speclist}
\title{Read and Combine Bruker MSI peaklist files.}
\usage{
load_peaklist(o.name, o.path = ".")

load_speclist(o.name, o.path = ".")

combine_peaklists(i.path, i.name = "peaklists", o.path = ".",
  o.name = NULL)
}
\arguments{
\item{o.name}{Name identifying the dataset of interest. If left as
\code{NULL} for \code{combine_peaklists}, this will default to
\code{basename(i.path)}, which is then what should be passed to the
\code{load_*} functions.}

\item{o.path}{Path to where output files should be written.}

\item{i.path}{Path to dataset of interest. This will usually be a folder with
a subfolder (\code{i.name}) that contains the peaklist files.}

\item{i.name}{Name of the subfolder in \code{i.path} that contains the
peaklist files.}
}
\value{
On successful completion returns the number of empty spectra found
  --- peaklist files with a header but no peaks. If no peaklist files are
  found at all returns -1 and a warning.
}
\description{
Peakpicking using Bruker mass spectrometry imaging software generally
produces a folder filled with many peaklists --- one for each spectrum.
\code{combine_peaklists} reads all such peaklist files in a given folder,
summarises and writes the relevant information to two tables.
The \code{load_*} functions are for reading the tables created by
\code{combine_peaklists}.
}
\details{
Typically, Bruker raw data will be stored all in a single folder (often with
the ".d" extension), with a folder name that identifies the relevant
information about the run. The only compulsory argument to
\code{combine_peaklists} is a filepath to that folder, \code{i.path} under
arguments. Typically the peaklist files themselves will be in a subfolder of
this folder, often called "peaklists" (can be specified via \code{i.name}).
This is the only information with which \code{combine_peaklists} is concerned
--- any other files in said folder are ignored. The output files will by
default be named using the folder name provided in \code{i.path}, but this
can be overridden by supplying \code{o.name}. These output files are by
default created in the current working directory but this can be changed via
\code{o.path}.

\code{combine_peaklists} creates and writes relevant information to two
tables: a peak-list, and a spectrum-list:
\itemize{
  \item Peak-list: A table in which each row corresponds to a peak. Consists
  of all original peaklists concatenated together with one additional
  column, \code{Acq}, identifying the spectrum from which a peak
  originated.
  \item Spectrum-list: A table in which each row corresponds to a spectrum.
  Columns contain information relevant at a spectrum level, specifically
  this table contains five columns total:
  \itemize{
    \item \code{fname}: The original peaklist filename.
    \item \code{R}:     The region numbers of the corresponding spectra.
    \item \code{X}:     The X-coordinates of the corresponding spectra.
    \item \code{Y}:     The Y-coordinates of the corresponding spectra.
    \item \code{Acq}:   An integer identifier used to cross-reference to
    the peak-list table.
  }
}

Note that \code{Acq} will correspond to the order of acquisition of the
spectra if the spectra where acquired in increasing order of first region
number, second Y-coordinate and third X-coordinate, in that order of
priority. This order or acquisition is a common default on Bruker Flex
instruments.
}
\section{Functions}{
\itemize{
\item \code{load_peaklist}: Reads Peak-list Produced by
\code{combine_peaklists}.

\item \code{load_speclist}: Reads Spectrum-list Produced by
\code{combine_peaklists}.
}}
\section{Assumptions}{

\code{combine_peaklists} assumes that the data is in a particular format --
that produced by Bruker software peakpicking at the time of writing this
vignette, specifically that:
\itemize{
  \item Each spectrum is represented by a single peaklist text-file
  containing any peaks detected in that spectrum.
  \item All peaklist files (for any given run) are stored in the same folder.
  \item The filename of each peaklist file contains a match to the perl
  regular expression \code{R(?P<r>\\d{2,3})X(?P<x>\\d{3,4})Y(?P<y>\\d{3,4})}
  in which the named capture groups `r`, `x`, and `y`, correspond to the
  region number, x-coordinate, and y-coordinate of the corresponding
  spectrum.
  \item The contents of each peaklist file is formatted such that it can be
  correctly read with the use \code{utils::read.table(. , header = TRUE)}.
  TODO: Update this to be more specific and instead use \code{base::scan},
  and remove the \code{load_*} functions instead replacing them with
  suggestions for direct use of a function such as \code{utils::read.csv}.
}

Note that \code{combine_peaklists} checks for duplicate spectra (peaklist
files with the same region number, X-coordinate, and Y-coordinate) and will
throw an error if it finds any.
}
\examples{
combine_peaklists("A1")

}
\seealso{
\code{\link{basename}}
}

