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



#' Difference in ProPortions Statistics (DIPPS)
#'
#' Calculates the DIPPS for the given subset. The argument descriptions are
#' generic as DIPPS can be applied to any binary (``occurrence'') data in which
#' each variable has two values (``occurrence'' and ``absence''). In the MSI
#' context, an occurrence is generally taken to be a peak, an observation is
#' generally taken to be a spectrum and a variable is generally taken to be
#' a mass range or peakgroup, possibly grouped via some clustering method such
#' as that offered by \code{dbscan}.
#'
#' \code{obs}, \code{var}, and \code{subset} must be equal length, and can be
#' taken from the output of \code{combine_peaklists} with relative ease -- see
#' example below. It is also assumed that equal entries in \code{obs} should
#' have equal entries in \code{subset} as well. TODO: I should add a check for
#' that.
#'
#' Note that from the perspective of treating occurrence in each variable
#' (seperately) being used as a binary classifier for membership in the subset,
#' the DIPPS can be thought of as the Informedness of these classifiers, i.e.
#' the DIPPS = sensitivity + specificity - 1.
#'
#' @param obs   A vector identifying the observation from which an occurrence
#'   originated.
#' @param var  A vector identifying the variable of which an occurrence is a
#'   realisation.
#' @param subset A vector identifying occurrences belonging to the subset of
#'   observations of interest.
#'
#' @return Successful completion will return a data.frame in which rows
#' represent variables (as identified by \code{var}), ordered in decreasing
#' order of DIPPS, and with seven columns:
#' \itemize{
#'   \item \code{var}.
#'   \item p.u: proportions of occurrence in the \code{subset == TRUE} subset of
#'     observations.
#'   \item p.d: proportions of occurrence in the \code{subset == FALSE} subset
#'     of observations.
#'   \item d:   p.u - p.d (DIPPS).
#'   \item c.u: the cosine distance centroid of the \code{subset == TRUE} subset
#'     of observations.
#'   \item cos: the cosine distance between c.u and the `template' vector t
#'     which contains ones in each peakgroup with a DIPPS equal to or greater
#'     than the DIPPS of the peakgroup the corresponding row represents.
#'   \item t:   the `template' vector for the heuristically chosen `optimal'
#'     DIPPS cutoff -- i.e. selecting a number of the highest DIPPS variables
#'     such that the cosine distance as described above is minimised, under the
#'     contraint that the dipps cutoff should be positive.
#' }
#'
#' @seealso \code{\link{combine_peaklists}}, \code{\link{dbscan}},
#'
#' Winderbaum, L. J. et al. Feature extraction for proteomics imaging
#' mass spectrometry data. The Annals of Applied Statistics.
#' 2015;9(4):1973-1996. doi: 10.1214/15-AOAS870.
#'
#' @examples
#' i.path = system.file("extdata", "test1", package = "dipps")
#' n.empty = combine_peaklists(i.path)
#' o.name = basename(i.path)
#' df.spec = load_speclist(o.name)
#' df.peak = load_peaklist(o.name)
#'
#' # Construct peakgroups
#' df.peak$group = dbscan(df.peak$m.z, eps = 0.1, mnpts = 1)
#'
#' # Select a subset of spectra expected to be overexpressed. In this case
#' # spectra with Y-coordinate greater than or equal to 170.
#' df.spec$subset = df.spec$Y >= 170
#' df.peak = merge(df.peak, df.spec[, c("Acq", "subset")])
#'
#' # Calculate DIPPS
#' df.dipps = dipps(df.peak$Acq, df.peak$group, df.peak$subset)
#'
#' @export
dipps <- function(obs, var, subset) {
  # Add checks on input

  df.peak = data.frame(obs = obs,
                       var = var,
                       subset = subset)

  # Calculate Proportions of Occurence
  nSpec_d = length(unique(obs[!subset]))
  nSpec_u = length(unique(obs[subset]))
  prop = plyr::ddply(df.peak,
                     c("var", "subset"),
                     plyr::summarise,
                     p = length(obs))
  prop[!prop$subset, "p"] = prop[!prop$subset, "p"]/nSpec_d
  prop[prop$subset,  "p"] = prop[prop$subset,  "p"]/nSpec_u

  # Reshape into vector form
  prop <- stats::reshape(prop,
                         timevar = "subset",
                         idvar = "var",
                         direction="wide")
  prop = replace(prop, is.na(prop), 0)
  names(prop)[names(prop) == "p.FALSE"] = "p.d"
  names(prop)[names(prop) == "p.TRUE"]  = "p.u"

  # Calculate DIPPS
  prop$d = prop$p.u - prop$p.d

  # Calculate subset centroid
  u.m = reshape2::dcast(subset(df.peak, subset),
                        obs ~ var,
                        value.var = "subset")
  acq = u.m[,1]
  u.m = as.matrix(u.m[,-1])
  u.m[!is.na(u.m)] = 1
  u.m[is.na(u.m)] = 0
  acq = data.frame(obs = acq,
                   n.peaks = rowSums(u.m))
  u.m = u.m/sqrt(acq$n.peaks)
  c.u = colMeans(u.m)
  c.u = c.u/norm(c.u,"2")

  tmp = match(colnames(u.m), prop$var)
  prop$c.u = 0
  prop[tmp,"c.u"] = c.u

  # Calculate Heuristic cutoff.
  p = nrow(prop)
  curMinCosD = 10
  prop = prop[order(prop$d, decreasing = TRUE), ]
  vN = 1:p
  prop$cos = NA
  for (n in vN){
    prop$t = 0
    prop[1:n, "t"] = 1/sqrt(n)
    prop[n, "cos"] = 1 - sum(prop$t * prop$c.u)
  }
  prop$t = 0
  prop[1:which.min(prop[prop$d > 0, "cos"]), "t"] = 1

  return(prop)
}










