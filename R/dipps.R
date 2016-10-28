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

#' Difference in ProPortions Statistic (DIPPS)
#'
#' Calculates the DIPPS for the given subset.
#'
#' @param df.peak Parameter
#' @param id.var Integer variable in \code{df.peak} that identifies unique
#'   spectra.
#' @param group.var Integer variable in \code{df.peak} that identifies peaks
#'   assigned to the same `peakgroup' -- typically grouped by m.z clustering.
#' @param subset.var Boolean variable in \code{df.peak} in which values of
#'   TRUE identify the subset of spectra expected to be `upregulated'.
#'
#' @seealso Winderbaum, L. J. et al. Feature extraction for proteomics imaging
#'   mass spectrometry data. The Annals of Applied Statistics.
#'   2015;9(4):1973-1996. doi: 10.1214/15-AOAS870.")
#' @export
dipps <- function(df.peak,
                  id.var = "Acq",
                  group.var = "group",
                  subset.var = "subset") {
  # Add checks on input

  # Calculate Proportions of Occurence
  nSpec_d = length(unique(df.peak[!df.peak[, subset.var], id.var]))
  nSpec_u = length(unique(df.peak[df.peak[, subset.var],  id.var]))
  prop = plyr::ddply(df.peak,
                     c(group.var, subset.var),
                     plyr::summarise,
                     p = length(eval(as.symbol(id.var))))
  prop[!prop[, subset.var], "p"] = prop[!prop[, subset.var], "p"]/nSpec_d
  prop[prop[, subset.var],  "p"] = prop[prop[, subset.var],  "p"]/nSpec_u

  # Reshape into vector form
  prop <- stats::reshape(prop,
                         timevar = subset.var,
                         idvar = group.var,
                         direction="wide")
  names(prop) = c(group.var,"p.d","p.u")

  prop = replace(prop, is.na(prop),0)
  # Calculate DIPPS
  prop$d = prop$p.u - prop$p.d

  return(prop)
}




#' Heuristic DIPPS cutoff
#'
#' Calculates a heuristic for the optimal DIPPS cutoff.
#'
#' @inheritParams dipps
#' @param prop.m The output from a call to \code{dipps}.
#'
#' @seealso Winderbaum, L. J. et al. Feature extraction for proteomics imaging
#'   mass spectrometry data. The Annals of Applied Statistics.
#'   2015;9(4):1973-1996. doi: 10.1214/15-AOAS870.")
#' @export
dipps_cutoff <- function(df.peak,  prop.m,
                      id.var = "Acq",
                      group.var = "group",
                      subset.var = "subset") {
  # Add checks for inputs.

  u.m = reshape2::dcast(subset(df.peak, eval(as.symbol(subset.var))),
                        eval(as.symbol(id.var)) ~ eval(as.symbol(group.var)),
                        value.var = subset.var)
  acq = u.m[,1]
  u.m = as.matrix(u.m[,-1])
  u.m[!is.na(u.m)] = 1
  u.m[is.na(u.m)] = 0
  acq = data.frame(Acq = acq,
                   n.peaks = rowSums(u.m))
  u.m = u.m/sqrt(acq$n.peaks)
  c.u = colMeans(u.m)
  c.u = c.u/norm(c.u,"2")

  tmp = match(colnames(u.m), prop.m[, group.var])
  dsum_in$c.u = 0
  dsum_in[tmp,"c.u"] = c.u

  # Find data-driven cutoff according to the DIPPS method using cosine distance.
  p = nrow(prop.m)
  curMinCosD = 10
  sortedDIPPS = sort(prop.m$d, index.return=TRUE)
  # vN = 1:floor(nrow(Summary_merged)/2)
  vN = 1:nrow(prop.m)
  cosD = vN
  for (n in vN){
    prop.m$t = 0
    prop.m[sortedDIPPS[(p - n + 1):p, "ix"], "t"] = 1/sqrt(n)
    cosD[n] = 1 - sum(prop.m$t * prop.m$c.u)
  }
  nStar = vN[which.min(cosD)]

  return(nStar)
}






