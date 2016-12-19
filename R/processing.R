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



#' Extract known masses from a peak-list.
#'
#' Typically used on a peak-list as produced by \code{combine_peaklists} to
#' extract peaks within a given tolerance of a set of known masses.
#'
#' @param df.peak A data.frame with column named \code{var}.
#' @param masses  A numeric vector of masses to extract.
#' @param margin  A positive number: the tolerance around \code{masses} to
#'   extract.
#' @param use_ppm A logical value indicating if \code{margin} is given in ppm,
#'   if not \code{margin} is assumed to be in the same units as \code{var}.
#' @param var     A character string representing the variable/ column in
#'   \code{df.peak} that contain the values to be subsetted.
#'
#' @seealso \code{\link{combine_peaklists}}
#'
#' @export
extract_masses <- function(df.peak,
                           masses,
                           margin = 0.3,
                           use_ppm = FALSE,
                           var = "m.z") {
  # TODO: Check inputs

  df.sub.dne = TRUE
  if (use_ppm) margin.ppm = margin
  for (i in 1:length(masses)){
    if (use_ppm) margin = margin.ppm*masses[i]*1e-6
    l = abs(df.peak[,var] - masses[i]) <= margin
    if (all(!l)) next
    if (df.sub.dne){
      df.sub <- transform(df.peak[l,], group = masses[i])
      df.sub.dne = FALSE
    } else {
      df.sub <- rbind(df.sub, transform(df.peak[l,], group = masses[i]))
    }
  }
  return(df.sub)
}


#' Density Based Clustering
#'
#' This density based clustering is a univariate optimisation of the DBSCAN*
#' algorithm as in section 3 of Campello et al. (2013).
#'
#' Note that a \code{mnpts} of 1 will produce groups that correspond to the
#' equivalence classes of the relation ``are within \code{eps} of each other''.
#' This can be useful for low-density or low-noise data.
#'
#' This is included in this package due to its usefulness in clustering peaks
#' by mass, see example below.
#'
#' @param x     A numeric vector of values to be clustered
#' @param eps   A positive number representing the window radius in which to
#'   group points / measure density.
#' @param mnpts A positive number representing the minimum density for a point
#'   to be included -- the minimum number of points within +/- \code{eps}.
#' @param pp    A logical value indicating whether progress should be printed to
#'   the console.
#'
#' @seealso \code{\link{combine_peaklists}}
#'
#' Section 3 of
#' Campello, Ricardo JGB, Davoud Moulavi, and Joerg Sander.
#' "Density-based clustering based on hierarchical density estimates." In
#' Advances in Knowledge Discovery and Data Mining, pp. 160-172. Springer Berlin
#' Heidelberg, 2013.
#'
#' @examples
#' i.path = system.file("extdata", "test1", package = "dipps")
#' n.emp  = combine_peaklists(i.path)
#' o.name = basename(i.path)
#' df.peak = load_peaklist(o.name)
#' df.peak$group = dbscan(df.peak$m.z)
#'
#' @export
dbscan <- function(x, eps = 0.05, mnpts = 100, pp = FALSE){

  if (!is.vector(x, mode="numeric")){
    stop(paste("dipps::dbscan: invalid input: x is not a numeric vector."))
  }

  if (!is.vector(eps, mode="numeric")){
    stop(paste("dipps::dbscan: invalid input: eps is not a numeric vector."))
  } else if (length(eps) != 1) {
    stop(paste("dipps::dbscan: invalid input: multiple eps values provided."))
  } else if (eps <= 0) {
    stop(paste("dipps::dbscan: invalid input: non-positive eps provided."))
  }
  
  if (!is.vector(mnpts, mode="numeric")){
    stop(paste("dipps::dbscan: invalid input: mnpts is not a numeric vector."))
  } else if (length(mnpts) != 1) {
    stop(paste("dipps::dbscan: invalid input: multiple mnpts values provided."))
  } else if (mnpts <= 0) {
    stop(paste("dipps::dbscan: invalid input: non-positive mnpts provided."))
  }
  
  if (!is.vector(pp, mode="logical")){
    stop(paste("dipps::dbscan: invalid input: pp is not a logical vector."))
  } else if (length(pp) != 1) {
    stop(paste("dipps::dbscan: invalid input: multiple pp values provided."))
  } 
  
  n = length(x)
  if (mnpts > n){
    warning("dipps::dbscan: No core points.")
    return(rep(0, n))
  }

  # sort
  if (pp) {
    print("Sorting...")
  }
  o = order(x)
  x = x[o]
  o = match(1:n, o)

  if (mnpts > 1) {
    d <- rep(1, n)
    cur.per = 0
    if (pp) {
      print("Counting Neighbours...")
    }
    for (i in 1:(mnpts-1)) {
      tmp = (x[(1+i):n] - x[1:(n-i)]) <= eps
      d[(1+i):n] = d[(1+i):n] + tmp
      d[1:(n-i)] = d[1:(n-i)] + tmp
      if (pp & floor(20 * i / (mnpts + 1)) > cur.per) {
        cur.per = floor(20 * i / (mnpts + 1))
        print(paste(toString(5 * cur.per), "%", sep = ""))
      }
    }
    # Identify core points
    p_core <- d >= mnpts
    n_core <- sum(p_core)
  } else if (mnpts == 1) {
    p_core = rep(TRUE, n)
    n_core = n
  }

  if(n_core == 0){
    warning("dipps::dbscan: No core points.")
    return(rep(0, n))
  }
  # Check there is more than one core point
  if(n_core == 1){
    return(as.numeric(p_core)[o])
  }
  # calc pairwise (adjacent) distances between core points (only)
  if (pp) {
    print("Grouping core points...")
  }
  x_core <- x[p_core]
  dist <- x_core[2:n_core] - x_core[1:(n_core-1)]
  clus <- c(1, 1 + cumsum(dist > eps))

  g = rep(0, length(x))
  g[p_core] = clus
  return(g[o])
}















