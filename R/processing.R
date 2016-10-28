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

#' Subsets a peaklist around known masses.
#'
#' Subsets a peaklist around given masses with the given tolerance.
#'
#' @param df.peak Parameter
#' @param mzs Parameter
#' @param margin Parameter
#' @param use_ppm Parameter
#' @param var Parameter
#'
#' @export
extract_masses <- function(df.peak, mzs,
                           margin = 0.3,
                           use_ppm = FALSE,
                           var = "m.z") {
  # TODO: Check inputs

  df.sub.dne = FALSE
  if (use_ppm) margin.ppm = margin
  for (i in 1:length(mzs)){
    if (use_ppm) margin = margin.ppm*mzs[i]*1e-6
    l = abs(df.peak[,var] - mzs[i]) <= margin
    if (all(!l)) next
    if (!df.sub.exists){
      df.sub <- transform(df.peak[l,], group = mzs[i])
      df.sub.exists = TRUE
    } else {
      df.sub <- rbind(df.sub, transform(df.peak[l,], group = mzs[i]))
    }
  }
  return(df.sub)
}


#' Tolerance Clustering
#'
#' Groups peaks together according to the equivalence classes induced by the
#' relation `are within tol of each other in mass'.
#'
#' @param df.peak Parameter
#' @param tol Parameter
#' @param var Parameter
#'
#' @examples
#' \dontrun{
#' dataset.name = "A1" # Peaklist files should be in ./A1/peaklists/
#' n.empty = combine_peaklists(dataset.name)
#' df.peak = load_peaklist(dataset.name)
#' df.tol  = tol_clus(df.peak)
#' # Count how many peaks in each group
#' counts = table(df.tol$group)
#' # Remove peaks belonging to groups with less than 100 peaks total.
#' df.sub = subset(df.tol, group %in% as.numeric(names(counts)[counts >= 100]))
#' }
#' @export
tol_clus <- function(df.peak,
                     tol = 0.1,
                     var = "m.z") {
  # TODO: Check inputs

  df.peak = df.peak[order(df.peak[, var]),]
  n = nrow(df.peak)
  df.peak <- transform(df.peak, group = 1)
  for (i in which(df.peak[2:n, var] - df.peak[1:(n-1), var] > tol)){
    df.peak[(i+1):n, "group"] <- df.peak[(i+1):n, "group"] + 1
  }
  return(df.peak)
}

#' Density Based Clustering
#'
#' Groups peaks together according to the equivalence classes induced by the
#' relation `are within tol of each other in mass'.
#'
#' @param df.peak Parameter
#' @param eps Parameter
#' @param mnpts Parameter
#' @param var Parameter
#' @param pp Parameter
#'
#' @seealso This is a univariate optimisation of the DBSCAN* algorithm as in
#' section 3 of:
#'
#' Campello, Ricardo JGB, Davoud Moulavi, and Joerg Sander.
#' "Density-based clustering based on hierarchical density estimates." In
#' Advances in Knowledge Discovery and Data Mining, pp. 160-172. Springer Berlin
#' Heidelberg, 2013.
#'
#' @examples
#' \dontrun{
#' dataset.name = "A1" # Peaklist files should be in ./A1/peaklists/
#' n.empty = combine_peaklists(dataset.name)
#' df.peak = load_peaklist(dataset.name)
#' df.tol  = dbscan(df.peak)
#' }
#' @export
dbscan <- function(df.peak, eps=0.05, mnpts=100, var="m.z", pp=TRUE){
  if (!any(var == names(df.peak))){
    stop("dipps::dbscan: Non-existent variable selected for clustering.")
  }
  # TODO: Add more checks on input

  n <- nrow(df.peak)
  # sort
  if(pp){
    print("Sorting...")
  }
  df.peak <- df.peak[order(df.peak[, var]), ]

  x <- df.peak[, var]
  count <- rep(0,n)
  if(pp){
    print("Counting Neighbours...")
    cur.per = 0
  }
  for(i in 1:(mnpts+1)){
    tmp = (x[(1+i):n] - x[1:(n-i)]) <= eps
    count[(1+i):n] = count[(1+i):n] + tmp
    count[1:(n-i)] = count[1:(n-i)] + tmp
    if (pp & floor(20 * i / (mnpts + 1)) > cur.per) {
      cur.per = floor(20 * i / (mnpts + 1))
      print(paste(toString(5 * cur.per), "%", sep = ""))
    }
  }
  # Identify core points
  p_core <- count >= mnpts
  n_core <- sum(p_core)
  if(n_core == 0){
    warning("dipps::dbscan: No core points")
    return(data.frame())
  }
  # Check there is more than one core point
  if(n_core == 1){
    df.peak$group <- as.numeric(p_core)
    return(df.peak)
  }
  # calc pairwise (adjacent) distances between core points (only)
  if (pp) {
    print("Grouping core points...")
  }
  x_core <- x[p_core]
  dist <- x_core[2:n_core] - x_core[1:(n_core-1)]
  clus <- c(1, 1 + cumsum(x_core > eps))

  df.peak$group = 0
  df.peak[p_core, "group"] = clus
  return(df.peak)
}















