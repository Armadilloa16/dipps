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


# This is for taking a subset of a peaklist centered around certain known m/z values,
# for example calibrants. fixed Da bins (use_ppm = FALSE) or ppm based tolerances
# (use_ppm = TRUE) are both supported.
mzMatch <- function(peaklist_in,mzList,binMargin=0.3,use_ppm=FALSE) {
  peaklist_subset_does_not_exist = TRUE
  if (use_ppm){
    binMargin_ppm <- binMargin
  }
  for (i in 1:length(mzList)){
    if (use_ppm){
      binMargin <- binMargin_ppm*mzList[i]/1000000
    }
    idx = which(abs(peaklist_in$m.z - mzList[i])<binMargin)
    if (length(idx) > 0){
      if (peaklist_subset_does_not_exist){
        peaklist_subset <- transform(peaklist_in[idx,],
                                     PeakGroup=mzList[i])
        peaklist_subset_does_not_exist = FALSE
      } else {
        peaklist_subset <- rbind(peaklist_subset,transform(peaklist_in[idx,],
                                                           PeakGroup=mzList[i]))
      }
    }
  }
  return(peaklist_subset)
}


# Does a peak-grouping, and annotates peaks by peakgroup, in case you want that.
# If you set a non-zero minGroupSize here any peaks not allocated to groups will be
# annotates peakgroup zero. This can always be done later though, with the table function
# in the localFunctions.R file for example, so I reccomend leaving minGroupSize = 0 here.
# You could modify tol (the tolerance used) if you wish however.
groupPeaks <- function(peaklist_in,tol = 0.1, minGroupSize = 0) {
  peaklist_in <- peaklist_in[order(peaklist_in$m.z),]
  nPeaks <- nrow(peaklist_in)
  peaklist_in <- transform(peaklist_in,PeakGroup = 1)
  for (i in which(peaklist_in[2:nPeaks,]$m.z - peaklist_in[1:(nPeaks-1),]$m.z > tol)){
    peaklist_in$PeakGroup[(i+1):nPeaks] <- peaklist_in$PeakGroup[(i+1):nPeaks] + 1
  }
  for (p in which(as.vector(table(peaklist_in$PeakGroup)) < minGroupSize)){
    peaklist_in[which(peaklist_in$PeakGroup == p),]$PeakGroup <- 0
  }
  return(peaklist_in)
}

dbscan_lw <- function(peaklist_in,eps=0.05,mnpts=100,cvar="m.z",pp=TRUE){
  # Implements DBSCAN* as in section 3 of:
  #
  # Campello, Ricardo JGB, Davoud Moulavi, and Joerg Sander.
  # "Density-based clustering based on hierarchical density estimates."
  # In Advances in Knowledge Discovery and Data Mining, pp. 160-172.
  # Springer Berlin Heidelberg, 2013.
  #
  # For one-dimensional objects only. Intended for clustering
  # peaks by m/z location in MALDI Imaging.

  if(sum(cvar==names(peaklist_in))==0){
    error("Non-existent variable selected for clustering.")
  }

  n <- nrow(peaklist_in)
  # sort
  if(pp){
    print("Sorting...")
  }
  peaklist_in <- peaklist_in[order(peaklist_in[,cvar]),]
  if(pp){
    print("Done")
  }
  p_locs <- peaklist_in[,cvar]
  p_core <- rep(0,n)
  if(pp){
    print("Counting neighbours")
  }
  for(i in 1:(mnpts+1)){
    temp = (p_locs[(1+i):n] - p_locs[1:(n-i)]) <= eps
    p_core[(1+i):n] = p_core[(1+i):n] + temp
    p_core[1:(n-i)] = p_core[1:(n-i)] + temp
    if(pp){
      print(paste(toString(i-1),"/",toString(mnpts)))
    }
  }
  # Identify core points
  p_core <- p_core >= mnpts
  n_core <- sum(p_core)
  if(n_core == 0){
    print("No core points")
    return(FALSE)
  }
  # Check there is more than one core point
  if(n_core == 1){
    peaklist_in$PeakGroup <- as.numeric(p_core)
    return(peaklist_in)
  }
  # calc pairwise (adjacent) distances between core points (only)
  p_locs <- p_locs[p_core]
  d_pair <- p_locs[2:n_core] - p_locs[1:(n_core-1)]
  clus <- c(1,1+cumsum(d_pair > eps))

  peaklist_in$PeakGroup = 0
  peaklist_in[p_core,"PeakGroup"] = clus
  return(peaklist_in)
}









DIPPS <- function(pl_in){
  # peaklist_all should be a data.frame with at least
  # three variables:
  #  - Acquisition (identifying unique spectra)
  #  - PeakGroup (identifying peakgroups)
  #  - Group (identifying regions to be compared by
  #     DIPPS), and should take three values:
  #      - 1 coding for the `downregulated' group, and
  #      - 2 coding for the `upregulated group.

  # Commented out -- meaning assume input is unique
  # Check for multiple peaks and remove.
  #   pl_in = unique(pl_in)

  nSpec_d = length(unique(subset(pl_in,Group == 1)$Acquisition))
  nSpec_u = length(unique(subset(pl_in,Group == 2)$Acquisition))

  prop = ddply(pl_in,
               c("PeakGroup","Group"),
               summarise,
               p = length(Acquisition)
  )
  prop[prop$Group == 1,"p"] = prop[prop$Group == 1,"p"]/nSpec_d
  prop[prop$Group == 2,"p"] = prop[prop$Group == 2,"p"]/nSpec_u

  prop <- reshape(prop,
                  timevar = "Group",
                  idvar="PeakGroup",
                  direction="wide")
  names(prop) = c("PeakGroup","p.d","p.u")

  prop = replace(prop,is.na(prop),0)
  prop$d = prop$p.u - prop$p.d

  return(prop)
}




dippsHeur <- function(pl_in,dsum_in){
  # Takes the output of DIPPS as input, and calculates
  # a heuristic cutoff for the `optimal' number of
  # variables with highest DIPPS.

  u.m = dcast(subset(pl_in,Group==2),
              Acquisition~PeakGroup,
              value.var="Group")
  acq = u.m[,1]
  u.m = as.matrix(u.m[,-1])
  u.m[!is.na(u.m)] = 1
  u.m[is.na(u.m)] = 0
  acq = data.frame(Acquisition = acq,
                   nPeaks = rowSums(u.m))
  u.m = u.m/sqrt(acq$nPeaks)
  c.u = colMeans(u.m)
  c.u = c.u/norm(c.u,"2")

  temp = match(colnames(u.m),dsum_in$PeakGroup)
  dsum_in$c.u = 0
  dsum_in[temp,"c.u"] = c.u

  # Find data-driven cutoff according to the DIPPS method using cosine distance.
  curMinCosD = 10
  sortedDIPPS = sort(dsum_in$d,index.return=TRUE)
  # vN = 1:floor(nrow(Summary_merged)/2)
  vN = 1:nrow(dsum_in)
  cosD = vN
  for (n in vN){
    dsum_in$t = 0
    dsum_in[tail(sortedDIPPS$ix,n),]$t = 1/sqrt(n)
    cosD[n] = 1 - sum(dsum_in$t * dsum_in$c.u)
  }
  nStar = vN[which.min(cosD)]

  return(nStar)
}














