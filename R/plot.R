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

#' Generate a \code{data.frame} for plotting with \code{ggplot2}
#'
#' @param df.peak Parameter
#' @param df.spec Parameter
#' @param plot.var Parameter
#' @param id.var Parameter
#'
#' @examples
#' \dontrun{
#' dataset.name = "A1"
#' combine_peaklists(dataset.name)
#' df.spec = load_speclist(dataset.name)
#' df.peak = load_peaklist(dataset.name)
#' # Select a m/z window of interest.
#' df.cal = subset(df.peak, abs(m.z - 1570.677) < 0.1)
#' if (sum(table(df.cal$Acq) > 1) > 0){
#'   stop(paste("multiple peaks from the same spectrum in given window.",
#'              "This should be dealt with before proceeding."))
#' }
#' df.cal$log.intensity = log1p(df.cal$intensity)
#' # Plot log-intensity of said m/z window.
#' df.plot = spatial_plot(df.cal, df.spec, plot.var = "log.intensity")
#' p = (ggplot(df.plot, aes(X, Y))
#'      + geom_tile(aes(fill=log.intensity,
#'                      alpha=as.numeric(!is.na(log.intensity))),
#'                  colour=NA)
#'      + geom_tile(data = df.plot, alpha=0.5*as.numeric(df.plot$empty))
#'      + coord_fixed()
#'      + guides(alpha = FALSE))
#' print(p)
#' }
#' @export
spatial_plot <- function(df.peak,
                         df.spec,
                         plot.var = "intensity",
                         id.var = "Acq"){
  # TODO: Add checks on input.

  # Check for multiple peaks per acq
  tmp = table(df.peak[, id.var])
  if (sum(tmp > 1) > 0){
    stop(paste("dipps::spatialPlot: Multiple peaks detected in ",
               id.var, ": ", toString(names(tmp[tmp>1])), sep = ""))
  }

  x.min = min(df.spec$X)
  x.max = max(df.spec$X)
  y.min = min(df.spec$Y)
  y.max = max(df.spec$Y)

  df.plot = merge(df.spec[, c("X", "Y", id.var)],
                  df.peak[, c(plot.var, id.var)],
                  all.x = TRUE)
  df.plot$empty = is.na(df.plot[, plot.var])
  df.plot[, id.var] = NULL

  df.plot = merge(data.frame(X = rep((x.min-1):(x.max+1),
                                     rep(y.max - y.min + 3,
                                         x.max - x.min + 3)),
                             Y = rep((y.min-1):(y.max+1),
                                     x.max - x.min + 3)),
                  df.plot, all.x = TRUE)
  df.plot[is.na(df.plot$empty), "empty"] = FALSE

  return(df.plot)
}






