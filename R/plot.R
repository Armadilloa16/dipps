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
#' @param x.var Parameter
#' @param y.var Parameter
#' @param return.plot Parameter
#' @param print.plot Parameter
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
#' }
#' @export
spatial_plot <- function(df.peak,
                         df.spec,
                         plot.var = "intensity",
                         id.var = "Acq",
                         x.var = "X",
                         y.var = "Y",
                         return.plot = FALSE,
                         print.plot = TRUE){
  # TODO: Add checks on input.

  # Check for multiple peaks per acq
  tmp = table(df.peak[, id.var])
  if (sum(tmp > 1) > 0){
    stop(paste("dipps::spatialPlot: Multiple peaks detected in ",
               id.var, ": ", toString(names(tmp[tmp>1])), sep = ""))
  }

  x.min = min(df.spec[, x.var])
  x.max = max(df.spec[, x.var])
  y.min = min(df.spec[, y.var])
  y.max = max(df.spec[, y.var])

  df.plot = merge(df.spec[, c(x.var, y.var, id.var)],
                  df.peak[, c(plot.var,     id.var)],
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

  if (return.plot | print.plot) {
    p = (ggplot2::ggplot(df.plot, ggplot2::aes(df.plot[, x.var],
                                               df.plot[, y.var]))
         + ggplot2::geom_tile(ggplot2::aes(fill = eval(as.symbol(plot.var)),
             alpha = as.numeric(!is.na(eval(as.symbol(plot.var))))),
                              colour = NA)
         + ggplot2::geom_tile(data = df.plot,
                              alpha = 0.5*as.numeric(df.plot$empty))
         + ggplot2::coord_fixed()
         + ggplot2::guides(alpha = FALSE,
                           fill = ggplot2::guide_colourbar(plot.var))
         + ggplot2::xlab(x.var)
         + ggplot2::ylab(y.var))
  }
  if (print.plot) {
    print(p)
  }

  if (return.plot){
    return(p)
  } else {
    return(df.plot)
  }
}






