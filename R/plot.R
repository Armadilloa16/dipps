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



#' Spatially Plot MSI Data
#'
#' Produces a spatial plot of the given MSI data, designed to be used with
#' combine_peaklists and associated load_* functions. Usually would be used on a
#' subset of the Peak-list produced by combine_peaklist, those near a particular
#' mass for example, see example below for an example of this example.
#'
#' Note \code{spatial_plot} assumes that there are no occurrences of multiple
#' peaks from the same spectrum in the input peak-list \code{df.peak}. Any such
#' occurrences should be dealt with before calling spatial_plot, as otherwise it
#' will throw an error.
#'
#' @param df.peak     A Peak-list data.frame as produced by
#'   \code{combine_peaklists}, or more commonly a subset thereof.
#' @param df.spec     A Spectra-list data.frame as produced by
#'   \code{combine_peaklists}.
#' @param plot.var    A string matching the variable name in \code{df.peak} to
#'   be plotted.
#' @param id.var      A string matching the variable name that uniquely
#'   identifies spectra in both \code{df.peak} and \code{df.spec}, thereby
#'   linking that information between them.
#' @param x.var       A string matching the variable name in \code{df.spec} that
#'   contains the X-coordinates of the spectra.
#' @param y.var       A string matching the variable name in \code{df.spec} that
#'   contains the Y-coordinates of the spectra.
#' @param return.plot A logical value that if TRUE will cause
#'   \code{spatial_plot} to return a ggplot object and if FALSE will instead
#'   return the data.frame underlying it.
#' @param print.plot  A logical value indicating whether the plot should be
#'   printed to the current graphics device.
#'
#' @return On successful completion returns either a ggplot object or the
#' data.frame used to produce said ggplot object, depending on the value of
#' \code{return.plot}.
#'
#' If \code{return.plot} is FALSE, the returned data.frame will contain rows
#' that represent unique X-Y coordinate pairs in a rectanglular region extending
#' from the minimum X-coordinate in \code{df.spec} minus one to the maximum plus
#' one (inclusive), and similarly from the minimum Y-coordinate in
#' \code{df.spec} minus one to the maximum plus one (inclusive), and for each
#' such unique X-Y coordinate pair four columns represent the relevant plotting
#' information:
#' \itemize{
#'   \item x.var: Columned named with the string \code{x.var}, containing the
#'     relevant X-coordinate values.
#'   \item y.var: Columned named with the string \code{y.var}, containing the
#'     relevant Y-coordinate values.
#'   \item plot.var: Column named with the string \code{plot.var}, containing
#'     the relevant values, including NA for missing values.
#'   \item empty: A logical column indicating X-Y coordinate pairs from which
#'     spectra were acquired, but no peaks are present in \code{df.peak}. This
#'     is used to distinguish such X-Y coordnate pairs from those in which no
#'     spectrum was acquired at all.
#'   }
#'
#' If \code{return.plot} is TRUE, the returned ggplot object is a geom_tile plot
#' with axes labelled with \code{x.var} and \code{y.var} respectively, a legend
#' for \code{plot.var}, X-Y coordinate pairs (or pixels) for which no spectrum
#' was acquired blank, and X-Y coordinate pairs (or pixels) for which a spectrum
#' was acquired but no peaks are present in \code{df.peak} in darkened grey. The
#' plot uses coord_fixed() to match X and Y scales.
#'
#' @seealso Intended to be used with \code{\link{combine_peaklists}}, associated
#'   \code{load_*} functions, and the \code{ggplot2} package
#'   (\url{https://cran.r-project.org/web/packages/ggplot2/index.html}),
#'   hadley's book on ggplot2 is a good resource, see his github
#'   (\url{https://github.com/hadley/ggplot2-book}) or you can access samples
#'   and/ or buy his book from the website (\url{http://ggplot2.org/book}).
#'
#' @examples
#' i.path = system.file("extdata", "test1", package = "dipps")
#' n.empty = combine_peaklists(i.path)
#' o.name = basename(i.path)
#' df.spec = load_speclist(o.name)
#' df.peak = load_peaklist(o.name)
#'
#' # Select a m/z window of interest, in this case m/z = 1570.677 +/- 0.1 Da
#' df.cal = subset(df.peak, abs(m.z - 1570.677) < 0.1)
#'
#' # Plot log-intensity of said m/z window.
#' df.cal$log.intensity = log1p(df.cal$intensity)
#' p = spatial_plot(df.cal, df.spec, plot.var = "log.intensity")
#'
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






