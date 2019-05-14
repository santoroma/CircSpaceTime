##' Rose diagram in ggplot2 inspired from rose.diag in package circular.
#'
#' @param x a vector of circular coordinates in radiants \eqn{[0,2 \pi)}.
#' @param bins number of bins
#' @param color color of the line and of the fill
#' @param alpha transparency
#' @param start the starting angle of the 0 (the North)
#' @param add add the rose_diag to an existing ggplot2 plot
#' @param template radiants or wind rose. the values are \code{c("rad","wind_rose")}..
#' default is \code{"rad"}
#' @param direction 1, clockwise; -1, anticlockwise. For template = "rad" direction is -1 while
#' for template = "wind_rose" direction is 1.
#' @return The plot in ggplot2 format.
#' @examples
#'
#'
#' library(CircSpaceTime)
#' x <- circular::rwrappedstable(200, index = 1.5, skewness = .5)
#' x1 <- circular::rwrappedstable(200, index = 2, skewness = .5)
#' x2 <- circular::rwrappedstable(200, index = 0.5, skewness = 1)
#' rose_diag(x, bins = 15, color = "green")
#' rose_diag(x1, bins = 15, color = "blue", alpha = .5, add = TRUE)
#' rose_diag(x2, bins = 15, color = "red", alpha = .5, add = TRUE)
#'
#' @export
rose_diag <- function(x, bins=15, color= "red", alpha = 1, start = 0, add = FALSE, template = "rad", direction = NULL) {
  # inspired from https://stackoverflow.com/questions/15125628/putting-mathematical-symbols-and-subscripts-mixed-with-regular-letters-in-r-ggpl
  x <- x %% (2*pi)
  if ( is.null(direction) & template == "wind_rose") direction = 1
  if ( is.null(direction) & template == "rad") direction = -1
  if (template == "rad") labels <- c(0, expression(pi/2), expression(pi), expression(3/2*pi))
  else if (template == "wind_rose") labels <- c("N", "E", "S", "W")
  else stop("You should specify the template 'rad' or 'wind_rose'")
  if (add == FALSE) {
    ggplot2::ggplot(mapping = ggplot2::aes(x = x/(2*pi))) +
      ggplot2::stat_bin(bins = bins, color = color, fill = color,alpha = alpha,boundary = 0, closed = "left") +
      ggplot2::scale_x_continuous(limits = c(0,1),
                         breaks = c(0,.25,.5,.75),
                         labels =  labels) +
      ggplot2::ylab(NULL) + ggplot2::xlab(NULL)  +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::coord_polar(start = start, direction = direction)
  } else{
    pp <- ggplot2::last_plot()
    print(pp +
            ggplot2::stat_bin(data = data.frame(x), mapping = ggplot2::aes(x = x/(2*pi)),bins = bins, color = color, fill = color,alpha = alpha,boundary = 0, closed = "left")  +
            ggplot2::coord_polar(start = pp$coordinates$start, direction = pp$coordinates$direction)
    )
  }
}
