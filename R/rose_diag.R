##' Rose diagram in ggplot2 inspired from rose.diag in package circular.
#'
#' @param x A vector of circular coordinates in radiants \eqn{[0,2 \pi)}.
#' @param bins Number of bins
#' @param color Color of the line and of the fill
#' @param alpha Transparency
#' @param start The starting angle of the 0 (the North)
#' @param add Add the rose diag to an existing ggplot2 plot
#' @param template radiants or wind rose
#' @param direction clockwise  or counterclockwise
#' @return The plot in ggplot2 format.
#' @examples
#' require(circular)
#' x<-rwrappedstable(200, index=1.5, skewness=.5)
#' x1<-rwrappedstable(200, index=2, skewness=.5)
#' x2<-rwrappedstable(200, index=0.5, skewness=1)
#' rose_diag(x,bins=15,color="green")
#' rose_diag(x1,bins=15,color="blue",alpha=.5,add=T)
#' rose_diag(x2,bins=15,color="red",alpha=.5,add=T)
#' @export
rose_diag <- function(x, bins=15, color= "red", alpha = 1, start = 0, add = FALSE, template = "rad", direction = NULL) {
  # inspired from https://stackoverflow.com/questions/15125628/putting-mathematical-symbols-and-subscripts-mixed-with-regular-letters-in-r-ggpl
  if ( is.null(direction) & template == "wind_rose") direction = 1
  if ( is.null(direction) & template == "rad") direction = -1
  if (template == "rad") labels <- c(0, expression(pi/2), expression(pi), expression(3/2*pi))
  else if (template == "wind_rose") labels <- c("N", "E", "S", "W")
  else stop

  require(ggplot2)
  if (add == FALSE) {
    ggplot(mapping = aes(x = x/(2*pi))) +
      stat_bin(bins = bins, color = color, fill = color,alpha = alpha,boundary = 0, closed = "left") +
      scale_x_continuous(limits = c(0,1),
                         breaks = c(0,.25,.5,.75),
                         labels =  labels) +
      ylab(NULL) + xlab(NULL)  + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      coord_polar(start = start, direction = direction)
  } else{
    pp <- last_plot()
    print(pp +
            stat_bin(data = data.frame(x), mapping = aes(x = x/(2*pi)),bins = bins, color = color, fill = color,alpha = alpha,boundary = 0, closed = "left")  +
            coord_polar(start = pp$coordinates$start, direction = pp$coordinates$direction)
    )
  }
}
