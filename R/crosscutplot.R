crosscutplot<-function (x, y, ct = 0.25, xlab = "", ylab = "", main = "", ylim = NULL)
{
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(length(x)==length(y))
  stopifnot((ct>0)&(ct<=.5))
  qx1 <- stats::quantile(x, ct)
  qx2 <- stats::quantile(x, 1 - ct)
  qy1 <- stats::quantile(y, ct)
  qy2 <- stats::quantile(y, 1 - ct)
  use <- ((x <= qx1) | (x >= qx2)) & ((y <= qy1) | (y >= qy2))
  if (is.null(ylim))
    graphics::plot(x, y, xlab = xlab, ylab = ylab,
                   main = main, type="n")
  else graphics::plot(x, y, xlab = xlab, ylab = ylab, ylim=ylim,
                      main = main, type="n")
  graphics::points(x[use], y[use], pch = 16)
  graphics::points(x[!use], y[!use], col = "gray", pch = 16)
  graphics::abline(h = c(qy1, qy2))
  graphics::abline(v = c(qx1, qx2))
}
