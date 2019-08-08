#' estimate_pcf_fast
#'
#' @description Fast estimation of the pair correlation function
#'
#' @param pattern Point pattern.
#' @param r A vector specifying the maximum radius of the sample circles and
#' the interval length between successive sample circles radii.
#' @param ... Arguments passed down to \code{kfun}.
#'
#' @details
#' The functions estimates the pair correlation function using the \code{ads}
#' package.
#'
#' @seealso
#' \code{\link{kfun}}
#'
#' @return fv.object
#'
#' @examples
#' pcf_species_a <- estimate_pcf_fast(species_a)
#'
#' @aliases estimate_pcf_fast
#' @rdname estimate_pcf_fast
#'
#' @references
#' Ripley, B.D. (1977) Modelling spatial patterns (with discussion). Journal of
#' the Royal Statistical Society, Series B, 39, 172-212.
#'
#' Stoyan, D, Kendall, W.S. and Mecke, J. (1995) Stochastic geometry and its
#' applications. 2nd edition. Springer Verlag.
#'
#' Stoyan, D. and Stoyan, H. (1994) Fractals, random shapes and point fields:
#' methods of geometrical statistics. John Wiley and Sons.

#' @export
estimate_pcf_fast <- function(pattern, r = NULL, ...){

  # unmark the pattern if marked
  if (spatstat::is.marked(pattern)) {

    pattern <- spatstat::unmark(pattern)
  }

  # calculat rmax if not provided
  # can be problematic depending on window
  if (is.null(r)) {

    # get maximum distance r
    upto <- spatstat::rmax.rule(W = pattern$window,
                                lambda = spatstat::intensity(pattern))

    # 25 steps to upto
    by <- upto / 25

    # combine into one vector
    r <- c(upto, by)
  }

  # error if length r
  else if (length(r) != 2 | !is.numeric(r)) {

    stop("r must be a vector with length.", call. = FALSE)
  }

  # convert to spp object
  pattern <- ads::ppp2spp(p = pattern)

  # calculate summary functions
  result <- ads::kfun(p = pattern,
                      upto = r[[1]], by = r[[2]],
                      nsim = 0)

  # get pcf only
  result <- cbind(r = result$r, result$g)

  # order identical to spatstat package
  result <- result[, c(1,3,2)]

  # set names identical to spatstat package
  names(result)[3] <- "iso"

  # construct fv spatstat object.
  result <- spatstat::fv(x = result,
                         argu = "r",
                         ylab = quote(g(r)),
                         valu = "iso",
                         fmla = . ~ r,
                         labl = c("r",
                                  "g[Pois](r)",
                                  "hat(g)[Ripley](r)"),
                         desc = c("distance argument r ",
                                  "theoretical Poisson g(r)",
                                  "isotropic-corrected estimate of g(r)"),
                         fname = "pcf")

  # return result
  return(result)
}
