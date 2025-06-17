#' @title Gaussian Kernel Function
#'
#' @description
#' Computes the Gaussian (or Radial Basis Function) kernel value for a given distance `h`.
#' This kernel is infinitely differentiable, implying very smooth functions and
#' a global influence of each data point.
#'
#' @param h A numeric value or vector representing the distance(s) between data points.
#' @return A numeric value or vector representing the Gaussian kernel similarity,
#'   which ranges from 0 to 1.
#'
#' @details
#' The Gaussian kernel is defined by the formula:
#' \ifelse{html}{\out{<i>K(h) = exp(-h<sup>2</sup>)</i>}}{\eqn{K(h) = \exp(-h^2)}}.
#' It is widely used due to its smoothness and tendency to produce
#' smoothly varying functions. As the distance `h` increases, the kernel value
#' rapidly approaches zero, indicating decreasing similarity.
#'
#' @keywords internal
#' # This function is not exported for direct user use.
#'
gaussianKernel <- function(h) {
  # The formula for the Gaussian kernel.
  # It takes the exponential of the negative square of the distance.
  base::exp(-h^2)
}
#' @title Matérn 5/2 Kernel Function
#'
#' @description
#' Computes the Matérn 5/2 kernel value for a given distance `h`. This kernel is
#' twice differentiable, leading to functions that are relatively smooth but can
#' still capture some local irregularities, offering a good balance between
#' smoothness and flexibility.
#'
#' @param h A numeric value or vector representing the distance(s) between data points.
#' @return A numeric value or vector representing the Matérn 5/2 kernel similarity.
#'
#' @details
#' The Matérn 5/2 kernel is defined by the formula:
#' \ifelse{html}{\out{<i>K(h) =\left(1 + \sqrt{5}h + \frac{5}{3}h^2\right) \exp(-\sqrt{5}h)</i>}}{\eqn{K(h) = \left(1 + \sqrt{5}h + \frac{5}{3}h^2\right) \exp(-\sqrt{5}h)}}
#' It provides a good compromise between the very smooth Gaussian kernel and
#' less smooth Matérn kernels, making it a popular choice in Gaussian Process
#' and kernel-based models.
#'
#' @keywords internal
#' # This function is not exported for direct user use.
#'
matern52Kernel <- function(h) {
  # The formula for the Matérn 5/2 kernel.
  # Combines a polynomial term with an exponential decay.
  (1 + base::sqrt(5) * h + (5 / 3) * h^2) * base::exp(-base::sqrt(5) * h)
}
#' @title Matérn 3/2 Kernel Function
#'
#' @description
#' Computes the Matérn 3/2 kernel value for a given distance `h`. This kernel is
#' once differentiable, producing functions that are less smooth than the Gaussian
#' but more flexible in capturing sharp changes or local variations in the data.
#'
#' @param h A numeric value or vector representing the distance(s) between data points.
#' @return A numeric value or vector representing the Matérn 3/2 kernel similarity.
#'
#' @details
#' The Matérn 3/2 kernel is defined by the formula:
#' \ifelse{html}{\out{<i>K(h) = (1 + \sqrt{3}h) \exp(-\sqrt{3}h)</i>}}{\eqn{K(h) = \left(1 + \sqrt{3}h\right) \exp\left(-\sqrt{3}h\right)}}.
#' It is often chosen when a model needs to capture more abrupt changes or
#' non-smooth behaviors in the underlying function.
#'
#' @keywords internal
#' # This function is not exported for direct user use.
#'
matern32Kernel <- function(h) {
  # The formula for the Matérn 3/2 kernel.
  # Combines a linear term with an exponential decay.
  (1 + base::sqrt(3) * h) * base::exp(-base::sqrt(3) * h)
}
