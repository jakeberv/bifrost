#' Generate Viridis Color Palette for Ranked Parameters
#'
#' Creates a named color mapping for a set of numeric parameters (e.g., evolutionary rates)
#' using the \pkg{viridis} color palette. Parameters are sorted in ascending order and
#' assigned evenly spaced colors by rank. Numeric magnitude and distance between parameter
#' values are not encoded in the colors.
#'
#' @param params A named numeric vector of parameter values (e.g., rates). The names will be
#'   preserved and used to label the resulting color mapping.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{NamedColors}}{A named character vector of hex color codes, with names
#'   corresponding to the input parameter names, ordered by increasing parameter value and
#'   colored at evenly spaced positions by rank.}
#'   \item{\code{ParamColorMapping}}{A named numeric vector of the sorted parameter values,
#'   maintaining the same order and names as \code{NamedColors}.}
#' }
#'
#' @details
#' This function is useful for plotting results where parameters should be visually
#' distinguished by their ordering (e.g., rate shifts across a phylogeny). By using the
#' perceptually uniform viridis palette, it avoids misleading color interpretations common
#' with rainbow scales. Colors encode sorted rank only; they do not represent the magnitude
#' of a parameter or the distance between parameter values.
#'
#' @examples
#' if (requireNamespace("viridis", quietly = TRUE)) {
#'   library(viridis)
#'   set.seed(1)
#'   rates <- c(A = 0.1, B = 0.5, C = 0.9)
#'   color_scale <- generateViridisColorScale(rates)
#'
#'   # View the color assignments
#'   color_scale$NamedColors
#'
#'   # Plot with colors
#'   barplot(color_scale$ParamColorMapping,
#'           col = color_scale$NamedColors,
#'           main = "Rates with Viridis Colors")
#' }
#'
#' @seealso \code{\link[viridis:viridis]{viridis::viridis()}} for details on
#'   the color palette.
#'
#' @importFrom viridis viridis
#' @export
generateViridisColorScale <- function(params) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.", call. = FALSE)
  }

  # Sort parameters and keep their names
  sorted_indices <- order(params)
  sorted_params <- params[sorted_indices]

  # Generate evenly spaced colors for the sorted parameter ranks
  colors <- viridis(length(sorted_params))

  # Associate each color with its state (name), using the sorted order
  named_sorted_colors <- setNames(colors, names(sorted_params))

  # Create a second list for the actual parameter values and their associated colors, also in sorted order
  param_color_mapping <- setNames(sorted_params, names(sorted_params))

  # Return both the original named sorted colors and the param-color mapping
  return(list("NamedColors" = named_sorted_colors, "ParamColorMapping" = param_color_mapping))
}
