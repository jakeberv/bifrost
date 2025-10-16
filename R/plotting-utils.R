#' Generate Scaled Viridis Color Palette for Rate Parameters
#'
#' Creates a named color mapping for a set of numeric parameters (e.g., evolutionary rates)
#' using the \pkg{viridis} color palette. Parameters are first sorted in ascending order and
#' normalized to the range \[0, 1\], then mapped to evenly spaced viridis colors for
#' intuitive visualization.
#'
#' @param params A named numeric vector of parameter values (e.g., rates). The names will be
#'   preserved and used to label the resulting color mapping.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{NamedColors}}{A named character vector of hex color codes, with names
#'   corresponding to the input parameter names, ordered by increasing parameter value.}
#'   \item{\code{ParamColorMapping}}{A named numeric vector of the sorted parameter values,
#'   maintaining the same order and names as \code{NamedColors}.}
#' }
#'
#' @details
#' This function is useful for plotting results where parameters should be visually
#' distinguished based on their magnitude (e.g., rate shifts across a phylogeny).
#' By using the perceptually uniform viridis palette, it avoids misleading color
#' interpretations common with rainbow scales.
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
#' @seealso [viridis::viridis()] for details on the color palette.
#'
#' @importFrom viridis viridis
#' @keywords internal
#' @noRd
generateViridisColorScale <- function(params) {
  # Sort parameters and keep their names
  sorted_indices <- order(params)
  sorted_params <- params[sorted_indices]

  # Normalize the sorted parameter values to a range from 0 to 1
  normalized_sorted_params <- (sorted_params - min(sorted_params)) / (max(sorted_params) - min(sorted_params))

  # Use the normalized values to get colors from the viridis palette
  colors <- viridis(length(normalized_sorted_params))

  # Associate each color with its state (name), using the sorted order
  named_sorted_colors <- setNames(colors, names(sorted_params))

  # Create a second list for the actual parameter values and their associated colors, also in sorted order
  param_color_mapping <- setNames(sorted_params, names(sorted_params))

  # Return both the original named sorted colors and the param-color mapping
  return(list("NamedColors" = named_sorted_colors, "ParamColorMapping" = param_color_mapping))
}
