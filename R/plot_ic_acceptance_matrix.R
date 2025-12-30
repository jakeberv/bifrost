#' Plot IC Acceptance Matrix with Optional Rate-of-Improvement Overlay
#'
#' @description
#' Create a two-layer base R plot that visualizes information criterion (IC) scores
#' across a sequence of sub-model evaluations, highlighting which steps were
#' \emph{accepted} vs \emph{rejected}. Optionally, a secondary y-axis overlays the
#' \strong{rate of improvement} (first difference of IC scores) as a line with markers.
#'
#' @details
#' The function expects a two-column object where:
#' \itemize{
#'   \item Column 1 contains the IC score at each step (numeric; lower is better).
#'   \item Column 2 contains an indicator for acceptance (0 = rejected, 1 = accepted).
#' }
#' The first IC value is treated as the \emph{baseline} and is plotted as a larger
#' black point and labeled. Accepted steps are drawn as blue filled points connected
#' by a thin line; rejected steps are drawn as small red crosses. When
#' \code{plot_rate_of_improvement = TRUE}, the function overlays a secondary y-axis on
#' the right that shows \code{diff(IC)} values (the per-step change in IC; more negative
#' implies improvement).
#'
#' The function uses only base graphics. It sets plot margins and \code{mgp} via
#' \code{par()}, and (when overlaying) uses \code{par(new = TRUE)} to layer the IC plot over the
#' rate-of-improvement axes. It does not restore previous graphical parameters.
#'
#' @param matrix_data A two-column \code{matrix} or \code{data.frame}. Column 1 must be
#'   numeric IC scores in evaluation order; Column 2 must be a logical or numeric flag
#'   (0/1) indicating whether the step was accepted.
#' @param plot_title \code{character(1)}. Title to draw above the plot.
#' @param plot_rate_of_improvement \code{logical(1)}. If \code{TRUE}, overlay the
#'   first differences of the IC series on a secondary (right) y-axis along with a
#'   horizontal reference line at zero.
#'
#' @return Invisibly returns \code{NULL}. Called for its plotting side effects.
#'
#' @details
#' **Axes and scaling.** Tick marks for the primary (IC) x/y axes are computed with
#' \code{pretty()} to give clean bounds. The secondary axis for the rate of improvement
#' uses fixed limits (\code{c(-400, 150)}) inside the function; adjust in source if your
#' expected \code{diff(IC)} range differs substantially.
#'
#' @examples
#' ic <- c(-1000, -1012, -1008, -1025, -1020, -1030)
#' accepted <- c(1, 0, 1, 0, 1)  # steps 2..6 relative to baseline
#' mat <- cbind(ic, c(1, accepted))  # mark baseline as accepted for plotting
#' plot_ic_acceptance_matrix(mat, plot_title = "IC Path")
#' # Avoid non-ASCII glyphs in titles on CRAN/CI:
#' plot_ic_acceptance_matrix(mat, plot_rate_of_improvement = TRUE)
#'
#' @seealso
#' \code{\link[graphics]{par}}, \code{\link[graphics]{plot}}, \code{\link[graphics]{axis}},
#' \code{\link[graphics]{lines}}, \code{\link[graphics]{points}}, \code{\link[graphics]{legend}},
#' \code{\link[graphics]{mtext}}, \code{\link[graphics]{title}}
#'
#' @importFrom graphics axis legend lines mtext par plot points text title
#' @importFrom grDevices rgb
#' @export
plot_ic_acceptance_matrix <- function(matrix_data,
                                      plot_title = "IC Acceptance Matrix Scatter Plot",
                                      plot_rate_of_improvement = TRUE) {
  # Adjust margins for balanced spacing
  par(mar = c(5, 5.5, 4, 6), mgp = c(3, 0.6, 0))  # Adjust mgp to move tick labels closer

  # Extract y-values and category values
  y_values <- matrix_data[, 1]
  categories <- matrix_data[, 2]

  # Calculate rate of improvement (differences between consecutive IC scores)
  rate_of_improvement <- diff(y_values)

  # Define x-axis and y-axis limits using pretty ticks with padding
  x_values <- seq_along(y_values)
  x_ticks <- pretty(x_values)
  y_ticks <- pretty(y_values)
  x_limits <- range(x_ticks)
  y_limits <- range(y_ticks)

  # Define limits for the rate of improvement (for the secondary y-axis)
  rate_ticks <- pretty(range(rate_of_improvement))
  rate_limits <- c(-400, 150) #range(rate_ticks)

  # Identify the baseline IC as the first IC score
  baseline_ic <- y_values[1]

  # Plot the rate of improvement optionally
  if (plot_rate_of_improvement) {
    plot(
      x_values[-1], rate_of_improvement,  # Rate of improvement (x values shifted for diff())
      col = NA,  # Suppress default plotting
      type = "n", lty = "solid", lwd = 1,  # Set up the plot environment
      xlab = "", ylab = "",  # Suppress axis labels for overlay
      xlim = x_limits, ylim = rate_limits,  # Secondary y-axis scaling
      xaxt = "n", yaxt = "n", bty = "n"  # Suppress axes for overlay
    )

    # Add a black horizontal line at y = 0 with restricted x-range
    lines(
      x = c(min(x_values), max(x_ticks) + 20),  # Extend from data limit to axis limit
      y = c(0, 0),  # Horizontal line at y = 0
      col = rgb(0, 0, 0, alpha = 0.5), lty = 1, lwd = 0.7)

    # Plot the rate of improvement curve
    lines(
      x_values[-1], rate_of_improvement,
      col = "grey",  # Semi-transparent black line
      lty = "solid", lwd = 0.8  # Thin solid line
    )

    # Add small black dots on the rate curve for accepted shifts
    accepted_x <- which(categories[-1] == 1)  # Accepted shifts correspond to categories == 1
    points(
      x_values[accepted_x + 1], rate_of_improvement[accepted_x],  # Offset by 1 for diff()
      col = rgb(0, 0, 0, alpha = 0.5), pch = 16, cex = 0.3  # Small black dots
    )

    # Add the secondary y-axis for rate of improvement with transparency
    axis(
      4, at = rate_ticks, labels = rate_ticks,
      las = 1, cex.axis = 0.75, tck = -0.02,
      col = rgb(0, 0, 0, alpha = 0.5),         # Color of ticks matches line transparency
      col.axis = rgb(0, 0, 0, alpha = 0.5)    # Color of tick labels matches line transparency
    )
  }

  # Plot the IC scores on top
  par(new = TRUE)  # Enable overlaying
  plot(
    x_values, y_values,
    col = NA,  # Suppress default point plotting; add manually
    xlab = "Sub-model evaluated",
    ylab = "",  # Leave blank; we'll use mtext for the y-axis label
    type = "n",  # Suppress plotting, just set up the environment
    xlim = x_limits,  # Adjust x-axis to tick-aligned limits
    ylim = y_limits,  # Adjust y-axis to IC score limits
    xaxt = "n",  # Suppress default x-axis
    yaxt = "n",  # Suppress default y-axis
    cex.lab = 1,  # Ensure consistent font size for axis labels
    bty = "n"  # Remove the box around the plot
  )

  # Custom title placement
  title(
    main = plot_title,
    line = 2,   # Adjust this value to raise/lower the title
    cex.main = 1.0 # Optional: Adjust title size
  )

  # Manually add x-axis
  axis(
    1, at = x_ticks, labels = x_ticks,
    cex.axis = 1, tck = -0.02  # Ticks and labels for x-axis
  )

  # Manually add y-axis for IC scores
  axis(
    2, at = y_ticks, labels = y_ticks,
    las = 1, cex.axis = 0.75, tck = -0.02  # Horizontal labels and ticks for y-axis
  )

  # Add the IC Score y-axis label using mtext for custom positioning
  mtext("IC Score", side = 2, line = 3.5, cex = 0.6)  # Customize 'line' as needed

  # Plot the rejected shifts (red dots, smaller size) first (behind accepted shifts)
  points(
    x_values[categories == 0], y_values[categories == 0],
    col = rgb(1, 0, 0, alpha = 1), pch = 3, cex = 0.4, lwd = 0.3  # Smaller red dots for rejected shifts
  )

  # Combine all blue dots (excluding the baseline IC) for a continuous line
  blue_x <- x_values[categories == 1]
  blue_y <- y_values[categories == 1]

  # Plot the blue line connecting all accepted shifts (BEHIND the dots)
  lines(
    blue_x, blue_y,
    col = "blue", type = "l", lwd = 1.1  # Thinner blue line
  )

  # Plot the blue dots for accepted shifts with a fine black outline (ON TOP of the line)
  points(
    x_values[-1][categories[-1] == 1], y_values[-1][categories[-1] == 1],
    col = "black", bg = "blue", pch = 21, cex = 0.8, lwd = 0.3  # Blue fill with hairline black outline
  )

  # Plot the baseline IC value as a larger black dot (after lines/dots for layering)
  points(
    x = x_values[1], y = baseline_ic, col = "black", pch = 19, cex = 1.0  # Black dot for baseline IC
  )

  # Add a label for the baseline IC
  text(
    x = x_values[1] + 2, y = baseline_ic, labels = paste0(round(baseline_ic, 2)),
    pos = 4, col = "black", cex = 0.6
  )

  # Add the label for the minimum accepted IC score
  min_accepted_index <- which(categories == 1 & y_values == min(y_values[categories == 1]))
  min_accepted_value <- y_values[min_accepted_index]

  text(
    x = min_accepted_index,
    y = min_accepted_value - diff(range(y_values)) * 0.02,  # Slight vertical offset
    labels = paste0(round(min_accepted_value, 2)),
    pos = 1, col = "black", cex = 0.6
  )

  # Add a clean legend for IC scores, rate of improvement, and baseline IC
  legend(
    "topright",
    inset = c(0.04, -0.10),
    legend = c("Rejected shift", "Accepted shift", "IC Change", "Baseline IC"),
    col = c("red", "blue", rgb(0, 0, 0, alpha = 0.5), "black"),
    lty = c(NA, 1, ifelse(plot_rate_of_improvement, 1, NA), NA),  # Line for "Accepted shift" and "IC Change"
    pch = c(3, 21, ifelse(plot_rate_of_improvement, NA, 19), 19), # Cross (3), Dot with line (21), "Baseline IC" (19)
    pt.bg = c(NA, "blue", NA, NA),  # Background for "Accepted shift"
    pt.lwd = c(NA, 0.5, NA, 0),     # Fine outline for "Accepted shift"
    cex = 0.65, bty = "n", xpd = TRUE  # Slightly smaller legend text and allow margin overlap
  )
}
