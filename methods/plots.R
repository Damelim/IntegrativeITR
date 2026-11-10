# Only the configured start seed is passed to this function.
plot_decision_boundaries <- function(model, specification, seed, methods, fits, rows, path) {
    axis <- seq(-1, 1, length.out = 301L)
    grid <- expand.grid(x1 = axis, x2 = axis)
    policy_basis <- model$b(grid)
    true_effect <- model$tau(grid)
    true_confounding <- model$c(grid)
    panel_count <- length(methods) + 2L
    panel_columns <- min(3L, panel_count)
    panel_rows <- ceiling(panel_count / panel_columns)
    device_type <- if (Sys.info()[["sysname"]] == "Darwin") "quartz" else "cairo"
    png(path, width = panel_columns * 1000L, height = panel_rows * 1000L,
        res = 220L, type = device_type)
    on.exit(dev.off())
    par(mfrow = c(panel_rows, panel_columns), mar = c(3, 3, 3, 1),
        mgp = c(1.8, 0.6, 0), oma = c(2, 0, 3, 0))

    draw_zero <- function(score, color = "black", line_type = 1L) {
        if (min(score) < 0 && max(score) > 0) {
            contour(axis, axis, matrix(score, 301L), levels = 0, add = TRUE,
                    drawlabels = FALSE, col = color, lty = line_type, lwd = 2)
        }
    }
    draw_region <- function(score, title) {
        image(axis, axis, matrix(as.integer(score > 0), 301L),
              breaks = c(-0.5, 0.5, 1.5), col = c("#FFF0DF", "#E0F0F5"),
              xlab = expression(X[1]), ylab = expression(X[2]), main = title, asp = 1)
        box()
    }
    draw_region(true_effect, "True treatment effect")
    draw_zero(true_effect)
    draw_zero(true_confounding, "#2267AC", 2L)
    draw_region(true_effect + true_confounding, "Observed target contrast")
    draw_zero(true_effect + true_confounding, "#C24128")
    draw_zero(true_effect, "black", 2L)

    for (method in methods) {
        label <- if (method == "CC_DRGO") "CC-DRGO" else method
        fit <- fits[[method]]
        result <- rows[rows$method == method & rows$seed == seed, , drop = FALSE]
        if (is.null(fit) || nrow(result) != 1L || !is.finite(result$regret)) {
            plot.new()
            title(label)
            text(0.5, 0.5, "Missing / failed")
            next
        }
        score <- drop(policy_basis %*% fit$theta)
        draw_region(score, sprintf("%s; regret = %.5f", label, result$regret))
        draw_zero(true_effect, "black", 2L)
        draw_zero(score, "#C24128")
    }
    mtext(sprintf("%s\n%s | seed %d | n0=%d, n1=%d | d_mu=%d",
                  paper_setting_labels[[specification$setting]], specification$condition, seed,
                  specification$n0, specification$n1, specification$d_mu),
          outer = TRUE, line = 1)
    mtext("Method panels: black dashed = true boundary; red = fitted boundary.",
          side = 1, outer = TRUE, line = 0.3)
    dev.off()
    on.exit(NULL)
    if (!file.exists(path) || file.info(path)$size == 0) stop("Plot was not created: ", path)
}
