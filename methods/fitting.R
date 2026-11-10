# Construct a fitter once per setting/condition and reuse its compiled backend.
make_fitter <- function(model, backend, cores) {
    control <- optimizer_control()
    number_of_starts <- 11L

    function(sample, method) {
        stopifnot(method %in% paper_methods)
        kernel_method <- if (method == "CC_DRGO") "CC" else method
        kernel <- make_method_kernel(kernel_method, model)
        kernel$itr_gradient <- backend$itr_gradient
        kernel$itr_loss <- backend$itr_loss
        last_solver <- NULL

        kernel$dc_policy <- function(Z, r, n_starts = number_of_starts,
                                     max_dc = 200L, seed = 1L, source_start = NULL) {
            # max_dc is a compatibility argument; control sets Newton/DC limits.
            solution <- solve_policy(
                Z, r, kernel, seed, n_starts, source_start,
                ctl = control, ls_start = TRUE, cores = cores
            )
            last_solver <<- solution
            kernel$diagnostics$starts <- n_starts
            kernel$diagnostics$optimizer_calls <- sum(vapply(
                solution$starts, function(start) start$bfgs_calls, integer(1)
            ))
            kernel$diagnostics$optimizer_nonconvergence <- sum(vapply(
                solution$starts, function(start) start$bfgs_failures, integer(1)
            ))
            kernel$diagnostics$dc_iterations <- vapply(
                solution$starts, function(start) nrow(start$trace), integer(1)
            )

            # Truth is used only after fitting, for an optimization diagnostic.
            grid <- expand.grid(x1 = seq(-1, 1, length.out = 11L),
                                x2 = seq(-1, 1, length.out = 11L))
            truth_coefficients <- qr.solve(model$b(grid), model$tau(grid))
            kernel$diagnostics$empirical_value_fit <- kernel$empirical_value(Z, r, solution$theta)
            kernel$diagnostics$empirical_value_truth <- kernel$empirical_value(Z, r, truth_coefficients)
            solution$theta
        }

        training <- list(source = sample$trial, target = sample$target,
                         n1 = nrow(sample$trial), n0 = nrow(sample$target))
        # Preserve the cluster study's fold and optimizer seed offsets.
        fit_seed <- 1000L + sample$seed
        arguments <- list(data = training, n_starts = number_of_starts,
                          max_dc = 200L, seed = fit_seed)

        if (method == "CC_DRGO") {
            scores <- kernel$fit_cross_fitted_scores(
                training, n_folds = 5L, seed = fit_seed,
                need_target = TRUE, need_weight = TRUE
            )
            correction_target <- rep(NA_real_, training$n0)
            coefficients_by_fold <- matrix(NA_real_, 5L, ncol(model$z(model$example)))
            for (fold in 1:5) {
                split <- scores$fold_fits[[fold]]
                fold_training <- list(
                    source = training$source[split$source_train_index, , drop = FALSE],
                    target = training$target[split$target_train_index, , drop = FALSE]
                )
                fold_training$n1 <- nrow(fold_training$source)
                fold_training$n0 <- nrow(fold_training$target)
                inner_scores <- kernel$fit_cross_fitted_scores(
                    fold_training, n_folds = 5L, seed = 40000L + fold,
                    need_target = TRUE, need_weight = TRUE
                )
                correction <- drgo_correction(inner_scores, model)
                coefficients_by_fold[fold, ] <- correction$gamma
                correction_target[split$target_eval_index] <- drop(
                    model$z(training$target[split$target_eval_index, , drop = FALSE]) %*%
                        correction$gamma
                )
                stopifnot(
                    !length(intersect(split$source_train_index, split$source_eval_index)),
                    !length(intersect(split$target_train_index, split$target_eval_index))
                )
            }
            stopifnot(all(is.finite(correction_target)))
            coefficients <- kernel$dc_policy(
                model$b(training$target),
                (scores$target$Gamma0_hat - correction_target) / training$n0,
                n_starts = number_of_starts, max_dc = 200L, seed = fit_seed + 1L
            )
            effect_rmse <- NA_real_
            confounding_rmse <- sqrt(mean((correction_target - model$c(training$target))^2))
        } else {
            details <- switch(
                method,
                IP = do.call(kernel$fit_integrative_plugin, arguments),
                WY = do.call(kernel$fit_wuyang, arguments),
                Target = , Pooled = do.call(kernel$fit_naive_pooled_itr, arguments),
                Wu2020 = do.call(kernel$fit_wu2020, c(arguments, list(ntree = 500L))),
                CC = do.call(kernel$fit_fixed_alpha_path, arguments)
            )
            coefficients <- if (method == "CC") details$eta_by_alpha[[1L]] else details$eta_hat
            effect_rmse <- if (method == "IP") {
                sqrt(mean((details$tau_target - model$tau(training$target))^2))
            } else NA_real_
            confounding_rmse <- if (method == "CC") {
                sqrt(mean((details$pointwise_c$c_target - model$c(training$target))^2))
            } else NA_real_
        }

        stopifnot(all(is.finite(coefficients)), !is.null(last_solver),
                  kernel$diagnostics$starts == number_of_starts)
        selected <- last_solver$starts[[last_solver$chosen]]
        list(
            theta = coefficients,
            converged_starts = sum(last_solver$converged),
            starts_total = number_of_starts,
            chosen_start = last_solver$chosen,
            chosen_status = selected$status,
            chosen_gradient = selected$terminal_gradient,
            retried_starts = sum(vapply(last_solver$starts, function(start) {
                !is.null(start$initial_attempt)
            }, logical(1))),
            total_dc_iterations = sum(vapply(last_solver$starts, function(start) {
                nrow(start$trace) + if (is.null(start$initial_attempt)) 0L else nrow(start$initial_attempt$trace)
            }, integer(1))),
            rank_deficient_fits = kernel$diagnostics$rank_deficient_fits,
            bfgs_calls = kernel$diagnostics$optimizer_calls,
            bfgs_nonconvergence = kernel$diagnostics$optimizer_nonconvergence,
            optimization_gap_lower_bound = max(
                0, kernel$diagnostics$empirical_value_truth - kernel$diagnostics$empirical_value_fit
            ),
            tau_rmse = effect_rmse,
            c_rmse = confounding_rmse
        )
    }
}
