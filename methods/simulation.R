# Checkpoints are internal; CSV tables and PNG figures are the public results.
validate_configuration <- function(configuration) {
    positive_integer <- function(value) {
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
            value >= 1 && value == floor(value)
    }
    stopifnot(positive_integer(configuration$replications),
              positive_integer(configuration$start_seed),
              positive_integer(configuration$cores))
    if (configuration$start_seed + configuration$replications - 1 > .Machine$integer.max - 900000L) {
        stop("The requested seed range exceeds the supported integer range.")
    }
    selections <- list(setting = paper_settings, condition = paper_conditions, methods = paper_methods)
    for (name in names(selections)) {
        values <- configuration[[name]]
        if (!is.character(values) || !length(values) || anyNA(values) ||
            anyDuplicated(values) || !all(values %in% selections[[name]])) {
            stop("Invalid ", name, " in the run configuration.")
        }
    }
    if (length(configuration$setting) != 1L || length(configuration$condition) != 1L) {
        stop("Each run.R executes exactly one setting and one condition.")
    }
    if (!is.character(configuration$output_directory) ||
        length(configuration$output_directory) != 1L ||
        is.na(configuration$output_directory) || !nzchar(configuration$output_directory)) {
        stop("output_directory must be a nonempty path.")
    }
    invisible(configuration)
}

save_checkpoint <- function(checkpoint, path) {
    temporary_path <- paste0(path, ".tmp")
    saveRDS(checkpoint, temporary_path)
    if (!file.rename(temporary_path, path)) {
        stop("Could not replace checkpoint: ", path)
    }
}

empty_results <- function() {
    data.frame(
        setting = character(), condition = character(), method = character(),
        seed = integer(), status = character(), regret = numeric(),
        test_MCSE = numeric(), seconds = numeric(), fit_seconds = numeric(), evaluation_seconds = numeric(),
        converged_starts = integer(), starts_total = integer(),
        chosen_start = integer(), chosen_status = character(), chosen_gradient = numeric(),
        retried_starts = integer(), total_dc_iterations = integer(),
        rank_deficient_fits = integer(), bfgs_calls = integer(), bfgs_nonconvergence = integer(),
        optimization_gap_lower_bound = numeric(), tau_rmse = numeric(), c_rmse = numeric(),
        error = character(), warnings = character(), stringsAsFactors = FALSE
    )
}

summarize_regret <- function(results, expected_replications) {
    keys <- unique(results[, c("setting", "condition", "method"), drop = FALSE])
    if (!nrow(keys)) return(data.frame())
    summaries <- lapply(seq_len(nrow(keys)), function(index) {
        key <- keys[index, , drop = FALSE]
        rows <- results[results$setting == key$setting & results$condition == key$condition &
                        results$method == key$method, , drop = FALSE]
        valid <- rows$status %in% c("completed", "optimizer_warning") & is.finite(rows$regret)
        values <- rows$regret[valid]
        count <- length(values)
        cbind(key, data.frame(
            n_expected = expected_replications, n_available = count,
            n_completed = sum(rows$status == "completed"),
            n_optimizer_warning = sum(rows$status == "optimizer_warning"),
            n_failed = sum(rows$status == "failed"),
            mean_regret = if (count) mean(values) else NA_real_,
            sd_regret = if (count > 1L) sd(values) else NA_real_,
            replication_MCSE = if (count > 1L) sd(values) / sqrt(count) else NA_real_,
            mean_seconds = if (any(valid)) mean(rows$seconds[valid]) else NA_real_
        ))
    })
    do.call(rbind, summaries)
}

paired_comparisons <- function(results, methods) {
    output <- data.frame(
        setting = character(), condition = character(), method_A = character(),
        method_B = character(), n_pairs = integer(), mean_A_minus_B = numeric(),
        paired_MCSE = numeric(), stringsAsFactors = FALSE
    )
    if (length(methods) < 2L || !nrow(results)) return(output)
    cases <- unique(results[, c("setting", "condition"), drop = FALSE])
    for (index in seq_len(nrow(cases))) {
        case <- cases[index, , drop = FALSE]
        rows <- results[results$setting == case$setting & results$condition == case$condition &
                        results$status %in% c("completed", "optimizer_warning"), , drop = FALSE]
        for (pair in combn(methods, 2L, simplify = FALSE)) {
            left <- rows[rows$method == pair[1L], c("seed", "regret"), drop = FALSE]
            right <- rows[rows$method == pair[2L], c("seed", "regret"), drop = FALSE]
            matched <- merge(left, right, by = "seed")
            difference <- matched$regret.x - matched$regret.y
            difference <- difference[is.finite(difference)]
            count <- length(difference)
            output <- rbind(output, data.frame(
                setting = case$setting, condition = case$condition,
                method_A = pair[1L], method_B = pair[2L], n_pairs = count,
                mean_A_minus_B = if (count) mean(difference) else NA_real_,
                paired_MCSE = if (count > 1L) sd(difference) / sqrt(count) else NA_real_
            ))
        }
    }
    output
}

write_result_tables <- function(results, directory, replications, methods) {
    write.csv(results, file.path(directory, "regret.csv"), row.names = FALSE)
    write.csv(summarize_regret(results, replications),
              file.path(directory, "summary.csv"), row.names = FALSE)
    write.csv(paired_comparisons(results, methods),
              file.path(directory, "paired_comparisons.csv"), row.names = FALSE)
}

fit_and_evaluate <- function(sample, specification, model, method, seed, fitter, evaluation = NULL) {
    started <- proc.time()[["elapsed"]]
    warning_messages <- character()
    fitted_policy <- NULL
    row <- data.frame(
        setting = specification$setting, condition = specification$condition, method = method,
        seed = as.integer(seed), status = "failed", regret = NA_real_, test_MCSE = NA_real_,
        seconds = NA_real_, fit_seconds = NA_real_, evaluation_seconds = NA_real_,
        converged_starts = NA_integer_, starts_total = 11L,
        chosen_start = NA_integer_, chosen_status = "", chosen_gradient = NA_real_,
        retried_starts = NA_integer_, total_dc_iterations = NA_integer_,
        rank_deficient_fits = NA_integer_, bfgs_calls = NA_integer_, bfgs_nonconvergence = NA_integer_,
        optimization_gap_lower_bound = NA_real_, tau_rmse = NA_real_, c_rmse = NA_real_,
        error = "", warnings = "", stringsAsFactors = FALSE
    )
    tryCatch(withCallingHandlers({
        # Match the independent cluster job's RNG and cache state before each fit.
        rm(list = ls(.policy_cache$entries, all.names = TRUE), envir = .policy_cache$entries)
        set.seed(900000L + seed)
        fit_started <- proc.time()[["elapsed"]]
        fit <- fitter(sample, method)
        row$fit_seconds <- unname(proc.time()[["elapsed"]] - fit_started)
        evaluation_started <- proc.time()[["elapsed"]]
        if (is.null(evaluation)) {
            evaluation <- list(effect = model$tau(sample$test), policy_basis = model$b(sample$test))
        }
        true_effect <- evaluation$effect
        policy_score <- drop(evaluation$policy_basis %*% fit$theta)
        loss <- abs(true_effect) * ((policy_score > 0) != (true_effect > 0))
        stopifnot(all(is.finite(loss)))
        row$status <- if (fit$converged_starts == fit$starts_total) "completed" else "optimizer_warning"
        row$regret <- mean(loss)
        row$test_MCSE <- sd(loss) / sqrt(length(loss))
        row$evaluation_seconds <- unname(proc.time()[["elapsed"]] - evaluation_started)
        diagnostic_names <- intersect(names(row), names(fit))
        for (name in diagnostic_names) row[[name]] <- fit[[name]]
        # Only coefficients are needed for the first seed's decision boundaries.
        fitted_policy <- list(theta = fit$theta)
    }, warning = function(warning) {
        warning_messages <<- c(warning_messages, conditionMessage(warning))
        invokeRestart("muffleWarning")
    }), error = function(error) {
        row$status <<- "failed"
        row$regret <<- NA_real_
        row$test_MCSE <<- NA_real_
        row$error <<- conditionMessage(error)
        fitted_policy <<- NULL
    })
    row$seconds <- unname(proc.time()[["elapsed"]] - started)
    row$warnings <- paste(unique(warning_messages), collapse = " | ")
    list(row = row, fit = fitted_policy)
}

run_simulation <- function(configuration, project_directory) {
    validate_configuration(configuration)
    # Each score is solved once; hashing a large matrix cannot yield a cache hit.
    previous_options <- options(itr.solver_cache = FALSE)
    on.exit(options(previous_options), add = TRUE)
    required_packages <- c("Rcpp", "randomForest")
    missing <- required_packages[!vapply(required_packages, requireNamespace,
                                        logical(1), quietly = TRUE)]
    if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))

    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1,
               VECLIB_MAXIMUM_THREADS = 1)
    cores <- if (.Platform$OS.type == "unix") as.integer(configuration$cores) else 1L
    seeds <- as.integer(configuration$start_seed + seq_len(configuration$replications) - 1L)
    specification <- setting_specification(configuration$setting, configuration$condition)
    case_name <- paste(configuration$setting, configuration$condition, sep = "/")
    output_directory <- configuration$output_directory
    code_directory <- file.path(project_directory, "methods")
    code_files <- sort(list.files(code_directory, pattern = "\\.(R|cpp)$", full.names = TRUE))
    fingerprint <- list(
        code = setNames(unname(tools::md5sum(code_files)), basename(code_files)),
        R = R.version.string,
        Rcpp = as.character(packageVersion("Rcpp")),
        randomForest = as.character(packageVersion("randomForest")),
        rng = RNGkind(), platform = R.version$platform,
        BLAS = sessionInfo()$BLAS, LAPACK = sessionInfo()$LAPACK,
        optimizer = optimizer_control(), starts = 11L, cores = cores
    )

    # Refuse to combine checkpoints from a different model, source, or runtime.
    checkpoint_path <- file.path(output_directory, "checkpoint.rds")
    signature <- list(specification = specification, fingerprint = fingerprint)
    if (file.exists(checkpoint_path)) {
        checkpoint <- readRDS(checkpoint_path)
        if (!identical(checkpoint$signature, signature)) {
            stop("Incompatible checkpoint in ", case_name,
                 ". Change output_folder in setup.R or move this condition's old output folder.")
        }
    } else {
        checkpoint <- list(signature = signature, results = empty_results(),
                           boundary_seed = configuration$start_seed, boundary_fits = list())
    }
    if (checkpoint$boundary_seed != configuration$start_seed) {
        checkpoint$boundary_seed <- configuration$start_seed
        checkpoint$boundary_fits <- list()
    }

    timings <- data.frame(phase = character(), seed = integer(), method = character(), seconds = numeric())
    record_timing <- function(phase, started, seed = NA_integer_, method = "") {
        elapsed <- unname(proc.time()[["elapsed"]] - started)
        timings <<- rbind(timings, data.frame(phase = phase, seed = seed, method = method, seconds = elapsed))
        message(sprintf("%s: %.2f s", phase, elapsed))
    }
    message(case_name, ": loading C++ backend (first build can take longer).")
    started <- proc.time()[["elapsed"]]
    backend <- load_fast_backend(code_directory)
    record_timing("backend_load", started)
    started <- proc.time()[["elapsed"]]
    model <- make_model(specification)
    fitter <- make_fitter(model, backend, cores)
    record_timing("model_setup", started)
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    write.csv(specification, file.path(output_directory, "settings.csv"), row.names = FALSE)
    saveRDS(list(configuration = configuration, seeds = seeds, fingerprint = fingerprint),
            file.path(output_directory, "run_config.rds"))
    writeLines(c(capture.output(dput(configuration)), "", capture.output(sessionInfo())),
               file.path(output_directory, "session.txt"))
    expected <- length(configuration$methods) * length(seeds)
    message(case_name, ": ", length(configuration$methods), " methods x ",
            length(seeds), " seeds = ", expected, " fits")
    message("Decision boundaries use seed ", configuration$start_seed, " only.")
    writeLines("running", file.path(output_directory, "status.txt"))
    plot_errors <- character()

    for (seed in seeds) {
        sample <- NULL
        evaluation <- NULL
        for (method in configuration$methods) {
            previous <- which(checkpoint$results$seed == seed & checkpoint$results$method == method)
            reusable <- length(previous) == 1L && checkpoint$results$status[previous] == "completed"
            needs_boundary <- seed == configuration$start_seed && is.null(checkpoint$boundary_fits[[method]])
            if (reusable && !needs_boundary) {
                message(case_name, " | seed ", seed, " | ", method, ": reused checkpoint.")
                next
            }
            if (is.null(sample)) {
                message(case_name, " | seed ", seed, ": preparing data and test basis.")
                started <- proc.time()[["elapsed"]]
                sample <- generate_data(specification, model, seed)
                evaluation <- list(effect = model$tau(sample$test), policy_basis = model$b(sample$test))
                record_timing("data_and_test_basis", started, seed)
            }
            message(case_name, " | seed ", seed, " | ", method, ": fitting...")
            started <- proc.time()[["elapsed"]]
            result <- fit_and_evaluate(sample, specification, model, method, seed, fitter, evaluation)
            record_timing("fit_and_evaluate", started, seed, method)
            if (length(previous)) checkpoint$results <- checkpoint$results[-previous, , drop = FALSE]
            checkpoint$results <- rbind(checkpoint$results, result$row)
            if (seed == configuration$start_seed) checkpoint$boundary_fits[[method]] <- result$fit
            save_checkpoint(checkpoint, checkpoint_path)
            print(result$row[, c("setting", "condition", "method", "seed", "status", "regret",
                                "fit_seconds", "evaluation_seconds", "retried_starts")],
                  row.names = FALSE)
        }
        # Keep regret records for every seed; keep policies only for the first.
        sample <- NULL
        evaluation <- NULL
        selected <- checkpoint$results$seed %in% seeds & checkpoint$results$method %in% configuration$methods
        current_results <- checkpoint$results[selected, , drop = FALSE]
        write_result_tables(current_results, output_directory, length(seeds), configuration$methods)
        if (seed == configuration$start_seed) {
            message(case_name, ": drawing decision boundaries.")
            started <- proc.time()[["elapsed"]]
            tryCatch(plot_decision_boundaries(
                model, specification, seed, configuration$methods,
                checkpoint$boundary_fits, current_results,
                file.path(output_directory, "decision_boundary.png")
            ), error = function(error) {
                plot_errors <<- c(plot_errors, paste(case_name, conditionMessage(error)))
            })
            record_timing("decision_boundary_plot", started, seed)
        }
        write.csv(timings, file.path(output_directory, "timing.csv"), row.names = FALSE)
    }

    numeric_success <- sum(is.finite(current_results$regret) &
                           current_results$status %in% c("completed", "optimizer_warning"))
    status <- if (nrow(current_results) == expected && all(current_results$status == "completed") &&
                  !length(plot_errors)) "completed" else if (numeric_success == 0L) "failed" else "mixed"
    report <- c(
        paste("Status:", status), paste("Case:", case_name), paste("Expected fits:", expected),
        paste("Completed:", sum(current_results$status == "completed")),
        paste("Optimizer warnings:", sum(current_results$status == "optimizer_warning")),
        paste("Failed:", sum(current_results$status == "failed")),
        paste("Plot errors:", length(plot_errors)), plot_errors,
        paste("Boundary seed:", configuration$start_seed),
        paste("Output:", normalizePath(output_directory)),
        "Stationarity does not certify a global optimum."
    )
    writeLines(status, file.path(output_directory, "status.txt"))
    writeLines(report, file.path(output_directory, "completion_report.txt"))
    cat(paste(report, collapse = "\n"), "\n")
    if (status != "completed") stop("Run finished with issues; inspect completion_report.txt and regret.csv.")
    invisible(current_results)
}
