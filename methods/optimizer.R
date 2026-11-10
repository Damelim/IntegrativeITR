# Newton/DC policy optimization with the original stopping rules and 11 starts.
# Z: policy basis matrix; r: policy score; k: compiled loss/gradient functions.
# Damping stabilizes Newton steps; it is not a statistical penalty.

.policy_cache <- new.env(parent = emptyenv())

.policy_cache$entries <- new.env(parent = emptyenv())

.policy_cache$last_hit <- FALSE

optimizer_control <- function(max_outer = 1000L, max_inner = 200L, inner_tol = 1e-07, outer_tol = 1e-06, 
    objective_tol = 1e-09, bfgs_fallback = 2000L) {
    stopifnot(max_outer >= 1L, max_inner >= 1L, inner_tol > 0, outer_tol > 0, objective_tol > 0, bfgs_fallback >= 
        0L)
    list(max_outer = max_outer, max_inner = max_inner, inner_tol = inner_tol, outer_tol = outer_tol, 
        objective_tol = objective_tol, bfgs_fallback = bfgs_fallback)
}

convex_newton <- function(Q, w, label, beta, cg, k, ctl) {
    margin <- function(b) label * drop(Q %*% b)
    fn <- function(b) 0.5 * k$itr_loss(w, margin(b), FALSE) - sum(cg * b)
    gr <- function(b) 0.5 * k$itr_gradient(Q, w, label, margin(b), FALSE) - cg
    value <- fn(beta)
    fallback <- FALSE
    line_failure <- FALSE
    for (iter in seq_len(ctl$max_inner)) {
        u <- margin(beta)
        g <- 0.5 * k$itr_gradient(Q, w, label, u, FALSE) - cg
        if (max(abs(g)) <= ctl$inner_tol) 
            break
        active <- u >= 0 & u <= 1
        qa <- Q[active, , drop = FALSE] * sqrt(w[active])
        H <- crossprod(qa)
        damping <- 1e-08 * max(1, max(diag(H)))
        C <- tryCatch(chol(H + diag(damping, ncol(Q))), error = function(e) NULL)
        direction <- if (is.null(C)) 
            -g
        else -backsolve(C, forwardsolve(t(C), g))
        slope <- sum(g * direction)
        if (!is.finite(slope) || slope >= 0) {
            direction <- -g
            slope <- -sum(g * g)
        }
        step <- 1
        accepted <- FALSE
        for (ls in 1:40) {
            candidate <- beta + step * direction
            new_value <- fn(candidate)
            if (is.finite(new_value) && new_value <= value + 1e-04 * step * slope + 1e-14) {
                accepted <- TRUE
                break
            }
            step <- step/2
        }
        if (!accepted) {
            line_failure <- TRUE
            break
        }
        beta <- candidate
        value <- new_value
    }
    gradient <- max(abs(gr(beta)))
    fallback_code <- NA_integer_
    if (gradient > ctl$inner_tol && ctl$bfgs_fallback > 0L) {
        fallback <- TRUE
        fit <- stats::optim(beta, fn, gr, method = "BFGS", control = list(maxit = ctl$bfgs_fallback, 
            reltol = 1e-13))
        if (is.finite(fit$value) && fit$value <= value + 1e-12) {
            beta <- fit$par
            value <- fit$value
        }
        gradient <- max(abs(gr(beta)))
        fallback_code <- fit$convergence
    }
    list(beta = beta, value = value, gradient = gradient, converged = gradient <= ctl$inner_tol, iterations = iter, 
        fallback = fallback, fallback_code = fallback_code, line_failure = line_failure)
}

# These matrices are identical for all starts and all retries of one policy fit.
prepare_policy_problem <- function(Z, r) {
    keep <- is.finite(r) & abs(r) > 1e-12
    z <- as.matrix(Z[keep, , drop = FALSE])
    rr <- r[keep]
    if (!length(rr)) return(NULL)
    w <- abs(rr)
    w <- w/sum(w)
    label <- sign(rr)
    scale <- sqrt(colSums(z * z * w))
    scale[!is.finite(scale) | scale < 1e-12] <- 1
    Q <- sweep(z, 2, scale, "/")
    list(Q = Q, weight = w, label = label, scale = scale)
}

solve_policy_start <- function(Z, r, start, k, ctl = optimizer_control(),
                               prepared = prepare_policy_problem(Z, r)) {
    if (is.null(prepared))
        return(list(theta = start, trace = data.frame(), status = "zero_score", terminal_gradient = 0,
            inner_failures = 0L, bfgs_calls = 0L, bfgs_failures = 0L))
    Q <- prepared$Q
    w <- prepared$weight
    label <- prepared$label
    scale <- prepared$scale
    beta <- as.numeric(start) * scale
    margin <- function(b) label * drop(Q %*% b)
    full <- function(b) k$itr_loss(w, margin(b), TRUE)
    gradient <- function(b) {
        u <- margin(b)
        0.5 * (k$itr_gradient(Q, w, label, u, FALSE) - k$itr_gradient(Q, w, label, u, TRUE))
    }
    value <- full(beta)
    trace <- vector("list", ctl$max_outer)
    status <- "outer_limit"
    failures <- bf <- bff <- 0L
    for (iter in seq_len(ctl$max_outer)) {
        cg <- 0.5 * k$itr_gradient(Q, w, label, margin(beta), TRUE)
        inner <- convex_newton(Q, w, label, beta, cg, k, ctl)
        bf <- bf + as.integer(inner$fallback)
        bff <- bff + as.integer(inner$fallback && !isTRUE(inner$fallback_code == 0L))
        candidate <- inner$beta
        new_value <- full(candidate)
        if (new_value > value + 1e-10) {
            status <- "descent_failure"
            candidate <- beta
            new_value <- value
        }
        g <- max(abs(gradient(candidate)))
        trace[[iter]] <- data.frame(iteration = iter, ramp_loss = new_value, empirical_value = sum(r * 
            (drop(Z %*% (candidate/scale)) > 0)), gradient_inf = g, inner_gradient_inf = inner$gradient, 
            inner_iterations = inner$iterations, inner_converged = inner$converged, bfgs_fallback = inner$fallback, 
            coefficient_norm = sqrt(sum((candidate/scale)^2)))
        stable <- abs(new_value - value) <= ctl$objective_tol * max(1, abs(value))
        beta <- candidate
        value <- new_value
        if (status == "descent_failure") 
            break
        if (!inner$converged) {
            failures <- failures + 1L
            status <- "inner_not_converged"
            break
        }
        if (stable && g <= ctl$outer_tol) {
            status <- "stationary"
            break
        }
    }
    list(theta = beta/scale, trace = do.call(rbind, trace[seq_len(iter)]), status = status, terminal_gradient = g, 
        inner_failures = failures, bfgs_calls = bf, bfgs_failures = bff)
}

solve_policy <- function(Z, r, k, seed, n_starts = 11L, source_start = NULL, ctl = optimizer_control(), 
    ls_start = TRUE, cores = 1L) {
    stopifnot(all(is.finite(Z)), all(is.finite(r)), n_starts > as.integer(ls_start))
    set.seed(seed)
    starts <- list()
    if (!is.null(source_start)) 
        starts[[1L]] <- source_start
    starts[[length(starts) + 1L]] <- rep(0, ncol(Z))
    while (length(starts) < n_starts - as.integer(ls_start)) starts[[length(starts) + 1L]] <- rnorm(ncol(Z), 
        sd = 0.35)
    if (ls_start) {
        use <- abs(r) > 1e-12
        theta <- if (any(use)) 
            lm.fit(Z[use, , drop = FALSE], r[use])$coefficients
        else rep(0, ncol(Z))
        theta[!is.finite(theta)] <- 0
        magnitude <- if (any(use)) 
            mean(abs(drop(Z[use, , drop = FALSE] %*% theta)))
        else 0
        if (is.finite(magnitude) && magnitude > 0) 
            theta <- theta/magnitude
        starts[[length(starts) + 1L]] <- as.numeric(theta)
    }
    .policy_cache$last_hit <- FALSE
    cache_enabled <- isTRUE(getOption("itr.solver_cache", TRUE))
    if (cache_enabled) {
        path <- tempfile("policy_key_")
        on.exit(unlink(path), add = TRUE)
        signature <- unname(tools::md5sum(file.path(.opt_dir, c("optimizer.R", "fast_kernels.cpp", "fast_backend.R"))))
        writeBin(serialize(list(Z = Z, r = r, starts = starts, control = ctl, source = signature, R = R.version.string), 
            NULL, version = 3L), path)
        cache_key <- unname(tools::md5sum(path))
        if (exists(cache_key, .policy_cache$entries, inherits = FALSE)) {
            .policy_cache$last_hit <- TRUE
            return(get(cache_key, .policy_cache$entries, inherits = FALSE))
        }
    }
    cores <- suppressWarnings(as.integer(cores))
    if (length(cores) != 1L || is.na(cores) || cores < 1L) 
        cores <- 1L
    if (.Platform$OS.type != "unix") 
        cores <- 1L
    cores <- min(cores, n_starts)
    prepared <- prepare_policy_problem(Z, r)
    task <- function(start) {
        fit <- solve_policy_start(Z, r, start, k, ctl, prepared)
        if (fit$status %in% c("inner_not_converged", "outer_limit")) {
            extended <- ctl
            extended$max_inner <- max(1000L, 5L * ctl$max_inner)
            extended$bfgs_fallback <- max(10000L, 5L * ctl$bfgs_fallback)
            if (fit$status == "outer_limit") 
                extended$max_outer <- max(2000L, 2L * ctl$max_outer)
            retry <- solve_policy_start(Z, r, fit$theta, k, extended, prepared)
            retry$initial_attempt <- fit
            retry$bfgs_calls <- retry$bfgs_calls + fit$bfgs_calls
            retry$bfgs_failures <- retry$bfgs_failures + fit$bfgs_failures
            retry$inner_failures <- retry$inner_failures + fit$inner_failures
            fit <- retry
        }
        fit
    }
    fits <- if (cores > 1L) 
        parallel::mclapply(starts[seq_len(n_starts)], task, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = FALSE)
    else lapply(starts[seq_len(n_starts)], task)
    if (any(vapply(fits, function(x) inherits(x, "try-error") || is.null(x), logical(1)))) 
        stop("A policy-start worker failed; no partial result is selected.")
    values <- vapply(fits, function(f) sum(r * (drop(Z %*% f$theta) > 0)), numeric(1))
    chosen <- which.max(values)
    result <- list(theta = fits[[chosen]]$theta, chosen = chosen, empirical_values = values, starts = fits, 
        converged = vapply(fits, function(f) f$status %in% c("stationary", "zero_score"), logical(1)))
    if (cache_enabled) 
        assign(cache_key, result, .policy_cache$entries)
    result
}

