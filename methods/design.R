# All settings and perturbations reported in the paper.
paper_settings <- c("setting1", "setting2")
paper_setting_labels <- c(
    setting1 = "Setting 1 (constant confounding)",
    setting2 = "Setting 2 (quadratic confounding)"
)
paper_conditions <- c(
    "baseline", "target_large", "trial_large", "no_shift",
    "strong_shift", "tau_underfit", "c_full", "policy_full"
)
paper_methods <- c("WY", "Target", "Pooled", "Wu2020", "IP", "CC", "CC_DRGO")

setting_specification <- function(setting, condition) {
    stopifnot(setting %in% paper_settings, condition %in% paper_conditions)
    policy_degree <- 6L
    effect_degree <- 6L
    confounding_degree <- if (setting == "setting1") 0L else 2L
    target_size <- 8000L
    trial_size <- 800L
    shift <- 0.75

    if (condition == "target_large") target_size <- 16000L
    if (condition == "trial_large") trial_size <- 1600L
    if (condition == "no_shift") shift <- 0
    if (condition == "strong_shift") shift <- 2
    if (condition == "tau_underfit") effect_degree <- 2L
    if (condition == "c_full") confounding_degree <- 6L
    if (condition == "policy_full") policy_degree <- 7L

    basis_dimension <- function(degree) as.integer((degree + 1) * (degree + 2) / 2)
    data.frame(
        setting = setting, condition = condition,
        n0 = target_size, n1 = trial_size, shift = shift, sigma = 1 / 3,
        b_degree = policy_degree, p_degree = effect_degree,
        z_degree = confounding_degree,
        d_D = basis_dimension(policy_degree), d_tau = basis_dimension(effect_degree),
        d_c = basis_dimension(confounding_degree), d_mu = 28L, d_arm = 28L
    )
}

# Basis column order and arithmetic match the completed cluster study.

basis_degree <- function(degree) {
    force(degree)
    function(X) {
        one <- function(x) {
            a <- matrix(1, length(x), degree + 1L)
            if (degree >= 1L) 
                a[, 2L] <- x
            if (degree >= 2L) 
                for (j in 2:degree) a[, j + 1L] <- ((2 * j - 1) * x * a[, j] - (j - 1) * a[, j - 1L])/j
            sweep(a, 2L, sqrt(2 * (0:degree) + 1), "*")
        }
        a <- one(X$x1)
        b <- one(X$x2)
        ij <- do.call(rbind, lapply(0:degree, function(j) cbind(j:0, 0:j)))
        ans <- vapply(seq_len(nrow(ij)), function(j) a[, ij[j, 1] + 1L] * b[, ij[j, 2] + 1L], numeric(nrow(X)))
        matrix(ans, nrow(X), nrow(ij))
    }
}

raw_tau6 <- function(X) {
    r <- X$x1^2 + X$x2^2
    (r - 0.2) * (r - 0.65) * (r - 1.2)
}

generate_cell <- function(spec, model, replication, ntest = 100000L) {
    stream <- function(offset, n, kind = "uniform") {
        set.seed(as.integer(replication + offset))
        if (kind == "normal") 
            rnorm(n)
        else runif(n)
    }
    draw <- function(n, s, offset) {
        u1 <- stream(offset + 10000L, n)
        u2 <- stream(offset + 20000L, n)
        transform <- function(u) {
            if (s == 0L || spec$shift == 0) 
                2 * u - 1
            else {
                sd <- sqrt(1/(2 * spec$shift))
                sd * qnorm(pnorm(-1/sd) + u * (pnorm(1/sd) - pnorm(-1/sd)))
            }
        }
        X <- data.frame(x1 = transform(u1), x2 = transform(u2))
        U <- ifelse(stream(offset + 30000L, n) < 0.5, -1, 1)
        A <- if (s == 0L) 
            as.integer(U == 1)
        else as.integer(stream(offset + 40000L, n) < 0.5)
        tau <- model$tau(X)
        mu <- X$x1 + X$x2
        Y <- mu + (A - 0.5) * tau + (1 - s) * U * model$c(X)/2 + spec$sigma * stream(offset + 50000L, 
            n, "normal")
        data.frame(X, A, Y)
    }
    test <- data.frame(x1 = 2 * stream(700000L, ntest) - 1, x2 = 2 * stream(800000L, ntest) - 1)
    list(target = draw(spec$n0, 0L, 100000L), trial = draw(spec$n1, 1L, 300000L), test = test, setting = spec$setting, 
        condition = spec$perturbation, rep = replication, seed = as.integer(replication + 1L), sigma = spec$sigma)
}

# The calibration grid is fixed and independent of all simulation seeds.
make_model <- function(specification) {
    calibration_grid <- expand.grid(
        x1 = seq(-1, 1, length.out = 601),
        x2 = seq(-1, 1, length.out = 601)
    )
    effect_scale <- 1 / IQR(raw_tau6(calibration_grid))
    model <- list(
        example = data.frame(x1 = 0, x2 = 0),
        covariates = c("x1", "x2"),
        weight_formula = S ~ x1 + x2 + I(x1^2) + I(x1 * x2) + I(x2^2),
        marginal_mean_model = "full"
    )
    model$tau <- function(X) effect_scale * raw_tau6(X)
    model$c <- if (specification$setting == "setting1") {
        function(X) rep(2, nrow(X))
    } else {
        function(X) with(X, 2 + x1 - x2 + x1^2 - x2^2)
    }
    model$b <- basis_degree(specification$b_degree)
    model$p <- basis_degree(specification$p_degree)
    model$z <- basis_degree(specification$z_degree)
    model$marginal_basis <- basis_degree(6L)
    model$arm_basis <- basis_degree(6L)
    model$variance_basis <- basis_degree(2L)
    model$sigma <- specification$sigma
    model
}

make_plan <- function(settings, conditions) {
    rows <- lapply(settings, function(setting) {
        do.call(rbind, lapply(conditions, function(condition) {
            setting_specification(setting, condition)
        }))
    })
    plan <- do.call(rbind, rows)
    rownames(plan) <- NULL
    plan
}

# The cluster generator uses the name 'perturbation' internally.
generate_data <- function(specification, model, seed) {
    generator_specification <- specification
    generator_specification$perturbation <- specification$condition
    generate_cell(generator_specification, model, replication = seed, ntest = 100000L)
}
