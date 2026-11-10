# Policy scores and cross-fitted nuisance estimation.
# b: policy basis; p: effect basis; z: confounding basis.
# Active score calculations are retained from the 50-replication study.
# fitting.R supplies 11 starts and the current optimizer to every method.
# Equal-size folds are randomized; their historical helper name is retained.

make_method_kernel <- function(method, model) {
    make_basis <- model$p
    make_policy_basis <- model$b
    make_confounding_basis <- make_pointwise_confounding_basis <- model$z
    make_variance_basis <- model$variance_basis
    diagnostics <- new.env(parent = emptyenv())
    diagnostics$starts <- diagnostics$optimizer_calls <- diagnostics$optimizer_nonconvergence <- 0L
    diagnostics$rank_deficient_fits <- diagnostics$least_squares_fits <- 0L
    if (method %in% c("CC")) 
        fit_fixed_alpha_path <- function(data, alpha_grid = seq(0, 0, by = 0), lambda_tau = 0, lambda_c = 0, 
            n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) {
            scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = TRUE, 
                need_weight = any(alpha_grid != 0))
            source = scores$samples$source
            target = scores$samples$target
            pointwise_c = fit_cross_fitted_pointwise_confounding(scores, lambda_tau = lambda_tau, lambda_c = lambda_c)
            eta_source = NULL
            if (any(alpha_grid == 1)) {
                source_objective = make_alpha_objective(source, target, scores, pointwise_c, alpha = 1)
                eta_source = dc_policy(source_objective$Z, source_objective$r, n_starts = n_starts, max_dc = max_dc, 
                  seed = seed + 1)
            }
            eta_by_alpha = vector("list", length(alpha_grid))
            names(eta_by_alpha) = paste0("alpha_", sprintf("%.1f", alpha_grid))
            for (j in seq_along(alpha_grid)) {
                alpha = alpha_grid[j]
                if (abs(alpha - 1) < 1e-12) {
                  eta_by_alpha[[j]] = eta_source
                }
                else {
                  objective = make_alpha_objective(source, target, scores, pointwise_c, alpha)
                  eta_by_alpha[[j]] = dc_policy(objective$Z, objective$r, n_starts = n_starts, max_dc = max_dc, 
                    seed = seed + 1)
                }
            }
            list(alpha_grid = alpha_grid, eta_by_alpha = eta_by_alpha, eta_source = eta_source, pointwise_c = pointwise_c, 
                scores = scores)
        }
    fit_integrative_plugin <- function(data, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) {
        samples = make_samples(data)
        source = samples$source
        target = samples$target
        folds = make_crossfit_folds(source, target, n_folds, seed)
        tau_target = rep(NA_real_, nrow(target))
        tau_coef_by_fold = matrix(NA_real_, n_folds, ncol(make_basis(model$example)))
        c_coef_by_fold = matrix(NA_real_, n_folds, ncol(make_confounding_basis(model$example)))
        for (k in seq_len(n_folds)) {
            tr1 = source[folds$source != k, , drop = FALSE]
            tr0 = target[folds$target != k, , drop = FALSE]
            ev0 = which(folds$target == k)
            inner = fit_inner_r_predictions(tr1, tr0, n_inner_folds = 5, seed = 20000 + k)
            joint = fit_integrative_r_fold(tr1, tr0, nuisance = NULL, lambda_tau = 0, lambda_c = 0, inner_predictions = inner)
            tau_coef_by_fold[k, ] = joint$tau_coef
            c_coef_by_fold[k, ] = joint$c_coef
            tau_target[ev0] = as.vector(make_basis(target[ev0, , drop = FALSE]) %*% joint$tau_coef)
        }
        if (any(!is.finite(tau_target))) 
            stop("Nonfinite out-of-fold treatment-effect prediction.")
        Z = make_policy_basis(target)
        r = tau_target/nrow(target)
        eta_hat = dc_policy(Z, r, n_starts = n_starts, max_dc = max_dc, seed = seed + 1)
        list(eta_hat = eta_hat, tau_target = tau_target, tau_coef_by_fold = tau_coef_by_fold, c_coef_by_fold = c_coef_by_fold, 
            folds = folds, empirical_value = empirical_value(Z, r, eta_hat))
    }
    fit_naive_pooled_itr <- function(data, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) {
        score = fit_cross_fitted_pooled_score(data, n_folds, seed)
        Z = make_policy_basis(score$pooled)
        r = score$Gamma_pool/nrow(score$pooled)
        eta_hat = dc_policy(Z, r, n_starts = n_starts, max_dc = max_dc, seed = seed + 1)
        list(eta_hat = eta_hat, association_coef = score$association_coef, pooled_empirical_value = empirical_value(Z, 
            r, eta_hat))
    }
    fit_wu2020 <- function(data, ntree = 500L, n_folds = 5L, n_starts = 8L, max_dc = 100L, seed = 1L) {
        target <- wu2020_sample(data, "rwd")
        rct <- wu2020_sample(data, "rct")
        folds <- wu2020_target_folds(nrow(target), n_folds, seed + 1L)
        tau_score_oof <- rep(NA_real_, nrow(target))
        q_rank_by_fold <- matrix(NA_integer_, n_folds, 2L, dimnames = list(NULL, c("A0", "A1")))
        for (ell in seq_len(n_folds)) {
            train <- target[folds != ell, , drop = FALSE]
            eval <- target[folds == ell, , drop = FALSE]
            feature_fit <- wu2020_fit_features(train, ntree, seed + 1000L * ell)
            h_rct <- wu2020_predict_features(feature_fit, rct)
            z_rct <- cbind(model$arm_basis(rct), h_rct)
            q_fit <- wu2020_q_fit(z_rct, rct$A, rct$Y)
            h_eval <- wu2020_predict_features(feature_fit, eval)
            z_eval <- cbind(model$arm_basis(eval), h_eval)
            tau_score_oof[folds == ell] <- as.vector(z_eval %*% q_fit$contrast)
            q_rank_by_fold[ell, ] <- q_fit$rank
        }
        if (any(!is.finite(tau_score_oof))) 
            stop("Nonfinite cross-fitted Wu2020 score.")
        policy_basis <- make_policy_basis(target)
        eta_hat <- dc_policy(Z = policy_basis, r = tau_score_oof/nrow(target), n_starts = n_starts, max_dc = max_dc, 
            seed = seed + 90000L)
        list(eta_hat = eta_hat, basis_fun = make_policy_basis, tau_score_oof = tau_score_oof, target_folds = folds, 
            q_rank_by_fold = q_rank_by_fold, projected_empirical_value = empirical_value(policy_basis, 
                tau_score_oof/nrow(target), eta_hat))
    }
    fit_wuyang <- function(data, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) {
        scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = FALSE)
        source = scores$samples$source
        Z1 = make_policy_basis(source)
        r1 = scores$source$w * scores$source$Gamma1_hat/nrow(source)
        eta_hat = dc_policy(Z1, r1, n_starts = n_starts, max_dc = max_dc, seed = seed + 1)
        list(eta_hat = eta_hat, scores = scores, source_empirical_value = empirical_value(Z1, r1, eta_hat))
    }
    make_alpha_objective <- function(source, target, scores, pointwise_c, alpha) {
        Z0 = make_policy_basis(target)
        Z1 = make_policy_basis(source)
        r0 = (1 - alpha) * (scores$target$Gamma0_hat - pointwise_c$c_target)/nrow(target)
        r1 = alpha * scores$source$w * scores$source$Gamma1_hat/nrow(source)
        list(Z = rbind(Z0, Z1), r = c(r0, r1))
    }
    clip <- function(x, lower = 0.02, upper = 0.98) {
        pmin(pmax(x, lower), upper)
    }
    fit_basis_regression <- function(y, basis) {
        fit <- lm.fit(x = as.matrix(basis), y = y)
        diagnostics$least_squares_fits <- diagnostics$least_squares_fits + 1L
        if (fit$rank < ncol(basis)) 
            diagnostics$rank_deficient_fits <- diagnostics$rank_deficient_fits + 1L
        coefficient <- as.numeric(fit$coefficients)
        coefficient[!is.finite(coefficient)] <- 0
        list(coefficient = coefficient)
    }
    fit_cross_fitted_pooled_score <- function(data, n_folds = 5, seed = 1) {
        pooled = make_pooled_sample(data)
        fold = make_stratified_folds(pooled$A, n_folds, seed)
        gamma_pool = rep(NA_real_, nrow(pooled))
        association_coef_by_fold = matrix(NA_real_, nrow = n_folds, ncol = ncol(model$arm_basis(model$example)))
        for (k in seq_len(n_folds)) {
            train = pooled[fold != k, , drop = FALSE]
            test = pooled[fold == k, , drop = FALSE]
            nuisance = fit_pooled_nuisance(train)
            gamma_pool[fold == k] = predict_pooled_score(nuisance, test)
            association_coef_by_fold[k, ] = nuisance$mu1_coef - nuisance$mu0_coef
        }
        list(pooled = pooled, Gamma_pool = gamma_pool, association_coef = colMeans(association_coef_by_fold))
    }
    if (method %in% c("WY")) 
        fit_cross_fitted_scores <- function(data, n_folds = 5, seed = 1, need_target = TRUE) {
            samples = make_samples(data)
            source = samples$source
            target = samples$target
            folds = make_crossfit_folds(source, target, n_folds, seed)
            source_score = data.frame(Gamma1_hat = rep(NA_real_, nrow(source)), tau1_hat = rep(NA_real_, 
                nrow(source)), m1_hat = rep(NA_real_, nrow(source)), w = rep(NA_real_, nrow(source)))
            if (need_target) {
                source_score$omega0_source_hat = NA_real_
                target_score = data.frame(Gamma0_hat = rep(NA_real_, nrow(target)), omega0_hat = rep(NA_real_, 
                  nrow(target)), m0_hat = rep(NA_real_, nrow(target)), e0_hat = rep(NA_real_, nrow(target)), 
                  w = rep(NA_real_, nrow(target)))
            }
            else {
                target_score = NULL
            }
            fold_fits = vector("list", n_folds)
            for (k in seq_len(n_folds)) {
                source_train_index = which(folds$source != k)
                target_train_index = which(folds$target != k)
                source_eval_index = which(folds$source == k)
                target_eval_index = which(folds$target == k)
                source_train = source[source_train_index, , drop = FALSE]
                target_train = target[target_train_index, , drop = FALSE]
                source_eval = source[source_eval_index, , drop = FALSE]
                target_eval = target[target_eval_index, , drop = FALSE]
                weight_fit = fit_transfer_weight_model(source_train, target_train)
                nuisance_fit = fit_nuisance_models(source_train, target_train, need_target = need_target)
                evaluation = predict_nuisance_scores(nuisance_fit, weight_fit, source_eval, target_eval, 
                  need_target = need_target)
                source_score[source_eval_index, names(evaluation$source)] = evaluation$source
                if (need_target) 
                  target_score[target_eval_index, ] = evaluation$target
                fold_fits[[k]] = list(nuisance = nuisance_fit, weight = weight_fit, source_train_index = source_train_index, 
                  target_train_index = target_train_index, source_eval_index = source_eval_index, target_eval_index = target_eval_index)
            }
            source_score$fold = folds$source
            if (need_target) 
                target_score$fold = folds$target
            model_names = c("mu11", "mu10")
            if (need_target) 
                model_names = c(model_names, "mu01", "mu00")
            averaged_models = setNames(lapply(model_names, function(name) average_fold_model(fold_fits, 
                name)), model_names)
            list(source = source_score, target = target_score, models = averaged_models, fold_fits = fold_fits, 
                folds = folds, samples = samples, n_folds = n_folds)
        }
    if (method %in% c("CC", "Fusion")) 
        fit_cross_fitted_scores <- function(data, n_folds = 5, seed = 1, need_target = TRUE, need_weight = TRUE) {
            samples = make_samples(data)
            source = samples$source
            target = samples$target
            folds = make_crossfit_folds(source, target, n_folds, seed)
            source_score = data.frame(Gamma1_hat = rep(NA_real_, nrow(source)), tau1_hat = rep(NA_real_, 
                nrow(source)), m1_hat = rep(NA_real_, nrow(source)), w = rep(NA_real_, nrow(source)))
            if (need_target) {
                source_score$omega0_source_hat = NA_real_
                target_score = data.frame(Gamma0_hat = rep(NA_real_, nrow(target)), omega0_hat = rep(NA_real_, 
                  nrow(target)), m0_hat = rep(NA_real_, nrow(target)), e0_hat = rep(NA_real_, nrow(target)), 
                  w = rep(NA_real_, nrow(target)))
            }
            else {
                target_score = NULL
            }
            fold_fits = vector("list", n_folds)
            for (k in seq_len(n_folds)) {
                source_train_index = which(folds$source != k)
                target_train_index = which(folds$target != k)
                source_eval_index = which(folds$source == k)
                target_eval_index = which(folds$target == k)
                source_train = source[source_train_index, , drop = FALSE]
                target_train = target[target_train_index, , drop = FALSE]
                source_eval = source[source_eval_index, , drop = FALSE]
                target_eval = target[target_eval_index, , drop = FALSE]
                weight_fit = if (need_weight) 
                  fit_transfer_weight_model(source_train, target_train)
                else NULL
                nuisance_fit = fit_nuisance_models(source_train, target_train, need_target = need_target)
                evaluation = predict_nuisance_scores(nuisance_fit, weight_fit, source_eval, target_eval, 
                  need_target = need_target)
                source_score[source_eval_index, names(evaluation$source)] = evaluation$source
                if (need_target) 
                  target_score[target_eval_index, ] = evaluation$target
                fold_fits[[k]] = list(nuisance = nuisance_fit, weight = weight_fit, source_train_index = source_train_index, 
                  target_train_index = target_train_index, source_eval_index = source_eval_index, target_eval_index = target_eval_index)
            }
            source_score$fold = folds$source
            if (need_target) 
                target_score$fold = folds$target
            model_names = c("mu11", "mu10")
            if (need_target) 
                model_names = c(model_names, "mu01", "mu00")
            averaged_models = setNames(lapply(model_names, function(name) average_fold_model(fold_fits, 
                name)), model_names)
            list(source = source_score, target = target_score, models = averaged_models, fold_fits = fold_fits, 
                folds = folds, samples = samples, n_folds = n_folds)
        }
    if (method %in% c("IP", "CC", "Fusion")) 
        fit_integrative_r_fold <- function(source, target, nuisance, lambda_tau = 0, lambda_c = 0, inner_predictions = NULL) {
            if (is.null(inner_predictions)) 
                stop("R-loss requires inner out-of-fold nuisance predictions.")
            e1 = rep(0.5, nrow(source))
            e0 = inner_predictions$e0
            source_residual_A = source$A - e1
            target_residual_A = target$A - e0
            response = c(source$Y - inner_predictions$m1, target$Y - inner_predictions$m0)
            tau_basis = rbind(make_basis(source), make_basis(target))
            c_basis_source = matrix(0, nrow(source), ncol(make_pointwise_confounding_basis(model$example)))
            c_basis_target = make_pointwise_confounding_basis(target)
            c_basis = rbind(c_basis_source, c_basis_target)
            residual_A = c(source_residual_A, target_residual_A)
            design = cbind(tau_basis * residual_A, c_basis * residual_A)
            p_tau = ncol(tau_basis)
            p_c = ncol(c_basis)
            penalty = diag(c(rep(lambda_tau, p_tau), rep(lambda_c, p_c)))
            lhs = crossprod(design)/nrow(design) + penalty
            rhs = crossprod(design, response)/nrow(design)
            coefficient = tryCatch(as.vector(solve(lhs, rhs)), error = function(e) as.vector(qr.solve(lhs, 
                rhs)))
            coefficient[!is.finite(coefficient)] = 0
            list(tau_coef = coefficient[seq_len(p_tau)], c_coef = coefficient[p_tau + seq_len(p_c)])
        }
    fit_linear_basis <- function(y, basis) {
        fit = lm.fit(x = as.matrix(basis), y = y)
        coefficient = as.numeric(fit$coefficients)
        coefficient[!is.finite(coefficient)] = 0
        coefficient
    }
    fit_nuisance_models <- function(source, target, need_target = TRUE) {
        models = list(mu11 = fit_arm_model(source, 1), mu10 = fit_arm_model(source, 0), m1 = fit_marginal_outcome_model(source))
        if (need_target) {
            models$mu01 = fit_arm_model(target, 1)
            models$mu00 = fit_arm_model(target, 0)
            models$m0 = fit_marginal_outcome_model(target)
            models$e0 = fit_target_propensity(target)
        }
        models
    }
    fit_pooled_nuisance <- function(train) {
        basis = model$arm_basis(train)
        use1 = train$A == 1
        use0 = train$A == 0
        mu1_coef = fit_linear_basis(train$Y[use1], basis[use1, , drop = FALSE])
        mu0_coef = fit_linear_basis(train$Y[use0], basis[use0, , drop = FALSE])
        propensity_fit = glm.fit(x = as.matrix(basis), y = train$A, family = binomial())
        propensity_coef = as.numeric(propensity_fit$coefficients)
        propensity_coef[!is.finite(propensity_coef)] = 0
        list(mu1_coef = mu1_coef, mu0_coef = mu0_coef, propensity_coef = propensity_coef)
    }
    fit_target_propensity <- function(target) {
        glm(reformulate(model$covariates, response = "A"), family = binomial(), data = target)
    }
    fit_transfer_weight_model <- function(source, target) {
        pooled <- rbind(transform(source, S = 1), transform(target, S = 0))
        model <- glm(model$weight_formula, family = binomial(), data = pooled)
        if (!isTRUE(model$converged) || any(!is.finite(coef(model)))) 
            stop("Weight model failed")
        ratio <- nrow(source)/nrow(target)
        norm <- mean(ratio * exp(-predict(model, newdata = source, type = "link")))
        stopifnot(is.finite(norm), norm > 0)
        list(model = model, sample_ratio = ratio, normalizer = norm)
    }
    make_crossfit_folds <- function(source, target, n_folds = 5, seed = 1) {
        list(source = make_stratified_folds(source$A, n_folds, seed + 1), target = make_stratified_folds(target$A, 
            n_folds, seed + 2), n_folds = n_folds)
    }
    make_pooled_sample <- function(data) {
        if (method == "Target") 
            data$target
        else rbind(data$source, data$target)
    }
    make_samples <- function(data) list(source = data$source, target = data$target)
    make_stratified_folds <- function(A, n_folds = 5, seed = 1) {
        if (length(n_folds) != 1L || !is.finite(n_folds) || n_folds < 2 || n_folds != as.integer(n_folds)) 
            stop("n_folds must be an integer >= 2.")
        if (length(A)%%n_folds != 0L) 
            stop("Equal folds require the sample size to be divisible by n_folds.")
        set.seed(seed)
        sample(rep(seq_len(n_folds), each = length(A)/n_folds))
    }
    if (method %in% c("WY")) 
        predict_nuisance_scores <- function(models, weight_fit, source, target = NULL, need_target = TRUE) {
            mu11 = predict_basis_model(models$mu11, source)
            mu10 = predict_basis_model(models$mu10, source)
            tau1_hat = mu11 - mu10
            Gamma1_hat = tau1_hat + source$A * (source$Y - mu11)/0.5 - (1 - source$A) * (source$Y - mu10)/0.5
            source_result = data.frame(Gamma1_hat = Gamma1_hat, tau1_hat = tau1_hat, m1_hat = predict_basis_model(models$m1, 
                source), w = predict_transfer_weight(weight_fit, source))
            if (!need_target) {
                return(list(source = source_result, target = NULL))
            }
            mu01 = predict_basis_model(models$mu01, target)
            mu00 = predict_basis_model(models$mu00, target)
            omega0_hat = mu01 - mu00
            e0_hat = predict_target_propensity(models$e0, target)
            Gamma0_hat = omega0_hat + target$A * (target$Y - mu01)/e0_hat - (1 - target$A) * (target$Y - 
                mu00)/(1 - e0_hat)
            source_result$omega0_source_hat = predict_basis_model(models$mu01, source) - predict_basis_model(models$mu00, 
                source)
            target_result = data.frame(Gamma0_hat = Gamma0_hat, omega0_hat = omega0_hat, m0_hat = predict_basis_model(models$m0, 
                target), e0_hat = e0_hat, w = predict_transfer_weight(weight_fit, target))
            list(source = source_result, target = target_result)
        }
    if (method %in% c("CC", "Fusion")) 
        predict_nuisance_scores <- function(models, weight_fit, source, target = NULL, need_target = TRUE) {
            mu11 = predict_basis_model(models$mu11, source)
            mu10 = predict_basis_model(models$mu10, source)
            tau1_hat = mu11 - mu10
            Gamma1_hat = tau1_hat + source$A * (source$Y - mu11)/0.5 - (1 - source$A) * (source$Y - mu10)/0.5
            source_result = data.frame(Gamma1_hat = Gamma1_hat, tau1_hat = tau1_hat, m1_hat = predict_basis_model(models$m1, 
                source), w = if (is.null(weight_fit)) 
                rep(1, nrow(source))
            else predict_transfer_weight(weight_fit, source))
            if (!need_target) {
                return(list(source = source_result, target = NULL))
            }
            mu01 = predict_basis_model(models$mu01, target)
            mu00 = predict_basis_model(models$mu00, target)
            omega0_hat = mu01 - mu00
            e0_hat = predict_target_propensity(models$e0, target)
            Gamma0_hat = omega0_hat + target$A * (target$Y - mu01)/e0_hat - (1 - target$A) * (target$Y - 
                mu00)/(1 - e0_hat)
            source_result$omega0_source_hat = predict_basis_model(models$mu01, source) - predict_basis_model(models$mu00, 
                source)
            target_result = data.frame(Gamma0_hat = Gamma0_hat, omega0_hat = omega0_hat, m0_hat = predict_basis_model(models$m0, 
                target), e0_hat = e0_hat, w = if (is.null(weight_fit)) 
                rep(1, nrow(target))
            else predict_transfer_weight(weight_fit, target))
            list(source = source_result, target = target_result)
        }
    predict_pooled_score <- function(fit, test) {
        basis = model$arm_basis(test)
        mu1 = as.vector(basis %*% fit$mu1_coef)
        mu0 = as.vector(basis %*% fit$mu0_coef)
        propensity = clip(plogis(as.vector(basis %*% fit$propensity_coef)))
        mu1 - mu0 + test$A * (test$Y - mu1)/propensity - (1 - test$A) * (test$Y - mu0)/(1 - propensity)
    }
    predict_target_propensity <- function(fit, data) {
        clip(as.numeric(predict(fit, newdata = data, type = "response")))
    }
    predict_transfer_weight <- function(fit, data) {
        weight = fit$sample_ratio * exp(-as.numeric(predict(fit$model, newdata = data, type = "link")))/fit$normalizer
        if (any(!is.finite(weight)) || any(weight <= 0)) 
            stop("Nonfinite/underflowed density ratios: inspect overlap/model fit.")
        weight
    }
    wu2020_sample <- function(data, domain = c("rwd", "rct")) {
        if (match.arg(domain) == "rwd") 
            data$target
        else data$source
    }
    wu2020_target_folds <- function(n, n_folds = 5L, seed = 1L) {
        if (n_folds < 2L || n%%n_folds != 0L) 
            stop("Wu2020 common-D projection requires equal target folds.")
        set.seed(seed)
        sample(rep(seq_len(n_folds), each = n/n_folds))
    }
    if (method %in% c("CC", "Fusion")) 
        fit_cross_fitted_pointwise_confounding <- function(crossfit, lambda_tau = 0, lambda_c = 0) {
            source = crossfit$samples$source
            target = crossfit$samples$target
            p_tau = ncol(make_basis(model$example))
            p_c = ncol(make_pointwise_confounding_basis(model$example))
            tau_coef_by_fold = matrix(NA_real_, crossfit$n_folds, p_tau)
            c_coef_by_fold = matrix(NA_real_, crossfit$n_folds, p_c)
            c_source = rep(NA_real_, nrow(source))
            c_target = rep(NA_real_, nrow(target))
            for (k in seq_len(crossfit$n_folds)) {
                fold_fit = crossfit$fold_fits[[k]]
                source_train = source[fold_fit$source_train_index, , drop = FALSE]
                target_train = target[fold_fit$target_train_index, , drop = FALSE]
                inner_predictions = fit_inner_r_predictions(source_train, target_train, n_inner_folds = 5, 
                  seed = 20000 + k)
                pointwise_fit = fit_integrative_r_fold(source_train, target_train, fold_fit$nuisance, 
                  lambda_tau = lambda_tau, lambda_c = lambda_c, inner_predictions = inner_predictions)
                tau_coef_by_fold[k, ] = pointwise_fit$tau_coef
                c_coef_by_fold[k, ] = pointwise_fit$c_coef
                c_source[fold_fit$source_eval_index] = as.vector(make_pointwise_confounding_basis(source[fold_fit$source_eval_index, 
                  , drop = FALSE]) %*% pointwise_fit$c_coef)
                c_target[fold_fit$target_eval_index] = as.vector(make_pointwise_confounding_basis(target[fold_fit$target_eval_index, 
                  , drop = FALSE]) %*% pointwise_fit$c_coef)
            }
            list(c_coef = colMeans(c_coef_by_fold), tau_coef = colMeans(tau_coef_by_fold), c_coef_by_fold = c_coef_by_fold, 
                tau_coef_by_fold = tau_coef_by_fold, c_source = c_source, c_target = c_target)
        }
    fit_inner_r_predictions <- function(source, target, n_inner_folds = 5, seed = 1) {
        folds = make_crossfit_folds(source, target, n_inner_folds, seed)
        m1 = rep(NA_real_, nrow(source))
        m0 = e0 = rep(NA_real_, nrow(target))
        for (j in seq_len(n_inner_folds)) {
            tr1 = which(folds$source != j)
            ev1 = which(folds$source == j)
            tr0 = which(folds$target != j)
            ev0 = which(folds$target == j)
            fit1 = fit_marginal_outcome_model(source[tr1, , drop = FALSE])
            fit0 = fit_marginal_outcome_model(target[tr0, , drop = FALSE])
            prop0 = fit_target_propensity(target[tr0, , drop = FALSE])
            m1[ev1] = predict_basis_model(fit1, source[ev1, , drop = FALSE])
            m0[ev0] = predict_basis_model(fit0, target[ev0, , drop = FALSE])
            e0[ev0] = predict_target_propensity(prop0, target[ev0, , drop = FALSE])
        }
        if (any(!is.finite(c(m1, m0, e0)))) 
            stop("Nonfinite inner out-of-fold nuisance predictions.")
        list(m1 = m1, m0 = m0, e0 = e0, folds = folds)
    }
    wu2020_q_fit <- function(z, A, Y) {
        z <- as.matrix(z)
        coefs <- matrix(0, nrow = ncol(z), ncol = 2L, dimnames = list(colnames(z), c("A0", "A1")))
        ranks <- counts <- integer(2)
        for (a in 0:1) {
            use <- A == a
            if (!any(use)) 
                stop("No RCT observations in an arm; cannot fit Q regression.")
            fit <- lm.fit(z[use, , drop = FALSE], Y[use])
            beta <- fit$coefficients
            beta[is.na(beta)] <- 0
            if (any(!is.finite(beta))) 
                stop("Nonfinite Q regression coefficients.")
            coefs[, a + 1L] <- beta
            ranks[a + 1L] <- fit$rank
            counts[a + 1L] <- sum(use)
        }
        list(coef = coefs, contrast = coefs[, 2L] - coefs[, 1L], rank = ranks, counts = counts)
    }
    wu2020_fit_features <- function(target, ntree = 500L, seed = 1L) {
        d0 <- wu2020_pretrain(target, target, ntree, seed)
        agree <- target$A == d0
        if (sum(agree) < 2L || sum(!agree) < 2L) 
            stop("Need at least two agreeing and two disagreeing target observations.")
        x <- wu2020_x(target)
        list(h1 = wu2020_rf(x[agree, ], target$A[agree], ntree, seed + 101L, classification = TRUE), 
            benefit_agree = wu2020_rf(x[agree, ], target$Y[agree], ntree, seed + 102L), benefit_disagree = wu2020_rf(x[!agree, 
                ], target$Y[!agree], ntree, seed + 103L))
    }
    wu2020_predict_features <- function(fit, data) {
        x <- wu2020_x(data)
        cbind(H1 = wu2020_rf_predict(fit$h1, x), H2 = wu2020_rf_predict(fit$benefit_agree, x) - wu2020_rf_predict(fit$benefit_disagree, 
            x))
    }
    wu2020_pretrain <- function(train, test, ntree, seed) {
        if (length(unique(train$A)) != 2L) 
            stop("Step-0 training sample must contain both treatments.")
        q <- wu2020_rf(cbind(wu2020_x(train), A = train$A), train$Y, ntree, seed)
        q1 <- wu2020_rf_predict(q, cbind(wu2020_x(test), A = 1))
        q0 <- wu2020_rf_predict(q, cbind(wu2020_x(test), A = 0))
        as.integer(q1 > q0)
    }
    wu2020_rf <- function(x, y, ntree, seed, classification = FALSE) {
        if (nrow(x) < 2L) 
            stop("Too few observations for a super-feature model.")
        if (length(unique(y)) == 1L) 
            return(list(constant = as.numeric(as.character(y[1]))))
        wu2020_with_seed(seed, function() {
            list(model = randomForest::randomForest(x = x, y = if (classification) factor(y, levels = 0:1) else y, 
                ntree = ntree, mtry = if (classification) max(1L, floor(sqrt(ncol(x)))) else max(1L, 
                  floor(ncol(x)/3)), nodesize = if (classification) 1L else 5L), classification = classification)
        })
    }
    wu2020_rf_predict <- function(fit, x) {
        if (!is.null(fit$constant)) 
            return(rep(fit$constant, nrow(x)))
        if (isTRUE(fit$classification)) 
            return(as.numeric(predict(fit$model, x, type = "prob")[, "1"]))
        as.numeric(predict(fit$model, x))
    }
    wu2020_with_seed <- function(seed, fun) {
        had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had) 
            old <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (had) assign(".Random.seed", old, envir = .GlobalEnv) else if (exists(".Random.seed", 
                envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        })
        set.seed(seed)
        fun()
    }
    wu2020_x <- function(z) z[, model$covariates, drop = FALSE]
    empirical_value <- function(Z, r, eta) {
        sum(r * as.numeric(as.vector(Z %*% eta) > 0))
    }
    fit_arm_model <- function(data, arm) {
        use <- data$A == arm
        fit <- fit_basis_regression(data$Y[use], model$arm_basis(data[use, , drop = FALSE]))
        fit$basis <- model$arm_basis
        fit
    }
    fit_marginal_outcome_model <- function(data) {
        fit <- fit_basis_regression(data$Y, model$marginal_basis(data))
        fit$basis <- model$marginal_basis
        fit
    }
    predict_basis_model <- function(fit, data) {
        basis <- if (is.null(fit$basis)) 
            model$p
        else fit$basis
        as.numeric(basis(data) %*% fit$coefficient)
    }
    average_fold_model <- function(fold_fits, model_name) {
        coefficient <- vapply(fold_fits, function(fit) fit$nuisance[[model_name]]$coefficient, numeric(ncol(model$arm_basis(model$example))))
        list(coefficient = rowMeans(coefficient), basis = model$arm_basis)
    }
    environment()
}

drgo_correction <- function(scores, model) {
    rct <- scores$samples$source
    tgt <- scores$samples$target
    n1 <- nrow(rct)
    n0 <- nrow(tgt)
    rho0 <- n0/(n0 + n1)
    rho1 <- n1/(n0 + n1)
    w_rct <- scores$source$w
    w_tgt <- scores$target$w
    g_rct <- 1/(rho0 * w_rct + rho1)
    g_tgt <- 1/(rho0 * w_tgt + rho1)
    phi1_rct <- scores$source$Gamma1_hat
    z_rct <- model$z(rct)
    z_tgt <- model$z(tgt)
    H <- crossprod(z_rct * (w_rct * g_rct), z_rct)/n1
    u_rct <- colSums(z_rct * (w_rct * g_rct * (scores$source$omega0_source_hat - phi1_rct)))/n1
    u_tgt <- colSums(z_tgt * (g_tgt * (scores$target$Gamma0_hat - scores$target$omega0_hat)))/n0
    gamma <- drop(solve(H, u_rct + u_tgt))
    list(gamma = gamma, c_target = drop(z_tgt %*% gamma))
}
