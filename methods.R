# Public interface: fit_itr(sample, method). Returned theta defines
# d(x) = 1{theta[1] + theta[2]*x1 > 0}. All nine learners use the same
# 16-start DC/BFGS policy search as the published run.
method_names <- c("WY","Target","Pooled","Wu2020","IP","CC","Fusion","CC_SCAD","Fusion_SCAD")

fit_itr <- function(sample, method = "CC") {
  stopifnot(method %in% method_names)
  model <- paper_setting(sample$setting)
  kernel <- fitting_functions(method, model)
  # Only observed training records enter the fitting routines.
  tr <- sample$trial; tg <- sample$target
  training <- list(n1=nrow(tr), n0=nrow(tg), x1_rct=tr$x1, x2_rct=tr$x2,
                   A_rct=tr$A, Y_rct=tr$Y, x1_rwd=tg$x1, x2_rwd=tg$x2,
                   A_rwd=tg$A, Y_rwd=tg$Y)
  control <- list(data=training, n_starts=16L, max_dc=200L, seed=1000+sample$seed)
  fit <- switch(method,
    IP         = do.call(kernel$fit_integrative_plugin, control),
    WY         = do.call(kernel$fit_wuyang, control),
    Target     = , Pooled = do.call(kernel$fit_naive_pooled_itr, control),
    Wu2020     = do.call(kernel$fit_wu2020, c(control,list(ntree=500L))),
    CC         = , CC_SCAD = do.call(kernel$fit_fixed_alpha_path, control),
    Fusion     = , Fusion_SCAD = do.call(kernel$fit_proposed_q, control))
  theta <- if (method %in% c("CC","CC_SCAD")) fit$eta_by_alpha[[1]] else fit$eta_hat
  audits <- kernel$diagnostics
  stopifnot(length(theta)==2L, all(is.finite(theta)), audits$starts>=16L,
            audits$starts%%16L==0L, audits$optimizer_calls>=audits$starts)
  list(theta=theta, q_hat=if (method=="Fusion") fit$q_fit$q_target else NULL,
       dc_starts=audits$starts, bfgs_calls=audits$optimizer_calls,
       rank_deficient_fits=audits$rank_deficient_fits)
}

# One lexical scope per fit keeps estimator-specific functions separate.
# Shared regression/optimization kernels are defined once. Where the archived
# implementations differ, explicit method branches retain their exact behavior.
fitting_functions <- function(method, model) {
  make_basis <- function(x1,x2) polynomial_basis(x1,x2,model$degree)
  make_policy_basis <- function(x1,x2) cbind(1,x1)
  make_confounding_basis <- make_pointwise_confounding_basis <- model$z
  make_variance_basis <- function(x1,x2) cbind(polynomial_basis(x1,x2,2L),x1^4)
  diagnostics <- new.env(parent=emptyenv())
  diagnostics$starts <- diagnostics$optimizer_calls <- diagnostics$optimizer_nonconvergence <- 0L
  diagnostics$rank_deficient_fits <- diagnostics$least_squares_fits <- 0L
  # --- Policy scores and fitting entry points ---
  if (method %in% c("CC")) fit_fixed_alpha_path <- function (data, alpha_grid = seq(0, 0, by = 0), lambda_tau = 0, lambda_c = 0, n_folds = 5, n_starts = 8, 
      max_dc = 100, seed = 1) 
  {
      scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = TRUE, need_weight = any(alpha_grid != 
          0))
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

  if (method %in% c("CC_SCAD")) fit_fixed_alpha_path <- function (data, alpha_grid = seq(0, 0, by = 0), lambda_tau = NULL, lambda_c = NULL, n_folds = 5, 
      n_starts = 8, max_dc = 100, seed = 1) 
  {
      scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = TRUE, need_weight = any(alpha_grid != 
          0))
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

  fit_integrative_plugin <- function (data, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) 
  {
      samples = make_samples(data)
      source = samples$source
      target = samples$target
      folds = make_crossfit_folds(source, target, n_folds, seed)
      tau_target = rep(NA_real_, nrow(target))
      tau_coef_by_fold = matrix(NA_real_, n_folds, ncol(make_basis(0, 0)))
      c_coef_by_fold = matrix(NA_real_, n_folds, ncol(make_confounding_basis(0, 0)))
      for (k in seq_len(n_folds)) {
          tr1 = source[folds$source != k, , drop = FALSE]
          tr0 = target[folds$target != k, , drop = FALSE]
          ev0 = which(folds$target == k)
          inner = fit_inner_r_predictions(tr1, tr0, n_inner_folds = 5, seed = 20000 + k)
          joint = fit_integrative_r_fold(tr1, tr0, nuisance = NULL, lambda_tau = 0, lambda_c = 0, 
              inner_predictions = inner)
          tau_coef_by_fold[k, ] = joint$tau_coef
          c_coef_by_fold[k, ] = joint$c_coef
          tau_target[ev0] = as.vector(make_basis(target$x1[ev0], target$x2[ev0]) %*% joint$tau_coef)
      }
      if (any(!is.finite(tau_target))) 
          stop("Nonfinite out-of-fold treatment-effect prediction.")
      Z = make_policy_basis(target$x1, target$x2)
      r = tau_target/nrow(target)
      eta_hat = dc_policy(Z, r, n_starts = n_starts, max_dc = max_dc, seed = seed + 1)
      list(eta_hat = eta_hat, tau_target = tau_target, tau_coef_by_fold = tau_coef_by_fold, c_coef_by_fold = c_coef_by_fold, 
          folds = folds, empirical_value = empirical_value(Z, r, eta_hat))
  }

  fit_naive_pooled_itr <- function (data, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) 
  {
      score = fit_cross_fitted_pooled_score(data, n_folds, seed)
      Z = make_policy_basis(score$pooled$x1, score$pooled$x2)
      r = score$Gamma_pool/nrow(score$pooled)
      eta_hat = dc_policy(Z, r, n_starts = n_starts, max_dc = max_dc, seed = seed + 1)
      list(eta_hat = eta_hat, association_coef = score$association_coef, pooled_empirical_value = empirical_value(Z, 
          r, eta_hat))
  }

  if (method %in% c("Fusion")) fit_proposed_q <- function (data, lambda_tau = 0, lambda_c = 0, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) 
  {
      scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = TRUE, need_weight = TRUE)
      pointwise_c = fit_cross_fitted_pointwise_confounding(scores, lambda_tau = lambda_tau, lambda_c = lambda_c)
      q_fit = fit_cross_fitted_q(scores, n_inner_folds = 5)
      objective = make_q_objective(scores, pointwise_c, q_fit)
      if (any(!is.finite(objective$r))) 
          stop("Nonfinite residual-fusion objective.")
      eta_hat = dc_policy(objective$Z, objective$r, n_starts = n_starts, max_dc = max_dc, seed = seed + 
          1)
      list(eta_hat = eta_hat, scores = scores, pointwise_c = pointwise_c, q_fit = q_fit, q_empirical_value = empirical_value(objective$Z, 
          objective$r, eta_hat))
  }

  if (method %in% c("Fusion_SCAD")) fit_proposed_q <- function (data, lambda_tau = NULL, lambda_c = NULL, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) 
  {
      scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = TRUE, need_weight = TRUE)
      pointwise_c = fit_cross_fitted_pointwise_confounding(scores, lambda_tau = lambda_tau, lambda_c = lambda_c)
      q_fit = fit_cross_fitted_q(scores, n_inner_folds = 5)
      objective = make_q_objective(scores, pointwise_c, q_fit)
      if (any(!is.finite(objective$r))) 
          stop("Nonfinite residual-fusion objective.")
      eta_hat = dc_policy(objective$Z, objective$r, n_starts = n_starts, max_dc = max_dc, seed = seed + 
          1)
      list(eta_hat = eta_hat, scores = scores, pointwise_c = pointwise_c, q_fit = q_fit, q_empirical_value = empirical_value(objective$Z, 
          objective$r, eta_hat))
  }

  fit_wu2020 <- function (data, ntree = 500L, n_folds = 5L, n_starts = 8L, max_dc = 100L, seed = 1L) 
  {
      target <- wu2020_sample(data, "rwd")
      rct <- wu2020_sample(data, "rct")
      folds <- wu2020_target_folds(nrow(target), n_folds, seed + 1L)
      tau_score_oof <- rep(NA_real_, nrow(target))
      q_rank_by_fold <- matrix(NA_integer_, n_folds, 2L, dimnames = list(NULL, c("A0", "A1")))
      for (ell in seq_len(n_folds)) {
          train <- target[folds != ell, , drop = FALSE]
          eval <- target[folds == ell, , drop = FALSE]
          feature_fit <- wu2020_fit_features(train, ntree, seed + 1000L * ell)
          h_rct <- wu2020_predict_features(feature_fit, rct$x1, rct$x2)
          z_rct <- cbind(make_basis(rct$x1, rct$x2), h_rct)
          q_fit <- wu2020_q_fit(z_rct, rct$A, rct$Y)
          h_eval <- wu2020_predict_features(feature_fit, eval$x1, eval$x2)
          z_eval <- cbind(make_basis(eval$x1, eval$x2), h_eval)
          tau_score_oof[folds == ell] <- as.vector(z_eval %*% q_fit$contrast)
          q_rank_by_fold[ell, ] <- q_fit$rank
      }
      if (any(!is.finite(tau_score_oof))) 
          stop("Nonfinite cross-fitted Wu2020 score.")
      policy_basis <- make_policy_basis(target$x1, target$x2)
      eta_hat <- dc_policy(Z = policy_basis, r = tau_score_oof/nrow(target), n_starts = n_starts, 
          max_dc = max_dc, seed = seed + 90000L)
      list(eta_hat = eta_hat, basis_fun = make_policy_basis, tau_score_oof = tau_score_oof, target_folds = folds, 
          q_rank_by_fold = q_rank_by_fold, projected_empirical_value = empirical_value(policy_basis, 
              tau_score_oof/nrow(target), eta_hat))
  }

  fit_wuyang <- function (data, n_folds = 5, n_starts = 8, max_dc = 100, seed = 1) 
  {
      scores = fit_cross_fitted_scores(data, n_folds = n_folds, seed = seed, need_target = FALSE)
      source = scores$samples$source
      Z1 = make_policy_basis(source$x1, source$x2)
      r1 = scores$source$w * scores$source$Gamma1_hat/nrow(source)
      eta_hat = dc_policy(Z1, r1, n_starts = n_starts, max_dc = max_dc, seed = seed + 1)
      list(eta_hat = eta_hat, scores = scores, source_empirical_value = empirical_value(Z1, r1, eta_hat))
  }

  make_alpha_objective <- function (source, target, scores, pointwise_c, alpha) 
  {
      Z0 = make_policy_basis(target$x1, target$x2)
      Z1 = make_policy_basis(source$x1, source$x2)
      r0 = (1 - alpha) * (scores$target$Gamma0_hat - pointwise_c$c_target)/nrow(target)
      r1 = alpha * scores$source$w * scores$source$Gamma1_hat/nrow(source)
      list(Z = rbind(Z0, Z1), r = c(r0, r1))
  }

  make_q_objective <- function (scores, pointwise_c, q_fit) 
  {
      source = scores$samples$source
      target = scores$samples$target
      corrected_target = scores$target$Gamma0_hat - pointwise_c$c_target
      r0 = q_fit$q_target * corrected_target/nrow(target)
      r1 = scores$source$w * (1 - q_fit$q_source) * scores$source$Gamma1_hat/nrow(source)
      list(Z = rbind(make_policy_basis(target$x1, target$x2), make_policy_basis(source$x1, source$x2)), 
          r = c(r0, r1))
  }

  # --- Cross-fitting, outcome regressions, and density ratios ---
  average_fold_model <- function (fold_fits, model_name) 
  {
      coefficient = vapply(fold_fits, function(fit) fit$nuisance[[model_name]]$coefficient, numeric(ncol(make_basis(0, 
          0))))
      list(coefficient = rowMeans(coefficient))
  }

  clip <- function (x, lower = 0.02, upper = 0.98) 
  {
      pmin(pmax(x, lower), upper)
  }

  fit_arm_model <- function (data, arm) 
  {
      use = data$A == arm
      fit_basis_regression(y = data$Y[use], basis = make_basis(data$x1[use], data$x2[use]))
  }

  fit_basis_regression <- function (y, basis) 
  {
      fit <- lm.fit(x = as.matrix(basis), y = y)
      diagnostics$least_squares_fits <- diagnostics$least_squares_fits + 1L
      if (fit$rank < ncol(basis)) 
          diagnostics$rank_deficient_fits <- diagnostics$rank_deficient_fits + 1L
      coefficient <- as.numeric(fit$coefficients)
      coefficient[!is.finite(coefficient)] <- 0
      list(coefficient = coefficient)
  }

  fit_cross_fitted_pooled_score <- function (data, n_folds = 5, seed = 1) 
  {
      pooled = make_pooled_sample(data)
      fold = make_stratified_folds(pooled$A, n_folds, seed)
      gamma_pool = rep(NA_real_, nrow(pooled))
      association_coef_by_fold = matrix(NA_real_, nrow = n_folds, ncol = ncol(make_basis(0, 0)))
      for (k in seq_len(n_folds)) {
          train = pooled[fold != k, , drop = FALSE]
          test = pooled[fold == k, , drop = FALSE]
          nuisance = fit_pooled_nuisance(train)
          gamma_pool[fold == k] = predict_pooled_score(nuisance, test)
          association_coef_by_fold[k, ] = nuisance$mu1_coef - nuisance$mu0_coef
      }
      list(pooled = pooled, Gamma_pool = gamma_pool, association_coef = colMeans(association_coef_by_fold))
  }

  if (method %in% c("WY")) fit_cross_fitted_scores <- function (data, n_folds = 5, seed = 1, need_target = TRUE) 
  {
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

  if (method %in% c("CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) fit_cross_fitted_scores <- function (data, n_folds = 5, seed = 1, need_target = TRUE, need_weight = TRUE) 
  {
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

  if (method %in% c("IP", "CC", "Fusion")) fit_integrative_r_fold <- function (source, target, nuisance, lambda_tau = 0, lambda_c = 0, inner_predictions = NULL) 
  {
      if (is.null(inner_predictions)) 
          stop("R-loss requires inner out-of-fold nuisance predictions.")
      e1 = rep(0.5, nrow(source))
      e0 = inner_predictions$e0
      source_residual_A = source$A - e1
      target_residual_A = target$A - e0
      response = c(source$Y - inner_predictions$m1, target$Y - inner_predictions$m0)
      tau_basis = rbind(make_basis(source$x1, source$x2), make_basis(target$x1, target$x2))
      c_basis_source = matrix(0, nrow(source), ncol(make_pointwise_confounding_basis(0, 0)))
      c_basis_target = make_pointwise_confounding_basis(target$x1, target$x2)
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

  if (method %in% c("CC_SCAD", "Fusion_SCAD")) fit_integrative_r_fold <- function (source, target, nuisance, lambda_tau = NULL, lambda_c = NULL, inner_predictions = NULL, 
      sc_grid = c(0, 0.5, 1, 1.5), n_penalty_folds = 5, nlambda = 50, seed = 1) 
  {
      if (is.null(inner_predictions)) 
          stop("R-loss requires inner out-of-fold nuisance predictions.")
      e1 = rep(0.5, nrow(source))
      e0 = inner_predictions$e0
      source_residual_A = source$A - e1
      target_residual_A = target$A - e0
      response = c(source$Y - inner_predictions$m1, target$Y - inner_predictions$m0)
      tau_basis = rbind(make_basis(source$x1, source$x2), make_basis(target$x1, target$x2))
      c_basis_source = matrix(0, nrow(source), ncol(make_pointwise_confounding_basis(0, 0)))
      c_basis_target = make_pointwise_confounding_basis(target$x1, target$x2)
      c_basis = rbind(c_basis_source, c_basis_target)
      residual_A = c(source_residual_A, target_residual_A)
      design = cbind(tau_basis * residual_A, c_basis * residual_A)
      p_tau = ncol(tau_basis)
      p_c = ncol(c_basis)
      if (!is.null(lambda_tau) || !is.null(lambda_c)) {
          if (is.null(lambda_tau) || is.null(lambda_c)) 
              stop("Specify both lambda_tau and lambda_c, or leave both NULL for SCAD CV.")
          penalty = diag(c(rep(lambda_tau, p_tau), rep(lambda_c, p_c)))
          lhs = crossprod(design)/nrow(design) + penalty
          rhs = crossprod(design, response)/nrow(design)
          coefficient = tryCatch(as.vector(solve(lhs, rhs)), error = function(e) as.vector(qr.solve(lhs, 
              rhs)))
          tuning = list(method = "ridge", sc = NA_real_, lambda_c = lambda_c, lambda_tau = lambda_tau, 
              cv_error = NA_real_)
      }
      else {
          if (!requireNamespace("ncvreg", quietly = TRUE))
              stop("Package 'ncvreg' is required for the two-penalty SCAD fit.")
          if (any(!is.finite(design)) || any(!is.finite(response))) 
              stop("Nonfinite R-loss design or response before SCAD fitting.")
          n_penalty_folds = min(n_penalty_folds, nrow(source), nrow(target))
          if (n_penalty_folds < 2) 
              stop("At least two observations per study are needed for SCAD CV.")
          set.seed(seed)
          source_fold = sample(rep(seq_len(n_penalty_folds), length.out = nrow(source)))
          target_fold = sample(rep(seq_len(n_penalty_folds), length.out = nrow(target)))
          penalty_fold = c(source_fold, target_fold)
          cv_fits = lapply(sc_grid, function(sc) {
              ncvreg::cv.ncvreg(X = design, y = response, family = "gaussian", penalty = "SCAD", gamma = 3.7, 
                  penalty.factor = c(rep(sc, p_tau), rep(1, p_c)), fold = penalty_fold, nfolds = n_penalty_folds, 
                  nlambda = nlambda, returnY = FALSE)
          })
          cv_error = vapply(cv_fits, function(fit) {
              if (any(is.finite(fit$cve))) 
                  min(fit$cve, na.rm = TRUE)
              else Inf
          }, numeric(1))
          if (!any(is.finite(cv_error))) 
              stop("All two-penalty SCAD fits failed.")
          best = which.min(cv_error)
          selected_fit = cv_fits[[best]]
          coefficient = as.numeric(stats::coef(selected_fit))
          if (length(coefficient) == ncol(design) + 1) 
              coefficient = coefficient[-1]
          if (length(coefficient) != ncol(design)) 
              stop("Unexpected coefficient dimension returned by ncvreg.")
          tuning = list(method = "SCAD-CV", sc = sc_grid[best], lambda_c = selected_fit$lambda.min, 
              lambda_tau = sc_grid[best] * selected_fit$lambda.min, cv_error = cv_error[best], cv_error_by_sc = cv_error, 
              sc_grid = sc_grid)
      }
      coefficient[!is.finite(coefficient)] = 0
      list(tau_coef = coefficient[seq_len(p_tau)], c_coef = coefficient[p_tau + seq_len(p_c)], tuning = tuning)
  }

  fit_linear_basis <- function (y, basis) 
  {
      fit = lm.fit(x = as.matrix(basis), y = y)
      coefficient = as.numeric(fit$coefficients)
      coefficient[!is.finite(coefficient)] = 0
      coefficient
  }

  fit_marginal_outcome_model <- function (data) 
  {
      fit_basis_regression(y = data$Y, basis = make_basis(data$x1, data$x2))
  }

  fit_nuisance_models <- function (source, target, need_target = TRUE) 
  {
      models = list(mu11 = fit_arm_model(source, 1), mu10 = fit_arm_model(source, 0), m1 = fit_marginal_outcome_model(source))
      if (need_target) {
          models$mu01 = fit_arm_model(target, 1)
          models$mu00 = fit_arm_model(target, 0)
          models$m0 = fit_marginal_outcome_model(target)
          models$e0 = fit_target_propensity(target)
      }
      models
  }

  fit_pooled_nuisance <- function (train) 
  {
      basis = make_basis(train$x1, train$x2)
      use1 = train$A == 1
      use0 = train$A == 0
      mu1_coef = fit_linear_basis(train$Y[use1], basis[use1, , drop = FALSE])
      mu0_coef = fit_linear_basis(train$Y[use0], basis[use0, , drop = FALSE])
      propensity_fit = glm.fit(x = as.matrix(basis), y = train$A, family = binomial())
      propensity_coef = as.numeric(propensity_fit$coefficients)
      propensity_coef[!is.finite(propensity_coef)] = 0
      list(mu1_coef = mu1_coef, mu0_coef = mu0_coef, propensity_coef = propensity_coef)
  }

  fit_target_propensity <- function (target) 
  {
      glm(A ~ x1 + x2, family = binomial(), data = target)
  }

  fit_transfer_weight_model <- function (source, target) 
  {
      pooled <- rbind(transform(source, S = 1), transform(target, S = 0))
      model <- glm(S ~ x1 + x2 + I(x1^2) + I(x1 * x2) + I(x2^2), family = binomial(), data = pooled)
      if (!isTRUE(model$converged) || any(!is.finite(coef(model)))) 
          stop("Weight model failed")
      ratio <- nrow(source)/nrow(target)
      norm <- mean(ratio * exp(-predict(model, newdata = source, type = "link")))
      stopifnot(is.finite(norm), norm > 0)
      list(model = model, sample_ratio = ratio, normalizer = norm)
  }

  make_crossfit_folds <- function (source, target, n_folds = 5, seed = 1) 
  {
      list(source = make_stratified_folds(source$A, n_folds, seed + 1), target = make_stratified_folds(target$A, 
          n_folds, seed + 2), n_folds = n_folds)
  }

  if (method %in% c("Target")) make_pooled_sample <- function (data) 
  {
      data.frame(x1 = data$x1_rwd, x2 = data$x2_rwd, A = data$A_rwd, Y = data$Y_rwd)
  }

  if (method %in% c("Pooled")) make_pooled_sample <- function (data) 
  {
      rbind(data.frame(x1 = data$x1_rct, x2 = data$x2_rct, A = data$A_rct, Y = data$Y_rct), data.frame(x1 = data$x1_rwd, 
          x2 = data$x2_rwd, A = data$A_rwd, Y = data$Y_rwd))
  }

  make_samples <- function (data) 
  {
      list(source = data.frame(x1 = data$x1_rct, x2 = data$x2_rct, A = data$A_rct, Y = data$Y_rct), 
          target = data.frame(x1 = data$x1_rwd, x2 = data$x2_rwd, A = data$A_rwd, Y = data$Y_rwd))
  }

  make_stratified_folds <- function (A, n_folds = 5, seed = 1) 
  {
      if (length(n_folds) != 1L || !is.finite(n_folds) || n_folds < 2 || n_folds != as.integer(n_folds)) 
          stop("n_folds must be an integer >= 2.")
      if (length(A)%%n_folds != 0L) 
          stop("Equal folds require the sample size to be divisible by n_folds.")
      set.seed(seed)
      sample(rep(seq_len(n_folds), each = length(A)/n_folds))
  }

  predict_basis_model <- function (fit, data) 
  {
      as.numeric(make_basis(data$x1, data$x2) %*% fit$coefficient)
  }

  if (method %in% c("WY")) predict_nuisance_scores <- function (models, weight_fit, source, target = NULL, need_target = TRUE) 
  {
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

  if (method %in% c("CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) predict_nuisance_scores <- function (models, weight_fit, source, target = NULL, need_target = TRUE) 
  {
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

  predict_pooled_score <- function (fit, test) 
  {
      basis = make_basis(test$x1, test$x2)
      mu1 = as.vector(basis %*% fit$mu1_coef)
      mu0 = as.vector(basis %*% fit$mu0_coef)
      propensity = clip(plogis(as.vector(basis %*% fit$propensity_coef)))
      mu1 - mu0 + test$A * (test$Y - mu1)/propensity - (1 - test$A) * (test$Y - mu0)/(1 - propensity)
  }

  predict_target_propensity <- function (fit, data) 
  {
      clip(as.numeric(predict(fit, newdata = data, type = "response")))
  }

  predict_transfer_weight <- function (fit, data) 
  {
      weight = fit$sample_ratio * exp(-as.numeric(predict(fit$model, newdata = data, type = "link")))/fit$normalizer
      if (any(!is.finite(weight)) || any(weight <= 0)) 
          stop("Nonfinite/underflowed density ratios: inspect overlap/model fit.")
      weight
  }

  wu2020_sample <- function (data, domain = c("rwd", "rct")) 
  {
      domain <- match.arg(domain)
      nms <- c("x1", "x2", "A", "Y")
      setNames(as.data.frame(data[paste0(nms, "_", domain)]), nms)
  }

  wu2020_target_folds <- function (n, n_folds = 5L, seed = 1L) 
  {
      if (n_folds < 2L || n%%n_folds != 0L) 
          stop("Wu2020 common-D projection requires equal target folds.")
      set.seed(seed)
      sample(rep(seq_len(n_folds), each = n/n_folds))
  }

  # --- Joint R-learning and confounding correction ---
  if (method %in% c("CC", "Fusion")) fit_cross_fitted_pointwise_confounding <- function (crossfit, lambda_tau = 0, lambda_c = 0) 
  {
      source = crossfit$samples$source
      target = crossfit$samples$target
      p_tau = ncol(make_basis(0, 0))
      p_c = ncol(make_pointwise_confounding_basis(0, 0))
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
          pointwise_fit = fit_integrative_r_fold(source_train, target_train, fold_fit$nuisance, lambda_tau = lambda_tau, 
              lambda_c = lambda_c, inner_predictions = inner_predictions)
          tau_coef_by_fold[k, ] = pointwise_fit$tau_coef
          c_coef_by_fold[k, ] = pointwise_fit$c_coef
          c_source[fold_fit$source_eval_index] = as.vector(make_pointwise_confounding_basis(source$x1[fold_fit$source_eval_index], 
              source$x2[fold_fit$source_eval_index]) %*% pointwise_fit$c_coef)
          c_target[fold_fit$target_eval_index] = as.vector(make_pointwise_confounding_basis(target$x1[fold_fit$target_eval_index], 
              target$x2[fold_fit$target_eval_index]) %*% pointwise_fit$c_coef)
      }
      list(c_coef = colMeans(c_coef_by_fold), tau_coef = colMeans(tau_coef_by_fold), c_coef_by_fold = c_coef_by_fold, 
          tau_coef_by_fold = tau_coef_by_fold, c_source = c_source, c_target = c_target)
  }

  if (method %in% c("CC_SCAD", "Fusion_SCAD")) fit_cross_fitted_pointwise_confounding <- function (crossfit, lambda_tau = NULL, lambda_c = NULL) 
  {
      source = crossfit$samples$source
      target = crossfit$samples$target
      p_tau = ncol(make_basis(0, 0))
      p_c = ncol(make_pointwise_confounding_basis(0, 0))
      tau_coef_by_fold = matrix(NA_real_, crossfit$n_folds, p_tau)
      c_coef_by_fold = matrix(NA_real_, crossfit$n_folds, p_c)
      c_source = rep(NA_real_, nrow(source))
      c_target = rep(NA_real_, nrow(target))
      tuning_by_fold = vector("list", crossfit$n_folds)
      for (k in seq_len(crossfit$n_folds)) {
          fold_fit = crossfit$fold_fits[[k]]
          source_train = source[fold_fit$source_train_index, , drop = FALSE]
          target_train = target[fold_fit$target_train_index, , drop = FALSE]
          inner_predictions = fit_inner_r_predictions(source_train, target_train, n_inner_folds = 5, 
              seed = 20000 + k)
          pointwise_fit = fit_integrative_r_fold(source_train, target_train, fold_fit$nuisance, lambda_tau = lambda_tau, 
              lambda_c = lambda_c, inner_predictions = inner_predictions, seed = 30000 + k)
          tuning_by_fold[[k]] = pointwise_fit$tuning
          tau_coef_by_fold[k, ] = pointwise_fit$tau_coef
          c_coef_by_fold[k, ] = pointwise_fit$c_coef
          c_source[fold_fit$source_eval_index] = as.vector(make_pointwise_confounding_basis(source$x1[fold_fit$source_eval_index], 
              source$x2[fold_fit$source_eval_index]) %*% pointwise_fit$c_coef)
          c_target[fold_fit$target_eval_index] = as.vector(make_pointwise_confounding_basis(target$x1[fold_fit$target_eval_index], 
              target$x2[fold_fit$target_eval_index]) %*% pointwise_fit$c_coef)
      }
      list(c_coef = colMeans(c_coef_by_fold), tau_coef = colMeans(tau_coef_by_fold), c_coef_by_fold = c_coef_by_fold, 
          tau_coef_by_fold = tau_coef_by_fold, tuning_by_fold = tuning_by_fold, c_source = c_source, 
          c_target = c_target)
  }

  fit_inner_r_predictions <- function (source, target, n_inner_folds = 5, seed = 1) 
  {
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

  # --- Fusion: estimate conditional noise and combine the two scores ---
  compute_q_weight <- function (v0, v1, w, n0, n1) 
  {
      if (any(!is.finite(c(v0, v1, w, n0, n1))) || any(c(v0, v1, w, n0, n1) <= 0)) 
          stop("q requires finite positive variances, weights, and sample sizes.")
      plogis(log(w) + log(v1) - log(n1) - log(v0) + log(n0))
  }

  fit_cross_fitted_q <- function (scores, n_inner_folds = 5) 
  {
      source = scores$samples$source
      target = scores$samples$target
      n1 = nrow(source)
      n0 = nrow(target)
      q_source = v0_source = v1_source = rep(NA_real_, n1)
      q_target = v0_target = v1_target = rep(NA_real_, n0)
      models = vector("list", scores$n_folds)
      for (k in seq_len(scores$n_folds)) {
          f = scores$fold_fits[[k]]
          tr1 = source[f$source_train_index, , drop = FALSE]
          tr0 = target[f$target_train_index, , drop = FALSE]
          i1 = f$source_eval_index
          i0 = f$target_eval_index
          ev1 = source[i1, , drop = FALSE]
          ev0 = target[i0, , drop = FALSE]
          vm = fit_q_variance_models(tr1, tr0, n_inner_folds, seed = 30000 + k)
          e0_1 = predict_target_propensity(f$nuisance$e0, ev1)
          e0_0 = predict_target_propensity(f$nuisance$e0, ev0)
          v0_source[i1] = predict_q_variance(vm$target, ev1, e0_1)
          v1_source[i1] = predict_q_variance(vm$source, ev1, rep(0.5, nrow(ev1)))
          v0_target[i0] = predict_q_variance(vm$target, ev0, e0_0)
          v1_target[i0] = predict_q_variance(vm$source, ev0, rep(0.5, nrow(ev0)))
          q_source[i1] = compute_q_weight(v0_source[i1], v1_source[i1], scores$source$w[i1], n0, n1)
          q_target[i0] = compute_q_weight(v0_target[i0], v1_target[i0], scores$target$w[i0], n0, n1)
          models[[k]] = vm
      }
      if (any(!is.finite(c(q_source, q_target)))) 
          stop("Incomplete q predictions.")
      list(q_source = q_source, q_target = q_target, v0_source = v0_source, v1_source = v1_source, 
          v0_target = v0_target, v1_target = v1_target, models = models)
  }

  fit_q_variance_models <- function (source, target, n_inner_folds = 5, seed = 1) 
  {
      folds <- make_crossfit_folds(source, target, n_inner_folds, seed)
      domain <- function(dat, fold) {
          residual <- rep(NA_real_, nrow(dat))
          for (j in seq_len(n_inner_folds)) for (a in 0:1) {
              tr <- dat[fold != j, , drop = FALSE]
              ev <- which(fold == j & dat$A == a)
              if (length(ev)) 
                  residual[ev] <- dat$Y[ev] - predict_basis_model(fit_arm_model(tr, a), dat[ev, , 
                    drop = FALSE])
          }
          stopifnot(all(is.finite(residual)))
          fits <- lapply(0:1, function(a) {
              use <- dat$A == a
              fit_basis_regression(residual[use]^2, make_variance_basis(dat$x1[use], dat$x2[use]))
          })
          list(arm0 = fits[[1]], arm1 = fits[[2]], variance_floor = max(1e-08, 1e-06 * mean(residual^2)), 
              inner_fold = fold)
      }
      list(source = domain(source, folds$source), target = domain(target, folds$target))
  }

  predict_q_variance <- function (fit, data, propensity) 
  {
      b <- make_variance_basis(data$x1, data$x2)
      s0 <- pmax(as.vector(b %*% fit$arm0$coefficient), fit$variance_floor)
      s1 <- pmax(as.vector(b %*% fit$arm1$coefficient), fit$variance_floor)
      s1/propensity + s0/(1 - propensity)
  }

  wu2020_q_fit <- function (z, A, Y) 
  {
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

  # --- Wu2020: auxiliary features ---
  wu2020_fit_features <- function (target, ntree = 500L, seed = 1L) 
  {
      d0 <- wu2020_pretrain(target, target, ntree, seed)
      agree <- target$A == d0
      if (sum(agree) < 2L || sum(!agree) < 2L) 
          stop("Need at least two agreeing and two disagreeing target observations.")
      x <- wu2020_x(target)
      list(h1 = wu2020_rf(x[agree, ], target$A[agree], ntree, seed + 101L, classification = TRUE), 
          benefit_agree = wu2020_rf(x[agree, ], target$Y[agree], ntree, seed + 102L), benefit_disagree = wu2020_rf(x[!agree, 
              ], target$Y[!agree], ntree, seed + 103L))
  }

  wu2020_predict_features <- function (fit, x1, x2) 
  {
      x <- data.frame(x1 = x1, x2 = x2)
      cbind(H1 = wu2020_rf_predict(fit$h1, x), H2 = wu2020_rf_predict(fit$benefit_agree, x) - wu2020_rf_predict(fit$benefit_disagree, 
          x))
  }

  wu2020_pretrain <- function (train, test, ntree, seed) 
  {
      if (length(unique(train$A)) != 2L) 
          stop("Step-0 training sample must contain both treatments.")
      q <- wu2020_rf(cbind(wu2020_x(train), A = train$A), train$Y, ntree, seed)
      q1 <- wu2020_rf_predict(q, cbind(wu2020_x(test), A = 1))
      q0 <- wu2020_rf_predict(q, cbind(wu2020_x(test), A = 0))
      as.integer(q1 > q0)
  }

  wu2020_rf <- function (x, y, ntree, seed, classification = FALSE) 
  {
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

  wu2020_rf_predict <- function (fit, x) 
  {
      if (!is.null(fit$constant)) 
          return(rep(fit$constant, nrow(x)))
      if (isTRUE(fit$classification)) 
          return(as.numeric(predict(fit$model, x, type = "prob")[, "1"]))
      as.numeric(predict(fit$model, x))
  }

  wu2020_with_seed <- function (seed, fun) 
  {
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

  wu2020_x <- function (z) 
  z[, c("x1", "x2"), drop = FALSE]

  # --- DC policy search and BFGS diagnostics ---
  if (method %in% c("WY", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) dc_policy <- function (Z, r, source_start = NULL, n_starts = 8, max_dc = 100, seed = 1) 
  {
      set.seed(seed)
      starts = list()
      if (!is.null(source_start)) 
          starts[[length(starts) + 1]] = source_start
      starts[[length(starts) + 1]] = rep(0, ncol(Z))
      while (length(starts) < n_starts) {
          starts[[length(starts) + 1]] = rnorm(ncol(Z), sd = 0.35)
      }
      fits = lapply(starts[seq_len(n_starts)], function(start) {
          dc_policy_once(Z, r, start = start, max_dc = max_dc)
      })
      values = vapply(fits, function(eta) empirical_value(Z, r, eta), numeric(1))
      fits[[which.max(values)]]
  }

  if (method %in% c("Target", "Pooled")) dc_policy <- function (Z, r, n_starts = 8, max_dc = 100, seed = 1) 
  {
      set.seed(seed)
      starts = list(rep(0, ncol(Z)))
      while (length(starts) < n_starts) {
          starts[[length(starts) + 1]] = rnorm(ncol(Z), sd = 0.35)
      }
      fits = lapply(starts, function(start) {
          dc_policy_once(Z, r, start = start, max_dc = max_dc)
      })
      values = vapply(fits, function(eta) empirical_value(Z, r, eta), numeric(1))
      fits[[which.max(values)]]
  }

  if (method %in% c("Wu2020")) dc_policy <- function (Z, r, n_starts = 8L, max_dc = 100L, seed = 1L) 
  {
      set.seed(seed)
      starts <- list(rep(0, ncol(Z)))
      while (length(starts) < n_starts) starts[[length(starts) + 1L]] <- rnorm(ncol(Z), sd = 0.35)
      fits <- lapply(starts, function(start) {
          dc_policy_once(Z, r, start, max_dc)
      })
      values <- vapply(fits, function(eta) empirical_value(Z, r, eta), numeric(1))
      fits[[which.max(values)]]
  }

  if (method %in% c("WY", "Target", "Pooled", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) dc_policy_once <- function (Z, r, start = NULL, max_dc = 200L, tol = 1e-08) 
  {
      .itr_last_beta <- NULL
      .itr_last_margin <- NULL
      .itr_margin <- function(beta) {
          if (!identical(beta, .itr_last_beta)) {
              .itr_last_margin <<- label * as.vector(Z %*% beta)
              .itr_last_beta <<- beta
          }
          .itr_last_margin
      }
      keep = is.finite(r) & abs(r) > 1e-12
      Z = as.matrix(Z[keep, , drop = FALSE])
      r = r[keep]
      if (length(r) == 0) 
          return(rep(0, ncol(Z)))
      label = sign(r)
      weight = abs(r)
      eta = if (is.null(start)) 
          rep(0, ncol(Z))
      else as.numeric(start)
      full_objective = function(beta) {
          margin = .itr_margin(beta)
          sum(weight * ramp_loss(margin))
      }
      for (iter in seq_len(max_dc)) {
          margin_old = .itr_margin(eta)
          concave_gradient = 0.5 * colSums(Z * (weight * l0_grad(margin_old) * label))
          convex_objective = function(beta) {
              margin = .itr_margin(beta)
              0.5 * sum(weight * l1_loss(margin)) - sum(concave_gradient * beta)
          }
          convex_gradient = function(beta) {
              margin = .itr_margin(beta)
              gradient = 0.5 * colSums(Z * (weight * l1_grad(margin) * label)) - concave_gradient
              gradient
          }
          opt = optim(eta, convex_objective, convex_gradient, method = "BFGS", control = list(maxit = 100, 
              reltol = 1e-10))
          eta_new = opt$par
          if (sum((eta_new - eta)^2) < tol || abs(full_objective(eta_new) - full_objective(eta)) < 
              tol) {
              eta = eta_new
              break
          }
          eta = eta_new
      }
      eta
  }

  if (method %in% c("Wu2020")) dc_policy_once <- function (Z, r, start = NULL, max_dc = 200L, tol = 1e-08) 
  {
      .itr_last_beta <- NULL
      .itr_last_margin <- NULL
      .itr_margin <- function(beta) {
          if (!identical(beta, .itr_last_beta)) {
              .itr_last_margin <<- label * as.vector(Z %*% beta)
              .itr_last_beta <<- beta
          }
          .itr_last_margin
      }
      keep <- is.finite(r) & abs(r) > 1e-12
      Z <- as.matrix(Z[keep, , drop = FALSE])
      r <- r[keep]
      if (!length(r)) 
          return(rep(0, ncol(Z)))
      label <- sign(r)
      weight <- abs(r)
      eta <- if (is.null(start)) 
          rep(0, ncol(Z))
      else as.numeric(start)
      full_objective <- function(beta) {
          margin <- .itr_margin(beta)
          sum(weight * ramp_loss(margin))
      }
      for (iter in seq_len(max_dc)) {
          margin_old <- .itr_margin(eta)
          concave_gradient <- 0.5 * colSums(Z * (weight * l0_grad(margin_old) * label))
          convex_objective <- function(beta) {
              margin <- .itr_margin(beta)
              0.5 * sum(weight * l1_loss(margin)) - sum(concave_gradient * beta)
          }
          convex_gradient <- function(beta) {
              margin <- .itr_margin(beta)
              gradient <- 0.5 * colSums(Z * (weight * l1_grad(margin) * label)) - concave_gradient
              gradient
          }
          opt <- optim(eta, convex_objective, convex_gradient, method = "BFGS", control = list(maxit = 100, 
              reltol = 1e-10))
          eta_new <- opt$par
          if (sum((eta_new - eta)^2) < tol || abs(full_objective(eta_new) - full_objective(eta)) < 
              tol) {
              eta <- eta_new
              break
          }
          eta <- eta_new
      }
      eta
  }

  empirical_value <- function (Z, r, eta) 
  {
      sum(r * as.numeric(as.vector(Z %*% eta) > 0))
  }

  if (method %in% c("WY", "Target", "Pooled", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) l0_grad <- function (u) 
  {
      ifelse(u >= 0, 0, ifelse(u >= -1, 2 * u, -2))
  }

  if (method %in% c("Wu2020")) l0_grad <- function (u) 
  ifelse(u >= 0, 0, ifelse(u >= -1, 2 * u, -2))

  if (method %in% c("WY", "Target", "Pooled", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) l0_loss <- function (u) 
  {
      ifelse(u >= 0, 0, ifelse(u >= -1, u^2, -2 * u - 1))
  }

  if (method %in% c("Wu2020")) l0_loss <- function (u) 
  ifelse(u >= 0, 0, ifelse(u >= -1, u^2, -2 * u - 1))

  if (method %in% c("WY", "Target", "Pooled", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) l1_grad <- function (u) 
  {
      ifelse(u >= 1, 0, ifelse(u >= 0, -2 * (1 - u), -2))
  }

  if (method %in% c("Wu2020")) l1_grad <- function (u) 
  ifelse(u >= 1, 0, ifelse(u >= 0, -2 * (1 - u), -2))

  if (method %in% c("WY", "Target", "Pooled", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) l1_loss <- function (u) 
  {
      ifelse(u >= 1, 0, ifelse(u >= 0, (1 - u)^2, 1 - 2 * u))
  }

  if (method %in% c("Wu2020")) l1_loss <- function (u) 
  ifelse(u >= 1, 0, ifelse(u >= 0, (1 - u)^2, 1 - 2 * u))

  optim <- function (...) 
  {
      args <- list(...)
      if (identical(args$method, "BFGS")) {
          args$control$maxit <- 200L
          args$control$reltol <- 1e-10
      }
      ans <- do.call(stats::optim, args)
      diagnostics$optimizer_calls <- diagnostics$optimizer_calls + 1L
      if (ans$convergence != 0) 
          diagnostics$optimizer_nonconvergence <- diagnostics$optimizer_nonconvergence + 1L
      ans
  }

  if (method %in% c("WY", "Target", "Pooled", "IP", "CC", "Fusion", "CC_SCAD", "Fusion_SCAD")) ramp_loss <- function (u) 
  {
      0.5 * (l1_loss(u) - l0_loss(u))
  }

  if (method %in% c("Wu2020")) ramp_loss <- function (u) 
  0.5 * (l1_loss(u) - l0_loss(u))

  # Count solver starts without changing arguments, RNG or fitted values.
  solve_once <- dc_policy_once
  dc_policy_once <- function(Z,r,start=NULL,max_dc=200L,...) {
    stopifnot(max_dc==200L)
    diagnostics$starts <- diagnostics$starts+1L
    solve_once(Z,r,start,max_dc,...)
  }
  environment()
}
