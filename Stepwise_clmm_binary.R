library(lme4)
library(dplyr)

stepwise_binomial <- function(Df, Y,
                              fixed_effect_variables,
                              random_effect_variables,
                              link = "logit",
                              nAGQ = 5,
                              min_improve = 1e-8,
                              max_steps = 100,
                              verbose = TRUE) {
  
  stopifnot(is.character(Y), length(Y) == 1)
  
  # ---- minimal hygiene ----
  # keep only columns we need, drop NAs, drop empty levels
  need <- unique(c(Y, fixed_effect_variables, random_effect_variables))
  Df <- Df[, need, drop = FALSE]
  Df <- Df %>% filter(complete.cases(across(everything())))
  Df[] <- lapply(Df, function(z) if (is.character(z)) factor(z) else z)
  Df <- droplevels(Df)
  
  # ensure binary response (0/1 numeric or 2-level factor)
  yv <- Df[[Y]]
  if (is.factor(yv)) {
    if (nlevels(yv) != 2) stop("Response must be binary (2 levels).")
    Df[[Y]] <- relevel(yv, ref = levels(yv)[1])
  } else if (!all(yv %in% c(0,1))) {
    stop("Response must be 0/1 or a 2-level factor.")
  }
  
  # drop fixed/random with <2 levels after droplevels
  fixed_effect_variables  <- fixed_effect_variables [vapply(fixed_effect_variables,  function(v) length(unique(Df[[v]])) >= 2, logical(1))]
  random_effect_variables <- random_effect_variables[vapply(random_effect_variables, function(v) length(unique(Df[[v]])) >= 2, logical(1))]
  
  options(contrasts = c("contr.treatment","contr.poly"))
  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  
  make_fixed_rhs <- function(vars) if (length(vars) == 0) "1" else paste(vars, collapse = " + ")
  make_rand_rhs  <- function(rand) if (length(rand) == 0) "" else paste0(" + ", paste(sprintf("(1 | %s)", rand), collapse = " + "))
  make_formula   <- function(y, fixed_vars, rand_vars) as.formula(paste0(y, " ~ ", make_fixed_rhs(fixed_vars), make_rand_rhs(rand_vars)))
  
  old_variables <- c(Y)
  fixed_in  <- character(0)
  random_in <- character(0)
  
  # start: intercept-only GLM
  old_model <- glm(make_formula(Y, fixed_in, character(0)), data = Df, family = binomial(link))
  old_aic   <- AIC(old_model)
  
  best_overall_improvement <- Inf
  step <- 0
  
  while ((length(fixed_effect_variables) + length(random_effect_variables)) > 0 &&
         (is.infinite(best_overall_improvement) || best_overall_improvement > min_improve) &&
         step < max_steps) {
    
    step <- step + 1
    if (verbose) message(sprintf("\nstep %d | current AIC: %.3f", step, old_aic))
    
    # try adding one fixed effect
    fixed_models       <- vector("list", length(fixed_effect_variables))
    fixed_improvements <- rep(-Inf, length(fixed_effect_variables))
    
    for (i in seq_along(fixed_effect_variables)) {
      f_fixed <- make_formula(Y, c(fixed_in, fixed_effect_variables[i]), random_in)
      fixed_model <- try(
        if (length(random_in) > 0)
          glmer(f_fixed, data = Df, family = binomial(link), nAGQ = nAGQ, control = ctrl)
        else
          glm  (f_fixed, data = Df, family = binomial(link)),
        silent = TRUE
      )
      if (!inherits(fixed_model, "try-error")) {
        fixed_improvements[i] <- old_aic - AIC(fixed_model)
        fixed_models[[i]]     <- fixed_model
      }
    }
    
    max_aic_improved_fixed <- ifelse(all(is.infinite(fixed_improvements)), -Inf, max(fixed_improvements, na.rm = TRUE))
    best_fixed_model <- if (is.finite(max_aic_improved_fixed)) fixed_models[[ which.max(fixed_improvements) ]] else NULL
    best_fixed_var   <- if (is.finite(max_aic_improved_fixed)) fixed_effect_variables[ which.max(fixed_improvements) ] else NA_character_
    
    # try adding one random intercept
    random_models       <- vector("list", length(random_effect_variables))
    random_improvements <- rep(-Inf, length(random_effect_variables))
    
    for (i in seq_along(random_effect_variables)) {
      f_rand <- make_formula(Y, fixed_in, c(random_in, random_effect_variables[i]))
      rand_model <- try(glmer(f_rand, data = Df, family = binomial(link), nAGQ = nAGQ, control = ctrl), silent = TRUE)
      if (!inherits(rand_model, "try-error")) {
        random_improvements[i] <- old_aic - AIC(rand_model)
        random_models[[i]]     <- rand_model
      }
    }
    
    max_aic_improved_random <- ifelse(all(is.infinite(random_improvements)), -Inf, max(random_improvements, na.rm = TRUE))
    best_random_model <- if (is.finite(max_aic_improved_random)) random_models[[ which.max(random_improvements) ]] else NULL
    best_random_var   <- if (is.finite(max_aic_improved_random)) random_effect_variables[ which.max(random_improvements) ] else NA_character_
    
    # choose best addition
    best_overall_improvement <- max(c(max_aic_improved_fixed, max_aic_improved_random), na.rm = TRUE)
    if (!is.finite(best_overall_improvement) || best_overall_improvement <= min_improve) break
    
    if (max_aic_improved_fixed >= max_aic_improved_random) {
      old_model <- best_fixed_model
      old_aic   <- AIC(old_model)
      fixed_in  <- c(fixed_in, best_fixed_var)
      fixed_effect_variables <- setdiff(fixed_effect_variables, best_fixed_var)
      if (verbose) message(sprintf("  + add FIXED: %s | ΔAIC = %.3f", best_fixed_var, best_overall_improvement))
    } else {
      old_model <- best_random_model
      old_aic   <- AIC(old_model)
      random_in <- c(random_in, best_random_var)
      random_effect_variables <- setdiff(random_effect_variables, best_random_var)
      if (verbose) message(sprintf("  + add RANDOM: (1|%s) | ΔAIC = %.3f", best_random_var, best_overall_improvement))
    }
  }
  
  list(
    final_fit = old_model,
    final_AIC = old_aic,
    formula   = formula(old_model),
    fixed_in  = fixed_in,
    random_in = random_in
  )
}
