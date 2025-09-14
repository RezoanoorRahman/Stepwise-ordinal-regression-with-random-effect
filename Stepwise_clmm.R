library(ordinal)
library(dplyr)

stepwise_ordinal <- function(Df, Y,
                                         fixed_effect_variables,
                                         random_effect_variables,
                                         link = "logit",
                                         nAGQ = 5,
                                         min_improve = 1e-8,
                                         max_steps = 100,
                                         verbose = TRUE) {
  
  stopifnot(is.character(Y), length(Y) == 1)
  if (!inherits(Df[[Y]], "ordered")) stop("Response must be an ordered factor.")
  
  make_fixed_rhs <- function(vars) if (length(vars) == 0) "1" else paste(vars, collapse = " + ")
  make_rand_rhs  <- function(rand) if (length(rand) == 0) "" else paste0(" + ", paste(sprintf("(1 | %s)", rand), collapse = " + "))
  make_formula   <- function(y, fixed_vars, rand_vars) as.formula(paste0(y, " ~ ", make_fixed_rhs(fixed_vars), make_rand_rhs(rand_vars)))
  
  old_variables <- c(Y)           
  fixed_in      <- character(0)
  random_in     <- character(0)
  
  old_model <- clm(make_formula(Y, fixed_in, random_in), data = Df, link = link)
  old_aic   <- AIC(old_model)
  
  best_overall_improvement <- Inf
  step <- 0
  
  while ( (length(fixed_effect_variables) + length(random_effect_variables)) > 0 &&
          (is.infinite(best_overall_improvement) || best_overall_improvement > min_improve) &&
          step < max_steps ) {
    
    step <- step + 1
    if (verbose) message(sprintf("\nstep %d | current AIC: %.3f", step, old_aic))
    
    fixed_models       <- vector("list", length(fixed_effect_variables))
    fixed_improvements <- rep(-10, length(fixed_effect_variables))
    
    for (i in seq_along(fixed_effect_variables)) {
      
      select_variables_fixed <- c(old_variables, fixed_effect_variables[i])
      
      f_fixed <- make_formula(Y, c(fixed_in, fixed_effect_variables[i]), random_in)
      
      if (length(random_in) > 0) {
        fixed_model <- try(clmm(f_fixed, data = Df, link = link, nAGQ = nAGQ), silent = TRUE)
      } else {
        fixed_model <- try(clm (f_fixed, data = Df, link = link), silent = TRUE)
      }
      
      if (!inherits(fixed_model, "try-error")) {
        fixed_improvements[i] <- old_aic - AIC(fixed_model)
        fixed_models[[i]]     <- fixed_model
      }
    }
    
    max_aic_improved_fixed <- ifelse(all(is.infinite(fixed_improvements)), -Inf, max(fixed_improvements, na.rm = TRUE))
    best_fixed_model <- if (is.finite(max_aic_improved_fixed)) fixed_models[[ which.max(fixed_improvements) ]] else NULL
    best_fixed_var   <- if (is.finite(max_aic_improved_fixed)) fixed_effect_variables[ which.max(fixed_improvements) ] else NA_character_
   
    random_models       <- vector("list", length(random_effect_variables))
    random_improvements <- rep(-10, length(random_effect_variables))
    
    for (i in seq_along(random_effect_variables)) {
      
      select_variables_random <- c(old_variables, random_effect_variables[i])
      
      f_rand <- make_formula(Y, fixed_in, c(random_in, random_effect_variables[i]))
      random_model <- try(clmm(f_rand, data = Df, link = link, nAGQ = nAGQ), silent = TRUE)
      
      if (!inherits(random_model, "try-error")) {
        random_improvements[i] <- old_aic - AIC(random_model)
        random_models[[i]]     <- random_model
      }
    }
    
    max_aic_improved_random <- ifelse(all(is.infinite(random_improvements)), -Inf, max(random_improvements, na.rm = TRUE))
    best_random_model <- if (is.finite(max_aic_improved_random)) random_models[[ which.max(random_improvements) ]] else NULL
    best_random_var   <- if (is.finite(max_aic_improved_random)) random_effect_variables[ which.max(random_improvements) ] else NA_character_
 
    diff_fixed_over_random   <- max_aic_improved_fixed - max_aic_improved_random
    best_overall_improvement <- max(c(max_aic_improved_fixed, max_aic_improved_random), na.rm = TRUE)
    
    if (!is.finite(best_overall_improvement) || best_overall_improvement <= min_improve) break
    
    if (diff_fixed_over_random > 0) {
      new_model <- best_fixed_model

      fixed_in          <- c(fixed_in, best_fixed_var)
      old_variables     <- c(old_variables, best_fixed_var)
      fixed_effect_variables <- setdiff(fixed_effect_variables, best_fixed_var)
      
      if (verbose) message(sprintf("  + add fixed: %s | Delta_AIC = %.3f", best_fixed_var, best_overall_improvement))
    } else {
      new_model <- best_random_model
      random_in          <- c(random_in, best_random_var)
      old_variables      <- c(old_variables, best_random_var)
      random_effect_variables <- setdiff(random_effect_variables, best_random_var)
      
      if (verbose) message(sprintf("  + add random: (1|%s) | Delta_AIC = %.3f", best_random_var, best_overall_improvement))
    }
    

    old_model <- new_model
    old_aic   <- AIC(old_model)
  }
  
  list(
    final_fit   = old_model,
    final_AIC   = old_aic,
    formula     = formula(old_model),
    fixed_in    = fixed_in,
    random_in   = random_in
  )
}

