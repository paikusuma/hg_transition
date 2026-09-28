library(tidyverse)
library(mvtnorm)
library(optimx)
library(furrr)
library(progressr)
library(edgeR)

# —— OU fitting function

fit_gene_limma <- function(gid,
                           limma_adjusted,
                           expr_adj,
                           meta,
                           tree,
                           n_restarts = 20) {
  
  # ── Extract gene row
  gene_lm <- limma_adjusted |> dplyr::filter(gene_id == gid)
  if (nrow(gene_lm) == 0) {
    return(tibble(gene_id = gid, converged = FALSE,
                  fail_reason = "gene not in limma_adjusted"))
  }
  if (!gid %in% rownames(expr_adj)) {
    return(tibble(gene_id = gid, converged = FALSE,
                  fail_reason = "gene not in expr_adj"))
  }
  
  # ── Per-sample residualised expression
  expr_vec <- expr_adj[gid, ]
  
  idx_pb <- which(meta$population == "PB")
  idx_pt <- which(meta$population == "PT")
  idx_ldy <- which(meta$population == "LDY")
  
  x_pb <- expr_vec[idx_pb]
  x_pt <- expr_vec[idx_pt]
  x_ldy <- expr_vec[idx_ldy]
  
  if (any(c(length(x_pb), length(x_pt), length(x_ldy)) < 2)) {
    return(tibble(gene_id = gid, converged = FALSE,
                  fail_reason = "too few samples in a population"))
  }
  if (!all(is.finite(c(x_pb, x_pt, x_ldy)))) {
    return(tibble(gene_id = gid, converged = FALSE,
                  fail_reason = "non-finite expression values"))
  }
  
  all_x <- c(x_pb, x_pt, x_ldy)
  data_mean <- mean(all_x)
  data_sd <- sd(all_x) + 1e-6
  n_obs <- length(all_x)
  T_root <- tree$T_root
  
  obs_means <- c(PB = mean(x_pb),
                 PT = mean(x_pt),
                 LDY = mean(x_ldy))
  
  # ── Alpha bounds scaled to tree
  # log_alpha_lo <- log(log(2) / (T_root * 100))  # half-life = 100 × T_root
  # log_alpha_hi <- log(log(2) / (T_root * 0.01)) # half-life = 0.01 × T_root
  
  # Replace the alpha bound calculation with fixed wider bounds
  log_alpha_lo <- log(log(2) / (T_root_gen * 100))  # slow OU
  log_alpha_hi <- 2                                 # exp(2)=7.39, half-life=0.09
  
  # ── OU tip variance helper
  ou_var <- function(alpha, sigma2, T) {
    (sigma2 / (2 * alpha)) * (1 - exp(-2 * alpha * T))
  }
  
  # ── Optimiser helper
  run_optimx <- function(init_fn, fn, lower, upper) {
    best_nll <- Inf
    best_par <- NULL
    for (i in seq_len(n_restarts)) {
      init <- init_fn()
      res <- tryCatch(
        optimx(
          par = init,
          fn = fn,
          method = c("L-BFGS-B", "Nelder-Mead"),
          lower = lower,
          upper = upper,
          control = list(maxit = 3000, kkt = FALSE)
        ),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        best_row <- res |> dplyr::arrange(value) |> dplyr::slice(1)
        if (isTRUE(best_row$value < best_nll)) {
          best_nll <- best_row$value
          best_par <- best_row
        }
      }
    }
    list(nll = best_nll, par = best_par)
  }
  
  # ── Model 1: BM
  # params: log_sigma2, ancestral,
  #         log_v_within_pb, log_v_within_pt, log_v_within_ldy
  bm_nll <- function(par) {
    sigma2 <- exp(par[["log_sigma2"]])
    ancestral <- par[["ancestral"]]
    vw_pb <- exp(par[["log_vw_pb"]])
    vw_pt <- exp(par[["log_vw_pt"]])
    vw_ldy <- exp(par[["log_vw_ldy"]])
    
    if (!all(is.finite(c(sigma2, ancestral,
                         vw_pb, vw_pt, vw_ldy)))) return(1e12)
    if (sigma2 < 1e-9) return(1e12)
    
    v_phy <- sigma2 * T_root
    ll <- sum(dnorm(x_pb,  ancestral, sqrt(v_phy + vw_pb),  log = TRUE)) +
      sum(dnorm(x_pt,  ancestral, sqrt(v_phy + vw_pt),  log = TRUE)) +
      sum(dnorm(x_ldy, ancestral, sqrt(v_phy + vw_ldy), log = TRUE))
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  # ── Model 2: OU shared alpha + shared theta
  # params: log_alpha, theta, log_sigma2,
  #         log_vw_pb, log_vw_pt, log_vw_ldy
  ou_shared_nll <- function(par) {
    alpha <- exp(par[["log_alpha"]])
    theta <- par[["theta"]]
    sigma2 <- exp(par[["log_sigma2"]])
    vw_pb <- exp(par[["log_vw_pb"]])
    vw_pt <- exp(par[["log_vw_pt"]])
    vw_ldy <- exp(par[["log_vw_ldy"]])
    
    if (!all(is.finite(c(alpha, theta, sigma2,
                         vw_pb, vw_pt, vw_ldy)))) return(1e12)
    if (alpha < 1e-9 || sigma2 < 1e-9) return(1e12)
    
    v_phy <- ou_var(alpha, sigma2, T_root)
    ll <- sum(dnorm(x_pb,  theta, sqrt(v_phy + vw_pb),  log = TRUE)) +
      sum(dnorm(x_pt,  theta, sqrt(v_phy + vw_pt),  log = TRUE)) +
      sum(dnorm(x_ldy, theta, sqrt(v_phy + vw_ldy), log = TRUE))
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  # ── Model 3: OU shared alpha + free theta
  # params: log_alpha, theta_pb, theta_pt, theta_ldy, log_sigma2,
  #         log_vw_pb, log_vw_pt, log_vw_ldy
  ou_ftheta_nll <- function(par) {
    alpha <- exp(par[["log_alpha"]])
    theta_pb <- par[["theta_pb"]]
    theta_pt <- par[["theta_pt"]]
    theta_ldy <- par[["theta_ldy"]]
    sigma2 <- exp(par[["log_sigma2"]])
    vw_pb <- exp(par[["log_vw_pb"]])
    vw_pt <- exp(par[["log_vw_pt"]])
    vw_ldy <- exp(par[["log_vw_ldy"]])
    
    if (!all(is.finite(c(alpha, theta_pb, theta_pt, theta_ldy,
                         sigma2, vw_pb, vw_pt, vw_ldy)))) return(1e12)
    if (alpha < 1e-9 || sigma2 < 1e-9) return(1e12)
    
    v_phy <- ou_var(alpha, sigma2, T_root)
    ll <- sum(dnorm(x_pb,  theta_pb,  sqrt(v_phy + vw_pb),  log = TRUE)) +
      sum(dnorm(x_pt,  theta_pt,  sqrt(v_phy + vw_pt),  log = TRUE)) +
      sum(dnorm(x_ldy, theta_ldy, sqrt(v_phy + vw_ldy), log = TRUE))
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  # ── Model 4: OU free alpha + shared theta
  # params: log_alpha_pb, log_alpha_pt, log_alpha_ldy,
  #         theta, log_sigma2,
  #         log_vw_pb, log_vw_pt, log_vw_ldy
  ou_falpha_nll <- function(par) {
    alpha_pb <- exp(par[["log_alpha_pb"]])
    alpha_pt <- exp(par[["log_alpha_pt"]])
    alpha_ldy <- exp(par[["log_alpha_ldy"]])
    theta <- par[["theta"]]
    sigma2 <- exp(par[["log_sigma2"]])
    vw_pb <- exp(par[["log_vw_pb"]])
    vw_pt <- exp(par[["log_vw_pt"]])
    vw_ldy <- exp(par[["log_vw_ldy"]])
    
    if (!all(is.finite(c(alpha_pb, alpha_pt, alpha_ldy,
                         theta, sigma2,
                         vw_pb, vw_pt, vw_ldy)))) return(1e12)
    if (any(c(alpha_pb, alpha_pt, alpha_ldy) < 1e-9) ||
        sigma2 < 1e-9) return(1e12)
    
    ll <- sum(dnorm(x_pb,  theta, sqrt(ou_var(alpha_pb,  sigma2, T_root) + vw_pb),  log = TRUE)) +
      sum(dnorm(x_pt,  theta, sqrt(ou_var(alpha_pt,  sigma2, T_root) + vw_pt),  log = TRUE)) +
      sum(dnorm(x_ldy, theta, sqrt(ou_var(alpha_ldy, sigma2, T_root) + vw_ldy), log = TRUE))
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  # ── Model 5: OU free alpha + free theta
  # params: log_alpha_pb, log_alpha_pt, log_alpha_ldy,
  #         theta_pb, theta_pt, theta_ldy, log_sigma2,
  #         log_vw_pb, log_vw_pt, log_vw_ldy
  ou_fboth_nll <- function(par) {
    alpha_pb <- exp(par[["log_alpha_pb"]])
    alpha_pt <- exp(par[["log_alpha_pt"]])
    alpha_ldy <- exp(par[["log_alpha_ldy"]])
    theta_pb <- par[["theta_pb"]]
    theta_pt <- par[["theta_pt"]]
    theta_ldy <- par[["theta_ldy"]]
    sigma2 <- exp(par[["log_sigma2"]])
    vw_pb <- exp(par[["log_vw_pb"]])
    vw_pt <- exp(par[["log_vw_pt"]])
    vw_ldy <- exp(par[["log_vw_ldy"]])
    
    if (!all(is.finite(c(alpha_pb, alpha_pt, alpha_ldy,
                         theta_pb, theta_pt, theta_ldy,
                         sigma2, vw_pb, vw_pt, vw_ldy)))) return(1e12)
    if (any(c(alpha_pb, alpha_pt, alpha_ldy) < 1e-9) ||
        sigma2 < 1e-9) return(1e12)
    
    ll <- sum(dnorm(x_pb,  theta_pb,  sqrt(ou_var(alpha_pb,  sigma2, T_root) + vw_pb),  log = TRUE)) +
      sum(dnorm(x_pt,  theta_pt,  sqrt(ou_var(alpha_pt,  sigma2, T_root) + vw_pt),  log = TRUE)) +
      sum(dnorm(x_ldy, theta_ldy, sqrt(ou_var(alpha_ldy, sigma2, T_root) + vw_ldy), log = TRUE))
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  # ── Initial value helpers
  # Use observed within-population variance as starting point for vw
  log_vw_init <- function() {
    c(log_vw_pb = log(var(x_pb)  + 1e-6) + rnorm(1, 0, 0.5),
      log_vw_pt = log(var(x_pt)  + 1e-6) + rnorm(1, 0, 0.5),
      log_vw_ldy = log(var(x_ldy) + 1e-6) + rnorm(1, 0, 0.5))
  }
  
  # ── Fit all 5 models
  bm_fit <- run_optimx(
    init_fn = function() c(
      log_sigma2 = runif(1, -4, 2),
      ancestral = rnorm(1, data_mean, data_sd),
      log_vw_init()
    ),
    fn = bm_nll,
    lower = c(-10, -Inf, -10, -10, -10),
    upper = c(  4,  Inf,   4,   4,   4)
  )
  
  ou_shared_fit <- run_optimx(
    init_fn = function() c(
      log_alpha = runif(1, log_alpha_lo, log_alpha_hi),
      theta = rnorm(1, data_mean, data_sd),
      log_sigma2 = runif(1, -4, 2),
      log_vw_init()
    ),
    fn = ou_shared_nll,
    lower = c(log_alpha_lo, -Inf, -10, -10, -10, -10),
    upper = c(log_alpha_hi,  Inf,   4,   4,   4,   4)
  )
  
  ou_ftheta_fit <- run_optimx(
    init_fn = function() c(
      log_alpha = runif(1, log_alpha_lo, log_alpha_hi),
      theta_pb = rnorm(1, obs_means["PB"],  data_sd),
      theta_pt = rnorm(1, obs_means["PT"],  data_sd),
      theta_ldy = rnorm(1, obs_means["LDY"], data_sd),
      log_sigma2 = runif(1, -4, 2),
      log_vw_init()
    ),
    fn = ou_ftheta_nll,
    lower = c(log_alpha_lo, -Inf, -Inf, -Inf, -10, -10, -10, -10),
    upper = c(log_alpha_hi,  Inf,  Inf,  Inf,   4,   4,   4,   4)
  )
  
  ou_falpha_fit <- run_optimx(
    init_fn = function() c(
      log_alpha_pb = runif(1, log_alpha_lo, log_alpha_hi),
      log_alpha_pt = runif(1, log_alpha_lo, log_alpha_hi),
      log_alpha_ldy = runif(1, log_alpha_lo, log_alpha_hi),
      theta = rnorm(1, data_mean, data_sd),
      log_sigma2 = runif(1, -4, 2),
      log_vw_init()
    ),
    fn = ou_falpha_nll,
    lower = c(log_alpha_lo, log_alpha_lo, log_alpha_lo, -Inf, -10, -10, -10, -10),
    upper = c(log_alpha_hi, log_alpha_hi, log_alpha_hi,  Inf,   4,   4,   4,   4)
  )
  
  ou_fboth_fit <- run_optimx(
    init_fn = function() c(
      log_alpha_pb = runif(1, log_alpha_lo, log_alpha_hi),
      log_alpha_pt = runif(1, log_alpha_lo, log_alpha_hi),
      log_alpha_ldy = runif(1, log_alpha_lo, log_alpha_hi),
      theta_pb = rnorm(1, obs_means["PB"],  data_sd),
      theta_pt = rnorm(1, obs_means["PT"],  data_sd),
      theta_ldy = rnorm(1, obs_means["LDY"], data_sd),
      log_sigma2 = runif(1, -4, 2),
      log_vw_init()
    ),
    fn = ou_fboth_nll,
    lower = c(log_alpha_lo, log_alpha_lo, log_alpha_lo,
              -Inf, -Inf, -Inf, -10, -10, -10, -10),
    upper = c(log_alpha_hi, log_alpha_hi, log_alpha_hi,
              Inf,  Inf,  Inf,   4,   4,   4,   4)
  )
  
  # ── Convergence check
  fits <- list(bm = bm_fit,
               ou_s = ou_shared_fit,
               ou_ft = ou_ftheta_fit,
               ou_fa = ou_falpha_fit,
               ou_fb = ou_fboth_fit)
  
  fits_ok <- !sapply(fits, function(f) is.null(f$par))
  if (!all(fits_ok)) {
    return(tibble(gene_id = gid, converged = FALSE,
                  fail_reason = paste("failed models:",
                                      paste(names(fits)[!fits_ok],
                                            collapse = ","))))
  }
  
  # ── Extract parameters
  bm_sigma2 <- exp(as.numeric(bm_fit$par[["log_sigma2"]]))
  bm_ancestral <- as.numeric(bm_fit$par[["ancestral"]])
  bm_vw_pb <- exp(as.numeric(bm_fit$par[["log_vw_pb"]]))
  bm_vw_pt <- exp(as.numeric(bm_fit$par[["log_vw_pt"]]))
  bm_vw_ldy <- exp(as.numeric(bm_fit$par[["log_vw_ldy"]]))
  bm_ll <- -bm_fit$nll
  
  ou_s_alpha <- exp(as.numeric(ou_shared_fit$par[["log_alpha"]]))
  ou_s_theta <- as.numeric(ou_shared_fit$par[["theta"]])
  ou_s_sigma2 <- exp(as.numeric(ou_shared_fit$par[["log_sigma2"]]))
  ou_s_vw_pb <- exp(as.numeric(ou_shared_fit$par[["log_vw_pb"]]))
  ou_s_vw_pt <- exp(as.numeric(ou_shared_fit$par[["log_vw_pt"]]))
  ou_s_vw_ldy <- exp(as.numeric(ou_shared_fit$par[["log_vw_ldy"]]))
  ou_s_ll <- -ou_shared_fit$nll
  
  ou_ft_alpha <- exp(as.numeric(ou_ftheta_fit$par[["log_alpha"]]))
  ou_ft_theta_pb <- as.numeric(ou_ftheta_fit$par[["theta_pb"]])
  ou_ft_theta_pt <- as.numeric(ou_ftheta_fit$par[["theta_pt"]])
  ou_ft_theta_ldy <- as.numeric(ou_ftheta_fit$par[["theta_ldy"]])
  ou_ft_sigma2 <- exp(as.numeric(ou_ftheta_fit$par[["log_sigma2"]]))
  ou_ft_vw_pb <- exp(as.numeric(ou_ftheta_fit$par[["log_vw_pb"]]))
  ou_ft_vw_pt <- exp(as.numeric(ou_ftheta_fit$par[["log_vw_pt"]]))
  ou_ft_vw_ldy <- exp(as.numeric(ou_ftheta_fit$par[["log_vw_ldy"]]))
  ou_ft_ll <- -ou_ftheta_fit$nll
  
  ou_fa_alpha_pb <- exp(as.numeric(ou_falpha_fit$par[["log_alpha_pb"]]))
  ou_fa_alpha_pt <- exp(as.numeric(ou_falpha_fit$par[["log_alpha_pt"]]))
  ou_fa_alpha_ldy <- exp(as.numeric(ou_falpha_fit$par[["log_alpha_ldy"]]))
  ou_fa_theta <- as.numeric(ou_falpha_fit$par[["theta"]])
  ou_fa_sigma2 <- exp(as.numeric(ou_falpha_fit$par[["log_sigma2"]]))
  ou_fa_vw_pb <- exp(as.numeric(ou_falpha_fit$par[["log_vw_pb"]]))
  ou_fa_vw_pt <- exp(as.numeric(ou_falpha_fit$par[["log_vw_pt"]]))
  ou_fa_vw_ldy <- exp(as.numeric(ou_falpha_fit$par[["log_vw_ldy"]]))
  ou_fa_ll <- -ou_falpha_fit$nll
  
  ou_fb_alpha_pb <- exp(as.numeric(ou_fboth_fit$par[["log_alpha_pb"]]))
  ou_fb_alpha_pt <- exp(as.numeric(ou_fboth_fit$par[["log_alpha_pt"]]))
  ou_fb_alpha_ldy <- exp(as.numeric(ou_fboth_fit$par[["log_alpha_ldy"]]))
  ou_fb_theta_pb <- as.numeric(ou_fboth_fit$par[["theta_pb"]])
  ou_fb_theta_pt <- as.numeric(ou_fboth_fit$par[["theta_pt"]])
  ou_fb_theta_ldy <- as.numeric(ou_fboth_fit$par[["theta_ldy"]])
  ou_fb_sigma2 <- exp(as.numeric(ou_fboth_fit$par[["log_sigma2"]]))
  ou_fb_vw_pb <- exp(as.numeric(ou_fboth_fit$par[["log_vw_pb"]]))
  ou_fb_vw_pt <- exp(as.numeric(ou_fboth_fit$par[["log_vw_pt"]]))
  ou_fb_vw_ldy <- exp(as.numeric(ou_fboth_fit$par[["log_vw_ldy"]]))
  ou_fb_ll <- -ou_fboth_fit$nll
  
  # ── AIC / BIC
  # BM        : 5 params (sigma2, ancestral, vw_pb, vw_pt, vw_ldy)
  # OU_shared : 6 params (alpha, theta, sigma2, vw_pb, vw_pt, vw_ldy)
  # OU_ftheta : 8 params (alpha, theta×3, sigma2, vw×3)
  # OU_falpha : 8 params (alpha×3, theta, sigma2, vw×3)
  # OU_fboth  : 10 params (alpha×3, theta×3, sigma2, vw×3)
  k <- c(bm = 5, ou_s = 6, ou_ft = 8, ou_fa = 8, ou_fb = 10)
  
  ll_vec <- c(bm = bm_ll,
              ou_s = ou_s_ll,
              ou_ft = ou_ft_ll,
              ou_fa = ou_fa_ll,
              ou_fb = ou_fb_ll)
  
  aic_vec <- 2 * k - 2 * ll_vec
  bic_vec <- k * log(n_obs) - 2 * ll_vec
  
  names(aic_vec) <- c("BM", "OU_shared", "OU_ftheta", "OU_falpha", "OU_fboth")
  names(bic_vec) <- names(aic_vec)
  
  best_model_aic <- names(which.min(aic_vec))
  best_model_bic <- names(which.min(bic_vec))
  
  # ── LRTs
  # BM → OU_shared  : df = 1 (add alpha)
  # OU_shared → OU_ftheta : df = 2 (free theta ×3 replaces theta ×1)
  # OU_shared → OU_falpha : df = 2 (free alpha ×3 replaces alpha ×1)
  # OU_ftheta → OU_fboth  : df = 2 (free alpha)
  # OU_falpha → OU_fboth  : df = 2 (free theta)
  # BM → OU_fboth         : df = 5
  
  lrt <- function(ll_full, ll_null, df) {
    stat <- max(2 * (ll_full - ll_null), 0)
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    c(stat = stat, pval = pval)
  }
  
  lrt_s_bm <- lrt(ou_s_ll,  bm_ll,    df = 1)
  lrt_ft_s <- lrt(ou_ft_ll, ou_s_ll,  df = 2)
  lrt_fa_s <- lrt(ou_fa_ll, ou_s_ll,  df = 2)
  lrt_fb_ft <- lrt(ou_fb_ll, ou_ft_ll, df = 2)
  lrt_fb_fa <- lrt(ou_fb_ll, ou_fa_ll, df = 2)
  lrt_fb_bm <- lrt(ou_fb_ll, bm_ll,    df = 5)
  
  # ── Effect sizes (from OU_fboth)
  pt_shift_abs <- ou_fb_theta_pt  - ou_fb_theta_pb
  ldy_shift_abs <- ou_fb_theta_ldy - ou_fb_theta_pb
  pt_progress <- if_else(abs(ldy_shift_abs) > 0.1,
                           pt_shift_abs / ldy_shift_abs,
                           NA_real_)
  
  t_half_pb_gen <- log(2) / ou_fb_alpha_pb
  t_half_pt_gen <- log(2) / ou_fb_alpha_pt
  t_half_ldy_gen <- log(2) / ou_fb_alpha_ldy
  t_half_pb_yr <- t_half_pb_gen  * generation_time
  t_half_pt_yr <- t_half_pt_gen  * generation_time
  t_half_ldy_yr <- t_half_ldy_gen * generation_time
  
  # ── Return
  tibble(
    gene_id = gid,
    converged = TRUE,
    n_pb = length(x_pb),
    n_pt = length(x_pt),
    n_ldy = length(x_ldy),
    n_obs = n_obs,
    # BM
    bm_sigma2 = bm_sigma2,
    bm_ancestral = bm_ancestral,
    bm_vw_pb = bm_vw_pb,
    bm_vw_pt = bm_vw_pt,
    bm_vw_ldy = bm_vw_ldy,
    bm_ll = bm_ll,
    bm_aic = aic_vec[["BM"]],
    bm_bic = bic_vec[["BM"]],
    # OU shared
    ou_s_alpha = ou_s_alpha,
    ou_s_theta = ou_s_theta,
    ou_s_sigma2 = ou_s_sigma2,
    ou_s_vw_pb = ou_s_vw_pb,
    ou_s_vw_pt = ou_s_vw_pt,
    ou_s_vw_ldy = ou_s_vw_ldy,
    ou_s_ll = ou_s_ll,
    ou_s_aic = aic_vec[["OU_shared"]],
    ou_s_bic = bic_vec[["OU_shared"]],
    # OU free theta
    ou_ft_alpha = ou_ft_alpha,
    ou_ft_theta_pb = ou_ft_theta_pb,
    ou_ft_theta_pt = ou_ft_theta_pt,
    ou_ft_theta_ldy = ou_ft_theta_ldy,
    ou_ft_sigma2 = ou_ft_sigma2,
    ou_ft_vw_pb = ou_ft_vw_pb,
    ou_ft_vw_pt = ou_ft_vw_pt,
    ou_ft_vw_ldy = ou_ft_vw_ldy,
    ou_ft_ll = ou_ft_ll,
    ou_ft_aic = aic_vec[["OU_ftheta"]],
    ou_ft_bic = bic_vec[["OU_ftheta"]],
    # OU free alpha
    ou_fa_alpha_pb = ou_fa_alpha_pb,
    ou_fa_alpha_pt = ou_fa_alpha_pt,
    ou_fa_alpha_ldy = ou_fa_alpha_ldy,
    ou_fa_theta = ou_fa_theta,
    ou_fa_sigma2 = ou_fa_sigma2,
    ou_fa_vw_pb = ou_fa_vw_pb,
    ou_fa_vw_pt = ou_fa_vw_pt,
    ou_fa_vw_ldy = ou_fa_vw_ldy,
    ou_fa_ll = ou_fa_ll,
    ou_fa_aic = aic_vec[["OU_falpha"]],
    ou_fa_bic = bic_vec[["OU_falpha"]],
    # OU free both
    ou_fb_alpha_pb = ou_fb_alpha_pb,
    ou_fb_alpha_pt = ou_fb_alpha_pt,
    ou_fb_alpha_ldy = ou_fb_alpha_ldy,
    ou_fb_theta_pb = ou_fb_theta_pb,
    ou_fb_theta_pt = ou_fb_theta_pt,
    ou_fb_theta_ldy = ou_fb_theta_ldy,
    ou_fb_sigma2 = ou_fb_sigma2,
    ou_fb_vw_pb = ou_fb_vw_pb,
    ou_fb_vw_pt = ou_fb_vw_pt,
    ou_fb_vw_ldy = ou_fb_vw_ldy,
    ou_fb_ll = ou_fb_ll,
    ou_fb_aic = aic_vec[["OU_fboth"]],
    ou_fb_bic = bic_vec[["OU_fboth"]],
    # Model selection
    best_model_aic = best_model_aic,
    best_model_bic = best_model_bic,
    # Delta AIC
    daic_ous_vs_bm = aic_vec[["BM"]] - aic_vec[["OU_shared"]],
    daic_ouft_vs_ous = aic_vec[["OU_shared"]] - aic_vec[["OU_ftheta"]],
    daic_oufa_vs_ous = aic_vec[["OU_shared"]] - aic_vec[["OU_falpha"]],
    daic_oufb_vs_bm = aic_vec[["BM"]] - aic_vec[["OU_fboth"]],
    # LRTs
    lrt_s_bm_stat = lrt_s_bm[["stat"]],
    lrt_s_bm_pval = lrt_s_bm[["pval"]],
    lrt_ft_s_stat = lrt_ft_s[["stat"]],
    lrt_ft_s_pval = lrt_ft_s[["pval"]],
    lrt_fa_s_stat = lrt_fa_s[["stat"]],
    lrt_fa_s_pval = lrt_fa_s[["pval"]],
    lrt_fb_ft_stat = lrt_fb_ft[["stat"]],
    lrt_fb_ft_pval = lrt_fb_ft[["pval"]],
    lrt_fb_fa_stat = lrt_fb_fa[["stat"]],
    lrt_fb_fa_pval = lrt_fb_fa[["pval"]],
    lrt_fb_bm_stat = lrt_fb_bm[["stat"]],
    lrt_fb_bm_pval = lrt_fb_bm[["pval"]],
    # Effect sizes
    pt_shift_abs = pt_shift_abs,
    ldy_shift_abs = ldy_shift_abs,
    pt_progress = pt_progress,
    t_half_pb_gen = t_half_pb_gen,
    t_half_pt_gen = t_half_pt_gen,
    t_half_ldy_gen = t_half_ldy_gen,
    t_half_pb_yr = t_half_pb_yr,
    t_half_pt_yr = t_half_pt_yr,
    t_half_ldy_yr = t_half_ldy_yr,
    # Observed means
    obs_mean_pb = obs_means[["PB"]],
    obs_mean_pt = obs_means[["PT"]],
    obs_mean_ldy = obs_means[["LDY"]]
  )
}
