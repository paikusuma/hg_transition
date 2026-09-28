library(tidyverse)
library(edgeR)

source("func_OU.R")

# ── MSMC-IM pairwise split times
generation_time <- 29

split_pb_pt <- 9195
split_pb_ldy <- 10677
split_pt_ldy <- 9037

# ── Root age = mean of all pairwise splits
T_root_years <- mean(c(split_pb_pt, split_pb_ldy, split_pt_ldy))
T_root_gen <- T_root_years / generation_time

cat("══ Star trifurcation tree ══════════════════════════════════\n")
cat(sprintf("  Root age (mean of MSMC-IM splits) : %.1f ya / %.1f gen\n",
            T_root_years, T_root_gen))
cat("  Topology:\n")
cat("        ┌──────────────── PB  (ancestral HG reference)\n")
cat("  ROOT ─┤──────────────── PT  (transitioning — focal population)\n")
cat("        └──────────────── LDY (agricultural endpoint)\n")

tree_params <- list(
  T_root = T_root_gen,
  T_root_years = T_root_years
)

# ── Load counts
dge <- DGEList(counts = counts_raw, group = meta$population)
cpm_dge <- cpm(dge)
cpm_PB <- cpm_dge[,grep("BLG",colnames(cpm_dge))]
cpm_PT <- cpm_dge[,grep("RPN",colnames(cpm_dge))]
cpm_LDY <- cpm_dge[,grep("SPI",colnames(cpm_dge))]

# set keep threshold to only retain genes that are expressed 
# at or over a CPM of one in at least half of the group
keep.PB <- rowSums(cpm_PB>1) >= (length(grep("BLG",rownames(dge$samples)))*0.5)
keep.PT <- rowSums(cpm_PT>1) >= (length(grep("RPN",rownames(dge$samples)))*0.5)
keep.LDY <- rowSums(cpm_LDY>1) >= (length(grep("SPI",rownames(dge$samples)))*0.5)

keep <- keep.PB | keep.PT | keep.LDY
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat(sprintf("Genes after filter: %d\n", nrow(dge)))

# ── Visualise batch x population balance
meta |>
  count(population, batch) |>
  pivot_wider(names_from = batch,
              values_from = n,
              values_fill = 0) |>
  print()

# Visual
meta |>
  ggplot(aes(x = batch, fill = population)) +
  geom_bar(position = "stack") +
  scale_fill_manual(
    values = c(PB = "steelblue", PT = "goldenrod3", LDY = "firebrick")
  ) +
  labs(title = "Batch × Population balance",
       subtitle = "Unbalanced = limma cannot cleanly separate batch from population",
       x = "Batch", y = "Number of samples") +
  theme_bw()


# ── Check cell type proportion collinearity
cell_cors <- meta |>
  select(Gran, CD19B, CD4T, CD8T, NK, CD14Mono) |>
  cor() |>
  round(2)
print(cell_cors)

meta |>
  select(Gran, CD19B, CD4T, CD8T, NK, CD14Mono) |>
  rowSums() |>
  summary()

# ── Build design matrix
design <- model.matrix(
  ~ 0 + population + batch + Age + sex + RIN + rRNA_ratio +
    Gran + CD19B + CD4T + CD8T + NK + CD14Mono,
  data = meta
)

# Check rank
cat(sprintf("Design matrix: %d rows × %d columns\n",
            nrow(design), ncol(design)))
cat(sprintf("Rank: %d  (full rank = %d)\n",
            qr(design)$rank, ncol(design)))
cat("Columns:\n")
print(colnames(design))

# ── TMM normalisation
dge <- calcNormFactors(dge, method = "TMM")

cat("Normalisation factors:\n")
print(summary(dge$samples$norm.factors))

# ── Voom transformation
png("voom_meanvariance.png", width = 600, height = 500)
v <- voom(dge, design, plot = TRUE)
dev.off()

cat(sprintf("Voom-transformed matrix: %d genes × %d samples\n",
            nrow(v$E), ncol(v$E)))

# ── Fit linear model
fit <- lmFit(v, design)

# ── Empirical Bayes moderation
fit <- eBayes(fit, robust = TRUE)

# ── Extract adjusted means directly
adj_mean_pb <- fit$coefficients[, "populationPB"]
adj_mean_pt <- fit$coefficients[, "populationPT"]
adj_mean_ldy <- fit$coefficients[, "populationLDY"]

# ── Standard errors for each population mean
sigma_post <- sqrt(fit$s2.post)

se_pb <- sigma_post * fit$stdev.unscaled[, "populationPB"]
se_pt <- sigma_post * fit$stdev.unscaled[, "populationPT"]
se_ldy <- sigma_post * fit$stdev.unscaled[, "populationLDY"]

# ── Coefficients for effect size (differences from PB)
coef_pt <- adj_mean_pt  - adj_mean_pb
coef_ldy <- adj_mean_ldy - adj_mean_pb

# ── Compile limma_adjusted
limma_adjusted <- tibble(
  gene_id = rownames(fit$coefficients),
  adj_mean_pb = adj_mean_pb,
  adj_mean_pt = adj_mean_pt,
  adj_mean_ldy = adj_mean_ldy,
  coef_pt = coef_pt,
  coef_ldy = coef_ldy,
  se_pb = se_pb,
  se_pt = se_pt,
  se_ldy = se_ldy,
  se2_pb = se_pb^2,
  se2_pt = se_pt^2,
  se2_ldy = se_ldy^2,
  # Moderated t and p for PT vs PB and LDY vs PB
  # With 0 + population, need contrasts for these
  sigma_post = sigma_post
)

cat(sprintf("limma_adjusted ready: %d genes\n", nrow(limma_adjusted)))

# ── Define contrasts
contrast_mat <- makeContrasts(
  PT_vs_PB = populationPT  - populationPB,
  LDY_vs_PB = populationLDY - populationPB,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2, robust = TRUE)

# ── Add DE results to limma_adjusted
limma_adjusted <- limma_adjusted |>
  mutate(
    t_pt = fit2$t[, "PT_vs_PB"],
    t_ldy = fit2$t[, "LDY_vs_PB"],
    pval_pt = fit2$p.value[, "PT_vs_PB"],
    pval_ldy = fit2$p.value[, "LDY_vs_PB"],
    fdr_pt = p.adjust(fit2$p.value[, "PT_vs_PB"],  method = "BH"),
    fdr_ldy = p.adjust(fit2$p.value[, "LDY_vs_PB"], method = "BH")
  )

limma_adjusted <- limma_adjusted |>
  mutate(
    mean_range = pmax(
      abs(adj_mean_pb - adj_mean_pt),
      abs(adj_mean_pb - adj_mean_ldy),
      abs(adj_mean_pt - adj_mean_ldy)
    ),
    max_se = pmax(se_pb, se_pt, se_ldy),
    se_ok = max_se < 10 & se_pb > 0 & se_pt > 0 & se_ldy > 0
  )

cat(sprintf("Genes passing SE filter: %d / %d\n",
            sum(limma_adjusted$se_ok),
            nrow(limma_adjusted)))

limma_adjusted <- limma_adjusted |> dplyr::filter(se_ok)

cat(sprintf("PT vs PB  FDR < 0.05: %d genes\n",
            sum(limma_adjusted$fdr_pt < 0.05)))
cat(sprintf("LDY vs PB FDR < 0.05: %d genes\n",
            sum(limma_adjusted$fdr_ldy < 0.05)))

# ── Residualise voom expression
nuisance_cols <- colnames(design)[!grepl("^population", colnames(design))]
cat("Nuisance columns being removed:\n")
print(nuisance_cols)

# Adjusted expression = observed - nuisance fitted values
# This preserves population signal + residual noise
nuisance_fit <- fit$coefficients[, nuisance_cols, drop = FALSE] %*%
  t(design[, nuisance_cols, drop = FALSE])

# v$E is genes x samples (voom log-CPM)
expr_adj <- v$E - nuisance_fit

cat(sprintf("expr_adj dimensions: %d genes × %d samples\n",
            nrow(expr_adj), ncol(expr_adj)))

# Verify: per-population means of expr_adj should match limma_adjusted means
pop_check <- tibble(
  gene_id = rownames(expr_adj),
  mean_pb_adj = rowMeans(expr_adj[, meta$population == "PB"]),
  mean_pt_adj = rowMeans(expr_adj[, meta$population == "PT"]),
  mean_ldy_adj= rowMeans(expr_adj[, meta$population == "LDY"])
) |>
  left_join(limma_adjusted |>
              select(gene_id, adj_mean_pb, adj_mean_pt, adj_mean_ldy),
            by = "gene_id") |>
  mutate(
    diff_pb = mean_pb_adj  - adj_mean_pb,
    diff_pt = mean_pt_adj  - adj_mean_pt,
    diff_ldy = mean_ldy_adj - adj_mean_ldy
  )

cat("\nMean difference between expr_adj means and limma coefficients:\n")
cat(sprintf("  PB : %.4f\n", mean(pop_check$diff_pb,  na.rm = TRUE)))
cat(sprintf("  PT : %.4f\n", mean(pop_check$diff_pt,  na.rm = TRUE)))
cat(sprintf("  LDY: %.4f\n", mean(pop_check$diff_ldy, na.rm = TRUE)))


# ── Stratified test sample across effect sizes
# Use limma coefficients to define effect size for stratification

effect_sizes_limma <- limma_adjusted |>
  mutate(
    max_diff = pmax(abs(coef_pt), abs(coef_ldy))
  )

set.seed(42)
test_genes_limma <- bind_rows(
  # Large effect
  effect_sizes_limma |>
    slice_max(max_diff, n = 10) |>
    mutate(effect_group = "large"),
  # Medium effect
  effect_sizes_limma |>
    filter(max_diff > quantile(max_diff, 0.6),
           max_diff < quantile(max_diff, 0.7)) |>
    slice_sample(n = 10) |>
    mutate(effect_group = "medium"),
  # Small effect
  effect_sizes_limma |>
    filter(max_diff > quantile(max_diff, 0.4),
           max_diff < quantile(max_diff, 0.5)) |>
    slice_sample(n = 10) |>
    mutate(effect_group = "small"),
  # Tiny effect
  effect_sizes_limma |>
    slice_min(max_diff, n = 10) |>
    mutate(effect_group = "tiny")
) |>
  distinct(gene_id, .keep_all = TRUE)

cat(sprintf("Test genes: %d\n", nrow(test_genes_limma)))

# ── Run test without parallelisation
test_results_limma <- map_dfr(
  test_genes_limma$gene_id,
  function(g) {
    cat(sprintf("Fitting: %s\n", g))
    fit_gene_limma(
      g,
      limma_adjusted = limma_adjusted,
      tree = tree_params,
      n_restarts = 15
    )
  }
)

# ── Diagnostics
cat("\n── Convergence ──\n")
cat(sprintf("  Converged : %d / %d\n",
            sum(test_results_limma$converged, na.rm = TRUE),
            nrow(test_results_limma)))

cat("\n── Alpha (should NOT be pinned at 403) ──\n")
summary(test_results_limma$ou_f_alpha)

cat("\n── ou_f_ll (should VARY across genes) ──\n")
summary(test_results_limma$ou_f_ll)

cat("\n── pt_progress ──\n")
summary(test_results_limma$pt_progress)

cat("\n── Model selection by effect group ──\n")
test_results_limma |>
  left_join(
    test_genes_limma |> select(gene_id, effect_group),
    by = "gene_id"
  ) |>
  count(effect_group, best_model) |>
  pivot_wider(names_from = best_model,
              values_from = n,
              values_fill = 0) |>
  print()

test_results_limma |>
  left_join(test_genes_limma, by = "gene_id") |>
  mutate(
    mean_range = pmax(
      abs(obs_mean_pb  - obs_mean_pt),
      abs(obs_mean_pb  - obs_mean_ldy),
      abs(obs_mean_pt  - obs_mean_ldy)
    ),
    snr = mean_range / sqrt((se2_pb + se2_pt + se2_ldy) / 3)
  ) |>
  select(gene_id, effect_group, mean_range, snr,
         se2_pb, ou_f_alpha, best_model) |>
  arrange(desc(snr)) |>
  print(n = 40)


limma_adjusted_nose <- limma_adjusted |>
  mutate(
    se2_pb = 1e-6,
    se2_pt = 1e-6,
    se2_ldy = 1e-6
  )

test_results_nose <- map_dfr(
  test_genes_limma$gene_id,
  function(g) {
    cat(sprintf("Fitting: %s\n", g))
    fit_gene_limma(
      g,
      limma_adjusted = limma_adjusted_nose,
      tree = tree_params,
      n_restarts = 15
    )
  }
)

test_results_nose |> count(best_model)
summary(test_results_nose$ou_f_alpha)

# ── Permutation calibration of LRT2
# Use a small-effect gene for the null distribution

set.seed(42)
n_perm <- 200

calibration_gene <- effect_sizes_limma |>
  filter(
    max_diff > quantile(max_diff, 0.4),
    max_diff < quantile(max_diff, 0.5)
  ) |>
  slice_sample(n = 1) |>
  pull(gene_id)

cat(sprintf("Calibration gene: %s\n", calibration_gene))

# Real fit
real_fit_limma <- fit_gene_limma(
  gene_id = calibration_gene,
  limma_adjusted = limma_adjusted,
  tree = tree_params,
  n_restarts = 15
)

cat(sprintf("Real LRT2 stat : %.3f\n", real_fit_limma$lrt2_stat))
cat(sprintf("Real LRT2 pval : %.4f\n", real_fit_limma$lrt2_pval))

# ── Permutation — shuffle population labels in limma_adjusted
# Shuffle which adjusted mean belongs to which population

gene_row <- limma_adjusted |>
  filter(gene_id == calibration_gene)

perm_lrt2 <- map_dbl(seq_len(n_perm), function(i) {
  
  # Shuffle the three adjusted means across population labels
  shuffled_means <- sample(c(
    gene_row$adj_mean_pb,
    gene_row$adj_mean_pt,
    gene_row$adj_mean_ldy
  ))
  
  # Shuffled SE2 values (keep paired with shuffled means)
  shuffled_se2 <- sample(c(
    gene_row$se2_pb,
    gene_row$se2_pt,
    gene_row$se2_ldy
  ))
  
  # Build permuted limma_adjusted row
  perm_row <- gene_row |>
    mutate(
      adj_mean_pb = shuffled_means[1],
      adj_mean_pt = shuffled_means[2],
      adj_mean_ldy = shuffled_means[3],
      se2_pb = shuffled_se2[1],
      se2_pt = shuffled_se2[2],
      se2_ldy = shuffled_se2[3]
    )
  
  # Refit with permuted data
  perm_limma <- bind_rows(
    perm_row,
    limma_adjusted |> filter(gene_id != calibration_gene)
  )
  
  fit <- tryCatch(
    fit_gene_limma(
      gene_id = calibration_gene,
      limma_adjusted = perm_limma,
      tree = tree_params,
      n_restarts = 5
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !fit$converged) return(NA_real_)
  return(fit$lrt2_stat)
}) |>
  na.omit()

# ── Report calibration
cat(sprintf("\nPermutation LRT2 distribution (n=%d):\n",
            length(perm_lrt2)))
print(summary(perm_lrt2))

cat(sprintf("\nTheoretical chi-sq(df=2) 95th pct : %.3f\n",
            qchisq(0.95, df = 2)))
cat(sprintf("Permutation 95th percentile        : %.3f\n",
            quantile(perm_lrt2, 0.95, na.rm = TRUE)))
cat(sprintf("Inflation factor                   : %.3f\n",
            quantile(perm_lrt2, 0.95, na.rm = TRUE) /
              qchisq(0.95, df = 2)))

# —— FULL RUN

library(furrr)
plan(multisession, workers = 30)

full_results <- future_map_dfr(
  limma_adjusted$gene_id,
  function(g) {
    fit_gene_limma(
      g,
      limma_adjusted = limma_adjusted,
      expr_adj = expr_adj,
      meta = meta,
      tree = tree_params,
      n_restarts = 15
    )
  },
  .progress = TRUE
)

write_csv(full_results, "ou_model_results_full_allmodels.csv")

cat("══ Convergence ══════════════════════════════════════════\n")
cat(sprintf("  Total genes     : %d\n", nrow(full_results)))
cat(sprintf("  Converged       : %d (%.1f%%)\n",
            sum(full_results$converged),
            100 * mean(full_results$converged)))
cat(sprintf("  Failed          : %d (%.1f%%)\n",
            sum(!full_results$converged),
            100 * mean(!full_results$converged)))

cat("\n══ Model selection (AIC) ════════════════════════════════\n")
full_results |>
  filter(converged) |>
  count(best_model_aic) |>
  mutate(pct = 100 * n / sum(n)) |>
  arrange(desc(n)) |>
  print()


cat("\n══ Alpha distributions (converged genes) ════════════════\n")

full_results |>
  filter(converged) |>
  select(ou_s_alpha, ou_ft_alpha,
         ou_fb_alpha_pb, ou_fb_alpha_pt, ou_fb_alpha_ldy) |>
  summary() |>
  print()

cat("\n══ Half-lives (OU_fboth, generations) ═══════════════════\n")
full_results |>
  filter(converged) |>
  select(t_half_pb_gen, t_half_pt_gen, t_half_ldy_gen) |>
  summary() |>
  print()

cat("\n══ Half-lives (OU_fboth, years) ═════════════════════════\n")
full_results |>
  filter(converged) |>
  select(t_half_pb_yr, t_half_pt_yr, t_half_ldy_yr) |>
  summary() |>
  print()

# Check: are any alpha values still pinned at bounds?
T_root_gen <- mean(c(9195, 10677, 9037)) / 29
log_alpha_lo <- log(log(2) / (T_root_gen * 100))
log_alpha_hi <- log(log(2) / (T_root_gen * 0.01))
alpha_lo <- exp(log_alpha_lo)
alpha_hi <- exp(log_alpha_hi)

cat(sprintf("\n  Alpha bounds: [%.6f, %.4f]\n", alpha_lo, alpha_hi))

n_pinned_hi <- full_results |>
  filter(converged) |>
  summarise(
    pb = sum(abs(ou_fb_alpha_pb  - alpha_hi) < 1e-4),
    pt = sum(abs(ou_fb_alpha_pt  - alpha_hi) < 1e-4),
    ldy = sum(abs(ou_fb_alpha_ldy - alpha_hi) < 1e-4)
  )

n_pinned_lo <- full_results |>
  filter(converged) |>
  summarise(
    pb = sum(abs(ou_fb_alpha_pb  - alpha_lo) < 1e-8),
    pt = sum(abs(ou_fb_alpha_pt  - alpha_lo) < 1e-8),
    ldy = sum(abs(ou_fb_alpha_ldy - alpha_lo) < 1e-8)
  )

cat("\n  Genes pinned at UPPER alpha bound:\n")
print(n_pinned_hi)
cat("\n  Genes pinned at LOWER alpha bound:\n")
print(n_pinned_lo)


cat("\n══ Log-likelihood distributions ════════════════════════\n")
full_results |>
  filter(converged) |>
  select(bm_ll, ou_s_ll, ou_ft_ll, ou_fa_ll, ou_fb_ll) |>
  summary() |>
  print()

cat("\n══ Delta AIC (positive = complex model better) ══════════\n")
full_results |>
  filter(converged) |>
  select(daic_ous_vs_bm, daic_ouft_vs_ous,
         daic_oufa_vs_ous, daic_oufb_vs_bm) |>
  summary() |>
  print()

# Sanity: OU_fboth should NEVER have lower ll than simpler models
# (it is the most complex — ll should be >= all others)
bad_ll <- full_results |>
  filter(converged) |>
  filter(ou_fb_ll < ou_ft_ll - 0.01 |
           ou_fb_ll < ou_fa_ll - 0.01 |
           ou_fb_ll < ou_s_ll  - 0.01 |
           ou_fb_ll < bm_ll    - 0.01)

cat(sprintf("\n  Genes where OU_fboth ll < simpler model: %d\n",
            nrow(bad_ll)))
cat("  (should be 0 — any non-zero indicates optimiser failure)\n")


cat("\n══ Effect sizes (from OU_ftheta — most common winner) ═══\n")
full_results |>
  filter(converged, best_model_aic == "OU_ftheta") |>
  select(pt_shift_abs, ldy_shift_abs, pt_progress) |>
  summary() |>
  print()

cat("\n══ PT progress distribution ═════════════════════════════\n")
cat("  (0 = same as PB, 1 = same as LDY, intermediate = transitioning)\n")
full_results |>
  filter(converged, !is.na(pt_progress)) |>
  summarise(
    n = n(),
    below_0 = sum(pt_progress < 0),
    btw_0_05 = sum(pt_progress >= 0   & pt_progress < 0.5),
    btw_05_1 = sum(pt_progress >= 0.5 & pt_progress < 1),
    above_1 = sum(pt_progress >= 1),
    median_ptp = median(pt_progress),
    mean_ptp = mean(pt_progress)
  ) |>
  print()

cat("\n══ LRT significance summary ═════════════════════════════\n")

full_results |>
  filter(converged) |>
  summarise(
    # OU signal at all (shared OU vs BM)
    n_ou_vs_bm_p05 = sum(lrt_s_bm_pval < 0.05),
    n_ou_vs_bm_p01 = sum(lrt_s_bm_pval < 0.01),
    # Different optima (free theta vs shared)
    n_ftheta_p05 = sum(lrt_ft_s_pval < 0.05),
    n_ftheta_p01 = sum(lrt_ft_s_pval < 0.01),
    # Different constraint (free alpha vs shared)
    n_falpha_p05 = sum(lrt_fa_s_pval < 0.05),
    n_falpha_p01 = sum(lrt_fa_s_pval < 0.01),
    # Full model vs BM
    n_fboth_vs_bm_p05= sum(lrt_fb_bm_pval < 0.05),
    n_fboth_vs_bm_p01= sum(lrt_fb_bm_pval < 0.01)
  ) |>
  pivot_longer(everything(),
               names_to = "test",
               values_to = "n_significant") |>
  mutate(pct = 100 * n_significant / sum(full_results$converged)) |>
  print()

# FDR correction on the key tests
full_results <- full_results |>
  mutate(
    fdr_ou_vs_bm = p.adjust(lrt_s_bm_pval,  method = "BH"),
    fdr_ftheta = p.adjust(lrt_ft_s_pval,   method = "BH"),
    fdr_falpha = p.adjust(lrt_fa_s_pval,   method = "BH"),
    fdr_fboth_bm = p.adjust(lrt_fb_bm_pval,  method = "BH")
  )

cat("\n══ FDR < 0.05 summary ═══════════════════════════════════\n")
full_results |>
  filter(converged) |>
  summarise(
    n_ou_vs_bm = sum(fdr_ou_vs_bm < 0.05, na.rm = TRUE),
    n_ftheta = sum(fdr_ftheta < 0.05, na.rm = TRUE),
    n_falpha = sum(fdr_falpha < 0.05, na.rm = TRUE),
    n_fboth_bm = sum(fdr_fboth_bm < 0.05, na.rm = TRUE)
  ) |>
  pivot_longer(everything(),
               names_to = "test",
               values_to = "n_fdr_significant") |>
  mutate(pct = 100 * n_fdr_significant / sum(full_results$converged)) |>
  print()


cat("\n══ Top 10 genes by PT shift (most changed from PB) ══════\n")
full_results |>
  filter(converged, best_model_aic %in% c("OU_ftheta", "OU_fboth")) |>
  arrange(desc(abs(pt_shift_abs))) |>
  select(gene_id, best_model_aic,
         ou_ft_theta_pb, ou_ft_theta_pt, ou_ft_theta_ldy,
         pt_shift_abs, ldy_shift_abs, pt_progress,
         fdr_ftheta) |>
  head(10) |>
  print() |> 
  pull(gene_id)

cat("\n══ Top 10 genes by LDY shift (most changed from PB) ════\n")
full_results |>
  filter(converged, best_model_aic %in% c("OU_ftheta", "OU_fboth")) |>
  arrange(desc(abs(ldy_shift_abs))) |>
  select(gene_id, best_model_aic,
         ou_ft_theta_pb, ou_ft_theta_pt, ou_ft_theta_ldy,
         pt_shift_abs, ldy_shift_abs, pt_progress,
         fdr_ftheta) |>
  head(10) |>
  print()


# Is daic_oufa_vs_ous really exactly -4 for every gene?
full_results |>
  filter(converged) |>
  pull(daic_oufa_vs_ous) |>
  table() |>
  head(20)

# Or is it approximately -4 with some variation?
full_results |>
  filter(converged) |>
  summarise(
    min = min(daic_oufa_vs_ous),
    p25 = quantile(daic_oufa_vs_ous, 0.25),
    p50 = median(daic_oufa_vs_ous),
    p75 = quantile(daic_oufa_vs_ous, 0.75),
    max = max(daic_oufa_vs_ous),
    n_gt0 = sum(daic_oufa_vs_ous > 0)   # how many genes where falpha beats shared?
  ) |>
  print()

# Check the actual likelihood gain
full_results |>
  filter(converged) |>
  mutate(ll_gain_fa = ou_fa_ll - ou_s_ll) |>
  summarise(
    min_gain = min(ll_gain_fa),
    median_gain = median(ll_gain_fa),
    max_gain = max(ll_gain_fa),
    n_gain_gt2 = sum(ll_gain_fa > 2),   # would overcome AIC penalty
    n_gain_gt1 = sum(ll_gain_fa > 1)
  ) |>
  print()


# Where are the NAs coming from?
full_results |>
  filter(converged) |>
  mutate(
    na_reason = case_when(
      is.na(pt_progress) & abs(ldy_shift_abs) <= 0.3 ~ "ldy_diff_too_small_0.3",
      is.na(pt_progress) & abs(ldy_shift_abs) <= 0.1 ~ "ldy_diff_too_small_0.1",
      is.na(pt_progress) ~ "other_na",
      TRUE ~ "has_value"
    )
  ) |>
  count(na_reason) |>
  print()

# Distribution of ldy_shift_abs across ALL converged genes
full_results |>
  filter(converged) |>
  summarise(
    n_total = n(),
    n_ldy_gt_03 = sum(abs(ldy_shift_abs) > 0.3, na.rm = TRUE),
    n_ldy_gt_01 = sum(abs(ldy_shift_abs) > 0.1, na.rm = TRUE),
    n_ldy_gt_005 = sum(abs(ldy_shift_abs) > 0.05, na.rm = TRUE),
    median_ldy_shift = median(abs(ldy_shift_abs), na.rm = TRUE)
  ) |>
  print()


# Use thetas from the winning model rather than always from OU_fboth
full_results <- full_results |>
  mutate(
    # Use OU_ftheta thetas when that model wins (most informative)
    theta_pb_best = case_when(
      best_model_aic == "OU_ftheta" ~ ou_ft_theta_pb,
      best_model_aic == "OU_fboth" ~ ou_fb_theta_pb,
      best_model_aic == "OU_shared" ~ ou_s_theta,
      TRUE ~ obs_mean_pb      # BM: use observed mean
    ),
    theta_pt_best = case_when(
      best_model_aic == "OU_ftheta" ~ ou_ft_theta_pt,
      best_model_aic == "OU_fboth" ~ ou_fb_theta_pt,
      best_model_aic == "OU_shared" ~ ou_s_theta,
      TRUE ~ obs_mean_pt
    ),
    theta_ldy_best = case_when(
      best_model_aic == "OU_ftheta" ~ ou_ft_theta_ldy,
      best_model_aic == "OU_fboth" ~ ou_fb_theta_ldy,
      best_model_aic == "OU_shared" ~ ou_s_theta,
      TRUE ~ obs_mean_ldy
    ),
    pt_shift_best = theta_pt_best  - theta_pb_best,
    ldy_shift_best = theta_ldy_best - theta_pb_best,
    pt_progress_best = if_else(
      abs(ldy_shift_best) > 0.1,
      pt_shift_best / ldy_shift_best,
      NA_real_
    )
  )

# How many genes now have pt_progress?
full_results |>
  filter(converged) |>
  summarise(
    n_with_progress = sum(!is.na(pt_progress_best)),
    pct = 100 * mean(!is.na(pt_progress_best)),
    median_progress = median(pt_progress_best, na.rm = TRUE),
    mean_progress = mean(pt_progress_best,   na.rm = TRUE)
  ) |>
  print()


# ── Add all final annotations
full_results <- full_results |>
  mutate(
    # Flag identifiability issues
    alpha_pinned_hi = abs(ou_ft_alpha - alpha_hi) < 1e-4,
    
    # Simplified model label (drop unidentifiable models)
    best_model_final = case_when(
      best_model_aic == "OU_ftheta" ~ "OU_ftheta",
      best_model_aic == "OU_shared" ~ "OU_shared",
      best_model_aic == "OU_falpha" ~ "OU_ftheta",  # reassign — unidentifiable
      best_model_aic == "OU_fboth" ~ "OU_ftheta",  # reassign — unidentifiable
      TRUE ~ "BM"
    ),
    
    # Direction of PT shift relative to LDY
    pt_direction = case_when(
      is.na(pt_progress_best) ~ "undefined",
      pt_progress_best < 0 ~ "opposite",
      pt_progress_best < 0.5 ~ "partial_early",
      pt_progress_best < 1 ~ "partial_late",
      pt_progress_best >= 1 ~ "overshoot"
    ),
    
    # Key significance flags
    sig_ou = fdr_ftheta < 0.05 & best_model_final == "OU_ftheta",
    sig_ou_bm = fdr_fboth_bm < 0.05
  )

# ── Summary of final gene sets ────────────────────────────────────────────────
cat("══ Final gene set summary ═══════════════════════════════════\n")
cat(sprintf("  Total converged              : %d\n",
            sum(full_results$converged)))
cat(sprintf("  Best model = OU_ftheta       : %d (%.1f%%)\n",
            sum(full_results$best_model_final == "OU_ftheta"),
            100 * mean(full_results$best_model_final == "OU_ftheta")))
cat(sprintf("  Best model = BM              : %d (%.1f%%)\n",
            sum(full_results$best_model_final == "BM"),
            100 * mean(full_results$best_model_final == "BM")))
cat(sprintf("  Significant OU (FDR<0.05)    : %d (%.1f%%)\n",
            sum(full_results$sig_ou, na.rm = TRUE),
            100 * mean(full_results$sig_ou, na.rm = TRUE)))
cat(sprintf("  Alpha pinned at upper bound  : %d (%.1f%%)\n",
            sum(full_results$alpha_pinned_hi, na.rm = TRUE),
            100 * mean(full_results$alpha_pinned_hi, na.rm = TRUE)))

cat("\n══ PT direction (genes with pt_progress_best) ═══════════════\n")
full_results |>
  filter(!is.na(pt_progress_best)) |>
  count(pt_direction) |>
  mutate(pct = 100 * n / sum(n)) |>
  arrange(desc(n)) |>
  print()

cat("\n══ Key gene sets ════════════════════════════════════════════\n")

# Set 1: Strong OU signal, PT transitioning toward LDY
genes_transitioning <- full_results |>
  filter(
    sig_ou,
    !is.na(pt_progress_best),
    pt_progress_best > 0,
    pt_progress_best < 1
  )
cat(sprintf("  PT transitioning (0 < progress < 1, FDR<0.05): %d genes\n",
            nrow(genes_transitioning)))

# Set 2: PT overshooting LDY
genes_overshoot <- full_results |>
  filter(sig_ou, pt_progress_best >= 1)
cat(sprintf("  PT overshooting LDY (progress >= 1, FDR<0.05): %d genes\n",
            nrow(genes_overshoot)))

# Set 3: PT moving opposite to LDY
genes_opposite <- full_results |>
  filter(sig_ou, pt_progress_best < 0)
cat(sprintf("  PT opposite to LDY (progress < 0, FDR<0.05)  : %d genes\n",
            nrow(genes_opposite)))

# Set 4: Strong BM genes (no directional selection)
genes_bm <- full_results |>
  filter(best_model_final == "BM")
cat(sprintf("  BM genes (drift only)                        : %d genes\n",
            nrow(genes_bm)))

# ── Save final annotated results
write_csv(full_results, "ou_model_results_full_final.csv")
