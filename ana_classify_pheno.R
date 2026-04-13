# =================================================================
# Phenotype Classification: Lifestyle Group Prediction (HG vs AGR)
# Methods: Univariate tests, PCA, LDA, QDA, Random Forest
# =================================================================

library(here)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(paletteer)
library(car)
library(rstatix)
library(glue)
library(MASS)
library(caret)
library(randomForest)
library(pROC)


# =================
# 1. Prepare Data
# =================

# Recode three-group variable into binary lifestyle label (HG vs AGR)
dat_all <- dat_all %>%
  mutate(Group2 = factor(
    ifelse(Group == "PunanBatu", "HG", "AGR"),
    levels = c("HG", "AGR")
  ))

df <- dat_all

# Phenotypic variables to test
vars_to_test <- c(
  "Bio_BMI", "Bio_BMRI", "Bio_BodyFat", "Bio_VisceralFat",
  "Bio_SubcutFat_Overall", "Bio_MuscleMass_Overall",
  "Blood_Glucose", "Blood_Cholesterol", "Blood_HDL",
  "Blood_Triglyceride", "Blood_LDL", "Blood_UricAcid"
)


# ========================================
# 2. Univariate Tests (Pairwise Wilcoxon)
# ========================================

# Run pairwise Wilcoxon tests with Bonferroni correction for each variable
results_pairwise_wilcox <- map_dfr(vars_to_test, function(var) {
  formula <- as.formula(glue("{var} ~ Group2"))
  pairwise_wilcox_test(df, formula, p.adjust.method = "bonferroni")
})

# write_tsv(results_pairwise_wilcox, here("result_pairwise_pheno_HGvAGR.txt"))


# ==========================
# 3. Effect Size Estimation
# ==========================

# Candidate variables for effect size analysis
cands <- c(
  "Bio_BMI", "Bio_BMRI", "Bio_BodyFat", "Bio_VisceralFat",
  "Bio_SubcutFat_Overall", "Blood_Cholesterol",
  "Blood_Glucose", "Blood_LDL", "Blood_Triglyceride"
)

# Reshape to long format for grouped testing
df_long <- df %>%
  dplyr::select(Group2, all_of(cands)) %>%
  pivot_longer(-Group2, names_to = "variable", values_to = "value")

# Wilcoxon test with Bonferroni correction
wilcox_results <- df_long %>%
  group_by(variable) %>%
  wilcox_test(value ~ Group2) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni"))

# Effect sizes (rank-biserial correlation)
eff_results <- df_long %>%
  group_by(variable) %>%
  wilcox_effsize(value ~ Group2) %>%
  ungroup()

# Merge test results and effect sizes, sorted by effect size
results_combined <- wilcox_results %>%
  left_join(eff_results, by = "variable") %>%
  arrange(desc(effsize))

print(results_combined)
# write_tsv(results_combined, here("result_effectsize_HGvAGR.txt"))


# ===========================
# 4. Visualise Top Variables
# ===========================

# Select top 6 variables by effect size
top_vars <- results_combined %>%
  slice_max(effsize, n = 6) %>%
  pull(variable)

df_long2 <- df %>%
  pivot_longer(all_of(top_vars), names_to = "variable", values_to = "value")

pl_top6vars <- ggplot(df_long2, aes(x = Group2, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal()


# ========================
# 5. PCA on Top Variables
# ========================

# Scale top variables and run PCA
mat <- df %>%
  dplyr::select(all_of(top_vars)) %>%
  scale()

pca <- prcomp(mat, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x) %>%
  bind_cols(Group = df$Group)

# Variance explained per component
pca_var <- pca$sdev^2 / sum(pca$sdev^2) * 100

pl_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  scale_color_paletteer_d("ggsci::default_aaas") +
  labs(
    x = paste0("PC1 (", round(pca_var[1], 2), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 2), "%)")
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(3, "mm"),
    axis.title = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 8)
  )

# ggsave(here("pl_pca_top6vars_HGvAGR.pdf"), pl_pca,
#        width = 114, height = 78, units = "mm", dpi = 300)


# =======================================
# 6. Prepare Scaled Data for Classifiers
# =======================================

# Fixed set of top predictors for all classifiers
top_vars <- c(
  "Bio_VisceralFat", "Bio_BMI", "Bio_SubcutFat_Overall",
  "Bio_BodyFat", "Blood_Cholesterol", "Blood_Triglyceride"
)

# Scale predictors to mean 0, sd 1
df_lda <- df %>%
  dplyr::select(Group2, all_of(top_vars)) %>%
  mutate(across(all_of(top_vars), ~scale(.)[, 1])) %>%
  mutate(Group2 = factor(Group2, levels = c("HG", "AGR")))

# Shared 5-fold cross-validation control
set.seed(2025)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  savePredictions = "final"
)


# =====================================================
# 7. Linear Discriminant Analysis (LDA) with 5-Fold CV
# =====================================================

lda_cv <- train(
  Group2 ~ .,
  data = df_lda,
  method = "lda",
  trControl = ctrl,
  preProcess = NULL  # already scaled
)

print(lda_cv)

# Confusion matrix on out-of-fold predictions
cv_conf_lda <- confusionMatrix(lda_cv$pred$pred, lda_cv$pred$obs)
print(cv_conf_lda)


# ========================================================
# 8. Quadratic Discriminant Analysis (QDA) with 5-Fold CV
# ========================================================

set.seed(2025)
qda_cv <- train(
  Group2 ~ .,
  data = df_lda,
  method = "qda",
  trControl = ctrl
)

print(qda_cv)

# Confusion matrix on out-of-fold predictions
cv_conf_qda <- confusionMatrix(qda_cv$pred$pred, qda_cv$pred$obs)
print(cv_conf_qda)


# ================================
# 9. Random Forest with 5-Fold CV
# ================================

# Tune grid: explore mtry values around the default sqrt(p)
mtry_vals <- floor(sqrt(ncol(df_lda) - 1)) + (-1:1)
tuneGrid <- expand.grid(mtry = unique(pmax(1, mtry_vals)))

set.seed(2025)
rf_cv <- train(
  Group2 ~ .,
  data = df_lda,
  method = "rf",
  trControl = ctrl,
  tuneGrid = tuneGrid,
  importance = TRUE
)

print(rf_cv)

# Extract out-of-fold predictions for the best mtry
preds_rf <- rf_cv$pred %>%
  filter(mtry == rf_cv$bestTune$mtry)

# Confusion matrix
rf_conf <- confusionMatrix(preds_rf$pred, preds_rf$obs)
print(rf_conf)


# ========================================
# 10. Variable Importance (Random Forest)
# ========================================

var_imp <- varImp(rf_cv, scale = FALSE)
print(var_imp)

# Tidy variable importance data for plotting
var_imp_df <- var_imp$importance %>%
  rownames_to_column("Variable") %>%
  dplyr::rename(Importance = Overall) %>%
  arrange(Importance)

pl_varimp <- ggplot(var_imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Variable Importance",
    x = "Variable",
    y = "Mean Decrease in Accuracy"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 8),
    title = element_text(size = 10)
  )


# ==============================
# 11. ROC Curve (Random Forest)
# ==============================

# Extract out-of-fold predicted probabilities for the positive class (HG)
rf_probs <- rf_cv$pred %>%
  dplyr::filter(mtry == rf_cv$bestTune$mtry) %>%
  dplyr::select(obs, HG)

# Compute ROC curve (AGR = negative class, HG = positive class)
roc_obj <- roc(
  response = rf_probs$obs,
  predictor = rf_probs$HG,
  levels = c("AGR", "HG")
)
print(roc_obj)

# Tidy ROC data for ggplot
roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

pl_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = paste0("ROC Curve — Random Forest (AUC = ", round(auc(roc_obj), 3), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 10)
  )

