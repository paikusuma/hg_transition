# =====================================================
# Weighted Gene Co-expression Network Analysis (WGCNA)
# =====================================================

library(here)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(WGCNA)
library(CorLevelPlot)


# ====================================
# 1. Load and Prepare Expression Data
# ====================================

# Load TMM-normalised log-CPM matrix and transpose for WGCNA (samples x genes)
norm.counts <- read_tsv(here("out3_log_cpm_dge_filt_norm_ftcnts.txt")) %>%
  column_to_rownames("Geneid") %>%
  t()

# Ensure all values are numeric
norm.counts[] <- sapply(norm.counts, as.numeric)


# =====================================
# 2. Soft-Thresholding Power Selection
# =====================================

power <- c(1:10, seq(12, 50, 2))
sft <- pickSoftThreshold(norm.counts, powerVector = power,
                         networkType = "signed", verbose = 5)
sft.data <- sft$fitIndices

# Plot scale-free topology fit and mean connectivity vs power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  labs(x = "Power", y = "Scale-free topology model fit (signed R2)") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean connectivity") +
  theme_classic()

pl_sft_power <- ggarrange(a1, a2, nrow = 2, ncol = 1, align = "hv")

# ggsave(here("pl_sft_power.pdf"), pl_sft_power, width = 174, height = 174,
#        units = "mm", dpi = 300)


# ===========================================
# 3. Network Construction (blockwiseModules)
# ===========================================

soft_power <- 14

# Temporarily override cor() with WGCNA's implementation to avoid conflicts
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(
  norm.counts,
  maxBlockSize = 14000,
  TOMType = "signed",
  power = soft_power,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  randomSeed = 2025,
  verbose = 3
)

# Restore base cor() after network construction
cor <- temp_cor


# =====================
# 4. Module Eigengenes
# =====================

module_eigengenes <- bwnet$MEs

# Clean up sample names (remove "_merged" suffix if present)
rownames(module_eigengenes) <- str_remove(rownames(module_eigengenes), "_merged")

# Remove outlier samples identified during QC
module_eigengenes_2 <- module_eigengenes[
  !rownames(module_eigengenes) %in% c("MAL_RPNII_030", "MAL_RPNII_042", "MAL_SPIII_032"), ]

# Summarise module sizes and plot dendrogram with merged/unmerged colours
head(module_eigengenes)
table(bwnet$colors)

plotDendroAndColors(
  bwnet$dendrograms[[1]],
  cbind(bwnet$unmergedColors, bwnet$colors),
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)


# ==============================
# 5. Phenotype Data Preparation
# ==============================

annot_pheno <- read_tsv(here("annot_pheno.txt"))

# Select traits of interest and encode categorical variables numerically
annot_pheno_int <- annot_pheno %>%
  dplyr::select(
    SampleID, Gender, Group,
    Blood_Glucose:Blood_LDL, Blood_UricAcid,
    Bio_BodyFat:Bio_BMI, Bio_SubcutFat_Overall, Bio_MuscleMass_Overall
  ) %>%
  mutate(
    Sex = if_else(Gender == "M", 0, 1),
    Ancestry = if_else(Group == "Lundayeh", 1, 0),
    Lifestyle = if_else(Group == "PunanBatu", 0, 1)
  ) %>%
  dplyr::select(-Group, -Gender) %>%
  filter(SampleID %in% rownames(module_eigengenes_2))

# Align sample order to match module eigengene matrix
annot_pheno_int <- annot_pheno_int[
  match(rownames(module_eigengenes_2), annot_pheno_int$SampleID), ] %>%
  column_to_rownames("SampleID")

# Load gene annotation for module membership labelling
annot_gene <- read_tsv(here("out3_topTable_voom_PBvLDY_annot.txt")) %>%
  dplyr::select(Geneid, gene_name)


# ============================
# 6. Module-Trait Correlation
# ============================

# Compute Pearson correlations between module eigengenes and phenotypes
module.trait.corr <- cor(module_eigengenes_2, annot_pheno_int, use = "p")
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, 75)

# Adjust p-values within each module using BH correction
module.trait.corr.adjpvals <- module.trait.corr.pvals %>%
  as.data.frame() %>%
  rownames_to_column("colors") %>%
  gather(-colors, key = "pheno", value = "p_vals") %>%
  group_by(colors) %>%
  mutate(adj_p = p.adjust(p_vals, method = "BH")) %>%
  ungroup() %>%
  dplyr::select(colors, pheno, adj_p) %>%
  spread(key = "pheno", value = "adj_p")

# Scatter plot: lifestyle vs BMI module associations
ggplot(module.trait.corr.adjpvals,
       aes(-log10(Lifestyle), -log10(Bio_BMI), label = colors)) +
  geom_point() +
  geom_text_repel(max.overlaps = Inf, size = 10 / .pt,
                  box.padding = 0.5, min.segment.length = 0, color = "black") +
  geom_vline(xintercept = c(-log10(0.01), -log10(0.05)),
             linetype = "dotted", color = "red") +
  geom_hline(yintercept = c(-log10(0.01), -log10(0.05)),
             linetype = "dotted", color = "red") +
  theme_bw()


# ===================================================
# 7. Module-Trait Correlation Heatmap (CorLevelPlot)
# ===================================================

heatmap.data <- merge(module_eigengenes_2, annot_pheno_int, by = "row.names")

pl_corlevel <- CorLevelPlot(
  heatmap.data,
  x = names(heatmap.data)[21:35],
  y = names(heatmap.data)[2:20],
  rotLabX = 45,
  cexLabX = 0.6,
  cexLabY = 0.6,
  cexCorval = 0.4,
  cexLabColKey = 0.6,
  col = c("blue1", "skyblue", "white", "pink", "red1")
)

# ggsave(here("pl_corlevelplot.pdf"), pl_corlevel,
#        width = 174, height = 174, units = "mm", dpi = 300)


# =======================
# 8. Module-Gene Mapping
# =======================

# Map each gene to its assigned module colour and join gene name annotation
module.gene.mapping <- as.data.frame(bwnet$colors) %>%
  rownames_to_column("Geneid") %>%
  left_join(annot_gene, by = "Geneid")

# write.table(module.gene.mapping, file = here("module_gene_mapping.tsv"),
#             sep = "\t", row.names = TRUE, col.names = NA)


# ===========================================
# 9. Module Membership and Gene Significance
# ===========================================

# Compute module membership (correlation of each gene with each eigengene)
module.membership.measure <- cor(module_eigengenes, norm.counts, use = "p")

# Compute and tidy associated p-values
module.membership.measure.pvals.df <- corPvalueStudent(
    module.membership.measure, nrow(norm.counts)
  ) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Geneid") %>%
  left_join(annot_gene, by = "Geneid")

# write_tsv(module.membership.measure.pvals.df,
#           here("module_membership_pvals.txt"))

