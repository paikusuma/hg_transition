# =============================================================
# DNA Methylation Differential Analysis
# Groups: PunanBatu (BLG), PunanTubu (RPN/SPI), Lundayeh (SPI)
# Platform: Illumina EPIC array
# =============================================================

library(here)
library(tidyverse)
library(magrittr)
library(sesame)
library(limma)
library(missMethyl)
library(clusterProfiler)
library(ggpubr)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(DMRcate)


# =========================================
# 1. Load Beta Values and Probe Annotation
# =========================================

# Load probe annotation (chromosome, position, strand)
anno <- read_tsv(here("Result1.txt")) %>%
  dplyr::select(Name, CHR, MAPINFO, Strand_FR) %>%
  arrange(Name) %>%
  distinct() %>%
  mutate(Strand_FR = if_else(Strand_FR == "F", "+", "-"))

# Load beta values and average across duplicate probes
dat <- read_tsv(here("Result1.txt")) %>%
  dplyr::select(Name, BLGII_001:MAL_SPIII_039, CHR, MAPINFO)

dat1 <- dat %>%
  group_by(Name, CHR, MAPINFO) %>%
  summarise_at(vars(BLGII_001:MAL_SPIII_039), mean, na.rm = TRUE) %>%
  column_to_rownames(var = "Name")

# write_tsv(dat1 %>% rownames_to_column(var = "cg"), here("Result1_averaged.txt"))


# ==============================================================
# 2. Cell-Type Deconvolution (sesame / estimateCellComposition)
# ==============================================================

# Load reference cell-type methylation profiles
betas_celltype <- getRefSet(
  c("CD4T", "CD8T", "CD19B", "CD56NK", "CD14Monocytes", "granulocytes"),
  platform = "EPIC"
)

# Subset reference to probes present in data
betas_decon <- betas_celltype[rownames(betas_celltype) %in% rownames(dat1), ]

# Load cell composition estimation function
source(here("cell_composition.R"))

# Estimate cell-type fractions per sample
inds <- colnames(dat1)
dats <- list()

for (i in inds) {
  dat_decon <- dat1[rownames(dat1) %in% rownames(betas_celltype), i]
  est_decon <- estimateCellComposition(betas_decon, dat_decon, refine = TRUE, dichotomize = TRUE)
  est_decon_frac <- as.matrix(est_decon$frac)
  colnames(est_decon_frac) <- i
  dats[[i]] <- est_decon_frac
}

dats_all <- do.call(cbind, dats) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

# write_tsv(dats_all, here("out_blooddeconv_meth.txt"))


# ================================
# 3. Merge with Sample Annotation
# ================================

# Load sample annotation; correct one known mislabelled sample
annot <- read_table(here("data_annot.txt")) %>%
  mutate(Group = ifelse(Sample == "MAL_SPIII_035", "PunanTubu", Group))

# Join deconvolution results with annotation
all <- dats_all %>%
  left_join(annot, by = "Sample")


# ===================================================
# 4. Design Matrix and Model Selection — Full Cohort
# ===================================================

design <- model.matrix(
  ~ 0 + Group + Age + Sex +
    granulocytes + CD19B + CD4T + CD8T + CD56NK + CD14Monocytes,
  data = all
)
colnames(design) <- gsub("Group", "", colnames(design))

# Candidate model list for AIC/BIC selection
designlist <- list(
  Group = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3]),
  Group_Bio = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,5], Age = design[,4]),
  Group_Blood = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3],
    Gran = design[,6], Bcell = design[,7], CD4 = design[,8], CD8 = design[,9], NK = design[,10], Mono = design[,11]),
  Group_Bio_Blood = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,5], Age = design[,4],
    Gran = design[,6], Bcell = design[,7], CD4 = design[,8], CD8 = design[,9], NK = design[,10], Mono = design[,11])
)

out_aic <- selectModel(dat1, designlist, criterion = "aic")
out_bic <- selectModel(dat1, designlist, criterion = "bic")
table(out_aic$pref)
table(out_bic$pref)


# ==================================================================
# 5. Differential Methylation — Full Cohort (Three-Group Contrasts)
# ==================================================================

# Define pairwise contrasts
contr.matrix <- makeContrasts(
  PBvPT  = PunanBatu - PunanTubu,
  PBvLDY = PunanBatu - Lundayeh,
  PTvLDY = PunanTubu - Lundayeh,
  levels = colnames(design)
)

# Convert beta values to M-values for linear modelling
Mval <- log2(dat1 / (1 - dat1))

fit <- lmFit(Mval, design)
vfit <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(vfit, robust = TRUE)

# Summarise DE CpG counts at various logFC thresholds
summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0))
summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5))
summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 1))

# Extract top tables for each contrast
topTable_PB_PT  <- topTable(efit, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")
topTable_PB_LDY <- topTable(efit, coef = 2, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")
topTable_PT_LDY <- topTable(efit, coef = 3, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")

# write_tsv(topTable_PB_PT,  here("out_topTable_cg_PBvPT.txt"))
# write_tsv(topTable_PB_LDY, here("out
# write_tsv(topTable_PB_LDY, here("out_topTable_cg_PBvLDY.txt"))
# write_tsv(topTable_PT_LDY, here("out_topTable_cg_PTvLDY.txt"))
# write.table(efit$p.value,  here("out_efit_pvalue.txt"))


# ===============================================
# 6. GO, KEGG, and EWAS Enrichment — Full Cohort
# ===============================================

# Load EPIC array annotation for gometh
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

# Background CpG set (all tested probes)
back <- topTable_PB_LDY$cg

# Load EWAS catalogue reference for trait enrichment
ref_ewas <- read_tsv(here("ref_epic_wholeblood.txt")) %>%
  dplyr::select(Trait, CpG)

# Loop over pairwise contrasts for GO, KEGG, and EWAS enrichment
for (i in c("PB_PT", "PB_LDY", "PT_LDY")) {

  dat_sig <- get(paste0("topTable_", i)) %>%
    filter(abs(logFC) > 0.5 & adj.P.Val < 0.05)

  # GO enrichment (promoter-proximal probes only)
  go_enrich <- gometh(
      sig.cpg = dat_sig$cg,
      all.cpg = back,
      array.type = "EPIC",
      genomic.features = c("TSS200", "TSS1500", "5UTR", "1stExon"),
      collection = "GO",
      anno = epic_anno,
      sig.genes = TRUE
    ) %>%
    filter(FDR < 0.05)

  # KEGG enrichment (promoter-proximal probes only)
  kegg_enrich <- gometh(
      sig.cpg = dat_sig$cg,
      all.cpg = back,
      array.type = "EPIC",
      genomic.features = c("TSS200", "TSS1500", "5UTR", "1stExon"),
      collection = "KEGG",
      anno = epic_anno,
      sig.genes = TRUE
    ) %>%
    filter(FDR < 0.05)

  # EWAS catalogue trait enrichment
  enrich_ewas <- enricher(
    gene = dat_sig$cg,
    universe = back,
    maxGSSize = 47000,
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    TERM2GENE = ref_ewas
  )

  # write_tsv(go_enrich,   here(paste0("out_cg_promoter_GO_adjp005_lfc05_",   i, ".txt")))
  # write_tsv(kegg_enrich, here(paste0("out_cg_promoter_KEGG_adjp005_lfc05_", i, ".txt")))
  # write_tsv(enrich_ewas@result %>% filter(p.adjust < 0.05),
  #           here(paste0("out_cg_EWAS_adjp005_lfc01_", i, ".txt")))
}


# =====================================
# 7. Subset to RNA-seq Matched Samples 
# =====================================

# Load RNA-seq deconvolution results and clean sample names
deconv_rna <- read_tsv(here("..", "rna", "counts", "sanity_check3",
                             "out3_deconcell2_predicted.txt")) %>%
  mutate(Sample = gsub("_merged.*", "", Sample))

# Subset methylation data to samples present in RNA-seq data
names_use <- colnames(dat1)[colnames(dat1) %in% deconv_rna$Sample]
dat2 <- dat1[, names_use]

# write_tsv(as.data.frame(dat2) %>% rownames_to_column("TargetID"),
#           here("Result1_subset.txt"))

# Fix known sample rename (MAL_SPIII_035 -> MAL_RPNII_035)
dat2 <- dat2 %>%
  dplyr::rename(MAL_RPNII_035 = MAL_SPIII_035)

# Build covariate table merging annotation and RNA-seq cell fractions
annot_sub <- read_table(here("data_annot.txt"))

covar <- tibble(Sample = colnames(dat2)) %>%
  left_join(annot_sub,   by = "Sample") %>%
  left_join(deconv_rna,  by = "Sample") %>%
  mutate(Group = ifelse(Group == "PunanBatu", "HG", "AGR"))

# Sanitise RNA-seq cell-type column names to match methylation naming
colnames(covar) %<>%
  gsub("B cells \\(CD19\\+\\)", "CD19B", .) %>%
  gsub("CD4\\+ T cells", "CD4T", .) %>%
  gsub("CD8\\+ T cells", "CD8T", .) %>%
  gsub("NK cells \\(CD3\\- CD56\\+\\)", "CD56NK", .) %>%
  gsub("Monocytes \\(CD14\\+\\)", "CD14Mono", .) %>%
  gsub("Granulocytes", "Gran", .)


# ======================================================
# 8. Design Matrix and Model Selection — Matched Subset
# ======================================================

design_sub <- model.matrix(
  ~ 0 + Group + Age + Sex +
    Gran + CD19B + CD4T + CD8T + CD56NK + CD14Mono,
  data = covar
)
colnames(design_sub) <- gsub("Group", "", colnames(design_sub))

designlist_sub <- list(
  Group = cbind(
    LDY = design_sub[,1], PB = design_sub[,2], PT = design_sub[,3]),
  Group_Bio = cbind(
    LDY = design_sub[,1], PB = design_sub[,2], PT = design_sub[,3], Sex = design_sub[,5], Age = design_sub[,4]),
  Group_Blood = cbind(
    LDY = design_sub[,1], PB = design_sub[,2], PT = design_sub[,3],
    Gran = design_sub[,6], Bcell = design_sub[,7], CD4 = design_sub[,8], CD8 = design_sub[,9], NK = design_sub[,10], Mono = design_sub[,11]),
  Group_Bio_Blood = cbind(
    LDY = design_sub[,1], PB = design_sub[,2], PT = design_sub[,3], Sex = design_sub[,5], Age = design_sub[,4],
    Gran = design_sub[,6], Bcell = design_sub[,7], CD4 = design_sub[,8], CD8 = design_sub[,9], NK = design_sub[,10], Mono = design_sub[,11])
)

out_aic_sub <- selectModel(dat2, designlist_sub, criterion = "aic")
out_bic_sub <- selectModel(dat2, designlist_sub, criterion = "bic")
table(out_aic_sub$pref)
table(out_bic_sub$pref)


# ==================================================================
# 9. Differential Methylation — Matched Subset (Multiple Contrasts)
# ==================================================================

# --- 9a. Three-group pairwise contrasts ---

contr.matrix.3grp <- makeContrasts(
  PBvPT  = PunanBatu - PunanTubu,
  PBvLDY = PunanBatu - Lundayeh,
  PTvLDY = PunanTubu - Lundayeh,
  levels = colnames(design_sub)
)

# --- 9b. Ancestry contrast (Punan vs Lundayeh) ---

contr.matrix.anc <- makeContrasts(
  PUNvLDY = PUN - LDY,
  levels = colnames(design_sub)
)

# --- 9c. Lifestyle contrast (Hunter-Gatherer vs Agriculturalist) ---

contr.matrix.lst <- makeContrasts(
  HGvAGR = HG - AGR,
  levels = colnames(design_sub)
)

# Convert beta values to M-values for linear modelling
Mval_sub <- log2(dat2 / (1 - dat2))

# Fit models and apply each contrast in turn
fit_sub <- lmFit(Mval_sub, design_sub)

# --- Three-group fit ---
vfit_3grp <- contrasts.fit(fit_sub, contrasts = contr.matrix.3grp)
efit_3grp <- eBayes(vfit_3grp, robust = TRUE)

summary(decideTests(efit_3grp, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0))
summary(decideTests(efit_3grp, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5))
summary(decideTests(efit_3grp, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 1))

topTable_PB_PT  <- topTable(efit_3grp, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")
topTable_PB_LDY <- topTable(efit_3grp, coef = 2, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")
topTable_PT_LDY <- topTable(efit_3grp, coef = 3, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")

# write_tsv(topTable_PB_PT,  here("out3_topTable_cg_PBvPT.txt"))
# write_tsv(topTable_PB_LDY, here("out3_topTable_cg_PBvLDY.txt"))
# write_tsv(topTable_PT_LDY, here("out3_topTable_cg_PTvLDY.txt"))
# write.table(efit_3grp$p.value, here("out3_cg_efit_pvalue.txt"))

# --- Ancestry fit ---
vfit_anc <- contrasts.fit(fit_sub, contrasts = contr.matrix.anc)
efit_anc <- eBayes(vfit_anc, robust = TRUE)

topTable_PUN_LDY <- topTable(efit_anc, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")

# write_tsv(topTable_PUN_LDY, here("out3_topTable_cg_PUNvLDY.txt"))

# --- Lifestyle fit ---
vfit_lst <- contrasts.fit(fit_sub, contrasts = contr.matrix.lst)
efit_lst <- eBayes(vfit_lst, robust = TRUE)

summary(decideTests(efit_lst, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0))
summary(decideTests(efit_lst, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5))
summary(decideTests(efit_lst, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 1))

topTable_HG_AGR <- topTable(efit_lst, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("cg")

# write_tsv(topTable_HG_AGR, here("out3_topTable_cg_HGvAGR.txt"))


# ===================================================
# 10. GO, KEGG, and EWAS Enrichment — Matched Subset
# ===================================================

# Load EPIC annotation (if not already loaded)
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

# Background CpG set from the matched subset
back_sub <- topTable_PB_LDY$cg

# Load EWAS catalogue reference
ref_ewas <- read_tsv(here("ref_epic_wholeblood.txt")) %>%
  dplyr::select(Trait, CpG)

# Loop over pairwise contrasts for GO, KEGG, and EWAS enrichment
for (i in c("PB_PT", "PB_LDY", "PT_LDY")) {

  dat_sig <- get(paste0("topTable_", i)) %>%
    filter(abs(logFC) > 0.5 & adj.P.Val < 0.05)

  # GO enrichment (all genomic features)
  go_enrich <- gometh(
      sig.cpg = dat_sig$cg,
      all.cpg = back_sub,
      array.type = "EPIC",
      collection = "GO",
      anno = epic_anno,
      sig.genes = TRUE
    ) %>%
    filter(FDR < 0.05)

  # KEGG enrichment (all genomic features)
  kegg_enrich <- gometh(
      sig.cpg = dat_sig$cg,
      all.cpg = back_sub,
      array.type = "EPIC",
      collection = "KEGG",
      anno = epic_anno,
      sig.genes = TRUE
    ) %>%
    filter(FDR < 0.05)

  # EWAS catalogue trait enrichment
  enrich_ewas <- enricher(
    gene = dat_sig$cg,
    universe = back_sub,
    maxGSSize = 47000,
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    TERM2GENE = ref_ewas
  )

  # write_tsv(go_enrich,   here(paste0("out3_cg_GO_adjp005_lfc05_",   i, ".txt")))
  # write_tsv(kegg_enrich, here(paste0("out3_cg_KEGG_adjp005_lfc05_", i, ".txt")))
  # write_tsv(enrich_ewas@result %>% filter(p.adjust < 0.05),
  #           here(paste0("out3_cg_EWAS_adjp005_lfc05_", i, ".txt")))
}


# =================================================
# 11. Differentially Methylated Regions (DMRcate)
# =================================================

# Remove SNP-associated and cross-hybridising probes
ALLMs_noSNPs <- rmSNPandCH(Mval_sub, mafcut = 0.05, rmcrosshyb = FALSE)

# --- 11a. DMR analysis: PunanTubu vs Lundayeh ---

cpg_annot_PTvLDY <- cpg.annotate(
  datatype = "array",
  object = ALLMs_noSNPs,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  contrasts = TRUE,
  cont.matrix = contr.matrix.3grp,
  design = design_sub,
  coef = "PTvLDY",
  fdr = 0.05
)

DMR_PTvLDY <- dmrcate(cpg_annot_PTvLDY, lambda = 1000, C = 2)
DMR_range_PTvLDY <- extractRanges(DMR_PTvLDY, genome = "hg38")

GO_DMR_PTvLDY <- goregion(
  DMR_range_PTvLDY,
  all.cpg    = rownames(ALLMs_noSNPs),
  collection = "GO",
  array.type = "EPIC",
  plot.bias  = TRUE
)

# write_tsv(as.data.frame(DMR_range_PTvLDY), here("out3_DMR_PTvLDY.txt"))
# write_tsv(as.data.frame(GO_DMR_PTvLDY),    here("out3_DMR_GO_PTvLDY.txt"))

# --- 11b. DMR analysis: PunanBatu vs Lundayeh ---

cpg_annot_PBvLDY <- cpg.annotate(
  datatype = "array",
  object = ALLMs_noSNPs,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  contrasts = TRUE,
  cont.matrix = contr.matrix.3grp,
  design = design_sub,
  coef = "PBvLDY",
  fdr = 0.05
)

DMR_PBvLDY <- dmrcate(cpg_annot_PBvLDY, lambda = 1000, C = 2)
DMR_range_PBvLDY <- extractRanges(DMR_PBvLDY, genome = "hg38")

GO_DMR_PBvLDY <- goregion(
  DMR_range_PBvLDY,
  all.cpg = rownames(ALLMs_noSNPs),
  collection = "GO",
  array.type = "EPIC",
  plot.bias = TRUE
)

# write_tsv(as.data.frame(DMR_range_PBvLDY), here("out3_DMR_PBvLDY.txt"))
# write_tsv(as.data.frame(GO_DMR_PBvLDY),    here("out3_DMR_GO_PBvLDY.txt"))

# --- 11c. DMR analysis: PunanBatu vs PunanTubu ---

cpg_annot_PBvPT <- cpg.annotate(
  datatype = "array",
  object = ALLMs_noSNPs,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  contrasts = TRUE,
  cont.matrix = contr.matrix.3grp,
  design = design_sub,
  coef = "PBvPT",
  fdr = 0.05
)

DMR_PBvPT <- dmrcate(cpg_annot_PBvPT, lambda = 1000, C = 2)
DMR_range_PBvPT <- extractRanges(DMR_PBvPT, genome = "hg38")

GO_DMR_PBvPT <- goregion(
  DMR_range_PBvPT,
  all.cpg = rownames(ALLMs_noSNPs),
  collection = "GO",
  array.type = "EPIC",
  plot.bias = TRUE
)

# write_tsv(as.data.frame(DMR_range_PBvPT), here("out3_DMR_PBvPT.txt"))
# write_tsv(as.data.frame(GO_DMR_PBvPT),    here("out3_DMR_GO_PBvPT.txt"))

# --- 11d. DMR analysis: Lifestyle (HG vs AGR) ---

cpg_annot_HGvAGR <- cpg.annotate(
  datatype = "array",
  object = ALLMs_noSNPs,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  contrasts = TRUE,
  cont.matrix = contr.matrix.lst,
  design = design_sub,
  coef = "HGvAGR",
  fdr = 0.05
)

DMR_HGvAGR <- dmrcate(cpg_annot_HGvAGR, lambda = 1000, C = 2)
DMR_range_HGvAGR <- extractRanges(DMR_HGvAGR, genome = "hg38")

GO_DMR_HGvAGR <- goregion(
  DMR_range_HGvAGR,
  all.cpg = rownames(ALLMs_noSNPs),
  collection = "GO",
  array.type = "EPIC",
  plot.bias = TRUE
)

# write_tsv(as.data.frame(DMR_range_HGvAGR), here("out3_DMR_LST_HGvAGR.txt"))
# write_tsv(as.data.frame(GO_DMR_HGvAGR),    here("out3_DMR_GO_LST_HGvAGR.txt"))

# --- 11e. DMR analysis: Ancestry (Punan vs Lundayeh) ---

cpg_annot_PUNvLDY <- cpg.annotate(
  datatype = "array",
  object = ALLMs_noSNPs,
  what = "M",
  arraytype = "EPICv1",
  analysis.type = "differential",
  contrasts = TRUE,
  cont.matrix = contr.matrix.anc,
  design = design_sub,
  coef = "PUNvLDY",
  fdr = 0.05
)

DMR_PUNvLDY <- dmrcate(cpg_annot_PUNvLDY, lambda = 1000, C = 2)
DMR_range_PUNvLDY <- extractRanges(DMR_PUNvLDY, genome = "hg38")

GO_DMR_PUNvLDY <- goregion(
  DMR_range_PUNvLDY,
  all.cpg = rownames(ALLMs_noSNPs),
  collection = "GO",
  array.type = "EPIC",
  plot.bias = TRUE
)

# write_tsv(as.data.frame(DMR_range_PUNvLDY), here("out3_DMR_ANC_PUNvLDY.txt"))
# write_tsv(as.data.frame(GO_DMR_PUNvLDY),    here("out3_DMR_GO_ANC_PUNvLDY.txt"))
