# =========================================================
# RNA-seq Differential Expression Analysis
# Groups: PunanBatu (BLG), PunanTubu (RPN), Lundayeh (SPI)
# =========================================================

library(here)
library(tidyverse)
library(edgeR)
library(ggpubr)
library(scales)
library(GGally)
library(DeconCell)
library(limma)
library(variancePartition)
library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(enrichplot)


# ===============================
# 1. Load and Prepare Count Data
# ===============================

# Read featureCounts output and remove low-quality/excluded samples
fc <- read.table(
    here("counts", "ALL.basic.tsl1_3.revised.ftcnts.out"),
    header = TRUE, check.names = FALSE
  ) %>%
  dplyr::select(
    -(MPI_381:SMB_WNG_021),
    -BLGII_007, -BLGII_015,
    -MAL_RPNII_008, -MAL_RPNII_033, -MAL_RPNII_036,
    -MAL_SPIII_003, -MAL_SPIII_006, -MAL_SPIII_014,
    -MAL_SPIII_030, -MAL_SPIII_018
  )
# Note: MAL_RPNII_007, MAL_RPNII_021, MAL_RPNII_036 absent from DNA data

# Strip Ensembl version suffixes from gene IDs
fc$Geneid <- sub("\\.[^.]*$", "", fc$Geneid)
rownames(fc) <- fc$Geneid
fc <- fc[, -1]

# Assign group labels based on sample ID prefixes
group <- data.frame(id = colnames(fc)) %>%
  mutate(group = case_when(
    grepl("BLG", id) ~ "PunanBatu",
    grepl("RPN", id) ~ "PunanTubu",
    grepl("SPI", id) ~ "Lundayeh",
    TRUE             ~ "Ctrl"
  )) %>%
  pull(group)

# Load sample annotation, excluding the same removed samples
annot <- read.table(
    here("data_annot.txt"),
    sep = "\t", header = TRUE
  ) %>%
  dplyr::select(Group:Lane) %>%
  dplyr::filter(!Sample %in% c(
    "MPI_381", "MTW_TLL_013", "SMB_WNG_021",
    "BLGII_007", "BLGII_015",
    "MAL_RPNII_008", "MAL_RPNII_033", "MAL_RPNII_036",
    "MAL_SPIII_003", "MAL_SPIII_006", "MAL_SPIII_014",
    "MAL_SPIII_030", "MAL_SPIII_018"
  ))


# ======================================
# 2. DGEList Construction and Filtering
# ======================================

dge <- DGEList(counts = fc, group = factor(group))

# Report mean and median library sizes (in millions)
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
message(sprintf("Mean library size: %.2f M | Median: %.2f M", L, M))

# Compute CPM per group for group-aware filtering
cpm_dge <- cpm(dge)
cpm_PB <- cpm_dge[, grep("BLG", colnames(cpm_dge))]
cpm_PT <- cpm_dge[, grep("RPN", colnames(cpm_dge))]
cpm_LDY <- cpm_dge[, grep("SPI", colnames(cpm_dge))]

# Retain genes with CPM > 1 in at least 50% of samples within any group
keep.PB <- rowSums(cpm_PB > 1) >= (length(grep("BLG", rownames(dge$samples))) * 0.5)
keep.PT <- rowSums(cpm_PT > 1) >= (length(grep("RPN", rownames(dge$samples))) * 0.5)
keep.LDY <- rowSums(cpm_LDY > 1) >= (length(grep("SPI", rownames(dge$samples))) * 0.5)

keep <- keep.PB | keep.PT | keep.LDY
dge_filt <- dge[keep, , keep.lib.sizes = FALSE]

# Append additional sample covariates (excluding columns already present)
dge_filt$samples <- bind_cols(
  dge_filt$samples,
  dplyr::select(annot, -(Group:Sample))
)

# Compute log-CPM matrices (pre- and post-filtering, before normalisation)
log_cpm_dge <- cpm(dge, log = TRUE)
log_cpm_dge_filt <- cpm(dge_filt, log = TRUE)

# write.table(log_cpm_dge_filt, file = here("out3_log_cpm_dge_filt_ftcnts.txt"),
#             sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(log_cpm_dge, file = here("out3_log_cpm_dge_ftcnts.txt"),
#             sep = "\t", col.names = TRUE, row.names = TRUE)


# ============================
# 3. Sample Correlation Check
# ============================

corr_all <- ggcorr(
  log_cpm_dge_filt,
  method = c("pairwise", "spearman"),
  limits = c(0.9, 1),
  midpoint = 0.95,
  size = 2,
  hjust = 1,
  low = "white",
  mid = "yellow",
  high = "red",
  legend.position = "bottom"
)

# ggsave(plot = corr_all, filename = here("pl3_logcpm_sancheck_all.pdf"),
#        width = 174, height = 174, units = "mm", dpi = 300)
# ggsave(plot = corr_all, filename = here("pl3_logcpm_sancheck_all.png"),
#        width = 174, height = 174, units = "mm", dpi = 300)


# =====================
# 4. TMM Normalisation
# =====================

dge_filt_norm <- calcNormFactors(dge_filt, method = "TMM")
log_cpm_dge_filt_norm <- cpm(dge_filt_norm, log = TRUE)

# write.table(log_cpm_dge_filt_norm, file = here("out3_log_cpm_dge_filt_norm_ftcnts.txt"),
#             sep = "\t", col.names = TRUE, row.names = TRUE)


# =======================================
# 5. Cell-Type Deconvolution (DeconCell)
# =======================================

# Prepare expression matrix and run DeconCell prediction
dCell.exp <- dCell.expProcessing(dge_filt$counts, trim = TRUE)
data("dCell.models")
prediction <- dCell.predict(dCell.exp, dCell.models, res.type = "median")

all.predicted.cellcounts <- prediction$dCell.prediction
predicted.cellcounts <- prediction$dCell.prediction[, c(
  "Granulocytes", "B cells (CD19+)", "CD4+ T cells",
  "CD8+ T cells", "NK cells (CD3- CD56+)", "Monocytes (CD14+)"
)]

# write.table(predicted.cellcounts, file = here("out3_deconcell2_predicted.txt"),
#             sep = "\t", col.names = TRUE, row.names = TRUE)
# write.table(all.predicted.cellcounts, file = here("out3_deconcell2_all_predicted.txt"),
#             sep = "\t", col.names = TRUE, row.names = TRUE)

# Append predicted cell-type proportions to sample metadata
dge_filt_norm$samples <- bind_cols(dge_filt_norm$samples, predicted.cellcounts)

# write_tsv(rownames_to_column(dge_filt_norm$samples), here("out3_covariates.txt"))

# Add ancestry and lifestyle covariates
dge_filt_norm$samples <- dge_filt_norm$samples %>%
  mutate(
    anc = if_else(grepl("Punan", group), "PUN", "LDY"),
    lst = if_else(grepl("PunanBatu", group), "HG", "AGR")
  )


# ======================
# 6. Variance Partition
# ======================

# Sanitise cell-type column names for use in formulae
colnames(dge_filt_norm$samples) %<>%
  gsub("B cells \\(CD19\\+\\)", "CD19B", .) %>%
  gsub("CD4\\+ T cells", "CD4T", .) %>%
  gsub("CD8\\+ T cells", "CD8T", .) %>%
  gsub("NK cells \\(CD3\\- CD56\\+\\)", "NK", .) %>%
  gsub("Monocytes \\(CD14\\+\\)", "Mono", .)

# Ensure batch and lane are treated as categorical
dge_filt_norm$samples$Batch_Extraction <- as.character(dge_filt_norm$samples$Batch_Extraction)
dge_filt_norm$samples$Lane <- as.character(dge_filt_norm$samples$Lane)

# Build info data frame with sample IDs for variancePartition
info2 <- dge_filt_norm$samples %>%
  mutate(Sample = rownames(dge_filt_norm$samples))

# Fit a simple group design for voom
design_vp <- model.matrix(~ 0 + group, info2)
colnames(design_vp) <- gsub("group", "", colnames(design_vp))
vobjGenes <- voom(dge_filt_norm, design_vp)

# Define variance partition formula
form <- ~ (1 | group) + (1 | Batch_Extraction) + RIN + rRNA_ratio +
          Age + (1 | Sex) + (1 | Lane) +
          Granulocytes + CD19B + CD4T + CD8T + NK + Mono

varPart <- fitExtractVarPartModel(vobjGenes, form, info2)
plotVarPart(sortCols(varPart))


# =====================================
# 7. Design Matrix and Model Selection
# =====================================

# Build full design matrix with all covariates for model selection
design <- model.matrix(
  ~ 0 + group + RIN + rRNA_ratio + Age + Sex + Granulocytes + CD19B + CD4T + CD8T + NK + Mono,
  data = dge_filt_norm$samples
)

# Tidy up column names
colnames(design) <- gsub("lst", "", colnames(design))
colnames(design) %<>%
  gsub("`B cells \\(CD19\\+\\)`", "CD19B", .) %>%
  gsub("`CD4\\+ T cells`", "CD4T", .) %>%
  gsub("`CD8\\+ T cells`", "CD8T", .) %>%
  gsub("`NK cells \\(CD3\\- CD56\\+\\)`", "NK",   .) %>%
  gsub("`Monocytes \\(CD14\\+\\)`", "Mono", .)

# Define candidate model list for AIC-based selection
designlist <- list(
  Group = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3]),
  Group_Bio = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,6], Age = design[,7]),
  Group_RIN = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], RIN = design[,4]),
  Group_RIN_rat = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], RIN = design[,4], rRNA_rat = design[,5]),
  Group_Bio_RIN = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,6], Age = design[,7], RIN = design[,4]),
  Group_Bio_RIN_rat = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,6], Age = design[,7], RIN = design[,4], rRNA_rat = design[,5]),
  Group_Blood = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3],
    Gran = design[,8], Bcell = design[,9], CD4 = design[,10], CD8 = design[,11], NK = design[,12], Mono = design[,13]),
  Group_Bio_Blood = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,6], Age = design[,7],
    Gran = design[,8], Bcell = design[,9], CD4 = design[,10], CD8 = design[,11], NK = design[,12], Mono = design[,13]),
  Group_Bio_RIN_Blood = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,6], Age = design[,7], RIN = design[,4],
    Gran = design[,8], Bcell = design[,9], CD4 = design[,10], CD8 = design[,11], NK = design[,12], Mono = design[,13]),
  Group_Bio_RIN_rat_Blood = cbind(
    LDY = design[,1], PB = design[,2], PT = design[,3], Sex = design[,6], Age = design[,7], RIN = design[,4], rRNA_rat = design[,5], 
    Gran = design[,8], Bcell = design[,9], CD4 = design[,10], CD8 = design[,11], NK = design[,12], Mono = design[,13])
)

# Select best-fitting model per gene using AIC
out <- selectModel(dge_filt_norm$counts, designlist, criterion = "aic")
table(out$pref)


# ==========================================
# 8. Voom + Differential Expression (limma)
# ==========================================

# --- 8a. Three-group comparison (PunanBatu vs PunanTubu vs Lundayeh) ---

# Build three-group design matrix
design_3grp <- model.matrix(~ 0 + group, info2)
colnames(design_3grp) <- gsub("group", "", colnames(design_3grp))

# Define pairwise contrasts
contr.matrix.3grp <- makeContrasts(
  PBvPT  = PunanBatu - PunanTubu,
  PBvLDY = PunanBatu - Lundayeh,
  PTvLDY = PunanTubu - Lundayeh,
  levels = colnames(design_3grp)
)

# Voom transformation and model fitting
par(mfrow = c(1, 2))
v_3grp  <- voom(dge_filt_norm, design_3grp, plot = TRUE, save.plot = TRUE)
vfit    <- lmFit(v_3grp, design_3grp)
vfit2   <- contrasts.fit(vfit, contrasts = contr.matrix.3grp)
efit    <- eBayes(vfit2, robust = TRUE)
plotSA(efit, main = "Final model: Mean-variance trend")

# Summarise DE gene counts at various logFC thresholds
summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0))
summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5))
summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 1))

# Load gene annotation for joining to results
gene_annot <- read_tsv(here("gene_annot.txt")) %>%
  dplyr::select(Geneid, Length, gene_name) %>%
  mutate(Geneid = sub("\\.[^.]*$", "", Geneid))

# Extract top tables for each pairwise contrast
topTable_PB_PT  <- topTable(efit, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("Geneid") %>%
  left_join(gene_annot, by = "Geneid")

topTable_PB_LDY <- topTable(efit, coef = 2, n = Inf, sort.by = "p") %>%
  rownames_to_column("Geneid") %>%
  left_join(gene_annot, by = "Geneid")

topTable_PT_LDY <- topTable(efit, coef = 3, n = Inf, sort.by = "p") %>%
  rownames_to_column("Geneid") %>%
  left_join(gene_annot, by = "Geneid")

# write_tsv(topTable_PB_PT,  here("out3_topTable_voom_PBvPT.txt"))
# write_tsv(topTable_PB_LDY, here("out3_topTable_voom_PBvLDY.txt"))
# write_tsv(topTable_PT_LDY, here("out3_topTable_voom_PTvLDY.txt"))
# write.table(efit$p.value,  here("out3_efit_pvalue.txt"))
# write.table(vfit$coefficients, here("out3_vfit_coeff.txt"))


# --- 8b. Ancestry comparison (Punan vs Lundayeh) ---

design_anc <- model.matrix(~ 0 + anc, dge_filt_norm$samples)
colnames(design_anc) <- gsub("anc", "", colnames(design_anc))

contr.matrix.anc <- makeContrasts(
  PUNvLDY = PUN - LDY,
  levels  = colnames(design_anc)
)

v_anc   <- voom(dge_filt_norm, design_anc, plot = FALSE)
vfit_anc  <- lmFit(v_anc, design_anc)
vfit2_anc <- contrasts.fit(vfit_anc, contrasts = contr.matrix.anc)
efit_anc  <- eBayes(vfit2_anc, robust = TRUE)

topTable_PUN_LDY <- topTable(efit_anc, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("Geneid") %>%
  left_join(gene_annot, by = "Geneid")

# write_tsv(topTable_PUN_LDY, here("out3_topTable_voom_ANC_PUNvLDY.txt"))


# --- 8c. Lifestyle comparison (Hunter-Gatherer vs Agriculturalist) ---

design_lst <- model.matrix(~ 0 + lst, dge_filt_norm$samples)
colnames(design_lst) <- gsub("lst", "", colnames(design_lst))

contr.matrix.lst <- makeContrasts(
  HGvAGR = HG - AGR,
  levels = colnames(design_lst)
)

v_lst <- voom(dge_filt_norm, design_lst, plot = FALSE)
vfit_lst <- lmFit(v_lst, design_lst)
vfit2_lst <- contrasts.fit(vfit_lst, contrasts = contr.matrix.lst)
efit_lst <- eBayes(vfit2_lst, robust = TRUE)

topTable_HG_AGR <- topTable(efit_lst, coef = 1, n = Inf, sort.by = "p") %>%
  rownames_to_column("Geneid") %>%
  left_join(gene_annot, by = "Geneid")

# write_tsv(topTable_HG_AGR, here("out3_topTable_voom_LST_HGvAGR.txt"))


# ==================================
# 9. BiomaRt Annotation of DE Genes
# ==================================

attribute_names <- c(
  "ensembl_gene_id", "entrezgene_id", "description",
  "gene_biotype", "chromosome_name",
  "start_position", "end_position", "strand"
)

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)

# Retrieve annotations for DE genes and background gene set
annot_DEG_pop <- getBM(
  attributes = attribute_names,
  filters = "ensembl_gene_id",
  values = filter(topTable_PUN_LDY, adj.P.Val < 0.05) %>% pull(Geneid),
  mart = mart
)

annot_back_pop <- getBM(
  attributes = attribute_names,
  filters = "ensembl_gene_id",
  values = topTable_PUN_LDY$Geneid,
  mart = mart
)

# Join annotations back to DE results
tab_DEG_pop <- left_join(
    filter(topTable_PUN_LDY, adj.P.Val < 0.05 & (
    logFC < -0.1 | logFC > 0.1)),
    annot_DEG_pop,
    by = c("Geneid" = "ensembl_gene_id")
  ) %>%
  distinct(Geneid, .keep_all = TRUE)

tab_back_pop <- left_join(
    topTable_PUN_LDY,
    annot_back_pop,
    by = c("Geneid" = "ensembl_gene_id")
  ) %>%
  distinct(Geneid, .keep_all = TRUE)


# =============================================================================
# 10. Gene Ontology and KEGG Enrichment Analysis
# =============================================================================

# --- 10a. GO over-representation analysis (Punan vs Lundayeh) ---

ego_bp_PUNvLDY <- enrichGO(
  gene = tab_DEG_pop$Geneid,
  universe = tab_back_pop$Geneid,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE,
  pool = TRUE
)

# write_tsv(ego_bp_PUNvLDY@result, here("out3_GO_ANC_PUNvLDY.txt"))


# --- 10b. KEGG over-representation analysis (Punan vs Lundayeh) ---

kegg_bp_PUNvLDY <- enrichKEGG(
    gene = as.character(na.omit(tab_DEG_pop$entrezgene_id)),
    universe = as.character(na.omit(tab_back_pop$entrezgene_id)),
    organism = "hsa",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  ) %>%
  setReadable(., OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# write_tsv(kegg_bp_PUNvLDY@result, here("out3_KEGG_ANC_PUNvLDY.txt"))
# write_tsv(tab_back_pop,           here("out3_topTable_voom_ANC_PUNvLDY_annot.txt"))


# ================================================================
# 11. GSEA (Gene Set Enrichment Analysis) — Lifestyle (HG vs AGR)
# ================================================================

# Build ranked gene list using signed -log10(adj.P.Val)
d <- topTable_HG_AGR %>%
  dplyr::select(Geneid, logFC, adj.P.Val) %>%
  mutate(
    rnk = case_when(
      logFC > 0 ~ -log10(adj.P.Val),
      logFC < 0 ~ -(-log10(adj.P.Val))
    )
  ) %>%
  dplyr::select(Geneid, rnk) %>%
  dplyr::filter(!is.na(Geneid)) %>%
  group_by(Geneid) %>%
  filter(n() == 1) %>%   # retain only unique gene IDs
  ungroup()

# Format as named, sorted vector required by gseGO
geneList <- setNames(d$rnk, as.character(d$Geneid))
geneList <- sort(geneList, decreasing = TRUE)

# Run GSEA on GO Biological Process terms
ego_gsea_bp_HGvAGR <- gseGO(
    geneList = geneList,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05,
    verbose = FALSE,
    keyType = "ENSEMBL"
  ) %>%
  setReadable(., OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")

# write_tsv(ego_gsea_bp_HGvAGR@result, here("out3_gse_GO_LST_HGvAGR.txt"))

