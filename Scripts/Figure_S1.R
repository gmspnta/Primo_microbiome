# PRIMO microbiome - Figure S1. PERMANOVA analysis
# Overall and interval-specific PERMANOVA of clinical and nutritional associations with microbiome composition, together with group-dispersion tests.

library(dplyr)
library(tibble)
library(purrr)
library(vegan)
library(Maaslin2)

# ============================================================
# 1) Load tables and prepare relative abundances

abundance <- data.table::fread(file.path(my_path, "abundance.tsv"), data.table = FALSE) %>%
  column_to_rownames("clade_name") %>% as.matrix() %>% t()  # samples × taxa

meta <- data.table::fread(file.path(my_path, "metadata.tsv"), data.table = FALSE)

# Check sample identifiers and align abundance rows to metadata
stopifnot(
  !anyNA(meta$lane_id),
  !anyDuplicated(meta$lane_id),
  !anyDuplicated(rownames(abundance)),
  all(meta$lane_id %in% rownames(abundance)))

abundance <- abundance[meta$lane_id, , drop = FALSE]
rownames(meta) <- meta$lane_id

# Convert each sample to relative abundance
abundance <- sweep(abundance, 1, rowSums(abundance), "/")

bray_dist <- vegdist(abundance, method = "bray")

# ============================================================
# 2) PERMANOVA

# Define clin_variables and nutri_variables explicitly above this section.
# Apply the original categorical coding before fitting these models.

run_permanova <- function(variables, block_by_infant = TRUE) {
  
  required <- unique(c(
    variables,
    if (block_by_infant) "ID"
  ))
  
  # Subset metadata and distances to the same complete-case samples
  keep <- complete.cases(meta[, required, drop = FALSE])
  meta_sub <- droplevels(meta[keep, , drop = FALSE])
  distance <- as.dist(as.matrix(bray_dist)[keep, keep, drop = FALSE])
  
  model_formula <- reformulate(variables, response = "distance")
  environment(model_formula) <- environment()
  
  vegan::adonis2(
    model_formula,
    data = meta_sub,
    by = "margin",
    permutations = 999,
    strata = if (block_by_infant) meta_sub$ID else NULL,
    parallel = 4
  )
}

set.seed(123)

# Clinical and nutritional variables: within-infant permutations
permanova_clin <- run_permanova(
  unique(c(clin_variables, "log10_read_count_postqc")))

permanova_nutri <- run_permanova(
  unique(c(nutri_variables, "log10_read_count_postqc")))

# ============================================================
# 3) Group dispersion

# Helper function to run betadisper + permutest
run_betadisper_test <- function(var, dist_obj, meta_df, id_var = "ID", nperm = 999) {
  
  # Pull grouping variable
  g <- meta_df[[var]]
  
  # Keep complete cases only (betadisper can't handle NA groups)
  keep <- !is.na(g)
  
  # Subset metadata and group
  meta_sub <- meta_df[keep, , drop = FALSE]
  g_sub <- factor(g[keep])
  
  # Need at least 2 groups
  if (nlevels(g_sub) < 2) {
    return(tibble(
      variable = var, n = sum(keep), n_groups = nlevels(g_sub),
      min_group_n = NA_integer_, max_group_n = NA_integer_,
      F_perm = NA_real_, p_perm = NA_real_,
      note = "Skipped: <2 non-missing groups"))
  }
  
  # Subset distance matrix to same samples
  dmat <- as.matrix(dist_obj)
  dmat_sub <- dmat[keep, keep, drop = FALSE]
  dist_sub <- as.dist(dmat_sub)
  
  # Group sizes
  group_sizes <- table(g_sub)
  perm_ctrl <- nperm
  
  # Run betadisper and permutation test
  bd <- betadisper(dist_sub, group = g_sub)
  pt <- permutest(bd, permutations = perm_ctrl)
  
  # Extract permutation F and p
  # (permutest(betadisper) stores them in $tab)
  F_perm <- unname(pt$tab[1, "F"])
  p_perm <- unname(pt$tab[1, "Pr(>F)"])
  
  tibble(
    variable = var,
    n = length(g_sub),
    n_groups = nlevels(g_sub),
    min_group_n = min(group_sizes),
    max_group_n = max(group_sizes),
    F_perm = F_perm,
    p_perm = p_perm,
    note = ifelse(any(group_sizes < 2), "Some groups have <2 samples", "")
  )
}

betadisper_results <- map_dfr(
  dispersion_vars,
  ~ run_betadisper_test(
    var = .x,
    dist_obj = bray_dist,
    meta_df = meta,
    nperm = 999)) %>%
  mutate(p_adj_BH = p.adjust(p_perm, method = "BH")) %>%
  arrange(p_perm)

# ============================================================
# 4) Maaslin associations

set.seed(123)

fit_growth <- Maaslin2::Maaslin2(
  input_data      = as.data.frame(abundance),
  input_metadata  = meta,
  output          = file.path("outputs", "maaslin2_growth"),
  fixed_effects   = c("weightZchange", "dol", "read_count_postqc"),
  random_effects  = c("ID"),
  min_prevalence  = 0.01,
  min_abundance   = 0.01,
  normalization  = "TSS",
  transform      = "LOG",
  standardize    = TRUE,
  analysis_method = "LM",
  correction     = "BH",
  max_significance = 0.1
)

# Clinical and nutritional variables identified through PERMANOVA
# were evaluated using analogous MaAsLin2 models.