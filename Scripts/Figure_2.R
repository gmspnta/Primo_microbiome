# PRIMO microbiome - Random Forest models for Figure 2
# Random forest prediction of growth trajectory using clinical variables, microbiome features, and combined predictors, including time-window and selected-taxa analyses.

# ========================
# 1. Load packages quietly
suppressPackageStartupMessages({
  library(tidymodels)  
  library(dplyr)      
  library(tibble)     
  library(purrr)       
  library(stringr)     
  library(readr)
  library(compositions)
  library(data.table)
})

tidymodels_prefer()

# ========================
# 2. Import data
my_path <- "Data/"

strain.abundance.all <- fread(file.path(my_path, "abundance.tsv")) %>% tibble::column_to_rownames(var = "clade_name")

meta.data <- fread(file.path(my_path, "metadata.tsv"))

# Define factors for outcomes
yesno_vars <- c(
  grep("^abx[0-9]*_yn$", names(meta.data), value = TRUE),
  "prom", "aph", "mathtn", "matdm", "bsi", "nec", "mortality",
  "chorio", "abxstatus", "mctoil", "FTT_2y", "IVH", "ROP", "PDA",
  "CLDstatus", "MVstatus", "probx_yn", "TPN_yn", "HMF_yn",
  "Hkcal_yn", "MCT_yn"
)

# Set variables as factors
meta.data <- meta.data %>% mutate(
  across(all_of(yesno_vars),     # Recode Yes/No variables to 0/1
         ~ dplyr::case_when(.x == "No" ~ 0, .x == "Yes" ~ 1, TRUE ~ NA_real_)),
  gender = case_when(gender == "Male" ~ 0, gender == "Female" ~ 1, TRUE ~ NA_real_),
  deliverymode = case_when(deliverymode == "NVD" ~ 0, deliverymode == "LSCS" ~ 1, TRUE ~ NA_real_),
  race = factor(race, levels = c("Malay", "Chinese", "Indian", "Others")),
  CLDseverity = case_when(
    CLDseverity == "No" ~ 0,
    CLDseverity == "Mild" ~ 1,
    CLDseverity == "Moderate" ~ 2,
    CLDseverity == "Severe" ~ 3,
    TRUE ~ NA_real_
  ),
  classification_bsid = case_when(
    is.na(classification_bsid) ~ NA_real_,
    classification_bsid == "None" ~ 0,
    classification_bsid == "Moderate" ~ 1,
    classification_bsid == "Severe" ~ 2,
    TRUE ~ NA_real_
  ),
  
  # Convert to factors
  across(all_of(yesno_vars),
         ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))),
  gender = factor(gender, levels = c(0, 1), labels = c("Male", "Female")),
  deliverymode = factor(deliverymode, levels = c(0, 1), labels = c("NVD", "LSCS")),
  classification_bsid = factor(classification_bsid, levels = c(0, 1, 2), labels = c("None", "Moderate", "Severe"), ordered = TRUE),
  CLDseverity = factor(CLDseverity, levels = c(0, 1, 2, 3), labels = c("No", "Mild", "Moderate", "Severe"), ordered = TRUE)
)

# Define traj_outcome as binary 0/1 where outcome of interest is Normal growth
meta.data <- meta.data %>%
  mutate(
    traj_outcome = if_else(trajectory2 == 1, "1", "0"), # Normal growth is 1, Poor growth is 0
    traj_outcome = factor(traj_outcome, levels = c("0", "1")))

# ID as factor 
meta.data$ID <- as.factor(meta.data$ID)

# ========================
# 3. Feature filtering and CLR transform:
meta.data.sel <- meta.data %>% select(c("lane_id", "ID", "dol", "corrGA",
                                        "GA", "EBM", "PDHM", "formula",
                                        "lengthofstay", "mortality", "deliverymode", "race", "gender",
                                        "prom", "bsi", "nec", "chorio", "aph", "mathtn", "matdm", 
                                        "IVH", "ROP", "PDA", "MVstatus", "MVduration", "CLDstatus", "CLDseverity",
                                        "dolstartfeeds", "dolfullfeeds", "tpnduration", "hmfduration", "27kcalduration",
                                        "mctoil",  "weightZchange",
                                        "abxstatus", "abx_duration", "abx_yn", "probx_yn",
                                        "TPN_yn", "HMF_yn", "Hkcal_yn", "MCT_yn",
                                        "traj_outcome"))

# Clinical data
dat_clin <- meta.data.sel 
dat_clin$ID <- as.factor(dat_clin$ID) # ensure ID as factor

# CLR-transform microbiome abundances
strain.clr.mat <- clr(t(strain.abundance.all) + 1e-6)

# Convert to data frame and keep lane_id as explicit column
strain.clr.mat <- as.data.frame(strain.clr.mat) %>% rownames_to_column(var = "lane_id")

# Microbiome data table
dat_micro <- strain.clr.mat %>% inner_join(meta.data.sel %>% select(lane_id, ID, traj_outcome), by = "lane_id")

# Clinical data + microbiome data
dat <- strain.clr.mat %>% inner_join(meta.data.sel, by = "lane_id")

## Filter microbiome data to top predictive taxa only
tree <- ape::read.tree(file.path(my_path, "top_taxa_tree.nwk"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

## Select 30 taxa (Include B. breve commensal and exclude B. breve M16V
strain.abundance.top30 <- strain.abundance.all[rownames(strain.abundance.all) %in% c(setdiff(tree$tip.label, "Bifidobacterium_breve_M16V"), "Bifidobacterium breve_C12"), , drop = FALSE]

# Repeat CLR-transformation and join with clinical data
strain.clr.top30 <- clr(t(strain.abundance.top30) + 1e-6)

strain.clr.top30 <- as.data.frame(strain.clr.top30) %>% rownames_to_column(var = "lane_id")

dat_micro_top30 <- strain.clr.top30 %>% inner_join(meta.data.sel %>% select(lane_id, ID, traj_outcome), by = "lane_id")

# Combined clin + microbiome data table with top 30 taxa
dat_top30 <- strain.clr.top30 %>% inner_join(meta.data.sel, by = "lane_id")

# Two-week windows based on day of life (dol):
# weeks 1–2 = days 0–13; weeks 2–3 = days 7–20.
# Subsequent windows shift both limits forward by 7 days:
# weeks 3–4 = days 14–27, weeks 4–5 = days 21–34, etc.

# Weeks 1–2
lane_ids_w12 <- dat_clin %>% filter(dol >= 0, dol <= 13) %>% arrange(ID, dol) %>% pull(lane_id)

strain_abundance_w12 <- strain.abundance.all[ , colnames(strain.abundance.all) %in% lane_ids_w12, drop = FALSE]
strain_abundance_w12 <- strain_abundance_w12[rowMeans(strain_abundance_w12) > 0, , drop = FALSE]

dat_w12 <- compositions::clr(t(strain_abundance_w12) + 1e-6) %>%
  as.data.frame() %>%
  rownames_to_column(var = "lane_id") %>%
  inner_join(meta.data.sel, by = "lane_id")

# Weeks 2–3
lane_ids_w23 <- dat_clin %>% filter(dol >= 7, dol <= 20) %>% arrange(ID, dol) %>% pull(lane_id)

strain_abundance_w23 <- strain.abundance.all[ , colnames(strain.abundance.all) %in% lane_ids_w23, drop = FALSE]
strain_abundance_w23 <- strain_abundance_w23[rowMeans(strain_abundance_w23) > 0, , drop = FALSE]

dat_w23 <- compositions::clr(t(strain_abundance_w23) + 1e-6) %>%
  as.data.frame() %>%
  rownames_to_column(var = "lane_id") %>%
  inner_join(meta.data.sel, by = "lane_id")

# For each subsequent window, repeat the same steps:
# select samples, remove taxa with zero mean abundance within the window,
# add the pseudocount, apply CLR, and join the clinical metadata.

# ========================
# 4. Repeated 5-fold cross-validation (10 repeats)

# Clinical data: keep all samples from each infant in the same fold
set.seed(1)
folds_clin <- group_vfold_cv(dat_clin, v = 5, repeats = 10, group = ID)

# Microbiome data: stratify folds by growth outcome
set.seed(1)
folds_micro <- vfold_cv(dat_micro, v = 5, repeats = 10, strata = traj_outcome)

# Apply the same outcome-stratified procedure to combined data,
# each two-week window, and datasets restricted to the selected taxa.
# Reset set.seed(1) before each call, replacing the input dataset
# and output object name as appropriate.

# =======================
# 5. Select co-variates to include in model
# Select co-variates available at baseline and exclude all co-morbidity outcomes
covars <- c("dolstartfeeds", "dol", "GA", "abx_yn", "EBM", "deliverymode", "race", "gender", "traj_outcome")

# Columns remaining that are NOT to be treated as microbiome features
meta_cols <- c("corrGA", "formula",
               "bsi", "nec" , "IVH", "ROP", "PDA", "MVstatus", "MVduration", "CLDstatus", "lengthofstay",
               "mortality", "CLDseverity", "hmfduration",  "PDHM",
               "27kcalduration", "mctoil", "abxstatus", "abx_duration", 
               "prom", "chorio", "aph", "mathtn", "matdm",
               "probx_yn", "TPN_yn", "HMF_yn", "Hkcal_yn", "MCT_yn",
               "lane_id", "ID",  "weightZchange", "dolfullfeeds", "tpnduration",
               covars)

# =======================
# 6. Select numeric microbiome features
get_micro_features <- function(data) {
  candidates <- setdiff(names(data), meta_cols)
  candidates[vapply(data[, candidates, drop = FALSE], is.numeric, logical(1))]
}

feature_cols     <- get_micro_features(dat_micro)
feature_w12_cols <- get_micro_features(dat_w12)
feature_w23_cols <- get_micro_features(dat_w23)

# Apply the same function to subsequent windows and top-taxa datasets.
# Example:
# feature_top30_cols <- get_micro_features(dat_micro_top30)

# =======================
# 7. Build recipes
# All models use variance filtering, median imputation, and normalization.
# Clinical and combined models additionally impute and encode categorical predictors.
make_recipe <- function(data, include_clinical = FALSE) {
  rec <- recipe(traj_outcome ~ ., data = data)
  
  if ("ID" %in% names(data)) {rec <- rec %>% update_role(ID, new_role = "id")}
  
  rec <- rec %>%
    step_zv(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors())
  
  if (include_clinical) {
    rec <- rec %>%
      step_impute_mode(all_nominal_predictors()) %>%
      step_dummy(all_nominal_predictors())}
  
  rec %>% step_normalize(all_numeric_predictors())}

# Clinical model
rec_clin <- make_recipe(
  dat_clin %>% select(-lane_id, -ID, -weightZchange),
  include_clinical = TRUE)

# Microbiome model
rec_micro <- make_recipe(
  dat_micro %>% select(-lane_id, -ID))

# Combined model with selected clinical covariates
rec_combined <- make_recipe(dat %>% select(all_of(c(feature_cols, covars, "ID"))), include_clinical = TRUE)

# Combined models for the first two time windows
rec_w12 <- make_recipe(dat_w12 %>% select(all_of(c(feature_w12_cols, covars, "ID"))), include_clinical = TRUE)

rec_w23 <- make_recipe(dat_w23 %>% select(all_of(c(feature_w23_cols, covars, "ID"))), include_clinical = TRUE)

# For subsequent windows and top-taxa datasets, use the same functions
# with the corresponding dataset and microbiome feature list.
# Microbiome-only models use the default include_clinical = FALSE.

# =======================
# 8. Random forest classification model
mod_rf <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
  set_engine("ranger", importance = "permutation") %>% set_mode("classification")

# Combine each recipe with the random forest model
make_rf_workflow <- function(rec) {workflow() %>% add_recipe(rec) %>% add_model(mod_rf)}

# Clinical, microbiome, and combined models
wf_rf_clin     <- make_rf_workflow(rec_clin)
wf_rf_micro    <- make_rf_workflow(rec_micro)
wf_rf_combined <- make_rf_workflow(rec_combined)

# First two time windows
wf_rf_w12 <- make_rf_workflow(rec_w12)
wf_rf_w23 <- make_rf_workflow(rec_w23)

# Top-taxa models
wf_rf_micro_top30    <- make_rf_workflow(rec_micro_top30)
wf_rf_combined_top30 <- make_rf_workflow(rec_combined_top30)

# Repeat for subsequent windows using the corresponding recipe.

# =========================
# 9. Metric and tuning control
conflicted::conflicts_prefer(yardstick::spec)

metrics <- metric_set(roc_auc, pr_auc, yardstick::accuracy, sens, spec)

ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything")

grid_rf <- grid_regular(
  mtry(range = c(5, min(20, length(feature_cols)))), # mtry bounded by #features
  min_n(range = c(2, 20)),
  levels = 6)

# =========================
# 10. Tune random forests with repeated cross-validation

tune_rf <- function(wf, folds, grid = grid_rf) {
  set.seed(123)
  tune_grid(wf, resamples = folds, grid = grid, metrics = metrics, control = ctrl)}

# Clinical, microbiome, and combined datasets
res_rf_clin     <- tune_rf(wf_rf_clin, folds_clin)
res_rf_micro    <- tune_rf(wf_rf_micro, folds_micro)
res_rf_combined <- tune_rf(wf_rf_combined, folds_combined)

# Combined datasets for two-week windows
res_rf_w12 <- tune_rf(wf_rf_w12, folds_w12)
res_rf_w23 <- tune_rf(wf_rf_w23, folds_w23)
# Repeat for subsequent windows using their workflows and matching folds.

# Top-30-taxa datasets
res_rf_micro_top30 <- tune_rf(wf_rf_micro_top30, folds_micro_top30)

res_rf_combined_top30 <- tune_rf(wf_rf_combined_top30, folds_combined_top30)

# Top-30-taxa datasets within selected windows e.g.
res_rf_micro_top30_w34 <- tune_rf(wf_rf_micro_top30_w34, folds_micro_top30_w34)

res_rf_combined_top30_w34 <- tune_rf(wf_rf_combined_top30_w34, folds_combined_top30_w34)
# Repeat for subsequent windows using their workflows and matching folds.

# =======================================================
# 11. Select best parameters, extract importance and pooled out-of-fold ROC-AUC values

extract_rf <- function(wf, res, data) {
  
  # 1. Select parameters using mean cross-validation PR-AUC
  best <- select_best(res, metric = "pr_auc")
  
  # 2. Fit the final model and extract permutation importance
  set.seed(123)
  fit <- finalize_workflow(wf, best) %>% fit(data)
  
  importance <- extract_fit_parsnip(fit)$fit$variable.importance
  vip <- tibble(
    feature = names(importance),
    importance = unname(importance)
  ) %>%
    arrange(desc(importance))
  
  # 3. Pool held-out predictions across all folds and repeats
  pred <- collect_predictions(res, parameters = best)
  
  roc <- pROC::roc(
    response = pred$traj_outcome, predictor = pred$.pred_1, levels = c("0", "1"), direction = "<", quiet = TRUE)
  
  list(best = best, fit = fit, vip = vip, pred = pred, roc = roc, auc = as.numeric(pROC::auc(roc)))
}

# Clinical, microbiome, and combined models
rf_clin     <- extract_rf(wf_rf_clin, res_rf_clin, dat_clin)
rf_micro    <- extract_rf(wf_rf_micro, res_rf_micro, dat_micro)
rf_combined <- extract_rf(wf_rf_combined, res_rf_combined, dat)

# Top-30-taxa models
rf_micro_top30 <- extract_rf(
  wf_rf_micro_top30, res_rf_micro_top30, dat_micro_top30
)
rf_combined_top30 <- extract_rf(
  wf_rf_combined_top30, res_rf_combined_top30, dat_top30
)

# Repeat for each time window using its matching workflow,
# tuning result, and analysis dataset.

# Access results
rf_clin$best                         # Selected hyperparameters
rf_clin$vip                          # Ranked feature importance
rf_clin$auc                          # Pooled out-of-fold ROC-AUC
rf_micro_top30$vip %>% slice_head(n = 60)

# ============================================================
# Select top common predictive taxa
# Filtered by presence in ≥4 of 8 windows and ranked in the overall top 60

importance_cols <- paste0("importance_", c("w12", "w23", "w34", "w45", "w56", "w67", "w78", "w810"))

top_taxa <- rf_micro$vip %>%
  slice_head(n = 60) %>%
  pull(feature)

rf_features_filt <- rf_features %>%
  mutate(n_windows = rowSums(!is.na(pick(all_of(importance_cols))))) %>%
  filter((n_windows >= 4 & feature %in% top_taxa))

# ============================================================
# 12. Plot pooled out-of-fold ROC curves

plot_rf_roc <- function(models, title = NULL) {
  labels <- sprintf("%s (AUC %.2f)", names(models), vapply(models, function(x) x$auc, numeric(1)))
  
  pROC::ggroc(purrr::map(models, "roc"), legacy.axes = TRUE, linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey70") +
    scale_colour_manual(values = c("grey25", "#3BB24C", "#da6310", "#0073C2FF"),
      breaks = names(models), labels = labels) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = title, x = "1 − Specificity", y = "Sensitivity", colour = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "vertical")
}

# Supply models in the desired display order
p_roc <- plot_rf_roc(
  list("Clinical" = rf_clin,
    "Top 30 taxa" = rf_micro_top30,
    "Top 30 taxa (week 4–5)" = rf_micro_top30_w45,
    "Clin+top 30 taxa (week 4–5)" = rf_combined_top30_w45),
  title = "Discharge outcome")

p_roc

p_roc_2y <- plot_rf_roc(
  list("Clinical" = rf_clin_2y,
    "Top 30 taxa" = rf_micro_top30_2y,
    "Top 30 taxa (week 4–5)" = rf_micro_top30_w45_2y,
    "Clin+top 30 taxa (week 4–5)" = rf_combined_top30_w45_2y),
  title = "2-year outcome")

library(patchwork)
(p_roc + p_roc_2y) + plot_layout(axes = "collect", axis_titles = "collect")
