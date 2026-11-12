# This repository provides R scripts supporting the manuscript "Gut microbiome assembly in very preterm infants predicts long-term growth trajectories".

## Software requirements

Analyses were performed using R version 4.6.1.

Required packages are listed at the beginning of each script. Major dependencies include `tidymodels`, `ranger`, `pROC`, `vegan`, `permute`, `Maaslin2`, `data.table`, `dplyr`, `tidyr`, and `ggplot2`.

## Analysis overview

### Figure S1 — PERMANOVA and MaAsLin2

Bray–Curtis dissimilarities are calculated from sample-level relative abundances. PERMANOVA evaluates clinical and nutritional associations using marginal tests, with sequencing depth included as specified in each model.
The script documents the permutation restrictions used for repeated infant sampling, the time intervals examined, and the group-dispersion tests. MaAsLin2 code illustrates the growth association model and the framework used to examine selected clinical and nutritional variables.

### Figure 2 — Random forest prediction

Random forest models use clinical predictors, microbiome features, or combined predictor sets. Models contain 1,000 trees, with `mtry` and `min_n` tuned using five-fold cross-validation repeated 10 times.
Hyperparameters are selected by maximizing mean cross-validation precision–recall area under the curve (PR-AUC). For the selected parameters, assessment predictions from all folds and repeats are pooled to calculate ROC-AUC and construct ROC curves. Each sample contributes one held-out prediction per repeat.
These ROC-AUC values are internal cross-validation estimates calculated from the resampling results used for tuning. Final full-data fits are used to extract permutation variable importance. Dataset-specific resampling and taxon-selection procedures are provided in the script.

### Figure 4 — Taxonomic contributions to HMO-related GHs

For each GH family, sample–strain annotated gene counts are expressed as percentages of that family’s total across samples. Contributions are summed across the selected families, averaged across samples within infants and then across infants within each growth group and interval.
A final global normalization makes the plotted contributions sum to 100% across all growth groups, time intervals, and taxonomic groups. These contributions describe annotated genomic functional potential.
