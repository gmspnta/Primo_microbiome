# # PRIMO microbiome - Figure 4. Species/strain-level contribution to HMO-related CAZymes
#
# Annotated gene counts are summed across six selected GH families.
# Counts are averaged across samples within each infant/interval,
# then across infants within each growth group/interval.
# A single final normalization makes contributions sum to 100%
# across the entire figure.
#
# Inputs:
# - hqmags_CAZymes.tsv: gene counts with sample and taxonomic labels
# - metadata.tsv: Run, ID, week, and trajectory2
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)})

# 1) Load inputs
my_path <- "Data/"

df <- fread(file.path(my_path, "hqmags_CAZymes.tsv"))

metadata_select <- fread(file.path(my_path, "metadata.tsv")) %>%
  select(lane_id, ID, week, trajectory2)

stopifnot(!anyDuplicated(metadata_select$Run))

sig_cazy <- c("GH123", "GH98", "GH29", "GH95", "GH89", "GH112")

# 2) Calculate contributions within each GH family
merged_df_filt <- df %>% filter(CAZy_family %in% sig_cazy)

contrib_overall <- merged_df_filt %>%
  group_by(Sample, CAZy_family, Strain, Species, Genus, Family, Order, Class, Phylum) %>%
  summarise(count_sum = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(CAZy_family) %>%
  mutate(family_total = sum(count_sum, na.rm = TRUE),
    pct_contrib = if_else(
      family_total > 0,
      100 * count_sum / family_total,
      NA_real_)) %>% ungroup()

# Each nonzero GH family should total 100%
contrib_overall %>% group_by(CAZy_family) %>%
  summarise(total_pct = sum(pct_contrib, na.rm = TRUE))

# 3) Add metadata and define plotting groups
# Use plot_cols matching these group names, in the desired order.
plot_group_levels <- c("Bifidobacterium breve M16V",
                       "Bifidobacterium breve commensal",
                       "Bifidobacterium longum",
                       "Bifidobacterium infantis",
                       "Other Bifidobacterium spp",
                       "Sarcina perfringens", 
                       "Clostridium neonatale",
                       "Clostridium butyricum",
                       "Clostridium paraputrificum",
                       "Other Clostridium spp",
                       "Cutibacterium avidum",
                       "Blautia hominis",
                       "Bacteroides uniformis",
                       "Other")  

contrib_plot_base <- contrib_overall %>%
  left_join(metadata_select, by = c("Sample" = "lane_id")) %>%
  mutate(
    week_interval = case_when(
      week %in% 1:2 ~ "1–2",
      week %in% 3:4 ~ "3–4",
      week %in% 5:6 ~ "5–6",
      week %in% 7:9 ~ "7–9",
      TRUE ~ NA_character_
    ),
    week_interval = factor(
      week_interval, levels = c("1–2", "3–4", "5–6", "7–9")
    ),
    trajectory2 = factor(
      trajectory2,
      levels = c("1", "2"),
      labels = c("Normal growth", "Poor growth")
    ),
    plot_group = case_when(
      Genus == "Bifidobacterium" & Strain == "Bifidobacterium breve M16V" ~ "Bifidobacterium breve M16V",
      Genus == "Bifidobacterium" & Strain == "Bifidobacterium breve commensal" ~ "Bifidobacterium breve commensal",
      Genus == "Bifidobacterium" & Strain == "Bifidobacterium infantis" ~ "Bifidobacterium infantis",
      Genus == "Bifidobacterium" & Strain == "Bifidobacterium longum" ~ "Bifidobacterium longum",
      Genus == "Bifidobacterium" ~ "Other Bifidobacterium spp",
      
      Genus == "Sarcina" & Strain == "Sarcina perfringens" ~ "Sarcina perfringens",
      Genus == "Clostridium" & Strain == "Clostridium neonatale" ~ "Clostridium neonatale",
      Genus == "Clostridium" & Strain == "Clostridium butyricum" ~ "Clostridium butyricum",
      Genus == "Clostridium" & Strain == "Clostridium paraputrificum" ~ "Clostridium paraputrificum",
      Genus == "Clostridium" ~ "Other Clostridium spp",
      Genus == "Clostridioides" ~ "Other Clostridium spp",
      
      Genus == "Blautia" & Strain == "Blautia hominis" ~ "Blautia hominis",
      
      Genus == "Cutibacterium" & Strain == "Cutibacterium avidum" ~ "Cutibacterium avidum",
      # Genus == "Cutibacterium" ~ "Cutibacterium_other",
      
      Genus == "Bacteroides" & Strain == "Bacteroides uniformis" ~ "Bacteroides uniformis",
      
      TRUE ~ "Other"
    ),
    plot_group = factor(plot_group, levels = plot_group_levels)
  ) %>%
  filter(!is.na(week_interval), !is.na(trajectory2))

# 4) Sum within samples, then average within and across infants
plot_dat_traj <- contrib_plot_base %>%
  group_by(trajectory2, week_interval, ID, Sample, plot_group) %>%
  summarise(
    sample_contrib = sum(pct_contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(trajectory2, week_interval, ID, plot_group) %>%
  summarise(
    infant_contrib = mean(sample_contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(trajectory2, week_interval, plot_group) %>%
  summarise(
    mean_infant_contrib = mean(infant_contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_contrib_sum = 100 * mean_infant_contrib /
      sum(mean_infant_contrib, na.rm = TRUE)
  )

stopifnot(near(sum(plot_dat_traj$pct_contrib_sum), 100))

# 5) Plot using your original named colour vector
ggplot(plot_dat_traj, aes(week_interval, pct_contrib_sum, fill = plot_group)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = plot_cols, breaks = plot_group_levels) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  facet_wrap(~ trajectory2, nrow = 1) +
  labs(
    x = "Week",
    y = "% contribution",
    fill = "Taxon",
    title = paste(sig_cazy, collapse = ", ")
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 12)
  )
