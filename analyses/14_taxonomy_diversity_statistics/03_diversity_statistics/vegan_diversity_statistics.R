library(tidyverse)
library(vegan)
library(FSA)
library(ggpubr)


# modify the summarise_to_order() function from the taxonomy_abundance_heatmap.R script to keep the replicate TPMs separate for calculating diversity statistics
process_replicates <- function(tpm_df, c2b_path, timepoint_label) {
  c2b <- read_tsv(c2b_path, col_names = c("Contig", "bin"), show_col_types = FALSE) %>%
    distinct(Contig, .keep_all = TRUE)
  
  tpm_df %>%
    inner_join(c2b, by = "Contig") %>%
    left_join(master_taxonomy, by = "bin") %>%
    group_by(order) %>%
    summarise(across(c(Rep1, Rep2, Rep3), sum, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(cols = c(Rep1, Rep2, Rep3), names_to = "Replicate", values_to = "TPM") %>%
    mutate(Timepoint = timepoint_label, 
           SampleID = paste(timepoint_label, Replicate, sep = "_"))
}

# combine all replicates across timepoints
master_replicates <- bind_rows(
  process_replicates(tpm_T0, "./rel_abundance_tpm/DM_C_T0/DM_C_T0_metabat2_contig2bin.tsv", "T0"),
  process_replicates(tpm_T3, "./rel_abundance_tpm/DM_T3/DM_T3_metabat_contig2bin.tsv", "T3"),
  process_replicates(tpm_T7, "./rel_abundance_tpm/DM_T7/DM_T7_metabat2_contig2bin.tsv", "T7")
)

# convert the master taxonomy object (combined MAGs and abundances) to wide format for vegan
# Each Timepoint becomes a row, each Order becomes a column
diversity_replicate_wide <- master_replicates %>%
  select(SampleID, order, TPM, Timepoint) %>%
  filter(!is.na(order), order != "",!str_detect(order, "UB|SZ|Unassigned")) %>% 
  pivot_wider(names_from = order, values_from = TPM, values_fill = 0) 


# create community data matrix for vegan
replicate_matrix <- diversity_replicate_wide %>%
  column_to_rownames("SampleID") %>%
  select(-Timepoint) %>%
  as.matrix()

# create metadata for plotting
metadata <- diversity_replicate_wide %>% select(Timepoint)

# calculate alpha diversity
alpha_metrics <- data.frame(
  SampleID = rownames(replicate_matrix),
  Shannon  = diversity(replicate_matrix, index = "shannon"),
  Richness = specnumber(replicate_matrix),
  Timepoint = metadata$Timepoint
)

# Kruskal-Wallis test to determine if the median Shannon diversity differs significantly across timepoints.
kw_test <- kruskal.test(Shannon ~ Timepoint, data = alpha_metrics)
print(kw_test)

# Dunn's posthoc test to determine which group comparisons were significantly different
dunn_res <- dunnTest(Shannon ~ Timepoint, 
                     data = alpha_metrics, 
                     method = "hochberg")
print(dunn_res)

# plot the alpha diversity
my_colors <- c("T0" = "#008B8B", "T3" = "#E3B448", "T7" = "#D85A2B")

ggboxplot(alpha_metrics, x = "Timepoint", y = "Shannon",
          fill = "Timepoint", palette = my_colors,
          add = "jitter", order = c("T0", "T3", "T7")) +
  # manually add the significant comparison between T3 and T7
  geom_signif(comparisons = list(c("T3", "T7")), 
              annotations = "*", 
              y_position = max(alpha_metrics$Shannon) * 1.05, 
              tip_length = 0.03) +
  # Add the global krusal wallis value
  stat_compare_means(method = "kruskal.test", 
                     label.x = 0.6, 
                     label.y = max(alpha_metrics$Shannon) * 1.2) + 
  labs(title = "Alpha Diversity Across Timepoints",
       y = "Shannon Index (H')") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# Beta Diversity
dist_bc <- vegdist(replicate_matrix, method = "bray")

# run PERMANOVA
perm_test <- adonis2(dist_bc ~ Timepoint, data = metadata, permutations = 999)
print(perm_test)

# Run NMDS
set.seed(42) # For reproducibility
nmds_res <- metaMDS(dist_bc, k = 2, trymax = 100)

# Extract NMDS scores for plotting
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$Timepoint <- metadata$Timepoint

# Extract p-value, R2 from the PERMANOVA and create a stat label for ggpubr
perm_p  <- perm_test$`Pr(>F)`[1]
perm_R2 <- round(perm_test$R2[1], 3)
stat_label<- paste0("PERMANOVA: R² = ", round(perm_test$R2[1], 3), 
                           ", p = ", perm_test$`Pr(>F)`[1])
                           
# Plot NMDS
ggscatter(nmds_scores, x = "NMDS1", y = "NMDS2",
          color = "Timepoint", 
          fill = "Timepoint",
          palette = c("#008B8B", "#E3B448", "#D85A2B"),
          shape = 21, size = 4,
          ellipse = TRUE, 
          ellipse.type = "convex",
          mean.point = FALSE) +
  # Stress Value
  annotate("text", x = Inf, y = -Inf, 
           label = paste("Stress:", round(nmds_res$stress, 3)), 
           hjust = 1.1, vjust = -2.8, size = 3, fontface = "italic") +
  # PERMANOVA Stats
  annotate("text", x = Inf, y = -Inf, 
           label = stat_label, 
           hjust = 1.05, vjust = -1.2, size = 3) +
  labs(title = "Beta Diversity",
       subtitle = "NMDS based on Bray-Curtis Dissimilarity") +
  theme_pubr() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "right"
  )