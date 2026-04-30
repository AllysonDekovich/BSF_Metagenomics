library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)


# define path at the beginning
setwd('/path/to/files')


# create functions to aid in data cleaning and manipulation

# 1. Function to isolate GTDB taxonomic strings 

clean_gtdb <- function(df) {
  df %>%
    select(user_genome, classification) %>% 
    # select specific columns, ignore the rest
    separate(classification, 
             into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
             sep=";") %>%
    mutate(across(everything(), ~str_replace(.x, ".*__", ""))) %>%
    dplyr::rename(bin = user_genome)
}

# 2. Read in relative abundance (TPM) files for all replicates
# Load relative abundance files and merge all replicates into a single TPM file

load_timepoint_tpm <- function(file1, file2, file3) {

  r1 <- read_tsv(file1, show_col_types = FALSE) %>% select(Contig, Rep1 = 2)
  r2 <- read_tsv(file2, show_col_types = FALSE) %>% select(Contig, Rep2 = 2)
  r3 <- read_tsv(file3, show_col_types = FALSE) %>% select(Contig, Rep3 = 2)
  r1 %>%
    inner_join(r2, by = "Contig") %>%
    inner_join(r3, by = "Contig")
}

# 3. Taxonomic Aggregation and Normalization
# This function reconstructs the total abundance of taxonomic Orders by mapping individual sequence contigs back to their respective bins and then to GTDBTK-assigned taxonomic classification.

summarise_to_order <- function(tpm_df, c2b_path, timepoint_label) {
  # read in contig2bin mapping file (each bin is mapped to all associated contigs)
  c2b <- read_tsv(c2b_path, col_names = c("Contig", "bin"), show_col_types = FALSE) %>%
    distinct(Contig, .keep_all = TRUE)
    
  tpm_df %>%
  # join the contig2bin assignments to the bins
    inner_join(c2b, by = "Contig") %>%
    # join those bins to GTDBTK taxonomic classifications
    left_join(master_taxonomy, by = "bin") %>%
    # group by Order level to consolidate all genomic content belonging to that taxon
    group_by(order) %>%
    # for each of the three replicates, sum the TPM of all contigs assigned to each Order.
    summarise(across(c(Rep1, Rep2, Rep3), sum, na.rm = TRUE), .groups = "drop") %>%
    # create a new column that takes the average of the 3 replicates to get a single value per timepoint for each Order
    mutate(final_avg_tpm = rowMeans(across(c(Rep1, Rep2, Rep3)), na.rm = TRUE),
           Timepoint = timepoint_label) %>%
    select(order, Timepoint, final_avg_tpm)
}

# Load all taxonomic summary files from GTDBTk and use clean_gtdb() function to isolate the appropriate information.
# Loading all taxonomic summaries and merging into one master reference
T0_bac_tax <- bind_rows(read_tsv('./gtdbtk.bac120_DM_C_T0.summary.tsv')) %>% clean_gtdb()
T0_arc_tax <- bind_rows(read_tsv('././gtdbtk.ar53_DM_C_T0.summary.tsv')) %>% clean_gtdb()
T3_tax <- read_tsv('./gtdbtk.bac120_DM_T3.summary.tsv') %>% clean_gtdb()
T7_tax <- read_tsv('./DM_T7_bac120.summary.tsv') %>% clean_gtdb()

# merge all taxonomic summaries across timepoints into a single file
master_taxonomy <- bind_rows(T0_bac_tax, T0_arc_tax, T3_tax, T7_tax)

# Use the load_timepoint_tpm() function to create a single, merged TPM file for each timepoint
tpm_T0 <- load_timepoint_tpm(
  "./rel_abundance_tpm/DM_C_T0/DM_C_T0_r1_Owings_S13_L001_tpm.tsv",
  "./rel_abundance_tpm/DM_C_T0/DM_C_T0_r2_Owings_S14_L001_tpm.tsv",
  "./rel_abundance_tpm/DM_C_T0/DM_C_T0_r3_Owings_S15_L001_tpm.tsv"
)

tpm_T3 <- load_timepoint_tpm(
  "./rel_abundance_tpm/DM_T3/DM_T3_r1_Owings_S16_L001_tpm.tsv",
  "./rel_abundance_tpm/DM_T3/DM_T3_r3_Owings_S17_L001_tpm.tsv",
  "./rel_abundance_tpm/DM_T3/DM_T3_r4_Owings_S18_L001_tpm.tsv"
)

tpm_T7 <- load_timepoint_tpm(
  "./rel_abundance_tpm/DM_T7/DM_T7_r1_Owings_S19_L001_tpm.tsv",
  "./rel_abundance_tpm/DM_T7/DM_T7_r2_Owings_S20_L001_tpm.tsv",
  "./rel_abundance_tpm/DM_T7/DM_T7_r4_Owings_S21_L001_tpm.tsv"
)


# Use the summarise_to_order() function to to get a single TPM value per timepoint for all identified orders 

master_long <- bind_rows(
  summarise_to_order(tpm_T0, "./rel_abundance_tpm/DM_C_T0/DM_C_T0_metabat2_contig2bin.tsv", "T0"),
  summarise_to_order(tpm_T3, "./rel_abundance_tpm/DM_T3/DM_T3_metabat_contig2bin.tsv", "T3"),
  summarise_to_order(tpm_T7, "./rel_abundance_tpm/DM_T7/DM_T7_metabat2_contig2bin.tsv", "T7")
) %>% 
# remove all Orders that are unknown or NA to clean up the final figure
  filter(!is.na(order), order != "", !str_detect(order, "UB|SZ"))


# Heatmap Construction

# manipulate the data to be conducive for heatmap format, where Timepoints (T0, T3, T7) are now column headers and each Order is now a single entry
# for any 'empty' values (i.e., an Order present in all timepoints except one) will be filled with a '0' instead of NA
heatmap_data <- master_long %>%
  pivot_wider(names_from = Timepoint, values_from = final_avg_tpm, values_fill = 0)

# convert heatmap to matrix form
p_matrix <- heatmap_data %>% column_to_rownames("order") %>% as.matrix()

# log transformation for the color heatmap
p_matrix_log <- log10(p_matrix + 1)

# Calculate the maximum log-abundance in the data to anchor the color scale's upper limit
max_log <- max(p_matrix_log, na.rm = TRUE)
# Fallback to prevent heatmap function from crashing; assigns max log to 1 if data is empty or all zeros
if(is.infinite(max_log) | is.na(max_log)) max_log <- 1

# create a specialized color palette, 0 (white) fir absence, 1 (Yellow) for low abundance, and max_log (Teal) for the highest observed abundance
col_fun <- colorRamp2(c(0, 1, max_log), c("#FFFFFF", "#F3C010", "#006D77"))

#Run the Heatmap
Heatmap(p_matrix_log,
        name = "Log10(TPM)", 
        col = col_fun,
        row_names_side = "left", 
        cluster_columns = FALSE, 
        column_order = c("T0", "T3", "T7"), 
        cluster_rows = TRUE,
        row_names_gp = gpar(fontsize = 9),
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- p_matrix[i, j]
          # 1. Determine the label text
          label <- if(val == 0) "0.0" else sprintf("%.1f", val)
          # 2. Determine the text color: 
          # Black only for exactly 0, White for everything else 
          txt_col <- if(val == 0) "black" else "white"
          # 3. Draw the text
          grid.text(label, x, y, gp = gpar(fontsize = 7, col = txt_col))
        })