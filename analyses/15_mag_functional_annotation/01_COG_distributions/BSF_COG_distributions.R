library(tidyverse)


# set working directory prior to running analyses
setwd("/path/to/files")

# List taxonomy files from GTDBTk
tax_files <- c("./gtdbtk.bac120_DM_C_T0.summary.tsv","./gtdbtk.ar53_DM_C_T0.summary.tsv","./gtdbtk.bac120_DM_T3.summary.tsv","./DM_T7_bac120.summary.tsv")

# Define a function to read in the taxonomy files, extract Order, and create a custom label: bin (Order)
taxonomy_master <- map_df(tax_files, function(f) {
  read_tsv(f) %>%
    # Select only id, classification
    select(user_genome, classification) 
}) %>%
  # Now extract the Order and create the label
  mutate(Order = str_extract(classification, "(?<=o__)[^;]+")) %>%
  mutate(Clean_Label = paste0(user_genome, " (", Order, ")"))

# Define the entire filepath to the eggNOG .annotation files and store them in an R variable
file_path <- "path/to/annotation/files"
eggnog_files <- list.files(path = file_path, pattern = "*.annotations", full.names = TRUE)

# Custom function to read in eggNOG mapper files and extract necessary data
eggnog_master <- map_df(eggnog_files, function(f) {
  read_tsv(f, comment = "##") %>% # Skip header rows
    mutate(Bin = basename(f)) %>% # Tag each row with the filename
    mutate(Bin = str_remove(Bin, ".emapper.annotations")) # Clean filename to match user_id
})

# COG Category Lookup Table; only kept most common COGs for clarity
cog_lookup <- c(
  "C" = "Energy production", "G" = "Carbohydrate transport",
  "E" = "Amino acid transport", "F" = "Nucleotide transport",
  "H" = "Coenzyme transport", "I" = "Lipid metabolism",
  "P" = "Inorganic ion transport", "Q" = "Secondary metabolites",
  "J" = "Translation", "L" = "Replication/Repair", 
  "K" = "Transcription", "M" = "Cell wall biogenesis",
  "O" = "Post-translational mod", "T" = "Signal transduction",
  "U" = "Intracellular trafficking", "V" = "Defense mechanisms"
)

# Merge functional annotations with taxonomic metadata and normalize COG categories
# NOTE: Many genes within a MAG are assigned to multiple COG categories (e.g. "JK"), so these assignments are expanded so each functional category is represented appropriately in the bin's total profile.
plot_data <- eggnog_master %>%
  # clean the eggNOG bin column immediately so it matches 'user_genome'
  mutate(Bin = str_remove(Bin, "_eggnog")) %>% 
  # isolate the COG data and expand multi-letter COG codes into individual rows
  select(Bin, COG_category) %>% 
  separate_rows(COG_category, sep = "") %>% 
  filter(COG_category != "" & COG_category != "-") %>%
  # join with GTDBTK taxonomy
  left_join(taxonomy_master, by = c("Bin" = "user_genome")) %>%
  # aggregate counts per COG category
  count(Bin, Order, COG_category)

# Additional filtering to removing unknown or NA taxonomic orders, extract timepoints from the metadata, and rename lettered COG categories to the descriptions in the COG lookup table
plot_data_filtered <- plot_data %>%
  # remove unknown/NA bins
  mutate(Bin = str_remove(Bin, "_eggnog")) %>% 
  filter(!str_detect(Order, "^SZ|^U")) %>%
  # Extract Timepoint (T0, T3, T7)
  mutate(Timepoint = str_extract(Bin, "T[037]")) %>%
  # Apply full COG names
  mutate(COG_Name = cog_lookup[COG_category]) %>%
  filter(!is.na(COG_Name)) %>%
  filter(!is.na(Order))

# Collapse the data by Order and calculate relative abundance of each COG per Order and Tiempoint
collapsed_data <- plot_data_filtered %>%
  group_by(Order, Timepoint, COG_Name) %>%
  summarise(Total_n = sum(n), .groups = "drop") %>%
  # Calculate relative abundance per Order/Timepoint
  group_by(Order, Timepoint) %>%
  mutate(RelAbund = (Total_n / sum(Total_n)) * 100) %>%
  ungroup()

# Ensure every Order has a spot in every Timepoint (for the empty spaces)
collapsed_complete <- collapsed_data %>%
  complete(Order, Timepoint, COG_Name, fill = list(RelAbund = 0))

# Cluster Order by similar "guilds". 
# Noticed there were a few different scenarios that the MAGs could cluster by:
# Sentinel: present at the start, absent at the end
# Core: present throughout (start -> end)
# Recruit: absent at the start, present at the end
# Transient: present in T3 stage only
occupancy <- collapsed_data %>%
  group_by(Order) %>%
  summarise(
    has_T0 = any(Timepoint == "T0"),
    has_T3 = any(Timepoint == "T3"),
    has_T7 = any(Timepoint == "T7"),
    .groups = "drop"
  ) %>%
  mutate(Strategy = case_when(
    has_T0 & !has_T7           ~ "Sentinel",      # Present at start, gone at end
    has_T0 & has_T7            ~ "Core",          # Present at start and end
    !has_T0 & has_T7           ~ "Recruit",       # Appeared later
    TRUE                       ~ "Transient"
  ))

# Function to cluster MAGs by similar patterns in COG profile across timepoints
get_clustered_suborder <- function(target_strategy) {
  orders_in_strat <- occupancy %>% filter(Strategy == target_strategy) %>% pull(Order)
  if(length(orders_in_strat) == 0) return(NULL)
  if(length(orders_in_strat) <= 2) return(orders_in_strat)
  
  mat <- collapsed_complete %>%
    filter(Order %in% orders_in_strat) %>%
    group_by(Order, COG_Name) %>%
    summarise(m = mean(RelAbund), .groups = "drop") %>%
    pivot_wider(names_from = COG_Name, values_from = m, values_fill = 0) %>%
    column_to_rownames("Order")
  
  hc <- hclust(dist(mat))
  return(labels(dist(mat))[hc$order])
}

# Assign factor levels to MAGs 
sentinel_order <- get_clustered_suborder("Sentinel")
core_order     <- get_clustered_suborder("Core")
recruit_order  <- get_clustered_suborder("Recruit")
transient_order <- get_clustered_suborder("Transient")

final_order <- c(sentinel_order, core_order, recruit_order, transient_order)
collapsed_complete$Order <- factor(collapsed_complete$Order, levels = rev(final_order))

# Bubble plot
ggplot(collapsed_complete, aes(x = COG_Name, y = Order)) +
  geom_point(aes(color = RelAbund), size = 5, shape = 16) +
  facet_grid(. ~ Timepoint) + 
  scale_color_viridis_c(option = "magma", name = "Relative Abundance (%)", direction = -1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 11, face = "italic"), 
    panel.grid.major = element_line(color = "grey95"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0.6, "lines"),
    legend.position = "right",
    plot.title = element_text(hjust=0.5)
  ) +
  labs(
    x = "Functional Category (COG)",
    y = "Taxonomic Order"
  )
