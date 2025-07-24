##### Set up
# Set working directory: setwd("/Users/hannahbazin/Desktop/Cambridge/Academics/Han_Lab/MPhil/mphil-project")
library(tidyverse)

#### Make list of in vitro validated (or similar for aducanumab and lecanemab) drugs and targets
# Define file paths
files <- c(
  "results/humanPVATsn/network_analysis/proximity_significant_drugs/significant_targets_whole_1000_step1.csv",
  "results/humanPVATsn/network_analysis/proximity_significant_drugs/significant_targets_whole_1000_step2.csv",
  "results/humanPVATsn/network_analysis/proximity_significant_drugs/significant_targets_whole_1000_step3.csv",
  "results/humanPVATsn/network_analysis/proximity_significant_drugs/significant_targets_whole_1000_full_diff.csv"
)

# Read and combine all files
all_data <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))

# Define the list of drugs to keep (lowercase for matching)
drug_list <- c(
  "dextromethorphan", "acyclovir", "minocycline", "clevidipine", "sacubitril",
  "fostamatinib", "phenserine", "somatotropin", "aducanumab", "lecanemab"
)

# Filter for selected drugs
filtered <- all_data[all_data$Drug_Name %in% drug_list, ]

# Remove duplicate rows
filtered_unique <- unique(filtered)

# Write to file
write.csv(
  filtered_unique,
  "results/humanPVATsn/network_analysis/proximity_significant_drugs/sign_targets_whole_1000_all_steps.csv",
  row.names = FALSE
)

### Save special version for cytoscape
write.csv(data.frame(Gene = unique(filtered_unique$Drug_Target), Validation_Target = 1), "results/humanPVATsn/cytoscape/validation_targets.csv", row.names = FALSE)


### See overlap between validation targets and DEGs
# Load the two files
val <- read.csv("results/humanPVATsn/cytoscape/validation_targets.csv", stringsAsFactors = FALSE)
deg <- read.csv("results/humanPVATsn/cytoscape/combined_deg.csv", stringsAsFactors = FALSE)

# Get the number of overlapping genes
overlap_count <- length(intersect(val$Gene, deg$Gene))
print(overlap_count)



