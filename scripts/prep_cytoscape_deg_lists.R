##### Set up
# Set working directory: setwd("/Users/hannahbazin/Desktop/Cambridge/Academics/Han_Lab/MPhil/mphil-project")
library(tidyverse)

############# STEP 1 -----

##### Complete DEG list
# Load 1-column CSV file
step1 <- read.csv("data/networks/step1_deg.csv", header = FALSE, col.names = "Gene")

# Add a dummy column for cytoscape import
step1$Dummy <- 1

# Save to new CSV
write.csv(step1, "results/humanPVATsn/cytoscape/step1_deg_with_dummy.csv", row.names = FALSE)

cat("Saved as step1_deg_with_dummy.csv\n")


##### Key gene list
step1_key <- read.csv("results/humanPVATsn/network_analysis/step1_key_genes_only.csv", header = FALSE, col.names = "Gene")
step1_key$Key <- 1
write.csv(step1_key, "results/humanPVATsn/cytoscape/step1_key_genes_only_with_dummy.csv", row.names = FALSE)

##### Pathway associations
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)

# Gene–pathway table
gene_pathway <- tribble(
  ~Pathway, ~Gene,
  # Notch signaling pathway
  "Notch signaling pathway", "ADAM10",
  "Notch signaling pathway", "JAG1",
  "Notch signaling pathway", "MAML3",
  "Notch signaling pathway", "MAML2",
  
  # substrate adhesion-dependent cell spreading
  "substrate adhesion-dependent cell spreading", "RHOA",
  "substrate adhesion-dependent cell spreading", "FN1",
  "substrate adhesion-dependent cell spreading", "FERMT2",
  "substrate adhesion-dependent cell spreading", "SRGAP2",
  "substrate adhesion-dependent cell spreading", "RADIL",
  "substrate adhesion-dependent cell spreading", "ANTXR1",
  "substrate adhesion-dependent cell spreading", "LAMB1",
  "substrate adhesion-dependent cell spreading", "LAMC1",
  
  # negative regulation of cell migration
  "negative regulation of cell migration", "JAG1",
  "negative regulation of cell migration", "DPYSL3",
  "negative regulation of cell migration", "NEDD9",
  "negative regulation of cell migration", "PTPRJ",
  "negative regulation of cell migration", "ROBO1",
  "negative regulation of cell migration", "TGFBR1",
  "negative regulation of cell migration", "TGFBR3",
  "negative regulation of cell migration", "TPM1",
  "negative regulation of cell migration", "VCL",
  "negative regulation of cell migration", "RECK",
  "negative regulation of cell migration", "LIMCH1",
  "negative regulation of cell migration", "KANK1",
  "negative regulation of cell migration", "SULF1",
  "negative regulation of cell migration", "EMILIN2",
  "negative regulation of cell migration", "CDH11",
  "negative regulation of cell migration", "FOXO3",
  "negative regulation of cell migration", "MITF",
  "negative regulation of cell migration", "PTEN",
  "negative regulation of cell migration", "PTPRK",
  "negative regulation of cell migration", "SLIT2",
  "negative regulation of cell migration", "DLC1",
  "negative regulation of cell migration", "NAV3",
  "negative regulation of cell migration", "PODN",
  
  # insulin-like growth factor receptor signaling pathway
  "insulin-like growth factor receptor signaling pathway", "PLCB1",
  "insulin-like growth factor receptor signaling pathway", "GHR",
  "insulin-like growth factor receptor signaling pathway", "IGF1",
  
  # cytoskeleton organization
  "cytoskeleton organization", "DST",
  "cytoskeleton organization", "CAPZB",
  "cytoskeleton organization", "DIAPH1",
  "cytoskeleton organization", "TPM1",
  "cytoskeleton organization", "LIMD1",
  "cytoskeleton organization", "CDC42BPB",
  "cytoskeleton organization", "SIPA1L3",
  "cytoskeleton organization", "FGD5",
  
  # negative regulation of canonical Wnt signaling pathway
  "negative regulation of canonical Wnt signaling pathway", "GLI3",
  "negative regulation of canonical Wnt signaling pathway", "IGFBP6",
  "negative regulation of canonical Wnt signaling pathway", "STK3",
  "negative regulation of canonical Wnt signaling pathway", "LIMD1",
  "negative regulation of canonical Wnt signaling pathway", "SOX13",
  "negative regulation of canonical Wnt signaling pathway", "KREMEN1",
  "negative regulation of canonical Wnt signaling pathway", "SHISA6",
  "negative regulation of canonical Wnt signaling pathway", "IGFBP4",
  "negative regulation of canonical Wnt signaling pathway", "BICC1",
  "negative regulation of canonical Wnt signaling pathway", "PRICKLE1",
  
  # semaphorin-plexin signaling pathway
  "semaphorin-plexin signaling pathway", "RHOA",
  "semaphorin-plexin signaling pathway", "PLXNA2",
  "semaphorin-plexin signaling pathway", "SEMA3C",
  "semaphorin-plexin signaling pathway", "NRP1",
  "semaphorin-plexin signaling pathway", "SEMA3A",
  
  # regulation of transcription by RNA polymerase II
  "regulation of transcription by RNA polymerase II", "ANXA4",
  "regulation of transcription by RNA polymerase II", "ATP2B4",
  "regulation of transcription by RNA polymerase II", "TRAK1",
  "regulation of transcription by RNA polymerase II", "CAMK2D",
  "regulation of transcription by RNA polymerase II", "FOXO3",
  "regulation of transcription by RNA polymerase II", "FOS",
  "regulation of transcription by RNA polymerase II", "JUN",
  "regulation of transcription by RNA polymerase II", "PPARG",
  
  # positive regulation of ERK1 and ERK2 cascade
  "positive regulation of ERK1 and ERK2 cascade", "APP",
  "positive regulation of ERK1 and ERK2 cascade", "CD44",
  "positive regulation of ERK1 and ERK2 cascade", "FERMT2",
  "positive regulation of ERK1 and ERK2 cascade", "ACKR3",
  "positive regulation of ERK1 and ERK2 cascade", "GLIPR2",
  "positive regulation of ERK1 and ERK2 cascade", "ANGPT1",
  "positive regulation of ERK1 and ERK2 cascade", "IGF1",
  "positive regulation of ERK1 and ERK2 cascade", "NRP1",
  "positive regulation of ERK1 and ERK2 cascade", "BMPER"
)

gene_pathway <- gene_pathway[, c("Gene", "Pathway")]

write.csv(gene_pathway, "results/humanPVATsn/cytoscape/step1_pathway_cytoscape.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

##### Combine key genes into csv
# Load necessary library
library(dplyr)

# Set working directory or use full paths
dir_path <- "results/humanPVATsn/network_analysis/"

# Define file names
files <- c("step1_key_genes_only.csv", 
           "step2_key_genes_only.csv", 
           "step3_key_genes_only.csv", 
           "full_key_genes_only.csv")

# Read and combine all files
all_key_genes <- files %>%
  lapply(function(file) {
    read.csv(file.path(dir_path, file), header = FALSE, col.names = "Gene")
  }) %>%
  bind_rows() %>%
  distinct() %>%                # Remove duplicates
  mutate(Key = 1)              # Add indicator column

# Save as new CSV for Cytoscape import
write.csv(all_key_genes, "results/humanPVATsn/cytoscape/combined_key_genes.csv", row.names = FALSE)

# Confirm
cat("Combined key genes saved to: results/humanPVATsn/cytoscape/combined_key_genes.csv\n")

#### Combine DEGs into one csv
# Define file paths
files <- c("data/networks/step1_deg.csv",
           "data/networks/step2_deg.csv",
           "data/networks/step3_deg.csv",
           "data/networks/full_diff_deg.csv")

# Read and combine all files
gene_list <- unique(unlist(lapply(files, function(file) read.csv(file, header = FALSE)[,1])))

# Create dataframe with required structure
combined_df <- data.frame(Gene = gene_list, DEG = 1)

# Save as CSV
write.csv(combined_df,
          file = "results/humanPVATsn/cytoscape/combined_deg.csv",
          row.names = FALSE)


############# STEP 2 -----
##### Complete DEG list
step2 <- read.csv("data/networks/step2_deg.csv", header = FALSE, col.names = "Gene")
step2$Step2 <- 1
write.csv(step2, "results/humanPVATsn/cytoscape/step2_deg_with_dummy.csv", row.names = FALSE)

##### Key gene list
step2_key <- read.csv("results/humanPVATsn/network_analysis/step2_key_genes_only.csv", header = FALSE, col.names = "Gene")
step2_key$Key2 <- 1
write.csv(step2_key, "results/humanPVATsn/cytoscape/step2_key_genes_only_with_dummy.csv", row.names = FALSE)

##### Pathway associations
# These are the original associations, they were grouped for visualisation - see below
gene_pathway2 <- tribble(
  ~Pathway, ~Gene,
  "negative regulation of protein phosphorylation", "FBLN1",
  "negative regulation of protein phosphorylation", "SLIT2",
  "negative regulation of protein phosphorylation", "PID1",
  "negative regulation of protein phosphorylation", "CHP1",
  "negative regulation of protein phosphorylation", "CORO1C",
  
  "fatty acid beta-oxidation", "ABCD2",
  "fatty acid beta-oxidation", "DECR1",
  "fatty acid beta-oxidation", "EHHADH",
  "fatty acid beta-oxidation", "HADHB",
  "fatty acid beta-oxidation", "HADH",
  "fatty acid beta-oxidation", "ADIPOQ",
  
  "negative regulation of gene expression", "APP",
  "negative regulation of gene expression", "CD34",
  "negative regulation of gene expression", "LRP1",
  "negative regulation of gene expression", "RBMS3",
  "negative regulation of gene expression", "MAPT",
  "negative regulation of gene expression", "PC",
  "negative regulation of gene expression", "PPARG",
  "negative regulation of gene expression", "PICALM",
  "negative regulation of gene expression", "YAP1",
  "negative regulation of gene expression", "WWP2",
  "negative regulation of gene expression", "CCDC3",
  
  "negative regulation of macrophage derived foam cell differentiation", "ABCA1",
  "negative regulation of macrophage derived foam cell differentiation", "PPARA",
  "negative regulation of macrophage derived foam cell differentiation", "PPARG",
  "negative regulation of macrophage derived foam cell differentiation", "ADIPOQ",
  
  "positive regulation of cholesterol efflux", "LRP1",
  "positive regulation of cholesterol efflux", "ABCA8",
  "positive regulation of cholesterol efflux", "CAV1",
  "positive regulation of cholesterol efflux", "PPARG",
  "positive regulation of cholesterol efflux", "ADIPOQ",
  "positive regulation of cholesterol efflux", "EEPD1",
  
  "negative regulation of cholesterol storage", "ABCA1",
  "negative regulation of cholesterol storage", "PPARA",
  "negative regulation of cholesterol storage", "PPARG",
  
  "negative regulation of sequestering of triglyceride", "PPARA",
  "negative regulation of sequestering of triglyceride", "PPARG",
  "negative regulation of sequestering of triglyceride", "ABHD5",
  "negative regulation of sequestering of triglyceride", "PNPLA2",
  
  "peptidyl-tyrosine phosphorylation", "FGFR1",
  "peptidyl-tyrosine phosphorylation", "FYN",
  "peptidyl-tyrosine phosphorylation", "PDGFRA",
  "peptidyl-tyrosine phosphorylation", "PTK2",
  "peptidyl-tyrosine phosphorylation", "INSR",
  
  "lipid storage", "CAV1",
  "lipid storage", "CD36",
  "lipid storage", "CIDEA",
  "lipid storage", "DGAT1",
  "lipid storage", "CIDEC",
  "lipid storage", "DGAT2",
  
  "negative regulation of transforming growth factor beta receptor signaling pathway", "CIDEA",
  "negative regulation of transforming growth factor beta receptor signaling pathway", "PPARA",
  "negative regulation of transforming growth factor beta receptor signaling pathway", "PPARG",
  "negative regulation of transforming growth factor beta receptor signaling pathway", "SMURF1",
  
  "hippo signaling", "FAT4",
  "hippo signaling", "TEAD1",
  "hippo signaling", "YAP1",
  
  "positive regulation of GTPase activity", "ITGB1",
  "positive regulation of GTPase activity", "LIMS1",
  "positive regulation of GTPase activity", "BCAR3",
  "positive regulation of GTPase activity", "USP6NL",
  "positive regulation of GTPase activity", "NET1",
  "positive regulation of GTPase activity", "FERMT2",
  "positive regulation of GTPase activity", "DOCK8",
  "positive regulation of GTPase activity", "DOCK11",
  
  "negative regulation of fat cell differentiation", "RUNX1T1",
  "negative regulation of fat cell differentiation", "RORA",
  "negative regulation of fat cell differentiation", "FOXO1",
  "negative regulation of fat cell differentiation", "TRIO",
  "negative regulation of fat cell differentiation", "ADIPOQ",
  "negative regulation of fat cell differentiation", "YAP1",
  "negative regulation of fat cell differentiation", "FERMT2",
  
  "protein-containing complex assembly", "LAMC1",
  "protein-containing complex assembly", "PARD3",
  "protein-containing complex assembly", "DMD",
  "protein-containing complex assembly", "FANCC",
  "protein-containing complex assembly", "PEX14",
  "protein-containing complex assembly", "TEAD1",
  "protein-containing complex assembly", "PICALM",
  "protein-containing complex assembly", "YAP1",
  "protein-containing complex assembly", "RBPMS",
  
  "amyloid-beta clearance", "C3",
  "amyloid-beta clearance", "IGF1R",
  "amyloid-beta clearance", "LRP1",
  "amyloid-beta clearance", "INSR",
  "amyloid-beta clearance", "MME",
  
  "regulation of ventricular cardiac muscle cell action potential", "CACNA1C",
  "regulation of ventricular cardiac muscle cell action potential", "CAV1",
  "regulation of ventricular cardiac muscle cell action potential", "DLG1",
  
  "positive regulation of cold-induced thermogenesis", "IGF1R",
  "positive regulation of cold-induced thermogenesis", "EBF2",
  "positive regulation of cold-induced thermogenesis", "CAV1",
  "positive regulation of cold-induced thermogenesis", "CD36",
  "positive regulation of cold-induced thermogenesis", "DECR1",
  "positive regulation of cold-induced thermogenesis", "ESRRG",
  "positive regulation of cold-induced thermogenesis", "FABP4",
  "positive regulation of cold-induced thermogenesis", "ACSL1",
  "positive regulation of cold-induced thermogenesis", "HADH",
  "positive regulation of cold-induced thermogenesis", "ADIPOQ",
  "positive regulation of cold-induced thermogenesis", "LPIN1",
  "positive regulation of cold-induced thermogenesis", "PDGFC",
  "positive regulation of cold-induced thermogenesis", "ADIPOR2",
  "positive regulation of cold-induced thermogenesis", "PPARGC1B",
  
  "positive regulation of glycolytic process", "APP",
  "positive regulation of glycolytic process", "ZBTB20",
  "positive regulation of glycolytic process", "INSR",
  "positive regulation of glycolytic process", "SLC4A4",
  "positive regulation of glycolytic process", "MLXIPL",
  
  "protein phosphorylation", "APP",
  "protein phosphorylation", "CAMK2D",
  "protein phosphorylation", "PRKG1",
  "protein phosphorylation", "DCLK1",
  "protein phosphorylation", "AKT3",
  "protein phosphorylation", "INSR",
  "protein phosphorylation", "IRAK2",
  "protein phosphorylation", "LIPE",
  "protein phosphorylation", "PRKCH",
  "protein phosphorylation", "PKN2",
  "protein phosphorylation", "MAPK10",
  "protein phosphorylation", "STK38L",
  "protein phosphorylation", "SIK2",
  "protein phosphorylation", "DAPK2",
  "protein phosphorylation", "SGK3",
  "protein phosphorylation", "HIPK2",
  "protein phosphorylation", "SNRK",
  "protein phosphorylation", "FNIP2",
  
  "positive regulation of cell migration", "IGF1R",
  "positive regulation of cell migration", "PDGFRA",
  "positive regulation of cell migration", "PPP3CA",
  "positive regulation of cell migration", "PTK2",
  "positive regulation of cell migration", "AKT2",
  "positive regulation of cell migration", "CAV1",
  "positive regulation of cell migration", "INSR",
  "positive regulation of cell migration", "ITGB1",
  "positive regulation of cell migration", "MYO1C",
  "positive regulation of cell migration", "TWIST1",
  "positive regulation of cell migration", "FERMT2",
  "positive regulation of cell migration", "SH3RF2",
  
  "positive regulation of fatty acid beta-oxidation", "AKT2",
  "positive regulation of fatty acid beta-oxidation", "ABCD2",
  "positive regulation of fatty acid beta-oxidation", "PPARA",
  "positive regulation of fatty acid beta-oxidation", "TWIST1",
  "positive regulation of fatty acid beta-oxidation", "IRS2",
  "positive regulation of fatty acid beta-oxidation", "PLIN5"
)

# These are grouped pathways for visualisation
gene_pathway2_grouped <- tribble(
  ~Pathway, ~Gene,
  
  # Group 1 – Lipid metabolism & storage
  "Lipid metabolism & storage", "ABCD2",
  "Lipid metabolism & storage", "DECR1",
  "Lipid metabolism & storage", "EHHADH",
  "Lipid metabolism & storage", "HADHB",
  "Lipid metabolism & storage", "HADH",
  "Lipid metabolism & storage", "ADIPOQ",
  "Lipid metabolism & storage", "AKT2",
  "Lipid metabolism & storage", "PPARA",
  "Lipid metabolism & storage", "TWIST1",
  "Lipid metabolism & storage", "IRS2",
  "Lipid metabolism & storage", "PLIN5",
  "Lipid metabolism & storage", "CAV1",
  "Lipid metabolism & storage", "CD36",
  "Lipid metabolism & storage", "CIDEA",
  "Lipid metabolism & storage", "DGAT1",
  "Lipid metabolism & storage", "CIDEC",
  "Lipid metabolism & storage", "DGAT2",
  "Lipid metabolism & storage", "PPARG",
  "Lipid metabolism & storage", "ABHD5",
  "Lipid metabolism & storage", "PNPLA2",
  "Lipid metabolism & storage", "ABCA1",
  "Lipid metabolism & storage", "ABCA8",
  "Lipid metabolism & storage", "EEPD1",
  "Lipid metabolism & storage", "LRP1",
  
  # Group 2 – Adipocyte differentiation & thermogenesis
  "Adipocyte differentiation & thermogenesis", "IGF1R",
  "Adipocyte differentiation & thermogenesis", "EBF2",
  "Adipocyte differentiation & thermogenesis", "CAV1",
  "Adipocyte differentiation & thermogenesis", "CD36",
  "Adipocyte differentiation & thermogenesis", "DECR1",
  "Adipocyte differentiation & thermogenesis", "ESRRG",
  "Adipocyte differentiation & thermogenesis", "FABP4",
  "Adipocyte differentiation & thermogenesis", "ACSL1",
  "Adipocyte differentiation & thermogenesis", "HADH",
  "Adipocyte differentiation & thermogenesis", "ADIPOQ",
  "Adipocyte differentiation & thermogenesis", "LPIN1",
  "Adipocyte differentiation & thermogenesis", "PDGFC",
  "Adipocyte differentiation & thermogenesis", "ADIPOR2",
  "Adipocyte differentiation & thermogenesis", "PPARGC1B",
  "Adipocyte differentiation & thermogenesis", "RUNX1T1",
  "Adipocyte differentiation & thermogenesis", "RORA",
  "Adipocyte differentiation & thermogenesis", "FOXO1",
  "Adipocyte differentiation & thermogenesis", "TRIO",
  "Adipocyte differentiation & thermogenesis", "YAP1",
  "Adipocyte differentiation & thermogenesis", "FERMT2",
  "Adipocyte differentiation & thermogenesis", "CIDEA",
  "Adipocyte differentiation & thermogenesis", "PPARA",
  "Adipocyte differentiation & thermogenesis", "PPARG",
  "Adipocyte differentiation & thermogenesis", "SMURF1",
  
  # Group 3 – Nuclear & transcriptional regulation
  "Nuclear & transcriptional regulation", "APP",
  "Nuclear & transcriptional regulation", "CD34",
  "Nuclear & transcriptional regulation", "LRP1",
  "Nuclear & transcriptional regulation", "RBMS3",
  "Nuclear & transcriptional regulation", "MAPT",
  "Nuclear & transcriptional regulation", "PC",
  "Nuclear & transcriptional regulation", "PPARG",
  "Nuclear & transcriptional regulation", "PICALM",
  "Nuclear & transcriptional regulation", "YAP1",
  "Nuclear & transcriptional regulation", "WWP2",
  "Nuclear & transcriptional regulation", "CCDC3",
  "Nuclear & transcriptional regulation", "LAMC1",
  "Nuclear & transcriptional regulation", "PARD3",
  "Nuclear & transcriptional regulation", "DMD",
  "Nuclear & transcriptional regulation", "FANCC",
  "Nuclear & transcriptional regulation", "PEX14",
  "Nuclear & transcriptional regulation", "TEAD1",
  "Nuclear & transcriptional regulation", "RBPMS",
  "Nuclear & transcriptional regulation", "FAT4",
  
  # Group 4 – Signal transduction & phosphorylation
  "Signal transduction & phosphorylation", "APP",
  "Signal transduction & phosphorylation", "CAMK2D",
  "Signal transduction & phosphorylation", "PRKG1",
  "Signal transduction & phosphorylation", "DCLK1",
  "Signal transduction & phosphorylation", "AKT3",
  "Signal transduction & phosphorylation", "INSR",
  "Signal transduction & phosphorylation", "IRAK2",
  "Signal transduction & phosphorylation", "LIPE",
  "Signal transduction & phosphorylation", "PRKCH",
  "Signal transduction & phosphorylation", "PKN2",
  "Signal transduction & phosphorylation", "MAPK10",
  "Signal transduction & phosphorylation", "STK38L",
  "Signal transduction & phosphorylation", "SIK2",
  "Signal transduction & phosphorylation", "DAPK2",
  "Signal transduction & phosphorylation", "SGK3",
  "Signal transduction & phosphorylation", "HIPK2",
  "Signal transduction & phosphorylation", "SNRK",
  "Signal transduction & phosphorylation", "FNIP2",
  "Signal transduction & phosphorylation", "FBLN1",
  "Signal transduction & phosphorylation", "SLIT2",
  "Signal transduction & phosphorylation", "PID1",
  "Signal transduction & phosphorylation", "CHP1",
  "Signal transduction & phosphorylation", "CORO1C",
  "Signal transduction & phosphorylation", "FGFR1",
  "Signal transduction & phosphorylation", "FYN",
  "Signal transduction & phosphorylation", "PDGFRA",
  "Signal transduction & phosphorylation", "PTK2",
  "Signal transduction & phosphorylation", "ZBTB20",
  "Signal transduction & phosphorylation", "SLC4A4",
  "Signal transduction & phosphorylation", "MLXIPL",
  
  # Group 5 – Migration & cytoskeleton
  "Migration & cytoskeleton", "IGF1R",
  "Migration & cytoskeleton", "PDGFRA",
  "Migration & cytoskeleton", "PPP3CA",
  "Migration & cytoskeleton", "PTK2",
  "Migration & cytoskeleton", "AKT2",
  "Migration & cytoskeleton", "CAV1",
  "Migration & cytoskeleton", "INSR",
  "Migration & cytoskeleton", "ITGB1",
  "Migration & cytoskeleton", "MYO1C",
  "Migration & cytoskeleton", "TWIST1",
  "Migration & cytoskeleton", "FERMT2",
  "Migration & cytoskeleton", "SH3RF2",
  "Migration & cytoskeleton", "LIMS1",
  "Migration & cytoskeleton", "BCAR3",
  "Migration & cytoskeleton", "USP6NL",
  "Migration & cytoskeleton", "NET1",
  "Migration & cytoskeleton", "DOCK8",
  "Migration & cytoskeleton", "DOCK11",
  "Migration & cytoskeleton", "CACNA1C",
  "Migration & cytoskeleton", "DLG1",
  
  # Group 6 – Lipid–immune cross-talk
  "Lipid–immune cross-talk", "ABCA1",
  "Lipid–immune cross-talk", "PPARA",
  "Lipid–immune cross-talk", "PPARG",
  "Lipid–immune cross-talk", "ADIPOQ",
  "Lipid–immune cross-talk", "C3",
  "Lipid–immune cross-talk", "IGF1R",
  "Lipid–immune cross-talk", "INSR",
  "Lipid–immune cross-talk", "LRP1",
  "Lipid–immune cross-talk", "MME"
)

gene_pathway2_grouped <- gene_pathway2_grouped[, c("Gene", "Pathway")]

# Save to file
write.csv(gene_pathway2_grouped, "results/humanPVATsn/cytoscape/step2_pathway_cytoscape.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")


####### STEP 3 -----
##### Complete DEG list
step3 <- read.csv("data/networks/step3_deg.csv", header = FALSE, col.names = "Gene")
step3$Step3 <- 1
write.csv(step3, "results/humanPVATsn/cytoscape/step3_deg_with_dummy.csv", row.names = FALSE)

##### Key gene list
step3_key <- read.csv("results/humanPVATsn/network_analysis/step3_key_genes_only.csv", header = FALSE, col.names = "Gene")
step3_key$Key3 <- 1
write.csv(step3_key, "results/humanPVATsn/cytoscape/step3_key_genes_only_with_dummy.csv", row.names = FALSE)

##### Pathway associations
# These are the original associations, they were grouped for visualisation - see below
gene_pathway_step3 <- tribble(
  ~Pathway, ~Gene,
  
  # skeletal system development
  "skeletal system development", "CDH11",
  "skeletal system development", "COL1A1",
  "skeletal system development", "COL1A2",
  "skeletal system development", "FBN1",
  "skeletal system development", "FGFR1",
  "skeletal system development", "GLI2",
  "skeletal system development", "IGF1",
  "skeletal system development", "RPS6KA3",
  "skeletal system development", "SH3PXD2B",
  
  # actin cytoskeleton organization
  "actin cytoskeleton organization", "FGF7",
  "actin cytoskeleton organization", "PRKG1",
  "actin cytoskeleton organization", "TNXB",
  "actin cytoskeleton organization", "EZR",
  "actin cytoskeleton organization", "ELMO1",
  "actin cytoskeleton organization", "ANTXR1",
  "actin cytoskeleton organization", "FGD5",
  "actin cytoskeleton organization", "SPTBN1",
  
  # positive regulation of ERK1 and ERK2 cascade
  "positive regulation of ERK1 and ERK2 cascade", "APP",
  "positive regulation of ERK1 and ERK2 cascade", "CD44",
  "positive regulation of ERK1 and ERK2 cascade", "IGF1",
  "positive regulation of ERK1 and ERK2 cascade", "PDGFRA",
  "positive regulation of ERK1 and ERK2 cascade", "ACKR3",
  "positive regulation of ERK1 and ERK2 cascade", "BMPER",
  "positive regulation of ERK1 and ERK2 cascade", "ADRA1A",
  "positive regulation of ERK1 and ERK2 cascade", "CD36",
  
  # collagen fibril organization
  "collagen fibril organization", "COL1A1",
  "collagen fibril organization", "COL1A2",
  "collagen fibril organization", "COL3A1",
  "collagen fibril organization", "COL5A1",
  "collagen fibril organization", "PXDN",
  
  # lncRNA-mediated post-transcriptional gene silencing
  "lncRNA-mediated post-transcriptional gene silencing", "NEAT1",
  "lncRNA-mediated post-transcriptional gene silencing", "MALAT1",
  
  # protein polyubiquitination
  "protein polyubiquitination", "CTNNB1",
  "protein polyubiquitination", "DDB2",
  "protein polyubiquitination", "FBXL7",
  "protein polyubiquitination", "SASH1",
  "protein polyubiquitination", "BCL2",
  "protein polyubiquitination", "LNPEP",
  
  # neuron projection development
  "neuron projection development", "APP",
  "neuron projection development", "NEDD4",
  "neuron projection development", "OPHN1",
  "neuron projection development", "CNTN4",
  "neuron projection development", "MAPT",
  
  # heart development
  "heart development", "CACNA1C",
  "heart development", "COL3A1",
  "heart development", "FBN1",
  "heart development", "GLI2",
  "heart development", "RBPJ",
  "heart development", "AKAP13",
  "heart development", "SH3PXD2B",
  "heart development", "PCSK5",
  "heart development", "PTEN",
  
  # receptor internalization
  "receptor internalization", "NEDD4",
  "receptor internalization", "ACKR3",
  "receptor internalization", "CD36",
  "receptor internalization", "ITGB1",
  "receptor internalization", "PICALM",
  
  # cellular response to fibroblast growth factor stimulus
  "cellular response to fibroblast growth factor stimulus", "ZFP36L1",
  "cellular response to fibroblast growth factor stimulus", "ZFP36L2",
  "cellular response to fibroblast growth factor stimulus", "CD44",
  "cellular response to fibroblast growth factor stimulus", "FGFR1",
  "cellular response to fibroblast growth factor stimulus", "NR4A1",
  
  # cytoskeleton organization
  "cytoskeleton organization", "DST",
  "cytoskeleton organization", "DPYSL2",
  "cytoskeleton organization", "TPM1",
  "cytoskeleton organization", "LIMD1",
  "cytoskeleton organization", "SIPA1L3",
  "cytoskeleton organization", "FGD5",
  "cytoskeleton organization", "SH3KBP1",
  "cytoskeleton organization", "SH3D19"
)

# These are grouped pathways for visualisation
gene_pathway3_grouped <- tribble(
  ~Pathway, ~Gene,
  "Structural remodelling", "ANTXR1",
  "Structural remodelling", "COL1A1",
  "Structural remodelling", "COL1A2",
  "Structural remodelling", "COL3A1",
  "Structural remodelling", "COL5A1",
  "Structural remodelling", "DPYSL2",
  "Structural remodelling", "DST",
  "Structural remodelling", "ELMO1",
  "Structural remodelling", "EZR",
  "Structural remodelling", "FGD5",
  "Structural remodelling", "FGF7",
  "Structural remodelling", "LIMD1",
  "Structural remodelling", "PRKG1",
  "Structural remodelling", "PXDN",
  "Structural remodelling", "SH3D19",
  "Structural remodelling", "SH3KBP1",
  "Structural remodelling", "SIPA1L3",
  "Structural remodelling", "SPTBN1",
  "Structural remodelling", "TNXB",
  "Structural remodelling", "TPM1",
  
  "Developmental signalling", "AKAP13",
  "Developmental signalling", "CACNA1C",
  "Developmental signalling", "CD44",
  "Developmental signalling", "CDH11",
  "Developmental signalling", "COL1A1",
  "Developmental signalling", "COL1A2",
  "Developmental signalling", "COL3A1",
  "Developmental signalling", "FBN1",
  "Developmental signalling", "FGFR1",
  "Developmental signalling", "GLI2",
  "Developmental signalling", "IGF1",
  "Developmental signalling", "NR4A1",
  "Developmental signalling", "PCSK5",
  "Developmental signalling", "PTEN",
  "Developmental signalling", "RBPJ",
  "Developmental signalling", "RPS6KA3",
  "Developmental signalling", "SH3PXD2B",
  "Developmental signalling", "ZFP36L1",
  "Developmental signalling", "ZFP36L2",
  
  "Neurogenic processes", "APP",
  "Neurogenic processes", "CNTN4",
  "Neurogenic processes", "MAPT",
  "Neurogenic processes", "NEDD4",
  "Neurogenic processes", "OPHN1",
  
  "Signal transduction and trafficking", "ACKR3",
  "Signal transduction and trafficking", "ADRA1A",
  "Signal transduction and trafficking", "APP",
  "Signal transduction and trafficking", "BMPER",
  "Signal transduction and trafficking", "CD36",
  "Signal transduction and trafficking", "CD44",
  "Signal transduction and trafficking", "IGF1",
  "Signal transduction and trafficking", "ITGB1",
  "Signal transduction and trafficking", "NEDD4",
  "Signal transduction and trafficking", "PDGFRA",
  "Signal transduction and trafficking", "PICALM",
  
  "Post-transcriptional and protein regulation", "BCL2",
  "Post-transcriptional and protein regulation", "CTNNB1",
  "Post-transcriptional and protein regulation", "DDB2",
  "Post-transcriptional and protein regulation", "FBXL7",
  "Post-transcriptional and protein regulation", "LNPEP",
  "Post-transcriptional and protein regulation", "MALAT1",
  "Post-transcriptional and protein regulation", "NEAT1",
  "Post-transcriptional and protein regulation", "SASH1"
)

gene_pathway3_grouped <- gene_pathway3_grouped[, c("Gene", "Pathway")]

# Save to file
write.csv(gene_pathway3_grouped, "results/humanPVATsn/cytoscape/step3_pathway_cytoscape.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")


########## FULL DIFF -----
##### Complete DEG list
full <- read.csv("data/networks/full_diff_deg.csv", header = FALSE, col.names = "Gene")
full$Full <- 1
write.csv(full, "results/humanPVATsn/cytoscape/full_deg_with_dummy.csv", row.names = FALSE)

##### Key gene list
full_key <- read.csv("results/humanPVATsn/network_analysis/full_key_genes_only.csv", header = FALSE, col.names = "Gene")
full_key$KeyFull <- 1
write.csv(full_key, "results/humanPVATsn/cytoscape/full_key_genes_only_with_dummy.csv", row.names = FALSE)

##### Pathway associations
# These are the original associations, they were grouped for visualisation - see below
gene_pathway_full <- tribble(
  ~Pathway, ~Gene,
  
  # positive regulation of fatty acid beta-oxidation
  "positive regulation of fatty acid beta-oxidation", "AKT2",
  "positive regulation of fatty acid beta-oxidation", "ABCD2",
  "positive regulation of fatty acid beta-oxidation", "PPARA",
  "positive regulation of fatty acid beta-oxidation", "TWIST1",
  "positive regulation of fatty acid beta-oxidation", "IRS2",
  "positive regulation of fatty acid beta-oxidation", "PLIN5",
  
  # NADP metabolic process
  "NADP metabolic process", "FMO2",
  "NADP metabolic process", "MDH1",
  "NADP metabolic process", "ME1",
  "NADP metabolic process", "PC",
  
  # caveola assembly
  "caveola assembly", "PACSIN2",
  "caveola assembly", "CAV1",
  "caveola assembly", "CAV2",
  
  # negative regulation of macrophage derived foam cell differentiation
  "negative regulation of macrophage derived foam cell differentiation", "ITGAV",
  "negative regulation of macrophage derived foam cell differentiation", "ABCA1",
  "negative regulation of macrophage derived foam cell differentiation", "PPARA",
  "negative regulation of macrophage derived foam cell differentiation", "PPARG",
  "negative regulation of macrophage derived foam cell differentiation", "ADIPOQ",
  "negative regulation of macrophage derived foam cell differentiation", "ABCA5",
  
  # Roundabout signaling pathway
  "Roundabout signaling pathway", "MYO9B",
  "Roundabout signaling pathway", "ROBO1",
  "Roundabout signaling pathway", "SLIT2",
  
  # regulation of bicellular tight junction assembly
  "regulation of bicellular tight junction assembly", "MYO1C",
  "regulation of bicellular tight junction assembly", "PRKACA",
  "regulation of bicellular tight junction assembly", "PRKCH",
  "regulation of bicellular tight junction assembly", "F11R",
  
  # fatty acid beta-oxidation using acyl-CoA dehydrogenase
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "ACADL",
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "ACADM",
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "ACADVL",
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "ETFA",
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "ETFDH",
  
  # positive regulation of phospholipase C activity
  "positive regulation of phospholipase C activity", "FGFR1",
  "positive regulation of phospholipase C activity", "PDGFRA",
  "positive regulation of phospholipase C activity", "PDGFRB",
  "positive regulation of phospholipase C activity", "ESR1",
  
  # negative regulation of epidermal growth factor receptor signaling pathway
  "negative regulation of epidermal growth factor receptor signaling pathway", "EGFR",
  "negative regulation of epidermal growth factor receptor signaling pathway", "PTPRJ",
  "negative regulation of epidermal growth factor receptor signaling pathway", "MVP",
  "negative regulation of epidermal growth factor receptor signaling pathway", "RNF115",
  "negative regulation of epidermal growth factor receptor signaling pathway", "ITGA1",
  "negative regulation of epidermal growth factor receptor signaling pathway", "PTPN3",
  "negative regulation of epidermal growth factor receptor signaling pathway", "DAB2IP",
  
  # cell motility
  "cell motility", "ACTB",
  "cell motility", "DST",
  "cell motility", "CD34",
  "cell motility", "CTTN",
  "cell motility", "PTK2",
  "cell motility", "TGFBR1",
  "cell motility", "IER2",
  "cell motility", "ELMO1",
  "cell motility", "ENPP2",
  
  # regulation of potassium ion transmembrane transporter activity
  "regulation of potassium ion transmembrane transporter activity", "NEDD4",
  "regulation of potassium ion transmembrane transporter activity", "FHL1",
  "regulation of potassium ion transmembrane transporter activity", "NEDD4L",
  
  # respiratory electron transport chain
  "respiratory electron transport chain", "ETFA",
  "respiratory electron transport chain", "ETFDH",
  "respiratory electron transport chain", "PPARGC1A",
  
  # cell migration involved in sprouting angiogenesis
  "cell migration involved in sprouting angiogenesis", "NR4A1",
  "cell migration involved in sprouting angiogenesis", "ROBO1",
  "cell migration involved in sprouting angiogenesis", "SLIT2",
  "cell migration involved in sprouting angiogenesis", "FGF2",
  "cell migration involved in sprouting angiogenesis", "GPLD1",
  "cell migration involved in sprouting angiogenesis", "VEGFA",
  "cell migration involved in sprouting angiogenesis", "NRP1",
  "cell migration involved in sprouting angiogenesis", "MIA3",
  
  # adipose tissue development
  "adipose tissue development", "VPS13B",
  "adipose tissue development", "SH3PXD2B",
  "adipose tissue development", "LRP5",
  "adipose tissue development", "NAMPT",
  "adipose tissue development", "PGRMC2",
  "adipose tissue development", "ARID5B",
  
  # negative regulation of anoikis
  "negative regulation of anoikis", "MCL1",
  "negative regulation of anoikis", "PTK2",
  "negative regulation of anoikis", "BCL2",
  "negative regulation of anoikis", "CAV1",
  "negative regulation of anoikis", "ITGB1",
  "negative regulation of anoikis", "PDK4",
  "negative regulation of anoikis", "TLE1",
  
  # wound healing
  "wound healing", "COL3A1",
  "wound healing", "SMAD3",
  "wound healing", "NF1",
  "wound healing", "PDGFRA",
  "wound healing", "PPP3CA",
  "wound healing", "TGFBR1",
  "wound healing", "TPM1",
  "wound healing", "DCBLD2",
  "wound healing", "FGF2",
  "wound healing", "FGF10",
  "wound healing", "CFLAR",
  "wound healing", "TSKU",
  "wound healing", "EPB41L4B",
  "wound healing", "MIA3",
  
  # fatty acid beta-oxidation
  "fatty acid beta-oxidation", "ACADM",
  "fatty acid beta-oxidation", "ABCD2",
  "fatty acid beta-oxidation", "DECR1",
  "fatty acid beta-oxidation", "EHHADH",
  "fatty acid beta-oxidation", "HADHA",
  "fatty acid beta-oxidation", "HADHB",
  "fatty acid beta-oxidation", "HADH",
  "fatty acid beta-oxidation", "SCP2",
  "fatty acid beta-oxidation", "ADIPOQ",
  "fatty acid beta-oxidation", "ACAD10",
  
  # insulin receptor signaling pathway
  "insulin receptor signaling pathway", "FER",
  "insulin receptor signaling pathway", "IGF1R",
  "insulin receptor signaling pathway", "AKT2",
  "insulin receptor signaling pathway", "CAV2",
  "insulin receptor signaling pathway", "FOXO1",
  "insulin receptor signaling pathway", "GPLD1",
  "insulin receptor signaling pathway", "INSR",
  "insulin receptor signaling pathway", "PDK4",
  "insulin receptor signaling pathway", "BCAR3",
  "insulin receptor signaling pathway", "IRS2",
  "insulin receptor signaling pathway", "NAMPT",
  "insulin receptor signaling pathway", "SORBS1",
  
  # regulation of small GTPase mediated signal transduction
  "regulation of small GTPase mediated signal transduction", "ABR",
  "regulation of small GTPase mediated signal transduction", "ARHGAP6",
  "regulation of small GTPase mediated signal transduction", "MYO9A",
  "regulation of small GTPase mediated signal transduction", "MYO9B",
  "regulation of small GTPase mediated signal transduction", "TRIO",
  "regulation of small GTPase mediated signal transduction", "ARHGAP29",
  "regulation of small GTPase mediated signal transduction", "ARHGEF10",
  "regulation of small GTPase mediated signal transduction", "AKAP13",
  "regulation of small GTPase mediated signal transduction", "SWAP70",
  "regulation of small GTPase mediated signal transduction", "ARHGAP26",
  "regulation of small GTPase mediated signal transduction", "DNMBP",
  "regulation of small GTPase mediated signal transduction", "ARHGEF12",
  "regulation of small GTPase mediated signal transduction", "SRGAP2",
  "regulation of small GTPase mediated signal transduction", "ARHGEF3",
  "regulation of small GTPase mediated signal transduction", "ARHGEF10L",
  "regulation of small GTPase mediated signal transduction", "PLEKHG1",
  "regulation of small GTPase mediated signal transduction", "ARHGAP21",
  "regulation of small GTPase mediated signal transduction", "FGD5",
  "regulation of small GTPase mediated signal transduction", "SPATA13",
  "regulation of small GTPase mediated signal transduction", "RASGRF2",
  "regulation of small GTPase mediated signal transduction", "ITSN1",
  "regulation of small GTPase mediated signal transduction", "KALRN",
  "regulation of small GTPase mediated signal transduction", "ARHGAP32",
  "regulation of small GTPase mediated signal transduction", "SRGAP3",
  "regulation of small GTPase mediated signal transduction", "FAM13A",
  "regulation of small GTPase mediated signal transduction", "NET1",
  "regulation of small GTPase mediated signal transduction", "DLC1",
  "regulation of small GTPase mediated signal transduction", "ARHGAP31",
  "regulation of small GTPase mediated signal transduction", "ARHGAP20",
  "regulation of small GTPase mediated signal transduction", "ARHGEF28",
  "regulation of small GTPase mediated signal transduction", "PREX2",
  "regulation of small GTPase mediated signal transduction", "FGD4",
  
  # cellular response to insulin stimulus
  "cellular response to insulin stimulus", "PIK3R1",
  "cellular response to insulin stimulus", "DENND4C",
  "cellular response to insulin stimulus", "AKT2",
  "cellular response to insulin stimulus", "FOXO1",
  "cellular response to insulin stimulus", "GPLD1",
  "cellular response to insulin stimulus", "INSR",
  "cellular response to insulin stimulus", "PCK1",
  "cellular response to insulin stimulus", "PDE3B",
  "cellular response to insulin stimulus", "PPARG",
  "cellular response to insulin stimulus", "IRS2",
  "cellular response to insulin stimulus", "ADIPOQ",
  "cellular response to insulin stimulus", "TBC1D4",
  "cellular response to insulin stimulus", "SORBS1",
  "cellular response to insulin stimulus", "SLC25A33",
  
  # positive regulation of cold-induced thermogenesis
  "positive regulation of cold-induced thermogenesis", "DYNC1H1",
  "positive regulation of cold-induced thermogenesis", "GRB10",
  "positive regulation of cold-induced thermogenesis", "IGF1R",
  "positive regulation of cold-induced thermogenesis", "IL4R",
  "positive regulation of cold-induced thermogenesis", "APPL2",
  "positive regulation of cold-induced thermogenesis", "EBF2",
  "positive regulation of cold-induced thermogenesis", "ACADL",
  "positive regulation of cold-induced thermogenesis", "CAV1",
  "positive regulation of cold-induced thermogenesis", "CD36",
  "positive regulation of cold-induced thermogenesis", "CEBPB",
  "positive regulation of cold-induced thermogenesis", "DECR1",
  "positive regulation of cold-induced thermogenesis", "ESRRG",
  "positive regulation of cold-induced thermogenesis", "FABP4",
  "positive regulation of cold-induced thermogenesis", "ACSL1",
  "positive regulation of cold-induced thermogenesis", "HADH",
  "positive regulation of cold-induced thermogenesis", "VEGFA",
  "positive regulation of cold-induced thermogenesis", "ADIPOQ",
  "positive regulation of cold-induced thermogenesis", "PPARGC1A",
  "positive regulation of cold-induced thermogenesis", "LPIN1",
  "positive regulation of cold-induced thermogenesis", "G0S2",
  "positive regulation of cold-induced thermogenesis", "PDGFC",
  "positive regulation of cold-induced thermogenesis", "ADIPOR2",
  "positive regulation of cold-induced thermogenesis", "PPARGC1B",
  
  # regulation of DNA-templated transcription
  "regulation of DNA-templated transcription", "ABL1",
  "regulation of DNA-templated transcription", "AHR",
  "regulation of DNA-templated transcription", "ZFHX3",
  "regulation of DNA-templated transcription", "ATRX",
  "regulation of DNA-templated transcription", "BACH1",
  "regulation of DNA-templated transcription", "BPTF",
  "regulation of DNA-templated transcription", "EFEMP1",
  "regulation of DNA-templated transcription", "FUS",
  "regulation of DNA-templated transcription", "GLI2",
  "regulation of DNA-templated transcription", "GLI3",
  "regulation of DNA-templated transcription", "SMAD3",
  "regulation of DNA-templated transcription", "RORA",
  "regulation of DNA-templated transcription", "TGFBR1",
  "regulation of DNA-templated transcription", "RPS6KA5",
  "regulation of DNA-templated transcription", "SOX13",
  "regulation of DNA-templated transcription", "ERC1",
  "regulation of DNA-templated transcription", "CEBPB",
  "regulation of DNA-templated transcription", "ATF2",
  "regulation of DNA-templated transcription", "ESR2",
  "regulation of DNA-templated transcription", "ESRRG",
  "regulation of DNA-templated transcription", "NR3C1",
  "regulation of DNA-templated transcription", "INSR",
  "regulation of DNA-templated transcription", "RREB1",
  "regulation of DNA-templated transcription", "SMARCA2",
  "regulation of DNA-templated transcription", "SP1",
  "regulation of DNA-templated transcription", "PPARGC1A",
  "regulation of DNA-templated transcription", "ATF7",
  "regulation of DNA-templated transcription", "JADE1",
  "regulation of DNA-templated transcription", "PPARGC1B"
)

# Define pathway-to-group mapping
pathway_to_group <- tribble(
  ~Pathway, ~Group,
  
  # Group 1: Insulin signalling
  "insulin receptor signaling pathway", "Insulin signalling",
  "cellular response to insulin stimulus", "Insulin signalling", 
  "positive regulation of phospholipase C activity", "Insulin signalling",
  
  # Group 2: Fatty acid beta-oxidation
  "positive regulation of fatty acid beta-oxidation", "Fatty acid beta-oxidation",
  "fatty acid beta-oxidation", "Fatty acid beta-oxidation",
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "Fatty acid beta-oxidation",
  
  # Group 3: Cell migration
  "regulation of bicellular tight junction assembly", "Cell migration",
  "cell motility", "Cell migration",
  "cell migration involved in sprouting angiogenesis", "Cell migration",
  "Roundabout signaling pathway", "Cell migration",
  "regulation of small GTPase mediated signal transduction", "Cell migration",
  
  # Group 4: Thermogenic processes
  "NADP metabolic process", "Thermogenic processes",
  "respiratory electron transport chain", "Thermogenic processes",
  "positive regulation of cold-induced thermogenesis", "Thermogenic processes",
  "regulation of DNA-templated transcription", "Thermogenic processes",
  
  # Group 5: Adipose tissue development
  "adipose tissue development", "Adipose tissue development",
  "negative regulation of macrophage derived foam cell differentiation", "Adipose tissue development",
  "negative regulation of anoikis", "Adipose tissue development", 
  "wound healing", "Adipose tissue development",
  
  # Group 6: Membrane organisation
  "caveola assembly", "Membrane organisation",
  "regulation of potassium ion transmembrane transporter activity", "Membrane organisation",
  
  # Group 7: Growth factor signaling
  "negative regulation of epidermal growth factor receptor signaling pathway", "Growth factor signaling"
)

# Create the gene-group mapping
gene_pathway_full_grouped <- gene_pathway_full %>%
  left_join(pathway_to_group, by = "Pathway") %>%
  select(Gene, Pathway = Group) %>%
  arrange(Pathway, Gene)

# Save to file
write.csv(gene_pathway_full_grouped, "results/humanPVATsn/cytoscape/full_diff_pathway_cytoscape.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  




















