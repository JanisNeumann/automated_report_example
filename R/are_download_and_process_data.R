#### Dependencies ####

library(GEOquery)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)



#### Settings ####

# Define GEO data set. Changing this will break this example code down the line as it assumes specific clinical variables.
# But technically, it can be changed.
GSE_id <- "GSE111629"

# Maximum download time.
options(timeout = max(99999, getOption("timeout")))

# Detection p-value threshold.
detection_pval_threshold <- 0.001



#### Ensure existence of data directory.
if (!dir.exists("data")) {
  dir.create("data")
}



#### IDAT File Processing ####

# Download and extract data.
getGEOSuppFiles(GSE_id, baseDir = "data")

untar(paste0("data/", GSE_id, "/", GSE_id, "_RAW.tar"),
  exdir = paste0("data/", GSE_id, "/idat")
)

sapply(
  list.files(paste0("data/", GSE_id, "/idat"),
    pattern = "idat.gz$", full = TRUE
  ),
  gunzip,
  overwrite = TRUE
)

# Read IDAT files.
idat_file_names <- list.files(paste0("data/", GSE_id, "/idat"),
  pattern = "Grn.idat", full.names = TRUE
)
rg_set <- read.metharray(idat_file_names, verbose = TRUE, force = TRUE)

# Calculate detection p-values.
# Then drop samples with mean or median detection p-value above threshold as well as
# CpG sites above the same value.
detection_pvals <- detectionP(rg_set)

failed_samples <- pmax(colMeans(detection_pvals), colMedians(detection_pvals)) > detection_pval_threshold

detection_pvals <- detection_pvals[, !failed_samples]

failed_sites <- pmax(rowMeans(detection_pvals), rowMedians(detection_pvals)) > detection_pval_threshold

rg_set <- rg_set[!failed_sites, !failed_samples]

rm(detection_pvals)
gc()

# Generate betas object.
betas <- t(getBeta(preprocessNoob(rg_set)))
rm(rg_set)
gc()



#### Metadata ####

# Download series matrix to get relevant data from.
GSE_series_matrix <- getGEO(GSE_id, GSEMatrix = TRUE)
saveRDS(GSE_series_matrix, paste0("data/", GSE_id, "series_matrix.rds"))

# Make clinical data frame retaining age, gender and diagnosis.
# For this simple example report, nothing else will be used.
clinical_data <- GSE_series_matrix$GSE111629_series_matrix.txt.gz@phenoData@data %>%
  dplyr::select("age:ch1", "disease state:ch1", "gender:ch1")

# Rename columns, make age a numeric variable, turn the other two into factors and add a logical Parkinson diagnosis column.
colnames(clinical_data) <- c("Age", "Disease_State", "Gender")
clinical_data$Age <- as.numeric(clinical_data$Age)
clinical_data$Disease_State <- as.factor(clinical_data$Disease_State)
clinical_data$Gender <- as.factor(clinical_data$Gender)
clinical_data <- data.frame(clinical_data,
  PD_positive = clinical_data$Disease_State == "Parkinson's disease (PD)"
)



#### Finalize & Save Files, Cleanup ####

# Filter out samples lacking either betas or clinical annotation.
rownames(betas) <- gsub("_.*", "", rownames(betas))
retained_samples <- Reduce(intersect, list(rownames(betas), rownames(clinical_data)))
betas <- betas[retained_samples, ]
clinical_data <- clinical_data[retained_samples, ]

saveRDS(betas, paste0("data/", GSE_id, "/betas.rds"))
saveRDS(clinical_data, paste0("data/", GSE_id, "/clinical_data.rds"))
rm(list = ls())
gc()
