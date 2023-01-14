#### Dependencies ####

# Packages
library(glmnet)

# Functions
source("R/are_split_data_classifier.R")



#### Preparation ####

# Run ID and output folder
run_id <- format(Sys.Date(), "%Y%m%d")
out_dir <- paste0("results/", run_id)
if (!dir.exists(out_dir)){
  dir.create(out_dir)
}

# Alphas to use.
alpha_min <- 0
alpha_max <- 1
alpha_n <- 11
alpha_vector <- seq(from = alpha_min, to = alpha_max, length.out = alpha_n)

# Load data
betas <- readRDS("data/GSE111629/betas.rds")
clinical_data <- readRDS("data/GSE111629/clinical_data.rds")

# Create training, development and holdout data split
modeling_data_split <- split_data_classifier(data = clinical_data,
                                             dev_ratio = 0.2,
                                             holdout_ratio = 0.2,
                                             y_var = "PD_positive")



#### Feature Selection ####

# Mean difference of more than 3% between Parkinson's patients and control cases.
# This is to increase stability to by increasing the chance of the model picking sites where the difference is due to biology and not just technical variance.
# (At least less than 2% is within the noise range of the arrays.)
# Of course, this is done using only samples that will not be part of the final test/holdout data.

betas_training_dev <- betas[c(rownames(modeling_data_split$train_set), rownames(modeling_data_split$dev_set)), ]

training_dev_diagnosis <- c(modeling_data_split$train_set$PD_positive,
                            modeling_data_split$dev_set$PD_positive)

colmeans_cases <- colMeans(betas_training_dev[training_dev_diagnosis,
                                              grep("^cg", colnames(betas_training_dev), value = TRUE, invert = FALSE)])
colmeans_controls <- colMeans(betas_training_dev[!training_dev_diagnosis,
                                                 grep("^cg", colnames(betas_training_dev), value = TRUE, invert = FALSE)])

colmeans_diff <- abs(colmeans_cases - colmeans_controls)

selected_cpgs <- names(colmeans_diff[colmeans_diff >= 0.03])

rm(betas_training_dev); gc()



#### Elastic Net Modeling ####

# Split betas into training, dev and holdout data while only retaining CpG sites selected in previous step.
betas_training <- betas[rownames(modeling_data_split$train_set), selected_cpgs]
betas_dev <- betas[rownames(modeling_data_split$dev_set), selected_cpgs]
betas_holdout <- betas[rownames(modeling_data_split$holdout_set), selected_cpgs]

# For each alpha to use, run cv.glmnet and glmnet and add them to lists of models.
elnet_cv_model_list <- list()
elnet_model_list <- list()

for (i in 1:length(alpha_vector)) {
  set.seed(16)
  i_cv <- cv.glmnet(betas_training, modeling_data_split$train_set$PD_positive,
                    family = "binomial", alpha = alpha_vector[i]
  )
  
  set.seed(16)
  i_model <- glmnet(betas_training, modeling_data_split$train_set$PD_positive,
                    family = "binomial", alpha = alpha_vector[i]
  )
  
  elnet_model_list[[i]] <- i_model
  elnet_cv_model_list[[i]] <- i_cv
  
  rm(i_cv, i_model)
  gc()
}
rm(i)

# Pick best alpha and lambda based on cross-validation results and dev set performance (balanced accuracy).
balanced_accuracy_vector <- c()
lambda_1se_vector <- c()

for (i in 1:length(elnet_model_list)) {
  
  lambda_1se_vector[i] <- elnet_cv_model_list[[i]]$lambda.1se
  
  i_pred <- predict(elnet_model_list[[i]], betas_dev, s = elnet_cv_model_list[[i]]$lambda.1se, type = "response")
  
  i_sensitivity <- caret::sensitivity(as.factor(round(i_pred, 0)), as.factor(modeling_data_split$dev_set$PD_positive + 0))
  i_specificity <- caret::specificity(as.factor(round(i_pred, 0)), as.factor(modeling_data_split$dev_set$PD_positive + 0))
  
  balanced_accuracy_vector[i] <- (i_sensitivity + i_specificity) / 2
  
}

alpha_best <- alpha_vector[which.max(balanced_accuracy_vector)]
lambda_best <- lambda_1se_vector[which.max(balanced_accuracy_vector)]
model_best <- elnet_model_list[[which.max(balanced_accuracy_vector)]]

print(paste0("Optimal balanced accuracy (",
             max(balanced_accuracy_vector),
             ") achieved using model ",
             which.max(balanced_accuracy_vector),
             " with an alpha of ",
             alpha_best,
             " and a lambda of ",
             lambda_best,
             "."))



#### Holdout Data Results ####

holdout_pred <- predict(model_best, betas_holdout, s = elnet_cv_model_list[[which.max(balanced_accuracy_vector)]]$lambda.1se,
                        type = "response")
holdout_votes <- round(holdout_pred, 0)
holdout_sensitivity <- caret::sensitivity(as.factor(holdout_votes), as.factor(modeling_data_split$holdout_set$PD_positive + 0))
holdout_specificity <- caret::specificity(as.factor(holdout_votes), as.factor(modeling_data_split$holdout_set$PD_positive + 0))



#### Define Report Inputs ####

# Cohort diagnosis table
report_cohort_table <- data.frame(Training = c(table(modeling_data_split$train_set$Disease_State)),
                                  Development = c(table(modeling_data_split$dev_set$Disease_State)),
                                  Holdout = c(table(modeling_data_split$holdout_set$Disease_State)))
report_cohort_table <- data.frame("Group" = rownames(report_cohort_table),
                                  report_cohort_table)

# Cohort plots for clinical data
report_clinical_data <- rbind.data.frame(data.frame(modeling_data_split$train_set,
                                                    Data = "Training"),
                                         data.frame(modeling_data_split$dev_set,
                                                    Data = "Development"),
                                         data.frame(modeling_data_split$holdout_set,
                                                    Data = "Holdout"))
report_clinical_data$Data <- factor(report_clinical_data$Data, levels = c("Training", "Development", "Holdout"))
colnames(report_clinical_data)[colnames(report_clinical_data) == "Disease_State"] <- "Group"

# PCA
report_pca <- prcomp(betas)
report_pca_imp <- data.frame(summary(report_pca)$importance)
report_pca_pd <- data.frame(V1 = report_pca$x[, 1],
                            V2 = report_pca$x[, 2],
                            Group = clinical_data$Disease_State,
                            Gender = clinical_data$Gender)

# Dev set ROC
report_roc_glmnet <- as.data.frame(roc.glmnet(model_best, betas_dev, modeling_data_split$dev_set$PD_positive,
                                              s = elnet_cv_model_list[[i]]$lambda.1se))

# Best model metrics
report_coef <- coef(model_best, s = lambda_best)
report_coef <- report_coef[report_coef[, 1] != 0, ]

report_metrics <- data.frame(Metric = c("Alpha", "Lambda", "CpGs Used",
                                        "Holdout Sensitivity", "Holdout Specificity", "Holdout Balanced Accuracy"),
                             Value = as.character(c(alpha_best, round(lambda_best, 3), length(report_coef) - 1,
                                                    round(holdout_sensitivity, 3), round(holdout_specificity, 3), round((holdout_sensitivity + holdout_specificity)/2, 3))))

# Holdout prediction plots

report_holdout_pred <- data.frame(Predicted <- holdout_pred,
                                  Observed <- modeling_data_split$holdout_set$Disease_State)



#### Save Results & Generate Report ####

saveRDS(model_best, paste0(out_dir, "/best_model.rds"))

report_input <- list(cohort_table = report_cohort_table,
                     clinical_data = report_clinical_data,
                     pca_imp = report_pca_imp,
                     pca_pd = report_pca_pd,
                     roc = report_roc_glmnet,
                     metrics = report_metrics,
                     holdout_pred = report_holdout_pred)

saveRDS(report_input,
        paste0(out_dir, "/report_input.rds"))

rmarkdown::render("R/are_analysis.Rmd",
                  knit_root_dir = getwd(),
                  output_file = paste0("../", out_dir, "/analysis_report.html"))
