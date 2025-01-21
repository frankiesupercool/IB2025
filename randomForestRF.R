library(randomForest)
library(snow)
library(parallel)
# Load data
load("./data/tcga_gene_expr_and_meth.RData")

ncores <- detectCores(logical = FALSE)
cl <- makeCluster(ncores)

PAM50 <- factor(PAM50)

print("small start")
small_sample <- brca.gene.expr[1:50, ]  # Taking only first 50 samples for testing
small_response <- PAM50[1:50]
small_response <- small_response[small_response != "Normal"]  # Remove "Normal" class
small_response <- factor(small_response)
small_result <- randomForest(x = small_sample, y = small_response, ntree = 500, mtry = 10)
print("small sample done")

print("Start randomForest")
# Perform normal Random Forest runs for each model
brca.gene.expr.rf <- randomForest(x = brca.gene.expr, y = PAM50, ntree = 5000, mtry = 10)
save(brca.gene.expr.rf, file = "./results/brca.gene.expr.rf.RData")

brca.meth.full.rf <- randomForest(x = brca.meth.full, y = PAM50, ntree = 5000, mtry = 10)
save(brca.meth.full.rf, file = "./results/brca.meth.full.rf.RData")

brca.combined.rf <- randomForest(x = combinedFeatureMatrix, y = PAM50, ntree = 5000, mtry = 10)
save(brca.combined.rf, file = "./results/brca.combined.rf.RData")
print("Finished randomForest")

print("Start bootstrapping")

# Define a custom bootstrap function for random forest
bootstrap_randomForest <- function(data, response, ntree = 500, mtry = 10, nboot = 10, cluster = NULL) {
  boot_models <- vector("list", nboot)

  # Perform bootstrapping with optional parallelization
  if (!is.null(cluster)) {
    clusterExport(cluster, varlist = c("data", "response", "randomForest", "ntree", "mtry"), envir = environment())
    boot_models <- parLapply(cluster, 1:nboot, function(i) {
      set.seed(i)  # Ensure reproducibility
      boot_idx <- sample(1:nrow(data), replace = TRUE)
      boot_sample <- data[boot_idx, ]
      boot_response <- response[boot_idx]
      randomForest(x = boot_sample, y = boot_response, ntree = ntree, mtry = mtry)
    })
  } else {
    for (i in 1:nboot) {
      set.seed(i)  # Ensure reproducibility
      boot_idx <- sample(1:nrow(data), replace = TRUE)
      boot_sample <- data[boot_idx, ]
      boot_response <- response[boot_idx]
      boot_models[[i]] <- randomForest(x = boot_sample, y = boot_response, ntree = ntree, mtry = mtry)
    }
  }

  return(boot_models)
}

# Apply the bootstrap function to each dataset
brca.gene.expr.boot <- bootstrap_randomForest(data = brca.gene.expr, response = PAM50, ntree = 500, mtry = 10, nboot = 10, cluster = cl)
save(brca.gene.expr.boot, file = "./results/brca.gene.expr.rf.boot.RData")

brca.meth.full.boot <- bootstrap_randomForest(data = brca.meth.full, response = PAM50, ntree = 500, mtry = 10, nboot = 10, cluster = cl)
save(brca.meth.full.boot, file = "./results/brca.meth.full.rf.boot.RData")

brca.combined.boot <- bootstrap_randomForest(data = combinedFeatureMatrix, response = PAM50, ntree = 500, mtry = 10, nboot = 10, cluster = cl)
save(brca.combined.boot, file = "./results/brca.combined.rf.boot.RData")

print("Finished bootstrapping")


# Stop the cluster after processing
stopCluster(cl)

# Evaluate and print results for each random forest model
# print("BRCA Gene Expression Model Results:")
# print(brca.gene.expr.rf)
print(paste("Out-of-bag error rate:", brca.gene.expr.rf$err.rate[nrow(brca.gene.expr.rf$err.rate), "OOB"]))
print("Variable importance:")
print(head(brca.gene.expr.rf$importance))

print("BRCA Methylation Model Results:")
print(brca.meth.full.rf)
print(paste("Out-of-bag error rate:", brca.meth.full.rf$err.rate[nrow(brca.meth.full.rf$err.rate), "OOB"]))
print("Variable importance:")
print(head(brca.meth.full.rf$importance))

print("BRCA Combined Model Results:")
print(brca.combined.rf)
print(paste("Out-of-bag error rate:", brca.combined.rf$err.rate[nrow(brca.combined.rf$err.rate), "OOB"]))
print("Variable importance:")
#print(brca.combined.rf$importance)
print(head(brca.combined.rf$importance))
