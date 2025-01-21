library(snow)
library(varSelRF)
library(ranger)

#load data
load("./data/tcga_gene_expr_and_meth.RData")

#start cluster, executing in parallel
#cl <- makeMPIcluster(16) #uses 16 cores. change to number of available cores.
ncores <- detectCores(logical = FALSE)
cl <- makeCluster(ncores)

PAM50 <- factor(PAM50)

print("small start")
small_sample <- brca.gene.expr[1:50,]
small_response <- PAM50[1:50]
small_response <- small_response[small_response != "Normal"]
small_response <- factor(small_response)
#small_result <- varSelRF(small_sample, small_response, keep.forest = TRUE)
small_result <- ranger(y = small_response, x = small_sample, num.trees = 500, mtry = 10, importance = "impurity", num.threads = ncores)
print("small sample done")

#perform normal varSelRF runs for each model
#brca.gene.expr.varSelRF <- varSelRF(brca.gene.expr, PAM50, keep.forest=T)
#save(brca.gene.expr.varSelRF, file="./results/brca.gene.expr.varSelRF.RData")
print("start gene exp")
brca.gene.expr.rf <- ranger(y = PAM50, x = brca.gene.expr, num.trees = 500, mtry = 10, importance = "impurity", num.threads = ncores)
save(brca.gene.expr.rf, file = "./results/brca.gene.expr.rf.ranger.RData")
print("gene exp done")

print("start boot strapping")
# Bootstrapping on each model using the above results as seed
# Ranger does not have a direct bootstrap method like varSelRFBoot, so you can implement bootstrapping manually.
bootnumber <- 10
boot_results <- list()

for (i in 1:bootnumber) {
  # Bootstrap sampling
  boot_sample <- brca.gene.expr[sample(nrow(brca.gene.expr), replace = TRUE), ]
  boot_response <- PAM50[sample(length(PAM50), replace = TRUE)]

  # Train a new Random Forest model on the bootstrapped sample
  boot_model <- ranger(y = boot_response, x = boot_sample, num.trees = 500, mtry = 10, importance = "impurity", num.threads = ncores)
  boot_results[[i]] <- boot_model
}

# Save the bootstrapped results
save(boot_results, file = "./results/brca.gene.expr.rf.boot.ranger.RData")
print("boot strapping done")


#stop cluster
stopCluster(cl)
