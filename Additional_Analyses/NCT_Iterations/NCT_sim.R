# jonashaslbeck@gmail.com; Feb 24, 2021

# --------------------------------------------------------------
# ---------- Get Iteration Number ------------------------------
# --------------------------------------------------------------

#!/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
print(iter)
iter <- as.numeric(iter)


# --------------------------------------------------------------
# ---------- Load Packages -------------------------------------
# --------------------------------------------------------------

timer_total <- proc.time()[3]

# Data generation
library(IsingSampler)
source("aux_functions.R")

# NCT
library(NetworkComparisonTest)

# Parallel
library(foreach)
library(parallel)
library(doParallel)



# Number of variables (same as in main simulation)
p <- 17
n_edges <- p*(p-1)/2

# Storage
iter_seq <- round(seq(50, 1000, length=11))
iter_seq
n_is <- length(iter_seq)
l_results <- list()


# ------ Generate Data ------

# Generate Graph
# We choose: delta=0.1, n=500

G_obj <- f_GenModel_Ising(iter = iter, 
                          p = 17, 
                          Pe = 0.2, 
                          delta = 0.3)
n <- 500

# Generate data from Ising models (copied from Ising model sim)
okay <- 0
counter <- 0
while(okay==0) {
  
  # Generate Data
  data1 <- IsingSampler(n = n, graph = G_obj$Graph1, thresholds = G_obj$thresh1, nIter = 200)
  data2 <- IsingSampler(n = n, graph = G_obj$Graph2, thresholds = G_obj$thresh2, nIter = 200)
  
  # Check variances: at least two instances of each category, in each variable, in both groups
  ind <- all(apply(data1, 2, function(x) min(c(sum(x==0), sum(x==1)))) > 1) & all(apply(data2, 2, function(x) min(c(sum(x==0), sum(x==1)))) > 1)
  if(!ind) counter <- counter + 1 else okay <- 1
  
  if(counter > 5000) stop("Ridiculous long resampling.")
  
} # end while



# ------ Run NCTs ------

cluster <- 11
cl <- makeCluster(cluster, outfile="")
registerDoParallel(cl)


out <- foreach(q = 1:11,
               .packages=c("NetworkComparisonTest", "IsingSampler"),
               .export=c("f_GenModel_Ising", "data1", "data2", "n_is"),
               .verbose=TRUE) %do% {
                 
                 m_results <- matrix(NA, n_is, 4)
                 
                 # ------ Run NCT ------
                 
                 NCT_obj <- NCT(data1 = data1,
                                data2 = data2,
                                gamma = 0.25,
                                it = iter_seq[q], 
                                test.edges = TRUE,
                                binary.data = TRUE,
                                edges = "all",
                                progressbar = FALSE)
                 
                 
                 # Compute significant parameter differences
                 alpha <- .05
                 m_pvals <- NCT_obj$einv.pval
                 m_pvals[, 1] <- as.numeric(gsub('var','', m_pvals[, 1]))
                 m_pvals[, 2] <- as.numeric(gsub('var','', m_pvals[, 2]))
                 
                 m_NCT_sgnf <- matrix(0, p, p)
                 for(j in 1:n_edges) if(m_pvals[j, 3] < alpha) m_NCT_sgnf[m_pvals[j, 1], m_pvals[j, 2]] <- m_NCT_sgnf[m_pvals[j, 2], m_pvals[j, 1]] <- 1
                 
                 EstDiffs <- (NCT_obj$nw1-NCT_obj$nw2)* m_NCT_sgnf
                 
                 # ------ Evaluate ------
                 
                 # Sensitivity & Specificity
                 true_diff <- (G_obj$G_diff!=0)*1
                 est_diff <- m_NCT_sgnf
                 sen <- mean(est_diff[true_diff==1])
                 pre <- mean(true_diff[est_diff==1])
                 
                 # Estimation error (present)
                 ee_p <- mean(abs(EstDiffs[true_diff==1] + G_obj$G_diff[true_diff==1]))
                 # Estimation error (present)
                 ee_a <- mean(abs(EstDiffs[true_diff==0] + G_obj$G_diff[true_diff==0]))
                 
                 # Save
                 m_results[q, 1] <- sen
                 m_results[q, 2] <- pre
                 m_results[q, 3] <- ee_p
                 m_results[q, 4] <- ee_a
                 
                 return(m_results)
                 
               } # end foreach


stopCluster(cl) # re-register cores


# Save
saveRDS(out, file= paste0("NCT_sim_iter", iter, ".RDS"))


print("Timing: ")
proc.time()[3] - timer_total