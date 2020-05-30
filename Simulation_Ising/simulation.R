# jonashaslbeck@gmail.com; May 2020

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

# Identify group differences
library(NetworkComparisonTest)
library(mgm)
library(BGGM)

# Parallel
library(foreach)
library(parallel)
library(doParallel)


# --------------------------------------------------------------
# ---------- Simulation Specs ----------------------------------
# --------------------------------------------------------------

# Number of nodes = median p in network predictability paper
p_NP <- c(14, 27, 17, 18, 12, 27, 12, 28, 63, 6, 16, 17, 11, 19, 13, 17, 20, 56)
p <- median(p_NP)
n_edges <- p*(p-1)/2

# n-variations
n_seq <- round(exp(seq(3, 8.5171, length=10)))

# --------------------------------------------------------------
# ---------- Aux function: Generate True Ising Models ----------
# --------------------------------------------------------------

source("aux_functions.R")

# --------------------------------------------------------------
# ---------- Generate Data -------------------------------------
# --------------------------------------------------------------

# Replace three pathological iterations 
# These took very long too complete; this is likely due to the BGGM sampling for low n
iter_set <- iter
if(iter == 96) iter_set <- 96 + 1000
if(iter == 113) iter_set <- 113 + 1000
if(iter == 199) iter_set <- 199 + 1000
if(iter == 68) iter_set <- 68 + 1000


cluster <- 10
cl <- makeCluster(cluster, outfile="")
registerDoParallel(cl)

out <- foreach(ni = 1:10,
               .packages=c("NetworkComparisonTest", "mgm", "IsingSampler"),
               .export=c("n_seq", "iter", "f_GenModel_Ising", "n_edges", "iter", "iter_set"),
               .verbose=TRUE) %dopar% {
                 
                 # Storage
                 a_out <- array(NA, dim = c(p, p, 2, 10, 3)) 
                 # p, p, 2 results (real: difference, logical: difference), 
                 # 8 estimation method; 8 n-variations; 3 for delta theta
                 # need 2 results
                 
                 m_time <- matrix(NA, 10, 3) # computation time
                 v_counter <- rep(NA, 3) # record number of resamples needed to obtain nonzero variances
                 
                 
                 n <- n_seq[ni] # sample size of each group
                 
                 d_seq <- c(.15, .3, 0.6) # updated May 2020
                 
                 for(del in 1:3) {
                   
                   set.seed(iter_set)
                   
                   # Generate Graph
                   G_obj <- f_GenModel_Ising(iter = iter_set, 
                                             p = 17, 
                                             Pe = 0.2, 
                                             delta = d_seq[del])
                   
                   # Sample data until all variables have nonzero variance
                   okay <- 0
                   counter <- 0
                   while(okay==0) {
                     
                     # Generate Data
                     data1 <- IsingSampler(n = n, graph = G_obj$Graph1, thresholds = G_obj$thresh1, nIter = 200)
                     data2 <- IsingSampler(n = n, graph = G_obj$Graph2, thresholds = G_obj$thresh2, nIter = 200)
                     
                     # ind <- any(apply(data1, 2, sd) < 0.01) | any(apply(data2, 2, sd) < 0.01)
                     
                     # Check variances: at least two instances of each category, in each variable, in both groups
                     ind <- all(apply(data1, 2, function(x) min(c(sum(x==0), sum(x==1)))) > 1) & all(apply(data2, 2, function(x) min(c(sum(x==0), sum(x==1)))) > 1)
                     if(!ind) counter <- counter + 1 else okay <- 1
                     
                     if(counter > 5000) stop("Ridiculous long resampling.")
                     
                   } # end while
                   
                   # Combine data sets
                   data <- cbind(rbind(data1, data2), c(rep(0, n), rep(1, n)))
                   
                   v_counter[del] <- counter
                   
                   
                   # --------------------------------------------------------------
                   # ---------- Estimation ----------------------------------------
                   # --------------------------------------------------------------
                   
                   
                   # ------ 1.1) NCT: alpha=.05 -----------------------------------
                   
                   timer <- proc.time()[3]
                   
                   NCT_obj <- NCT(data1 = data1,
                                  data2 = data2,
                                  gamma = 0.25, # Default value of NCT
                                  it = 250, # to save some cost
                                  test.edges = TRUE,
                                  binary.data = TRUE,
                                  edges = "all",
                                  progressbar = FALSE)
                   
                   m_time[1, del] <-  proc.time()[3] - timer
                   m_time[2, del] <-  proc.time()[3] - timer # copy for consistency in data structure
                   
                   # Compute signficiant parameter differences
                   alpha <- .05
                   m_pvals <- NCT_obj$einv.pval
                   m_pvals[, 1] <- as.numeric(gsub('var','', m_pvals[, 1]))
                   m_pvals[, 2] <- as.numeric(gsub('var','', m_pvals[, 2]))
                   
                   m_NCT_sgnf <- matrix(0, p, p)
                   for(i in 1:n_edges) if(m_pvals[i, 3] < alpha) m_NCT_sgnf[m_pvals[i, 1], m_pvals[i, 2]] <- m_NCT_sgnf[m_pvals[i, 2], m_pvals[i, 1]] <- 1
                   
                   # Save in output object
                   a_out[, , 1, 1, del] <- m_NCT_sgnf
                   a_out[, , 2, 1, del] <- (NCT_obj$nw1-NCT_obj$nw2)* m_NCT_sgnf
                   
                   
                   # ------ 1.2) NCT: alpha=.01 -----------------------------------
                   
                   # Compute signficiant parameter differences
                   alpha <- .01
                   m_pvals <- NCT_obj$einv.pval
                   m_pvals[, 1] <- as.numeric(gsub('var','',m_pvals[, 1]))
                   m_pvals[, 2] <- as.numeric(gsub('var','',m_pvals[, 2]))
                   
                   m_NCT_sgnf <- matrix(0, p, p)
                   for(i in 1:n_edges) if(m_pvals[i, 3] < alpha) m_NCT_sgnf[m_pvals[i, 1], m_pvals[i, 2]] <- m_NCT_sgnf[m_pvals[i, 2], m_pvals[i, 1]] <- 1
                   
                   # Save in output object
                   a_out[, , 1, 2, del] <- m_NCT_sgnf
                   a_out[, , 2, 2, del] <- (NCT_obj$nw1-NCT_obj$nw2)* m_NCT_sgnf
                   
                   
                   # ------ 2.1) mgm: BIC + AND -----------------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("c", p), "c"),
                                  level = c(rep(2, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = 0,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "AND", 
                                  threshold = "none")
                   
                   m_time[3, del] <-  proc.time()[3] - timer
                   
                   # Get differences out of mgm model object
                   out_MNM <- f_Diffs_from_mgm(MNM_obj)
                   
                   # Save in output object
                   a_out[, , 1, 3, del] <- out_MNM$diff_logic
                   a_out[, , 2, 3, del] <- out_MNM$diff_value
                   
                   
                   # ------ 2.2) mgm: EBIC, gamma=0.25 + AND ----------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("c", p), "c"),
                                  level = c(rep(2, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = .25,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "AND", 
                                  threshold = "none")
                   
                   m_time[4, del] <-  proc.time()[3] - timer
                   
                   # Get differences out of mgm model object
                   out_MNM <- f_Diffs_from_mgm(MNM_obj)
                   
                   # Save in output object
                   a_out[, , 1, 4, del] <- out_MNM$diff_logic
                   a_out[, , 2, 4, del] <- out_MNM$diff_value
                   
                   
                   # ------ 2.3) mgm: CV, 10fold + AND ---------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("c", p), "c"),
                                  level = c(rep(2, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "CV",
                                  pbar = TRUE,
                                  # pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "AND", 
                                  threshold = "none")
                   
                   m_time[5, del] <-  proc.time()[3] - timer
                   
                   # Get differences out of mgm model object
                   out_MNM <- f_Diffs_from_mgm(MNM_obj)
                   
                   # Save in output object
                   a_out[, , 1, 5, del] <- out_MNM$diff_logic
                   a_out[, , 2, 5, del] <- out_MNM$diff_value
                   
                   mean(out_MNM$diff_logic[G_obj$G_diff!=0])
                   
                   
                   # ------ 2.4) mgm: BIC + OR -----------------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("c", p), "c"),
                                  level = c(rep(2, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = 0,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "OR", 
                                  threshold = "none")
                   
                   m_time[6, del] <-  proc.time()[3] - timer
                   
                   # Get differences out of mgm model object
                   out_MNM <- f_Diffs_from_mgm(MNM_obj)
                   
                   # Save in output object
                   a_out[, , 1, 6, del] <- out_MNM$diff_logic
                   a_out[, , 2, 6, del] <- out_MNM$diff_value
                   
                   
                   # ------ 2.5) mgm: EBIC, gamma=0.25 + OR ----------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("c", p), "c"),
                                  level = c(rep(2, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = .25,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "OR", 
                                  threshold = "none")
                   
                   m_time[7, del] <-  proc.time()[3] - timer
                   
                   # Get differences out of mgm model object
                   out_MNM <- f_Diffs_from_mgm(MNM_obj)
                   
                   # Save in output object
                   a_out[, , 1, 7, del] <- out_MNM$diff_logic
                   a_out[, , 2, 7, del] <- out_MNM$diff_value
                   
                   
                   # ------ 2.6) mgm: CV, 10fold + OR ---------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("c", p), "c"),
                                  level = c(rep(2, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "CV",
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "OR", 
                                  threshold = "none")
                   
                   m_time[8, del] <-  proc.time()[3] - timer
                   
                   # Get differences out of mgm model object
                   out_MNM <- f_Diffs_from_mgm(MNM_obj)
                   
                   # Save in output object
                   a_out[, , 1, 8, del] <- out_MNM$diff_logic
                   a_out[, , 2, 8, del] <- out_MNM$diff_value
                   
                   
                   # ------ 3.1) BGGM: Bayes Factor ---------------------------------
                   
                   timer <- proc.time()[3]
                   
                   if(ni>2) {
                     fit <- BGGM::ggm_compare_explore(data1, data2, type="binary")
                     select_graph <- BGGM::select(fit, post_prob = 0.50)
                     
                     # Save in output object
                     a_out[, , 1, 9, del] <- select_graph$adj_10
                     a_out[, , 2, 9, del] <- select_graph$pcor_mat_10  
                   } else {
                     
                     # No estimation for ni=1 (i.e., n=20)
                     a_out[, , 1, 9, del] <- matrix(NA, p, p)
                     a_out[, , 2, 9, del] <- matrix(NA, p, p)
                   }
                   
                   m_time[9, del] <-  proc.time()[3] - timer
                   
                   
                   # ------ 3.2) BGGM: Posterior Thresholding ---------------------------------
                   
                   timer <- proc.time()[3]
                   
                   if(ni>2) {
                     fit <- BGGM::ggm_compare_estimate(data1, data2, analytic = FALSE, type = "binary")
                     select_graph <- BGGM::select(fit, cred = 0.95)
                     
                     # Save in output object
                     a_out[, , 1, 10, del] <- select_graph$adj[[1]]
                     a_out[, , 2, 10, del] <- select_graph$pcor_adj[[1]]
                   } else {
                     # No estimation for ni=1 (i.e., n=20)
                     a_out[, , 1, 10, del] <- matrix(NA, p, p)
                     a_out[, , 2, 10, del] <- matrix(NA, p, p)
                   }
                   
                   m_time[10, del] <-  proc.time()[3] - timer
                   
                   
                   
                   # Print Progress
                   print(paste("ni = ", ni, " del = ", del))
                   
                 } # end for: d
                 
                 
                 outlist <- list("array" = a_out, 
                                 "m_time" = m_time, 
                                 "counter" = v_counter)
                 
                 return(outlist)
                 
               } # end for: n


stopCluster(cl)


# --------------------------------------------------------------
# ---------- Postprocess & Save ---------------------------------
# --------------------------------------------------------------

# Turn list into another array dimensions
n_methods <- 10
out_final <- array(NA, dim = c(p, p, 2, n_methods, 3, 10))
a_time <-  array(NA, dim = c(10, n_methods, 3))
m_counter <- matrix(NA, nrow=10, ncol=3)

for(ni in 1:10) out_final[, , , , , ni] <- out[[ni]]$array
for(ni in 1:10) a_time[ni, , ] <- out[[ni]]$m_time
for(ni in 1:10) m_counter[ni, ] <- out[[ni]]$counter


# Save true model
d_seq <- c(.15, .3, 0.6) 

G_obj_d1 <- f_GenModel_Ising(iter = iter_set, 
                             p = 17, 
                             Pe = 0.2, 
                             delta = d_seq[1])
G_obj_d2 <- f_GenModel_Ising(iter = iter_set, 
                             p = 17, 
                             Pe = 0.2, 
                             delta = d_seq[2])
G_obj_d3 <- f_GenModel_Ising(iter = iter_set, 
                             p = 17, 
                             Pe = 0.2, 
                             delta = d_seq[3])

# Everything into outlist to save
outlist <- list("array" = out_final, 
                "timing" = a_time, 
                "counter" = m_counter,
                "truemod_d1"=G_obj_d1,
                "truemod_d2"=G_obj_d2,
                "truemod_d3"=G_obj_d3)


saveRDS(outlist, file=paste0("Simresults_Ising_Iter", iter, ".RDS"))

print("Timing: ")
proc.time()[3] - timer_total

