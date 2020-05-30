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
library(MASS)

# Identify group differences
library(NetworkComparisonTest)
library(EstimateGroupNetwork)
library(mgm)
library(BGGM)
library(psychonetrics)
library(dplyr)

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
# ---------- Aux functions -------------------------------------
# --------------------------------------------------------------

source("aux_functions.R") # source helper functions

# --------------------------------------------------------------
# ---------- Generate Data -------------------------------------
# --------------------------------------------------------------

cluster <- 10
cl <- makeCluster(cluster, outfile="")
registerDoParallel(cl)

out <- foreach(ni = 1:10,
               .packages=c("NetworkComparisonTest", "EstimateGroupNetwork", 
                           "mgm", "BGGM", "dplyr", "psychonetrics"),
               .export=c("n_seq", "iter", "f_GenModel", "n_edges", "mle_fisher", "f_Diffs_from_mgm", 
                         "psychonetricsGGM"),
               .verbose=TRUE) %dopar% {
                 
                 
                 # Storage
                 a_out <- array(NA, dim = c(p, p, 2, 18, 3)) 
                 # p, p, 2 results (real: difference, logical: difference), 
                 # 18 estimation methods; 10 n-variations; 3 for delta theta
                 # need 2 results, because BGGM only gives logical
                 
                 m_time <- matrix(NA, 18, 3) # computation time
                 
                 d_seq <- c(0.05, .1, .2)
                 
                 for(del in 1:3) {
                   
                   # Generate Graph
                   G_obj <- f_GenModel(iter = iter, 
                                       p = 17, 
                                       Pe = 0.2, 
                                       delta = d_seq[del])
                   
                   n <- n_seq[ni]
                   
                   data1 <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma = G_obj$Sigma1)
                   data2 <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma = G_obj$Sigma2)
                   data <- cbind(rbind(data1, data2), c(rep(0, n), rep(1, n)))
                   l_data <- list(data1, data2)
                   
                   
                   # --------------------------------------------------------------
                   # ---------- Estimation ----------------------------------------
                   # --------------------------------------------------------------
                   
                   
                   # ------ 1.1) NCT: alpha=.05 -----------------------------------
                   
                   timer <- proc.time()[3]
                   
                   NCT_obj <- NCT(data1 = data1,
                                  data2 = data2,
                                  gamma = .25, # Default value of NCT
                                  it = 1000, # run 1000 iterations for actual sim
                                  test.edges = TRUE,
                                  binary.data = FALSE,
                                  edges = "all",
                                  progressbar = FALSE, 
                                  make.positive.definite = TRUE)
                   
                   m_time[1, del] <-  proc.time()[3] - timer
                   m_time[2, del] <-  proc.time()[3] - timer
                   
                   # Compute signficiant parameter differences
                   alpha <- .05
                   m_pvals <- NCT_obj$einv.pval
                   m_pvals[, 1] <- as.numeric(gsub('var','',m_pvals[, 1]))
                   m_pvals[, 2] <- as.numeric(gsub('var','',m_pvals[, 2]))
                   
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
                   
                   
                   # ------ 2.1) FGL: BIC -----------------------------------------
                   
                   timer <- proc.time()[3]
                   
                   FGL_obj <- EstimateGroupNetwork(X = l_data,
                                                   inputType = "list.of.dataframes",
                                                   method = "InformationCriterion",
                                                   strategy = "sequential",
                                                   gamma = 0)
                   
                   m_time[3, del] <-  proc.time()[3] - timer
                   
                   # Save in output object
                   a_out[, , 1, 3, del] <- FGL_obj[[1]] != FGL_obj[[2]]
                   a_out[, , 2, 3, del] <- FGL_obj[[1]] - FGL_obj[[2]]
                   
                   
                   # ------ 2.2) FGL: EBIC; gamma=.25 -----------------------------
                   
                   timer <- proc.time()[3]
                   
                   FGL_obj <- EstimateGroupNetwork(X = l_data,
                                                   inputType = "list.of.dataframes",
                                                   method = "InformationCriterion",
                                                   strategy = "sequential",
                                                   gamma = .25)
                   
                   m_time[4, del] <-  proc.time()[3] - timer
                   
                   # Save in output object
                   a_out[, , 1, 4, del] <- FGL_obj[[1]] != FGL_obj[[2]]
                   a_out[, , 2, 4, del] <- FGL_obj[[1]] - FGL_obj[[2]]
                   
                   
                   # ------ 2.3) FGL: CV 10fold -----------------------------------------
                   
                   timer <- proc.time()[3]
                   
                   FGL_obj <- EstimateGroupNetwork(X = l_data,
                                                   inputType = "list.of.dataframes",
                                                   method = "crossvalidation",
                                                   strategy = "sequential",
                                                   gamma = 0)
                   
                   m_time[5, del] <-  proc.time()[3] - timer
                   
                   
                   # Save in output object
                   a_out[, , 1, 5, del] <- FGL_obj[[1]] != FGL_obj[[2]]
                   a_out[, , 2, 5, del] <- FGL_obj[[1]] - FGL_obj[[2]]
                   
                   
                   # ------ 3.1) BGGM with Bayes Factor; BF_cut=1 ------------------
                   
                   timer <- proc.time()[3]
                   
                   fit <- BGGM::ggm_compare_explore(data1, data2)
                   select_graph <- BGGM::select(fit, post_prob = 0.50)
                   
                   m_time[6, del] <-  proc.time()[3] - timer
                   
                   # Save in output object
                   a_out[, , 1, 6, del] <- select_graph$adj_10
                   a_out[, , 2, 6, del] <- select_graph$pcor_mat_10
                   
                   
                   # ------ 3.2) BGGM difference from posterior (analytical) ----------
                   
                   timer <- proc.time()[3]
                   
                   fit <- BGGM::ggm_compare_estimate(data1, data2, analytic = TRUE)
                   select_graph <- BGGM::select(fit, cred = 0.95)
                   
                   m_time[7, del] <-  proc.time()[3] - timer
                   
                   # Save in output object
                   a_out[, , 1, 7, del] <- select_graph$adj[[1]]
                   a_out[, , 2, 7, del] <- select_graph$pcor_adj[[1]]
                   
                   
                   # ------ 3.3) BGGM difference from posterior (sampling) ----------
                   
                   timer <- proc.time()[3]
                   
                   fit <- BGGM::ggm_compare_estimate(data1, data2, analytic = FALSE)
                   select_graph <- BGGM::select(fit, cred = 0.95)
                   
                   m_time[7, del] <-  proc.time()[3] - timer
                   
                   # Save in output object
                   a_out[, , 1, 8, del] <- select_graph$adj[[1]]
                   a_out[, , 2, 8, del] <- select_graph$pcor_adj[[1]]
                   
                   
                   # ------ 4.1) mgm: BIC + AND -----------------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("g", p), "c"),
                                  level = c(rep(1, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = 0,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "AND", 
                                  threshold = "none")
                   
                   m_time[9, del] <-  proc.time()[3] - timer
                   
                   # Get matrix of differences
                   m_3way <- MNM_obj$interactions$indicator[[2]]
                   m_3way <- matrix(m_3way, ncol=3)
                   
                   # Get differences out of mgm model object 
                   out_MNM <- f_Diffs_from_mgm(m_3way)
                   
                   # Save in output object
                   a_out[, , 1, 9, del] <- out_MNM$diff_logic
                   a_out[, , 2, 9, del] <- out_MNM$diff_value
                   
                   
                   # ------ 4.2) mgm: EBIC, gamma=0.25 + AND ----------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("g", p), "c"),
                                  level = c(rep(1, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = .25,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "AND", 
                                  threshold = "none")
                   
                   m_time[10, del] <-  proc.time()[3] - timer
                   
                   # Get matrix of differences
                   m_3way <- MNM_obj$interactions$indicator[[2]]
                   m_3way <- matrix(m_3way, ncol=3)
                   
                   # Get differences out of mgm model object 
                   out_MNM <- f_Diffs_from_mgm(m_3way)
                   
                   # Save in output object
                   a_out[, , 1, 10, del] <- out_MNM$diff_logic
                   a_out[, , 2, 10, del] <- out_MNM$diff_value
                   
                   
                   # ------ 4.3) mgm: CV, 10fold + AND ---------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("g", p), "c"),
                                  level = c(rep(1, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "CV",
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "AND", 
                                  threshold = "none")
                   
                   m_time[11, del] <-  proc.time()[3] - timer
                   
                   # Get matrix of differences
                   m_3way <- MNM_obj$interactions$indicator[[2]]
                   m_3way <- matrix(m_3way, ncol=3)
                   
                   # Get differences out of mgm model object 
                   out_MNM <- f_Diffs_from_mgm(m_3way)
                   
                   # Save in output object
                   a_out[, , 1, 11, del] <- out_MNM$diff_logic
                   a_out[, , 2, 11, del] <- out_MNM$diff_value
                   
                   
                   # ------ 4.4) mgm: BIC + OR -----------------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("g", p), "c"),
                                  level = c(rep(1, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = 0,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "OR", 
                                  threshold = "none")
                   
                   m_time[12, del] <-  proc.time()[3] - timer
                   
                   # Get matrix of differences
                   m_3way <- MNM_obj$interactions$indicator[[2]]
                   m_3way <- matrix(m_3way, ncol=3)
                   
                   # Get differences out of mgm model object 
                   out_MNM <- f_Diffs_from_mgm(m_3way)
                   
                   # Save in output object
                   a_out[, , 1, 12, del] <- out_MNM$diff_logic
                   a_out[, , 2, 12, del] <- out_MNM$diff_value
                   
                   
                   # ------ 4.5) mgm: EBIC, gamma=0.25 + OR ----------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("g", p), "c"),
                                  level = c(rep(1, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "EBIC",
                                  lambdaGam = .25,
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "OR", 
                                  threshold = "none")
                   
                   m_time[13, del] <-  proc.time()[3] - timer
                   
                   # Get matrix of differences
                   m_3way <- MNM_obj$interactions$indicator[[2]]
                   m_3way <- matrix(m_3way, ncol=3)
                   
                   # Get differences out of mgm model object 
                   out_MNM <- f_Diffs_from_mgm(m_3way)
                   
                   # Save in output object
                   a_out[, , 1, 13, del] <- out_MNM$diff_logic
                   a_out[, , 2, 13, del] <- out_MNM$diff_value
                   
                   
                   # ------ 4.6) mgm: CV, 10fold + OR ---------------------------------
                   
                   timer <- proc.time()[3]
                   
                   MNM_obj <- mgm(data = data,
                                  type = c(rep("g", p), "c"),
                                  level = c(rep(1, p), 2),
                                  moderators = p+1,
                                  lambdaSel = "CV",
                                  pbar = FALSE,
                                  signInfo = FALSE,
                                  ruleReg = "OR", 
                                  threshold = "none")
                   
                   m_time[14, del] <-  proc.time()[3] - timer
                   
                   # Get matrix of differences
                   m_3way <- MNM_obj$interactions$indicator[[2]]
                   m_3way <- matrix(m_3way, ncol=3)
                   
                   # Get differences out of mgm model object 
                   out_MNM <- f_Diffs_from_mgm(m_3way)
                   
                   # Save in output object
                   a_out[, , 1, 14, del] <- out_MNM$diff_logic
                   a_out[, , 2, 14, del] <- out_MNM$diff_value
                   
                   
                  
                   # ------ 5.1) Fisher Z-transform; alpha = 0.05 -------------------------
                   
                   timer <- proc.time()[3]
                   out_fisher <- mle_fisher(data1, data2, alpha = 0.05)
                   m_time[15, del] <-  proc.time()[3] - timer
                   
                   a_out[, , 1, 15, del] <- out_fisher$diff_logic
                   a_out[, , 2, 15, del] <- out_fisher$diff_value
                   
                   
                   # ------ 5.2) Fisher Z-transform; alpha = 0.01 -------------------------                 
                   
                   timer <- proc.time()[3]
                   out_fisher <- mle_fisher(data1, data2, alpha = 0.01)
                   m_time[16, del] <-  proc.time()[3] - timer
                   
                   a_out[, , 1, 16, del] <- out_fisher$diff_logic
                   a_out[, , 2, 16, del] <- out_fisher$diff_value
                   
                   
                   # ------ 6.1) Psychonetrics pruning approach, alpha = 0.05 -------------------------                 
                   
                   timer <- proc.time()[3]
                   
                   if(ni == 1) {
                     on_out <- matrix(NA, p, p) 
                   } else {
                     on_out <- psychonetricsGGM(data_0 = data1, data_1 = data2, alpha = 0.05)
                   }
                     
                   m_time[17, del] <-  proc.time()[3] - timer
                   
                   a_out[, , 1, 17, del] <- on_out != 0
                   a_out[, , 2, 17, del] <- on_out
                   
                   
                   # ------ 6.2) Psychonetrics pruning approach, alpha = 0.01 -------------------------                 
                   
                   timer <- proc.time()[3]
                   
                   if(ni == 1) {
                     on_out <- matrix(NA, p, p) 
                   } else {
                     on_out <- psychonetricsGGM(data_0 = data1, data_1 = data2, alpha = 0.01)
                   }
                   
                   m_time[17, del] <-  proc.time()[3] - timer
                   
                   a_out[, , 1, 18, del] <- on_out != 0
                   a_out[, , 2, 18, del] <- on_out
                   
                   
                   print(paste("ni = ", ni, " del = ", del))
                   
                 } # end for: d
                 
                 # Save all
                 outlist <- list("array" = a_out, 
                                 "time" = m_time)
                 
                 return(outlist)
                 
               } # end for: n


stopCluster(cl)


# --------------------------------------------------------------
# ---------- Postprocess & Save ---------------------------------
# --------------------------------------------------------------

# Turn list into another array dimensions
out_final <- array(NA, dim = c(p, p, 2, 18, 3, 10))
a_time <-  array(NA, dim = c(10, 18, 3))
for(ni in 1:10) out_final[, , , , , ni] <- out[[ni]]$array
for(ni in 1:10) a_time[ni, , ] <- out[[ni]]$time


# Save true model
d_seq <- c(0.05, .1, .2)

G_obj_d1 <- f_GenModel(iter = iter, 
                       p = 17, 
                       Pe = 0.2, 
                       delta = d_seq[1])
G_obj_d2 <- f_GenModel(iter = iter, 
                       p = 17, 
                       Pe = 0.2, 
                       delta = d_seq[2])
G_obj_d3 <- f_GenModel(iter = iter, 
                       p = 17, 
                       Pe = 0.2, 
                       delta = d_seq[3])

# Everything into outlist to save
outlist <- list("array" = out_final, 
                "timing" = a_time, 
                "truemod_d1"=G_obj_d1,
                "truemod_d2"=G_obj_d2,
                "truemod_d3"=G_obj_d3)


saveRDS(outlist, file=paste0("Simresults_Iter", iter, ".RDS"))

print("Timing: ")
proc.time()[3] - timer_total


