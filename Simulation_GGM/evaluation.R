# jonashaslbeck@gmail.com; May 20, 2021

# --------------------------------------------------------------
# ---------- Load & Prepare Data -------------------------------
# --------------------------------------------------------------

# Get Data

simResDir <- "Simulation_GGM/output/"

v_files <- list.files(simResDir)
n_files <- length(v_files)
l_files <- list()
i_files <- 1:n_files
nIter <- length(i_files)

for(i in i_files) l_files[[i]] <- readRDS(paste0(simResDir, v_files[i]))

# Some simulation settings
p <- 17
n_seq <- round(exp(seq(3, 8.5171, length=10)))
n_methods <- 18

# Create overall output arrays:
a_out <- array(NA, dim = c(p, p, 2, n_methods+1, 3, 10, nIter))
a_time <-  array(NA, dim = c(10, n_methods, 3, nIter))

# Loop in simulation results
for(i in 1:nIter) {
  
  # Estimates
  a_out[, , , 1:n_methods, , , i] <- l_files[[i]]$array
  a_time[, , , i] <- l_files[[i]]$timing
  
  # True 
  a_out[, , 1, n_methods+1, 1, , i] <- (l_files[[i]]$truemod_d1$G_diff != 0) * 1
  a_out[, , 1, n_methods+1, 2, , i] <- (l_files[[i]]$truemod_d2$G_diff != 0) * 1
  a_out[, , 1, n_methods+1, 3, , i] <- (l_files[[i]]$truemod_d3$G_diff != 0) * 1
  
  a_out[, , 2, n_methods+1, 1, , i] <- l_files[[i]]$truemod_d1$G_diff
  a_out[, , 2, n_methods+1, 2, , i] <- l_files[[i]]$truemod_d2$G_diff
  a_out[, , 2, n_methods+1, 3, , i] <- l_files[[i]]$truemod_d3$G_diff
  
} # end for




# --------------------------------------------------------------
# ---------- Count errors in psychonetrics ---------------------
# --------------------------------------------------------------

# example of NA in psychonetrics
a_out[, , 2, 17, 1, 1, 3] 
a_out[, , 2, 18, 1, 1, 3] 

# Get all NAs
tb <- table(is.na(a_out[, , 2, 17, , , ]))
tb / sum(tb)

# Get all NAs conditional on sample size
tb <- table(is.na(a_out[, , 2, 17, , 8, ]))
tb / sum(tb)
# OK: around 1% failure rate for ni > 1


# --------------------------------------------------------------
# ---------- Preprocess ----------------------------------------
# --------------------------------------------------------------

# Sensitivity, Precision, Estimation Errors, Specificity; probably: average estimation error (average across detected differences)

a_P <- array(NA, dim = c(6, n_methods, 3, 10, nIter))

for(m in 1:n_methods) { # method
  for(n in 1:10) { # sample size
    for(d in 1:3) { # delta theta
      for(i in 1:nIter) { # iterations
        
        # Sensitivity, Precision, Specificity
        true <- a_out[, , 1, n_methods+1, d, n, i]
        absent <- true != 1
        a_P[1, m, d, n, i] <- mean(a_out[, , 1, m, d, n, i][true==1]) # Sensitivity
        a_P[2, m, d, n, i] <- mean(true[a_out[, , 1, m, d, n, i]!=0]) # Precision
        a_P[6, m, d, n, i] <- mean((a_out[, , 1, m, d, n, i]==0)[absent]) # Specificity
        
        # --- Estimation Errors ---
        # Among present Edges
        true_est <- a_out[, , 2, n_methods+1, d, n, i][true==1]
        est_est <- - a_out[, , 2, m, d, n, i][true==1]
        a_P[3, m, d, n, i] <- mean(abs(true_est-est_est))
        
        # Among absent Edges
        true_est <- a_out[, , 2, n_methods+1, d, n, i][true==0]
        est_est <- - a_out[, , 2, m, d, n, i][true==0]
        a_P[4, m, d, n, i] <- mean(abs(true_est-est_est))
        
        # Combined
        true_est <- a_out[, , 2, n_methods+1, d, n, i]
        est_est <- - a_out[, , 2, m, d, n, i]
        a_P[5, m, d, n, i] <- mean(abs(true_est-est_est))
        
      }
    }
  }
}

# Mean aggregate across iterations
a_P_means <- apply(a_P, 1:4, function(x) mean(x, na.rm = TRUE))
a_P_freqNA <- apply(a_P, 1:4, function(x) sum(is.na(x)))
a_P_means[a_P_freqNA > 90] <- NA

# Save aggregate result; we need it for the Sparsity comparison in the appendix
saveRDS(a_P_means, "Simulation_GGM/GGM_results.RDS")


# --------------------------------------------------------------
# ---------- Plotting Figure 1 ---------------------------------
# --------------------------------------------------------------

library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
cols <- cols[-c(6:7)]
lwd=1.5
title.cex <- 1

sc <- 2
pdf("Simulation_GGM/figures/Fig_Sim_GGM.pdf", width = 4.3*sc, height = 3.5*sc)

# ----- Setup Layout -----

lmat <- rbind(
  rep(26,5),
  c(1:5), 
  c(6, 14:17),
  c(7, 18:21), 
  c(8, 22:25), 
  c(9:13))

lo <- layout(mat = lmat, 
             widths = c(.1, 1, 1, 1, 1), 
             heights = c(.5, .1, 1, 1, 1, .1))

# layout.show(lo)

# ------ Plot all the labels into the layout ------

par(mar=rep(0, 4))

plot.new() # empty

# titles horizontal (delta)
v_h_titles <- c("Sensitivity", "Precision", "Estimation Error (present)", "Estimation Error (absent)")
for(s in 1:4) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, v_h_titles[s], cex=1.4, srt = 0)  
}

# titles vertical (errors)
v_v_titles <- list()
v_v_titles[[1]] <- bquote(paste(Delta, theta, " = 0.05"))
v_v_titles[[2]] <- bquote(paste(Delta, theta, " = 0.10"))
v_v_titles[[3]] <- bquote(paste(Delta, theta, " = 0.20"))
for(s in 1:3) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, v_v_titles[[s]], srt = 90, cex=1.5)  
}


plot.new() # empty

# titles bottom
for(s in 1:4) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, "  Sample Size per Group", cex=1.1)  
}

# ------ Plot the data ------

par(mar=c(3,3,1,2))

true_par <- c(.05, .1, .2)

# ----- Plotting Sensitivity and Precision ------

for(d in 1:3) {
  
  for(s in c(1,2,3,4)) {
    
    if(s < 3) {
      # Set up plotting area: precision /sensitivity
      plot.new()
      plot.window(ylim=c(0,1), xlim=c(1,10))
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, las=2)
    } else if(s==3) {
      # Set up plotting area: estimation error (present)
      plot.new()
      plot.window(ylim=c(0,.25), xlim=c(1,10))
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, las=2)
      # abline(h=true_par[d], lty=1, color="lightgrey")
    } else {
      # Set up plotting area: estimation error (absent)
      plot.new()
      plot.window(ylim=c(0,.25), xlim=c(1,10)) 
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, las=2)
    }
    
    
    # Plot Data:
    # Plot Data:
    # NCT
    lines(a_P_means[s, 1, d, ], col=cols[1], lty=1, lwd=lwd) # 0.05
    # lines(a_P_means[s, 2, d, ], col=cols[1], lty=2, lwd=lwd) # 0.01
    
    # Fisher / NHST
    lines(a_P_means[s, 15, d, ], col=cols[5], lty=1, lwd=lwd) # 0.05
    # lines(a_P_means[s, 16, d, ], col=cols[5], lty=2, lwd=lwd) # 0.01
    
    # FGL
    # lines(a_P_means[s, 3, d, ], col=cols[2], lty=1, lwd=lwd) # BIC
    lines(a_P_means[s, 4, d, ], col=cols[2], lty=1, lwd=lwd) # EBIC
    lines(a_P_means[s, 5, d, ], col=cols[2], lty=2, lwd=lwd) # CV
    
    # BGGM
    lines(a_P_means[s, 6, d, ], col=cols[3], lty=1, lwd=lwd) # Bayes-Factor
    lines(a_P_means[s, 8, d, ], col=cols[3], lty=2, lwd=lwd) # posterior difference (sampling)
    # lines(a_P_means[s, 7, d, ], col=cols[3], lty=3, lwd=lwd) # posterior difference (analytical)
    
    # MGM
    # lines(a_P_means[s, 9, d, ], col=cols[4], lty=1, lwd=lwd) # BIC + AND
    lines(a_P_means[s, 10, d, ], col=cols[4], lty=1, lwd=lwd) # EBIC + AND
    lines(a_P_means[s, 11, d, ], col=cols[4], lty=2, lwd=lwd) # CV + AND
    
    # lines(a_P_means[s, 12, d, ], col=cols[4], lty=4, lwd=lwd) # BIC + OR
    lines(a_P_means[s, 13, d, ], col=cols[4], lty=3, lwd=lwd) # EBIC + OR
    lines(a_P_means[s, 14, d, ], col=cols[4], lty=4, lwd=lwd) # CV + OR
    
    # Psychonetrics
    psychonetrics_results <- a_P_means[s, 17, d, ]
    psychonetrics_results[1] <- NA
    lines(psychonetrics_results, col=cols[6], lty=1, lwd=lwd) # alpha = 0.05
    # lines(a_P_means[s, 18, d, ], col=cols[6], lty=4, lwd=lwd) # alpha = 0.01
    
    
  } # for: s
  
  
} # for: d


# Plot legend in lower panel (no 26)
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

legend_cex <- 1.15
legend_y <- 1

# Batch 1: NCT + Fisher
legend(0.1, legend_y, legend = c(as.expression(bquote(paste("NCT, ", alpha, " = 0.05"))),
                           # as.expression(bquote(paste("NCT, ", alpha, " = 0.01"))),
                           as.expression(bquote(paste("Fisher, ", alpha, " = 0.05"))),
                           # as.expression(bquote(paste("Fisher, ", alpha, " = 0.01")))
                           as.expression(bquote(paste("Partial Pruning, ", alpha, " = 0.05")))
                           ),
       col = c(cols[1], cols[5],
               cols[6]),
       lwd = lwd, lty=c(1,1,
                        1),
       bty = "n", cex = legend_cex)

# Batch 2: FGL + BGGM
legend(0.4, legend_y, legend = c(as.expression(bquote(paste("FGL EBIC, ", gamma, " = 0.25"))),
                           "FGL CV",
                           "BGGM, Bayes Factor",
                           "BGGM, Post Diff"),
       col = c(cols[2], cols[2],
               cols[3], cols[3]),
       lwd = lwd, lty=c(1,2,
                        1,2),
       bty = "n", cex = legend_cex)


# Batch 3: Moderation
legend(.7, legend_y, legend = c(#"MNM BIC",
  as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + AND"))),
  "MNM CV + AND",
  #"MNM BIC + OR",
  as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + OR"))),
  "MNM CV + OR"),
  col = c(cols[4], cols[4], cols[4], cols[4]),
  lwd = lwd, lty=c(1,2,3,4),
  bty = "n", cex = legend_cex)


dev.off()


# --------------------------------------------------------------
# ---------- Marginal Performance ------------------------------
# --------------------------------------------------------------

# Marginalize over n and delta theta

method_labels <- c("NCT1", "NCT2", 
                   "FGL_bic", "FGL_ebic", "FGL_cv",
                   "BGGM_BF", "BGGM_P_an", "BGGM_P_sa",
                   "mgm_bic", "mgm_ebic", "mgm_cv", "mgm_bic_or", "mgm_ebic_or", "mgm_cv_or", 
                   "Fisher1", "Fisher2",
                   "psychonetrics_al1", "psychonetrics_al2") # this is the order in the vector

methods_labels_proper <- c(c(as.expression(bquote(paste("NCT, ", alpha, " = 0.05"))),
                             as.expression(bquote(paste("Fisher, ", alpha, " = 0.05")))
), c(as.expression(bquote(paste("FGL EBIC, ", gamma, " = 0.25"))),
     "FGL CV",
     "BGGM, Bayes Factor",
     "BGGM, Post Diff"),
c(as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + AND"))),
  "MNM CV + AND",
  as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + OR"))),
  "MNM CV + OR"), 
as.expression(bquote(paste("Partial Pruning, ", alpha, " = 0.05")))) # this is the order in the above figure


methods_select <- c(1, 15, 4, 5, 6, 8, 10:11, 13:14, 17) # Select same as in graph above
method_labels[methods_select] # Sanity check

methods_labels_proper

# Aggregate over N
a_P_means_aggN <- apply(a_P_means, 1:2, function(x) mean(x, na.rm=TRUE))
a_P_qntls_aggN <- apply(a_P_means, 1:2, function(x) quantile(x, na.rm=TRUE, probs = c(.25, .75)))
dim(a_P_means_aggN)

# Overview matrix: est(present)
m_ov <- matrix(NA, 11, 2)
m_ov[, 2] <- round(a_P_means_aggN[3, methods_select], 4)
m_ov[, 1] <- method_labels[methods_select]
or_1 <- order(m_ov[, 2])
m_ov_ord <- m_ov[or_1, ]

# Overview matrix: est(absent)
m_ov2 <- matrix(NA, 11, 2)
m_ov2[, 2] <- round(a_P_means_aggN[4, methods_select], 4)
m_ov2[, 1] <- method_labels[methods_select]
m_ov2_ord <- m_ov2[order(as.numeric(m_ov2[, 2])), ]


# Proper figure

sc <- 1.2
pdf("Simulation_GGM/figures/Fig_Sim_marginal_Gauss.pdf", width=7*sc, height=4.5*sc)

par(mar=c(9.2,5,2,2))
plot.new()
plot.window(ylim=c(0, .15), xlim=c(1,11))
axis(2, las=2)
# axis(1, labels = methods_labels_proper[methods_select][or_1], at=1:11, las=2)
# cols_fix <- rep("black", 11)

axis(1, labels = rep("", 11), at=1:11, 
     las=2, cex.axis=.8)

cols_fix <- cols[c(1,5,2,2,3,3,4,4,4,4,6)[or_1]]
for(i in 1:11) axis(1, labels = methods_labels_proper[or_1][i], at=i, 
                    las=2, cex.axis=.79, col.axis=cols_fix[i])

title(ylab="Estimation error", line=3.5)

## error: present
# mean
points(m_ov[, 2][or_1], pch=20, cex=1.3)
# upper/lower quantile
segments(x0 = 1:11, y0 = a_P_qntls_aggN[1, 3, ][methods_select][or_1], 
         x1 = 1:11, y1 = a_P_qntls_aggN[2, 3, ][methods_select][or_1], 
         lwd=1.5)


## error: absent
# mean
points(m_ov2[, 2][or_1], pch=20, col="grey", cex=1.3)
# upper/lower quantile
segments(x0 = 1:11, y0 = a_P_qntls_aggN[1, 4, ][methods_select][or_1], 
         x1 = 1:11, y1 = a_P_qntls_aggN[2, 4, ][methods_select][or_1], 
         col="grey", 
         lwd=1.5)

legend("topright", legend=c("Mean & 25-75% quantiles (present)", "Mean & 25-75% quantiles (absent)"), 
       col=c("black", "grey"), bty="n", text.col=c("black", "grey"))


dev.off()

# --------------------------------------------------------------
# ---------- Plotting Figure 6 in Appendix (Extended Fig1) -----
# --------------------------------------------------------------

# This figure shows results of some additional (specifications of) methods
# that were not shown in the main text

library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
cols <- cols[-c(6:7)]
lwd=1.5
title.cex <- 1

sc <- 2
pdf("Simulation_GGM/figures/Fig_Sim_GGM_APP.pdf", width = 4.3*sc, height = 3.5*sc)

# ----- Setup Layout -----

lmat <- rbind(rep(26,5),
              c(1:5), 
              c(6, 14:17),
              c(7, 18:21), 
              c(8, 22:25), 
              c(9:13))

lo <- layout(mat = lmat, 
             widths = c(.1, 1, 1, 1, 1), 
             heights = c(.5, .1, 1, 1, 1, .1))

# layout.show(lo)

# ------ Plot all the labels into the layout ------

par(mar=rep(0, 4))

plot.new() # empty

# titles horizontal (delta)
v_h_titles <- c("Sensitivity", "Precision", "Estimation Error (present)", "Estimation Error (absent)")
for(s in 1:4) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, v_h_titles[s], cex=1.4, srt = 0)  
}

# titles vertical (errors)
v_v_titles <- list()
v_v_titles[[1]] <- bquote(paste(Delta, theta, " = 0.05"))
v_v_titles[[2]] <- bquote(paste(Delta, theta, " = 0.10"))
v_v_titles[[3]] <- bquote(paste(Delta, theta, " = 0.20"))
for(s in 1:3) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, v_v_titles[[s]], srt = 90, cex=1.5)  
}


plot.new() # empty

# titles bottom
for(s in 1:4) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5, .5, "  Sample Size per Group", cex=1.1)  
}

# ------ Plot the data ------

par(mar=c(3,3,1,2))

true_par <- c(.05, .1, .2)

# ----- Plotting Sensitivity and Precision ------

for(d in 1:3) {
  
  for(s in c(1,2,3,4)) {
    
    if(s < 3) {
      # Set up plotting area: precision /sensitivity
      plot.new()
      plot.window(ylim=c(0,1), xlim=c(1,10))
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, las=2)
    } else if(s==3) {
      # Set up plotting area: estimation error (present)
      plot.new()
      plot.window(ylim=c(0,.25), xlim=c(1,10))
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, las=2)
      # abline(h=true_par[d], lty=1, color="lightgrey")
    } else {
      # Set up plotting area: estimation error (absent)
      plot.new()
      plot.window(ylim=c(0,.25), xlim=c(1,10)) 
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, las=2)
    }
    
    
    # Plot Data:
    # Plot Data:
    # NCT
    lines(a_P_means[s, 1, d, ], col=cols[1], lty=1, lwd=lwd) # 0.05
    lines(a_P_means[s, 2, d, ], col=cols[1], lty=2, lwd=lwd) # 0.01
    
    # Fisher / NHST
    lines(a_P_means[s, 15, d, ], col=cols[5], lty=1, lwd=lwd) # 0.05
    lines(a_P_means[s, 16, d, ], col=cols[5], lty=2, lwd=lwd) # 0.01
    
    # FGL
    # lines(a_P_means[s, 3, d, ], col=cols[2], lty=1, lwd=lwd) # BIC
    lines(a_P_means[s, 4, d, ], col=cols[2], lty=1, lwd=lwd) # EBIC
    lines(a_P_means[s, 5, d, ], col=cols[2], lty=2, lwd=lwd) # CV
    
    # BGGM
    lines(a_P_means[s, 6, d, ], col=cols[3], lty=1, lwd=lwd) # Bayes-Factor
    lines(a_P_means[s, 8, d, ], col=cols[3], lty=2, lwd=lwd) # posterior difference (sampling)
    lines(a_P_means[s, 7, d, ], col=cols[3], lty=3, lwd=lwd) # posterior difference (analytical)
    
    # MGM
    # lines(a_P_means[s, 9, d, ], col=cols[4], lty=1, lwd=lwd) # BIC + AND
    lines(a_P_means[s, 10, d, ], col=cols[4], lty=1, lwd=lwd) # EBIC + AND
    lines(a_P_means[s, 11, d, ], col=cols[4], lty=2, lwd=lwd) # CV + AND
    
    # lines(a_P_means[s, 12, d, ], col=cols[4], lty=4, lwd=lwd) # BIC + OR
    lines(a_P_means[s, 13, d, ], col=cols[4], lty=3, lwd=lwd) # EBIC + OR
    lines(a_P_means[s, 14, d, ], col=cols[4], lty=4, lwd=lwd) # CV + OR
    
    # Psychonetrics
    psychonetrics_results <- a_P_means[s, 17, d, ]
    psychonetrics_results[1] <- NA
    lines(psychonetrics_results, col=cols[6], lty=1, lwd=lwd) # alpha = 0.05
    psychonetrics_results2 <- a_P_means[s, 18, d, ]
    psychonetrics_results2[1] <- NA
    lines(psychonetrics_results2, col=cols[6], lty=2, lwd=lwd) # alpha = 0.01
    
  } # for: s
  
} # for: d


# Plot legend in lower panel (no 26)
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

legend_cex <- 1.15
legend_y <- 1

# Batch 1: NCT + Fisher 
legend(0.02, legend_y, legend = c(as.expression(bquote(paste("NCT, ", alpha, " = 0.05"))),
                                  as.expression(bquote(paste("NCT, ", alpha, " = 0.01"))),
                                  as.expression(bquote(paste("Fisher, ", alpha, " = 0.05"))),
                                  as.expression(bquote(paste("Fisher, ", alpha, " = 0.01")))
),
col = c(cols[1], cols[1],
        cols[5], cols[5]),
lwd = lwd, lty=c(1,2,
                 1,2),
bty = "n", cex = legend_cex)

# Batch 2: Psychonetrics + FGL
legend(0.23, legend_y, legend = c(as.expression(bquote(paste("Partial Pruning, ", alpha, " = 0.05"))),
                                  as.expression(bquote(paste("Partial Pruning, ", alpha, " = 0.01"))),
                                  as.expression(bquote(paste("FGL EBIC, ", gamma, " = 0.25"))),
                                  "FGL CV"
),
col = c(cols[6], cols[6],
        cols[2], cols[2]),
lwd = lwd, lty=c(1,2,
                 1,2),
bty = "n", cex = legend_cex)


# Batch 3: BGGM
legend(0.48, legend_y, legend = c("BGGM, Bayes Factor",
                                  "BGGM, Post Diff (sampling)",
                                  "BGGM, Post Diff (analytical)"),
       col = c(cols[3], cols[3], cols[3]),
       lwd = lwd, lty=c(1,2,3),
       bty = "n", cex = legend_cex)


# Batch 4: Moderation
legend(.755, legend_y, legend = c(#"MNM BIC",
  as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + AND"))),
  "MNM CV + AND",
  #"MNM BIC + OR",
  as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + OR"))),
  "MNM CV + OR"),
  col = c(cols[4], cols[4], cols[4], cols[4]),
  lwd = lwd, lty=c(1,2,3,4),
  bty = "n", cex = legend_cex)


dev.off()





# --------------------------------------------------------------
# ---------- Plot Figure: Timing -------------------------------
# --------------------------------------------------------------

method_labels <- c("NCT1", "NCT2", 
                   "FGL_bic", "FGL_ebic", "FGL_cv",
                   "BGGM_BF", "BGGM_P_an", "BGGM_P_sa",
                   "mgm_bic", "mgm_ebic", "mgm_cv", "mgm_bic_or", "mgm_ebic_or", "mgm_cv_or", 
                   "Fisher1", "Fisher2",
                   "psychonetrics_al1", "psychonetrics_al2")

dim(a_time)
m_timing <- apply(a_time, 1:2, mean)
colnames(m_timing) <- method_labels

m_timing  <- m_timing / 60

library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
cols <- cols[-c(6:7)]


sc <- 1.1
pdf("Simulation_GGM/figures/Figure_Timing_GGM_cmb.pdf", width = 2*6*sc, height=4.5*sc)

par(mfrow=c(1,2))

plot.new()
plot.window(xlim=c(1,10), ylim=c(0, 20))
axis(1, n_seq, at=1:10)
axis(2, las=2)

lwd <- 2
lines(m_timing[, 2], col=cols[1], lwd=lwd) # NCT
lines(m_timing[, 10], col=cols[4], lwd=lwd) # mgm with (E)BIC
lines(m_timing[, 11], col=cols[4], lty=2, lwd=lwd) # mgm with CV

lines(m_timing[, 6], col=cols[3], lwd=lwd) # BGGM: BF
lines(m_timing[, 7], col=cols[3], lty=2, lwd=lwd) # BGGM: Posterior: analytical
lines(m_timing[, 8], col=cols[3], lty=3, lwd=lwd) # BGGM: Posterior: sampling


lines(m_timing[, 15], col=cols[5], lty=1, lwd=lwd) # Fisher

lines(m_timing[, 4], col=cols[2], lty=1, lwd=lwd) # FGL: (E)BIC
lines(m_timing[, 5], col=cols[2], lty=2, lwd=lwd) # FGL: CV

lines(m_timing[, 17], col=cols[6], lty=1, lwd=lwd) # FGL: CV

title(xlab="Sample size per Group", ylab="Runtime (minutes)")
mtext(text = "Running times: GGM", line=1.5, cex=1.3)

legend("topleft", 
       legend = c("NCT", 
                  "MNM with (E)BIC", 
                  "MNM with CV", 
                  "BGGM with BF", 
                  "BGGM with Post. Diff. (analytical)",
                  "BGGM with Post. Diff. (sampling)",
                  "Fisher", 
                  "FGL with (E)BIC", 
                  "FGL with CV", 
                  "Partial Pruning"), 
       col=cols[c(1,4,4,3,3,3,5,2,2,6)], 
       lty=c(1,1,2,1,2,3,1,1,2,1), bty="n", lwd=rep(lwd, 9))


# dev.off()
# sc <- 1.1
# pdf("Simulation_GGM/figures/Figure_Timing_GGM_noFGL.pdf", width = 6*sc, height=4.5*sc)

plot.new()
plot.window(xlim=c(1,10), ylim=c(0, 1))
axis(1, n_seq, at=1:10)
axis(2, las=2)

lwd <- 2
lines(m_timing[, 2], col=cols[1], lwd=lwd) # NCT
lines(m_timing[, 10], col=cols[4], lwd=lwd) # mgm with (E)BIC
lines(m_timing[, 11], col=cols[4], lty=2, lwd=lwd) # mgm with CV

lines(m_timing[, 6], col=cols[3], lwd=lwd) # BGGM: BF
lines(m_timing[, 7], col=cols[3], lty=2, lwd=lwd) # BGGM: Posterior: analytical
lines(m_timing[, 8], col=cols[3], lty=3, lwd=lwd) # BGGM: Posterior: sampling

lines(m_timing[, 15], col=cols[5], lty=1, lwd=lwd) # Fisher

lines(m_timing[, 4], col=cols[2], lty=1, lwd=lwd) # FGL: (E)BIC
# lines(m_timing[, 5], col=cols[2], lty=2, lwd=lwd) # FGL: CV

lines(m_timing[, 17], col=cols[6], lty=1, lwd=lwd) # psychonetrics

title(xlab="Sample size per Group", ylab="Runtime (minutes)")
mtext(text = "Running times: GGM (without FGL CV)", line=1.5, cex=1.3)

# legend("topleft", 
#        legend = c("NCT", 
#                   "MNM with (E)BIC", 
#                   "MNM with CV", 
#                   "BGGM with BF", 
#                   "BGGM with Post. Diff. (analytical)",
#                   "BGGM with Post. Diff. (sampling)",
#                   "Fisher", 
#                   "FGL with (E)BIC", 
#                   "Partial Pruning"), 
#        col=cols[c(1,4,4,3,3,3,5,2,2,6)], 
#        lty=c(1,1,2,1,2,3,1,1,2,1), bty="n", lwd=rep(lwd, 9))


dev.off()







