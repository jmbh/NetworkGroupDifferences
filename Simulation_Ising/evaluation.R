# jonashaslbeck@gmail.com; May 4, 2021

# --------------------------------------------------------------
# ---------- Load Data -----------------------------------------
# --------------------------------------------------------------

simResDir <- "Simulation_Ising/output/"

v_files <- list.files(simResDir)
n_files <- length(v_files)

l_files <- list()

i_files <- 1:n_files
nIter <- length(i_files)

for(i in i_files) l_files[[i]] <- readRDS(paste0(simResDir, v_files[i]))

# Create overall output arrays:
n_methods <- 10
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
# ---------- Preprocess ----------------------------------------
# --------------------------------------------------------------

# Sensitivity, Precision, Estimation Error; probably: average estimation error (average across detected differences)

a_P <- array(NA, dim = c(5, n_methods, 3, 10, nIter))

for(m in 1:n_methods) { # method
  for(n in 1:10) { # sample size
    for(d in 1:3) { # delta theta
      for(i in 1:nIter) { # iterations
        
        # Sensitivity and Precision
        true <- a_out[, , 1, n_methods+1, d, n, i]
        a_P[1, m, d, n, i] <- mean(a_out[, , 1, m, d, n, i][true==1]) # Sensitivity
        a_P[2, m, d, n, i] <- mean(true[a_out[, , 1, m, d, n, i]!=0]) # Precision
        
        # Absolute Estimation Error (among true differences)
        ## Estimation Errors
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


# --------------------------------------------------------------
# ---------- Plot Figure 2 -------------------------------------
# --------------------------------------------------------------

# n-variations (from simulation)
n_seq <- round(exp(seq(3, 8.5171, length=10)))

library(RColorBrewer)
cols <- brewer.pal(5, "Set1")
lwd=1.5
title.cex <- 1

sc <- 2
pdf("Simulation_Ising/figures/Fig_Sim_Ising.pdf", width = 4.3*sc, height = 3.5*sc)

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
v_v_titles[[1]] <- bquote(paste(Delta, theta, " = 0.15"))
v_v_titles[[2]] <- bquote(paste(Delta, theta, " = 0.30"))
v_v_titles[[3]] <- bquote(paste(Delta, theta, " = 0.60"))
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

for(d in 1:3) {
  
  # ----- Plotting Sensitivity and Precision ------
  for(s in 1:4) {
    
    if(s < 3) {
      # Set up plotting area
      plot.new()
      plot.window(ylim=c(0,1), xlim=c(1,10))
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, seq(0, 1, length=5), las=2)
    } else {
      # Set up plotting area
      plot.new()
      plot.window(ylim=c(0, 1), xlim=c(1,10))
      axis(1, at = 1:10, labels = n_seq, las=2)
      axis(2, seq(0, 1, length=5), las=2)
    }
    
    # Plot Data: 
    # NCT
    lines(a_P_means[s, 1, d, ], col=cols[1], lty=1, lwd=lwd) # alpha = 0.05
    lines(a_P_means[s, 2, d, ], col=cols[1], lty=2, lwd=lwd) # alpha = 0.01
    
    # MGM
    # lines(a_P_means[s, 3, d, ], col=cols[4], lty=1, lwd=lwd) # BIC + AND
    lines(a_P_means[s, 4, d, ], col=cols[4], lty=1, lwd=lwd) # EBIC + AND
    lines(a_P_means[s, 5, d, ], col=cols[4], lty=2, lwd=lwd) # CV + AND
    
    # lines(a_P_means[s, 6, d, ], col=cols[4], lty=4, lwd=lwd) # BIC + OR
    lines(a_P_means[s, 7, d, ], col=cols[4], lty=3, lwd=lwd) # EBIC + OR
    lines(a_P_means[s, 8, d, ], col=cols[4], lty=4, lwd=lwd) # CV + OR
    
    # BGGM
    lines(a_P_means[s, 9, d, ], col=cols[3], lty=1, lwd=lwd) # Bayes factor
    lines(a_P_means[s, 10, d, ], col=cols[3], lty=2, lwd=lwd) # Post Diff
    
    
  } # for: s
  
} # for: d


# Plot legend in lower panel (no 26)
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))

legend_cex <- 1.15
legend_y <- 1

# Batch 1: NCT
legend(0.1, legend_y, legend = c(as.expression(bquote(paste("NCT, ", alpha, " = 0.05"))), 
                                 as.expression(bquote(paste("NCT, ", alpha, " = 0.01")))),
       col = c(cols[1], cols[1], cols[4]),
       lwd = lwd, lty=c(1,2, 1), 
       bty = "n", cex= legend_cex)

# Batch 2: BGGM
legend(0.4, legend_y, legend = c("BGGM, Bayes factor",
                                 "BGGM, Post diff"),
       col = c(cols[3], cols[3]),
       lwd = lwd, lty=c(1,2), 
       bty = "n", cex= legend_cex)

# Batch 2: MNM
legend(0.7, legend_y, legend = c(as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + AND"))),
                                 "MNM CV + AND", 
                                 as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + OR"))), 
                                 "MNM CV + OR"), 
       col = c(cols[4],cols[4],cols[4],cols[4]), 
       lwd = lwd, lty=c(1:4), 
       bty = "n", cex= legend_cex)

dev.off()





# --------------------------------------------------------------
# ---------- Marginal Performance ------------------------------
# --------------------------------------------------------------

# Marginalize over n and delta theta

method_labels <- c("NCT1", "NCT2", 
                   "mgm_bic", "mgm_ebic", "mgm_cv", "mgm_bic_or", "mgm_ebic_or", "mgm_cv_or", 
                   "BGGM_BF", "BGGM_P") # this is the order in the vector

methods_labels_proper <- c(c(as.expression(bquote(paste("NCT, ", alpha, " = 0.05"))), 
                             as.expression(bquote(paste("NCT, ", alpha, " = 0.01")))),
                           c(as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + AND"))),
                             "MNM CV + AND", 
                             as.expression(bquote(paste("MNM EBIC, ", gamma, " = 0.25 + OR"))), 
                             "MNM CV + OR"), 
                           c("BGGM, Bayes factor",
                             "BGGM, Post diff")) # this is the order in the above figure

methods_select <- c(1,2,4,5,7,8,9,10) # Select same as in graph above
method_labels[methods_select] # Sanity check
n_methods <- length(methods_select)
methods_labels_proper

# Aggregate over N
a_P_means_aggN <- apply(a_P_means, 1:2, function(x) mean(x, na.rm=TRUE))
a_P_qntls_aggN <- apply(a_P_means, 1:2, function(x) quantile(x, na.rm=TRUE, probs = c(.25, .75)))
dim(a_P_means_aggN)

# Overview matrix: est(present)
m_ov <- matrix(NA, n_methods, 2)
m_ov[, 2] <- round(a_P_means_aggN[3, methods_select], 4)
m_ov[, 1] <- method_labels[methods_select]
or_1 <- order(m_ov[, 2])
m_ov_ord <- m_ov[or_1, ]

# Overview matrix: est(absent)
m_ov2 <- matrix(NA, n_methods, 2)
m_ov2[, 2] <- round(a_P_means_aggN[4, methods_select], 4)
m_ov2[, 1] <- method_labels[methods_select]
m_ov2_ord <- m_ov2[order(as.numeric(m_ov2[, 2])), ]


# Proper figure

sc <- 1.2
pdf("Simulation_Ising/figures/Fig_Sim_marginal_Ising.pdf", width=7*sc*0.9, height=4.8*sc)

par(mar=c(9.2,5,2,2))
plot.new()
plot.window(ylim=c(0, .5), xlim=c(1,n_methods))
axis(2, las=2)
axis(1, labels = rep("", n_methods), at=1:n_methods, 
     las=2, cex.axis=.8)
cols_fix <- cols[c(1,1,4,4,4,4,3,3)[or_1]]
# cols_fix <- rep("black", n_methods)
for(i in 1:n_methods) axis(1, labels = methods_labels_proper[or_1][i], at=i, 
                    las=2, cex.axis=.8, col.axis=cols_fix[i])

title(ylab="Estimation error", line=3.5)

## error: present
# mean
points(m_ov[, 2][or_1], pch=20, cex=1.3)
# upper/lower quantile
segments(x0 = 1:n_methods, y0 = a_P_qntls_aggN[1, 3, ][methods_select][or_1], 
         x1 = 1:n_methods, y1 = a_P_qntls_aggN[2, 3, ][methods_select][or_1], 
         lwd=1.5)


## error: absent
# mean
points(m_ov2[, 2][or_1], pch=20, col="grey", cex=1.3)
# upper/lower quantile
segments(x0 = 1:n_methods, y0 = a_P_qntls_aggN[1, 4, ][methods_select][or_1], 
         x1 = 1:n_methods, y1 = a_P_qntls_aggN[2, 4, ][methods_select][or_1], 
         col="grey", 
         lwd=1.5)

legend("topleft", legend=c("Mean error (present)", "Mean error (absent)"), 
       col=c("black", "grey"), pch=c(20,20), bty="n")


dev.off()



# --------------------------------------------------------------
# ---------- Plot Figure: Timing -------------------------------
# --------------------------------------------------------------

method_labels <- c("NCT1", "NCT2", 
                   "mgm_bic", "mgm_ebic", "mgm_cv", "mgm_bic_or", "mgm_ebic_or", "mgm_cv_or", 
                   "BGGM_BF", "BGGM_P")

dim(a_time)
m_timing <- apply(a_time, 1:2, mean)
colnames(m_timing) <- method_labels

m_timing  <- m_timing / 60

library(RColorBrewer)
cols <- brewer.pal(5, "Set1")


sc <- 1.1
pdf("Simulation_Ising/figures/Figure_Timing_Ising.pdf", width = 6*sc, height=4.5*sc)

plot.new()
plot.window(xlim=c(1,10), ylim=c(0, 40))
axis(1, n_seq, at=1:10)
axis(2, las=2)

lwd <- 2
lines(m_timing[, 2], col=cols[1], lwd=lwd) # NCT
lines(m_timing[, 3], col=cols[4], lwd=lwd) # mgm with (E)BIC
lines(m_timing[, 5], col=cols[4], lty=2, lwd=lwd) # mgm with CV
lines(m_timing[, 9], col=cols[3], lwd=lwd) # BGGM: BF
lines(m_timing[, 10], col=cols[3], lty=2, lwd=lwd) # BGGM: Posterior

title(xlab="Sample size per Group", ylab="Runtime (minutes)")
mtext(text = "Running times: Ising Model", line=1.5, cex=1.3)

legend("topleft", 
       legend = c("NCT", 
                  "MNM with (E)BIC", 
                  "MNM with CV", 
                  "BGGM with BF", 
                  "BGGM with Post. Diff."), 
       col=cols[c(1,4,4,3,3)], 
       lty=c(1,1,2,1,2), bty="n", lwd=rep(lwd, 5))


dev.off()







