# jonashaslbeck@gmail.com; May 2020

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
# ---------- Plotting Figure 6 in Appendix ---------------------
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

