# jonashaslbeck@gmail.com; May 4, 2021

# --------------------------------------------------------------
# ---------- Load Simulation Results ---------------------------
# --------------------------------------------------------------

# Load files
simResDir <- "Additional_Analyses/NCT_Iterations/output/"

v_files <- list.files(simResDir)
n_files <- length(v_files)
l_files <- list()
for(i in 1:n_files) l_files[[i]] <- readRDS(paste0(simResDir, v_files[i]))

# Collect in one object
a_out <- array(NA, dim=c(n_files, 11, 4))

for(i in 1:n_files) {
  for(j in 1:11) {
    a_out[i, j, ] <- l_files[[i]][[j]][j, ]
  }
}

# Aggregate over iterations
m_out <- apply(a_out, 2:3, function(x) mean(x, na.rm = TRUE))



# --------------------------------------------------------------
# ---------- Make Figure ---------------------------------------
# --------------------------------------------------------------

iter_seq <- round(seq(50, 1000, length=11))


pdf("Additional_Analyses/NCT_Iterations/figures/Fig_App_NCTiterationss.pdf", width = 8, height = 4)

par(mfrow=c(1,2))

# Sensitivity and Precision
plot.new()
plot.window(xlim=range(iter_seq), ylim=c(0,1))
axis(1, iter_seq)
axis(2, las=2)
title(xlab="Number of NCT iterations")

lines(iter_seq, m_out[, 1], lwd=2)
lines(iter_seq, m_out[, 2], lty=2, lwd=2)

legend("topright", 
       legend=c("Sensitivity", "Precision"), 
       lty=1:2, lwd=c(2,2), bty="n", cex=0.9)

# Estimation Errors absend & present
plot.new()
plot.window(xlim=range(iter_seq), ylim=c(0,.3))
axis(1, iter_seq)
axis(2, las=2)
title(xlab="Number of NCT iterations")

lines(iter_seq, m_out[, 3], lwd=2)
lines(iter_seq, m_out[, 4], lty=2, lwd=2)

legend("right", 
       legend=c("Abs. Estimation Error (present)", "Abs. Estimation Error (absent)"), 
       lty=1:2, lwd=c(2,2), bty="n", cex=0.9)

dev.off()











