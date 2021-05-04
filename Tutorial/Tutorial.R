# jonashaslbeck@gmail.com; May 2020

# ---------------------------------------------------------------------------------------
# ---------- Basic Data Info ------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# library(devtools)
# install_github("jmbh/mgm") # version 1.2-10
library(mgm)

data(dataGD)
dim(dataGD)
head(dataGD)

library(qgraph)

# ---------------------------------------------------------------------------------------
# ---------- Fit moderated MGM ----------------------------------------------------------
# ---------------------------------------------------------------------------------------

mgm_obj <- mgm(data = dataGD, 
               type = c("g", "g", "c", "g", "c", "g", "c"), 
               level = c(1, 1, 2, 1, 3, 1, 3), 
               moderators = 7, 
               lambdaSel = "EBIC", 
               lambdaGam = 0.25, 
               ruleReg = "AND")


# ---------------------------------------------------------------------------------------
# ---------- Condition mgm object on moderator (group variable) --------------------------
# ---------------------------------------------------------------------------------------

l_mgm_cond <- list()
for(g in 1:3) l_mgm_cond[[g]] <- condition(mgm_obj, values = list("7"=g))


# ---------------------------------------------------------------------------------------
# ---------- Visualize pairwise model for each group -------------------------------------
# ---------------------------------------------------------------------------------------

# Compute Maximum for each group
v_max <- rep(NA, 3)
for(g in 1:3) v_max[g] <- max(l_mgm_cond[[g]]$pairwise$wadj)

pdf("Tutorial/Tutorial_pw_group.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
for(g in 1:3) {
  qgraph(input = l_mgm_cond[[g]]$pairwise$wadj, 
         edge.color = l_mgm_cond[[g]]$pairwise$edgecolor, 
         layout = "circle", mar=c(2,3,5,3),
         maximum = max(v_max), vsize=16, esize=23, 
         edge.labels  = TRUE, edge.label.cex = 3)
  mtext(text = paste0("Group ", g), line=2.5)
}
dev.off()


# ---------------------------------------------------------------------------------------
# ---------- Inspect Categorical Interactions ------------------------------------------
# ---------------------------------------------------------------------------------------

int35_g1 <- showInteraction(l_mgm_cond[[3]], int = c(3,5))
int35_g1$parameters


# ---------------------------------------------------------------------------------------
# ---------- Make Difference Graph for Appendix -----------------------------------------
# ---------------------------------------------------------------------------------------

# Get adjacency matrices with signs (if defined)
l_sadj <- list()
for(g in 1:3) {
  sadj_g <- l_mgm_cond[[g]]$pairwise$wadj
  sadj_g[l_mgm_cond[[g]]$pairwise$edgecolor=="red"] <- sadj_g[l_mgm_cond[[g]]$pairwise$edgecolor=="red"] * -1
  l_sadj[[g]] <- sadj_g
}

# Make differences
GDiff_12 <- l_sadj[[1]] - l_sadj[[2]]
GDiff_13 <- l_sadj[[1]] - l_sadj[[3]]
GDiff_23 <- l_sadj[[2]] - l_sadj[[3]]

sc <- .8
pdf("Tutorial/Tutorial_Diffnetworks.pdf", width = 7*sc, height = 7*sc)

# Make layout
lmat <- rbind(c(0,1,2), 
              c(3,5,6),
              c(4,0,7))
lo <- layout(mat = lmat, widths = c(.1,1,1), heights = c(.1,1,1))
# layout.show(lo)

# Plot Labels
PlotLabels <- function(label, srt=0, cex=1) {
  
  par(mar=rep(0,4))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1,1))
  
  text(0, 0, label=label, cex=cex, srt=srt)
}

cex_fix <- 1.3
PlotLabels("Group 2", cex=cex_fix)
PlotLabels("Group 3", cex=cex_fix)
PlotLabels("Group 1", cex=cex_fix, srt=90)
PlotLabels("Group 2", cex=cex_fix, srt=90)

# Plot difference networks
maximum <- max(abs(GDiff_13), abs(GDiff_12), abs(GDiff_23))

mar_fix <- 5
qgraph(GDiff_12, mar=rep(mar_fix,4), vsize=13, esize=21)
qgraph(GDiff_13, mar=rep(mar_fix,4), vsize=13, esize=21)
qgraph(GDiff_23, mar=rep(mar_fix,4), vsize=13, esize=21)

dev.off()




