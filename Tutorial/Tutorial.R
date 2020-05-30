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

pdf("Tutorial_pw_group.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
library(qgraph)
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


