# jonashaslbeck@gmail.com; March 12, 2021

# ------------------------------------------------------------------
# ----------- Psychonetrics Wrapper: Ising model -------------------
# ------------------------------------------------------------------

psychonetricsIsing <- function(data_0, data_1, alpha =  0.05) {
  
  # Code below is from Sacha Epskamp
  
  vars <- paste0("V",1:ncol(data_1))
  
  data_1 <- as.data.frame(data_1)
  names(data_1) <- vars
  data_1$group <- 1
  
  data_0 <- as.data.frame(data_0)
  names(data_0) <- vars
  data_0$group <- 0
  
  data <- rbind(data_0, data_1)
  
  mod <- Ising(data, group = "group", vars = vars) %>% runmodel
  
  # Pruned:
  mod_prune <-  mod %>%
    prune(alpha = alpha)
  
  mod_union <- mod_prune%>%
    unionmodel %>% # Include edges in each model that were significant in any model
    groupequal("omega")  %>% # Constrain omega to be equal
    runmodel # Run the model
  
  # Now, iterativly free parameters along modification indices
  curMod <- mod_union
  
  repeat{
    pars <- curMod@parameters
    
    # Make a data frame in which the equality free parameters are summed:
    miDF <- pars %>% filter(!fixed) %>% group_by(group,row,col,matrix) %>%
      summarize(mi_free = sum(mi_free)) %>% 
      arrange(-mi_free)
    
    # Free the best parameter:
    propMod <- curMod %>% 
      groupfree(miDF$matrix[1],miDF$row[1],miDF$col[1]) %>% 
      runmodel
    
    # Test BIC:
    if (propMod@fitmeasures$bic < curMod@fitmeasures$bic){
      curMod <- propMod
    } else {
      break
    }
  }
  mod_partialpooled <- curMod %>% prune(alpha = alpha)
  
  # Select best model:
  mods <- list(mod, mod_prune, mod_union, mod_partialpooled)
  
  best <- which.min(c(
    mod@fitmeasures$bic,
    mod_prune@fitmeasures$bic,
    mod_union@fitmeasures$bic,
    mod_partialpooled@fitmeasures$bic
  ))
  
  bestmod <- mods[[best]]
  
  pcordiffs <- bestmod@modelmatrices[[1]]$omega - bestmod@modelmatrices[[2]]$omega
  
  return(pcordiffs)
  
} # eoF


# --------------------------------------------------------------
# ---------- Function to generate true Ising models ------------
# --------------------------------------------------------------

f_GenModel_Ising <- function(iter, 
                             p = 17, 
                             Pe = 0.2, 
                             delta = 1) {
  
  set.seed(iter) 
  
  n_edges <- p*(p-1)/2
  
  # Step 1: the base graph
  G <- matrix(0, p, p)
  G[upper.tri(G)] <- sample(0:1, size = n_edges, prob = c(1-Pe, Pe), replace = TRUE)
  G[upper.tri(G)] <- G[upper.tri(G)] * runif(n_edges, -1, 2)
  G <- G + t(G)
  
  # Step 2: generating differences
  # Let's say 20 are different; always delta=.2; 10 positive sign, 10 negative sign
  G_diff <- matrix(0, p, p)
  G_diff[upper.tri(G_diff)] <- sample(1:n_edges, size = n_edges, replace = FALSE)
  G_diff[G_diff %in% 1:20] <- delta
  G_diff[G_diff %in% 21:n_edges] <- 0
  G_diff <- G_diff + t(G_diff)
  
  # Get the two groups
  Graph1 <- G
  Graph2 <- G + G_diff
  
  thresh1 <- -colSums(Graph1) / 2
  thresh2 <- -colSums(Graph2) / 2
  
  # Define outlist  
  outlist <- list("Graph1"=Graph1, 
                  "Graph2"=Graph2, 
                  "thresh1"=thresh1,
                  "thresh2"=thresh2,
                  "G_diff"=G_diff)
  
  return(outlist)
  
} # eoF



# --------------------------------------------------------------
# ---------- Function to extract differences from mgm object ---
# --------------------------------------------------------------


f_Diffs_from_mgm <- function(MNM_obj) {
  
  p <- 17
  
  # Get 3-way interactions
  m_3way <- MNM_obj$interactions$indicator[[2]]
  m_3way <- matrix(m_3way, ncol=3)
  
  # Compute edges that are different
  mgm_diff <- matrix(0, p, p)
  if(nrow(m_3way) > 0) {
    m_3way <- as.matrix(m_3way)
    for(i in 1:nrow(m_3way)) {
      int_i <- m_3way[i, ]
      if(int_i[3] != 18) next else mgm_diff[int_i[1], int_i[2]] <- mgm_diff[int_i[2], int_i[1]] <- 1
    }
  }
  
  # Compute differences between edges as function of moderator variable 18
  # Get two group networks
  obj1 <- condition(MNM_obj, values = list("18"=0))
  ints1 <- obj1$pairwise$wadj
  ints1[obj1$pairwise$edgecolor == "red"] <- ints1[obj1$pairwise$edgecolor == "red"] * -1
  
  obj2 <- condition(MNM_obj, values = list("18"=1))
  ints2 <- obj2$pairwise$wadj
  ints2[obj2$pairwise$edgecolor == "red"] <- ints2[obj2$pairwise$edgecolor == "red"] * -1
  
  # Difference
  multi <- 2 # multiply with 2 to get to parameter scale from Ising Sampler
  G_1_wo_group <- ints1[-(p+1), -(p+1)] * multi # multiply with 2 to get to parameter scale from Ising Sampler
  G_2_wo_group <- ints2[-(p+1), -(p+1)] * multi
  
  G_diff <- G_1_wo_group - G_2_wo_group
  
  
  outlist <- list("diff_logic" = mgm_diff,
                  "diff_value" = G_diff)
  
} # eoF