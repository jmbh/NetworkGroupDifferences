# jonashaslbeck@gmail.com; Updated March 9, 2021

# ------------------------------------------------------------------
# ----------- Psychonetrics Wrapper: GGM ---------------------------
# ------------------------------------------------------------------

psychonetricsGGM <- function(data_0, data_1, alpha =  0.05) {
  
  # Code below in this function is from Sacha Epskamp; updated March 2021
  
  vars <- paste0("V",1:ncol(data_1))
  
  data_1 <- as.data.frame(data_1)
  names(data_1) <- vars
  data_1$group <- 1
  
  data_0 <- as.data.frame(data_0)
  names(data_0) <- vars
  data_0$group <- 0
  
  data <- rbind(data_0, data_1)
  
  
  mod <- ggm(data, group = "group", vars = vars, standardize = "z") %>% runmodel
  
  # Pruned:
  mod_prune <-  mod %>%
    partialprune(alpha = alpha)
  
  pcordiffs <- mod_prune@modelmatrices[[1]]$omega - mod_prune@modelmatrices[[2]]$omega
  
  return(pcordiffs)
  
} # eoF


# --------------------------------------------------------------
# ---------- Function to extract differences from mgm object ---
# --------------------------------------------------------------

f_Diffs_from_mgm <- function(m_3way) {
  
  p <- 17
  
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
  G_diff <- ints1[-(p+1), -(p+1)] - ints2[-(p+1), -(p+1)]
  
  
  outlist <- list("diff_logic" = mgm_diff,
                  "diff_value" = G_diff)
  
} # eoF




# --------------------------------------------------------------
# ---------- MLE hypothesis test on differences ----------------
# --------------------------------------------------------------

mle_fisher <- function(x, y, alpha = 0.01){
  
  # code in this function from Donny Williams
  
  # scale
  x <- scale(x, scale = F)
  y <- scale(y, scale = F)
  n <- nrow(x)
  p <- ncol(x)
  
  # prep
  k  <- p - 1
  mat <- mat_values <- matrix(0, p, p)
  se_diff <- sqrt((1 / (n - k - 3)) + (1 / (n - k - 3)))
  
  # get partial correlations
  x_pcs  <- cov2cor(solve(nrow(x)^-1 *   t(x) %*% x)) * -1
  y_pcs  <- cov2cor(solve(nrow(y)^-1 *   t(y) %*% y)) * -1
  
  # fishers z-transform
  pcs_diff <- BGGM:::fisher_z(x_pcs[upper.tri(x_pcs)]) - BGGM:::fisher_z(y_pcs[upper.tri(y_pcs)])
  z_score <- abs(pcs_diff / se_diff)
  
  # threshold at alpha level
  pvalues <- (1 - pnorm(z_score)) * 2
  
  # save results in matrices
  mat[upper.tri(mat)] <- ifelse(pvalues < alpha, 1, 0)
  mat_values[upper.tri(mat_values)] <- pcs_diff
  mat_values <- mat_values + t(mat_values)
  diff_logic <- mat + t(mat)
  mat_values[diff_logic!=1] <- 0
  
  outlist <- list("diff_logic" = diff_logic, 
                  "diff_value" = mat_values)
  
  return(outlist)
  
} # eoF


# --------------------------------------------------------------
# ---------- Aux function: Generate True Models ----------------
# --------------------------------------------------------------


f_GenModel <- function(iter, p = 17, Pe = 0.2, delta = .2) {
  
  set.seed(iter) 
  
  n_edges <- p*(p-1)/2
  
  okay <- FALSE
  counter <- 0
  while(okay == FALSE) {
    
    G <- matrix(0, p, p)
    G[upper.tri(G)] <- sample(0:1, size = n_edges, prob = c(1-Pe, Pe), replace = TRUE)
    G[upper.tri(G)] <- G[upper.tri(G)] * runif(n_edges, -.2, .4)
    G <- G + t(G)
    
    # Let's say 20 are different; always delta=.2; 10 positive sign, 10 negative sign
    G_diff <- matrix(0, p, p)
    G_diff[upper.tri(G_diff)] <- sample(1:n_edges, size = n_edges, replace = FALSE)
    
    G_diff[G_diff %in% 1:20] <- delta
    G_diff[G_diff %in% 21:n_edges] <- 0
    
    G_diff <- G_diff + t(G_diff)
    
    # Get the two groups
    Sigma1_inv <- G
    Sigma2_inv <- G + G_diff
    
    diag(Sigma1_inv) <- 1
    diag(Sigma2_inv) <- 1
    
    # Going to correlation
    # Sig1
    m1 = -Sigma1_inv
    diag(m1) = -diag(m1)
    m1 = solve(m1)
    eig1 <- eigen(m1)
    
    # Sig1
    m2 = -Sigma2_inv
    diag(m2) = -diag(m2)
    m2 = solve(m2)
    eig2 <- eigen(m2)
    
    
    if(all(eig1$values > 0) & all(eig2$values > 0)) {
      
      okay <- TRUE
      
      Sigma1 <- cov2cor(m1)
      Sigma2 <- cov2cor(m2)
      
    } else {
      counter <- counter + 1 # to keep track of how long resampling was
    }
    
  } # end: while: finding two pos def covariances
  
  outlist <- list("Sigma1"=Sigma1, 
                  "Sigma2"=Sigma2, 
                  "counter"=counter, 
                  "G_diff"=G_diff, 
                  "Sigma1_inv" = Sigma1_inv, 
                  "Sigma2_inv" = Sigma2_inv)
  
  return(outlist)
  
} # eoF

