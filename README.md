# Reproducibility Archive

The files in this reproducibility archive allow one to reproduce the simulation results and the tutorial in the paper "Estimating Group Differences in Network Models using Moderation Analysis", preprint: https://psyarxiv.com/926pv.


### Simulation_GGM

This folder contains files that allow one to reproduce the simulation study on estimating group differences in the GGM.

- `simulation.R` contains a simulation script parallelized over 12 cores, which runs a single iteration of the simulation
- `aux_functions.R` contains a number of helper functions that are sourced in `simulation.R`
- `submit_all.sh` and `submit_jobs.sh` are bash scripts that send `simulation.R` with different seeds to nodes of a cluster computer. However, the simulation could in principle be reproduced by running `simulation.R` with seeds 1, 2, ..., 200 on a local computer. Each node took 3-4 hours to run (parallelized on 12 cores)
- The folder /output contains the 200 output files that are produced by the above three files
- `evaluation.R` preprocesses the output in /output and produces Figure 1 in the paper
- The folder /figures contains the two figures plotted in `evaluation.R`


### Simulation_Ising

This folder contains files that allow one to reproduce the simulation study on estimating group differences in the Ising model.

- `simulation.R` contains a simulation script parallelized over 12 cores, which runs a single iteration of the simulation
- `aux_functions.R` contains a number of helper functions that are sourced in `simulation.R`
- `submit_all.sh` and `submit_jobs.sh` are bash scripts as above
- The folder /output contains the 200 output files that are produced by the above three files
- `evaluation.R` preprocesses the output in /output and produces Figure 2 in the paper
- The folder /figures contains the figure plotted in `evaluation.R`


### Additional_Analyses

#### Sparsity of Graph

- `sim_sparsityG.R` contains a single iteration of the simulation study analogous to the GGM simulation script above, except that the edge-probability of the base-graph is varied instead of the size of the group differences
- `submit_all.sh` and `submit_jobs.sh` are bash scripts as above
- `evaluation.R` creates the figure for this analysis

#### NCT Iterations

- `NCT_sim.R` contains a single iteration of a new simulation studying the extent to which the performance of the NCT depends on the number of iterations
- `submit_all.sh` and `submit_jobs.sh` are bash scripts as above
- `evaluation.R` creates the figure for this analysis

### Tutorial

- `Tutorial.R` contains the code to reproduce the tutorial shown in Section 4.

### SessionInfo()

Here is the session info from the environment in which all simulations were run:

sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS:   /sara/eb/AVX2/Debian10/EB_production/2020/software/R/4.0.2-intel-2020a/lib/R/lib/libR.so
LAPACK: /sara/eb/AVX2/Debian10/EB_production/2020/software/R/4.0.2-intel-2020a/lib/R/modules/lapack.so

locale:
 [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
 [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=en_US   
 [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] dplyr_0.8.5                 psychonetrics_0.9          
 [3] EstimateGroupNetwork_0.3.1  MASS_7.3-51.6              
 [5] doParallel_1.0.15           iterators_1.0.12           
 [7] foreach_1.5.0               BGGM_2.0.3                 
 [9] mgm_1.2-11                  NetworkComparisonTest_2.2.1

loaded via a namespace (and not attached):
  [1] minqa_1.2.4          colorspace_1.4-1     rjson_0.2.20        
  [4] ellipsis_0.3.0       ggridges_0.5.2       htmlTable_1.13.3    
  [7] corpcor_1.6.9        base64enc_0.1-3      rstudioapi_0.11     
 [10] lavaan_0.6-5         IsingFit_0.3.1       fansi_0.4.1         
 [13] mvtnorm_1.1-0        codetools_0.2-16     splines_4.0.2       
 [16] mnormt_1.5-6         knitr_1.28           glasso_1.11         
 [19] Formula_1.2-3        nloptr_1.2.2.1       cluster_2.1.0       
 [22] png_0.1-7            compiler_4.0.2       backports_1.1.6     
 [25] assertthat_0.2.1     Matrix_1.2-18        cli_2.0.2           
 [28] acepack_1.4.1        htmltools_0.4.0      tools_4.0.2         
 [31] igraph_1.2.5         coda_0.19-3          gtable_0.3.0        
 [34] glue_1.4.0           reshape2_1.4.4       Rcpp_1.0.4.6        
 [37] GA_3.2               statnet.common_4.3.0 vctrs_0.2.4         
 [40] nlme_3.1-147         gbRd_0.4-11          psych_1.9.12.31     
 [43] xfun_0.13            stringr_1.4.0        network_1.16.0      
 [46] lme4_1.1-23          lifecycle_0.2.0      gtools_3.8.2        
 [49] statmod_1.4.34       scales_1.1.0         BDgraph_2.62        
 [52] huge_1.3.4.1         RColorBrewer_1.1-2   BFpack_0.3.2        
 [55] VCA_1.4.3            pbapply_1.4-2        gridExtra_2.3       
 [58] ggplot2_3.3.0        IsingSampler_0.2.1   rpart_4.1-15        
 [61] reshape_0.8.8        latticeExtra_0.6-29  stringi_1.4.6       
 [64] checkmate_2.0.0      optimx_2020-4.2      boot_1.3-25         
 [67] bibtex_0.4.2.2       shape_1.4.4          Rdpack_0.11-1       
 [70] rlang_0.4.5          pkgconfig_2.0.3      d3Network_0.5.2.1   
 [73] pracma_2.2.9         lattice_0.20-41      purrr_0.3.4         
 [76] htmlwidgets_1.5.1    tidyselect_1.0.0     GGally_1.5.0        
 [79] plyr_1.8.6           magrittr_1.5         R6_2.4.1            
 [82] Hmisc_4.4-0          combinat_0.0-8       sna_2.5             
 [85] mgcv_1.8-31          pillar_1.4.3         whisker_0.4         
 [88] foreign_0.8-79       survival_3.1-12      abind_1.4-5         
 [91] nnet_7.3-14          tibble_3.0.1         crayon_1.3.4        
 [94] fdrtool_1.2.15       jpeg_0.1-8.1         grid_4.0.2          
 [97] qgraph_1.6.5         data.table_1.12.8    pbivnorm_0.6.0      
[100] bain_0.2.4           matrixcalc_1.0-3     digest_0.6.25       
[103] numDeriv_2016.8-1.1  tidyr_1.0.2          extraDistr_1.9.1    
[106] stats4_4.0.2         munsell_0.5.0        glmnet_3.0-2    