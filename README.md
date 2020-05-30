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
- `submit_all.sh` and `submit_jobs.sh` are bash scripts that send `simulation.R` with different seeds to nodes of a cluster computer. However, the simulation could in principle be reproduced by running `simulation.R` with seeds 1, 2, ..., 200 on a local computer. Each node took 5-6 hours to run (parallelized on 12 cores)
- The folder /output contains the 200 output files that are produced by the above three files
- `evaluation.R` preprocesses the output in /output and produces Figure 2 in the paper
- The folder /figures contains the figure plotted in `evaluation.R`


### Tutorial

- `Tutorial.R` contains the code to reproduce the tutorial shown in Section 4.

### SessionInfo()

Here is the session info from the environment in which all simulations were run:

R version 3.6.1 (2019-07-05)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Debian GNU/Linux 10 (buster)

Matrix products: default
BLAS/LAPACK: /home/haslbeck/.conda/envs/my_root/lib/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
 [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=en_US   
 [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] doParallel_1.0.15           iterators_1.0.12           
[3] foreach_1.5.0               dplyr_0.8.5                
[5] psychonetrics_0.7.1         BGGM_2.0.0                 
[7] mgm_1.2-9                   EstimateGroupNetwork_0.2.2 
[9] NetworkComparisonTest_2.2.1

loaded via a namespace (and not attached):
  [1] minqa_1.2.4          colorspace_1.4-1     rjson_0.2.20        
  [4] ellipsis_0.3.0       ggridges_0.5.2       htmlTable_1.13.3    
  [7] corpcor_1.6.9        base64enc_0.1-3      rstudioapi_0.11     
 [10] lavaan_0.6-5         IsingFit_0.3.1       fansi_0.4.1         
 [13] mvtnorm_1.1-0        codetools_0.2-16     splines_3.6.1       
 [16] mnormt_1.5-7         knitr_1.28           glasso_1.11         
 [19] Formula_1.2-3        nloptr_1.2.2.1       cluster_2.0.8       
 [22] png_0.1-7            compiler_3.6.1       backports_1.1.7     
 [25] assertthat_0.2.1     Matrix_1.2-17        cli_2.0.2           
 [28] acepack_1.4.1        htmltools_0.4.0      tools_3.6.1         
 [31] igraph_1.2.5         coda_0.19-3          gtable_0.3.0        
 [34] glue_1.4.1           reshape2_1.4.4       Rcpp_1.0.4.6        
 [37] GA_3.2               statnet.common_4.3.0 vctrs_0.3.0         
 [40] nlme_3.1-139         gbRd_0.4-11          psych_1.9.12.31     
 [43] xfun_0.13            stringr_1.4.0        network_1.16.0      
 [46] lme4_1.1-23          lifecycle_0.2.0      gtools_3.8.2        
 [49] statmod_1.4.34       MASS_7.3-51.3        scales_1.1.1        
 [52] BDgraph_2.62         huge_1.3.4.1         RColorBrewer_1.1-2  
 [55] BFpack_0.2.1         VCA_1.4.2            pbapply_1.4-2       
 [58] gridExtra_2.3        ggplot2_3.3.0        IsingSampler_0.2.1  
 [61] rpart_4.1-15         reshape_0.8.8        latticeExtra_0.6-29 
 [64] stringi_1.4.6        ucminf_1.1-4         checkmate_2.0.0     
 [67] optimx_2020-4.2      boot_1.3-20          bibtex_0.4.2.2      
 [70] shape_1.4.4          Rdpack_0.11-1        rlang_0.4.6         
 [73] pkgconfig_2.0.3      d3Network_0.5.2.1    pracma_2.2.9        
 [76] lattice_0.20-38      purrr_0.3.4          htmlwidgets_1.5.1   
 [79] tidyselect_1.1.0     GGally_1.5.0         plyr_1.8.6          
 [82] magrittr_1.5         R6_2.4.1             Hmisc_4.4-0         
 [85] combinat_0.0-8       sna_2.5              mgcv_1.8-28         
 [88] pillar_1.4.4         whisker_0.4          foreign_0.8-71      
 [91] survival_3.1-12      abind_1.4-5          nnet_7.3-12         
 [94] tibble_3.0.1         crayon_1.3.4         fdrtool_1.2.15      
 [97] jpeg_0.1-8.1         grid_3.6.1           qgraph_1.6.5        
[100] data.table_1.12.8    pbivnorm_0.6.0       bain_0.2.4          
[103] matrixcalc_1.0-3     digest_0.6.25        tidyr_1.0.3         
[106] numDeriv_2016.8-1.1  stats4_3.6.1         munsell_0.5.0       
[109] glmnet_4.0     

