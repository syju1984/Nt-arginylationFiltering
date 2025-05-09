# Nt-arginylationFiltering

code re-run was successful in the session:
R version 4.5.0 (2025-04-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: Asia/Seoul
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrepel_0.9.6          UniProtKeywords_0.99.7 rrvgo_1.20.0           ReactomePA_1.52.0     
 [5] org.Hs.eg.db_3.21.0    AnnotationDbi_1.70.0   Biobase_2.68.0         clusterProfiler_4.16.0
 [9] ggVennDiagram_1.5.2    viridis_0.6.5          viridisLite_0.4.2      ggpointdensity_0.1.0  
[13] ComplexHeatmap_2.24.0  yardstick_1.3.2        workflowsets_1.1.0     workflows_1.2.0       
[17] tune_1.3.0             rsample_1.3.0          recipes_1.3.0          parsnip_1.3.1         
[21] modeldata_1.4.0        infer_1.0.8            dials_1.4.0            tidymodels_1.3.0      
[25] scales_1.4.0           OrgMassSpecR_0.5-3     broom_1.0.8            conflicted_1.2.0      
[29] Rdisop_1.68.0          Rcpp_1.0.14            wesanderson_0.3.7      dagLogo_1.46.0        
[33] ggseqlogo_0.2          RColorBrewer_1.1-3     ggpubr_0.6.0           Peptides_2.4.6        
[37] Biostrings_2.76.0      GenomeInfoDb_1.44.0    XVector_0.48.0         IRanges_2.42.0        
[41] S4Vectors_0.46.0       BiocGenerics_0.54.0    generics_0.1.3         devtools_2.4.5        
[45] usethis_3.1.0          lubridate_1.9.4        forcats_1.0.0          stringr_1.5.1         
[49] dplyr_1.1.4            purrr_1.0.4            readr_2.1.5            tidyr_1.3.1           
[53] tibble_3.2.1           ggplot2_3.5.2          tidyverse_2.0.0       

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2           progress_1.2.3              urlchecker_1.0.1           
  [4] nnet_7.3-20                 vctrs_0.6.5                 ggtangle_0.0.6             
  [7] digest_0.6.37               png_0.1-8                   shape_1.4.6.1              
 [10] BiocBaseUtils_1.10.0        parallelly_1.44.0           MASS_7.3-65                
 [13] reshape2_1.4.4              httpuv_1.6.16               foreach_1.5.2              
 [16] qvalue_2.40.0               withr_3.0.2                 xfun_0.52                  
 [19] ggfun_0.1.8                 ellipsis_0.3.2              survival_3.8-3             
 [22] memoise_2.0.1               gson_0.1.0                  profvis_0.4.0              
 [25] tidytree_0.4.6              GlobalOptions_0.1.2         gtools_3.9.5               
 [28] R.oo_1.27.1                 Formula_1.2-5               prettyunits_1.2.0          
 [31] KEGGREST_1.48.0             promises_1.3.2              httr_1.4.7                 
 [34] rstatix_0.7.2               restfulr_0.0.15             globals_0.17.0             
 [37] rstudioapi_0.17.1           UCSC.utils_1.4.0            miniUI_0.1.2               
 [40] DOSE_4.2.0                  reactome.db_1.92.0          curl_6.2.2                 
 [43] ncdf4_1.24                  ggraph_2.2.1                polyclip_1.10-7            
 [46] GenomeInfoDbData_1.2.14     SparseArray_1.8.0           xtable_1.8-4               
 [49] ade4_1.7-23                 doParallel_1.0.17           evaluate_1.0.3             
 [52] S4Arrays_1.8.0              BiocFileCache_2.16.0        hms_1.1.3                  
 [55] GenomicRanges_1.60.0        colorspace_2.1-1            filelock_1.0.3             
 [58] NLP_0.3-2                   reticulate_1.42.0           treemap_2.4-4              
 [61] magrittr_2.0.3              later_1.4.2                 ggtree_3.16.0              
 [64] lattice_0.22-7              future.apply_1.11.3         lhs_1.2.0                  
 [67] XML_3.99-0.18               cowplot_1.1.3               matrixStats_1.5.0          
 [70] class_7.3-23                pillar_1.10.2               nlme_3.1-168               
 [73] iterators_1.0.14            pwalign_1.4.0               gridBase_0.4-7             
 [76] GPfit_1.0-9                 caTools_1.18.3              compiler_4.5.0             
 [79] RSpectra_0.16-2             stringi_1.8.7               gower_1.0.2                
 [82] SummarizedExperiment_1.38.1 GenomicAlignments_1.44.0    plyr_1.8.9                 
 [85] crayon_1.5.3                abind_1.4-8                 BiocIO_1.18.0              
 [88] gridGraphics_0.5-1          graphlayouts_1.2.2          bit_4.6.0                  
 [91] fastmatch_1.1-6             codetools_0.2-20            openssl_2.3.2              
 [94] AnVILBase_1.2.0             slam_0.1-55                 GetoptLong_1.0.5           
 [97] tm_0.7-16                   mime_0.13                   splines_4.5.0              
[100] circlize_0.4.16             dbplyr_2.5.0                DiceDesign_1.10            
[103] knitr_1.50                  blob_1.2.4                  clue_0.3-66                
[106] BiocVersion_3.21.1          rjsoncons_1.3.2             seqLogo_1.74.0             
[109] mzR_2.42.0                  fs_1.6.6                    listenv_0.9.1              
[112] pkgbuild_1.4.7              ggsignif_0.6.4              ggplotify_0.1.2            
[115] Matrix_1.7-3                tzdb_0.5.0                  tweenr_2.0.3               
[118] pkgconfig_2.0.3             pheatmap_1.0.12             tools_4.5.0                
[121] cachem_1.1.0                RSQLite_2.3.11              DBI_1.2.3                  
[124] graphite_1.54.0             fastmap_1.2.0               rmarkdown_2.29             
[127] Rsamtools_2.24.0            AnnotationHub_3.16.0        patchwork_1.3.0            
[130] BiocManager_1.30.25         graph_1.86.0                carData_3.0-5              
[133] rpart_4.1.24                farver_2.1.2                tidygraph_1.3.1            
[136] yaml_2.3.10                 MatrixGenerics_1.20.0       rtracklayer_1.68.0         
[139] cli_3.6.5                   lifecycle_1.0.4             askpass_1.2.1              
[142] lava_1.8.1                  sessioninfo_1.2.3           backports_1.5.0            
[145] BiocParallel_1.42.0         timechange_0.3.0            gtable_0.3.6               
[148] rjson_0.2.23                umap_0.2.10.0               parallel_4.5.0             
[151] ape_5.8-1                   jsonlite_2.0.0              TFBSTools_1.46.0           
[154] bitops_1.0-9                bit64_4.6.0-1               yulab.utils_0.2.0          
[157] GOSemSim_2.34.0             UniProt.ws_2.48.0           R.utils_2.13.0             
[160] timeDate_4041.110           lazyeval_0.2.2              shiny_1.10.0               
[163] htmltools_0.5.8.1           enrichplot_1.28.2           GO.db_3.21.0               
[166] rappdirs_0.3.3              glue_1.8.0                  TFMPvalue_0.0.9            
[169] httr2_1.1.2                 RCurl_1.98-1.17             treeio_1.32.0              
[172] BSgenome_1.76.0             motifStack_1.52.0           gridExtra_2.3              
[175] igraph_2.1.4                R6_2.6.1                    cluster_2.1.8.1            
[178] wordcloud_2.6               pkgload_1.4.0               aplot_0.2.5                
[181] ipred_0.9-15                DirichletMultinomial_1.50.0 DelayedArray_0.34.1        
[184] tidyselect_1.2.1            ProtGenerics_1.40.0         ggforce_0.4.2              
[187] xml2_1.3.8                  car_3.1-3                   future_1.40.0              
[190] furrr_0.3.1                 data.table_1.17.0           htmlwidgets_1.6.4          
[193] fgsea_1.34.0                biomaRt_2.64.0              rlang_1.1.6                
[196] remotes_2.5.0               hardhat_1.4.1               prodlim_2025.04.28         
