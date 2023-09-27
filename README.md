# Temporal trends in antibiotic resistance in Europe, 1998-2019

This repository organises the code for the upcoming publication

## Folder contents

- `Masterscript.R`: overall script which calls all the individual scripts. The idea here is to make the workflow clear to an external person using the scripts.

- `/data`: folder which contains aggregated datasets (in 'original' subfolder, note that raw unaggregated data is not shared); processed versions of these data;
`/data/results`: datasets resulting from model fitting and correlation calculations.

- `/code`: folder containing code to generate temporal fits and correlations.
code/preprocessing: code that generates the aggregated datasets from the raw data (not shared).
code/obsolete: older analyses no longer relevant to paper, but kept for reference (for now)

- `/analysis`: folder containing code to further analyse data sets resulting from model fitting and correlation calculations (i.e. generate figures, tables etc).
`/analysis/output`: write figures/tables here

- `supporting_analysis`: folder containing supporting analyses - either as part of our thinking or to go in supplementary materials
supporting_analysis/output: write figures/tables here

- `package management`: Packages are managed using `renv`. All package versions are specified in the `renv.lock` and can easily be installed with this framework in your system. Futhermore, we provide the `R sessionInfo()` that was used to generate the figures.

- `author contribution`: This package shows the state of the repository at publication. The history was truncated and therefore the contributions are not reflective of the project. 

```
sessionInfo()
R version 4.2.3 (2023-03-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] rworldmap_1.3-6   sp_1.6-0          wrapr_2.0.9       scales_1.2.1      jtools_2.2.1      broom_1.0.4       ggstance_0.3.6   
 [8] DescTools_0.99.49 gridExtra_2.3     lmtest_0.9-40     zoo_1.8-12        nnet_7.3-18       svglite_2.1.1     ggplot2_3.4.2    
[15] plyr_1.8.8        tidyr_1.3.0       dplyr_1.1.2       minpack.lm_1.2-3

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9         mvtnorm_1.1-3      lattice_0.21-8     class_7.3-21       digest_0.6.31      utf8_1.2.3         R6_2.5.1          
 [8] cellranger_1.1.0   backports_1.4.1    rootSolve_1.8.2.3  e1071_1.7-12       spam_2.9-1         httr_1.4.4         pillar_1.9.0      
[15] rlang_1.1.1        Exact_3.2          readxl_1.4.1       rstudioapi_0.14    data.table_1.14.6  MetBrewer_0.2.0    car_3.1-1         
[22] Matrix_1.5-3       labeling_0.4.2     foreign_0.8-84     pander_0.6.5       munsell_0.5.0      proxy_0.4-27       compiler_4.2.3    
[29] pkgconfig_2.0.3    systemfonts_1.0.4  tidyselect_1.2.0   tibble_3.2.1       lmom_2.9           expm_0.999-7       fansi_1.0.4       
[36] viridisLite_0.4.2  crayon_1.5.2       withr_2.5.0        MASS_7.3-60        gtable_0.3.3       lifecycle_1.0.3    magrittr_2.0.3    
[43] gld_2.6.6          carData_3.0-5      cli_3.6.1          mapproj_1.2.11     farver_2.1.1       renv_0.17.3        viridis_0.6.2     
[50] ellipsis_0.3.2     generics_0.1.3     vctrs_0.6.2        boot_1.3-28.1      RColorBrewer_1.1-3 tools_4.2.3        glue_1.6.2        
[57] purrr_1.0.1        maps_3.4.1         fields_14.1        abind_1.4-5        colorspace_2.1-0   maptools_1.1-6     dotCall64_1.0-2
```
