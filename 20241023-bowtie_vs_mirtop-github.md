20241023-miRNA-bowtie vs mirtop
================
yufu
2024-10-23

## Aim

This document aims to compare the difference/similarity of two type of
miRNA expression estimation method (bowtie and mirtop).  
bowtie: Its basically the ‘quick & dirty’ way to count, then the
pipeline turns this into CPM for visualisation. In contrast, mirtop
collapses counts across all iso-miRs for a gene, and reports counts per
miR.  
These were previously done in
<https://cgt.coh.org/haemdata/news/index.html#haemdata-0009010>
(nf-core/smrnaseq v2.1.0)

The comparison will be made at not only gene expression correlation but
also the state-space.  

## Load Required Libraries

``` r
# load the pinboard before render the file
setwd(r'(C:\Users\yufu\OneDrive - City of Hope National Medical Center\Documents\miRNA_AML)')
# library(haemdata)
# library(glmGamPoi)
# use_pinboard("onedrive") #when using local R
# use_pinboard("devel") #when using R on the apollo server


# loading package
library(haemdata)
library(glmGamPoi)
library(tidyr)
library(edgeR)
library(devtools)
library(ggbiplot)
library(dplyr)
library(compositions)
library(tidyHeatmap)
library(factoextra)
library(openxlsx)
```

## Data Loading and Processing

Load data from haemdata  
bowtie_cpm have been deposited in the haemdata, while mirtop is still
raw count and need to be converted here.

``` r
# loading data from haemdata
get_pin_list()
```

    ##  [1] "hsa_mrna_flt3_GENCODEr40_qc.rds"                          "hsa_mrna_flt3_GENCODEr40_qc_1tpm_in_5samples.csv"        
    ##  [3] "hsa_mrna_kim_GENCODEr40_qc.rds"                           "hsa_mrna_kim_GENCODEr40_qc_1tpm_in_5samples.csv"         
    ##  [5] "hsa_mrna_mds_GENCODEr40_qc.rds"                           "hsa_mrna_mds_GENCODEr40_qc_1tpm_in_5samples.csv"         
    ##  [7] "hsa_mrna_non_aml_GENCODEr40_qc.rds"                       "hsa_mrna_non_aml_GENCODEr40_qc_1tpm_in_5samples.csv"     
    ##  [9] "metadata_hsa.csv"                                         "metadata_mmu.csv"                                        
    ## [11] "mmu_10x_aml2022_GENCODEm28_HLT.h5ad"                      "mmu_10x_aml2022_GENCODEm28_HLT.rds"                      
    ## [13] "mmu_10x_blastcrisis_GENCODEm28_HLT.h5ad"                  "mmu_10x_blastcrisis_GENCODEm28_HLT.rds"                  
    ## [15] "mmu_10x_mir142ko_GENCODEm28_HLT.h5ad"                     "mmu_10x_mir142ko_GENCODEm28_HLT.rds"                     
    ## [17] "mmu_mirna_all_mice_miRBase22_bowtie_cpm.csv"              "mmu_mirna_all_mice_miRBase22_mirna_qc.csv"               
    ## [19] "mmu_mirna_all_mice_miRBase22_mirtop_counts.csv"           "mmu_mrna_all_mice_GENCODEm28_HLT_qc.rds"                 
    ## [21] "mmu_mrna_all_mice_GENCODEm28_HLT_qc_1tpm_in_5samples.csv" "mmu_mrna_techrep_GENCODEm28_HLT_qc.rds"                  
    ## [23] "mmu_mrna_techrep_GENCODEm28_HLT_qc_1tpm_in_5samples.csv"

``` r
metadata<-get_pin("metadata_mmu.csv", version = "20230907T012227Z-d3324")


cpm.bowtie<-get_pin("mmu_mirna_all_mice_miRBase22_bowtie_cpm.csv", version = "20230907T012024Z-a1f46")
colnames(cpm.bowtie)[1]<-"miRNA"
cpm.bowtie <- cpm.bowtie %>%
  dplyr::select(-miRNA) %>%
  magrittr::set_rownames(cpm.bowtie$miRNA)

counts.mirtop<-get_pin("mmu_mirna_all_mice_miRBase22_mirtop_counts.csv", version = "20230907T012005Z-6a6f3")
counts.mirtop <- counts.mirtop %>%
  dplyr::select(-miRNA) %>%
  magrittr::set_rownames(counts.mirtop$miRNA)


# replace NA to 0
sum(is.na(cpm.bowtie))
```

    ## [1] 3688

``` r
cpm.bowtie[is.na(cpm.bowtie)]<-0

sum(is.na(counts.mirtop))
```

    ## [1] 2936

``` r
counts.mirtop[is.na(counts.mirtop)]<-0
```

# Convert mirtop raw counts to cpm

Codes from
<https://github.com/drejom/smrnaseq/blob/master/bin/edgeR_miRBase.r>
(also cloned from nf-core/smrnaseq v2.1.0)

``` r
# make mirtop into cpm ================================================
raw_to_cpm<-function(data){
  # Remove genes with 0 reads in all samples
  row_sub = apply(data, 1, function(row) all(row ==0 ))
  # Only subset if at least one sample is remaining
  nr_keep <- sum(row_sub)
  if (nr_keep > 0){
    data<-data[!row_sub,]
  }
  # Also check for colSums > 0, otherwise DGEList will fail if samples have entirely colSum == 0 #Fixes #134
  drop_colsum_zero <- (colSums(data, na.rm=T) != 0) # T if colSum is not 0, F otherwise
  data <- data[, drop_colsum_zero] # all the non-zero columns

  # Normalization
  dataDGE<-DGEList(counts=data,genes=rownames(data))
  o <- order(rowSums(dataDGE$counts), decreasing=TRUE)
  dataDGE <- dataDGE[o,]

  # TMM
  dataNorm <- calcNormFactors(dataDGE)

  # Print normalized read counts to file
  dataNorm_df<-as.data.frame(cpm(dataNorm))
  return(dataNorm_df)

  }

cpm.mirtop<-raw_to_cpm(counts.mirtop)
```

# Extract samples from 2018 cohort

``` r
# arrange the time point and mouse ID in metadata

metadata$timepoint <- factor(metadata$timepoint, levels = c("T0","T1","T2","T3","T4","T5","T6","T6p5","T7","T8","T9","T10","W0",
                                                                     "W2","W3","W4","W5","W6","W7","W8","W9","W10","W11","W12","W13","W14",
                                                                     "W15","W16","W17","W18","W19","W20","W22","W26","W28","W34","W35","W36",
                                                                     "W40","W41","W48","END"))

metadata$mouse_id <- factor(metadata$mouse_id, levels = c("3334","3335","3336","3338","3339","3340","3341","3342","3346","3349",
                                                        "3350","3357","3368","3370","4309","4311","4321","4324","4329","4419",
                                                        "4433","4436","4443","4501","4506","4510","4534","4535"))

metadata<- metadata %>%
  arrange(mouse_id, timepoint)
# create 2018 cohort list

# miRNA129_ID<-metadata.local$ID[metadata.local$cohort=="AML.miRNA.2018"]
# saveRDS(miRNA129_ID,"miRNA129_ID.rds")
miRNA129_ID<-readRDS("miRNA129_ID.rds")

metadata.AML.miRNA.2018<-metadata[metadata$library_id %in% miRNA129_ID,]%>%
  drop_na(timepoint)

AML.miRNA.2018.library_id<-metadata.AML.miRNA.2018$library_id

# # intersect miRNA species
# miRNA.intersect<-intersect(row.names(counts.bowtie), row.names(counts.mirtop))
# counts.bowtie<-counts.bowtie[row.names(counts.bowtie) %in% miRNA.intersect,]
# counts.mirtop<-counts.mirtop[row.names(counts.mirtop) %in% miRNA.intersect,]

# extract 2018 from bowtie and re-arrange the order of samples in the matrix to fit the order in metadata
cpm.bowtie.2018<-cpm.bowtie[,which(colnames(cpm.bowtie) %in% AML.miRNA.2018.library_id)]
index<-match(metadata.AML.miRNA.2018$library_id, colnames(cpm.bowtie.2018))
cpm.bowtie.2018<-cpm.bowtie.2018[,index]

# extract 2018 from mirtop and re-arrange the order of samples in the matrix to fit the order in metadata
cpm.mirtop.2018<-cpm.mirtop[,which(colnames(cpm.mirtop) %in% AML.miRNA.2018.library_id)]
index<-match(metadata.AML.miRNA.2018$library_id, colnames(cpm.mirtop.2018))
cpm.mirtop.2018<-cpm.mirtop.2018[,index]

# check colnames are at the same order in bowtie and mirtop
all(colnames(cpm.mirtop.2018) == colnames(cpm.bowtie.2018))
```

    ## [1] TRUE

# Make correlation heatmaps between bowtie and mirtop

``` r
# make the combined_cpm matrix 
cpm.bowtie.2018.2<-cpm.bowtie.2018
cpm.mirtop.2018.2<-cpm.mirtop.2018

cpm.bowtie.2018.2$gene<-row.names(cpm.bowtie.2018.2)
cpm.mirtop.2018.2$gene<-row.names(cpm.mirtop.2018.2)

combined_cpm <- dplyr::left_join(cpm.bowtie.2018.2, cpm.mirtop.2018.2, by="gene", suffix = c(".bowtie", ".mirtop")) |>
  tidyr::drop_na()

mat <- combined_cpm |>
  dplyr::select(-gene) |>
  as.matrix() |>
  magrittr::set_rownames(combined_cpm$gene)


# draw correlation heatmap with dendrogram:
correlate_combined_runs <- function(mat, sample_sheet) {
  mat <- mat
  combined_sample_sheet <- sample_sheet
  
  corr_mat_long <- (cor(compositions::acomp(mat), method = "spearman"))^2 |>
    as_tibble(rownames = "sample1") |>
    pivot_longer(!sample1, names_to = "sample2", values_to = "corr_R^2")

  corr_mat_long<-corr_mat_long%>%
    separate(sample1, into = c(NA,"library_id"), remove = F)
  corr_mat_long$library_id <- paste0("COHP_",corr_mat_long$library_id)


  tidy_combined_data <- dplyr::left_join(
    corr_mat_long,
    combined_sample_sheet,
    by = c("library_id" = "library_id")
  )
  unique(tidy_combined_data$timepoint)
  p <- tidy_combined_data |>
    tidyHeatmap::heatmap(sample1,
                         sample2,
                         `corr_R^2`,
                         scale = "none",
                         palette_value = circlize::colorRamp2(c(-0.01, 0, 1), c("red", "white", "blue"))) |>
    tidyHeatmap::add_tile(mouse_id) |>
    tidyHeatmap::add_tile(timepoint) |>
    tidyHeatmap::add_tile(treatment) |>
    tidyHeatmap::add_tile(sex) |>
    tidyHeatmap::add_tile(sample_weeks)

  # tidyHeatmap::save_pdf(p, paste0("results/", matrix, "_corr_heatmap.pdf"), width = 16L, units = "in")
  
  return(p)
}

# correlation heatmap without dendrogram:
correlate_combined_runs_2 <- function(mat, sample_sheet) {
  mat <- mat
  combined_sample_sheet <- sample_sheet
  
  corr_mat_long <- (cor(compositions::acomp(mat), method = "spearman"))^2 |>
    as_tibble(rownames = "sample1") |>
    pivot_longer(!sample1, names_to = "sample2", values_to = "corr_R^2")
  
  corr_mat_long<-corr_mat_long%>%
    separate(sample1, into = c(NA,"library_id"), remove = F)
  corr_mat_long$library_id <- paste0("COHP_",corr_mat_long$library_id)
  
  
  tidy_combined_data <- dplyr::left_join(
    corr_mat_long,
    combined_sample_sheet,
    by = c("library_id" = "library_id")
  )
  
  p <- tidy_combined_data |>
    tidyHeatmap::heatmap(sample1,
                         sample2,
                         `corr_R^2`,
                         scale = "none",
                         cluster_rows = F,
                         cluster_columns = FALSE,
                         palette_value = circlize::colorRamp2(c(-0.01, 0, 1), c("red", "white", "blue"))) |>
    tidyHeatmap::add_tile(mouse_id) |>
    tidyHeatmap::add_tile(timepoint) |>
    tidyHeatmap::add_tile(treatment) |>
    tidyHeatmap::add_tile(sex) |>
    tidyHeatmap::add_tile(sample_weeks)
  
  # tidyHeatmap::save_pdf(p, paste0("results/", matrix, "_corr_heatmap.pdf"), width = 16L, units = "in")
  
  return(p)
 } 

correlate_combined_runs(mat,metadata.AML.miRNA.2018)
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 66564 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/make_correlation-1.png)<!-- -->

``` r
correlate_combined_runs_2(mat,metadata.AML.miRNA.2018)
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 66564 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/make_correlation-2.png)<!-- -->

## PCA

cpm were subjected to offset, log transform, then perform PCA by
{prcomp}

``` r
var.offset.t.log<-function(mat){
  # drop those rows with zero var
  mat<-mat[which(apply(mat, 1, var) != 0), ]
  # add offset to the dataset
  mat<-mat+0.01
  # transpose
  mat<-t(mat)
  # make log
  mat<-log(mat,2)
  return(mat)
}


cpm.bowtie.2018.pcocessed<-var.offset.t.log(cpm.bowtie.2018)
cpm.mirtop.2018.pcocessed<-var.offset.t.log(cpm.mirtop.2018)

# for test
a<-cpm.bowtie.2018.pcocessed
meta<-metadata.AML.miRNA.2018
suffix<-"bowtie"

# PCA
PCA_plot<-function(a,meta,suffix, negative){
# perfrom PCA
pc <- prcomp(a,
             center = TRUE,
             scale. = F)
# elbow plot 
elb<-fviz_eig(pc)+
  ggtitle(suffix)

# adjust the orientation
if (negative==T){
  pc$x=-pc$x
  pc$rotation=-pc$rotation
  
}

# PC1 PC2 plot
PC12<-ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = meta[["treatment"]],
              # colour = meta[["treatment"]],
              ellipse = F,
              circle = F,
              ellipse.prob = 0.68,
              var.axes = FALSE,)+
  geom_path(aes(colour = meta[["treatment"]]), arrow = arrow(), group= meta[["mouse_id"]])+
  ggtitle(suffix)


pc.sample.df<-as.data.frame(pc$x) #sample
pc.sample.df$ID<-row.names(pc.sample.df)

pc.mi.df<-as.data.frame(pc$rotation)  #miRNA

# miRNA-PC1 vs PC2 colored by sample
PC12.sample<-ggbiplot(pc,
               obs.scale = 1,
               var.scale = 1,
               groups = meta[["mouse_id"]],
               # colour = meta[["treatment"]],
               ellipse = F,
               circle = F,
               ellipse.prob = 0.68,
               var.axes = FALSE)+
  geom_path(aes(colour = meta[["mouse_id"]]), arrow = arrow(), group= meta[["mouse_id"]])+
  ggtitle(suffix)


print(elb)
print(PC12)
# print(PC12.sample)


aa<-list("pc" = pc,
     "elb" = elb,
     "PC12" = PC12,
     "pc.sample.df"=pc.sample.df,
     "pc.mi.df"=pc.mi.df)

return(aa)

}

set.seed(126)
PCA.bowtie<-PCA_plot(cpm.bowtie.2018.pcocessed,metadata.AML.miRNA.2018,"bowtie",negative = F)
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PCA-1.png)<!-- -->![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PCA-2.png)<!-- -->

``` r
PCA.mirtop<-PCA_plot(cpm.mirtop.2018.pcocessed,metadata.AML.miRNA.2018,"mirtop",negative = T)
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PCA-3.png)<!-- -->![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PCA-4.png)<!-- -->

## Merge PC values from bowtie and mirtop, make plots to demonstrate the PC1 and PC2

``` r
# merge
total.meta<-merge(metadata.AML.miRNA.2018, PCA.bowtie[["pc.sample.df"]],c(1:2,ncol(PCA.bowtie)), by.x="library_id", by.y = "ID")
total.meta<-merge(total.meta, PCA.mirtop[["pc.sample.df"]],c(1:2,ncol(PCA.mirtop)), by.x="library_id", by.y = "ID", suffixes = c(".bowtie",".mirtop"))
total.meta.sort <- total.meta %>%                 # Order table with dplyr
  arrange(timepoint)

# write table
cohort<-"miRNA-2018-bowtie_vs_mirtop"
file.name<-paste(Sys.Date(),cohort,".xlsx", sep = "-")
print(file.name)
```

    ## [1] "2024-10-23-miRNA-2018-bowtie_vs_mirtop-.xlsx"

``` r
openxlsx::write.xlsx(total.meta.sort,file =file.name)

# miRNA-PC1 vs PC2
ggplot(total.meta.sort, aes(x=PC1.bowtie, y=PC2.bowtie, group=mouse_id, color=mouse_id))+
  geom_point()+
  geom_path(arrow = arrow())+
  ggtitle("PC1 vs PC2-bowtie")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/other_PC_plot-1.png)<!-- -->

``` r
ggplot(total.meta.sort, aes(x=PC1.mirtop, y=PC2.mirtop, group=mouse_id, color=mouse_id))+
  geom_point()+
  geom_path(arrow = arrow())+
  ggtitle("PC1 vs PC2-mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/other_PC_plot-2.png)<!-- -->

``` r
# miRNA-PC1 vs time
ggplot(total.meta.sort, aes(x=timepoint, y=PC1.bowtie, group=mouse_id, color=mouse_id))+
  geom_point(aes(shape=treatment), size = 2)+
  scale_shape_manual(values=c(16, 1))+
  geom_path(aes(linetype=treatment))+
  scale_y_reverse()+
  theme(axis.text.x = element_text(size=8,angle=45, hjust = 1, vjust = 1))+
  ggtitle("PC1 vs time-bowtie")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/other_PC_plot-3.png)<!-- -->

``` r
ggplot(total.meta.sort, aes(x=timepoint, y=PC1.mirtop, group=mouse_id, color=mouse_id))+
  geom_point(aes(shape=treatment), size = 2)+
  scale_shape_manual(values=c(16, 1))+
  geom_path(aes(linetype=treatment))+
  scale_y_reverse()+
  theme(axis.text.x = element_text(size=8,angle=45, hjust = 1, vjust = 1))+
  ggtitle("PC1 vs time-mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/other_PC_plot-4.png)<!-- -->

``` r
# miRNA-PC2 vs time
ggplot(total.meta.sort, aes(x=timepoint, y=PC2.bowtie, group=mouse_id, color=mouse_id))+
  geom_point(aes(shape=treatment), size = 2)+
  scale_shape_manual(values=c(16, 1))+
  geom_path(aes(linetype=treatment))+
  scale_y_reverse()+
  theme(axis.text.x = element_text(size=8,angle=45, hjust = 1, vjust = 1))+
  ggtitle("PC2 vs time-bowtie")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/other_PC_plot-5.png)<!-- -->

``` r
ggplot(total.meta.sort, aes(x=timepoint, y=PC2.mirtop, group=mouse_id, color=mouse_id))+
  geom_point(aes(shape=treatment), size = 2)+
  scale_shape_manual(values=c(16, 1))+
  geom_path(aes(linetype=treatment))+
  scale_y_reverse()+
  theme(axis.text.x = element_text(size=8,angle=45, hjust = 1, vjust = 1))+
  ggtitle("PC2 vs time-mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/other_PC_plot-6.png)<!-- -->

## Make plots to compare bowtie and mirtop

``` r
ggplot(total.meta.sort, aes(x=PC1.bowtie, y=PC1.mirtop, group=mouse_id, color=mouse_id))+
  geom_point()+
  geom_path(arrow = arrow())+
  ggtitle("PC1.bowtie vs PC1.mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PC_plot_between_bowtiw_mirtop-1.png)<!-- -->

``` r
ggplot(total.meta.sort, aes(x=PC2.bowtie, y=PC2.mirtop, group=mouse_id, color=mouse_id))+
  geom_point()+
  geom_path(arrow = arrow())+
  ggtitle("PC2.bowtie vs PC2.mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PC_plot_between_bowtiw_mirtop-2.png)<!-- -->

``` r
# merge the miRNA eigenvalue from bowtie and mirtop
bowtie.eigenvalue<-PCA.bowtie[["pc.mi.df"]][,c(1:2)]
bowtie.eigenvalue$ID<-row.names(PCA.bowtie[["pc.mi.df"]])
mirtop.eigenvalue<-PCA.mirtop[["pc.mi.df"]][,c(1:2)]
mirtop.eigenvalue$ID<-row.names(PCA.mirtop[["pc.mi.df"]])

merge.eigenvalue<-merge(bowtie.eigenvalue, mirtop.eigenvalue, by="ID", suffixes = c(".bowtie",".mirtop"))
all(is.na(merge.eigenvalue))
```

    ## [1] FALSE

``` r
ggplot(merge.eigenvalue, aes(x=PC1.bowtie, y=PC1.mirtop))+
  geom_point()+
  ggtitle("eigenvalue-PC1.bowtie vs PC1.mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PC_plot_between_bowtiw_mirtop-3.png)<!-- -->

``` r
ggplot(merge.eigenvalue, aes(x=PC2.bowtie, y=PC2.mirtop))+
  geom_point()+
  ggtitle("eigenvalue-PC2.bowtie vs PC2.mirtop")
```

![](20241023-bowtie_vs_mirtop-github_files/figure-gfm/PC_plot_between_bowtiw_mirtop-4.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] openxlsx_4.2.5.2    factoextra_1.0.7    tidyHeatmap_1.8.1   compositions_2.0-8  dplyr_1.1.4         ggbiplot_0.55       scales_1.3.0       
    ##  [8] plyr_1.8.9          ggplot2_3.5.1       devtools_2.4.5      usethis_2.2.3       edgeR_3.36.0        limma_3.50.3        tidyr_1.3.0        
    ## [15] glmGamPoi_1.6.0     haemdata_0.0.0.9012
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_2.1-0            ggsignif_0.6.4              rjson_0.2.21                ellipsis_0.3.2              circlize_0.4.16            
    ##   [6] XVector_0.34.0              GenomicRanges_1.46.1        GlobalOptions_0.1.2         fs_1.6.4                    clue_0.3-65                
    ##  [11] rstudioapi_0.16.0           farver_2.1.2                ggpubr_0.6.0                remotes_2.5.0               ggrepel_0.9.4              
    ##  [16] fansi_1.0.5                 codetools_0.2-20            doParallel_1.0.17           cachem_1.0.8                robustbase_0.99-3          
    ##  [21] knitr_1.48                  pkgload_1.4.0               jsonlite_1.8.8              Cairo_1.6-2                 broom_1.0.6                
    ##  [26] cluster_2.1.6               Microsoft365R_2.4.0         png_0.1-8                   shiny_1.8.1.1               compiler_4.1.1             
    ##  [31] httr_1.4.7                  backports_1.5.0             Matrix_1.6-4                fastmap_1.1.1               cli_3.6.1                  
    ##  [36] later_1.3.1                 htmltools_0.5.7             tools_4.1.1                 gtable_0.3.5                glue_1.6.2                 
    ##  [41] GenomeInfoDbData_1.2.7      rappdirs_0.3.3              Rcpp_1.0.11                 carData_3.0-5               Biobase_2.54.0             
    ##  [46] vctrs_0.6.5                 iterators_1.0.14            tensorA_0.36.2.1            xfun_0.45                   stringr_1.5.1              
    ##  [51] mime_0.12                   miniUI_0.1.1.1              lifecycle_1.0.4             rstatix_0.7.2               dendextend_1.17.1          
    ##  [56] DEoptimR_1.1-3              zlibbioc_1.40.0             MASS_7.3-60                 promises_1.2.1              MatrixGenerics_1.6.0       
    ##  [61] parallel_4.1.1              SummarizedExperiment_1.24.0 RColorBrewer_1.1-3          ComplexHeatmap_2.20.0       yaml_2.3.9                 
    ##  [66] curl_5.2.1                  memoise_2.0.1               gridExtra_2.3               jose_1.2.0                  stringi_1.8.2              
    ##  [71] highr_0.11                  pins_1.3.0                  S4Vectors_0.32.4            foreach_1.5.2               BiocGenerics_0.40.0        
    ##  [76] zip_2.3.0                   AzureAuth_1.3.3             pkgbuild_1.4.4              AzureGraph_1.3.4            shape_1.4.6.1              
    ##  [81] GenomeInfoDb_1.30.1         rlang_1.1.2                 pkgconfig_2.0.3             matrixStats_1.1.0           bitops_1.0-7               
    ##  [86] evaluate_0.24.0             lattice_0.22-5              purrr_1.0.2                 labeling_0.4.3              patchwork_1.2.0            
    ##  [91] htmlwidgets_1.6.4           tidyselect_1.2.1            magrittr_2.0.3              R6_2.5.1                    IRanges_2.28.0             
    ##  [96] generics_0.1.3              profvis_0.3.8               DelayedArray_0.20.0         withr_3.0.0                 pillar_1.9.0               
    ## [101] abind_1.4-5                 RCurl_1.98-1.9              tibble_3.2.1                bayesm_3.1-6                car_3.1-2                  
    ## [106] crayon_1.5.3                utf8_1.2.4                  urlchecker_1.0.1            rmarkdown_2.27              viridis_0.6.5              
    ## [111] GetoptLong_1.0.5            locfit_1.5-9.10             digest_0.6.33               xtable_1.8-4                httpuv_1.6.12              
    ## [116] openssl_2.2.0               stats4_4.1.1                munsell_0.5.1               viridisLite_0.4.2           sessioninfo_1.2.2          
    ## [121] askpass_1.2.0
