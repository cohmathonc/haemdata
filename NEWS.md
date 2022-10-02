# Roadmap
* fix erroneous AML mouse metadata
* fix mislabelled FLT3 AML samples, add new samples
* simplify `metadata_hsa`
* include miRNA for AML samples
* revised colour scheme
* ggPCA and PC vs. time plots
* scRNAseq for CML; add AML scRNAseq h5a to metadata_mmu

# dev
* add `percent_ckit` to mmu_metadata.csv
  
# haemdata 0.0.0.9005
* add metadata_mmu template for new samples
* add pins for "all mice" `SummarisedExperiment` and csv expression matrix
* update AML sample and mouse metadata provided by Dr. Kuo: Email 2022-9-2 and 
  PSON Teams channel (`General|Copy of matadata_mmu_pivoted_AMLmice.YK.xlsx`)
* add Howto.Rmd
* add `assay` to `metadata_mmu`
* `sample_sheets` now 4-column nf-core format
* numerous small fixes

# haemdata 0.0.0.9004
* make `auth_type = "device_code"` the default for MS365R
* added CML mRNA 2022 part 2 (71 samples)
* updated vignettes with figs and examples for mRNAseq & scRNAseq
* added sample dates from PSON Teams channel (`General|AML.Seq.Samples_dates.xlsx`), along with computed columns `sample_weeks`, `age_at_end`, `age_at_start`, `age_at_sample`
* added `R-CMD-CHK` GitHub Action for CI testing
* added RNAseq reads from 90 libraries published by Kim et al 2020 [SciRep](https://www.nature.com/articles/s41598-020-76933-2); [PRJEB27973](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB27973)
* streamlined package installation by removing all but the minimal dependencies
* moved to {`Microsoft365R`} and {`pins`} for storing processed datasets and metadata. 
  
# haemdata 0.0.0.9000

* added a `NEWS.md` file to track changes to the package.
