# haemdata 0.0.0.9012
* add CML blast crisis scRNAseq experiment; 6 mice, 53 samples ("CML.blastcrisis.2023" cohort)
* add Notebook rendering to targets pipeline (via build_package() function)
* [isomiRs](https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html) object added to targets pipeline. Can publish pin if any interest
* fix FLT3 cohort metadata, include 2 rerun samples, remove irrelevant sample
* fix incorrect dod (6 CML mice)
* add functions to facilitate survival analyses: add_survival_columns(), plot_cohort_survival()
* use "code-link: true" in qmd documents

# haemdata 0.0.0.9011
* add notebooks
* renamed mirna & aml2022 10X pins
* add CML miR-142 KO scRNAseq (LMPP & T cells) 18 samples; RNAseq 25 samples ("mir142ko" cohort)
* update AML Cellplex analysis - Cellranger v7.0.0 (includes intron counts)
* change sample_id prefix to MHO from PSON
* add miRNA sample metadata from @yufu1120
* simplify Isilon paths (`/net/nfs-irwrsrchnas01/labs` & `/net/isi-dcnl/ifs/user_data` -> `/labs`)

# haemdata 0.0.0.9010
* move target pipeline to `./preprocessing`
* move _target.R script to _target.qmd
* add treatment mice metadata from @luechi
* make 2016 miRNA timepoints 0-offset to match mRNA (miRNA labels are 1-based)
* harmonise miRNA for all AML samples (nf-core/smrnaseq v2.1.0)
* miRNA fastqs gzip compressed in place, keeping originals
* fix mislabelled human FLT3 AML samples, add new samples
* retired GRCm38 mRNAseq pipeline; GENCODEm28 is the default mmu reference genome for mRNAseq
* add a MHO sample identifier for all mmu samples
* add a warning about reproducibility when pins are loaded without setting `version`
* rename `sample` to `library_id` in metadata_mmu
* rename `project` to `cohort` in metadata_mmu
* add `ref_dim1` & `ref_dim2` to metadata_mmu, for PC, UMAP or other coordinates from dimensionality reduction
* nf-core sample sheets use `library_id` column for `sample` column

# haemdata 0.0.0.9009
  * fix metadata_mmu 
  * fix pin versions on Onedrive

# haemdata 0.0.0.9008
  * freeze on updating GRCm38 pins
  * include AML scRNAseq 10X hda5 file paths in metadata_mmu
  * implement update_metadata_mmu() and retire make_metadata_mmu()
  * fix erroneous AML mouse sample metadata identified by @yufu1120 (DOD)
  * fix erroneous AML scRNAseq bulk sample labels (PBMC_CKIT)
  * fix mislabelled scRNAseq 10X samples (previous fix 321c166; now fixed in Cellranger sample sheet)
  * other small typos and spelling corrections

# haemdata 0.0.0.9007
  * add SingleR `cell_type` and `cell_type_fine` to Seurat objects/h5ad files
  
# haemdata 0.0.0.9006
* add ROADMAP.md
* add `percent_ckit`, `qc_pass_mapping` to metadata_mmu.csv
* flag rather than drop mouse samples that fail mapping (`qc_pass_mapping`)
* fix erroneous AML mouse sample metadata identified by Yu-Husan Fu (@yufu1120)
* other small fixes
  
# haemdata 0.0.0.9005
* add metadata_mmu template for new samples
* add pins for "all mice" `SummarisedExperiment` and csv expression matrix
* update AML sample and mouse metadata provided by @yahueikuo: Email 2022-9-2 and 
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
