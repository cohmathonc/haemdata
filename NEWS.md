# haemdata 0.0.0.9001

* updated vignettes with figs and examples for mRNAseq & scRNAseq
* Added sample dates from PSON Teams channel (`General|AML.Seq.Samples_dates.xlsx`), along with computed columns `sample_weeks`, `age_at_end`, `age_at_start`, `age_at_sample`
* Added `R-CMD-CHK` GitHub Action for CI testing
* Added RNAseq reads from 90 libraries published by Kim et al 2020 [SciRep](https://www.nature.com/articles/s41598-020-76933-2); [PRJEB27973](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB27973)
* Streamlined package installation by removing all but the minimal dependencies
* Moved to {`Microsoft365R`} and {`pins`} for storing processed datasets and metadata. 
  
# haemdata 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
