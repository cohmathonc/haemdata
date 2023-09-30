# Roadmap
* update {targets} pipeline to use {crew} for parallelisation
* add 'normal' human BM + PBMC samples
  http://cgt.coh.org/MHO/AML.human.data/results/rnaseq.PE.GENCODEv33/MultiQC/multiqc_report.html
  http://cgt.coh.org/MHO/AML.human.data/results/rnaseq.SE.GENCODEv33/MultiQC/multiqc_report.html
* use scavenger partition on Apollo
* add state_space labels 
* add notebooks, including "cohort report"
    * use langevitour in cohort report
* remove metadata from SummarisedExperiments
* add 10X metrics for all scRNAseq libraries
* simplify & rationalise `metadata_hsa`
* implement [`pointblank`](https://rich-iannone.github.io/pointblank/index.html) data dictionaries
* revised colour scheme
* ggPCA and PC vs. time plots
* fix Isilon paths for CML, AML 2016 & 2018 mRNAseq, "devel" pin
* update nf-core version to 3.10+ and remove rRNA
* use {babelwhale} to run cutadapt
* consolidate everything into an Arrow dataframe; use tidySummarizedExperiment framework
* improve scRNAseq QC
    * https://github.com/plger/scDblFinder
    * https://github.com/broadinstitute/CellBender
    * https://github.com/wmacnair/SampleQC
    * new [hexsticker](http://gradientdescending.com/how-to-generate-a-hex-sticker-with-openai-and-cropcircles/)?
* https://github.com/cxli233/SimpleTidy_GeneCoEx
* CoFrEE: Estimate DNA Copy Number from Genome-wide RNA Expression Data
    https://www.biorxiv.org/content/10.1101/2023.08.25.554898v1
* https://github.com/omnideconv/immunedeconv