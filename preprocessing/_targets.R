# Setup
library(targets)
library(tarchetypes)
library(hprcc)
#tar_unscript()

# Location of the nf-core cache and other intermediate analysis files
nf_core_cache <- "/labs/rrockne/MHO/haemdata-nf-core-cache"

# Load the R scripts & functions:
tar_source(here::here("preprocessing/scripts"))

tar_option_set(
    packages = c("tidyverse", "SummarizedExperiment", "haemdata", "tarchetypes", "hprcc"),
    error = "abridge" # continue with running targets on error
)

# Pipeline globals
# nf-core pipeline versions
rnaseq_release <- "3.7"
smrnaseq_release <- "2.1.0"
scrnaseq_release <- "2.1.0"

# Publish pins to onedrive or devel?
publish_location <- "devel"

# Cohort aliases
# Regex patterns to select cohorts using dplyr::filter() or similar functions.
all_mice <- "^AML.mRNA.2016$|^AML.mRNA.2018.all_samples$|^AML.mRNA.2020$|
        |^AML.mRNA.2021.RxGroup1$|^AML.mRNA.2021.RxGroup2$|^AML.mRNA.2021.RxGroup2_pt2$|
        |^AML.mRNA.2022.RxGroup3$|^CML.mRNA.2021$|^CML.mRNA.2022$|^CML.mRNA.2022_pt2$"

aml_mice <- "^AML.mRNA.2016$|^AML.mRNA.2018.all_samples$|^AML.mRNA.2020$|
        |^AML.mRNA.2021.RxGroup1$|^AML.mRNA.2021.RxGroup2$|^AML.mRNA.2021.RxGroup2_pt2$|
        |^AML.mRNA.2022.RxGroup3$"

chemo_mice <- "^AML.mRNA.2020$"

treatment_mice <- "^AML.mRNA.2021.RxGroup1$|^AML.mRNA.2021.RxGroup2$|
        |^AML.mRNA.2021.RxGroup2_pt2$|^AML.mRNA.2022.RxGroup3$"

cml_mice <- "^CML.mRNA.2021$|^CML.mRNA.2022$|^CML.mRNA.2022_pt2$"

old_mice <- "^senescence$"
list(
# Metadata
# Mouse metadata
tar_target(metadata_mmu_prepub, update_metadata_mmu()),

# Human metadata
tar_target(
    metadata_hsa,
    make_metadata_hsa(
        dplyr::full_join(sample_sheet_2022_3, sample_sheet_2022_2) |>
        dplyr::full_join(sample_sheet_2022_4))
    ),

# Assays
# mRNA
# Sample sheets
tar_target(sample_sheet_2016_2022, make_rnaseq_sample_sheet(all_mice, metadata_mmu_prepub)),
tar_target(sample_sheet_CML_3, make_rnaseq_sample_sheet("^CML.mRNA.2022_pt2$", metadata_mmu_prepub)),
tar_target(sample_sheet_mrna_mir142ko, make_rnaseq_sample_sheet("^CML.mRNA.mir142ko$", metadata_mmu_prepub)),
tar_target(sample_sheet_mrna_old, make_rnaseq_sample_sheet("^senescence.mRNA.2022$", metadata_mmu_prepub)),

# Symlink fastqs
tar_target(sample_sheet_mir142ko_mrna_cache, symlink_mrna_fastqs("mmu_mrna_mir142ko/fastq", sample_sheet_mrna_mir142ko)),
tar_target(sample_sheet_old_mrna_cache, symlink_mrna_fastqs("mmu_mrna_old_mice/fastq", sample_sheet_mrna_old)),

# AML validation samples and human datasets
tar_target(sample_sheet_2017_1, parse_metadata_AML.validation.2017()),
tar_target(sample_sheet_2020_2, parse_metadata_AML.mRNA.novaseq_validation.2020()),
tar_target(
    sample_sheet_techrep,
    rbind(
        sample_sheet_2017_1,
        sample_sheet_2020_2
    )
),
tar_target(sample_sheet_2022_2, parse_metadata_AML.mRNA.HSA_FLT3.2022()),
tar_target(sample_sheet_2022_3, parse_metadata_MDS.rnaseq.EGAD00001003891()),
tar_target(sample_sheet_2022_4, parse_metadata_AML.PRJEB27973()),
tar_target(sample_sheet_non_aml, dplyr::filter(metadata_hsa, grepl("AML.non_aml", cohort))),

# Run nf-core/rnaseq
tar_target(mmu_mrna_2016_2022_qc,
    run_nf_core_rnaseq("mmu_mrna_2016_2022", sample_sheet_2016_2022, "GENCODEm28_HLT"),
    format = "file"
),
tar_target(mmu_mrna_techrep_qc,
    run_nf_core_rnaseq("mmu_mrna_techrep", sample_sheet_techrep, "GENCODEm28_HLT"),
    format = "file"
),
tar_target(mmu_mrna_cml3_qc,
    run_nf_core_rnaseq("mmu_mrna_cml3_qc", sample_sheet_CML_3, "GENCODEm28_HLT"),
    format = "file"
),
tar_target(mmu_mrna_mir142ko_qc,
    run_nf_core_rnaseq("mmu_mrna_mir142ko", sample_sheet_mir142ko_mrna_cache, "GENCODEm28_HLT"),
    format = "file"
),
tar_target(mmu_mrna_old_qc,
    run_nf_core_rnaseq("mmu_mrna_old_mice", sample_sheet_old_mrna_cache, "GENCODEm28_HLT"),
    format = "file"
),
tar_target(hsa_mrna_flt3_qc,
    run_nf_core_rnaseq("hsa_mrna_flt3", sample_sheet_2022_2, "GENCODEr40"),
    format = "file"
),
tar_target(hsa_mrna_mds_qc,
    run_nf_core_rnaseq("hsa_mrna_mds", sample_sheet_2022_3, "GENCODEr40"),
    format = "file"
),
tar_target(hsa_mrna_kim_qc,
    run_nf_core_rnaseq("hsa_mrna_kim", sample_sheet_2022_4, "GENCODEr40"),
    format = "file"
),
tar_target(hsa_mrna_non_aml_qc,
    run_nf_core_rnaseq("hsa_mrna_non_aml", sample_sheet_non_aml, "GENCODEr40"),
    format = "file"
),

# Make SummarisedExperiments
# Mouse datasets
tar_target(mmu_mrna_2016_2022_qc_se,
    annotate_se(
        get_rnaseq_se(mmu_mrna_2016_2022_qc),
        metadata_mmu_prepub,
        mmu_mrna_2016_2022_qc
)),
tar_target(
    mmu_mrna_techrep_qc_se,
    annotate_se(
        get_rnaseq_se(mmu_mrna_techrep_qc),
        metadata_mmu_prepub,
        mmu_mrna_techrep_qc
)),
tar_target(
    mmu_mrna_cml3_qc_se,
    annotate_se(
        get_rnaseq_se(mmu_mrna_cml3_qc),
        metadata_mmu_prepub,
        mmu_mrna_cml3_qc
)),
tar_target(
    mmu_mrna_mir142ko_qc_se,
    annotate_se(
        get_rnaseq_se(mmu_mrna_mir142ko_qc),
        metadata_mmu_prepub,
        mmu_mrna_mir142ko_qc
)),
tar_target(
    mmu_mrna_old_qc_se,
    annotate_se(
        get_rnaseq_se(mmu_mrna_old_qc),
        metadata_mmu_prepub,
        mmu_mrna_old_qc
)),

# Human datasets
tar_target(hsa_mrna_flt3_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_flt3_qc), metadata_hsa, hsa_mrna_flt3_qc)),
tar_target(hsa_mrna_mds_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_mds_qc), metadata_hsa, hsa_mrna_mds_qc)),
tar_target(hsa_mrna_kim_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_kim_qc), metadata_hsa, hsa_mrna_kim_qc)),
tar_target(hsa_mrna_non_aml_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_non_aml_qc), metadata_hsa, hsa_mrna_non_aml_qc)),

# Merge SummarisedExperiments
tar_target(
    mmu_mrna_2016_2022_cml3_mice_qc_se, 
    merge_mrna_se(
        mmu_mrna_2016_2022_qc_se,
        mmu_mrna_cml3_qc_se,
        new_name = "mmu_mrna_all_mice_GENCODEm28_HLT_qc")
),
tar_target(
    mmu_mrna_2016_2022_cml3_mir142ko_mice_qc_se, 
    merge_mrna_se(
        mmu_mrna_2016_2022_cml3_mice_qc_se,
        mmu_mrna_mir142ko_qc_se,
        new_name = "mmu_mrna_all_mice_GENCODEm28_HLT_qc")
),
tar_target(
    mmu_mrna_all_mice_qc_se, 
    merge_mrna_se(
        mmu_mrna_2016_2022_cml3_mir142ko_mice_qc_se,
        mmu_mrna_old_qc_se,
        new_name = "mmu_mrna_all_mice_GENCODEm28_HLT_qc")
),

# Quality control
tar_target(
    metadata_mmu,
    flag_lowly_mapped_mmu_se(mmu_mrna_all_mice_qc_se, metadata_mmu_prepub)),

# Remove human samples that fail the STAR uniquely mapped reads threshold of 5%
tar_target(hsa_mrna_flt3_qc_se_flt, qc_filter_se(hsa_mrna_flt3_qc_se)),
tar_target(hsa_mrna_mds_qc_se_flt, qc_filter_se(hsa_mrna_mds_qc_se)),
tar_target(hsa_mrna_kim_qc_se_flt, qc_filter_se(hsa_mrna_kim_qc_se)),
tar_target(hsa_mrna_non_aml_qc_se_flt, qc_filter_se(hsa_mrna_non_aml_qc_se)),

# miRNA
# Sample sheets
tar_target(miRNA_sample_sheet_pretrimmed, pretrimmed_sample_sheet(metadata_mmu_prepub)),
tar_target(
    miRNA_sample_sheet_untrimmed, untrimmed_sample_sheet(metadata_mmu_prepub)),
tar_target(miRNA_sample_sheet_trimmed, run_cutadapt_smrna(miRNA_sample_sheet_untrimmed)),

# nf-core/smRNAseq pipeline
tar_target(mmu_mirna_pretrimmed_multiqc,
    run_nfcore_smrnaseq("mmu_mirna_pretrimmed", miRNA_sample_sheet_pretrimmed),
    format = "file"
    ),
tar_target(mmu_mirna_trimmed_multiqc,
    run_nfcore_smrnaseq("mmu_mirna_untrimmed", miRNA_sample_sheet_trimmed),
    format = "file"
),

# Bowtie counts
tar_target(mmu_mirna_pretrimmed_cpm, nfcore_smrna_cpm(mmu_mirna_pretrimmed_multiqc)),
tar_target(mmu_mirna_trimmed_cpm, nfcore_smrna_cpm(mmu_mirna_trimmed_multiqc)),
tar_target(mmu_mirna_all_mice_miRBase22_bowtie_cpm, merge(mmu_mirna_pretrimmed_cpm, mmu_mirna_trimmed_cpm, by = 'row.names', all = TRUE)),

# mirTop counts
tar_target(mmu_mirna_pretrimmed_mirtop_counts, nfcore_mirtop_counts(mmu_mirna_pretrimmed_multiqc)),
tar_target(mmu_mirna_trimmed_mirtop_counts, nfcore_mirtop_counts(mmu_mirna_trimmed_multiqc)),
tar_target(mmu_mirna_all_mice_miRBase22_mirtop_counts, 
    dplyr::left_join(
        mmu_mirna_pretrimmed_mirtop_counts, mmu_mirna_trimmed_mirtop_counts, by = "miRNA")),

# mirTop isomiRs
tar_target(mmu_mirna_pretrimmed_mirtop_isomir, nfcore_mirtop_isomir(mmu_mirna_pretrimmed_multiqc)),
tar_target(mmu_mirna_trimmed_mirtop_isomir, nfcore_mirtop_isomir(mmu_mirna_trimmed_multiqc)),
tar_target(mmu_mirna_mirtop_isomir,
    left_join(
        mmu_mirna_pretrimmed_mirtop_isomir,
        mmu_mirna_trimmed_mirtop_isomir,
            by = c("seq", "mir", "mism", "add", "t5", "t3")) |>
    mutate(across(starts_with("COHP"), ~ replace_na(., 0)))
    ),
tar_target(mmu_mirna_isomir, make_isomir(mmu_mirna_mirtop_isomir, metadata_mmu),
    resources = hprcc::medium
),

# mirTrace & samtools metrics
tar_target(mmu_mirna_pretrimmed_qc, nfcore_smrna_qc(mmu_mirna_pretrimmed_multiqc)),
tar_target(mmu_mirna_trimmed_qc, nfcore_smrna_qc(mmu_mirna_trimmed_multiqc)),
tar_target(mmu_mirna_all_mice_miRBase22_mirna_qc, rbind(mmu_mirna_pretrimmed_qc, mmu_mirna_trimmed_qc)),

# 10X
# Load data
# Multiplexed samples
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT,
    seurat_import_cellplex("^AML.scRNAseq.2022$", metadata_mmu_prepub),
    resources = hprcc::medium
),

# Single sample libraries
tar_target(sample_sheet_mir142ko_10x, make_scrnaseq_sample_sheet("^CML.10x.mir142ko$", metadata_mmu_prepub)),
tar_target(sample_sheet_blastcrisis_10x, parse_10x_blastcrisis()),

# Link fastqs
tar_target(sample_sheet_mir142ko_10x_cache, symlink_10x_fastqs("mmu_10X_mir142_ko/fastq", sample_sheet_mir142ko_10x)),
tar_target(sample_sheet_blastcrisis_10x_cache, symlink_10x_fastqs("mmu_10X_blastcrisis/fastq", sample_sheet_blastcrisis_10x)),

# Run nf-core/scrnaseq pipeline
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_multiqc,
    run_nfcore_scrnaseq("mmu_10X_mir142_ko", sample_sheet_mir142ko_10x_cache)
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_multiqc,
    run_nfcore_scrnaseq("mmu_10X_blastcrisis", sample_sheet_blastcrisis_10x_cache)
),

# Import Cellranger output
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT,
    seurat_import_cellranger(mmu_10x_mir142ko_GENCODEm28_HLT_multiqc),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT,
    seurat_import_cellranger(mmu_10x_blastcrisis_GENCODEm28_HLT_multiqc),
    resources = hprcc::medium
),

# QC filter
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT_qc, 
    seurat_perform_cell_qc(mmu_10x_aml2022_GENCODEm28_HLT),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_qc, 
    seurat_perform_cell_qc(mmu_10x_mir142ko_GENCODEm28_HLT),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_qc, 
    seurat_perform_cell_qc(mmu_10x_blastcrisis_GENCODEm28_HLT),
    resources = hprcc::medium
),

# Integrate samples with SCTransform
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT_sct,
    seurat_sctransform(mmu_10x_aml2022_GENCODEm28_HLT_qc),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_sct,
    seurat_sctransform(mmu_10x_mir142ko_GENCODEm28_HLT_qc),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_sct,
    seurat_sctransform(mmu_10x_blastcrisis_GENCODEm28_HLT_qc),
    resources = hprcc::medium
),

# Annotate clusters & UMAP
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT_sct_clust,
    seurat_annotate_clusters_and_umap(mmu_10x_aml2022_GENCODEm28_HLT_sct),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_sct_clust,
    seurat_annotate_clusters_and_umap(mmu_10x_mir142ko_GENCODEm28_HLT_sct),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_sct_clust,
    seurat_annotate_clusters_and_umap(mmu_10x_blastcrisis_GENCODEm28_HLT_sct),
    resources = hprcc::medium
),

# Annotate cell cycle
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT_sct_clust_cc,
    seurat_annotate_cell_cycle(mmu_10x_aml2022_GENCODEm28_HLT_sct_clust),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_sct_clust_cc,
    seurat_annotate_cell_cycle(mmu_10x_mir142ko_GENCODEm28_HLT_sct_clust),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_sct_clust_cc,
    seurat_annotate_cell_cycle(mmu_10x_blastcrisis_GENCODEm28_HLT_sct_clust),
    resources = hprcc::medium
),

# Annotate cell type
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT_sct_clust_cc_ct,
    seurat_annotate_cell_type(mmu_10x_aml2022_GENCODEm28_HLT_sct_clust_cc),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_sct_clust_cc_ct,
    seurat_annotate_cell_type(mmu_10x_mir142ko_GENCODEm28_HLT_sct_clust_cc),
    resources = hprcc::medium
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_sct_clust_cc_ct,
    seurat_annotate_cell_type(mmu_10x_blastcrisis_GENCODEm28_HLT_sct_clust_cc),
    resources = hprcc::medium
),

# Publish pins
# Metadata pins
tar_target(
    metadata_mmu_pins,
    publish_metadata(metadata_mmu)
),
tar_target(
    metadata_hsa_pins,
    publish_metadata(metadata_hsa)
),

# mRNA pins
tar_target(
    mmu_mrna_all_mice_GENCODEm28_pins, publish_se(mmu_mrna_all_mice_qc_se),
    resources = hprcc::small
),
tar_target(
    mmu_mrna_techrep_GENCODEm28_pins, publish_se(mmu_mrna_techrep_qc_se),
    resources = hprcc::small
),
tar_target(
    hsa_mrna_flt3_GENCODEm28_pins, publish_se(hsa_mrna_flt3_qc_se_flt),
    resources = hprcc::small
),
tar_target(
    hsa_mrna_mds_GENCODEm28_pins, publish_se(hsa_mrna_mds_qc_se_flt),
    resources = hprcc::small
),
tar_target(
    hsa_mrna_kim_GENCODEm28_pins, publish_se(hsa_mrna_kim_qc_se_flt),
    resources = hprcc::small
),
tar_target(
    hsa_mrna_non_aml_GENCODEm28_pins, publish_se(hsa_mrna_non_aml_qc_se_flt),
    resources = hprcc::small
),

# miRNA pins
tar_target(
    mmu_mirna_all_mice_miRBase22_mirtop_counts_pins,
    publish_mirtop_counts(mmu_mirna_all_mice_miRBase22_mirtop_counts),
    resources = hprcc::small
),
tar_target(
    mmu_mirna_all_mice_miRBase22_bowtie_cpm_pins,
    publish_bowtie_cpm(mmu_mirna_all_mice_miRBase22_bowtie_cpm),
    resources = hprcc::small
),
tar_target(
    mmu_mirna_all_mice_miRBase22_mirna_qc_pins,
    publish_mirna_qc(mmu_mirna_all_mice_miRBase22_mirna_qc),
    resources = hprcc::small
),

# 10X pins
tar_target(
    mmu_10x_aml2022_GENCODEm28_HLT_pins,
    publish_seurat(mmu_10x_aml2022_GENCODEm28_HLT_sct_clust_cc_ct),
    resources = hprcc::small
),
tar_target(
    mmu_10x_mir142ko_GENCODEm28_HLT_pins,
    publish_seurat(mmu_10x_mir142ko_GENCODEm28_HLT_sct_clust_cc_ct),
    resources = hprcc::small
),
tar_target(
    mmu_10x_blastcrisis_GENCODEm28_HLT_pins,
    publish_seurat(mmu_10x_blastcrisis_GENCODEm28_HLT_sct_clust_cc_ct),
    resources = hprcc::medium
),

# Collect latest pin versions
tar_target(
    latest_published_data, collect_pin_versions(
        mmu_10x_mir142ko_GENCODEm28_HLT_pins,
        mmu_10x_aml2022_GENCODEm28_HLT_pins,
        mmu_10x_blastcrisis_GENCODEm28_HLT_pins,
        # metadata
        metadata_mmu_pins,
        metadata_hsa_pins,
        # micro RNAseq
        mmu_mirna_all_mice_miRBase22_mirtop_counts_pins,
        mmu_mirna_all_mice_miRBase22_bowtie_cpm_pins,
        mmu_mirna_all_mice_miRBase22_mirna_qc_pins,
        # bulk RNAseq
        mmu_mrna_all_mice_GENCODEm28_pins,
        mmu_mrna_techrep_GENCODEm28_pins,
        hsa_mrna_flt3_GENCODEm28_pins,
        hsa_mrna_mds_GENCODEm28_pins,
        hsa_mrna_kim_GENCODEm28_pins,
        hsa_mrna_non_aml_GENCODEm28_pins
    ),
    format = "file"
),

# Build package
tar_target(
    built_package,
    build_package(latest_published_data),
    resources = hprcc::medium
)
)