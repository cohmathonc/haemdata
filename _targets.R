# Created by use_targets().
# Follow the manual to check and run the pipeline:
# https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

###### PIPELINE SETUP ######
# nf-core pipeline version
rnaseq_release <- "3.7"
###########################
# Load the R scripts & functions:
for (file in list.files("scripts", full.names = TRUE)) source(file)

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
# Set target options:
tar_option_set(
    packages = c("haemdata", "tidyverse", "SummarizedExperiment"), # packages that targets need to run
    # imports = c("haemdata"), # packages that targets need to run
    error = "continue", # continue or stop on error
    # format = "qs", # default storage format
    # # Set other options as needed.
    storage = "worker",
    retrieval = "worker"
    # garbage_collection = TRUE,
    # memory = "transient"
)

# check that we're running on Apollo
if (!file.exists(nf_core_cache)) {
    stop("This pipeline must be run on Apollo")
}

list(
    # make the package logo
    tar_target(logo, make_logo()),

    ###############################################################################################
    # Build sample sheets -------------------------------------------------
    # AML
    tar_target(sample_sheet_2016_1, parse_metadata_AML.mRNA.2016()),
    tar_target(sample_sheet_2017_1, parse_metadata_AML.validation.2017()),
    tar_target(sample_sheet_2018_1, parse_metadata_AML.mRNA.2018.all_samples()),
    tar_target(sample_sheet_2020_1, parse_metadata_AML.mRNA.2020()),
    tar_target(sample_sheet_2020_2, parse_metadata_AML.mRNA.novaseq_validation.2020()),
    tar_target(sample_sheet_2021_1, parse_metadata_AML.mRNA.2021.RxGroup1()),
    tar_target(sample_sheet_2021_2, parse_metadata_AML.mRNA.2021.RxGroup2()),
    tar_target(sample_sheet_2021_3, parse_metadata_AML.mRNA.2021.RxGroup2_pt2()),
    tar_target(sample_sheet_2021_4, parse_metadata_AML.mRNA.2022.RxGroup3()),
    # CML
    tar_target(sample_sheet_CML_1, parse_metadata_CML.mRNA.2021()),
    tar_target(sample_sheet_CML_2, parse_metadata_CML.mRNA.2022()),
    tar_target(sample_sheet_CML_3, parse_metadata_CML.mRNA.2022_pt2()),
    # AML scRNAseq
    tar_target(sample_sheet_2022_1, parse_metadata_AML.scRNAseq.2022()),
    # Patient data, mRNAseq
    # FLT3 AML patients (COH Biobank); MDS from EGAD00001003891; AML patients from PRJEB27973
    tar_target(sample_sheet_2022_2, parse_metadata_AML.mRNA.HSA_FLT3.2022()),
    tar_target(sample_sheet_2022_3, parse_metadata_MDS.rnaseq.EGAD00001003891()),
    tar_target(sample_sheet_2022_4, parse_metadata_AML.PRJEB27973()),

    # Make sample sheets for the nf-core/rnaseq pipeline as described in the manual:
    # https://nf-co.re/rnaseq/usage#samplesheet-input
    # CML mice
    tar_target(
        sample_sheet_CML,
        rbind(sample_sheet_CML_1, sample_sheet_CML_2)
    ),
    # 2016_2022 mice (no validation sets included)
    tar_target(
        sample_sheet_2016_2022,
        rbind(
            sample_sheet_2016_1,
            sample_sheet_2018_1,
            sample_sheet_2020_1,
            sample_sheet_2021_1,
            sample_sheet_2021_2,
            sample_sheet_2021_3,
            sample_sheet_2021_4,
            sample_sheet_CML_1,
            sample_sheet_CML_2,
            sample_sheet_2022_1
        )
    ),

    # Make a sample sheet for the validation sets 2017_1 & 2020_2
    # which are essentially technical replicates and so are not included with
    # the remaining samples
    tar_target(
        sample_sheet_techrep,
        rbind(
            sample_sheet_2017_1,
            sample_sheet_2020_2
        )
    ),

    ###############################################################################################
    # run the nf-core pipeline - mRNA -------------------------------------------
    # Names adhere to the format "species_protocol_cohort_ref-genome_workflow"
    # cohort: {2016_2022, validation, flt3, mds, 2022_1}
    # species: {mmu|hsa}
    # protocol: {mrna|10x|mirna}
    # ref-genome: {GENCODEr40|GENCODEm28_HLT|GRCm38_HLT}
    # workflow: {qc|salmon|seurat}}
    # Construct the `run_folder` parameter for run_nf_core_rnaseq() as "species_protocol_cohort"

    # GENCODEm28_HLT Full QC pipeline
    tar_target(mmu_mrna_2016_2022_qc,
        run_nf_core_rnaseq("mmu_mrna_2016_2022", sample_sheet_2016_2022, "GENCODEm28_HLT"),
        format = "file"
    ),
    tar_target(mmu_mrna_techrep_qc,
        run_nf_core_rnaseq("mmu_mrna_techrep", sample_sheet_techrep, "GENCODEm28_HLT"),
        format = "file"
    ),
    tar_target(mmu_mrna_aml2016_qc,
        run_nf_core_rnaseq("mmu_mrna_aml2016", sample_sheet_2016_1, "GENCODEm28_HLT"),
        format = "file"
    ),
    tar_target(mmu_mrna_cml3_qc,
        run_nf_core_rnaseq("mmu_mrna_cml3_qc", sample_sheet_CML_3, "GENCODEm28_HLT"),
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
    # Salmon only
    tar_target(CML.mRNA_salmon,
        run_nf_core_rnaseq("mmu_mrna_cml", sample_sheet_CML, "GENCODEm28_HLT", qc = FALSE),
        format = "file"
    ),

    # run the nf-core pipeline GRCm38_HLT
    # Salmon only
    tar_target(mmu_mrna_aml2016_GRCm38_HLT_salmon,
        run_nf_core_rnaseq("mmu_mrna_aml2016", sample_sheet_2016_1, "GRCm38_HLT", qc = FALSE),
        format = "file"
    ),
    tar_target(CML.mRNA_GRCm38_HLT_salmon,
        run_nf_core_rnaseq("mmu_mrna_cml", sample_sheet_CML, "GRCm38_HLT", qc = FALSE),
        format = "file"
    ),
    tar_target(mmu_mrna_2016_2022_GRCm38_HLT_salmon,
        run_nf_core_rnaseq("mmu_mrna_2016_2022", sample_sheet_2016_2022, "GRCm38_HLT", qc = FALSE),
        format = "file"
    ),

    ############################################################################
    # Consolidate metadata -----------------------------------------------------
    ############################################################################
    # consolidate metadata across all mice samples
    tar_target(metadata_mmu_prepub, make_metadata_mmu(rbind(sample_sheet_2016_2022, sample_sheet_CML_3))),
    # consolidate metadata across human samples
    tar_target(metadata_hsa, make_metadata_hsa(
        dplyr::full_join(sample_sheet_2022_3, sample_sheet_2022_2) |>
        dplyr::full_join(sample_sheet_2022_4)
    )),

    ###############################################################################################
    # make SummarisedExperiments -----------------------------------------------------
    # from each nf-core/rnaseq pipeline run, annotate the salmon/salmon.merged.gene_counts.rds
    # SummarizedExperiment with sample metadata and QC metrics
    # GENCODEm28_HLT
    tar_target(mmu_mrna_2016_2022_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_2016_2022_qc), metadata_mmu_prepub, mmu_mrna_2016_2022_qc
    )),
    tar_target(
        mmu_mrna_cml_qc_se,
        subset(mmu_mrna_2016_2022_qc_se,
            select = grepl("CML", SummarizedExperiment::colData(mmu_mrna_2016_2022_qc_se)$project)
        )
    ),
    tar_target(
        mmu_mrna_aml_qc_se,
        subset(mmu_mrna_2016_2022_qc_se,
            select = grepl("AML", SummarizedExperiment::colData(mmu_mrna_2016_2022_qc_se)$project)
        )
    ),
    tar_target(mmu_mrna_aml2016_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_aml2016_qc), metadata_mmu_prepub, mmu_mrna_aml2016_qc
    )),
    tar_target(mmu_mrna_techrep_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_techrep_qc), metadata_mmu_prepub, mmu_mrna_techrep_qc
    )),
    tar_target(mmu_mrna_cml3_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_cml3_qc), metadata_mmu_prepub, mmu_mrna_cml3_qc
    )),

    # GRCm38_HLT 
    tar_target(mmu_mrna_2016_2022_GRCm38_HLT_salmon_se, annotate_se(
        get_rnaseq_se(mmu_mrna_2016_2022_GRCm38_HLT_salmon), metadata_mmu_prepub, mmu_mrna_2016_2022_GRCm38_HLT_salmon
    )),
    tar_target(
        mmu_mrna_cml_GRCm38_HLT_salmon_se,
        subset(mmu_mrna_2016_2022_GRCm38_HLT_salmon_se,
            select = grepl("CML", SummarizedExperiment::colData(mmu_mrna_2016_2022_GRCm38_HLT_salmon_se)$project)
        )
    ),
    tar_target(
        mmu_mrna_aml_GRCm38_HLT_salmon_se,
        subset(mmu_mrna_2016_2022_GRCm38_HLT_salmon_se,
            select = grepl("AML", SummarizedExperiment::colData(mmu_mrna_2016_2022_GRCm38_HLT_salmon_se)$project)
        )
    ),

    # GENCODEr40
    tar_target(hsa_mrna_flt3_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_flt3_qc), metadata_hsa, hsa_mrna_flt3_qc)),
    tar_target(hsa_mrna_mds_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_mds_qc), metadata_hsa, hsa_mrna_mds_qc)),
    tar_target(hsa_mrna_kim_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_kim_qc), metadata_hsa, hsa_mrna_kim_qc)),

    ###############################################################################################
    # Merge SummarisedExperiments
    tar_target(mmu_mrna_all_mice_qc_se, merge_mrna_se(mmu_mrna_2016_2022_qc_se, mmu_mrna_cml3_qc_se,
        new_name = "mmu_mrna_all_mice_GENCODEm28_HLT_qc")),
    # ###############################################################################################
    # Amend metadata_mmu with qc_pass_mapping column
    tar_target(metadata_mmu, flag_lowly_mapped_mmu_se(mmu_mrna_all_mice_qc_se, metadata_mmu_prepub)),
    # ###############################################################################################
    # Drop hsa samples that fail mapping threshold  -- for SummarisedExperiments generated from qc runs
    tar_target(hsa_mrna_flt3_qc_se_flt, qc_filter_se(hsa_mrna_flt3_qc_se)),
    tar_target(hsa_mrna_mds_qc_se_flt, qc_filter_se(hsa_mrna_mds_qc_se)),
    tar_target(hsa_mrna_kim_qc_se_flt, qc_filter_se(hsa_mrna_kim_qc_se)),

    ###############################################################################################
    ### find outliers
    # tar_target(mmu_mrna_2016_2022_qc_se_outliers, find_outliers_se(mmu_mrna_2016_2022_qc_se),
    #     resources = apollo_large),
    # tar_target(mmu_mrna_cml_qc_se_outliers, find_outliers_se(mmu_mrna_cml_qc_se),
    #     resources = apollo_large),
    # tar_target(mmu_mrna_aml_qc_se_outliers, find_outliers_se(mmu_mrna_aml_qc_se),
    #     resources = apollo_large),
    # tar_target(mmu_mrna_aml2016_qc_se_outliers, find_outliers_se(mmu_mrna_aml2016_qc_se),
    #     resources = apollo_large),
    # tar_target(hsa_mrna_flt3_qc_se_outliers, find_outliers_se(hsa_mrna_flt3_qc_se),
    #     resources = apollo_large)
    ###############################################################################################
    # Publish mRNAseq metadata (with pins package) ----------------------------------------------
    ###############################################################################################
    # save SummarisedExperiments to disk
    tar_target(mmu_mrna_all_mice_GENCODEm28_pins, publish_se(mmu_mrna_all_mice_qc_se)),
    tar_target(mmu_mrna_techrep_GENCODEm28_pins, publish_se(mmu_mrna_techrep_qc_se)),
    tar_target(hsa_mrna_flt3_GENCODEm28_pins, publish_se(hsa_mrna_flt3_qc_se_flt)),
    tar_target(hsa_mrna_mds_GENCODEm28_pins, publish_se(hsa_mrna_mds_qc_se_flt)),
    tar_target(hsa_mrna_kim_GENCODEm28_pins, publish_se(hsa_mrna_kim_qc_se_flt)),
    tar_target(metadata_mmu_pins, publish_metadata(metadata_mmu)),
    tar_target(metadata_hsa_pins, publish_metadata(metadata_hsa)),

    ###############################################################################################
    # single cell RNA-seq ------------------------------------------------------------------------
    # load data
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT, seurat_import_objects("GENCODEm28_HLT"),
        resources = apollo_medium
    ),
    tar_target(mmu_10x_2022_1_GRCm38_HLT, seurat_import_objects("GRCm38_HLT"),
        resources = apollo_medium
    ),

    # # QC filter ---------------------------------------------------------------
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT_qc, seurat_perform_cell_qc(mmu_10x_2022_1_GENCODEm28_HLT),
        resources = apollo_medium
    ),
    tar_target(mmu_10x_2022_1_GRCm38_HLT_qc, seurat_perform_cell_qc(mmu_10x_2022_1_GRCm38_HLT),
        resources = apollo_medium
    ),

    # # # Integrate cells with SCTransform
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT_sct, seurat_sctransform(mmu_10x_2022_1_GENCODEm28_HLT_qc),
        resources = apollo_bigmem
    ),
    tar_target(mmu_10x_2022_1_GRCm38_HLT_sct, seurat_sctransform(mmu_10x_2022_1_GRCm38_HLT_qc),
        resources = apollo_bigmem
    ),
    # # Make clusters  ---------------------------------------------------------------
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust, seurat_annotate_clusters_and_umap(
        mmu_10x_2022_1_GENCODEm28_HLT_sct
    ),
    resources = apollo_large
    ),
    tar_target(mmu_10x_2022_1_GRCm38_HLT_sct_clust, seurat_annotate_clusters_and_umap(
        mmu_10x_2022_1_GRCm38_HLT_sct
    ),
    resources = apollo_large
    ),
    # # Annotate cell cycle ---------------------------------------------------------
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc, seurat_annotate_cell_cycle(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust),
        resources = apollo_medium
    ),
    tar_target(mmu_10x_2022_1_GRCm38_HLT_sct_clust_cc, seurat_annotate_cell_cycle(mmu_10x_2022_1_GRCm38_HLT_sct_clust),
        resources = apollo_medium
    ),

    # # Annotate cell type -----------------------------------------------------------
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc_ct, seurat_annotate_cell_cycle(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc),
        resources = apollo_medium
    ),
    tar_target(mmu_10x_2022_1_GRCm38_HLT_sct_clust_cc_ct, seurat_annotate_cell_cycle(mmu_10x_2022_1_GRCm38_HLT_sct_clust_cc),
        resources = apollo_medium
    ),

    # Publish Seurat objects -----------------------------------------------------------
    tar_target(mmu_10x_2022_1_GENCODEm28_HLT_pins, publish_seurat(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc_ct),
    resources = apollo_small),
    tar_target(mmu_10x_2022_1_GRCm38_HLT_pins, publish_seurat(mmu_10x_2022_1_GRCm38_HLT_sct_clust_cc_ct),
    resources = apollo_small),

    ######### Collect latest pin versions #########
    # TODO make a function to prune pins not in a release
    tar_target(
        latest_published_data,
        haemdata::write_data("published_pins", rbind(
            # 10X single cell RNA-seq
            mmu_10x_2022_1_GRCm38_HLT_pins,
            mmu_10x_2022_1_GENCODEm28_HLT_pins,
            # metadata
            metadata_mmu_pins,
            metadata_hsa_pins,
            # bulk RNAseq
            mmu_mrna_all_mice_GENCODEm28_pins,
            mmu_mrna_techrep_GENCODEm28_pins,
            hsa_mrna_flt3_GENCODEm28_pins,
            hsa_mrna_mds_GENCODEm28_pins,
            hsa_mrna_kim_GENCODEm28_pins
        ))
    ),
            tar_target(built_package, build_package(latest_published_data), resources = apollo_shortmem)
)
