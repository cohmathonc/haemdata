# Created by use_targets().
# Follow the manual to check and run the pipeline:
# https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

###### PIPELINE SETUP ######
# nf-core pipeline version
rnaseq_release <- "3.7"
# Post chat to teams with new data
teams <- FALSE
###########################

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
# Set target options:
tar_option_set(
    packages = c("tidyverse", "SummarizedExperiment"), # packages that targets need to run
    # format = "qs", # default storage format
    # # Set other options as needed.
    storage = "worker",
    retrieval = "worker"
    # garbage_collection = TRUE,
    # memory = "transient"
)

# Load the R scripts & functions:
for (file in list.files("R", full.names = TRUE)) source(file)

list(
    # make the package logo
    tar_target(logo, make_logo()),

    ###############################################################################################
    # Build sample sheets
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
    # AML scRNAseq
    tar_target(sample_sheet_2022_1, parse_metadata_AML.scRNAseq.2022()),
    # Human FLT3 patients - COH Biobank & MDS from EGAD00001003891
    tar_target(sample_sheet_2022_2, parse_metadata_AML.mRNA.HSA_FLT3.2022()),
    tar_target(sample_sheet_2022_3, parse_metadata_MDS.rnaseq.EGAD00001003891()),

    # combine sample sheets for the nf-core/rnaseq pipeline
    tar_target(
        sample_sheet_CML,
        rbind(sample_sheet_CML_1, sample_sheet_CML_2)
    ),
    # all mice (no validation sets included)
    tar_target(
        sample_sheet_all_mice,
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

    # Make a "validation" run with samples from the validation sets 2017, 2022
    # these are essentially technical replicates and so are not included with
    # the remaining samples
    tar_target(
        sample_sheet_validation_mice,
        rbind(
            sample_sheet_2017_1,
            sample_sheet_2020_2
        )
    ),

    ###############################################################################################
    # run the nf-core pipeline
    # Names adhere to the format "cohort_species_protocol_ref-genome_workflow"
    # cohort: {all, validation, flt3, mds}
    # species: {mmu|hsa}
    # protocol: {mrna|scrna|mirna}
    # ref-genome: {GENCODEr40|GENCODEm28_HLT|GRCm38_HLT}
    # workflow: {qc|salmon}}
    # Name the "run_folder" as "species_protocol_cohort"

    # GENCODEm28_HLT Full QC pipeline
    tar_target(all_mice.mRNA_qc,
        run_nf_core_rnaseq("mmu_mrna_all", sample_sheet_all_mice, "GENCODEm28_HLT"),
        format = "file"
    ),
    tar_target(validation_qc,
        run_nf_core_rnaseq("mmu_mrna_validation", sample_sheet_validation_mice, "GENCODEm28_HLT"),
        format = "file"
    ),
    tar_target(AML.mRNA.HSA_FLT3_qc,
        run_nf_core_rnaseq("hsa_mrna_flt3", sample_sheet_2022_2, "GENCODEr40"),
        format = "file"
    ),
    tar_target(hsa_mrna_mds_qc,
        run_nf_core_rnaseq("hsa_mrna_mds", sample_sheet_2022_3, "GENCODEr40"),
        format = "file"
    ),
    tar_target(AML.mRNA.2016_qc,
        run_nf_core_rnaseq("mmu_mrna_aml2016", sample_sheet_2016_1, "GENCODEm28_HLT"),
        format = "file"
    ),

    # Salmon only
    tar_target(CML.mRNA_salmon,
        run_nf_core_rnaseq("mmu_mrna_cml", sample_sheet_CML, "GENCODEm28_HLT", qc = FALSE),
        format = "file"
    ),

    # run the nf-core pipeline GRCm38_HLT
    # Salmon only
    tar_target(AML.mRNA.2016_salmon_GRCm38_HLT,
        run_nf_core_rnaseq("mmu_mrna_aml2016", sample_sheet_2016_1, "GRCm38_HLT", qc = FALSE),
        format = "file"
    ),
    tar_target(CML.mRNA_salmon_GRCm38_HLT,
        run_nf_core_rnaseq("mmu_mrna_cml", sample_sheet_CML, "GRCm38_HLT", qc = FALSE),
        format = "file"
    ),
    tar_target(all_mice.mRNA_salmon_GRCm38_HLT,
        run_nf_core_rnaseq("mmu_mrna_all", sample_sheet_all_mice, "GRCm38_HLT", qc = FALSE),
        format = "file"
    ),

    ############################################################################
    # Consolidate metadata -----------------------------------------------------
    ############################################################################
    # consolidate metadata across all mice samples
    tar_target(metadata_mmu, make_metadata_mmu(sample_sheet_all_mice)),
    # consolidate metadata across all human samples
    tar_target(metadata_hsa, make_metadata_hsa(
        dplyr::full_join(sample_sheet_2022_3, sample_sheet_2022_2)
    )),

    ###############################################################################################
    # make SummarisedExperiment from each pipeline run, annotate with metadata and QC metrics
    # GENCODEm28_HLT
    tar_target(all_mice.mRNA_qc_se, annotate_se(
        get_rnaseq_se(all_mice.mRNA_qc), metadata_mmu, all_mice.mRNA_qc
    )),
    tar_target(
        CML.all_mice.mRNA_qc_se,
        subset(all_mice.mRNA_qc_se,
            select = grepl("CML", SummarizedExperiment::colData(all_mice.mRNA_qc_se)$project)
        )
    ),
    tar_target(
        AML.all_mice.mRNA_qc_se,
        subset(all_mice.mRNA_qc_se,
            select = grepl("AML", SummarizedExperiment::colData(all_mice.mRNA_qc_se)$project)
        )
    ),
    tar_target(AML.mRNA.2016_qc_se, annotate_se(
        get_rnaseq_se(AML.mRNA.2016_qc), metadata_mmu, AML.mRNA.2016_qc
    )),
    tar_target(validation_qc_se, annotate_se(
        get_rnaseq_se(validation_qc), metadata_mmu, validation_qc
    )),

    # GRCm38_HLT
    tar_target(all_mice.mRNA_salmon_GRCm38_HLT_se, annotate_se(
        get_rnaseq_se(all_mice.mRNA_salmon_GRCm38_HLT), metadata_mmu, all_mice.mRNA_salmon_GRCm38_HLT
    )),
    tar_target(
        CML.all_mice.mRNA_salmon_GRCm38_HLT_se,
        subset(all_mice.mRNA_salmon_GRCm38_HLT_se,
            select = grepl("CML", SummarizedExperiment::colData(all_mice.mRNA_salmon_GRCm38_HLT_se)$project)
        )
    ),
    tar_target(
        AML.all_mice.mRNA_salmon_GRCm38_HLT_se,
        subset(all_mice.mRNA_salmon_GRCm38_HLT_se,
            select = grepl("AML", SummarizedExperiment::colData(all_mice.mRNA_salmon_GRCm38_HLT_se)$project)
        )
    ),

    # GENCODEr40
    tar_target(AML.mRNA.HSA_FLT3_qc_se, annotate_se(get_rnaseq_se(AML.mRNA.HSA_FLT3_qc), metadata_hsa, AML.mRNA.HSA_FLT3_qc)),
    tar_target(hsa_mrna_mds_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_mds_qc), metadata_hsa, hsa_mrna_mds_qc)),

    ###############################################################################################
    # Drop samples that fail mapping threshold  -- only for SummarisedExperiments generated from qc runs
    tar_target(AML.mRNA.2016_qc_se_flt, qc_filter_se(AML.mRNA.2016_qc_se)),
    tar_target(all_mice.mRNA_qc_se_flt, qc_filter_se(all_mice.mRNA_qc_se)),
    tar_target(validation_qc_se_flt, qc_filter_se(validation_qc_se)),
    tar_target(AML.mRNA.HSA_FLT3_qc_se_flt, qc_filter_se(AML.mRNA.HSA_FLT3_qc_se)),
    tar_target(hsa_mrna_mds_qc_se_flt, qc_filter_se(hsa_mrna_mds_qc_se)),
    ###############################################################################################
    ### find outliers
    # tar_target(all_mice.mRNA_qc_se_outliers, find_outliers_se(all_mice.mRNA_qc_se),
    #     resources = apollo_large),
    # tar_target(CML.all_mice.mRNA_qc_se_outliers, find_outliers_se(CML.all_mice.mRNA_qc_se),
    #     resources = apollo_large),
    # tar_target(AML.all_mice.mRNA_qc_se_outliers, find_outliers_se(AML.all_mice.mRNA_qc_se),
    #     resources = apollo_large),
    # tar_target(AML.mRNA.2016_qc_se_outliers, find_outliers_se(AML.mRNA.2016_qc_se),
    #     resources = apollo_large),
    # tar_target(AML.mRNA.HSA_FLT3_qc_se_outliers, find_outliers_se(AML.mRNA.HSA_FLT3_qc_se),
    #     resources = apollo_large)

    ###############################################################################################
    # Publish mRNAseq metadata (with pins package)
    ###############################################################################################
    # save SummarisedExperiments to disk
    tar_target(all_mice_GENCODEm28_pins, publish_se(all_mice.mRNA_qc_se_flt)),
    tar_target(validation_mice_GENCODEm28_pins, publish_se(validation_qc_se_flt)),
    tar_target(HSA_COH_FLT3_GENCODEm28_pins, publish_se(AML.mRNA.HSA_FLT3_qc_se_flt)),
    tar_target(HSA_EGA_MDS_GENCODEm28_pins, publish_se(hsa_mrna_mds_qc_se_flt)),

    tar_target(metadata_mmu_pins, publish_metadata(metadata_mmu)),
    tar_target(metadata_hsa_pins, publish_metadata(metadata_hsa)),

    ###############################################################################################
    # single cell RNA-seq
    # load data
    tar_target(GENCODEm28_HLT_seurat, seurat_import_objects("GENCODEm28_HLT"),
        resources = apollo_medium),
    tar_target(GRCm38_HLT_seurat, seurat_import_objects("GRCm38_HLT"),
        resources = apollo_medium),

    # # QC filter ---------------------------------------------------------------
    tar_target(GENCODEm28_HLT_seurat_qc, seurat_perform_cell_qc(GENCODEm28_HLT_seurat),
        resources = apollo_bigmem),
    tar_target(GRCm38_HLT_seurat_qc, seurat_perform_cell_qc(GRCm38_HLT_seurat),
        resources = apollo_bigmem),

    # # Integrate cells using Reciprocal PCA -------------------------------------
    tar_target(GENCODEm28_HLT_seurat_rpca, seurat_integrate_RPCA(GENCODEm28_HLT_seurat_qc),
        resources = apollo_large),
    tar_target(GRCm38_HLT_seurat_rpca, seurat_integrate_RPCA(GRCm38_HLT_seurat_qc),
        resources = apollo_large),

    # # # Integrate cells with SCTransform
    # # tar_target(GENCODEm28_HLT_seurat_sct, integrate_SCTransform(GENCODEm28_HLT_seurat_qc)),
    # # tar_target(GRCm38_HLT_seurat_sct, integrate_SCTransform(GRCm38_HLT_seurat_qc)),

    # # Make clusters  ---------------------------------------------------------------
    tar_target(GENCODEm28_HLT_seurat_rpca_clust, seurat_annotate_umap_and_clusters(GENCODEm28_HLT_seurat_rpca),
        resources = apollo_large),
    tar_target(GRCm38_HLT_seurat_rpca_clust, seurat_annotate_umap_and_clusters(GRCm38_HLT_seurat_rpca),
        resources = apollo_large),

    # # Annotate cell cycle ---------------------------------------------------------
    tar_target(GENCODEm28_HLT_seurat_rpca_clust_cc, seurat_annotate_cell_cycle(GENCODEm28_HLT_seurat_rpca_clust),
        resources = apollo_large),
    tar_target(GRCm38_HLT_seurat_rpca_clust_cc, seurat_annotate_cell_cycle(GRCm38_HLT_seurat_rpca_clust),
        resources = apollo_large),

    # # Export objects ---------------------------------------------------------------------
    # tar_target(GENCODEm28_HLT_seurat_rpca_h5, export_seurat_object_h5ad(GENCODEm28_HLT_seurat_rpca)),
    # tar_target(GRCm38_HLT_seurat_rpca_h5, export_seurat_object_h5ad(GRCm38_HLT_seurat_rpca))


    ######### Collect latest pin versions #########
    tar_target(
        latest_published_data,
        helpeRs::write_data("published_pins", rbind(
            metadata_mmu_pins,
            metadata_hsa_pins,
            all_mice_GENCODEm28_pins,
            validation_mice_GENCODEm28_pins,
            HSA_COH_FLT3_GENCODEm28_pins
        ))
    )
)