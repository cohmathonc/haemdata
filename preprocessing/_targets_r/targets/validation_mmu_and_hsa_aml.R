list(
# AML validation samples
    tar_target(sample_sheet_2017_1, parse_metadata_AML.validation.2017()),
    tar_target(sample_sheet_2020_2, parse_metadata_AML.mRNA.novaseq_validation.2020()),
    tar_target(
        sample_sheet_techrep,
        rbind(
            sample_sheet_2017_1,
            sample_sheet_2020_2
        )
    ),
# AML & MDS patients (COH Biobank); MDS from EGAD00001003891; AML patients from PRJEB27973
    tar_target(sample_sheet_2022_2, parse_metadata_AML.mRNA.HSA_FLT3.2022()),
    tar_target(sample_sheet_2022_3, parse_metadata_MDS.rnaseq.EGAD00001003891()),
    tar_target(sample_sheet_2022_4, parse_metadata_AML.PRJEB27973())
)
