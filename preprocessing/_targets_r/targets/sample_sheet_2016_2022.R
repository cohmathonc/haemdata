list(
    tar_target(pattern_2016_2022, "^AML.mRNA.2016$|^AML.mRNA.2018.all_samples$|^AML.mRNA.2020$|
        |^AML.mRNA.2021.RxGroup1$|^AML.mRNA.2021.RxGroup2$|^AML.mRNA.2021.RxGroup2_pt2$|
        |^AML.mRNA.2022.RxGroup3$|^CML.mRNA.2021$|^CML.mRNA.2022$"), 
    
    tar_target(sample_sheet_2016_2022, published_metadata_mmu |>
        dplyr::filter(str_detect(cohort, pattern_2016_2022)) |>
        dplyr::select(library_id, fastq_1, fastq_2, strandedness)),

    tar_target(sample_sheet_CML_3, published_metadata_mmu |>
        dplyr::filter(str_detect(cohort, "^CML.mRNA.2022_pt2$")) |>
        dplyr::select(library_id, fastq_1, fastq_2, strandedness))
)        
