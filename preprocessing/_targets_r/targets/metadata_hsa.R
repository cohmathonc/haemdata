tar_target(
    metadata_hsa, make_metadata_hsa(
        dplyr::full_join(sample_sheet_2022_3, sample_sheet_2022_2) |>
        dplyr::full_join(sample_sheet_2022_4))
    )
