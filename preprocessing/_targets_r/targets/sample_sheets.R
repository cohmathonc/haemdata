list(
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
    )
)
