list(
    tar_target(mmu_mrna_2016_2022_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_2016_2022_qc), metadata_mmu_prepub, mmu_mrna_2016_2022_qc
    )),
    tar_target(mmu_mrna_techrep_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_techrep_qc), metadata_mmu_prepub, mmu_mrna_techrep_qc
    )),
    tar_target(mmu_mrna_cml3_qc_se, annotate_se(
        get_rnaseq_se(mmu_mrna_cml3_qc), metadata_mmu_prepub, mmu_mrna_cml3_qc
    ))
)
