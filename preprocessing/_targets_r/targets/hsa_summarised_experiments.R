list(
    tar_target(hsa_mrna_flt3_qc_se, annotate_se(hsa_mrna_flt3_qc_se_bare, metadata_hsa, hsa_mrna_flt3_qc)),
    tar_target(hsa_mrna_mds_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_mds_qc), metadata_hsa, hsa_mrna_mds_qc)),
    tar_target(hsa_mrna_kim_qc_se, annotate_se(get_rnaseq_se(hsa_mrna_kim_qc), metadata_hsa, hsa_mrna_kim_qc))
)
