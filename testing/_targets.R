options(hprcc.slurm_logs = TRUE, hprcc.default_partition = "fast,all")
# _targets.R
library(targets)
library(tarchetypes)
library(tidyverse)
library(Seurat)
library(hprcc)

tar_option_set(
    error = "continue",
    packages = c("tidyverse", "Seurat"),
    storage = "main",
    retrieval = "main"
)

get_feature_abundance <- function(dir, feature_names) {
    # Initialize result at the start
    result <- data.frame(library_id = basename(dir) %>% str_remove("^sample-"))

    message("Sample:", result$library_id)
    matrix_path <- file.path(dir, "outs", "filtered_feature_bc_matrix")
    counts <- Read10X(data.dir = matrix_path)
    message("read...")
    seurat_obj <- CreateSeuratObject(counts = counts)
    message("created SeuratObj")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
    message("normalised...")

# Get all features at once
features_present <- feature_names[feature_names %in% rownames(seurat_obj)]
if (length(features_present) > 0) {
    # Get raw counts and colSums for normalization 
    feature_data <- GetAssayData(seurat_obj, slot = "counts")[features_present, , drop = FALSE]
    lib_sizes <- colSums(GetAssayData(seurat_obj, slot = "counts"))
    
    # Calculate pseudobulk (sum counts, then normalize)
    for (feature in feature_names) {
        result[[paste0(feature, "_abundance")]] <- 
            if (feature %in% features_present) {
                sum(feature_data[feature, ]) / sum(lib_sizes) * 1e6  # Scale to CPM
            } else {
                0L
            }
    }
} else {
    result[paste0(feature_names, "_abundance")] <- NA_real_
}

return(result)
}

# Use shortened base paths
BASE_PATHS <- list(
    blastcrisis_2024 = "/home/domeally/workspaces/haemdata/analysis_cache/mmu_10X_blastcrisis/nfcore-scrnaseq-v2.1.0-mmu/cellranger/count",
    blastcrisis_m28 = "/home/domeally/workspaces/haemdata/analysis_cache/mmu_10X_blastcrisis_backup/nfcore-scrnaseq-v2.1.0-mmu/cellranger/count",
    mir142ko_2024 = "/home/domeally/workspaces/haemdata/analysis_cache/mmu_10X_mir142_ko/nfcore-scrnaseq-v2.1.0-mmu/cellranger/count",
    mir142ko_m28 = "/home/domeally/workspaces/haemdata/analysis_cache/mmu_10X_mir142_ko_backup/nfcore-scrnaseq-v2.1.0-mmu/cellranger/count"
)

# Create all combinations of parameters with shorter names
configs <- expand.grid(
    cohort = c("blast_crisis", "mir142ko"),
    genome = c("2024", "m28"), # Shortened genome names
    feature_name = c("hsa-BCR-ABL-gene", "HSA_BCR_ABL1_gene", "HSA_BCR_gene", "HSA_ABL1_gene"),
    stringsAsFactors = FALSE
) %>%
    filter(
        (genome == "2024" & feature_name == "hsa-BCR-ABL-gene") |
            (genome == "m28" & feature_name %in% c("HSA_BCR_ABL1_gene", "HSA_BCR_gene", "HSA_ABL1_gene"))
    ) %>%
    mutate(
        path_key = paste0(
            if_else(cohort == "blast_crisis", "blastcrisis", "mir142ko"),
            "_",
            genome
        ),
        base_dir = unlist(BASE_PATHS[path_key])
    )

# Group configs by cohort and genome to process multiple features together
configs_grouped <- configs %>%
    group_by(cohort, genome, base_dir) %>%
    summarise(
        feature_names = list(feature_name),
        .groups = "drop"
    )

# Define pipeline with shorter target names
list(
    tar_map(
        values = configs_grouped,
        names = c("cohort", "genome"),
        tar_target(
            name = dirs,
            command = list.dirs(base_dir, full.names = TRUE, recursive = FALSE) %>%
                grep("sample-", ., value = TRUE),
            deployment = "main"
        ),
        tar_target(
            name = abundances,
            command = get_feature_abundance(dirs, feature_names),
            pattern = map(dirs),
            resources = small
        ),
        tar_target(
            name = df,
            command = bind_rows(abundances),
            deployment = "main"
        ),
        tar_target(
            name = results,
            command = write.csv(
                df,
                file = sprintf("%s_%s_abundances.csv", cohort, genome),
                row.names = FALSE
            ),
            deployment = "main",
            format = "file"
        ),
        # New target for generating plots
        tar_target(
            name = abundance_plot,
            command = {
                # Read metadata
                metadata <- readRDS("metadata.rds")

                # Join with metadata
                plot_data <- left_join(df, metadata)

                # Set up abundance column based on genome version
                if (genome == "2024") {
                    plot_data$abundance <- plot_data$`hsa-BCR-ABL-gene_abundance`
                } else {
                    plot_data$abundance <- plot_data$HSA_BCR_ABL1_gene_abundance
                }

                # Clean data and add small offset to zeros
                plot_data <- plot_data %>%
                    filter(!is.na(abundance)) %>%
                    mutate(abundance = ifelse(abundance == 0,
                        min(abundance[abundance > 0], na.rm = TRUE) / 10,
                        abundance
                    ))

                if (cohort == "mir142ko") {
                    # For mir142ko cohort - create boxplot with points by tissue and condition
                    p <- ggplot(plot_data, aes(
                        x = tissue,
                        y = abundance,
                        fill = treatment, # Using treatment (KO/WT) for fill
                        color = treatment # And for point colors
                    )) +
                        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
                        geom_jitter(
                            position = position_jitterdodge(dodge.width = 0.75),
                            size = 3, alpha = 0.7
                        ) +
                        theme_classic() +
                        theme(
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            panel.grid.major.y = element_line(color = "grey90"),
                            legend.position = "right"
                        ) +
                        #scale_y_log10() + # Log scale for better visualization of differences
                        labs(
                            x = "Tissue",
                            y = "BCR-ABL Gene Abundance (CPM)",
                            fill = "Treatment",
                            color = "Treatment",
                            title = sprintf("BCR-ABL Gene Abundance by Tissue\n%s - %s", cohort, genome)
                        )
                } else {
                    # Keep original time series plot for blast crisis cohort
                    p <- ggplot(plot_data, aes(
                        x = timepoint,
                        y = abundance,
                        color = treatment,
                        group = mouse_id
                    )) +
                        geom_line() +
                        geom_point(size = 3) +
                        geom_text(
                            data = plot_data %>%
                                group_by(mouse_id) %>%
                                slice_max(timepoint, n = 1),
                            aes(label = mouse_id),
                            hjust = -0.2,
                            size = 3
                        ) +
                        theme_classic() +
                        theme(
                            panel.grid.major.y = element_line(color = "grey90"),
                            legend.position = "right"
                        ) +
                        labs(
                            x = "Timepoint",
                            y = "BCR-ABL Gene Abundance (CPM)",
                            color = "Treatment",
                            title = sprintf("BCR-ABL Gene Abundance Over Time\n%s - %s", cohort, genome)
                        )
                }

                # Save plot with consistent dimensions
                ggsave(
                    sprintf("%s_%s_abundance_plot.png", cohort, genome),
                    p,
                    width = 10,
                    height = 8,
                    dpi = 300
                )
            },
            deployment = "main",
            format = "file"
        )
    )
)