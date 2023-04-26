library(Seurat)
library(scCustomize)

seurat_object <- get_pin("mmu_10x_blastcrisis_GENCODEm28_HLT.rds")
metadata <- readxl::read_xlsx("data-raw/metadata_mmu_blastcrisis.xlsx") |>
    rename(orig.ident = library_id) |>
    mutate(weeks = readr::parse_number(timepoint))


igc_metadata <- tibble::tribble(
     ~orig.ident, ~read_format, ~ilab_order, ~read_depth,
    "COHP_50358", "28+10+10+101", 4282L, 118.361248,
    "COHP_50359", "28+10+10+101", 4282L, 106.889749,
    "COHP_50360", "28+10+10+101", 4282L, 102.877744,
    "COHP_50361", "28+10+10+101", 4282L, 104.732072,
    "COHP_50362", "28+10+10+101", 4282L, 117.152266,
    "COHP_50363", "28+10+10+101", 4282L, 106.602778,
    "COHP_50546", "28+10+10+101", 4293L, 30.766402,
    "COHP_50547", "28+10+10+101", 4293L, 35.994608,
    "COHP_50548", "28+10+10+101", 4293L, 32.60602,
    "COHP_50549", "28+10+10+101", 4293L, 43.368361,
    "COHP_50550", "28+10+10+101", 4293L, 32.708344,
    "COHP_50551", "28+10+10+101", 4293L, 37.35342,
    "COHP_50620", "28+10+10+101", 4300L, 22.543784,
    "COHP_50621", "28+10+10+101", 4300L, 34.424001,
    "COHP_50622", "28+10+10+101", 4300L, 37.207438,
    "COHP_50623", "28+10+10+101", 4300L, 32.72098,
    "COHP_50624", "28+10+10+101", 4300L, 40.773399,
    "COHP_50625", "28+10+10+101", 4300L, 28.351739,
    "COHP_50680", "28+10+10+101", 4306L, 24.233,
    "COHP_50681", "28+10+10+101", 4306L, 30.07,
    "COHP_50682", "28+10+10+101", 4306L, 30.976,
    "COHP_50683", "28+10+10+101", 4306L, 27.263,
    "COHP_50684", "28+10+10+101", 4306L, 41.505,
    "COHP_50685", "28+10+10+101", 4306L, 27.038,
    "COHP_50686", "28+10+10+101", 4308L, 24.563,
    "COHP_50687", "28+10+10+101", 4308L, 29.554,
    "COHP_50688", "28+10+10+101", 4308L, 32.108,
    "COHP_50689", "28+10+10+101", 4308L, 28.182,
    "COHP_50690", "28+10+10+101", 4308L, 33.583,
    "COHP_50691", "28+10+10+101", 4308L, 30.536,
    "COHP_50738", "28+10+10+101", 4318L, 32.122,
    "COHP_50739", "28+10+10+101", 4318L, 35.145,
    "COHP_50740", "28+10+10+101", 4318L, 47.273,
    "COHP_50741", "28+10+10+101", 4318L, 36.336,
    "COHP_50742", "28+10+10+101", 4318L, 32.114,
    "COHP_50858", "28+10+10+101", 4330L, 21.628,
    "COHP_50859", "28+10+10+101", 4330L, 24.068,
    "COHP_50860", "101+8+8+101", 4330L, 22.821,
    "COHP_50861", "101+8+8+101", 4330L, 25.706,
    "COHP_50862", "101+8+8+101", 4330L, 25.354,
    "COHP_50976", "101+8+8+101", 4340L, 26.441,
    "COHP_50977", "101+8+8+101", 4340L, 27.259,
    "COHP_50978", "101+8+8+101", 4340L, 23.513,
    "COHP_50979", "101+8+8+101", 4340L, 19.358,
    "COHP_50980", "101+8+8+101", 4340L, 20.333,
    "COHP_51077", "101+8+8+101", 4348L, 34.33,
    "COHP_51078", "101+8+8+101", 4353L, 40.59,
    "COHP_51079", "101+8+8+101", 4353L, 36.69,
    "COHP_51110", "101+8+8+101", 4353L, 38.7,
    "COHP_51111", "101+8+8+101", 4364L, 40.35,
    "COHP_51112", "101+8+8+101", 4364L, 43.34,
    "COHP_51168", "101+8+8+101", 4364L, 31.94,
    "COHP_51169", "101+8+8+101", 4364L, 28.67
)

sample_median <- Median_Stats(seurat_object, group_by_var = "orig.ident")

sample_cellcount <- seurat_object@meta.data |> add_count(orig.ident) |> select(orig.ident, cell_count = n) |> distinct()

sample_metadata <- metadata |>
    left_join(igc_metadata) |>
    left_join(sample_median) |>
    left_join(sample_cellcount) |>
    select(-dod) # NAs not liked by Seurat

obj <- Add_Sample_Meta(
    seurat_object, sample_metadata, join_by_seurat = "orig.ident", join_by_meta = "orig.ident"
)

FeaturePlot_scCustom(seurat_object, features = "nFeature_RNA")
FeaturePlot_scCustom(obj, features = "cell_count")

FeaturePlot(obj, features = "nFeature_RNA", split.by = "read_format")

FeaturePlot(obj, features = "Median_nFeature_RNA", split.by = "read_format")


DimPlot_scCustom(obj, group.by = "read_format")
DimPlot(obj, group.by = "read_format")

DimPlot_scCustom(obj, group.by = "ilab_order")
DimPlot(obj, group.by = "ilab_order")


DimPlot_All_Samples(obj, meta_data_column = "ilab_order", num_col = 3, pt.size = 0.5)
DimPlot_All_Samples(obj, meta_data_column = "orig.ident", num_col = 3, pt.size = 0.5)

DimPlot_scCustom(obj, group.by = "ilab_order")
DimPlot_scCustom(obj, split.by = "ilab_order", raster = FALSE)
DimPlot_scCustom(obj, split.by = "read_format", raster = FALSE)


FeaturePlot_scCustom(obj, features = "read_depth", na_cutoff = 90, alpha_exp = 0.75)
FeaturePlot_scCustom(obj, features = "read_depth", na_cutoff = NULL, raster = FALSE)
FeaturePlot_scCustom(obj, features = "weeks", alpha_exp = 0.2)

Stacked_VlnPlot(
    obj, features = "timepoint", x_lab_rotate = TRUE
)

Idents(object = obj) <- "read_format"

Stacked_VlnPlot(
    obj, features = c("cell_count", "read_depth", "nCount_RNA", "nFeature_RNA", "percent_mt", "percent_ribo", "Median_percent_ribo"), x_lab_rotate = TRUE)


# Get cell names
high_rd <- WhichCells(obj, expression = read_depth > 100)

# Make into list
cells <- list(high_readdepth = high_rd)

# Plot
Cell_Highlight_Plot(obj, cells_highlight = cells)
