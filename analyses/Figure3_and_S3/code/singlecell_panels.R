# The embedding and bulk projections are preserved, rather than re-estimated.
singlecell <- readRDS("data/singlecell/reference_plot_data.rds")
cells <- singlecell$cells
cell_groups <- list(
  Fibroblast = c("Fibro_homeostatic", "Fibro_inflammatory_remodeling"),
  Epithelial = c("Epithelial_basal", "Epithelial_differentiated"),
  Endothelial = c("Endothelial_activated_ACKR1", "Endothelial_general", "Endothelial_capillary", "Endothelial_arterial", "Endothelial_lymphatic"),
  Immune = c("T_NK_cytotoxic", "T_cells_IL7R", "B_cells", "Plasma_cells", "Myeloid_inflammatory", "Mast_cells"),
  SMC = "Mural_pericyte_SMC"
)
cells$macro_group <- "Other"
for (group in names(cell_groups)) cells$macro_group[cells$celltype %in% cell_groups[[group]]] <- group
cell_palette <- c(Fibroblast = "#E51C23", Epithelial = "#F2C500", Endothelial = "#762A83", Immune = "#377EB8", SMC = "#B98245", Other = "#BDBDBD")
umap_plot <- ggplot(cells, aes(x, y)) +
  geom_point(colour = "#222222", size = .25) +
  geom_point(aes(colour = macro_group), size = .15) +
  scale_colour_manual(values = cell_palette) + coord_equal() + theme_void() +
  theme(legend.position = "bottom") + labs(colour = NULL)
bulk <- singlecell$bulk
projection_plot <- ggplot(cells, aes(x, y)) + geom_point(size = .15, colour = "#DDDDDD") +
  geom_point(data = bulk, aes(proj_UMAP_1, proj_UMAP_2, colour = treatment), size = 1.6) +
  scale_colour_manual(values = c(control = "grey45", inhibitor = "#328A37", mimic = "#923D8E")) +
  coord_equal() + theme_void() + theme(legend.position = "bottom") + labs(colour = NULL)

# Pooled macrocluster summaries retain cell-count weighting.
marker_checkpoint <- readRDS("data/singlecell/marker_dotplot_macro_identity.rds")
markers <- marker_checkpoint$data
markers$macro_group <- factor(markers$macro_group,
  levels = c("SMC / perivascular", "Immune", "Endothelial", "Epithelial", "Fibroblast"))
markers$pct_detected <- 100 * markers$n_detected / markers$n_cells
markers$expression_z <- ave(log1p(markers$mean_normalized_counts), markers$gene,
  FUN = function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else rep(0, length(x)))
stopifnot(max(abs(markers$expression_z - marker_checkpoint$data$expression_z)) < 1e-10)
marker_order <- unique(markers$gene[order(-markers$gene_y)])
markers$gene <- factor(markers$gene, levels = rev(marker_order))
marker_plot <- ggplot(markers, aes(gene, macro_group)) +
  geom_point(aes(size = pct_detected, colour = expression_z), shape = 16) +
  scale_size_area(max_size = 4, limits = c(0, 100), breaks = c(10, 50, 100), name = "Cells expressing (%)") +
  scale_colour_gradient2(low = "#3B6FA3", mid = "#F2F2F2", high = "#BD3C4A", limits = c(-2.5, 2.5), oob = scales::squish, name = "Expression z score") +
  labs(x = NULL, y = NULL) + theme_figure(9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
stopifnot(nrow(cells) == 133335L, nrow(bulk) == 96L,
          sum(unique(markers[, c("macro_group", "n_cells")])[, "n_cells"]) == 131475L)
