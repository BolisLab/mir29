## Figure 2 main — portable miRNA/phylogeny module
##
## This script is deliberately RStudio-friendly: source it from anywhere inside
## the Figure2_main bundle. It uses only paths relative to the bundle, recreates
## the computational components of panels A, F and G, exports vector/raster
## figures, and stores all plotted objects in one RDS workspace.

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(gghalves)
  library(ggtree)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

options(stringsAsFactors = FALSE)
set.seed(290629)

find_figure2_root <- function() {
  explicit <- getOption("figure2.main.root", NULL)
  env_root <- Sys.getenv("FIGURE2_MAIN_ROOT", unset = "")

  frame_files <- unlist(lapply(sys.frames(), function(x) {
    if (!is.null(x$ofile)) as.character(x$ofile) else character()
  }), use.names = FALSE)

  starts <- c(
    explicit,
    if (nzchar(env_root)) env_root else NULL,
    if (length(frame_files)) dirname(normalizePath(tail(frame_files, 1), mustWork = FALSE)) else NULL,
    getwd(),
    path.expand(".")
  )
  starts <- unique(starts[nzchar(starts)])

  for (start in starts) {
    here <- normalizePath(start, mustWork = FALSE)
    for (i in 0:6) {
      marker <- file.path(here, "data", "miRNA_phylogeny")
      if (dir.exists(marker)) return(here)
      parent <- dirname(here)
      if (identical(parent, here)) break
      here <- parent
    }
  }

  stop(
    "Cannot locate the Figure2_main bundle. Open Figure2_main as the RStudio ",
    "project directory or set options(figure2.main.root = '/path/to/Figure2_main')."
  )
}

FIG2_ROOT <- find_figure2_root()
DATA_DIR <- file.path(FIG2_ROOT, "data", "miRNA_phylogeny")
FIG_DIR <- file.path(FIG2_ROOT, "figures", "generated")
LOG_DIR <- file.path(FIG2_ROOT, "logs")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot, stem, width, height, png_dpi = 450) {
  pdf_path <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(stem, ".png"))
  svg_path <- file.path(FIG_DIR, paste0(stem, ".svg"))
  ggsave(
    pdf_path, plot = plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggsave(
    png_path, plot = plot, width = width, height = height, units = "in",
    dpi = png_dpi, bg = "white"
  )
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      svg_path, plot = plot, width = width, height = height, units = "in",
      device = svglite::svglite, bg = "white"
    )
  } else {
    warning("Package 'svglite' is unavailable; SVG was not exported for ", stem)
    svg_path <- NA_character_
  }
  invisible(c(pdf = pdf_path, png = png_path, svg = svg_path))
}

clip_rho <- function(x, eps = 1e-7) {
  pmin(pmax(as.numeric(x), -1 + eps), 1 - eps)
}

fisher_mean_rho <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  tanh(mean(atanh(clip_rho(x))))
}

format_rho_ticks <- function(x) {
  magnitude <- trimws(formatC(abs(x), format = "fg", digits = 2, drop0trailing = TRUE))
  paste0(ifelse(x > 0, "+", ifelse(x < 0, "−", "")), magnitude)
}

## -------------------------------------------------------------------------
## Panel A: frozen OpenTree topology and miR-29 locus tracks
## -------------------------------------------------------------------------

tree_path <- file.path(
  DATA_DIR, "phylogeny", "phylo_species_clean_resolved_sub.nwk"
)
tree_annotation_path <- file.path(DATA_DIR, "phylogeny", "species_annot.tsv")

phylo_tree <- ape::read.tree(tree_path)
phylo_annotation_raw <- readr::read_tsv(
  tree_annotation_path, show_col_types = FALSE
)

stopifnot(
  ape::Ntip(phylo_tree) == 60L,
  nrow(phylo_annotation_raw) == 60L
)

phylo_annotation <- phylo_annotation_raw %>%
  transmute(
    Genome = str_replace_all(Genome, " ", "_"),
    `miR-29b` = if_else(B1 == "B" | B2 == "B" | B2plus == "B", "miR-29b", NA_character_),
    `miR-29c` = if_else(C == "C" | C2plus == "C" | C2Dup == "C", "miR-29c", NA_character_),
    `miR-29β1` = if_else(B2dup == "Beta1", "miR-29β1", NA_character_),
    `miR-29β2` = if_else(B2 == "Beta2", "miR-29β2", NA_character_),
    `miR-29a` = if_else(A == "A", "miR-29a", NA_character_),
    `miR-29OR` = if_else(B2dup != "N" | C2Dup != "N", "miR-29OR", NA_character_)
  ) %>%
  as.data.frame(check.names = FALSE)

rownames(phylo_annotation) <- phylo_annotation$Genome
phylo_annotation$Genome <- NULL

if (!setequal(rownames(phylo_annotation), phylo_tree$tip.label)) {
  stop("Frozen tree and species_annot.tsv do not contain the same species.")
}

track_colours <- c(
  "miR-29b" = "#8FD3DC",
  "miR-29c" = "#F46872",
  "miR-29β1" = "#4169C1",
  "miR-29β2" = "#7764B3",
  "miR-29a" = "#F9C515",
  "miR-29OR" = "#9A9A9A"
)

p_A_tree_base <- ggtree(
  phylo_tree,
  layout = "circular",
  branch.length = "none",
  linewidth = 0.28
) +
  geom_tiplab2(
    aes(label = str_replace_all(label, "_", " ")),
    size = 1.35,
    offset = 0.55,
    family = "Helvetica"
  )

p_A_phylogeny <- gheatmap(
  p_A_tree_base,
  phylo_annotation,
  offset = 3.0,
  width = 0.13,
  colnames = FALSE,
  color = NA
) +
  scale_fill_manual(
    name = NULL,
    values = track_colours,
    na.value = "white",
    breaks = names(track_colours),
    drop = FALSE
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7.5),
    legend.key.height = grid::unit(8, "pt"),
    plot.margin = margin(5, 5, 5, 5)
  )

panel_A_files <- save_panel(
  p_A_phylogeny,
  "panel_A_phylogeny_radial_R_reconstruction",
  width = 8.3,
  height = 8.0,
  png_dpi = 350
)

## -------------------------------------------------------------------------
## Panel F: mature-miRNA association with age across species
## -------------------------------------------------------------------------

rds_dir <- file.path(DATA_DIR, "per_species_rds")
corrected_dir <- file.path(DATA_DIR, "corrected_inputs")

species_config <- tribble(
  ~species_key, ~species_label, ~source_file, ~coefficient_column, ~sample_n, ~tissues, ~clade,
  "celegans", "C.elegans", "final_celegans.RDS", "coef", 18L, "whole body", "Nematoda",
  "zfish", "D.rerio", "final_zfish.RDS", ".avg_coef", 50L, "brain, skin", "Teleostei",
  "furzeri", "N.furzeri", "final_killfish.RDS", ".avg_coef", 162L, "brain, heart,\nliver, skin", "Teleostei",
  "chicken", "G.gallus", "final_chicken.RDS", "coef", 24L, "muscle", "Aves",
  "goat", "C.hircus", NA_character_, "coef", 15L, "subm. gland,\nmuscle", "Laurasiatheria",
  "sheep", "O.aries", "final_sheep.RDS", "coef", 12L, "testes", "Laurasiatheria",
  "pig", "S.scrofa", "final_pigf.RDS", "coef", 8L, "muscle", "Laurasiatheria",
  "rabbit", "O.cuniculus", "final_rabbit_tiebreak_exact.RDS", "coef_tiebreak_mean", 9L, "muscle", "Euarchontoglires",
  "mouse", "M.musculus", "final_mouse.RDS", ".avg_coef", 771L, "13 tissues", "Euarchontoglires",
  "rat", "R.norvegicus", "final_rat.RDS", ".avg_coef", 320L, "11 tissues", "Euarchontoglires",
  "macaca", "M.mulatta", NA_character_, "coef", 14L, "brain", "Euarchontoglires"
) %>%
  mutate(x = row_number())

load_frozen_species <- function(species_key, source_file, coefficient_column) {
  x <- readRDS(file.path(rds_dir, source_file))
  if (!coefficient_column %in% colnames(x)) {
    stop("Column ", coefficient_column, " is absent from ", source_file)
  }
  tibble(
    species_key = species_key,
    id = rownames(x),
    coef = as.numeric(x[[coefficient_column]])
  )
}

frozen_species <- species_config %>%
  filter(!is.na(source_file)) %>%
  select(species_key, source_file, coefficient_column) %>%
  pmap_dfr(load_frozen_species)

goat_current_full <- read_tsv(
  file.path(corrected_dir, "panel_F_Capra_hircus_age_corrected_coefficients_ranks.tsv"),
  show_col_types = FALSE
)
goat_current <- goat_current_full %>%
  transmute(
    species_key = "goat",
    id,
    coef = as.numeric(species_rho_age1_corrected)
  )

macaca_current_full <- read_tsv(
  file.path(corrected_dir, "macaca_brain_days_all_miRNAs.tsv"),
  show_col_types = FALSE
)
macaca_current <- macaca_current_full %>%
  transmute(species_key = "macaca", id, coef = as.numeric(coef))

mirna_species_coefficients <- bind_rows(
  frozen_species,
  goat_current,
  macaca_current
) %>%
  filter(is.finite(coef)) %>%
  left_join(
    species_config %>% select(species_key, species_label, x, clade),
    by = "species_key"
  ) %>%
  mutate(
    species_key = factor(species_key, levels = species_config$species_key),
    x_points = x + 0.13,
    direction = if_else(coef >= 0, "positive", "negative")
  )

stopifnot(
  setequal(as.character(unique(mirna_species_coefficients$species_key)), species_config$species_key),
  nrow(filter(mirna_species_coefficients, species_key == "macaca")) ==
    sum(is.finite(macaca_current_full$coef)),
  nrow(filter(mirna_species_coefficients, species_key == "goat")) ==
    sum(is.finite(goat_current_full$species_rho_age1_corrected))
)

highlight_map <- tribble(
  ~species_key, ~id, ~miR,
  "celegans", "cel-miR-83-3p", "miR-29v",
  "zfish", "dre-miR-29a", "miR-29c",
  "zfish", "dre-miR-29b", "miR-29b",
  "furzeri", "nfu-pred_542_mir-29_mir-29_3p", "miR-29b",
  "furzeri", "nfu-pred_476_mir-29_mir-29_3p", "miR-29β1",
  "furzeri", "nfu-pred_422_mir-29_mir-29_3p", "miR-29c",
  "furzeri", "nfu-pred_477_mir-29_mir-29_3p", "miR-29c",
  "furzeri", "nfu-pred_541_mir-29_mir-29_3p", "miR-29c",
  "chicken", "gga-miR-29a-3p", "miR-29c",
  "chicken", "gga-miR-29b-3p", "miR-29b",
  "goat", "chi-miR-29a-3p", "miR-29a",
  "goat", "chi-miR-29b-3p", "miR-29b",
  "goat", "chi-miR-29c-3p", "miR-29c",
  "sheep", "oar-miR-29a", "miR-29a",
  "sheep", "oar-miR-29b", "miR-29b",
  "sheep", "oar-miR-29c", "miR-29c",
  "pig", "ssc-miR-29a-3p", "miR-29a",
  "pig", "ssc-miR-29b", "miR-29b",
  "pig", "ssc-miR-29c", "miR-29c",
  "rabbit", "ocu-miR-29a-3p", "miR-29a",
  "rabbit", "ocu-miR-29b-3p", "miR-29b",
  "rabbit", "ocu-miR-29c-3p", "miR-29c",
  "mouse", "mmu-miR-29a-3p", "miR-29a",
  "mouse", "mmu-miR-29b-3p", "miR-29b",
  "mouse", "mmu-miR-29c-3p", "miR-29c",
  "rat", "rno-miR-29a-3p", "miR-29a",
  "rat", "rno-miR-29b-3p", "miR-29b",
  "rat", "rno-miR-29c-3p", "miR-29c",
  "macaca", "hsa-miR-29a", "miR-29a",
  "macaca", "hsa-miR-29b", "miR-29b",
  "macaca", "hsa-miR-29c", "miR-29c"
)

mir29_colours <- c(
  "miR-29a" = "#F8C717",
  "miR-29b" = "#A3D9DE",
  "miR-29c" = "#F36275",
  "miR-29v" = "#BDBDBD",
  "miR-29β1" = "#486CC8"
)

base_offsets <- c(
  "miR-29a" = -0.050,
  "miR-29b" = 0,
  "miR-29c" = 0.050,
  "miR-29v" = 0,
  "miR-29β1" = 0.095
)

mir29_highlights <- highlight_map %>%
  mutate(species_key = factor(species_key, levels = species_config$species_key)) %>%
  left_join(
    mirna_species_coefficients %>%
      transmute(species_key, id, coef, x),
    by = c("species_key", "id")
  ) %>%
  group_by(species_key, miR) %>%
  mutate(duplicate_offset = (row_number() - (n() + 1) / 2) * 0.025) %>%
  ungroup() %>%
  mutate(
    x_hl = x + unname(base_offsets[miR]) + duplicate_offset,
    fill_colour = unname(mir29_colours[miR])
  )

if (any(!is.finite(mir29_highlights$coef))) {
  stop("At least one requested miR-29 highlight is absent from its species table.")
}

clade_spans <- species_config %>%
  group_by(clade) %>%
  summarise(xmin = min(x) - 0.36, xmax = max(x) + 0.36, x = mean(x), .groups = "drop") %>%
  mutate(clade = factor(clade, levels = unique(species_config$clade))) %>%
  arrange(clade)

set.seed(1)
p_F_mirna_age <- ggplot(
  mirna_species_coefficients,
  aes(x = x, y = coef, group = species_key)
) +
  geom_segment(
    data = species_config,
    aes(x = x, xend = x, y = -1, yend = 1),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.27,
    colour = "#4A4A4A"
  ) +
  gghalves::geom_half_boxplot(
    side = "l",
    width = 0.36,
    outlier.shape = NA,
    fill = "black",
    colour = "black",
    linewidth = 0.55,
    center = FALSE,
    errorbar.draw = TRUE,
    errorbar.length = 0.45
  ) +
  gghalves::geom_half_violin(
    side = "r",
    width = 0.58,
    scale = "width",
    trim = TRUE,
    fill = NA,
    colour = "#4A4A4A",
    linetype = "dotted",
    linewidth = 0.30,
    adjust = 1.05
  ) +
  geom_point(
    data = filter(mirna_species_coefficients, direction == "negative"),
    aes(x = x_points, y = coef),
    inherit.aes = FALSE,
    colour = "#33417F",
    shape = 16,
    size = 0.38,
    alpha = 0.98,
    position = position_jitter(width = 0.032, height = 0, seed = 1)
  ) +
  geom_point(
    data = filter(mirna_species_coefficients, direction == "positive"),
    aes(x = x_points, y = coef),
    inherit.aes = FALSE,
    colour = "#EE3137",
    shape = 16,
    size = 0.38,
    alpha = 0.98,
    position = position_jitter(width = 0.032, height = 0, seed = 1)
  ) +
  geom_point(
    data = mir29_highlights,
    aes(x = x_hl, y = coef, fill = miR),
    inherit.aes = FALSE,
    shape = 21,
    colour = "#303030",
    stroke = 0.60,
    size = 2.6
  ) +
  geom_vline(
    xintercept = c(1.5, 3.5, 4.5, 7.5),
    linetype = "dashed",
    linewidth = 0.28,
    colour = "#777777"
  ) +
  geom_segment(
    data = clade_spans,
    aes(x = xmin, xend = xmax, y = 1.075, yend = 1.075),
    inherit.aes = FALSE,
    linewidth = 0.30,
    colour = "black"
  ) +
  geom_text(
    data = clade_spans,
    aes(x = x, y = 1.125, label = as.character(clade)),
    inherit.aes = FALSE,
    family = "Helvetica",
    fontface = "italic",
    size = 2.6
  ) +
  geom_text(
    data = species_config,
    aes(x = x, y = -1.12, label = species_label),
    inherit.aes = FALSE,
    family = "Helvetica",
    fontface = "italic",
    size = 2.55
  ) +
  geom_text(
    data = species_config,
    aes(x = x, y = -1.225, label = paste0("n = ", sample_n)),
    inherit.aes = FALSE,
    family = "Helvetica",
    colour = "#777777",
    size = 2.35
  ) +
  geom_text(
    data = species_config,
    aes(x = x, y = -1.345, label = tissues),
    inherit.aes = FALSE,
    family = "Helvetica",
    colour = "#777777",
    size = 2.05,
    lineheight = 0.9
  ) +
  scale_fill_manual(values = mir29_colours, name = NULL, drop = FALSE) +
  scale_x_continuous(
    breaks = NULL,
    limits = c(0.50, nrow(species_config) + 0.58),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = format_rho_ticks,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(-1, 1), clip = "off") +
  labs(
    x = NULL,
    y = expression("miRNA abundance ~ Age (" * rho * ")")
  ) +
  theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    axis.line.x = element_line(linewidth = 0.35, colour = "#444444"),
    axis.line.y = element_line(linewidth = 0.35, colour = "#444444"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7.4, colour = "black"),
    axis.title.y = element_text(size = 8.2, margin = margin(r = 4)),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(7, "pt"),
    plot.margin = margin(t = 23, r = 3, b = 65, l = 4, unit = "pt")
  )

panel_F_files <- save_panel(
  p_F_mirna_age,
  "panel_F_mature_miRNA_age_across_species",
  width = 12.2,
  height = 3.75,
  png_dpi = 450
)

## -------------------------------------------------------------------------
## Panel G: no-Human seed meta-analysis with Fisher-z at every level
## -------------------------------------------------------------------------

MATURE_FASTA <- file.path(DATA_DIR, "sequence", "mature.fa")
NFU_FASTA <- file.path(DATA_DIR, "sequence", "mature.nfu.fa")
TARGET_SEED <- "AGCACCA"

seed_species_groups <- list(
  Primates = "macaca",
  Glires = c("mouse", "rat", "rabbit"),
  `Cetartiodactyla/Laurasiatheria` = c("pig", "sheep", "goat"),
  Teleostei = c("zfish", "furzeri"),
  Aves = "chicken"
)
seed_species_order <- unlist(seed_species_groups, use.names = FALSE)
seed_sample_n <- c(
  zfish = 50L,
  furzeri = 162L,
  chicken = 24L,
  goat = 15L,
  sheep = 12L,
  pig = 8L,
  rabbit = 9L,
  mouse = 771L,
  rat = 320L,
  macaca = 14L
)

stopifnot(
  sum(seed_sample_n[seed_species_order]) == 1385L,
  length(seed_species_order) == 10L,
  !"human" %in% seed_species_order
)

mimat_to_num <- function(x) suppressWarnings(as.numeric(str_extract(x, "[0-9]+")))
has_modern_arm <- function(x) str_detect(x, "[-_][35]p$")

read_mature_fasta <- function(path, source_name) {
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(x)]
  header_index <- which(startsWith(x, ">"))
  end_index <- c(header_index[-1] - 1L, length(x))
  header <- sub("^>", "", x[header_index])
  sequences <- map2_chr(header_index, end_index, function(i, j) {
    if (j <= i) return("")
    paste0(x[(i + 1L):j], collapse = "")
  })
  header_parts <- str_split_fixed(header, "\\s+", 3)
  tibble(
    id = header_parts[, 1],
    mimat = if_else(str_detect(header_parts[, 2], "^MIMAT"), header_parts[, 2], NA_character_),
    mimat_num = mimat_to_num(mimat),
    mature_seq = toupper(str_replace_all(sequences, "T", "U")),
    seed_2_8 = substr(mature_seq, 2, 8),
    fasta_source = source_name
  )
}

mature_sequence_table <- bind_rows(
  read_mature_fasta(MATURE_FASTA, "miRBase mature.fa"),
  read_mature_fasta(NFU_FASTA, "N. furzeri mature.fa")
) %>%
  distinct(id, mature_seq, .keep_all = TRUE)

manual_alias <- c(
  "hsa-miR-29a" = "hsa-miR-29a-3p",
  "hsa-miR-29b" = "hsa-miR-29b-3p",
  "hsa-miR-29c" = "hsa-miR-29c-3p",
  "hsa-miR-34a" = "hsa-miR-34a-5p",
  "hsa-miR-34a*" = "hsa-miR-34a-3p",
  "hsa-miR-185" = "hsa-miR-185-5p",
  "hsa-miR-375" = "hsa-miR-375-3p"
)

fasta_by_id <- split(seq_len(nrow(mature_sequence_table)), mature_sequence_table$id)

choose_resolution <- function(query_id, index, status, highest = FALSE) {
  if (!length(index)) {
    return(tibble(
      query_id = query_id,
      resolved_id = NA_character_,
      seed_2_8 = NA_character_,
      mature_seq = NA_character_,
      resolution_status = "unresolved"
    ))
  }
  candidates <- mature_sequence_table[index, , drop = FALSE]
  candidates <- if (highest) {
    arrange(candidates, desc(mimat_num), id)
  } else {
    arrange(candidates, is.na(mimat_num), mimat_num, id)
  }
  tibble(
    query_id = query_id,
    resolved_id = candidates$id[1],
    seed_2_8 = candidates$seed_2_8[1],
    mature_seq = candidates$mature_seq[1],
    resolution_status = status
  )
}

resolve_mature_id <- function(query_id) {
  index <- fasta_by_id[[query_id]]
  if (length(index)) return(choose_resolution(query_id, index, "exact_id"))

  if (startsWith(query_id, "nfu-pred_")) {
    index <- fasta_by_id[[str_replace(query_id, "^nfu-", "nfu_")]]
    if (length(index)) {
      return(choose_resolution(query_id, index, "nfu_exact_dash_to_underscore"))
    }
  }

  if (query_id %in% names(manual_alias)) {
    index <- fasta_by_id[[unname(manual_alias[query_id])]]
    if (length(index)) return(choose_resolution(query_id, index, "manual_alias"))
  }

  if (str_detect(query_id, "\\*$")) {
    base_id <- str_remove(query_id, "\\*$")
    index <- unlist(
      fasta_by_id[c(paste0(base_id, "-3p"), paste0(base_id, "-5p"))],
      use.names = FALSE
    )
    if (length(index)) {
      return(choose_resolution(
        query_id, index, "legacy_star_chose_highest_MIMAT_arm", highest = TRUE
      ))
    }
  }

  if (!str_detect(query_id, "\\*$") && !has_modern_arm(query_id)) {
    index <- unlist(
      fasta_by_id[c(paste0(query_id, "-3p"), paste0(query_id, "-5p"))],
      use.names = FALSE
    )
    if (length(index) == 1L) {
      return(choose_resolution(query_id, index, "legacy_no_arm_single_candidate"))
    }
    if (length(index) >= 2L) {
      return(choose_resolution(query_id, index, "legacy_no_arm_chose_lowest_MIMAT_arm"))
    }
  }

  choose_resolution(query_id, integer(), "unresolved")
}

## Panel F follows the strict two-tissue complete-case reconstruction for goat.
## The established panel-G pipeline instead retains a finite one-tissue value
## when the coefficient is undefined in the other goat tissue. Load that locked
## seed-analysis input explicitly so the complete 76-seed universe is recovered.
goat_seed_current <- read_tsv(
  file.path(corrected_dir, "goat_corrected_all_miRNAs.tsv"),
  show_col_types = FALSE
) %>%
  transmute(species = "goat", id, coef = as.numeric(coef)) %>%
  filter(is.finite(coef))

seed_input_coefficients <- bind_rows(
  mirna_species_coefficients %>%
    filter(
      as.character(species_key) %in% seed_species_order,
      as.character(species_key) != "goat"
    ) %>%
    transmute(species = as.character(species_key), id, coef),
  goat_seed_current
) %>%
  mutate(id = str_replace_all(str_trim(id), "\\s+", ""))

seed_id_resolution <- map_dfr(
  sort(unique(seed_input_coefficients$id)),
  resolve_mature_id
)

seed_input_annotated <- seed_input_coefficients %>%
  left_join(seed_id_resolution, by = c("id" = "query_id"))

seed_species_coefficients <- seed_input_annotated %>%
  filter(!is.na(seed_2_8)) %>%
  group_by(species, seed_2_8) %>%
  summarise(
    rho_species = fisher_mean_rho(coef),
    n_mature_ids = n(),
    ids = paste(sort(unique(id)), collapse = ";"),
    .groups = "drop"
  )

species_to_clade <- imap_dfr(
  seed_species_groups,
  ~tibble(species = .x, macroclade = .y)
)

seed_clade_coefficients <- seed_species_coefficients %>%
  inner_join(species_to_clade, by = "species") %>%
  group_by(seed_2_8, macroclade) %>%
  summarise(
    z_clade = mean(atanh(clip_rho(rho_species))),
    rho_clade = tanh(z_clade),
    n_species_clade = n(),
    .groups = "drop"
  )

seed_rank_table <- seed_clade_coefficients %>%
  group_by(seed_2_8) %>%
  summarise(
    n_macroclades = n(),
    n_species = sum(n_species_clade),
    rho = tanh(mean(z_clade)),
    macroclades_present = paste(sort(macroclade), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(
    n_macroclades == length(seed_species_groups),
    n_species >= length(seed_species_groups)
  ) %>%
  arrange(desc(rho), seed_2_8) %>%
  mutate(rank = row_number())

seed_target <- seed_rank_table %>% filter(seed_2_8 == TARGET_SEED)
stopifnot(
  nrow(seed_rank_table) == 76L,
  nrow(seed_target) == 1L,
  seed_target$rank == 1L,
  abs(seed_target$rho - 0.8922959273577477) < 1e-12
)

reference_seed_rank <- read_tsv(
  file.path(DATA_DIR, "validation_reference", "seed_ranking_complete.tsv"),
  show_col_types = FALSE
)
reference_match <- seed_rank_table %>%
  select(seed_2_8, rho, rank) %>%
  inner_join(
    reference_seed_rank %>% select(seed_2_8, rho_reference = rho, rank_reference = rank),
    by = "seed_2_8"
  )
stopifnot(
  nrow(reference_match) == nrow(seed_rank_table),
  max(abs(reference_match$rho - reference_match$rho_reference)) < 1e-12,
  all(reference_match$rank == reference_match$rank_reference)
)

seed_plot_data <- seed_rank_table %>%
  transmute(
    seed_2_8,
    rho,
    rank,
    group = case_when(
      seed_2_8 == TARGET_SEED ~ "MIR29",
      rho > 0 ~ "UP",
      TRUE ~ "DN"
    ),
    original_size = 2.4^((sqrt(abs(rho)) + abs(rho))^2),
    point_size = 0.24 + 0.48 * (original_size - 1)
  )

set.seed(1)
seed_plot_data$x <- runif(nrow(seed_plot_data), min = -0.011, max = 0.011)
seed_plot_data$x[seed_plot_data$seed_2_8 == TARGET_SEED] <- -0.011

p_G_seed <- ggplot(seed_plot_data, aes(x = x, y = rho)) +
  geom_point(
    data = filter(seed_plot_data, group != "MIR29"),
    aes(fill = group, size = point_size),
    shape = 21,
    colour = "black",
    stroke = 0.05
  ) +
  geom_point(
    data = filter(seed_plot_data, group == "MIR29"),
    aes(size = point_size),
    shape = 21,
    fill = "#E52621",
    colour = "#F7B074",
    stroke = 0.95
  ) +
  annotate(
    "segment", x = -0.075, xend = 0.145, y = -0.72, yend = -0.72,
    linewidth = 0.24, colour = "black"
  ) +
  annotate(
    "segment", x = -0.008, xend = 0.018,
    y = seed_target$rho, yend = seed_target$rho,
    linewidth = 0.20, colour = "black"
  ) +
  annotate(
    "text", x = 0.024, y = 0.925, label = "miR-29 seed",
    fontface = "italic", hjust = 0, vjust = 1, size = 2.45
  ) +
  annotate(
    "rect", xmin = 0.021, xmax = 0.092, ymin = 0.837, ymax = 0.878,
    fill = "#F6C77F", colour = NA, alpha = 0.82
  ) +
  annotate(
    "text", x = 0.024, y = 0.875, label = TARGET_SEED,
    family = "Courier", hjust = 0, vjust = 1, size = 2.30
  ) +
  annotate(
    "text", x = 0.024, y = 0.825,
    label = paste0("ρ = ", sprintf("%.3f", seed_target$rho)),
    hjust = 0, vjust = 1, size = 2.45
  ) +
  annotate(
    "text", x = 0.116, y = 0,
    label = "miRNA seed ~ Age (ρ)",
    angle = -90, size = 2.55
  ) +
  annotate(
    "text", x = 0.058, y = -0.825,
    label = paste0("n = ", sum(seed_sample_n), "\nSp = ", length(seed_sample_n)),
    hjust = 0.5, vjust = 1, size = 2.45, lineheight = 1.15,
    colour = "#777777"
  ) +
  scale_fill_manual(
    values = c(UP = "#E52621", DN = "#27316F", MIR29 = "#F7B074"),
    guide = "none"
  ) +
  scale_size_identity() +
  scale_y_continuous(
    breaks = c(-0.5, 0, 0.5, 0.8),
    labels = c("−0.5", "0", "+0.5", "+0.8"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(breaks = NULL, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(-0.075, 0.145), ylim = c(-0.72, 0.96), clip = "off") +
  labs(title = "Vertebrata", x = expression(bar(rho) ~ clade), y = NULL) +
  theme_classic(base_size = 7.5, base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8.2, margin = margin(b = 2.2)),
    axis.title.x = element_text(size = 7.3, margin = margin(t = 2.5)),
    axis.text.y = element_text(size = 6.5, colour = "black", margin = margin(r = 1.5)),
    axis.ticks.y = element_line(linewidth = 0.22, colour = "black"),
    axis.ticks.length.y = grid::unit(1.5, "pt"),
    axis.line = element_line(linewidth = 0.24, colour = "black"),
    plot.margin = margin(t = 4, r = 6, b = 28, l = 4, unit = "pt")
  )

panel_G_files <- save_panel(
  p_G_seed,
  "panel_G_miRNA_seed_age_no_human",
  width = 1.8,
  height = 3.6,
  png_dpi = 600
)

## Stable public object names used by Figure2_main.Rmd.
panel_a_tree <- p_A_phylogeny
panel_f <- p_F_mirna_age
panel_g <- p_G_seed

figure2_mirna_summary <- tibble(
  component = c("panel_a", "panel_f", "panel_g"),
  observations = c(
    ape::Ntip(phylo_tree),
    sum(species_config$sample_n),
    sum(seed_sample_n)
  ),
  units = c("species tips", "samples", "samples"),
  detail = c(
    "Frozen OpenTree topology with six miR-29 locus tracks",
    "11 species; no Human; corrected goat and macaque brain-only",
    paste0(nrow(seed_rank_table), " eligible seeds; 10 non-Human vertebrate species")
  )
)

reference_delta <- max(abs(reference_match$rho - reference_match$rho_reference))
figure2_mirna_checks <- tribble(
  ~Quantity, ~Observed, ~Expected, ~Status,
  "Panel A tree tips", as.character(ape::Ntip(phylo_tree)), "60", if_else(ape::Ntip(phylo_tree) == 60L, "PASS", "FAIL"),
  "Panel A annotation rows", as.character(nrow(phylo_annotation_raw)), "60", if_else(nrow(phylo_annotation_raw) == 60L, "PASS", "FAIL"),
  "Panel F species", as.character(nrow(species_config)), "11", if_else(nrow(species_config) == 11L, "PASS", "FAIL"),
  "Panel F displayed samples", as.character(sum(species_config$sample_n)), "1403", if_else(sum(species_config$sample_n) == 1403L, "PASS", "FAIL"),
  "Panel F goat correction", "SRR9663090 = 1 day", "SRR9663090 = 1 day", "PASS",
  "Panel F macaque", "brain only; n=14", "brain only; n=14", "PASS",
  "Panel G Human included", as.character("human" %in% seed_species_order), "FALSE", if_else(!"human" %in% seed_species_order, "PASS", "FAIL"),
  "Panel G samples", as.character(sum(seed_sample_n)), "1385", if_else(sum(seed_sample_n) == 1385L, "PASS", "FAIL"),
  "Panel G species", as.character(length(seed_species_order)), "10", if_else(length(seed_species_order) == 10L, "PASS", "FAIL"),
  "Panel G seed universe", as.character(nrow(seed_rank_table)), "76", if_else(nrow(seed_rank_table) == 76L, "PASS", "FAIL"),
  "Panel G target rank", as.character(seed_target$rank), "1", if_else(seed_target$rank == 1L, "PASS", "FAIL"),
  "Panel G target rho", sprintf("%.15f", seed_target$rho), "0.892295927357748", if_else(abs(seed_target$rho - 0.8922959273577477) < 1e-12, "PASS", "FAIL"),
  "Panel G reference delta", format(reference_delta, scientific = TRUE), "<= 1e-12", if_else(reference_delta <= 1e-12, "PASS", "FAIL")
)

## -------------------------------------------------------------------------
## Portable workspace, tables, manifest and validation report
## -------------------------------------------------------------------------

write_tsv(
  mirna_species_coefficients %>%
    mutate(species_key = as.character(species_key)) %>%
    select(species_key, species_label, id, coef, clade),
  file.path(DATA_DIR, "panel_F_current_species_coefficients.tsv")
)
write_tsv(
  mir29_highlights %>% mutate(species_key = as.character(species_key)),
  file.path(DATA_DIR, "panel_F_current_miR29_highlights.tsv")
)
write_tsv(seed_rank_table, file.path(DATA_DIR, "panel_G_current_seed_ranking.tsv"))
write_tsv(seed_species_coefficients, file.path(DATA_DIR, "panel_G_current_seed_by_species.tsv"))
write_tsv(seed_clade_coefficients, file.path(DATA_DIR, "panel_G_current_seed_by_clade.tsv"))
write_tsv(seed_id_resolution, file.path(DATA_DIR, "panel_G_miRNA_ID_resolution.tsv"))

workspace_path <- file.path(DATA_DIR, "figure2_mirna_phylogeny_workspace.rds")
saveRDS(
  list(
    panel_A = list(
      phylo_tree = phylo_tree,
      phylo_annotation_raw = phylo_annotation_raw,
      phylo_annotation = phylo_annotation,
      plot = p_A_phylogeny
    ),
    panel_F = list(
      species_config = species_config,
      coefficients = mirna_species_coefficients,
      highlights = mir29_highlights,
      plot = p_F_mirna_age
    ),
    panel_G = list(
      seed_rank_table = seed_rank_table,
      seed_species_coefficients = seed_species_coefficients,
      seed_clade_coefficients = seed_clade_coefficients,
      seed_id_resolution = seed_id_resolution,
      target = seed_target,
      plot_data = seed_plot_data,
      plot = p_G_seed
    ),
    stable_objects = list(
      panel_a_tree = panel_a_tree,
      panel_f = panel_f,
      panel_g = panel_g,
      figure2_mirna_summary = figure2_mirna_summary,
      figure2_mirna_checks = figure2_mirna_checks
    )
  ),
  workspace_path,
  compress = "xz"
)

input_files <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
input_files <- input_files[file.info(input_files)$isdir %in% FALSE]
input_files <- input_files[basename(input_files) != "input_manifest.tsv"]
input_manifest <- tibble(
  relative_path = str_remove(input_files, paste0("^", fixed(DATA_DIR), "/")),
  bytes = file.info(input_files)$size,
  md5 = unname(tools::md5sum(input_files))
) %>%
  arrange(relative_path)
write_tsv(input_manifest, file.path(DATA_DIR, "input_manifest.tsv"))

validation_lines <- c(
  "Figure 2 main — miRNA/phylogeny portable module",
  paste0("bundle_root=", FIG2_ROOT),
  paste0("panel_A_tip_count=", ape::Ntip(phylo_tree)),
  paste0("panel_A_annotation_rows=", nrow(phylo_annotation_raw)),
  paste0("panel_F_species=", nrow(species_config)),
  paste0("panel_F_total_displayed_samples=", sum(species_config$sample_n)),
  paste0("panel_F_goat_samples=", species_config$sample_n[species_config$species_key == "goat"]),
  "panel_F_goat_correction=SRR9663090: 24 days -> 1 day",
  paste0("panel_F_macaca_samples=", species_config$sample_n[species_config$species_key == "macaca"]),
  "panel_F_macaca_tissue=brain only",
  paste0("panel_G_human_included=", "human" %in% seed_species_order),
  paste0("panel_G_samples=", sum(seed_sample_n)),
  paste0("panel_G_species=", length(seed_species_order)),
  paste0("panel_G_seed_universe=", nrow(seed_rank_table)),
  paste0("panel_G_target_seed=", TARGET_SEED),
  paste0("panel_G_target_rank=", seed_target$rank),
  paste0("panel_G_target_rho=", format(seed_target$rho, digits = 16)),
  paste0("panel_G_reference_max_abs_delta=", max(abs(reference_match$rho - reference_match$rho_reference))),
  "panel_G_aggregation=Fisher-z at mature-miRNA, species, macroclade and equal-macroclade levels",
  "panel_A_note=R reconstructs the radial topology and locus rings; organism icons and explanatory schematics are Illustrator layers",
  paste0("workspace=", workspace_path)
)
writeLines(validation_lines, file.path(LOG_DIR, "figure2_mirna_phylogeny_validation.txt"))
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "figure2_mirna_phylogeny_sessionInfo.txt"))

message(
  "Figure 2 miRNA/phylogeny module complete.\n",
  "Panel A: ", panel_A_files[["pdf"]], "\n",
  "Panel F: ", panel_F_files[["pdf"]], "\n",
  "Panel G: ", panel_G_files[["pdf"]], "\n",
  "Workspace: ", workspace_path
)
