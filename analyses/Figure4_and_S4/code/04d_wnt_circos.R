# S4g: three effect-size-selected regions among four native CHIR-vs-CTRL calls.
# Marker sizes are visual pointers, not genomic region widths or assay coverage.
suppressPackageStartupMessages(library(grid))
wnt <- input("wnt/dmrcate_display_data.rds")
calls <- as.data.frame(input("wnt/historical_DMRCATE_calls.rds")[["CHIR.vs.CTRL"]])
calls$retained <- abs(calls$meandiff) >= .10
calls$locus <- c("PHLDB2", "DKK2", "C1QTNF3 / AMACR", "HOXB3 / HOXB-AS3")[match(calls$seqnames, c("chr3", "chr4", "chr5", "chr17"))]
stopifnot(nrow(calls) == 4, sum(calls$retained) == 3)
save_table(calls, "S4g_complete_CHIR_discovery_and_selection")
selected <- subset(calls, retained)
cyto <- as.data.frame(input("wnt/hg38_cytobands.rds"))
names(cyto) <- c("chromosome", "start0", "end0", "band", "stain")
ci <- read.delim("data/wnt/chromosome_geometry.tsv")
selected$angle <- ci$start_degrees[match(selected$seqnames, ci$chromosome)] -
  ((selected$start - 1 + selected$end) / 2) / ci$length_bp[match(selected$seqnames, ci$chromosome)] *
  ci$span_degrees[match(selected$seqnames, ci$chromosome)]
draw_circos <- function() {
  xy <- function(angle, radius) cbind(.5 + radius * cos(angle * pi / 180),
                                     .5 + radius * sin(angle * pi / 180))
  text <- function(label, x, y, size = 10, col = "#263D4C", just = "center", face = "plain")
    grid.text(label, x, y, just = just, gp = gpar(fontfamily = "Helvetica", fontsize = size, col = col, fontface = face))
  arc <- function(a1, a2, fill, border = NA) {
    angles <- seq(a1, a2, length.out = max(3, ceiling(abs(a2 - a1) * 4)))
    pts <- rbind(xy(angles, .30), xy(rev(angles), .266))
    grid.polygon(pts[, 1], pts[, 2], gp = gpar(fill = fill, col = border, lwd = .4))
  }
  band_colors <- c(gneg = "#F1F4F5", gpos25 = "#D6DFE3", gpos50 = "#B9C7CE", gpos75 = "#9BAFB9",
                   gpos100 = "#829AA7", acen = "#A8B5BE", gvar = "#DDE4E8", stalk = "#E7ECEF")
  for (j in seq_len(nrow(ci))) {
    chromosome <- ci[j, ]
    bands <- subset(cyto, chromosome == ci$chromosome[j])
    for (k in seq_len(nrow(bands))) {
      a1 <- chromosome$start_degrees - bands$start0[k] / chromosome$length_bp * chromosome$span_degrees
      a2 <- chromosome$start_degrees - bands$end0[k] / chromosome$length_bp * chromosome$span_degrees
      arc(a1, a2, band_colors[bands$stain[k]])
    }
    arc(chromosome$start_degrees, chromosome$start_degrees - chromosome$span_degrees, NA, "#A7B8C1")
    p <- xy(chromosome$start_degrees - chromosome$span_degrees / 2, .33)
    text(sub("chr", "", chromosome$chromosome), p[1], p[2], 8, "#7F95A0")
  }
  positions <- data.frame(locus = c("DKK2", "C1QTNF3 / AMACR", "HOXB3 / HOXB-AS3"),
                          x = c(.86, .84, .17), y = c(.63, .17, .88), just = c("center", "center", "center"))
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]; loc <- positions[positions$locus == row$locus, ]
    color <- ifelse(row$locus == "DKK2", "#B17629", "#4E8CA6")
    p <- xy(row$angle, .306)
    q <- xy(row$angle, .355)
    endpoint_y <- loc$y + ifelse(loc$y < .5, .11, -.10)
    grid.lines(c(p[1], q[1], loc$x), c(p[2], q[2], endpoint_y), gp = gpar(col = color, lwd = .7))
    grid.points(p[1], p[2], pch = 21, size = unit(3.4, "mm"), gp = gpar(fill = color, col = "white"))
    label <- gsub(" / ", " /\n", row$locus, fixed = TRUE)
    text(label, loc$x, loc$y + .035, 10, color, loc$just, "bold")
    text(sprintf("min. smoothed FDR\n%.2e", row$min_smoothed_fdr), loc$x, loc$y - .035, 7.5, color, loc$just)
  }
  text("3", .5, .55, 30, face = "bold")
  text("CHIR-responsive DMRs", .5, .49, 10)
  text("|mean delta beta| >= 0.10", .5, .45, 8)
  text("hg38", .5, .415, 8)
}
circos_grob <- grid.grabExpr(draw_circos(), width = 7.5, height = 6.5, wrap = TRUE)
save_panel(cowplot::ggdraw(circos_grob), "S4g_CHIR_responsive_loci", 7.5, 6.5)
