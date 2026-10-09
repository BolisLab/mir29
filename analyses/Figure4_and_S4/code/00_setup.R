# Shared configuration. Run from the bundle root, never from the source archive.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(cowplot)
})
COLORS <- c(CTRL = "#686868", INH = "#248636", MIM = "#923d8e",
            CHIR = "#427DB5", COMBO = "#C57B2D")
CS <- c(INH = "#248636", inhCTRL = "#A8CDAE", MIM = "#923D8E",
        mimCTRL = "#D4ADCE", CTRL = "#686868")
TP_COLORS <- c("#e47c64", "#d02d18", "#2fb5b5", "#00736d", "#9881ba", "#56338a")
ORDINAL_X <- c(1, 2, 4, 5, 7, 8)
HOURS <- c(72, 96, 168, 192, 264, 288)
panels <- list()
for (d in c("figures/reconstructed", "tables", "validation"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
theme_figure <- function(size = 10) {
  theme_classic(base_size = size, base_family = "Helvetica") +
    theme(plot.title = element_text(size = size, face = "bold"),
          legend.title = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), legend.position = "bottom")
}
theme_set(theme_figure())
base_probe <- function(ids) sub("_[^_]+$", "", ids)
input <- function(file) readRDS(file.path("data", file))
metadata <- function() input("primary/sample_metadata.rds")
# Optional upstream regeneration requires the original large matrices.
# Normal rendering uses the documented frozen analytical checkpoints.
read_input <- function(file) {
  folder <- Sys.getenv("FIG4_RAW_DATA", unset = "")
  if (!nzchar(folder)) stop("Set FIG4_RAW_DATA to the original matrix directory for upstream recalculation.")
  readRDS(file.path(folder, file))
}
cache <- function(name, calculation) {
  path <- file.path("data/primary", paste0(name, ".rds"))
  if (!file.exists(path)) stop("Missing frozen checkpoint: ", path,
                             ". Restore it or run the documented upstream procedure explicitly.")
  readRDS(path)
}
save_table <- function(x, name) {
  write.table(as.data.frame(x), file.path("tables", paste0(name, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
save_panel <- function(plot, name, width = 6, height = 4) {
  panels[[name]] <<- plot
  path <- file.path("figures/reconstructed", name)
  ggsave(paste0(path, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(paste0(path, ".png"), plot, width = width, height = height,
         dpi = 160, bg = "white", limitsize = FALSE)
  message("Reconstructed ", name)
  invisible(plot)
}
assert_close <- function(actual, expected, label, tolerance = 1e-8) {
  error <- max(abs(actual - expected), na.rm = TRUE)
  if (!is.finite(error) || error > tolerance) stop(label, ": validation failed (", error, ")")
  cat(label, "\tPASS\t", format(error, digits = 12), "\n",
      file = "validation/numeric_checks.tsv", append = TRUE)
  invisible(TRUE)
}
# Quadratic curves are descriptive. The ribbon interpolates pointwise mean +/-
# SEM, not a regression confidence band. Controls have solid markers/error bars.
trajectory <- function(d, ylab, xpos = ORDINAL_X, title = NULL) {
  d <- as.data.frame(d)
  d$x <- xpos[as.integer(d$TP2)]
  sm <- d %>% group_by(series, x, TP2) %>%
    summarise(mean = mean(score), sem = sd(score) / sqrt(n()), n = n(), .groups = "drop")
  stopifnot(all(sm$n >= 3), !anyNA(sm$sem))
  fitted <- sm %>% group_by(series) %>% group_modify(~ {
    degree <- min(2, nrow(.x) - 1)
    model <- lm(mean ~ poly(x, degree, raw = TRUE), .x)
    grid <- seq(min(.x$x), max(.x$x), length.out = 151)
    data.frame(x = grid, fit = predict(model, data.frame(x = grid)),
               lower = approx(.x$x, .x$mean - .x$sem, grid)$y,
               upper = approx(.x$x, .x$mean + .x$sem, grid)$y)
  }) %>% ungroup()
  ggplot(sm, aes(x, mean, color = series, group = series)) +
    geom_ribbon(data = fitted, aes(x, ymin = lower, ymax = upper, fill = series),
                inherit.aes = FALSE, alpha = .15) +
    geom_line(aes(linetype = series), linewidth = .3, alpha = .6) +
    geom_line(data = fitted, aes(y = fit, linetype = series), linewidth = .85) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .14, linewidth = .4) +
    geom_point(aes(shape = series), size = 2.8, stroke = .7) +
    scale_color_manual(values = CS, labels = c(INH = "INH", inhCTRL = "CTRL (INH)",
                                              MIM = "MIM", mimCTRL = "CTRL (MIM)", CTRL = "CTRL")) +
    scale_fill_manual(values = CS, guide = "none") +
    scale_linetype_manual(values = c(INH = 1, MIM = 1, inhCTRL = 2, mimCTRL = 2, CTRL = 2), guide = "none") +
    scale_shape_manual(values = c(INH = 16, MIM = 16, inhCTRL = 1, mimCTRL = 1, CTRL = 1), guide = "none") +
    scale_x_continuous(breaks = xpos, labels = paste0("T", seq_along(xpos))) +
    labs(x = NULL, y = ylab, title = title) + theme_figure()
}
