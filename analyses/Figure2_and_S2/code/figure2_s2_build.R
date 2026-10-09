# Assemble the approved Supplementary S2 modules in isolated environments.
# This script also saves directly loadable plot/data objects. No terminal
# arguments are used; the working directory is the portable project folder.
if (!exists("figure2_s2_root")) figure2_s2_root <- "."
if (!exists("s2_export")) s2_export <- FALSE
figure2_s2_root <- normalizePath(figure2_s2_root, mustWork = TRUE)
s2_destination_env <- environment()
s2_module_envs <- lapply(c("abc", "defg", "hi"), function(module) {
  e <- new.env(parent = globalenv())
  e$figure2_s2_root <- figure2_s2_root
  e$s2_export <- s2_export
  source(file.path(figure2_s2_root, "code", paste0("figure2_s2_", module, ".R")), local = e)
  check_name <- paste0("s2_", module, "_checks")
  checks <- get(check_name, envir = e)
  if (is.logical(checks$status)) checks$status <- ifelse(checks$status, "PASS", "FAIL")
  assign(check_name, checks, envir = e)
  e
})
names(s2_module_envs) <- c("abc", "defg", "hi")
for (module in names(s2_module_envs)) {
  e <- s2_module_envs[[module]]
  wanted <- c(ls(e, pattern = "^p_S2"), paste0("s2_", module, "_data"),
              paste0("s2_", module, "_checks"))
  for (name in wanted) assign(name, get(name, envir = e), envir = s2_destination_env)
}
s2_plot_names <- paste0("p_S2", letters[1:9])
stopifnot(all(vapply(s2_plot_names, exists, logical(1), envir = s2_destination_env)))
s2_plot_objects <- mget(c(s2_plot_names, "p_S2i_clades", "p_S2i_species"), envir = s2_destination_env)
s2_checks <- dplyr::bind_rows(
  abc = s2_abc_checks, defg = s2_defg_checks, hi = s2_hi_checks, .id = "module")
stopifnot(all(s2_checks$status == "PASS"))
s2_data <- list(abc = s2_abc_data, defg = s2_defg_data, hi = s2_hi_data)
s2_out <- file.path(figure2_s2_root, "figures/supplementary_s2")
dir.create(s2_out, recursive = TRUE, showWarnings = FALSE)
write.table(s2_checks, file.path(s2_out, "S2_reproducibility_checks.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
if (isTRUE(s2_export)) {
  # Include main plots when called by the combined notebook/RStudio runner.
  if (!exists("figure2_plot_objects")) figure2_plot_objects <- list()
  combined_plot_objects <- list(Figure2 = figure2_plot_objects, FigureS2 = s2_plot_objects)
  saveRDS(combined_plot_objects, file.path(s2_out, "Figure2_and_S2_plot_objects.rds"), compress = "gzip")
  saveRDS(s2_data, file.path(s2_out, "S2_analytical_data.rds"), compress = "gzip")
  workspace_names <- c(s2_plot_names, "p_S2i_clades", "p_S2i_species", "s2_plot_objects",
                       "figure2_plot_objects", "s2_data", "s2_checks")
  save(list = workspace_names, envir = s2_destination_env,
       file = file.path(s2_out, "Figure2_and_S2_workspace.RData"), compress = "gzip")
}
message("S2 a-i reconstructed; ", nrow(s2_checks), " checks passed.")
