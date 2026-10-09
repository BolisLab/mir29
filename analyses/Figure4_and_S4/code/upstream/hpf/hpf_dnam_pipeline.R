# Reconstruct the historical HPF DNAm analytical checkpoints for Figure 4 r/s.
# The calculation functions never read historical result files. Those are used
# only by validate_hpf_dnam.R after regeneration.

hpf_require_packages <- function() {
  for (package in c("limma", "tibble", "dplyr", "tidyr")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  }
}

hpf_time_info <- function() {
  # Retained verbatim for historical RDS compatibility. T1-T6 hours below are
  # historical nominal labels, not elapsed calendar time; display T1-T8 in figures.
  tibble::tibble(
    timepoint = 1:8, timepoint_name = paste0("T", 1:8),
    display_label = c(paste0("T", 1:5, "\n", (1:5)*96, " h"),
                      "T6\n576 h\nend treatment", "T7\n+168 h", "T8\n+240 h"),
    phase = c(rep("ACTIVE_TREATMENT", 6), rep("WITHDRAWAL", 2)))
}

hpf_pair_samples <- function(samples) {
  stopifnot(is.data.frame(samples), !anyDuplicated(rownames(samples)),
            all(c("GROUP", "TP2") %in% names(samples)))
  tp <- as.integer(as.character(samples$TP2))
  group <- as.character(samples$GROUP)
  if (anyNA(tp) || nrow(samples) != 16L) stop("Expected the 16 historical HPF samples.")
  sample_id <- rownames(samples)
  select_one <- function(timepoint, condition) {
    ids <- sample_id[tp == timepoint & group == condition]
    if (length(ids) != 1L) stop("Expected one ", condition, " sample at T", timepoint)
    ids
  }
  ti <- hpf_time_info()
  tibble::tibble(timepoint = ti$timepoint, timepoint_name = ti$timepoint_name,
                ctrl_sample = vapply(1:8, select_one, character(1), condition="CTRL"),
                combo_sample = vapply(1:8, select_one, character(1), condition="COMBO"),
                display_label = ti$display_label, phase = ti$phase)
}

hpf_paired_delta <- function(matrix, pairs) {
  stopifnot(all(c(pairs$ctrl_sample, pairs$combo_sample) %in% colnames(matrix)))
  result <- matrix[, pairs$combo_sample, drop=FALSE] - matrix[, pairs$ctrl_sample, drop=FALSE]
  colnames(result) <- paste0("T", 1:8)
  result
}

hpf_collapse_mean <- function(matrix) {
  # Arithmetic mean across EPIC suffix probes, separately per sample, with
  # missing values omitted. An all-missing CpG/sample remains NA.
  prefix <- sub("_.*$", "", rownames(matrix))
  numerator <- rowsum(matrix, group=prefix, reorder=TRUE, na.rm=TRUE)
  denominator <- rowsum(1L * !is.na(matrix), group=prefix, reorder=TRUE)
  result <- numerator / denominator
  result[denominator == 0] <- NA_real_
  result
}

hpf_matrix_check <- function(actual, expected, label, tolerance=1e-12) {
  stopifnot(identical(dimnames(actual), dimnames(expected)))
  error <- 0
  for (j in seq_len(ncol(actual))) {
    if (!identical(is.na(actual[,j]), is.na(expected[,j]))) stop(label, ": missing-value mismatch")
    delta <- abs(actual[,j] - expected[,j])
    error <- max(error, delta, na.rm=TRUE)
  }
  if (!is.finite(error) || error > tolerance) stop(label, ": maximum error ", error)
  error
}

hpf_load_inputs <- function(input_dir) {
  samples <- readRDS(file.path(input_dir, "hpf_samples_all.RDS"))
  pairs <- hpf_pair_samples(samples)
  bfile <- file.path(input_dir, "hpf_beta_all.collapsed.RDS")
  mfile <- file.path(input_dir, "hpf_mval_all.collapsed.RDS")
  raw <- readRDS(file.path(input_dir, "hpf_beta_all.RDS"))
  stopifnot(is.matrix(raw), is.numeric(raw), !anyDuplicated(rownames(raw)),
            !anyDuplicated(colnames(raw)), all(is.na(raw) | (raw>0 & raw<1)))
  # Use the supplied original intermediates; raw beta is used to audit their
  # provenance. If they are absent, both can be reconstructed from raw beta.
  b <- if (file.exists(bfile)) readRDS(bfile) else hpf_collapse_mean(raw)
  m <- if (file.exists(mfile)) readRDS(mfile) else hpf_collapse_mean(log2(raw/(1-raw)))
  stopifnot(identical(dimnames(b), dimnames(m)), !anyDuplicated(rownames(b)))
  beta_error <- hpf_matrix_check(hpf_collapse_mean(raw), b, "Collapsed beta")
  m_error <- hpf_matrix_check(hpf_collapse_mean(log2(raw/(1-raw))), m, "Collapsed M")
  summary <- tibble::tibble(n_original_rows=as.integer(nrow(raw)),
    n_collapsed_CpGs=as.integer(nrow(b)), n_removed_rows=as.integer(nrow(raw)-nrow(b)),
    collapse_required=nrow(raw)!=nrow(b))
  list(beta=b, historical_M=m, pairs=pairs, collapse_summary=summary,
       import_check=data.frame(matrix=c("beta_collapsed", "M_collapsed"),
                               maximum_absolute_error=c(beta_error,m_error)))
}

hpf_probe_table <- function(beta, historical_M, pairs) {
  delta <- hpf_paired_delta(historical_M, pairs)
  dbeta <- hpf_paired_delta(beta, pairs)
  keep <- rowSums(!is.finite(delta))==0 & rowSums(!is.finite(dbeta))==0
  delta <- delta[keep,,drop=FALSE]
  dbeta <- dbeta[keep,,drop=FALSE]
  design <- cbind(T6_offset=1, ActiveSlope=c(-5:0,0,0),
                  PostSlope=c(rep(0,6),168/96,240/96))
  rownames(design) <- colnames(delta)
  fit <- limma::lmFit(delta, design)
  contrasts <- cbind(ActiveSlope=c(0,1,0), T6_offset=c(1,0,0),
    PostSlope=c(0,0,1), T7_offset=c(1,0,168/96), T8_offset=c(1,0,240/96),
    T8_minus_T6=c(0,0,240/96))
  rownames(contrasts) <- colnames(design)
  tested <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend=FALSE, robust=FALSE)
  q <- apply(tested$p.value, 2, p.adjust, method="BH")
  c <- tested$coefficients
  p <- tested$p.value
  out <- tibble::tibble(probe=rownames(delta),
    active_slope_M=c[,"ActiveSlope"], active_slope_pvalue=p[,"ActiveSlope"],
    active_slope_FDR=q[,"ActiveSlope"], T6_offset_M=c[,"T6_offset"],
    T6_offset_pvalue=p[,"T6_offset"], T6_offset_FDR=q[,"T6_offset"],
    post_slope_M_per_96h=c[,"PostSlope"], post_slope_pvalue=p[,"PostSlope"],
    post_slope_FDR=q[,"PostSlope"], T7_offset_M=c[,"T7_offset"],
    T7_offset_FDR=q[,"T7_offset"], T8_offset_M=c[,"T8_offset"],
    T8_offset_FDR=q[,"T8_offset"], T8_minus_T6_M=c[,"T8_minus_T6"],
    T8_minus_T6_FDR=q[,"T8_minus_T6"], raw_T6_delta_M=delta[,6],
    raw_T7_delta_M=delta[,7], raw_T8_delta_M=delta[,8],
    raw_T6_delta_beta=dbeta[,6], raw_T7_delta_beta=dbeta[,7], raw_T8_delta_beta=dbeta[,8])
  out$T6_significant <- out$T6_offset_FDR < 0.05
  out$active_significant <- out$active_slope_FDR < 0.05
  out$post_significant <- out$post_slope_FDR < 0.05
  out$same_sign_T6_T7 <- sign(out$T6_offset_M)==sign(out$T7_offset_M)
  out$same_sign_T6_T8 <- sign(out$T6_offset_M)==sign(out$T8_offset_M)
  out$post_change_same_direction_as_T6 <- sign(out$T8_minus_T6_M)==sign(out$T6_offset_M)
  out$post_change_opposite_to_T6 <- sign(out$T8_minus_T6_M)!=sign(out$T6_offset_M)
  out$retention_ratio_T7 <- ifelse(abs(out$T6_offset_M)>0.05, out$T7_offset_M/out$T6_offset_M, NA_real_)
  out$retention_ratio_T8 <- ifelse(abs(out$T6_offset_M)>0.05, out$T8_offset_M/out$T6_offset_M, NA_real_)
  # Only these three classes occur in the historical data. Do not infer an
  # undocumented taxonomy for new datasets with significant post slopes.
  if (any(out$post_significant)) stop("Post-slope significant CpGs: historical classification no longer applies.")
  out$withdrawal_class <- ifelse(!out$T6_significant, "NO_SIGNIFICANT_T6_OFFSET",
    ifelse(out$same_sign_T6_T8, "MAINTAINED_NO_DETECTABLE_POST_CHANGE", "UNRESOLVED"))
  # The historical tibble stores probe IDs in their own column, not as names
  # attached to individual numeric/logical columns.
  out[] <- lapply(out, unname)
  list(table=out, design=design, contrasts=contrasts)
}

hpf_projection_curve <- function(delta, ids, ref, weights, method, universe) {
  idx <- match(ids, rownames(delta))
  stopifnot(!anyNA(idx))
  r <- ref[idx]; w <- weights[idx]
  stopifnot(all(is.finite(r)), all(is.finite(w)), all(w>0))
  projection <- as.numeric(crossprod(w*r, delta[idx,,drop=FALSE])/sum(w*r^2))
  if (!is.finite(projection[6]) || projection[6]==0) stop("Invalid T6 projection")
  ti <- hpf_time_info()
  tibble::tibble(method=method, universe=universe, timepoint_name=ti$timepoint_name,
    timepoint=ti$timepoint, raw_projection=projection, retention=projection/projection[6],
    n_features=as.integer(length(ids)), display_label=ti$display_label, phase=ti$phase)
}

hpf_decision_trials <- function(beta, pairs, probe_table, collapse_summary) {
  # This later historical stage converts ALREADY COLLAPSED BETA to M.
  # The exact core remains the core selected by the earlier T1-T8 piecewise fit.
  m <- log2(beta/(1-beta))
  delta <- hpf_paired_delta(m, pairs)
  delta <- delta[rowSums(!is.finite(delta))==0,,drop=FALSE]
  core_original <- probe_table$probe[probe_table$T6_significant]
  core <- intersect(core_original, rownames(delta))
  design <- cbind(T6_effect=1, ActiveSlope=-5:0)
  rownames(design) <- paste0("T",1:6)
  fit <- limma::eBayes(limma::lmFit(delta[,1:6,drop=FALSE],design), trend=TRUE, robust=FALSE)
  ranking <- rownames(delta)[order(abs(fit$t[,"T6_effect"]), decreasing=TRUE)]
  remainder <- setdiff(ranking,core)
  universes <- list(PANEL_P_CORE_EXACT=core)
  for (size in c(5000L,10000L,25000L,50000L)) {
    stopifnot(length(core)<size, nrow(delta)>=size)
    universes[[paste0("CORE_PLUS_TO_",size)]] <- c(core,head(remainder,size-length(core)))
  }
  methods <- c("RAW_T6_UNWEIGHTED", "MODEL_T6_UNWEIGHTED", "MODEL_T6_PRECISION_WEIGHTED")
  refs <- list(delta[,6],fit$coefficients[,"T6_effect"],fit$coefficients[,"T6_effect"])
  precision <- 1/(fit$s2.post * fit$stdev.unscaled[,"T6_effect"]^2)
  # Cap weights across the complete universe BEFORE subsetting. Type 8
  # reproduces the cap recorded in the historical weighted checkpoint.
  weight_cap <- as.numeric(quantile(precision, 0.995, type=8))
  weights <- list(rep(1,nrow(delta)),rep(1,nrow(delta)),
    pmin(precision, weight_cap))
  nested <- dplyr::bind_rows(lapply(seq_along(methods), function(k) {
    dplyr::bind_rows(lapply(names(universes),function(nm) {
      hpf_projection_curve(delta,universes[[nm]],refs[[k]],weights[[k]],methods[k],nm)
    }))
  }))
  comparators <- dplyr::bind_rows(lapply(seq_along(methods),function(k) {
    hpf_projection_curve(delta,rownames(delta),refs[[k]],weights[[k]],methods[k],"ALL_COMPLETE")
  }))
  all_curves <- dplyr::bind_rows(nested,comparators)
  envelope <- dplyr::summarise(
    dplyr::group_by(nested,method,timepoint,timepoint_name,display_label,phase),
    median_retention=median(retention), minimum_retention=min(retention),
    maximum_retention=max(retention), lower_quartile=as.numeric(quantile(retention,.25)),
    upper_quartile=as.numeric(quantile(retention,.75)), n_universes=dplyr::n(), .groups="drop")
  endpoints <- tidyr::pivot_wider(
    dplyr::select(dplyr::filter(all_curves,timepoint %in% c(7L,8L)),
                  method,universe,n_features,timepoint_name,retention),
    names_from=timepoint_name,values_from=retention)
  endpoints <- dplyr::arrange(endpoints,method,dplyr::desc(n_features))
  envelope_endpoints <- dplyr::select(dplyr::filter(envelope,timepoint %in% c(7L,8L)),
    method,timepoint_name,median_retention,minimum_retention,maximum_retention,
    lower_quartile,upper_quartile)
  core_source <- "HPF_DNAm_WITHDRAWAL_ANALYSIS/HPF_probe_retention_table_full.rds: T6_significant"
  out <- list(
    source=list(preprocessing="hpf_beta_all -> collapse beta by prefix -> convert collapsed beta to M-value",
      panel_p_core_source=core_source,
      exact_core_definition="historical T6_offset_FDR < 0.05 CpGs used for panel p",
      expanded_universe_ranking="|moderated t for modeled T6 effect estimated from T1-T6 only|",
      note=paste("The inherited exact core is selected by the historical T1-T8 piecewise fit.",
                 "Only the expansion ranking, modeled T6 reference and precision weights use T1-T6 alone.")),
    preprocessing=list(collapse_summary=collapse_summary,n_complete_CpGs=as.integer(nrow(delta))),
    time_info=hpf_time_info(), pair_table=pairs,
    core_overlap_summary=tibble::tibble(source=core_source,
      n_original_panel_p_core=length(core_original),n_core_in_corrected_complete_matrix=length(core),
      n_core_missing_after_corrected_preprocessing=length(setdiff(core_original,core)),
      overlap_fraction=length(core)/length(core_original)),
    missing_core_ids=setdiff(core_original,core),
    universe_summary=tibble::tibble(universe=names(universes),
      n_features=vapply(universes,length,integer(1)),
      contains_full_exact_core=vapply(universes,function(x) all(core %in% x),logical(1))),
    universe_ids=universes,
    curves=list(nested_responsive=nested,all_complete_comparators=comparators,
                all=all_curves,responsive_envelope=envelope),
    endpoint_summary=endpoints,envelope_endpoint_summary=envelope_endpoints)
  list(object=out,design=design,weight_cap=weight_cap)
}

hpf_final_export <- function(decision, probe_table, probe_path) {
  ids <- decision$universe_ids$PANEL_P_CORE_EXACT
  nested <- dplyr::arrange(dplyr::filter(decision$curves$nested_responsive,
                                       method=="RAW_T6_UNWEIGHTED"),universe,timepoint)
  scatter <- as.data.frame(probe_table[probe_table$probe %in% ids,
                                      c("probe","raw_T6_delta_M","raw_T8_delta_M")])
  names(scatter) <- c("feature_id","effect_T6","effect_T8")
  rownames(scatter) <- NULL
  stats <- tibble::tibble(
    retention=sum(scatter$effect_T6*scatter$effect_T8)/sum(scatter$effect_T6^2),
    pearson_r=cor(scatter$effect_T6,scatter$effect_T8),
    spearman_r=cor(scatter$effect_T6,scatter$effect_T8,method="spearman"),n=nrow(scatter))
  list(source=list(decision_object="HPF_DNAm_retention_decision_trials from workspace",
    probe_table=normalizePath(probe_path,mustWork=FALSE),core_definition="PANEL_P_CORE_EXACT",
    curve_method="RAW_T6_UNWEIGHTED",effect_definition="raw paired COMBO minus CTRL M-value"),
    exact_core_ids=ids,missing_core_ids=setdiff(ids,scatter$feature_id),nested_curves=nested,
    exact_core_curve=dplyr::filter(nested,universe=="PANEL_P_CORE_EXACT"),
    envelope_curve=dplyr::filter(decision$curves$responsive_envelope,method=="RAW_T6_UNWEIGHTED"),
    scatter_data=scatter,scatter_stats=stats,universe_summary=decision$universe_summary)
}

run_hpf_dnam <- function(input_dir, output_dir) {
  hpf_require_packages()
  input_dir <- normalizePath(input_dir,mustWork=TRUE)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  filenames <- c("HPF_probe_retention_table_full.rds", "HPF_DNAm_retention_decision_trials.rds",
                 "HPF_DNAm_FINAL_CORE_TRIALS.rds")
  if (any(file.exists(file.path(output_dir,filenames)))) stop("Choose a fresh output directory; existing results are preserved.")
  message("Loading and validating the original imported HPF matrices.")
  input <- hpf_load_inputs(input_dir)
  message("Rebuilding the historical T1-T8 piecewise probe table.")
  probes <- hpf_probe_table(input$beta,input$historical_M,input$pairs)
  message("Rebuilding T1-T6 ranked universes and all retention methods.")
  decision <- hpf_decision_trials(input$beta,input$pairs,probes$table,input$collapse_summary)
  final <- hpf_final_export(decision$object,probes$table,file.path(output_dir,filenames[1]))
  objects <- list(probes$table,decision$object,final)
  for (i in seq_along(objects)) saveRDS(objects[[i]],file.path(output_dir,filenames[i]),compress="gzip")
  input_names <- c("hpf_beta_all.RDS","hpf_beta_all.collapsed.RDS",
                   "hpf_mval_all.collapsed.RDS","hpf_samples_all.RDS")
  input_paths <- file.path(input_dir,input_names)
  input_paths <- input_paths[file.exists(input_paths)]
  manifest <- data.frame(path=input_paths,bytes=file.info(input_paths)$size,
                         md5=unname(tools::md5sum(input_paths)))
  write.table(manifest,file.path(output_dir,"input_manifest.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(input$import_check,file.path(output_dir,"import_verification.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  saveRDS(list(piecewise_design=probes$design,piecewise_contrasts=probes$contrasts,
    piecewise_eBayes=list(trend=FALSE,robust=FALSE),
    expansion_design=decision$design,expansion_eBayes=list(trend=TRUE,robust=FALSE),
    precision_weight_cap=list(quantile=0.995,quantile_type=8,
                              universe="ALL_COMPLETE",value=decision$weight_cap),
    historical_core_uses_T7_T8=TRUE,expansion_ranking_uses_T7_T8=FALSE,
    historical_M="mean(log2(beta / (1-beta))) across suffix probes",
    decision_M="log2(mean(beta) / (1-mean(beta))) across suffix probes",
    ratio_denominator_minimum_absolute_M=0.05),file.path(output_dir,"analysis_specification.rds"))
  capture.output(sessionInfo(),file=file.path(output_dir,"sessionInfo.txt"))
  message("Written: ",paste(filenames,collapse=", "))
  print(final$scatter_stats)
  invisible(setNames(objects,filenames))
}
