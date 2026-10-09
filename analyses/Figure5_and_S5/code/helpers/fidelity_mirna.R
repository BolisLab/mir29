tcga_sample_id_local <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    paste(z[1:4], collapse = "-")
  }, character(1))
}

tcga_patient_id_local <- function(x) {
  substr(tcga_sample_id_local(x), 1, 12)
}

tcga_type_code_local <- function(x) {
  sid <- tcga_sample_id_local(x)
  vapply(strsplit(sid, "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    substr(z[4], 1, 2)
  }, character(1))
}

prepare_exact_clock_overlap <- function(cancer) {
  P <- allPred450K[[cancer]]
  if (is.null(P)) {
    stop("No allPred450K object for ", cancer)
  }
  P <- as.data.frame(P)
  age_idx <- which(colnames(P) == "AGE")
  hor_idx <- which(colnames(P) == "horvath")
  if (length(age_idx) == 0) {
    stop(cancer, ": AGE column not found in allPred450K")
  }
  if (length(hor_idx) == 0) {
    stop(cancer, ": horvath column not found in allPred450K")
  }
  age_idx <- age_idx[1]
  hor_idx <- hor_idx[1]
  A <- data.frame(MethID = rownames(P), SampleID = tcga_sample_id_local(rownames(P)), PatientID = tcga_patient_id_local(rownames(P)),
    TypeCode = tcga_type_code_local(rownames(P)), AGE = suppressWarnings(as.numeric(P[, age_idx])),
    Horvath = suppressWarnings(as.numeric(P[, hor_idx])), stringsAsFactors = FALSE)
  A <- A[A$TypeCode == "01" & is.finite(A$AGE) & is.finite(A$Horvath) & !is.na(A$SampleID), , drop = FALSE]
  A <- A[!duplicated(A$SampleID), , drop = FALSE]
  A <- A[order(A$PatientID, A$SampleID), , drop = FALSE]
  A <- A[!duplicated(A$PatientID), , drop = FALSE]
  A
}

build_analysis_data <- function(cancer) {
  Fobj <- TCGA_fidelity_main15[[cancer]]
  Mobj <- TCGA_miR29[[cancer]]
  if (is.null(Fobj) || is.null(Fobj$status) || Fobj$status != "OK") {
    stop("Invalid fidelity object for ", cancer)
  }
  if (is.null(Mobj) || is.null(Mobj$status) || Mobj$status != "OK") {
    stop("Invalid miRNA object for ", cancer)
  }
  F <- Fobj$scores[Fobj$scores$SampleType == "PrimaryTumor", c("SampleID", "Stable_FidelityLoss_RMS_M"),
    drop = FALSE]
  names(F)[names(F) == "Stable_FidelityLoss_RMS_M"] <- "FidelityLoss"
  M <- Mobj$table[, c("SampleID", "pan_miR29"), drop = FALSE]
  A <- prepare_exact_clock_overlap(cancer)
  A <- A[, c("SampleID", "PatientID", "AGE", "Horvath"), drop = FALSE]
  D <- merge(F, M, by = "SampleID")
  D <- merge(D, A, by = "SampleID")
  D <- D[complete.cases(D[, c("pan_miR29", "FidelityLoss", "AGE", "Horvath")]), , drop = FALSE]
  D <- D[order(D$PatientID, D$SampleID), , drop = FALSE]
  D <- D[!duplicated(D$PatientID), , drop = FALSE]
  rownames(D) <- NULL
  D
}

run_spline_residualized_spearman <- function(cancer) {
  D <- build_analysis_data(cancer)
  N <- nrow(D)
  if (N <= 6) {
    stop(cancer, ": insufficient N for nonlinear residualized correlation.")
  }
  fit_pan <- lm(pan_miR29 ~ ns(AGE, df = 3), data = D)
  fit_fid <- lm(FidelityLoss ~ ns(AGE, df = 3), data = D)
  pan_resid <- residuals(fit_pan)
  fid_resid <- residuals(fit_fid)
  ct <- suppressWarnings(cor.test(pan_resid, fid_resid, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  rho_clip <- min(max(rho, -0.999999), 0.999999)
  z <- atanh(rho_clip)
  se_z <- 1/sqrt(N - 6)
  ci_low <- tanh(z - qnorm(0.975) * se_z)
  ci_high <- tanh(z + qnorm(0.975) * se_z)
  data.frame(Cancer = cancer, N = N, rho = rho, CI_low = ci_low, CI_high = ci_high, pvalue = ct$p.value,
    Fisher_z = z, Fisher_SE = se_z, stringsAsFactors = FALSE)
}

