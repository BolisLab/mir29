tcga_sample_id_mrc <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    paste(z[1:4], collapse = "-")
  }, character(1))
}

tcga_patient_id_mrc <- function(x) {
  substr(tcga_sample_id_mrc(x), 1, 12)
}

tcga_type_code_mrc <- function(x) {
  sid <- tcga_sample_id_mrc(x)
  vapply(strsplit(sid, "-", fixed = TRUE), function(z) {
    if (length(z) < 4) {
      return(NA_character_)
    }
    substr(z[4], 1, 2)
  }, character(1))
}

get_age_table_mrc <- function(cancer) {
  P <- allPred450K[[cancer]]
  if (is.null(P)) {
    stop(cancer, ": allPred450K missing.")
  }
  P <- as.data.frame(P)
  age_idx <- which(colnames(P) == "AGE")
  if (length(age_idx) == 0) {
    stop(cancer, ": AGE column not found.")
  }
  age_idx <- age_idx[1]
  D <- data.frame(SampleID = tcga_sample_id_mrc(rownames(P)), PatientID = tcga_patient_id_mrc(rownames(P)),
    TypeCode = tcga_type_code_mrc(rownames(P)), AGE = suppressWarnings(as.numeric(P[, age_idx])),
    stringsAsFactors = FALSE)
  D <- D[D$TypeCode == "01" & is.finite(D$AGE) & !is.na(D$PatientID), , drop = FALSE]
  D <- D[order(D$PatientID, D$SampleID), , drop = FALSE]
  D <- D[!duplicated(D$PatientID), , drop = FALSE]
  D[, c("PatientID", "AGE"), drop = FALSE]
}

get_mrc_table <- function(cancer) {
  PC <- FIG6E_FINAL$perCancer[[cancer]]
  if (is.null(PC)) {
    stop(cancer, ": absent from FIG6E_FINAL$perCancer.")
  }
  X <- if (!is.null(PC$mrc)) {
    PC$mrc
  } else {
    PC$erc
  }
  if (is.null(X)) {
    stop(cancer, ": MRC/ERC table missing.")
  }
  X <- as.data.frame(X)
  value_candidates <- c("MRC_Z", "MRC_score", "MRC_Score", "MRC", "ERC_Z", "ERC_score", "ERC_Score",
    "ERC_raw")
  value_col <- value_candidates[value_candidates %in% names(X)][1]
  if (length(value_col) == 0 || is.na(value_col)) {
    stop(cancer, ": MRC/ERC score column not identified.")
  }
  if ("Patient" %in% names(X)) {
    PatientID <- as.character(X$Patient)
  } else if ("PatientID" %in% names(X)) {
    PatientID <- as.character(X$PatientID)
  } else if ("Sample" %in% names(X)) {
    PatientID <- substr(as.character(X$Sample), 1, 12)
  } else if ("SampleID" %in% names(X)) {
    PatientID <- substr(as.character(X$SampleID), 1, 12)
  } else {
    stop(cancer, ": no patient/sample identifier in MRC table.")
  }
  if ("Type" %in% names(X)) {
    TypeCode <- as.character(X$Type)
  } else if ("TypeCode" %in% names(X)) {
    TypeCode <- as.character(X$TypeCode)
  } else if ("Sample" %in% names(X)) {
    TypeCode <- tcga_type_code_mrc(X$Sample)
  } else if ("SampleID" %in% names(X)) {
    TypeCode <- tcga_type_code_mrc(X$SampleID)
  } else {
    stop(cancer, ": cannot identify tumor rows in MRC table.")
  }
  OUT <- data.frame(PatientID = PatientID, TypeCode = TypeCode, MRC = suppressWarnings(as.numeric(X[[value_col]])),
    stringsAsFactors = FALSE)
  OUT <- OUT[OUT$TypeCode == "01" & is.finite(OUT$MRC) & !is.na(OUT$PatientID), , drop = FALSE]
  OUT <- OUT[order(OUT$PatientID), , drop = FALSE]
  OUT <- OUT[!duplicated(OUT$PatientID), , drop = FALSE]
  OUT[, c("PatientID", "MRC"), drop = FALSE]
}

get_fidelity_table_mrc <- function(cancer) {
  Fobj <- TCGA_fidelity_main15[[cancer]]
  if (is.null(Fobj) || is.null(Fobj$status) || Fobj$status != "OK") {
    stop(cancer, ": invalid fidelity object.")
  }
  S <- as.data.frame(Fobj$scores)
  S <- S[S$SampleType == "PrimaryTumor", , drop = FALSE]
  if (!"Stable_FidelityLoss_RMS_M" %in% names(S)) {
    stop(cancer, ": Stable_FidelityLoss_RMS_M missing.")
  }
  if ("PatientID" %in% names(S)) {
    PatientID <- as.character(S$PatientID)
  } else {
    PatientID <- substr(as.character(S$SampleID), 1, 12)
  }
  OUT <- data.frame(PatientID = PatientID, FidelityLoss = suppressWarnings(as.numeric(S$Stable_FidelityLoss_RMS_M)),
    stringsAsFactors = FALSE)
  OUT <- OUT[is.finite(OUT$FidelityLoss) & !is.na(OUT$PatientID), , drop = FALSE]
  OUT <- OUT[!duplicated(OUT$PatientID), , drop = FALSE]
  OUT
}

build_mrc_fidelity_data <- function(cancer) {
  M <- get_mrc_table(cancer)
  F <- get_fidelity_table_mrc(cancer)
  A <- get_age_table_mrc(cancer)
  D <- merge(M, F, by = "PatientID")
  D <- merge(D, A, by = "PatientID")
  D <- D[complete.cases(D[, c("MRC", "FidelityLoss", "AGE")]), , drop = FALSE]
  rownames(D) <- NULL
  D
}

run_mrc_fidelity_cancer <- function(cancer) {
  D <- build_mrc_fidelity_data(cancer)
  N <- nrow(D)
  if (N <= 6) {
    stop(cancer, ": insufficient sample size.")
  }
  fit_mrc <- lm(MRC ~ ns(AGE, df = 3), data = D)
  fit_fid <- lm(FidelityLoss ~ ns(AGE, df = 3), data = D)
  mrc_res <- residuals(fit_mrc)
  fid_res <- residuals(fit_fid)
  ct <- suppressWarnings(cor.test(mrc_res, fid_res, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  rho_clip <- min(max(rho, -0.999999), 0.999999)
  fisher_z <- atanh(rho_clip)
  fisher_se <- 1/sqrt(N - 6)
  CI_low <- tanh(fisher_z - qnorm(0.975) * fisher_se)
  CI_high <- tanh(fisher_z + qnorm(0.975) * fisher_se)
  data.frame(Cancer = cancer, N = N, rho = rho, CI_low = CI_low, CI_high = CI_high, pvalue = ct$p.value,
    Fisher_z = fisher_z, Fisher_SE = fisher_se, stringsAsFactors = FALSE)
}

