# Optional matrix-level reference derivation; not run during document rendering.  BN must contain
# prefiltered, imputed autosomal beta values from normal samples.
MAX_NORMAL_SD_BETA <- 0.05
MAX_ABS_AGE_COR <- 0.25
MAX_AGE_DBETA_DECADE <- 0.01
MIN_STABLE_CPGS <- 5000
ESCAPE_BETA_FLOOR <- 0.05
DBETA_DIRECTIONAL <- 0.1
BLOCK_SIZE <- 50
beta_to_M <- function(x, eps = 1e-06) {
  x <- pmin(pmax(x, eps), 1 - eps)
  log2(x/(1 - x))
}

derive_stable_reference <- function(BN, normal_age) {
  if (ncol(BN) != length(normal_age)) {
    stop("Normal methylation matrix and age vector have different lengths.")
  }

  normal_sd_beta <- matrixStats::rowSds(BN)
  is_stable <- normal_sd_beta <= MAX_NORMAL_SD_BETA

  age_centered <- normal_age - mean(normal_age)
  age_ss <- sum(age_centered^2)
  if (!is.finite(age_ss) || age_ss <= 0)
    stop("Normal age has no usable variance.")

  age_num <- as.vector(BN %*% age_centered)
  age_delta_beta_decade <- 10 * age_num/age_ss

  denom_cor <- sqrt(((ncol(BN) - 1) * normal_sd_beta^2) * age_ss)

  age_cor <- age_num/denom_cor
  age_cor[!is.finite(age_cor)] <- 0

  is_age_independent <- abs(age_cor) <= MAX_ABS_AGE_COR & abs(age_delta_beta_decade) <= MAX_AGE_DBETA_DECADE

  stable <- is_stable & is_age_independent
  stable_cpgs <- rownames(BN)[stable]

  if (length(stable_cpgs) < MIN_STABLE_CPGS) {
    stop("Too few stable age-independent CpGs: ", length(stable_cpgs))
  }

  BNs <- BN[stable_cpgs, , drop = FALSE]

  ref_beta <- matrixStats::rowMedians(BNs)
  ref_beta_mad <- matrixStats::rowMads(BNs, center = ref_beta, constant = 1.4826)
  ref_M <- matrixStats::rowMedians(beta_to_M(BNs))

  names(ref_beta) <- stable_cpgs
  names(ref_beta_mad) <- stable_cpgs
  names(ref_M) <- stable_cpgs

  list(stable_CpGs = stable_cpgs, ref_beta = ref_beta, ref_beta_mad = ref_beta_mad, ref_M = ref_M,
    normal_sd_beta = normal_sd_beta[stable], age_cor = age_cor[stable], age_delta_beta_decade = age_delta_beta_decade[stable])
}


score_against_reference <- function(B, ann, reference, block_size = BLOCK_SIZE) {
  cpgs <- intersect(reference$stable_CpGs, rownames(B))
  if (length(cpgs) < MIN_STABLE_CPGS) {
    stop("Too few reference CpGs available in target matrix: ", length(cpgs))
  }

  B <- B[cpgs, , drop = FALSE]
  ref_beta <- reference$ref_beta[cpgs]
  ref_beta_mad <- reference$ref_beta_mad[cpgs]
  ref_M <- reference$ref_M[cpgs]
  escape_width <- pmax(3 * ref_beta_mad, ESCAPE_BETA_FLOOR)

  n <- ncol(B)
  RMS_M <- MeanAbsM <- MeanAbsBeta <- MeanSignedBeta <- EscapeFraction <- Hyper10 <- Hypo10 <- rep(NA_real_,
    n)

  starts <- seq(1, n, by = block_size)

  for (s in starts) {
    jj <- s:min(s + block_size - 1, n)
    BB <- B[, jj, drop = FALSE]

    dBeta <- sweep(BB, 1, ref_beta, "-")
    dM <- sweep(beta_to_M(BB), 1, ref_M, "-")

    RMS_M[jj] <- sqrt(colMeans(dM^2))
    MeanAbsM[jj] <- colMeans(abs(dM))
    MeanAbsBeta[jj] <- colMeans(abs(dBeta))
    MeanSignedBeta[jj] <- colMeans(dBeta)
    EscapeFraction[jj] <- colMeans(abs(dBeta) > escape_width)
    Hyper10[jj] <- colMeans(dBeta >= DBETA_DIRECTIONAL)
    Hypo10[jj] <- colMeans(dBeta <= -DBETA_DIRECTIONAL)

    rm(BB, dBeta, dM)
  }

  score_df <- data.frame(Stable_FidelityLoss_RMS_M = RMS_M, Stable_EscapeFraction = EscapeFraction,
    Stable_MeanAbsDeltaM = MeanAbsM, Stable_MeanAbsDeltaBeta = MeanAbsBeta, Stable_MeanSignedDeltaBeta = MeanSignedBeta,
    Stable_HyperFraction10 = Hyper10, Stable_HypoFraction10 = Hypo10, row.names = colnames(B))

  cbind(ann, score_df[ann$MethID, , drop = FALSE])
}


