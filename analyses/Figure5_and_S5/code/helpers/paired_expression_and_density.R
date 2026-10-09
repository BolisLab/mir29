one_sample_per_patient <- function(D, sample_type) {
  X <- D[D$Type == sample_type, , drop = FALSE]
  if (nrow(X) == 0) {
    return(X)
  }
  X <- X[order(X$Patient, X$Sample), , drop = FALSE]
  X <- X[!duplicated(X$Patient), , drop = FALSE]
  X
}

run_histotype_meta <- function(D, readout) {
  X <- D[D$Readout == readout, , drop = FALSE]
  X <- X[match(MASTER_CANCERS16, X$Cancer), , drop = FALSE]
  m <- metafor::rma(yi = Effect, sei = SE, data = X, method = "REML", test = "knha", slab = Cancer)
  stopifnot(m$k == 16L)
  m
}

build_organ_effects <- function(D, readout) {
  X <- D[D$Readout == readout, , drop = FALSE]
  X$Organ <- unname(ORGAN_MAP[X$Cancer])
  if (anyNA(X$Organ)) {
    stop("Missing organ mapping.")
  }
  SPLIT <- split(X, X$Organ)
  O <- do.call(rbind, lapply(names(SPLIT), function(org) {
    Z <- SPLIT[[org]]
    if (nrow(Z) == 1) {
      return(data.frame(Organ = org, Histotypes = Z$Cancer, N_histotypes = 1L, Effect = Z$Effect,
        SE = Z$SE, CI_low = Z$Effect - 1.96 * Z$SE, CI_high = Z$Effect + 1.96 * Z$SE, tau2_within = 0,
        Q_p_within = NA_real_, stringsAsFactors = FALSE))
    }
    m_org <- metafor::rma(yi = Effect, sei = SE, data = Z, method = "REML", test = "z")
    data.frame(Organ = org, Histotypes = paste(Z$Cancer, collapse = "/"), N_histotypes = nrow(Z),
      Effect = as.numeric(m_org$b[1]), SE = m_org$se, CI_low = m_org$ci.lb, CI_high = m_org$ci.ub,
      tau2_within = m_org$tau2, Q_p_within = m_org$QEp, stringsAsFactors = FALSE)
  }))
  rownames(O) <- NULL
  O
}

run_organ_meta <- function(D, readout) {
  O <- build_organ_effects(D, readout)
  m <- metafor::rma(yi = Effect, sei = SE, data = O, method = "REML", test = "knha", slab = Organ)
  stopifnot(m$k == 13L)
  list(effects = O, meta = m)
}

weighted_kde2d <- function(x, y, weights, n = 200, h_x = NULL, h_y = NULL, expand = 0.1) {
  ok <- is.finite(x) & is.finite(y) & is.finite(weights) & weights > 0
  x <- x[ok]
  y <- y[ok]
  weights <- weights[ok]
  weights <- weights/sum(weights)
  if (is.null(h_x)) {
    h_x <- 1.1 * stats::bw.nrd0(x)
  }
  if (is.null(h_y)) {
    h_y <- 1.1 * stats::bw.nrd0(y)
  }
  if (!is.finite(h_x) || h_x <= 0) {
    h_x <- sd(x) * 0.3
  }
  if (!is.finite(h_y) || h_y <= 0) {
    h_y <- sd(y) * 0.3
  }
  xr <- range(x)
  yr <- range(y)
  xpad <- diff(xr) * expand
  ypad <- diff(yr) * expand
  gx <- seq(xr[1] - xpad, xr[2] + xpad, length.out = n)
  gy <- seq(yr[1] - ypad, yr[2] + ypad, length.out = n)
  Z <- matrix(0, nrow = length(gx), ncol = length(gy))
  norm_const <- 1/(2 * pi * h_x * h_y)
  for (i in seq_along(x)) {
    kx <- exp(-0.5 * ((gx - x[i])/h_x)^2)
    ky <- exp(-0.5 * ((gy - y[i])/h_y)^2)
    Z <- Z + weights[i] * outer(kx, ky)
  }
  Z <- Z * norm_const
  GRID <- expand.grid(x = gx, y = gy)
  GRID$density <- as.vector(Z)
  GRID$density_scaled <- GRID$density/max(GRID$density, na.rm = TRUE)
  list(grid = GRID, h_x = h_x, h_y = h_y)
}

