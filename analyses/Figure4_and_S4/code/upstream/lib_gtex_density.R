# Historical chronological-age-window density calculation; formatting only.
# Original script did not save a random seed. The driver records a new fixed seed.
meta_curves_fast_absolute <- function(meth_mat, annot_list, tissues = names(annot_list), age_col = "AGE", 
    win_size = 15, step_size = 5, target_n = 20000, min_vals = 1000, nbins_beta = 256, n_age_bins = 15, 
    smooth_sigma = 0.01, min_windows_per_bin = 6) {
    stopifnot(is.matrix(meth_mat))
    storage.mode(meth_mat) <- "double"
    beta_breaks <- seq(0, 1, length.out = nbins_beta + 1)
    beta_centers <- (beta_breaks[-1] + beta_breaks[-length(beta_breaks)])/2
    bin_w <- 1/nbins_beta
    all_ages <- unlist(lapply(annot_list, function(a) {
        if (!is.null(a) && age_col %in% colnames(a)) 
            as.numeric(a[[age_col]])
        else NULL
    }))
    all_ages <- all_ages[is.finite(all_ages)]
    global_min <- min(all_ages, na.rm = TRUE)
    global_max <- max(all_ages, na.rm = TRUE)
    if (!is.finite(global_min) || !is.finite(global_max) || global_max <= global_min) 
        stop("Invalid age range across tissues")
    age_breaks <- seq(global_min, global_max, length.out = n_age_bins + 1)
    age_centers <- (age_breaks[-1] + age_breaks[-length(age_breaks)])/2
    meta_counts <- matrix(0, nrow = nbins_beta, ncol = n_age_bins)
    win_counts <- integer(n_age_bins)
    bin_index_beta <- function(x, nb) {
        i <- floor(x * nb) + 1L
        i[i < 1L] <- 1L
        i[i > nb] <- nb
        i
    }
    smooth_reflect <- function(y, sigma_beta) {
        if (!is.finite(sigma_beta) || sigma_beta <= 0) 
            return(y)
        sigma_bins <- max(0.5, sigma_beta/bin_w)
        half <- max(3L, ceiling(3 * sigma_bins))
        kx <- (-half):half
        ker <- exp(-0.5 * (kx/sigma_bins)^2)
        ker <- ker/sum(ker)
        n <- length(y)
        ypad <- c(rev(y[1:half]), y, rev(y[(n - half + 1):n]))
        yf <- stats::filter(ypad, ker, sides = 2)
        out <- as.numeric(yf[(half + 1):(half + n)])
        s <- sum(out) * bin_w
        if (s > 0) 
            out <- out/s
        out
    }
    for (t in tissues) {
        ann <- annot_list[[t]]
        if (is.null(ann) || is.null(rownames(ann)) || !(age_col %in% colnames(ann))) 
            next
        samp <- intersect(rownames(ann), colnames(meth_mat))
        if (length(samp) < win_size) 
            next
        ages <- ann[samp, age_col]
        o <- order(ages, na.last = NA)
        if (length(o) < win_size) 
            next
        B <- meth_mat[, samp[o], drop = FALSE]
        A <- as.numeric(ages[o])
        starts <- seq(1, ncol(B) - win_size + 1, by = step_size)
        if (length(starts) == 0) 
            next
        for (i in starts) {
            idx <- i:(i + win_size - 1)
            beta_vec <- as.numeric(B[, idx, drop = FALSE])
            beta_vec <- beta_vec[is.finite(beta_vec)]
            n_avail <- length(beta_vec)
            if (n_avail < min_vals) 
                next
            take <- min(target_n, n_avail)
            if (n_avail > take) 
                beta_vec <- sample(beta_vec, size = take)
            ids <- bin_index_beta(beta_vec, nbins_beta)
            counts <- tabulate(ids, nbins = nbins_beta)
            age_c <- mean(A[idx], na.rm = TRUE)
            k <- findInterval(age_c, age_breaks, rightmost.closed = TRUE)
            if (k < 1L) 
                k <- 1L
            if (k > n_age_bins) 
                k <- n_age_bins
            meta_counts[, k] <- meta_counts[, k] + counts
            win_counts[k] <- win_counts[k] + 1L
        }
    }
    keep <- (colSums(meta_counts) > 0) & (win_counts >= min_windows_per_bin)
    if (!any(keep)) 
        stop("All age bins dropped; lower 'min_windows_per_bin' or reduce 'n_age_bins'.")
    meta_counts <- meta_counts[, keep, drop = FALSE]
    age_centers <- age_centers[keep]
    dens_mat <- apply(meta_counts, 2, function(y) {
        s <- sum(y) * bin_w
        if (s == 0) 
            return(y)
        y <- y/s
        smooth_reflect(y, smooth_sigma)
    })
    data.frame(beta = rep(beta_centers, times = ncol(dens_mat)), density = as.vector(dens_mat), age_mid = rep(age_centers, 
        each = nbins_beta), age_bin = factor(seq_len(ncol(dens_mat))[rep(seq_len(ncol(dens_mat)), each = nbins_beta)]))
}
