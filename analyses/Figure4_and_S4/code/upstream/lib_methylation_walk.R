# Historical walk algorithms; mechanically reformatted, no analytical changes.
# Source: methylation_walk_omd.R. Only reusable functions are retained.
distance_label <- function(mode = c("RMS", "MS")) {
    mode <- match.arg(mode)
    if (mode == "RMS") 
        "RMS ΔM"
    else "MS ΔM²"
} 

vector_distance <- function(x, mode = c("RMS", "MS")) {
    mode <- match.arg(mode)
    if (mode == "RMS") {
        sqrt(mean(x^2, na.rm = TRUE))
    }
    else {
        mean(x^2, na.rm = TRUE)
    }
} 

step_distances <- function(x, mode = c("RMS", "MS")) {
    mode <- match.arg(mode)
    if (mode == "RMS") {
        sqrt(colMeans(x^2, na.rm = TRUE))
    }
    else {
        colMeans(x^2, na.rm = TRUE)
    }
} 

cosine_similarity <- function(a, b) {
    den <- sqrt(sum(a^2, na.rm = TRUE) * sum(b^2, na.rm = TRUE))
    if (!is.finite(den) || den <= 0) 
        return(NA_real_)
    sum(a * b, na.rm = TRUE)/den
} 

rotate_xy <- function(xy, theta) {
    R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE)
    as.matrix(xy) %*% t(R)
} 

make_spline_path <- function(df, n = SMOOTH_CURVE_N) {
    df <- df %>% dplyr::arrange(time)
    cumulative_arc_length <- c(0, cumsum(sqrt(diff(df$x_plot)^2 + diff(df$y_plot)^2)))
    if (nrow(df) < 3 || length(unique(cumulative_arc_length)) < 3) {
        return(df %>% dplyr::transmute(condition, experiment, arm, curve_order = dplyr::row_number(), 
            x_plot, y_plot))
    }
    interpolation_axis <- seq(min(cumulative_arc_length), max(cumulative_arc_length), length.out = n)
    spline_x <- stats::spline(x = cumulative_arc_length, y = df$x_plot, xout = interpolation_axis, method = "natural")
    spline_y <- stats::spline(x = cumulative_arc_length, y = df$y_plot, xout = interpolation_axis, method = "natural")
    tibble::tibble(condition = df$condition[1], experiment = df$experiment[1], arm = df$arm[1], curve_order = seq_along(interpolation_axis), 
        x_plot = spline_x$y, y_plot = spline_y$y)
} 

probe_base <- function(x) {
    sub("_.*$", "", as.character(x))
} 

resolve_probe_ids <- function(requested_ids, matrix_ids) {
    requested_ids <- unique(as.character(requested_ids))
    requested_ids <- requested_ids[!is.na(requested_ids)]
    exact <- intersect(requested_ids, matrix_ids)
    requested_base <- unique(probe_base(requested_ids))
    matrix_base <- probe_base(matrix_ids)
    unique(c(exact, matrix_ids[matrix_base %in% requested_base]))
} 

finite_probe_ids <- function(mval, sample_ids = colnames(mval)) {
    sample_ids <- intersect(sample_ids, colnames(mval))
    keep <- rowSums(is.finite(mval[, sample_ids, drop = FALSE])) == length(sample_ids)
    rownames(mval)[keep]
} 

select_top_variable <- function(mval, sample_ids, n_top) {
    sample_ids <- intersect(sample_ids, colnames(mval))
    finite_ids <- finite_probe_ids(mval, sample_ids)
    x <- mval[finite_ids, sample_ids, drop = FALSE]
    rv <- matrixStats::rowVars(x)
    ord <- order(rv, decreasing = TRUE, na.last = NA)
    ord <- ord[seq_len(min(length(ord), as.integer(n_top)))]
    selected <- rownames(x)[ord]
    intersect(selected, finite_probe_ids(mval))
} 

prepare_matrix <- function(mval, probes) {
    probes <- intersect(unique(probes), rownames(mval))
    if (length(probes) == 0) 
        stop("No requested probes found.")
    x <- mval[probes, , drop = FALSE]
    keep <- rowSums(is.finite(x)) == ncol(x)
    x <- x[keep, , drop = FALSE]
    if (nrow(x) < 10) 
        stop("Too few finite probes remain.")
    x
} 

condition_to_experiment <- function(x) {
    dplyr::case_when(grepl("^inh_", x) ~ "Inhibitor experiment", grepl("^mim_", x) ~ "Mimic experiment", 
        TRUE ~ NA_character_)
} 

condition_to_arm <- function(x) {
    dplyr::case_when(x %in% c("inh_CTRL", "mim_CTRL") ~ "Matched CTRL", x == "inh_INH" ~ "miR-29 inhibitor", 
        x == "mim_MIM" ~ "miR-29 mimic", TRUE ~ x)
} 

build_centroids <- function(mval, anno, conditions = CONDITIONS, times = TIME_LEVELS) {
    grid <- expand.grid(condition = conditions, time = times, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid <- grid[order(match(grid$condition, conditions), match(grid$time, times)), , drop = FALSE]
    centroids <- matrix(NA_real_, nrow = nrow(mval), ncol = nrow(grid), dimnames = list(rownames(mval), 
        paste(grid$condition, paste0("T", grid$time), sep = "__")))
    for (i in seq_len(nrow(grid))) {
        idx <- anno$condition_walk == grid$condition[i] & anno$time_walk == grid$time[i]
        centroids[, i] <- rowMeans(mval[, idx, drop = FALSE])
    }
    centroids
} 

build_steps <- function(centroids, conditions = CONDITIONS, times = TIME_LEVELS) {
    step_meta <- dplyr::bind_rows(lapply(conditions, function(cc) {
        tibble::tibble(condition = cc, from = times[-length(times)], to = times[-1])
    }))
    step_meta$step_id <- paste0(step_meta$condition, "__T", step_meta$from, "_to_T", step_meta$to)
    step_vectors <- matrix(NA_real_, nrow = nrow(centroids), ncol = nrow(step_meta), dimnames = list(rownames(centroids), 
        step_meta$step_id))
    for (i in seq_len(nrow(step_meta))) {
        from_col <- paste(step_meta$condition[i], paste0("T", step_meta$from[i]), sep = "__")
        to_col <- paste(step_meta$condition[i], paste0("T", step_meta$to[i]), sep = "__")
        step_vectors[, i] <- centroids[, to_col] - centroids[, from_col]
    }
    list(step_vectors = step_vectors, step_meta = step_meta)
} 

embed_steps <- function(step_vectors, step_lengths) {
    n_cpg <- nrow(step_vectors)
    gram <- crossprod(step_vectors)/n_cpg
    eig <- eigen(gram, symmetric = TRUE)
    vals <- pmax(eig$values[1:2], 0)
    xy <- eig$vectors[, 1:2, drop = FALSE] %*% diag(sqrt(vals), nrow = 2)
    raw_len <- sqrt(rowSums(xy^2))
    good <- raw_len > 1e-12 & is.finite(raw_len)
    xy[good, ] <- xy[good, , drop = FALSE]/raw_len[good] * step_lengths[good]
    if (any(!good)) {
        theta <- seq(0, 2 * pi, length.out = sum(!good) + 1)[-1]
        xy[!good, 1] <- cos(theta) * step_lengths[!good]
        xy[!good, 2] <- sin(theta) * step_lengths[!good]
    }
    colnames(xy) <- c("dx", "dy")
    xy
} 

make_walk <- function(mval, anno, feature_id, feature_label, distance_mode = c("RMS", "MS")) {
    distance_mode <- match.arg(distance_mode)
    mval <- prepare_matrix(mval, rownames(mval))
    centroids <- build_centroids(mval, anno)
    st <- build_steps(centroids)
    step_lengths <- step_distances(st$step_vectors, distance_mode)
    st$step_meta$step_length <- step_lengths
    step_xy <- embed_steps(st$step_vectors, step_lengths)
    path_df <- dplyr::bind_rows(lapply(CONDITIONS, function(cc) {
        xy <- matrix(0, nrow = length(TIME_LEVELS), ncol = 2)
        for (k in seq_len(length(TIME_LEVELS) - 1)) {
            idx <- which(st$step_meta$condition == cc & st$step_meta$from == TIME_LEVELS[k] & st$step_meta$to == 
                TIME_LEVELS[k + 1])
            xy[k + 1, ] <- xy[k, ] + step_xy[idx, ]
        }
        tibble::tibble(condition = cc, time = TIME_LEVELS, time_label = paste0("T", TIME_LEVELS), x = xy[, 
            1], y = xy[, 2])
    }))
    ref <- path_df %>% dplyr::filter(condition == "inh_CTRL", time == max(TIME_LEVELS))
    angle <- atan2(ref$y, ref$x)
    path_df[, c("x", "y")] <- rotate_xy(path_df[, c("x", "y")], -angle)
    path_df <- path_df %>% dplyr::mutate(x_plot = x * PLOT_SCALE, y_plot = y * PLOT_SCALE, experiment = condition_to_experiment(condition), 
        arm = condition_to_arm(condition))
    ref_net <- centroids[, paste("inh_CTRL", paste0("T", max(TIME_LEVELS)), sep = "__")] - centroids[, 
        paste("inh_CTRL", paste0("T", min(TIME_LEVELS)), sep = "__")]
    metrics <- dplyr::bind_rows(lapply(CONDITIONS, function(cc) {
        idx_steps <- st$step_meta$condition == cc
        path_length <- sum(st$step_meta$step_length[idx_steps])
        net <- centroids[, paste(cc, paste0("T", max(TIME_LEVELS)), sep = "__")] - centroids[, paste(cc, 
            paste0("T", min(TIME_LEVELS)), sep = "__")]
        net_displacement <- vector_distance(net, distance_mode)
        tibble::tibble(feature_id = feature_id, feature_space = feature_label, distance_mode = distance_mode, 
            distance_label = distance_label(distance_mode), condition = cc, n_cpg = nrow(mval), path_length = path_length, 
            net_displacement = net_displacement, straightness = net_displacement/path_length, cosine_to_inhibitor_ctrl = cosine_similarity(net, 
                ref_net))
    })) %>% dplyr::mutate(experiment = condition_to_experiment(condition), arm = condition_to_arm(condition)) %>% 
        dplyr::group_by(experiment) %>% dplyr::mutate(path_length_ratio_to_ctrl = path_length/path_length[arm == 
        "Matched CTRL"], net_displacement_ratio_to_ctrl = net_displacement/net_displacement[arm == "Matched CTRL"]) %>% 
        dplyr::ungroup()
    curve_df <- dplyr::bind_rows(lapply(split(path_df, path_df$condition), make_spline_path, n = SMOOTH_CURVE_N))
    rm(centroids, st)
    invisible(gc())
    list(feature_id = feature_id, feature_space = feature_label, distance_mode = distance_mode, distance_label = distance_label(distance_mode), 
        n_cpg = nrow(mval), path_df = path_df, curve_df = curve_df, metrics = metrics)
}
