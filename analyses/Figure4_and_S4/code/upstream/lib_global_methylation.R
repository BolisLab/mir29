# Historical mixed-model functions, unchanged apart from formatting.
build_delta_contrasts <- function(emm_grid, time_values, baseline_time_hours) {
    emm_grid_df <- as.data.frame(emm_grid@grid)
    emm_grid_df$mGROUP <- as.character(emm_grid_df$mGROUP)
    emm_grid_df$TIME_factor <- as.character(emm_grid_df$TIME_factor)
    baseline_label <- as.character(baseline_time_hours)
    contrast_list <- list()
    add_weight <- function(weight_vector, group_name, time_label, weight_value) {
        matching_index <- which(emm_grid_df$mGROUP == group_name & emm_grid_df$TIME_factor == time_label)
        if (length(matching_index) != 1) {
            stop("Could not uniquely identify ", group_name, " at time ", time_label, ".")
        }
        weight_vector[matching_index] <- weight_vector[matching_index] + weight_value
        weight_vector
    }
    nonbaseline_times <- time_values[time_values != baseline_time_hours]
    for (current_time in nonbaseline_times) {
        current_label <- as.character(current_time)
        ctrl_weights <- rep(0, nrow(emm_grid_df))
        ctrl_weights <- add_weight(ctrl_weights, "inh_CTRL", current_label, 0.5)
        ctrl_weights <- add_weight(ctrl_weights, "inh_CTRL", baseline_label, -0.5)
        ctrl_weights <- add_weight(ctrl_weights, "mim_CTRL", current_label, 0.5)
        ctrl_weights <- add_weight(ctrl_weights, "mim_CTRL", baseline_label, -0.5)
        contrast_list[[paste("CTRL", current_time, sep = "__")]] <- ctrl_weights
        inhibitor_weights <- rep(0, nrow(emm_grid_df))
        inhibitor_weights <- add_weight(inhibitor_weights, "inh_INH", current_label, 1)
        inhibitor_weights <- add_weight(inhibitor_weights, "inh_INH", baseline_label, -1)
        contrast_list[[paste("Inhibitor", current_time, sep = "__")]] <- inhibitor_weights
        mimic_weights <- rep(0, nrow(emm_grid_df))
        mimic_weights <- add_weight(mimic_weights, "mim_MIM", current_label, 1)
        mimic_weights <- add_weight(mimic_weights, "mim_MIM", baseline_label, -1)
        contrast_list[[paste("Mimic", current_time, sep = "__")]] <- mimic_weights
    }
    contrast_list
} 

fit_own_baseline_metric <- function(data, metric_name) {
    data_metric <- data[is.finite(data[[metric_name]]), , drop = FALSE]
    observed_cells <- table(data_metric$mGROUP, data_metric$TIME_factor)
    if (any(observed_cells == 0)) {
        stop("At least one mGROUP × TIME combination is missing for ", metric_name, ".")
    }
    factor_formula <- stats::as.formula(paste0(metric_name, " ~ mGROUP * TIME_factor"))
    numeric_formula <- stats::as.formula(paste0(metric_name, " ~ mGROUP * TIME_FROM_T1_DAYS"))
    fit_factor <- nlme::lme(fixed = factor_formula, random = ~1 | BLOCK, data = data_metric, method = "REML", 
        na.action = stats::na.omit, control = nlme::lmeControl(opt = "optim", maxIter = 300, msMaxIter = 300, 
            returnObject = TRUE))
    fit_numeric <- nlme::lme(fixed = numeric_formula, random = ~1 | BLOCK, data = data_metric, method = "ML", 
        na.action = stats::na.omit, control = nlme::lmeControl(opt = "optim", maxIter = 300, msMaxIter = 300, 
            returnObject = TRUE))
    emm_full <- emmeans::emmeans(fit_factor, specs = ~mGROUP * TIME_factor, weights = "equal")
    delta_contrast_list <- build_delta_contrasts(emm_grid = emm_full, time_values = time_values, baseline_time_hours = baseline_time_hours)
    delta_contrast_object <- emmeans::contrast(emm_full, method = delta_contrast_list, adjust = "none")
    delta_df <- as.data.frame(confint(delta_contrast_object, level = confidence_level))
    lower_column <- intersect(c("lower.CL", "asymp.LCL"), colnames(delta_df))[1]
    upper_column <- intersect(c("upper.CL", "asymp.UCL"), colnames(delta_df))[1]
    if (is.na(lower_column) || is.na(upper_column)) {
        stop("Could not identify confidence interval columns for ", metric_name, ".")
    }
    delta_df$lower.CL <- delta_df[[lower_column]]
    delta_df$upper.CL <- delta_df[[upper_column]]
    parsed_names <- strsplit(as.character(delta_df$contrast), "__", fixed = TRUE)
    delta_df$GROUP <- vapply(parsed_names, function(x) {
        x[1]
    }, character(1))
    delta_df$TIME_HOURS <- as.numeric(vapply(parsed_names, function(x) {
        x[2]
    }, character(1)))
    delta_df <- delta_df %>% dplyr::mutate(GROUP = factor(GROUP, levels = plot_group_levels), metric = metric_name) %>% 
        dplyr::select(GROUP, TIME_HOURS, estimate, SE, df, lower.CL, upper.CL, metric)
    baseline_rows <- tibble::tibble(GROUP = factor(plot_group_levels, levels = plot_group_levels), TIME_HOURS = baseline_time_hours, 
        estimate = 0, SE = 0, df = NA_real_, lower.CL = 0, upper.CL = 0, metric = metric_name)
    trajectory_df <- dplyr::bind_rows(baseline_rows, delta_df) %>% dplyr::arrange(GROUP, TIME_HOURS)
    slope_emmeans <- emmeans::emtrends(fit_numeric, specs = ~mGROUP, var = "TIME_FROM_T1_DAYS")
    slope_methods <- list(`Inhibitor vs matched CTRL` = c(-1, 1, 0, 0), `Mimic vs matched CTRL` = c(0, 
        0, -1, 1))
    slope_contrasts <- emmeans::contrast(slope_emmeans, method = slope_methods, adjust = "none")
    slope_contrast_df <- as.data.frame(summary(slope_contrasts, infer = c(TRUE, TRUE), level = confidence_level, 
        adjust = "none"))
    list(metric = metric_name, fit_factor = fit_factor, fit_numeric = fit_numeric, emm_full = emm_full, 
        delta_contrast_object = delta_contrast_object, trajectory = trajectory_df, slope_emmeans = slope_emmeans, 
        slope_contrasts = slope_contrast_df)
}
