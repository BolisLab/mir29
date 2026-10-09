# S4i: paired beta-scale factorial contrast in each complete experiment.
# The extended region includes the core. BH adjustment covers these two tests.
w <- input("wnt/dmrcate_display_data.rds")
profile <- as.data.frame(w$profiles$DKK2)
profile <- subset(profile, BATCH %in% c("EXP66", "EXP67", "EXP68"))
core <- subset(w$cores, locus == "DKK2" & anchor_source == "CHIR.vs.CTRL")
all_sites <- sort(unique(profile$position))
regions <- list(`3-CpG core` = all_sites[all_sites >= core$start & all_sites <= core$end],
                `16-CpG extended` = all_sites)
interaction_runs <- bind_rows(lapply(names(regions), function(region) {
  d <- profile %>% filter(position %in% regions[[region]]) %>%
    group_by(BATCH, GROUP) %>% summarise(beta = mean(beta), .groups = "drop") %>%
    pivot_wider(names_from = GROUP, values_from = beta)
  data.frame(region, BATCH = d$BATCH,
             interaction_pp = 100 * ((d$COMBO - d$CHIR) - (d$INH - d$CTRL)))
}))
interaction_tests <- interaction_runs %>% group_by(region) %>% group_modify(~ {
  tt <- t.test(.x$interaction_pp, mu = 0)
  data.frame(n = nrow(.x), estimate_pp = mean(.x$interaction_pp),
             lower95 = tt$conf.int[1], upper95 = tt$conf.int[2], p = tt$p.value)
}) %>% ungroup() %>% mutate(q = p.adjust(p, "BH"),
                           region = factor(region, levels = rev(names(regions))))
p <- ggplot(interaction_tests, aes(estimate_pp, region)) +
  geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), orientation = "y", width = .15,
                color = "#C57B2D", linewidth = .7) +
  geom_point(shape = 23, fill = "#C57B2D", color = "white", size = 4) +
  geom_text(aes(x = upper95, label = sprintf("q = %.4f", q)), hjust = -.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(.1, .38))) +
  labs(x = "Interaction (Delta beta, percentage points)", y = NULL,
       title = "DKK2: paired interaction") + theme_figure()
save_panel(p, "S4i_DKK2_interaction_forest", 6.6, 2.6)
save_table(interaction_runs, "S4i_interactions_per_experiment")
save_table(interaction_tests, "S4i_interaction_tests")
expected <- as.data.frame(input("wnt/interaction_statistics_by_gene.rds"))
expected <- subset(expected, gene == "DKK2")
assert_close(sort(interaction_tests$p), sort(expected$paired_t_p), "DKK2 paired interaction p values")
