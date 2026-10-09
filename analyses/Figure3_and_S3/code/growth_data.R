# Original cell counts; columns identify treatment and independent experiment.
growth_raw <- list(
  INH = data.frame(time = c(0, 72, 96, 168, 192, 264, 288),
    ctrl_EXP15 = c(23360,46600,53200,75800,80200,81200,78800),
    ctrl_EXP25 = c(39000,47000,44000,59200,68200,78200,73400),
    ctrl_EXP27 = c(50200,109400,115800,148800,148600,142800,120500),
    inh_EXP15 = c(23360,55200,64000,76800,87800,88200,82200),
    inh_EXP25 = c(39000,43400,44600,60600,65200,71800,66400),
    inh_EXP27 = c(50200,95400,102800,121000,124800,126400,130600)),
  MIM = data.frame(time = c(0, 72, 96, 168, 192, 264, 288),
    ctrl_EXP45 = c(40800,63600,62600,68000,66400,89000,87200),
    ctrl_EXP49 = c(74600,152400,110000,157400,132600,111000,131200),
    ctrl_EXP52 = c(54000,69800,78400,83200,82200,90000,80800),
    mim_EXP45 = c(40800,49400,55800,59200,50400,55400,71400),
    mim_EXP49 = c(74600,142000,100800,120400,114800,NA,148800),
    mim_EXP52 = c(54000,70200,75800,76000,73000,75000,79200)))
normalize_growth <- function(raw) {
  do.call(rbind, lapply(setdiff(names(raw), "time"), function(id) {
    value <- raw[[id]]
    data.frame(time = raw$time, group = sub("_.*", "", id), experiment = sub(".*_", "", id),
      y = 100 * (value - value[1]) / (max(value, na.rm = TRUE) - value[1]))
  }))
}
fit_growth <- function(d) {
  do.call(rbind, lapply(split(d, d$group), function(z) {
    fit <- minpack.lm::nlsLM(y ~ plateau * (1 - exp(-k * time)), data = na.omit(z),
      start = list(plateau = max(z$y, na.rm = TRUE), k = .01), lower = c(0, 0), upper = c(200, 1))
    out <- data.frame(time = seq(0, 288, length.out = 250), group = z$group[1])
    out$yhat <- predict(fit, out)
    out
  }))
}
