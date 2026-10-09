suppressPackageStartupMessages(library(lme4))
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(script)) dirname(normalizePath(script)) else getwd()
out <- file.path(root, "replayed")
dir.create(out, showWarnings = FALSE)
d <- read.csv(file.path(root, "native_clock_predictions.csv"), check.names = FALSE)
stopifnot(nrow(d) == 192L, all(is.finite(d$score)))
results <- list(); diagnostics <- list()
for (clock_name in unique(d$clock)) for (control_set in c("matched", "pooled")) {
  rows <- list()
  for (arm_name in c("INH", "MIM")) {
    controls <- if (control_set == "pooled") c("inhCTRL", "mimCTRL") else if (arm_name == "INH") "inhCTRL" else "mimCTRL"
    dat <- subset(d, clock == clock_name & GROUP %in% c(controls, arm_name))
    stopifnot(!anyDuplicated(dat$SHORT_ID))
    dat$group <- factor(ifelse(dat$GROUP == arm_name, "treated", "ctrl"), levels=c("ctrl", "treated"))
    dat$time_f <- factor(dat$TIMEF)
    dat$BATCH <- factor(dat$BATCH)
    fits <- list(
      no_group=lmer(score ~ time_f + (1|BATCH), dat, REML=FALSE),
      shift=lmer(score ~ group + time_f + (1|BATCH), dat, REML=FALSE),
      poly_shift=lmer(score ~ group + poly(TIME_H,2) + (1|BATCH), dat, REML=FALSE),
      poly_interaction=lmer(score ~ group * poly(TIME_H,2) + (1|BATCH), dat, REML=FALSE))
    key <- paste(clock_name, control_set, arm_name, sep="_")
    rows[[arm_name]] <- data.frame(clock=clock_name, controls=control_set, arm=arm_name,
      n=nrow(dat), n_runs=nlevels(dat$BATCH), n_control=sum(dat$group=="ctrl"), n_treated=sum(dat$group=="treated"),
      estimate=unname(fixef(fits$shift)["grouptreated"]),
      p_shift=anova(fits$no_group,fits$shift)$`Pr(>Chisq)`[2],
      p_poly=anova(fits$poly_shift,fits$poly_interaction)$`Pr(>Chisq)`[2])
    for (nm in names(fits)) diagnostics[[paste(key,nm)]] <- data.frame(clock=clock_name,controls=control_set,arm=arm_name,fit=nm,
      singular=isSingular(fits[[nm]]), optimizer_code=paste(fits[[nm]]@optinfo$conv$opt,collapse=";"),
      message=paste(fits[[nm]]@optinfo$conv$lme4$messages,collapse=";"))
  }
  z <- do.call(rbind,rows)
  q <- p.adjust(c(z$p_shift,z$p_poly),"BH")
  z$q_shift <- q[1:2]; z$q_poly <- q[3:4]
  results[[paste(clock_name,control_set)]] <- z
}
results <- do.call(rbind,results); rownames(results) <- NULL
frozen <- read.csv(file.path(root, "global_tests.csv"))
check <- merge(results, frozen, by=c("clock","controls","arm"), suffixes=c("_replayed","_frozen"))
check$log_q_error <- abs(log(check$q_shift_replayed/check$q_shift_frozen))
stopifnot(nrow(check)==8L, max(check$log_q_error)<1e-4)
write.csv(results,file.path(out,"global_tests.csv"),row.names=FALSE)
write.csv(check,file.path(out,"replay_validation.csv"),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,"model_diagnostics.csv"),row.names=FALSE)

# Preserve original DESeq2 inference: extract genes only after the genome-wide BH check.
genes <- c("EDA2R", "PTCHD4", "DUSP6", "EREG")
de <- lapply(c("INHvsCTRL", "MIMvsCTRL"), function(contrast) {
  z <- read.csv(file.path(root,paste0("DESeq2_",contrast,"_original.csv")))
  keep <- is.finite(z$padj)
  stopifnot(max(abs(p.adjust(z$pvalue[keep], "BH")-z$padj[keep]))<1e-14)
  a <- z[match(genes,z$SYMBOL), ]
  stopifnot(identical(a$SYMBOL,genes))
  data.frame(gene=a$SYMBOL,Ensembl_gene_ID=a$Ensembl_gene_ID,contrast=contrast,
    log2_fold_change=a$log2FoldChange,P_value=a$pvalue,q_global_BH=a$padj,
    genes_retained_for_BH=sum(keep))
})
write.csv(do.call(rbind,de),file.path(out,"Fig3i_global_DESeq2_matched.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(results,digits=12,row.names=FALSE)
cat("All clock values and genome-wide DESeq2 BH families reproduced.\n")
