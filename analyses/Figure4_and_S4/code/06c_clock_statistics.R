suppressPackageStartupMessages({library(data.table);library(lme4);library(pbkrtest);library(emmeans)})
d <- as.data.table(input("clocks/primary_clocks.rds")$scores)
# Use own controls and a biological-run random intercept. The primary clock
# test is the average treatment offset adjusted for categorical TP, as in main.
# A separate interaction test is retained, not substituted according to P.
rows <- tp_rows <- diag_rows <- list()
for(n in unique(d$clock)) for(arm in c('INH','MIM')) {
 dd<-copy(d[clock==n & tolower(as.character(EXPERIMENT))==tolower(arm)])
 dd[,`:=`(condition=factor(ifelse(GROUP=='CTRL','control','treated'),c('control','treated')),
          time=factor(TP2,levels=1:6),run=droplevels(run))]
 stopifnot(nrow(dd)==36,all(table(dd$run,dd$time,dd$condition)==1))
 dd[,y:=as.numeric(scale(score))]
 fit<-lmer(y~condition+time+(1|run),dd,REML=TRUE)
 null<-lmer(y~time+(1|run),dd,REML=TRUE)
 full<-lmer(y~condition*time+(1|run),dd,REML=TRUE)
 global<-KRmodcomp(fit,null)$stats; inter<-KRmodcomp(full,fit)$stats
 rows[[paste(n,arm)]]<-data.table(clock=n,arm=arm,n_experiments=3,
   p_average_shift=global$p.value,p_time_interaction=inter$p.value,
   estimate_shift=fixef(fit)['conditiontreated']*sd(dd$score),
   df_average=global$ddf,df_interaction=inter$ddf,
   singular_average=isSingular(fit),singular_interaction=isSingular(full))
 em<-emmeans(full,~condition|time,lmer.df='kenward-roger')
 ct<-as.data.frame(contrast(em,'revpairwise'))
 ct$clock<-n;ct$arm<-arm;ct$estimate<-ct$estimate*sd(dd$score);ct$SE<-ct$SE*sd(dd$score)
 tp_rows[[paste(n,arm)]]<-ct
}
stats<-rbindlist(rows); stats[,q_BH_two_arms:=p.adjust(p_average_shift,'BH'),by=clock]
stats[,q_interaction_two_arms:=p.adjust(p_time_interaction,'BH'),by=clock]
stats[,q_BH_all_clocks:=p.adjust(p_average_shift,'BH')]
tp<-rbindlist(tp_rows);tp[,q_BH_12_TP_tests:=p.adjust(p.value,'BH'),by=clock]

save_table(stats, "S4_epigenetic_clock_global_statistics")
save_table(tp, "S4_epigenetic_clock_timepoint_statistics")
