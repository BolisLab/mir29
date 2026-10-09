suppressPackageStartupMessages({library(data.table);library(lme4);library(emmeans)})
htert_bundle <- input("clocks/htert_clocks.rds")
htert_meta <- htert_bundle$predictions
# New inference on the NEW hTERT scores; no inherited primary-HGF P values.
rows<-posts<-paired<-diagnostics<-list()
for(clock in c('RNAAgeCalc','tAge','PASTA')) for(armname in c('INH','MIM')) {
 dd<-htert_meta[htert_meta$arm==armname,,drop=FALSE]
 dd$score<-dd[[clock]];dd$y<-as.numeric(scale(dd$score));sdy<-sd(dd$score)
 dd$condition<-factor(ifelse(dd$is_control=='CTRL','control','treated'),c('control','treated'))
 dd$tp<-factor(dd$TIME);dd$run<-factor(dd$BATCH)
 stopifnot(all(table(dd$run,dd$condition,dd$tp)==1))
 fits<-list(null=lmer(y~tp+(1|run),dd,REML=FALSE),
            shift=lmer(y~condition+tp+(1|run),dd,REML=FALSE),
            full=lmer(y~condition*tp+(1|run),dd,REML=FALSE))
 # Three TPs: a quadratic interaction spans exactly the categorical interaction.
 pshift<-anova(fits$null,fits$shift)$`Pr(>Chisq)`[2]
 pint<-anova(fits$shift,fits$full)$`Pr(>Chisq)`[2]
 kr<-pbkrtest::KRmodcomp(update(fits$shift,REML=TRUE),update(fits$null,REML=TRUE))$stats$p.value
 means<-as.data.table(dd)[,.(score=mean(score)),by=.(run,condition)]
 bw<-dcast(means,run~condition,value.var='score');bw[,delta:=treated-control]
 tt<-t.test(bw$delta)
 rows[[paste(clock,armname)]]<-data.table(clock,arm=armname,n_runs=nlevels(dd$run),n_samples=nrow(dd),
   estimate_shift=fixef(fits$shift)['conditiontreated']*sdy,p_global_shift=pshift,
   p_interaction=pint,p_KR_average_shift=kr,p_runmean_paired=tt$p.value,
   runmean_difference=mean(bw$delta),runmean_lower95=tt$conf.int[1],runmean_upper95=tt$conf.int[2])
 post<-as.data.table(as.data.frame(contrast(emmeans(fits$full,~condition|tp,lmer.df='kenward-roger'),'revpairwise')))
 post[,`:=`(clock=clock,arm=armname,estimate=estimate*sdy,SE=SE*sdy,p_Holm_3TP=p.adjust(p.value,'holm'))]
 posts[[paste(clock,armname)]]<-post
 bw[,`:=`(clock=clock,arm=armname)];paired[[paste(clock,armname)]]<-bw
 diagnostics[[paste(clock,armname)]]<-rbindlist(lapply(names(fits),function(n)
  data.table(clock,arm=armname,fit=n,singular=isSingular(fits[[n]]),
   messages=paste(fits[[n]]@optinfo$conv$lme4$messages,collapse='; '))))
}
global<-rbindlist(rows);post<-rbindlist(posts)
global[,`:=`(q_global_BH4=NA_real_,q_interaction_BH4=NA_real_)]
for(n in unique(global$clock)) {
 ix<-which(global$clock==n);q<-p.adjust(c(global$p_global_shift[ix],global$p_interaction[ix]),'BH')
 global$q_global_BH4[ix]<-q[1:2];global$q_interaction_BH4[ix]<-q[3:4]
}
global[,q_runmean_BH2:=p.adjust(p_runmean_paired,'BH'),by=clock]
post[,q_BH_6TP:=p.adjust(p.value,'BH'),by=clock]

save_table(global, "S4k_clock_global_statistics")
save_table(post, "S4k_clock_timepoint_statistics")
save_table(rbindlist(diagnostics), "S4k_model_diagnostics")
old <- as.data.table(htert_bundle$global)
ch <- merge(global, old, by = c("clock", "arm"), suffixes = c("_new", "_old"))
assert_close(ch$q_global_BH4_new, ch$q_global_BH4_old, "hTERT global clock q values", 1e-5)
display <- as.data.frame(htert_bundle$display)
display$TP2 <- display$TIME
display$series <- as.character(display$GROUP)
clock_plots <- lapply(c("RNAAgeCalc", "tAge", "PASTA"), function(nm) trajectory(subset(display, clock == nm), ifelse(nm == "tAge", "Relative tAge (model-scaled units)", paste(nm, "score")), xpos = 1:3, title = nm))
save_panel(plot_grid(plotlist = clock_plots, nrow = 1), "S4k_hTERT_transcriptomic_clocks", 12, 4.5)
save_table(display, "S4k_clock_display_values")
