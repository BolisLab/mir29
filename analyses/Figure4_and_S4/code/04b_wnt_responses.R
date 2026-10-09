# S4i: same probes in every sample, three complete biological experiments.
w <- input("wnt/dmrcate_display_data.rds")
groups <- c('INH','CHIR','COMBO')
batches <- c('EXP66','EXP67','EXP68')
p <- as.data.table(w$profiles[['DKK2']])[BATCH %in% batches & GROUP %in% c('CTRL',groups)]
core <- as.data.table(w$cores)[locus=='DKK2' & anchor_source=='CHIR.vs.CTRL']
stopifnot(nrow(core)==1,!anyDuplicated(p[,.(BATCH,GROUP,position)]),all(is.finite(p$beta)))
sites <- sort(unique(p$position))
regions <- list('Selected core'=sites[sites>=core$start & sites<=core$end],
                'Entire region'=sites)
stopifnot(length(regions[[1]])==3L,length(regions[[2]])==16L,
          all(regions[[1]] %in% regions[[2]]))
dat <- rbindlist(lapply(names(regions),function(region){
  a <- p[position %in% regions[[region]],
         .(mean_beta=mean(beta),n_cpg=.N),by=.(BATCH,GROUP)]
  stopifnot(nrow(a)==12L,all(a$n_cpg==length(regions[[region]])))
  cc <- a[GROUP=='CTRL',.(BATCH,control_beta=mean_beta)]
  d <- merge(a[GROUP!='CTRL'],cc,by='BATCH')
  d[,`:=`(delta_pp=100*(mean_beta-control_beta),feature=region)]
  d
}))
sm <- dat[,.(mean_pp=mean(delta_pp),sem_pp=sd(delta_pp)/sqrt(.N),n=.N),by=.(feature,GROUP)]
stopifnot(nrow(dat)==18L,all(sm$n==3L))
dat[,feature:=factor(feature,levels=names(regions))]
sm[,feature:=factor(feature,levels=names(regions))]
dat[,xgroup:=match(GROUP,groups)]
dat[,xpoint:=xgroup+c(EXP66=-.055,EXP67=0,EXP68=.055)[BATCH]]
sm[,xgroup:=match(GROUP,groups)]

cols <- c(INH='#248636',CHIR='#3F7FB5',COMBO='#D28A32')
labels <- c(INH='INH',CHIR='CHIR',COMBO='INH+CHIR')
ink <- '#263D4C'
facet_labels <- c('Selected core'='DKK2\n3-CpG core',
                  'Entire region'='DKK2\n16-CpG extended region')
upper <- ceiling(max(c(dat$delta_pp,sm$mean_pp+sm$sem_pp))/5)*5
plot <- ggplot(dat,aes(x=xpoint,y=delta_pp))+
  geom_hline(yintercept=0,linetype='dashed',colour='#A5B2B9',linewidth=.4)+
  geom_line(aes(group=BATCH),colour='#B5C3CA',linewidth=.6,alpha=.85)+
  geom_point(aes(fill=GROUP),colour='white',shape=21,size=2.8,stroke=.35)+
  geom_errorbar(data=sm,aes(x=xgroup,y=mean_pp,ymin=mean_pp-sem_pp,ymax=mean_pp+sem_pp,
                          colour=GROUP),width=.14,linewidth=.7,inherit.aes=FALSE)+
  geom_point(data=sm,aes(x=xgroup,y=mean_pp,fill=GROUP),shape=23,colour='white',
             stroke=.6,size=4.6,inherit.aes=FALSE)+
  facet_wrap(~feature,nrow=1,labeller=as_labeller(facet_labels))+
  scale_fill_manual(values=cols,breaks=groups,labels=labels)+
  scale_colour_manual(values=cols,breaks=groups,labels=labels)+
  scale_x_continuous(limits=c(.65,3.4),breaks=1:3,labels=labels[groups],expand=c(0,0))+
  scale_y_continuous(limits=c(-1,upper),breaks=seq(0,upper,5),expand=c(0,0))+
  labs(x=NULL,y='Methylation change vs matched CTRL\n(percentage points)')+
  theme_classic(base_size=11,base_family='Helvetica')+
  theme(text=element_text(colour=ink),axis.text=element_text(colour=ink),
        axis.line=element_line(colour='#8C9BA4',linewidth=.35),
        axis.ticks=element_line(colour='#8C9BA4',linewidth=.35),
        legend.position='none',plot.margin=margin(12,16,10,10),
        strip.background=element_rect(fill='#F0F5F6',colour=NA),
        strip.text=element_text(size=10,colour=ink,lineheight=1.15),
        panel.spacing=grid::unit(6,'mm'),axis.text.x=element_text(size=9.5))

save_panel(plot, "S4i_DKK2_paired_responses", 8.8, 5.3)
save_table(dat, "S4i_DKK2_paired_responses")
