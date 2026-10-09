# S4h: original DMR calls, matched response and aligned CpG heatmap.
COL <- COLORS
th <- theme_figure()
mini <- th + theme(text=element_text(size=8), axis.title=element_text(size=8),
  axis.text=element_text(size=7), plot.title=element_text(size=9,face='bold'),
  legend.text=element_text(size=7), legend.key.size=unit(.27,'cm'),
  legend.title=element_text(size=7), plot.margin=margin(5,7,5,4))

wd <- input("wnt/dmrcate_display_data.rds")
cores <- as.data.table(wd$cores)
LOC <- 'DKK2'
original_orders <- c('CHIR.vs.CTRL','INH.vs.CTRL','COMBO.vs.INH','COMBO.vs.CHIR','COMBO.vs.CTRL')
orders <- c('INH.vs.CTRL', 'CHIR.vs.CTRL', 'COMBO.vs.CTRL')
labels <- c('INH vs CTRL', 'CHIR vs CTRL', 'INH+CHIR vs CTRL')
stacks <- list()
audit <- list()

for (i in seq_along(LOC)) {
  gene <- LOC[i]
  core <- cores[locus==gene & anchor_source=='CHIR.vs.CTRL']
  stopifnot(nrow(core)==1)
  pp <- as.data.table(wd$profiles[[gene]])[BATCH %in% c('EXP66','EXP67','EXP68')]
  cc <- pp[GROUP=='CTRL', .(BATCH, position, control=beta)]
  dd <- merge(pp[GROUP!='CTRL'], cc, by=c('BATCH','position'))
  dd[, delta:=100*(beta-control)]
  dd[, GROUP:=factor(GROUP,c('INH','CHIR','COMBO'))]
  original_calls <- as.data.table(wd$calls)[locus==gene & contrast %in% original_orders]
  cl <- copy(original_calls[contrast %in% orders])
  stopifnot(setequal(cl$contrast, orders))
  # Retain the original full genomic window even though two tracks are hidden.
  xr <- range(c(pp$position,original_calls$start,original_calls$end),na.rm=TRUE)
  normx <- function(v) (v-xr[1])/diff(xr)
  sites <- data.table(position=sort(unique(pp$position)))
  sites[, `:=`(x=normx(position), hx=(seq_len(.N)-.5)/.N,
               core=position>=core$start & position<=core$end)]
  dd <- merge(dd,sites,by='position')
  cl[, `:=`(x=normx(start), xe=normx(end), y=length(orders)+1-match(contrast,orders))]
  br <- pretty(xr,n=3); br <- br[br>=xr[1] & br<=xr[2]]
  xs <- scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=normx(br),
    labels=format(br,scientific=FALSE,big.mark=',',trim=TRUE))
  a <- ggplot()+annotate('rect',xmin=normx(core$start),xmax=normx(core$end),
    ymin=-Inf,ymax=Inf,fill='#DAE5EB',alpha=.8)+
    geom_segment(data=cl,aes(x=0,xend=1,y=y,yend=y),colour='#DCE3E7',linewidth=.3)+
    geom_segment(data=cl[called==TRUE],aes(x=x,xend=xe,y=y,yend=y),colour='#CD6253',linewidth=2.5)+
    geom_text(data=cl[called==FALSE],aes(x=.52,y=y,label='No overlapping DMR'),size=2,colour='#6F8190')+
    xs+scale_y_continuous(breaks=seq_along(orders),labels=rev(labels),
      limits=c(.4,length(orders)+.6),expand=c(0,0))+
    labs(x=paste0(core$seqnames,' position (hg38)'),y=NULL,title='Original DMRcate calls')+
    mini+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.8))
  sm <- dd[,.(delta=mean(delta),sem=sd(delta)/sqrt(.N)),by=.(position,x,GROUP)]
  b <- ggplot(sm,aes(x,delta,colour=GROUP,group=GROUP))+
    annotate('rect',xmin=normx(core$start),xmax=normx(core$end),ymin=-Inf,ymax=Inf,fill='#DAE5EB',alpha=.5)+
    geom_hline(yintercept=0,colour='grey65',linetype='dashed')+
    geom_errorbar(aes(ymin=delta-sem,ymax=delta+sem),width=.009,linewidth=.3)+
    geom_line(linewidth=.5)+geom_point(size=1.4)+xs+
    scale_colour_manual(values=COL,labels=c(INH='INH',CHIR='CHIR',COMBO='INH+CHIR'))+
    labs(x=NULL,y='Delta beta (pp)',colour=NULL,title='Matched-control profile: mean +/- SEM')+
    mini+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='bottom')
  con <- ggplot(sites)+geom_segment(aes(x=x,xend=hx,y=1,yend=0),colour='#91A7B4',linewidth=.3)+
    scale_x_continuous(limits=c(0,1),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1),expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(0,7,0,4))
  dd[,row:=paste(ifelse(GROUP=='COMBO','INH+CHIR',as.character(GROUP)),BATCH)]
  ro <- as.vector(t(outer(c('INH','CHIR','INH+CHIR'),c('EXP66','EXP67','EXP68'),paste)))
  dd[,row:=factor(row,rev(ro))]
  lim <- max(5,ceiling(max(abs(dd$delta))/5)*5)
  h <- ggplot(dd,aes(hx,row,fill=delta))+geom_tile(width=.97/nrow(sites),height=.94)+
    geom_tile(data=dd[core==TRUE],width=.97/nrow(sites),height=.94,fill=NA,colour='#7B91A0',linewidth=.3)+
    scale_fill_gradient2(low='#427DB5',mid='white',high='#CF654E',limits=c(-lim,lim),
      breaks=c(-lim,0,lim),name='Delta beta (pp)')+
    scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=sites$hx,labels=seq_len(nrow(sites)))+
    labs(x='CpGs in genomic order; CHIR core outlined',y=NULL)+mini+
    theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=6.8),
      axis.text.x=element_text(size=5.5))
  body <- plot_grid(a,b,con,h,ncol=1,rel_heights=c(1.08,1.3,.27,1.4),align='v',axis='lr')
  title <- ggdraw()+draw_label(gene,size=13,fontface='bold',colour='#263D4C')
  stacks[[gene]] <- plot_grid(title,body,ncol=1,rel_heights=c(.26,4.05))
  audit[[gene]] <- list(genomic_window=xr,sites=sites,profile=sm,heatmap=dd,tracks=cl,
    lower_plots=list(profile=b,connectors=con,heatmap=h))
}

save_panel(stacks[["DKK2"]], "S4h_DKK2_aligned_locus", 6.8, 8.0)
save_table(audit[["DKK2"]]$profile, "S4h_DKK2_profile")
