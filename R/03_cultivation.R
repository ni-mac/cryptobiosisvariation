# Cultivation analysis

run_cultivation <- function(project_root, bootstrap_iterations=1000L) {
 suppressPackageStartupMessages({library(readxl);library(dplyr);library(tidyr);library(stringr);library(purrr);library(survival);library(ggplot2);library(patchwork)})
 out <- file.path(project_root,'Results/statistics/cultivation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
 cache<-file.path(project_root,'Results/cache/cultivation');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
 write_out<-function(x,nm)write.csv(x,file.path(out,nm),row.names=FALSE,na='')
 num<-function(x)suppressWarnings(as.numeric(as.character(x)))
 key<-function(df,cols)do.call(paste,c(df[cols],sep='|'))
 raw<-read.csv(file.path(project_root,'Results/inputs/cultivation_observations.csv'),colClasses='character',check.names=FALSE,na.strings='') %>% mutate(Source_snapshot_row=row_number()+1L,.before=1)
 line_keys<-c('Experiment_id','Strain','Salinity_PSU','Replicate');plate_keys<-c(line_keys,'Generation');obs_keys<-c(plate_keys,'Day')
 parse_status<-function(x)map_dbl(as.character(x),function(v){h<-str_extract_all(ifelse(is.na(v),'',v),'(?<![0-9])[0-7](?![0-9])')[[1]];if(length(h)==1)as.numeric(h) else NA_real_})
 mal<-raw %>% mutate(Day=num(Day),Salinity_PSU=num(Salinity_PSU),Replicate=num(Replicate),Status_code=coalesce(parse_status(Status_raw),num(Status_code))) %>% filter(!is.na(Status_code),Status_code %in% 0:7)
 conflicts<-mal %>% group_by(across(all_of(obs_keys))) %>% summarise(nstatus=n_distinct(Status_code),.groups='drop') %>% filter(nstatus>1) %>% dplyr::select(all_of(obs_keys))
 daily<-mal %>% anti_join(conflicts,by=obs_keys) %>% group_by(across(all_of(obs_keys))) %>% arrange(Source_snapshot_row,.by_group=TRUE) %>% slice(1) %>% ungroup() %>% mutate(Generation_numeric=num(sub('^[Ff]','',trimws(Generation))),analysis_line_id=key(pick(everything()),line_keys),analysis_plate_id=key(pick(everything()),plate_keys),Is_egg=as.integer(Status_code==2),Is_juvenile=as.integer(Status_code==3),Is_dormant=as.integer(Status_code==4)) %>% arrange(across(all_of(obs_keys)))
 endpoints<-daily %>% group_by(across(all_of(line_keys))) %>% summarise(ever_offspring=as.integer(any(Status_code %in% c(2,3))),ever_dormant=as.integer(any(Status_code==4)),.groups='drop') %>% mutate(CE=sub('MAL','CE',Experiment_id),parent=floor(Replicate))
 lines<-endpoints %>% group_by(CE,Salinity_PSU,parent) %>% summarise(ever_offspring=max(ever_offspring),ever_dormant=max(ever_dormant),.groups='drop')
 counts<-lines %>% group_by(CE,PSU=Salinity_PSU) %>% summarise(n=n(),offspring=sum(ever_offspring),cryptobiosis=sum(ever_dormant),.groups='drop') %>% mutate(proportion=offspring/n)
 stopifnot(all(counts$n[counts$CE %in% c('CE2','CE3')]==30))
 write_out(lines,'line_endpoints_merged.csv');write_out(counts,'observed_counts.csv')
 wilson<-function(k,n){p<-k/n;z<-qnorm(.975);den<-1+z^2/n;cbind(lower=(p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/den,upper=(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/den)}
 crypto<-cbind(counts,cryptobiosis_percent=100*counts$cryptobiosis/counts$n,100*wilson(counts$cryptobiosis,counts$n));write_out(crypto,'cryptobiosis_points_Wilson95.csv');write_out(cbind(counts,wilson(counts$offspring,counts$n)),'offspring_points_Wilson95.csv')
 tests<-preds<-mods<-list()
 for(g in unique(counts$CE)){
  z<-subset(counts,CE==g);z$s<-z$PSU/10
  b<-glm(cbind(offspring,n-offspring)~s,data=z,family=binomial());q<-update(b,family=quasibinomial());m<-if(g=='CE2')q else b
  tt<-anova(update(m,.~1),m,test=if(g=='CE2')'F' else 'Chisq')
  tests[[g]]<-data.frame(CE=g,model=if(g=='CE2')'quasibinomial' else 'binomial',dispersion=summary(q)$dispersion,statistic=if(g=='CE2')tt$F[2] else tt$Deviance[2],df=tt$Df[2],p=if(g=='CE2')tt$`Pr(>F)`[2] else tt$`Pr(>Chi)`[2],OR_per10=exp(coef(m)[2]))
  gr<-data.frame(s=seq(0,max(z$s),length.out=151));pp<-predict(m,gr,type='link',se.fit=TRUE)
  preds[[g]]<-data.frame(CE=g,PSU=gr$s*10,fit=plogis(pp$fit),lower=plogis(pp$fit-1.96*pp$se.fit),upper=plogis(pp$fit+1.96*pp$se.fit));mods[[g]]<-m
 }
 occurrence_tests<-do.call(rbind,tests);occurrence_tests$p_Holm<-p.adjust(occurrence_tests$p,'holm');pd<-do.call(rbind,preds)
 write_out(occurrence_tests,'occurrence_tests.csv');write_out(pd,'occurrence_predictions.csv');saveRDS(mods,file.path(cache,'occurrence_models.rds'))
 parse_date<-function(v){
  if(is.na(v)||!nzchar(v)||v=='NA')return(as.Date(NA));n<-num(v)
  if(is.finite(n)&&n>40000&&n<50000)return(as.Date(n,origin='1899-12-30'))
  for(fmt in c('%d.%m.%Y','%Y-%m-%d')){a<-suppressWarnings(as.Date(v,format=fmt));if(!is.na(a))return(a)};as.Date(NA)
 }
 dates<-list();origins<-list();k<-0
 books<-list.files(file.path(project_root,'Data/raw/cultivation'),pattern='_plate_records[.]xlsx$',recursive=TRUE,full.names=TRUE)
 stopifnot(length(books)==3)
 for(book in books){exp<-sub('^CE','MAL',sub('_.*$','',basename(book)))
  for(sh in excel_sheets(book)){
   cells<-suppressMessages(read_excel(book,sheet=sh,col_names=FALSE,col_types='text',range=cell_limits(c(1,1),c(NA,NA)),.name_repair='minimal'))
   if(sh=='Days'&&ncol(cells)>=2){for(i in seq_len(nrow(cells))){if(!is.na(cells[[2]][i])&&cells[[2]][i]=='0'){a<-parse_date(cells[[1]][i]);if(!is.na(a))origins[[exp]]<-a}}}
   hits<-which(as.matrix(cells) %in% c('Date of Plate','Date_Plate'));if(!length(hits))next
   rows_hit<-row(as.matrix(cells))[hits];hr<-min(rows_hit);header<-as.character(cells[hr,])
   find_col<-function(v){j<-which(header %in% v);if(length(j))j[1] else NA_integer_}
   sc<-find_col(c('PSU','Treatment(PSU)'));rc<-find_col('Replicate');gc<-find_col('Generation');dc<-find_col(c('Date of Plate','Date_Plate'))
   if(anyNA(c(sc,rc,gc,dc))||hr>=nrow(cells))next
   for(i in (hr+1):nrow(cells)){dt<-parse_date(cells[[dc]][i]);if(is.na(dt))next
    k<-k+1;dates[[k]]<-data.frame(Experiment_id=exp,Salinity_PSU=num(cells[[sc]][i]),Replicate=num(cells[[rc]][i]),Generation=cells[[gc]][i],plate_date=as.character(dt),source=paste0(basename(book),':',sh,'!',readxl_col_letter(dc),i))
   }
  }
 }
 date_entries<-if(length(dates))do.call(rbind,dates) else stop('No plate dates found')
 write_out(date_entries,'plate_date_sources.csv')
 date_summary<-date_entries %>% group_by(Experiment_id,Salinity_PSU,Replicate,Generation) %>% summarise(plate_dates=paste(sort(unique(plate_date)),collapse='|'),date_sources=paste(source,collapse='|'),date_conflict=n_distinct(plate_date)>1,.groups='drop')
 generations<-daily %>% group_by(across(all_of(c(plate_keys,'analysis_line_id','analysis_plate_id','Generation_numeric')))) %>% summarise(generation_start_day=min(Day),first_offspring_day=if(any(Status_code %in% c(2,3)))min(Day[Status_code %in% c(2,3)]) else NA_real_,last_observation_day=max(Day),offspring_recorded=as.integer(any(Status_code %in% c(2,3))),.groups='drop') %>% arrange(across(all_of(plate_keys)))
 generations<-left_join(generations,date_summary,by=c('Experiment_id','Salinity_PSU','Replicate','Generation'));generations$date_conflict[is.na(generations$date_conflict)]<-FALSE
 generations$documented_start_day<-vapply(seq_len(nrow(generations)),function(i){z<-generations[i,];o<-origins[[z$Experiment_id]];if(!z$date_conflict&&!is.na(z$plate_dates)&&!is.null(o))as.numeric(as.Date(z$plate_dates)-o) else NA_real_},numeric(1))
 generations$start_used<-ifelse(is.finite(generations$documented_start_day),generations$documented_start_day,generations$generation_start_day)
 generations$time<-generations$first_offspring_day-generations$start_used
 generations$exclusion<-ifelse(generations$date_conflict,'conflicting plate dates',ifelse(generations$offspring_recorded==1 & generations$time<=0,'no documented start before offspring',''))
 write_out(generations,'generation_timing_audit.csv')
 d<-subset(generations,offspring_recorded==1 & exclusion=='');stopifnot(nrow(d)==611)
 d<-d[order(d$Experiment_id,d$Strain,d$Salinity_PSU,d$Replicate,d$Generation),]
 d$source_line_id<-d$analysis_line_id;d$analysis_line_id<-paste(d$Experiment_id,d$Strain,d$Salinity_PSU,floor(d$Replicate),sep='|')
 d$MAL<-factor(d$Experiment_id,levels=c('MAL1','MAL2','MAL3'));d$s<-d$Salinity_PSU/10;d$s2<-d$s^2;write_out(d,'analysis_data.csv')
 mm<-survreg(Surv(time,rep(1,nrow(d)))~MAL*s,data=d,dist='loglogistic',robust=TRUE,cluster=analysis_line_id);bc<-coef(mm);V<-vcov(mm)[names(bc),names(bc)]
 tt<-do.call(rbind,lapply(levels(d$MAL),function(g){L<-setNames(rep(0,length(bc)),names(bc));L['s']<-1;if(g!='MAL1')L[paste0('MAL',g,':s')]<-1;e<-sum(L*bc);se<-sqrt(drop(t(L)%*%V%*%L));data.frame(CE=sub('MAL','CE',g),ratio_per10=exp(e),lower=exp(e-1.96*se),upper=exp(e+1.96*se),p=2*pnorm(-abs(e/se)))}));tt$p_Holm<-p.adjust(tt$p,'holm');write_out(tt,'timing_tests.csv');saveRDS(mm,file.path(cache,'timing_model_robust.rds'))
 specs<-expand.grid(distribution=c('exponential','weibull','lognormal','loglogistic'),salinity_form=c('linear','quadratic'),stringsAsFactors=FALSE);candidates<-list();models<-list()
 for(i in seq_len(nrow(specs))){nm<-paste(specs$distribution[i],specs$salinity_form[i],sep='_');f<-if(specs$salinity_form[i]=='linear')Surv(time,rep(1,nrow(d)))~MAL*s else Surv(time,rep(1,nrow(d)))~MAL*(s+s2);m<-survreg(f,data=d,dist=specs$distribution[i]);models[[nm]]<-m;candidates[[nm]]<-data.frame(model=nm,AIC=AIC(m))}
 cmp<-do.call(rbind,candidates);cmp<-cmp[order(cmp$AIC),];cmp$delta_AIC<-cmp$AIC-min(cmp$AIC);write_out(cmp,'model_comparison.csv');fit<-models$loglogistic_linear
 grid<-do.call(rbind,lapply(levels(d$MAL),function(g)data.frame(MAL=factor(g,levels=levels(d$MAL)),Salinity_PSU=seq(0,max(d$Salinity_PSU[d$MAL==g]),length.out=151))))
 grid$s<-grid$Salinity_PSU/10;grid$s2<-grid$s^2;grid$median_days<-as.numeric(predict(fit,newdata=grid,type='quantile',p=.5))
 strata<-split(d,interaction(d$MAL,d$Salinity_PSU,drop=TRUE));boot<-matrix(NA_real_,bootstrap_iterations,nrow(grid));success<-logical(bootstrap_iterations)
 set.seed(20260915)
 for(b in seq_len(bootstrap_iterations)){
  if(b%%100==0)message('Cultivation bootstrap ',b,'/',bootstrap_iterations)
  chunks<-list();k<-1
  for(z in strata){ids<-unique(z$analysis_line_id);chosen<-ids[sample.int(length(ids),length(ids),replace=TRUE)];for(j in seq_along(chosen)){w<-z[z$analysis_line_id==chosen[j],];w$boot_id<-paste(b,k,j,sep='_');chunks[[k]]<-w;k<-k+1}}
  db<-do.call(rbind,chunks);mb<-try(survreg(Surv(time,rep(1,nrow(db)))~MAL*s,data=db,dist='loglogistic'),silent=TRUE)
  if(!inherits(mb,'try-error')){boot[b,]<-predict(mb,newdata=grid,type='quantile',p=.5);success[b]<-all(is.finite(boot[b,]))}
 }
 stopifnot(sum(success)>0);grid$lower_95<-apply(boot[success,,drop=FALSE],2,quantile,.025);grid$upper_95<-apply(boot[success,,drop=FALSE],2,quantile,.975)
 ratio<-grid;ratio$ratio_vs_0<-NA_real_;rb<-boot
 for(g in levels(d$MAL)){ii<-which(grid$MAL==g);j<-ii[which.min(grid$Salinity_PSU[ii])];ratio$ratio_vs_0[ii]<-grid$median_days[ii]/grid$median_days[j];rb[,ii]<-boot[,ii,drop=FALSE]/boot[,j]}
 ratio$ratio_lower_95<-apply(rb[success,,drop=FALSE],2,quantile,.025);ratio$ratio_upper_95<-apply(rb[success,,drop=FALSE],2,quantile,.975)
 write_out(ratio,'relative_prolongation_predictions.csv');write_out(data.frame(requested=bootstrap_iterations,successful=sum(success)),'bootstrap_summary.csv');saveRDS(list(model=fit,predictions=grid,bootstrap=boot,success=success),file.path(cache,'prediction_bootstrap.rds'))
 render_cultivation_figure(project_root,d,ratio,pd,counts,crypto,tt)
 invisible(list(occurrence=occurrence_tests,timing=tt))
}
readxl_col_letter<-function(n){s<-'';while(n>0){n<-n-1L;s<-paste0(LETTERS[n%%26L+1L],s);n<-n%/%26L};s}

render_cultivation_figure<-function(project_root,t,a,pd,d,crypto,effects){
 sys.source(file.path(project_root,'R/01_figure_style.R'),envir=environment())
 t$CE<-sub('MAL','CE',t$MAL);a$CE<-sub('MAL','CE',a$MAL)
 cols<-c(CE1='#4477AA',CE2='#228B75',CE3='#C47732');labels<-c('CE1','CE2','CE3')
 annotations<-data.frame(CE=c('CE1','CE2','CE3'),x=34.6,y=c(74,69,64));annotations$percent<-round(100*(effects$ratio_per10[match(annotations$CE,effects$CE)]-1));annotations$label<-sprintf('%s: +%d%% per +10 PSU',annotations$CE,annotations$percent)
 aplot<-ggplot()+geom_ribbon(data=pd,aes(PSU,ymin=100*lower,ymax=100*upper,fill=CE),alpha=.08,show.legend=FALSE)+geom_ribbon(data=a,aes(Salinity_PSU,ymin=lower_95,ymax=upper_95,fill=CE),alpha=.2,show.legend=FALSE)+geom_point(data=t,aes(Salinity_PSU,time,colour=CE),position=position_jitter(width=.3,height=0,seed=17),size=1.15,alpha=.26,show.legend=FALSE)+geom_line(data=a,aes(Salinity_PSU,median_days,colour=CE,linetype='Time to first offspring'),linewidth=1.1,show.legend=FALSE)+geom_line(data=pd,aes(PSU,100*fit,colour=CE,linetype='Offspring occurrence'),linewidth=1,show.legend=FALSE)+geom_point(data=d,aes(PSU,100*proportion,colour=CE),shape=21,fill='white',size=3,stroke=1,show.legend=FALSE)+geom_text(data=annotations,aes(x=x,y=y,label=label,colour=CE),hjust=1,size=3.5,fontface='plain',show.legend=FALSE)+geom_point(data=data.frame(CE=c('CE1','CE2','CE3'),x=NA_real_,y=NA_real_),aes(x,y,colour=CE),shape=21,fill='white',size=3,stroke=1,na.rm=TRUE,show.legend=TRUE)+geom_line(data=data.frame(x=NA_real_,y=NA_real_,kind=c('Time to first offspring','Offspring occurrence')),aes(x,y,linetype=kind),colour='black',linewidth=.9,na.rm=TRUE,show.legend=TRUE)+scale_colour_manual(values=cols,labels=labels)+scale_fill_manual(values=cols,guide='none')+scale_linetype_manual(values=c('Time to first offspring'='solid','Offspring occurrence'='22'))+scale_y_continuous(name='Offspring occurrence (%)',limits=c(0,100),breaks=seq(0,100,20),sec.axis=sec_axis(~.,name='Time to first offspring (days)',breaks=seq(0,100,20)))+scale_x_continuous(limits=c(-.4,35.4),breaks=seq(0,35,5),expand=expansion(mult=0))+labs(x='Salinity (PSU)',colour=NULL,linetype=NULL,tag='A')+guides(colour=guide_legend(order=1,nrow=1,byrow=TRUE,override.aes=list(shape=21,fill='white',size=3,linetype=0)),linetype=guide_legend(order=2,override.aes=list(colour='black',linewidth=1,shape=NA)))+theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position='bottom',legend.box='vertical',plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=9),legend.key.width=grid::unit(.45,'cm'),legend.key.height=grid::unit(.55,'cm'),legend.spacing.y=grid::unit(.12,'cm'),legend.margin=margin(5,2,2,2),axis.text.x=element_text(size=11,colour='black'),axis.ticks.x=element_line(colour='black',linewidth=.7),axis.ticks.length=grid::unit(2.5,'mm'),axis.text.y=element_text(size=11,colour='black'),axis.title.y=element_text(size=14,face='plain'),axis.title.y.right=element_text(size=14,face='plain'),plot.tag=element_text(face='bold',size=18),legend.text=ggtext::element_markdown(size=11,lineheight=1.1,margin=margin(l=1)),legend.spacing.x=grid::unit(2,'mm'),plot.margin=margin(12,14,8,12))
 key_labels <- c('Time to first offspring', 'Offspring occurrence')
 key_y <- c(.72, .28)
 key_grobs <- lapply(seq_along(key_labels), function(i) {
  txt <- grid::textGrob(key_labels[i], x=grid::unit(1,'npc'),
   y=grid::unit(key_y[i],'npc'), just='right',
   gp=grid::gpar(fontsize=11))
  key_end <- grid::unit(1,'npc') - grid::grobWidth(txt) - grid::unit(2,'mm')
  grid::grobTree(
   grid::segmentsGrob(x0=key_end-grid::unit(11,'mm'), x1=key_end,
    y0=grid::unit(key_y[i],'npc'), y1=grid::unit(key_y[i],'npc'),
    gp=grid::gpar(lwd=.9*72.27/25.4,lty=if(i==1) 'solid' else '22')),
   txt)
 })
 line_key <- patchwork::wrap_elements(full=do.call(grid::grobTree,key_grobs))
 aplot <- aplot + guides(linetype='none') + theme(axis.title.x=element_text(size=14,face='plain',margin=margin(t=5))) +
  inset_element(line_key, left=.58, bottom=.81, right=(34.6+.4)/35.8, top=.94, align_to='panel')
 bplot<-ggplot(crypto,aes(PSU,cryptobiosis_percent,colour=CE,group=CE))+geom_errorbar(aes(ymin=lower,ymax=upper),width=.35,linewidth=.75,position=position_dodge(width=.65))+geom_point(shape=21,fill='white',stroke=1.1,size=3.4,position=position_dodge(width=.65))+scale_colour_manual(values=cols,labels=labels)+scale_x_continuous(limits=c(-.4,35.4),breaks=seq(0,35,5),expand=expansion(mult=0))+scale_y_continuous(limits=c(0,40),breaks=seq(0,40,10),expand=expansion(mult=c(.02,.04)))+labs(x='Salinity (PSU)',y='Lineages with\ncryptobiosis (%)',colour=NULL,tag='B')+guides(colour='none')+theme_bw(base_size=16)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),plot.tag=element_text(face='bold',size=18),axis.title=element_text(size=14,face='plain'),axis.title.x=element_text(margin=margin(t=5)),axis.text=element_text(size=11,colour='black'),axis.ticks.length=grid::unit(2.5,'mm'),axis.ticks=element_line(linewidth=.7),plot.margin=margin(10,14,8,12))
 aplot<-aplot & theme(axis.title.x=element_text(size=14,face='plain',hjust=.5,margin=margin(t=5)),axis.title.y=element_text(size=14,face='plain',hjust=.5,vjust=1,margin=margin(r=6)),axis.title.y.right=element_text(size=14,face='plain',hjust=.5,margin=margin(l=6)),axis.text.y.left=element_text(margin=margin(r=2.2)),axis.text.y.right=element_text(margin=margin(l=2.2)))
 bplot<-bplot+theme(axis.title.x=element_text(size=14,face='plain',hjust=.5,margin=margin(t=5)),axis.title.y=element_text(size=14,face='plain',hjust=.5,margin=margin(r=6)))
 fig<-(aplot/bplot)+plot_layout(heights=c(2,1),guides='collect') & theme(legend.position='bottom',legend.box='vertical')
 out<-file.path(project_root,'Manuscript/figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
 figure_temp<-file.path(project_root,'.tmp','Figure_3_cultivation_render.png')
 figure_grob<-patchwork::patchworkGrob(fig)
 shift_left_title<-function(g){
  if(inherits(g,'gtable')){
   for(i in seq_along(g$grobs)){
    if(g$layout$name[i]=='ylab-l-1'){
     title_grob<-g$grobs[[i]]
     for(j in seq_along(title_grob$children)) if(inherits(title_grob$children[[j]],'text')) title_grob$children[[j]]$x<-title_grob$children[[j]]$x+grid::unit(12,'pt')
     g$grobs[[i]]<-title_grob
    }else g$grobs[[i]]<-shift_left_title(g$grobs[[i]])
   }
  }
  g
 }
 figure_grob<-shift_left_title(figure_grob)
 ggsave(figure_temp,figure_grob,width=9,height=8.2,dpi=300,bg='white')
 if(!file.copy(figure_temp,file.path(out,'Figure_3_cultivation.png'),overwrite=TRUE)) stop('Could not replace cultivation figure')
 saveRDS(fig,file.path(project_root,'Results/cache/cultivation/figure.rds'))
}
