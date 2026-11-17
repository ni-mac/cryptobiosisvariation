# Supplemental visuals drawn from documented protocol and analysed recovery counts.
run_supplement_figures <- function(project_root) {
 suppressPackageStartupMessages({library(grid);library(ggplot2);library(patchwork);library(dplyr);library(cowplot)})
 root<-normalizePath(project_root,winslash='/',mustWork=TRUE)
 sys.source(file.path(root,'R/01_figure_style.R'),envir=environment())
 out<-file.path(root,'Manuscript/Supplement/Supplement_figures')
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 # S1 is a table only; no flow diagram and no D3.
 png(file.path(out,'Figure_S1_desiccation_design.png'),width=2700,height=810,res=300,bg='white')
 grid.newpage()
 widths<-c(.19,.19,.23,.19,.20);edges<-c(0,cumsum(widths))
 headers<-c('Experiment','Preconditioning\ncomparison','Desiccation\ncomparison','Desiccation\nconditions','Rehydration')
 rows<-list(c('Desiccation\nexperiment 1 (D1)','None','Non-desiccated control\nvs desiccated treatment','8 h\n<10% RH','M9 buffer\n17 h'),
            c('Desiccation\nexperiment 2 (D2)','0 h vs 96 h\n97–98% RH','Both groups\ndesiccated','48 h\n<10% RH','dH2O, 24 h\nshaking at 50 rpm'))
 fills<-c('#eeeeee','#fff2cc','#f4cccc','#cfe8f3','#d9d2e9')
 for(j in 1:5)for(i in 0:2){
  yc<-if(i==0).83 else if(i==1).515 else .185;hh<-if(i==0).28 else .33
  grid.rect(x=edges[j]+widths[j]/2,y=yc,width=widths[j],height=hh,gp=gpar(fill=if(i==0)fills[j] else adjustcolor(fills[j],alpha.f=.25),col='#888',lwd=.7))
  grid.text(if(i==0)headers[j] else rows[[i]][j],x=edges[j]+widths[j]/2,y=yc,
    gp=gpar(fontsize=if(i==0)11.5 else 10.5,fontface=if(i==0)'bold' else 'plain',fontfamily='sans'))
 }
 dev.off()
 # S2: preserve recovery counts, pooled denominators, individual points and crosses.
 st<-file.path(root,'Results/statistics/salinity')
 points<-read.csv(file.path(st,'m9_recovery_analysis_ready.csv'))
 pooled<-read.csv(file.path(st,'m9_recovery_by_salinity.csv'))
 tested<-read.csv(file.path(st,'m9_recovery_tested_salinity_grid.csv'))
 labels<-c('ALT.22.08'='bolditalic(Panagrolaimus)~bold("sp.")~plain("(WWL002)")',
           'PAP.22.17'='bolditalic(Acrobeloides~tricornis)~plain("(WWL709)")',
           'A. max'='bolditalic(Acrobeloides~maximus)~plain("(A. max)")')
 panel<-function(id,strain){
  p<-filter(points,Experiment_id==id);q<-filter(pooled,Experiment_id==id)
  missing<-filter(tested,Experiment_id==id,!has_cryptobiotic_candidates)
  q<-mutate(q,label=ifelse(recovered_animals>0,paste0(recovered_animals,'/',cryptobiotic_candidates),''),
            label_y=ifelse(recovered_fraction_candidate>.95,100*recovered_fraction_candidate-9,100*recovered_fraction_candidate+8))
  ggplot()+
   geom_point(data=p,aes(salinity_PSU,100*recovered_fraction_candidate),colour='#999999',alpha=.5,size=1.7,position=position_jitter(width=.15,height=0,seed=1211))+
   geom_point(data=q,aes(salinity_PSU,100*recovered_fraction_candidate),shape=21,fill=figure_style$blue,colour='#333333',size=2.8,stroke=.65)+
   geom_point(data=missing,aes(salinity_PSU,-1.5),shape=4,colour='#666666',size=2.8,stroke=.9)+
   geom_text(data=q,aes(salinity_PSU,label_y,label=label),size=3.4,colour='#333333')+
   facet_wrap(~factor(labels[[strain]]),labeller=label_parsed)+
   scale_x_continuous(breaks=seq(0,40,10))+
   scale_y_continuous(breaks=seq(0,100,25))+
   coord_cartesian(xlim=c(0,40),ylim=c(-4,104))+
   labs(x=NULL,y=NULL,title=if(grepl('ALT',id)) paste0(unique(p$temperature_C),' °C') else NULL)+figure_theme()+theme(plot.margin=margin(4,5,4,5),plot.title=element_text(hjust=.5,size=16,face='bold'),plot.title.position='panel')
 }
 core<-((panel('ALT.22.08_22C','ALT.22.08')+panel('ALT.22.08_16C','ALT.22.08'))/
          (panel('A.max_22C','A. max')+panel('PAP.22.17_16C','PAP.22.17')))+plot_layout(heights=c(1,1))
 plot<-ggdraw()+draw_plot(core,x=.065,y=.065,width=.935,height=.935)+
  draw_label('Salinity (PSU)',x=.55,y=.060,size=14)+
  draw_label('Recovery (%)',x=.060,y=.54,angle=90,size=14)
 recovery_temp<-file.path(root,'.tmp','Figure_S2_recovery_render.png')
 ggsave(recovery_temp,plot,width=9,height=7.5,dpi=300,bg='white')
 if(!file.copy(recovery_temp,file.path(out,'Figure_S2_recovery.png'),overwrite=TRUE)) warning('Rendered figure must be copied to final destination')
 # Explicit plotting inputs enable checking that only style, not recovery results, changed.
 dir.create(file.path(root,'Results/cache/supplement'),recursive=TRUE,showWarnings=FALSE)
 saveRDS(list(points=points,pooled=pooled,tested=tested),file.path(root,'Results/cache/supplement/recovery_plot_inputs.rds'))
 invisible(out)
}
