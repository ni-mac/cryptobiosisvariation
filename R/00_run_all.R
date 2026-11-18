# Entry point. Run from any working directory: Rscript R/00_run_all.R
# No package installation, downloads or files outside the project are requested.
run_all <- function(project_root, modules=c('salinity','cultivation','desiccation','supplement_figures'), bootstrap_mode='bootstrap_reference_replay', cultivation_bootstraps=1000L) {
 root<-normalizePath(project_root,winslash='/',mustWork=TRUE)
 old_wd<-getwd();setwd(root);on.exit(setwd(old_wd),add=TRUE)
 needed<-c('readxl','dplyr','tidyr','stringr','purrr','survival','ggplot2','patchwork','glmmTMB','MASS','readr','scales','cowplot','ggtext','emmeans','DHARMa','png','jsonlite','digest')
 missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
 if(length(missing))stop('Install required R packages before running: ',paste(missing,collapse=', '))
 dir.create(file.path(root,'.tmp'),showWarnings=FALSE)
 dir.create(file.path(root,'Results'),showWarnings=FALSE)
 dir.create(file.path(root,'Results','run_info'),recursive=TRUE,showWarnings=FALSE)
 started<-Sys.time();timings<-list()
 sources<-c(salinity='02_salinity.R',cultivation='03_cultivation.R',desiccation='04_desiccation.R',supplement_figures='05_supplement_figures.R')
 for(nm in modules){
  if(!nm %in% names(sources))stop('Unknown module: ',nm)
  cat('\nRunning ',nm,'\n',sep='');t0<-Sys.time()
  e<-new.env(parent=globalenv());sys.source(file.path(root,'R',sources[[nm]]),envir=e)
  fun<-e[[paste0('run_',nm)]]
  if(nm=='salinity')fun(root,bootstrap_mode=bootstrap_mode) else if(nm=='cultivation')fun(root,bootstrap_iterations=cultivation_bootstraps) else fun(root)
  timings[[nm]]<-data.frame(module=nm,seconds=as.numeric(difftime(Sys.time(),t0,units='secs')))
 }
 write.csv(do.call(rbind,timings),file.path(root,'Results/run_info/run_times.csv'),row.names=FALSE)
 versions<-data.frame(package=needed,version=vapply(needed,function(x)as.character(packageVersion(x)),character(1)))
 write.csv(versions,file.path(root,'Results/run_info/package_versions.csv'),row.names=FALSE)
 capture.output(sessionInfo(),file=file.path(root,'Results/run_info/sessionInfo.txt'))
 writeLines(c(paste('Started:',started),paste('Finished:',Sys.time()),paste('Salinity bootstrap mode:',bootstrap_mode),paste('Cultivation bootstrap iterations:',cultivation_bootstraps)),file.path(root,'Results/run_info/last_run.txt'))
 invisible(timings)
}
args<-commandArgs(FALSE);script<-sub('^--file=','',grep('^--file=',args,value=TRUE)[1])
if(!is.na(script)&&basename(script)=='00_run_all.R'){
 root<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
 extra<-commandArgs(trailingOnly=TRUE)
 modules<-if(length(extra))strsplit(extra[1],',',fixed=TRUE)[[1]] else c('salinity','cultivation','desiccation','supplement_figures')
 run_all(root,modules=modules)
}
