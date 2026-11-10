#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
file_arg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())])
repo<-normalizePath(file.path(dirname(file_arg),'../..'))
suppressMessages(devtools::load_all(repo,quiet=TRUE))
directory<-if(length(args))args[1] else Sys.getenv('RIBOCRYPT_TRANSLATION_DIR','')
if(!nzchar(directory))stop('Supply a prepared Translation dataset directory.')
port<-if(length(args)>1L)as.integer(args[2]) else 9284L
shiny::runApp(ribocrypt_translation_app(directory),host='127.0.0.1',port=port,launch.browser=FALSE)
