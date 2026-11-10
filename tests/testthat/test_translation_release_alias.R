test_that('Translation sessions pin current before lazy reads', {
  root<-tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
  old<-file.path(root,'v1');new<-file.path(root,'v2');dir.create(old);dir.create(new)
  for(path in c(old,new)) {
    saveRDS(list(),file.path(path,'atlas.rds'))
    saveRDS(list(species='Homo sapiens',marker=basename(path)),file.path(path,'dataset.rds'))
  }
  current<-file.path(root,'current');expect_true(file.symlink(old,current))
  shiny::testServer(translation_server,args=list(browser_options=c(translation_data_dir=current)),{
    session$setInputs(species='Homo sapiens')
    expect_identical(base(),normalizePath(old))
    unlink(current);expect_true(file.symlink(new,current));session$flushReact()
    expect_identical(translation_cached_read(file.path(base(),'dataset.rds'))$marker,'v1')
  })
  shiny::testServer(translation_server,args=list(browser_options=c(translation_data_dir=current)),{
    session$setInputs(species='Homo sapiens')
    expect_identical(base(),normalizePath(new))
    expect_identical(translation_cached_read(file.path(base(),'dataset.rds'))$marker,'v2')
  })
})
