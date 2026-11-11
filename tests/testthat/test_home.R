test_that("home groups modalities, merged tracks, collections and studies", {
  n <- c("all_merged-Homo_sapiens_modalities", "all_merged-Homo_sapiens_RNA",
         "all_samples-Homo_sapiens_RNA", "PRJ123-homo_sapiens")
  expect_equal(rc_home_groups(n),c("Modalities","Merged tracks","Collections","Single experiments"))
})

test_that("metadata search is literal, case insensitive and preserves Run IDs", {
  m <- data.table::data.table(Run=c("SRR1","SRR2","SRR3"), GENE=c("EIF4A2","eIf4G2",NA_character_), sample_title=c("control","treated","eIF4 control"))
  expect_equal(rc_home_search(m,"eif4"),c("SRR1","SRR2","SRR3"))
  expect_equal(rc_home_search(m,"EIF4A2"),"SRR1")
  expect_length(rc_home_search(m,"."),0)
  expect_length(rc_home_search(m,"not found"),0)
})

test_that("home catalog is limited to the supplied experiment catalog", {
  e <- data.table::data.table(name=c("a","b"),organism=c("Homo sapiens","Homo sapiens x virus"))
  m <- data.table::data.table(ScientificName=c("Homo sapiens","Mus musculus"),Run=c("A","B"),study_accession=c("X","Y"))
  a <- rc_home_catalog(e,m)
  expect_equal(a$name,"Homo sapiens")
  expect_equal(a$libraries,1L)
  expect_equal(a$studies,1L)
})

test_that("logo navigation is enabled only with a home page", {
  expect_match(as.character(rc_title(TRUE)),"#home",fixed=TRUE)
  expect_false(grepl("#home",as.character(rc_title(FALSE)),fixed=TRUE))
})

test_that("every species clade is a wheel clade and viral references are Viruses", {
  x <- utils::read.csv(system.file("home", "species_clades.csv", package = "RiboCrypt"), stringsAsFactors = FALSE)
  wheel <- c("Mammals","Other vertebrates","Invertebrates","Plants","Fungi","Protists","Bacteria","Archaea","Viruses")
  expect_true(all(x$clade %in% wheel))
  expect_false(anyDuplicated(x$species) > 0)
  expect_equal(x$clade[x$species %in% c("Sars cov2", "Sars rozman2016")], c("Viruses", "Viruses"))
})
