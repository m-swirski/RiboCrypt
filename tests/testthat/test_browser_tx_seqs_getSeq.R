df_browser_tx_seqs <- ORFik::ORFik.template.experiment()[9:10, ]
tx_browser_tx_seqs <- ORFik::loadRegion(df_browser_tx_seqs, "mrna")
fa_browser_tx_seqs <- ORFik::findFa(df_browser_tx_seqs)

test_that("browser_tx_seqs_getSeq matches extractTranscriptSeqs for single-exon plus transcripts", {
  tx <- tx_browser_tx_seqs[1]
  tx[[1]] <- tx[[1]][1]

  expected <- GenomicFeatures::extractTranscriptSeqs(fa_browser_tx_seqs, tx)
  observed <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx)
  observed_no_names <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx, keep.names = FALSE)

  expect_identical(observed, expected)
  expect_null(names(observed_no_names))
  expect_identical(as.character(observed_no_names), unname(as.character(expected)))
})

test_that("browser_tx_seqs_getSeq matches extractTranscriptSeqs for single-exon minus transcripts", {
  tx <- tx_browser_tx_seqs[5]
  tx[[1]] <- tx[[1]][1]

  expected <- GenomicFeatures::extractTranscriptSeqs(fa_browser_tx_seqs, tx)
  observed <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx)
  observed_no_names <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx, keep.names = FALSE)

  expect_identical(observed, expected)
  expect_null(names(observed_no_names))
  expect_identical(as.character(observed_no_names), unname(as.character(expected)))
})

test_that("browser_tx_seqs_getSeq matches extractTranscriptSeqs for multi-exon plus transcripts", {
  tx <- tx_browser_tx_seqs[1]

  expected <- GenomicFeatures::extractTranscriptSeqs(fa_browser_tx_seqs, tx)
  observed <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx)
  observed_no_names <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx, keep.names = FALSE)

  expect_identical(observed, expected)
  expect_null(names(observed_no_names))
  expect_identical(as.character(observed_no_names), unname(as.character(expected)))
})

test_that("browser_tx_seqs_getSeq matches extractTranscriptSeqs for multi-exon minus transcripts", {
  tx <- tx_browser_tx_seqs[5]

  expected <- GenomicFeatures::extractTranscriptSeqs(fa_browser_tx_seqs, tx)
  observed <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx)
  observed_no_names <- browser_tx_seqs_getSeq(fa_browser_tx_seqs, tx, keep.names = FALSE)

  expect_identical(observed, expected)
  expect_null(names(observed_no_names))
  expect_identical(as.character(observed_no_names), unname(as.character(expected)))
})
