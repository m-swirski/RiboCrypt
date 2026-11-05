test_that("compact sequence payload reconstructs the original six traces", {
  for (sequence in c(paste(rep("ATGTAA", 23), collapse = ""),
                     paste(rep("ATGNNNTAG", 15), collapse = ""))) {
    for (style in c("one_letter", "three_letters")) {
      dna <- Biostrings::DNAString(sequence)
      payload <- RiboCrypt:::fetch_JS_seq(dna, nplots = 1, distance = 250,
                                           aa_letter_code = style, input_id = "x")
      letters <- strsplit(sequence, "", fixed = TRUE)[[1]]
      ir <- IRanges::IRanges(start = seq.int(3),
                              width = length(letters) - seq.int(0, 2))
      aas <- suppressWarnings(Biostrings::translate(Biostrings::extractAt(dna, ir),
        if.fuzzy.codon = if (grepl("[^ACGT]", sequence)) "X" else "error"))
      aas <- lapply(as.character(aas), function(x) strsplit(x, "", fixed = TRUE)[[1]])
      for (frame in seq_len(3)) {
        nt <- payload$traces[[frame]]
        positions <- seq(frame, length(letters), 3)
        expect_equal(nt$length, length(positions))
        expect_identical(letters[positions], substring(payload$sequence, positions, positions))
        aa <- payload$traces[[frame + 3]]
        expect_equal(aa$length, length(aas[[frame]]))
        observed <- strsplit(payload$translations[[frame]], "", fixed = TRUE)[[1]]
        if (style == "three_letters") {
          observed <- unlist(payload$aa_codes[observed], use.names = FALSE)
          codes <- Biostrings::AMINO_ACID_CODE; codes["*"] <- "*"
          expected <- unname(codes[match(aas[[frame]], names(codes))])
        } else expected <- aas[[frame]]
        expect_identical(unname(observed), unname(expected))
      }
    }
  }
})
