context("Mapping file")

minimal_map_lines <- c(
  "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription",
  "A.1\tAACCGGTT\tACGT\tTest sample")
minimal_map_df <- data.frame(
  BarcodeSequence="AACCGGTT",
  LinkerPrimerSequence="ACGT",
  Description="Test sample")
rownames(minimal_map_df) <- "A.1"

test_that("Minimal QIIME mapping file is loaded correctly", {
  map_fp <- tempfile()
  writeLines(minimal_map_lines, map_fp)

  expect_equal(load.qiime.mapping.file(map_fp), minimal_map_df)

  unlink(map_fp)
})