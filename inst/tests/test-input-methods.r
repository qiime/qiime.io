context("OTU table")

test_that("Minimal OTU table is loaded correctly", {
  minimal_otu_lines <- c(
    "# Constructed from biom file",
    "#OTU ID\tA.1\tB.1\tConsensus Lineage",
    "0\t2.0\t1.0\tBacteria; Bacteroidetes",
    "1\t0.0\t0.0\tBacteria",
    "2\t0.0\t0.0\tBacteria; Actinobacteria")
  otu_fp <- tempfile()
  writeLines(minimal_otu_lines, otu_fp)
  
  minimal_otus <- matrix(c(2.0, 1.0, 0.0, 0.0, 0.0, 0.0), nrow=2)
  rownames(minimal_otus) <- c("A.1", "B.1")
  colnames(minimal_otus) <- c("0", "1", "2")
  minimal_lineage <- c(
    "Bacteria; Bacteroidetes", "Bacteria", "Bacteria; Actinobacteria")
  
  obs_nolineage <- load.qiime.otu.table(otu_fp)
  expect_equal(obs_nolineage, minimal_otus)
  
  obs_lineage <- load.qiime.otu.table(otu_fp, include.lineages=T)
  expect_equal(obs_lineage, list(otus=minimal_otus, lineages=minimal_lineage))
  
  unlink(otu_fp)
})


context("Mapping file")

test_that("Minimal mapping file is loaded correctly", {
  minimal_map_lines <- c(
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription",
    "A.1\tAACCGGTT\tACGT\tTest sample")
  map_fp <- tempfile()
  writeLines(minimal_map_lines, map_fp)
  
  minimal_map <- data.frame(
    BarcodeSequence="AACCGGTT",
    LinkerPrimerSequence="ACGT",
    Description="Test sample")
  rownames(minimal_map) <- "A.1"
  
  expect_equal(load.qiime.mapping.file(map_fp), minimal_map)
  
  unlink(map_fp)
})

test_that("Whitespace is stripped from fields in mapping file", {
  map_fp <- tempfile()
  writeLines(c("#SampleID\tA\tB", "A.1 \t val1\tval2 "), map_fp)

  expected_df <- data.frame(A="val1", B="val2")
  rownames(expected_df) <- "A.1"

  expect_equal(load.qiime.mapping.file(map_fp), expected_df)

  unlink(map_fp)
})

test_that("Comments are ignored in mapping file", {
  map_fp <- tempfile()
  writeLines(c("#SampleID\tA\tB", "# Comment 1", "w\tx\ty", "#C2"), map_fp)

  expected_df <- data.frame(A="x", B="y")
  rownames(expected_df) <- "w"

  expect_equal(load.qiime.mapping.file(map_fp), expected_df)

  unlink(map_fp)
})

test_that("Double quotes are stripped from values in mapping file", {
  map_fp <- tempfile()
  writeLines(c("#SampleID\tA\tB", "w\t\"x\"\ty"), map_fp)

  expected_df <- data.frame(A="x", B="y")
  rownames(expected_df) <- "w"

  expect_equal(load.qiime.mapping.file(map_fp), expected_df)

  unlink(map_fp)
})

test_that("Empty mapping file produces error", {
  empty_fp <- tempfile()
  file.create(empty_fp)
  expect_error(load.qiime.mapping.file(empty_fp), "empty")
  unlink(empty_fp)
})

test_that("No leading # in header produces warning", {
  nohead_fp <- tempfile()
  writeLines(c("SampleID\tA\tB", "w\tx\ty"), nohead_fp)
  expect_warning(load.qiime.mapping.file(nohead_fp), "Header")
  unlink(nohead_fp)
})

test_that("No header line produces warning", {
  nohead_fp <- tempfile()
  writeLines(c("a\tb\tc", "w\tx\ty"), nohead_fp)
  expect_warning(load.qiime.mapping.file(nohead_fp), "Header")
  unlink(nohead_fp)
})

test_that("Header but no values produces error", {
  norows_fp <- tempfile()
  writeLines(c("#SampleID\tA\tB", "#Comment"), norows_fp)
  expect_error(load.qiime.mapping.file(norows_fp), "line")
  unlink(norows_fp)
})


context("Taxon table")

test_that("Minimal distance matrix is loaded correctly", {
  minimal_taxon_lines <- c(
    "Taxon\tA.1\tB.1",
    "Bacteria;Actinobacteria\t0.13\t0.03",
    "Bacteria;Bacteroidetes\t0.17\t0.22")
  taxon_fp <- tempfile()
  writeLines(minimal_taxon_lines, taxon_fp)
  
  taxon_table <- matrix(c(0.13, 0.03, 0.17, 0.22), nrow=2)
  rownames(taxon_table) <- c("A.1", "B.1")
  colnames(taxon_table) <- c(
    "Bacteria;Actinobacteria", "Bacteria;Bacteroidetes")
  
  expect_equal(load.qiime.taxon.table(taxon_fp), taxon_table)
  
  unlink(taxon_fp)
})


context("Distance matrix")

test_that("Minimal distance matrix is loaded correctly", {
  minimal_dm_lines <- c(
    "\tA.1\tB.1\tC.1",
    "A.1\t0.0\t0.8\t0.7",
    "B.1\t0.8\t0.0\t0.6",
    "C.1\t0.7\t0.6\t0.0")
  dm_fp <- tempfile()
  writeLines(minimal_dm_lines, dm_fp)
  
  dm <- matrix(c(0.0, 0.8, 0.7, 0.8, 0.0, 0.6, 0.7, 0.6, 0.0), nrow=3)
  rownames(dm) <- c("A.1", "B.1", "C.1")
  colnames(dm) <- c("A.1", "B.1", "C.1")

  expect_equal(load.qiime.distance.matrix(dm_fp), dm)
  
  unlink(dm_fp)
})


context("Non-overlapping samples")

test_that("Non-overlapping samples are properly removed", {
  map <- data.frame(
    BarcodeSequence="AACCGGTT",
    LinkerPrimerSequence="ACGT",
    Description="Test sample")
  rownames(map) <- "A.1"

  otus <- matrix(c(2.0, 1.0, 0.0, 0.0, 0.0, 0.0), nrow=2)
  rownames(otus) <- c("A.1", "B.1")
  
  taxa <- matrix(c(0.13, 0.03, 0.17, 0.22), nrow=2)
  rownames(taxa) <- c("A.1", "B.1")
  
  distmat <- matrix(c(0.0, 0.8, 0.7, 0.8, 0.0, 0.6, 0.7, 0.6, 0.0), nrow=3)
  rownames(distmat) <- c("A.1", "B.1", "C.1")
  colnames(distmat) <- c("A.1", "B.1", "C.1")
  
  observed_nomap <- remove.nonoverlapping.samples(
    otus=otus, taxa=taxa, distmat=distmat)
  expect_equal(observed_nomap$map, NULL)
  expect_equal(rownames(observed_nomap$otus), c("A.1", "B.1"))
  expect_equal(rownames(observed_nomap$taxa), c("A.1", "B.1"))
  expect_equal(rownames(observed_nomap$distmat), c("A.1", "B.1"))
  
  observed_map <- remove.nonoverlapping.samples(
    map=map, otus=otus, taxa=taxa, distmat=distmat)
  expect_equal(rownames(observed_map$map), "A.1")
  expect_equal(rownames(observed_map$otus), "A.1")
  expect_equal(rownames(observed_map$taxa), "A.1")
  expect_equal(rownames(observed_map$distmat), "A.1")
})
