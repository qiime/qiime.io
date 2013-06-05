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


