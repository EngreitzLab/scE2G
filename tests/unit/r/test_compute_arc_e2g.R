suppressMessages({
  library(data.table)
})

## compute_arc_e2g.R is a Snakemake `script:` that reads its inputs from the
## `snakemake` S4 object, so run the real script end to end on small synthetic
## inputs with a stand-in for that object.
## testthat::test_dir() runs each test file with its own directory as the
## working directory, so anchor back up to the repo root.
script_path <- normalizePath(file.path("..", "..", "..", "workflow", "scripts",
                                       "feature_computation", "compute_arc_e2g.R"))
setClass("MockSnakemake", representation(input = "list", output = "list", params = "list"))

## One row per gene. tpm = gene-level RNA metric in Pairs.Kendall (NA: gene not
## in Pairs.Kendall). kendall = "ok" (all pairs non-NA), "allNA", "mixed".
genes <- data.table(
  TargetGene = c("OXTR", "OXT",   "PCDHB10", "PCDHB11", "PCDHB1", "ZZZ1",
                 "ARC",  "ARCN1", "MDH1B",   "MDH1",    "MIX",    "GENE1", "GENE2", "GENE3"),
  tpm        = c(3.894,  0,       2,         3,         0,        0,
                 5,      25,      1.1,       NA,        4,        10,      20,      30),
  kendall    = c("ok",   "allNA", "ok",      "ok",      "allNA",  "allNA",
                 "ok",   "ok",    "ok",      "absent",  "mixed",  "ok",    "ok",    "ok"))
n_pairs <- 3L
pairs <- genes[rep(seq_len(.N), each = n_pairs)]
pairs[, pair := rep(seq_len(n_pairs), nrow(genes))]
pairs[, start := 1000000L + 10000L * rep(seq_len(nrow(genes)), each = n_pairs) + 1000L * pair]
pairs[, end := start + 500L]
pairs[, ABC.Score := round(0.01 * seq_len(.N) + 0.005 * pair, 4)]
pairs[, Kendall := fcase(kendall == "ok", round(0.02 * seq_len(.N) - 0.1, 4),
                         kendall == "mixed", c(NA, 0.15, NA)[pair],
                         default = NA_real_)]

abc <- pairs[, .(chr = "chr1", start, end, name = paste0("e", seq_len(.N)),
                 class = "intergenic", TargetGene, ABC.Score)]
kendall <- pairs[kendall != "absent",
                 .(chr = "chr1", start, end, TargetGene,
                   PeakName = paste("chr1", start, end, sep = "-"),
                   PairName = paste0("chr1-", start, "-", end, "_", TargetGene),
                   mean_log_normalized_rna = tpm / 100,
                   RnaDetectedPercent = pmin(1, tpm / 30),
                   RnaPseudobulkTPM = tpm,
                   Kendall)]

run_dir <- tempfile("compute_arc_e2g_test_")
dir.create(run_dir)
abc_path <- file.path(run_dir, "EnhancerPredictionsAllPutative.tsv.gz")
kendall_path <- file.path(run_dir, "Pairs.Kendall.tsv.gz")
arc_path <- file.path(run_dir, "EnhancerPredictionsAllPutative_ARC.tsv.gz")
fwrite(abc, abc_path, sep = "\t")
fwrite(kendall, kendall_path, sep = "\t")

script_env <- new.env()
script_env$snakemake <- new("MockSnakemake",
                            input = list(abc_predictions = abc_path, kendall_predictions = kendall_path),
                            output = list(arc_predictions = arc_path),
                            params = list(abc_score_col = "ABC.Score"))
script_messages <- character(0)
withCallingHandlers(
  sys.source(script_path, envir = script_env),
  message = function(m) {
    script_messages <<- c(script_messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  },
  warning = function(w) {
    # genomation's package-load notice is not about this script
    if (grepl("replacing previous import", conditionMessage(w))) invokeRestart("muffleWarning")
  })
arc <- fread(arc_path)
setkey(arc, start)

rna_of <- function(gene) unique(arc[TargetGene == gene, .(mean_log_normalized_rna, RnaDetectedPercent, RnaPseudobulkTPM)])

test_that("every ABC pair is written, with one set of RNA metrics per gene", {
  expect_equal(nrow(arc), nrow(abc))
  expect_equal(arc[, uniqueN(.SD), by = TargetGene, .SDcols = c("RnaPseudobulkTPM")]$V1, rep(1L, nrow(genes)))
})

test_that("genes present in Pairs.Kendall keep their own RNA metrics", {
  for (g in genes[kendall %in% c("ok", "mixed", "allNA"), TargetGene]) {
    tpm <- genes[TargetGene == g, tpm]
    expect_equal(rna_of(g)$RnaPseudobulkTPM, tpm, info = g)
    expect_equal(rna_of(g)$mean_log_normalized_rna, tpm / 100, info = g)
    expect_equal(rna_of(g)$RnaDetectedPercent, min(1, tpm / 30), info = g)
  }
})

test_that("a gene whose Kendall is NA on every pair is not given a prefix-matching gene's metrics", {
  # OXT is a unique prefix of OXTR; it must keep its own zero expression
  expect_equal(rna_of("OXT")$RnaPseudobulkTPM, 0)
  expect_false(isTRUE(all.equal(rna_of("OXT")$RnaPseudobulkTPM, rna_of("OXTR")$RnaPseudobulkTPM)))
  # PCDHB1 is an ambiguous prefix (PCDHB10, PCDHB11); ZZZ1 is no prefix at all
  expect_equal(rna_of("PCDHB1")$RnaPseudobulkTPM, 0)
  expect_equal(rna_of("ZZZ1")$RnaPseudobulkTPM, 0)
  # an exact match wins over a longer symbol it prefixes
  expect_equal(rna_of("ARC")$RnaPseudobulkTPM, 5)
})

test_that("a gene absent from Pairs.Kendall gets NA metrics, not a prefix-matching gene's, and is reported", {
  # MDH1 is a unique prefix of MDH1B (TPM 1.1)
  expect_true(all(is.na(unlist(rna_of("MDH1")))))
  expect_true(any(grepl("^1 of 14 target genes have no RNA metrics.*First genes: MDH1$", trimws(script_messages))))
})

test_that("Kendall and ARC.E2G.Score are computed from the non-NA Kendall pairs only", {
  expected <- pairs[, .(start, Kendall = fifelse(kendall == "absent", NA_real_, Kendall), ABC.Score)]
  setkey(expected, start)
  expect_equal(arc$Kendall, expected$Kendall)
  # pairs without a Kendall value keep their ABC score as the ARC-E2G score
  no_kendall <- is.na(expected$Kendall)
  expect_equal(arc$ARC.E2G.Score[no_kendall], expected$ABC.Score[no_kendall])
  expect_true(all(arc$ARC.E2G.Score[!no_kendall] > 0))
})
