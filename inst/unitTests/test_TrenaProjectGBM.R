library(TrenaProjectGBM)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("gbm"))
   gbm <- TrenaProjectGBM();
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basic()
   test_getVariants()
   test_geneInfoTable()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basic <- function()
{
   printf("---- test_basic")

   expected <-  c("ABCA7", "APOE", "BIN1", "CASS4", "CD2AP", "CD33", "CELF1", "CLU", "CR1", "DSG2", "EPHA1",
                  "FERMT2", "HLA-DRB1", "HLA-DRB5", "INPP5D", "MEF2C", "MS4A6A", "NME8", "PICALM", "PTK2B",
                  "RIN3", "SLC24A4", "SORL1", "TOMM40", "ZCWPW1", "TREM2")

   checkTrue(all(expected %in% getSupportedGenes(gbm)))

   expected <- c("brain_wellington_16", "brain_wellington_20", "brain_hint_16", "brain_hint_20")
   checkTrue(all(expected %in% getFootprintDatabaseNames(gbm)))
   checkEquals(getFootprintDatabaseHost(gbm), "khaleesi.systemsbiology.net")

   expected <- c("cerebellum.15167x263", "rosmap.14235x632", "temporalCortex.15167x264")
   checkTrue(all(expected %in% getExpressionMatrixNames(gbm)))
   mtx <- getExpressionMatrix(gbm, "temporalCortex.15167x264")
   checkEquals(dim(mtx), c(15167, 264))

      # note this next test only works wtih tcx.
      # TODO: get covariate data from Mayo cerebellum and rosmap expression matrices

   tbl.covar <- getCovariatesTable(gbm)
   checkEquals(dim(tbl.covar), c(278, 8))   # includes all samples in mtx retrieved above
   checkTrue(all(colnames(mtx) %in% tbl.covar$ID))


   setTargetGene(gbm, "INPP5D")
   checkEquals(getTargetGene(gbm), "INPP5D")

   tbl.transcripts <- getTranscriptsTable(gbm)
   checkTrue(nrow(tbl.transcripts) >= 3)

   tbl.enhancers <- getEnhancers(gbm)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) > 30)

   tbl.dhs <- getEncodeDHS(gbm)
   checkTrue(nrow(tbl.dhs) > 700)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))

   chromosome <- unique(c(tbl.dhs$chrom, tbl.enhancers$chrom))
   checkEquals(length(chromosome), 1)
   start <- min(tbl.dhs$chromStart)
   end   <- max(tbl.dhs$chromEnd)

   tbl.chipSeq <- getChipSeq(gbm, chrom=chromosome, start=start, end=end, tfs=NA)
   checkTrue(nrow(tbl.chipSeq) > 20000)
   checkEquals(colnames(tbl.chipSeq), c("chrom", "start", "endpos", "tf", "name", "strand", "peakStart", "peakEnd"))

   tbl.chipSeq.stat1 <- getChipSeq(gbm, chrom=chromosome, start=start, end=end, tfs="STAT1")
   checkTrue(nrow(tbl.chipSeq.stat1) > 80)
   checkTrue(nrow(tbl.chipSeq.stat1) < 105)

} # test_basic
#------------------------------------------------------------------------------------------------------------------------
test_getVariants <- function()
{
   message(sprintf("--- test_getVariants"))
   names <- getVariantDatasetNames(gbm)
   dataset.name <- "GWAS.snps"
   checkTrue(dataset.name %in% names)
   tbl.snps <- getVariantDataset(gbm, dataset.name)
   checkTrue(nrow(tbl.snps) > 40000)
   checkEquals(colnames(tbl.snps), c("chrom", "start", "end", "rsid", "variant", "reference", "pval", "pScore"))

} # test_getVariants
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancerRegionsAllGenes <- function()
{
   printf("--- test_getEnhancerRegionsAllGenes")

   goi <- getSupportedGenes(gbm)
   x <- lapply(goi, function(gene) {
             setTargetGene(gbm, gene)
             getGeneEnhancersRegion(gbm, 5)
             })
   regions <- unlist(x)
   names(regions) <- goi
   write(regions, file="regions.txt", sep="\n")
    # https://genome.ucsc.edu/cgi-bin/hgLiftOver
    # created regions-hg19.txt
   # hand edit into enhancerRegions.csv
   # tbl.enhancerRegions <- read.table("enhancerRegions.csv", sep=",", as.is=TRUE)
   # save(tbl.eh, file="../inst/extdata/misc/tbl.enhancerRegions.RData)
   # mv that file to inst/extdata/misc/en

} # test_getEnhancerRegionsAllGenes
#------------------------------------------------------------------------------------------------------------------------
test_geneInfoTable <- function()
{
   printf("--- test_geneInfoTable")
   checkTrue("data.frame" %in% is(getGeneInfoTable(gbm)))

} # test_geneInfoTable
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
