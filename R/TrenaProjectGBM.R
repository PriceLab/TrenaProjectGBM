#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#'
#' @title TrenaProjectGBM-class
#'
#' @name TrenaProjectGBM-class
#' @rdname TrenaProjectGBM-class
#' @aliases TrenaProjectGBM
#' @exportClass TrenaProjectGBM
#'

.TrenaProjectGBM <- setClass("TrenaProjectGBM",
                               contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectPlacneta
#'
#' @description
#' Expression, variant and covariate data for the 21 GBM loci
#'
#' @rdname TrenaProject-class
#'
#' @export
#'
#' @return An object of the TrenaProjectGBM class
#'

TrenaProjectGBM <- function(quiet=TRUE)

{
   genomeName <- "hg38"

   gbm.ad.genes <- sort(c("CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM",
                           "ABCA7", "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1",
                           "SLC24A4", "RIN3", "DSG2", "INPP5D", "MEF2C", "NME8", "ZCWPW1",
                           "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40", "TREM2"))

   footprintDatabaseNames <- c("brain_wellington_16",
                               "brain_wellington_20",
                               "brain_hint_16",
                               "brain_hint_20")

     # very temporarily
   geneInfoTable.path <- system.file(package="TrenaProject", "extdata", "geneInfoTable.RData")

   expressionDirectory <- system.file(package="TrenaProjectGBM", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProjectGBM", "extdata", "variants")
   footprintDatabaseHost <- "khaleesi.systemsbiology.net"

   covariatesFile <- system.file(package="TrenaProjectGBM", "extdata", "covariates", "covariates.RData")
   stopifnot(file.exists(expressionDirectory))
   stopifnot(file.exists(variantsDirectory))
   stopifnot(file.exists(covariatesFile))

   .TrenaProjectGBM(TrenaProject(supportedGenes=gbm.ad.genes,
                                  genomeName=genomeName,
                                  geneInfoTable.path=geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  expressionDirectory=expressionDirectory,
                                  variantsDirectory=variantsDirectory,
                                  covariatesFile=covariatesFile,
                                  quiet=quiet
                                  ))

} # TrenaProjectGBM, the constructor
#----------------------------------------------------------------------------------------------------
