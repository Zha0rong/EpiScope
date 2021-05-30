##' Given an AnnotationDB obtained from AnnotationHub, retrieve gene symbol, ensembl ID and Entrez ID
##'
##'
##' @title ExtractGeneInfo
##' @param annoDb AnnotationDB
##' @return data.frame
##' @importFrom AnnotationDbi select


ExtractGeneInfo <- function(annoDb){
    annoDb=annoDb
    genetable=select(annoDb,keys = keys(annoDb),columns = c("ENSEMBL",'GENENAME','SYMBOL'))
    return(genetable)
}
