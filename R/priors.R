####################################################################################
# Functions for estimating or simulating effect size and expression rank priors
#################################################################################

#' Calculation of gene ranks from a normalized count matrix
#'
#' @param countMatrix normalized count matrix with gene names as row names
#' @param diff.expr.genes Vector of significant genes (DE or eQTLs), the gene names
#' in the vector need to match the row names of the count matrix
#'
#' @return Data frame with expression rank for each DE gene
#'
#' @export
#'
gene.rank.calculation<-function(countMatrix,diff.expr.genes){

  #Calculate mean expression of all genes
  gene.rank.calculation.vector(rowMeans(countMatrix),rownames(countMatrix),diff.expr.genes)

}

#' Calculation of gene ranks using a vector of mean expression per gene
#'
#' @param meanVector Vector with mean value per gene in the sample
#' (use normalized counts)
#' @param geneNames Gene name corresponding to the mean value
#' (need to be in the same order as the mean vector)
#' @param diff.expr.genes Vector of significant genes (DE or eQTLs), the gene names
#' in the vector need to match the row names of the count matrix
#'
#' @return Data frame with expression rank for each DE gene
#'
#' @export
#'
gene.rank.calculation.vector<-function(meanVector,geneNames,diff.expr.genes){

  #Calculate mean expression of all genes
  gene.expr<-data.frame(gene_symbol=geneNames, meanExpr=meanVector,
                        stringsAsFactors = FALSE)
  gene.expr$diff.expressed<-ifelse(gene.expr$gene_symbol %in% diff.expr.genes,1,0)

  #Order genes according to the expression values
  gene.expr<-gene.expr[order(gene.expr$meanExpr, decreasing = T),]
  #Calculate for each gene which fraction of all genes has an expression value equal or smaller (always all   genes before it in the data frame)
  gene.expr$cumFraction<-cumsum(rep(1/nrow(gene.expr),nrow(gene.expr)))
  gene.expr$numberGenes<-cumsum(rep(1,nrow(gene.expr)))
  gene.expr$rank<-rank(-gene.expr$meanExpr,ties.method="min")

  #Delete all not differentially expressed genes
  gene.expr<-gene.expr[gene.expr$diff.expressed==1,]

  return(gene.expr[,c("gene_symbol","cumFraction","rank")])
}

#' Simulation of gene ranks uniformly distributed in a certain interval
#'
#' If the given interval is not completely dividable by the number
#' of genes to simulate (numGenes), the number of simulated ranks might
#' be slightly different from the numGenes.
#'
#' @param start Start of the rank interval (and value of the first simulated rank)
#' @param end End of the rank interval
#' @param numGenes Number of ranks to simulate
#'
#' @return Vector with uniforml distributed gene ranks
#'
#' @export
#'
uniform.ranks.interval<-function(start,end,numGenes){
  return(seq(start,end,by=round((end-start)/numGenes)))
}

#' Simulation of effect sizes for DE genes (FoldChanges)
#'
#' @param mean of the normal distribution of the log fold changes
#' @param sd of the normal distribution of the log fold changes
#' @param numGenes Number of effect sizes to simulate
#'
#' @return Vector with uniforml distributed gene ranks
#'
#' @export
#'
effectSize.DE.simulation<-function(mean,sd,numGenes){
  effectSizes.logFold<-rnorm(numGenes,mean=mean,sd=sd)
  effectSizes<-2^effectSizes.logFold
  return(effectSizes)
}

#' Simulation of effect sizes for eQTL genes (FoldChanges)
#'
#' @param mean of the normal distribution of the log fold changes
#' @param sd of the normal distribution of the log fold changes
#' @param numGenes Number of effect sizes to simulate
#'
#' @return Vector with uniforml distributed gene ranks
#'
#' @export
#'
effectSize.eQTL.simulation<-function(mean,sd,numGenes){

  require(HardyWeinberg)

  sampled.Zscore<-rnorm(numGenes,mean=mean,sd=sd)
  #Set all values below the mean to values larger than the mean ...
  sampled.Zscore[sampled.Zscore<mean]<-2*mean-sampled.Zscore[sampled.Zscore<mean]
  Rsq<-HardyWeinberg::ifisherz(sampled.Zscore)^2
}
