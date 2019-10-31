#' Reformat count matrix to 3D pseudocount matrix
#'
#' Sum up the counts per cell type and person for each gene and create a
#' 3 dimensional matrix of genes x individuals x cell types
#'
#' @param expr.singlets gene count matrix (already filtered for multiplets etc.)
#' @param annotation Data frame with annotation of each cell to a cell type and donor,
#' the data frame must be sorted in the same way as the expression matrix (required columns: individual)
#' @param colName Column name of the cell type annotation in the annotation data frame
#'
#' @return 3d matrix of genes x individuals x cell type
#'
#' @export
#'
create.pseudobulk<-function(expr.singlets, annotation, colName="cell.type"){
  #Get dimensions
  individuals<-unique(annotation$individual)
  celltypes<-levels(annotation[,colName])

  expr.array<-array(0, dim=c(nrow(expr.singlets),length(individuals),
                             length(celltypes)),
                    dimnames=c(list(gene.id=rownames(expr.singlets),indiv.id=individuals,
                                    ct.id=celltypes)))

  for(indiv in individuals){
    associated.barcodes<-(annotation$individual==indiv)

    expr.person <- expr.singlets[,associated.barcodes,drop=FALSE]
    cell.annotation.person <- annotation[associated.barcodes,colName]
    expr.array[,as.character(indiv),] <- t(apply(expr.person, 1, tapply, cell.annotation.person,
                                                 sum, na.rm=T))
  }

  expr.array[is.na(expr.array)] <- 0

  return(expr.array)
}

#' Get expressed genes from pseudobulk
#'
#' Return a data frame with all genes, which are defined as expressed per cell type (at least min.counts
#' in perc.indiv fraction of all individuals)
#'
#' @param expr.array 3d pesudobulk matrix of genes x individuals x cell type (output of create.pseudobulk)
#' @param min.counts  More than is number of UMI counts for each gene per person and cell type is required
#' to defined it as expressed in one individual
#' @param perc.indiv  Percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#'
#' @return Data frame with all expressed genes per cell type (columns: gene name, cell type,
#' fraction in how many samples the gene is expressed)
#'
#' @import reshape2
#'
#' @export
#'
calculate.gene.counts<-function(expr.array,min.counts=10,perc.indiv=0.5){

  require(reshape2)

  ## get the percentage of persons with counts > min.count per gene and cell type
  pct.expr <- t(apply(expr.array > min.counts, 1, apply, 2, sum) / dim(expr.array)[2])

  pct.expr.reformated <- melt(pct.expr)
  colnames(pct.expr.reformated) <- c("gene", "cell.type", "percent.expressed")

  #Delete all, which are not expressed in at least half of the individuals
  pct.expr.reformated <- pct.expr.reformated[pct.expr.reformated$percent.expressed > perc.indiv,]

  return(pct.expr.reformated)

}

#' Estimate gene expression probability
#'
#' This function estimates the expression probability of each gene in pseudobulk
#' with a certain cutoff of more than min.counts UMI counts
#'
#' @param mu Estimated mean value in sc data per gene (vector)
#' @param size Estimated size parameter in sc data from the negative binomial fit (1/dispersion parameter), also per gene (vector)
#' @param cells.per.person Mean number of cells per person and cell type
#' @param min.counts  More than is number of UMI counts for each gene per person and cell type is required to defined it as expressed in one individual
#' @param num.indivs  Total sample size
#' @param perc.indiv  Percentage of individuals that need to have this gene expressed to define it as globally expressed
#'
#' @return Vector with expression probabilites for each gene
#'
#' @export
#'
estimate.gene.counts<-function(mu,size,cells.per.person,num.indivs,
                               min.counts=10,perc.indiv=0.5){

  #Generate basis data.frame
  fits.allIndivs<-data.frame(mu=mu,size=size,cells.per.person=cells.per.person)

  #Calculate summed NegBinomial distribution
  fits.allIndivs$size<-fits.allIndivs$size*fits.allIndivs$cells.per.person
  fits.allIndivs$mu<-fits.allIndivs$mu*fits.allIndivs$cells.per.person

  #Calculate probability to observe > n counts in one individual
  fits.allIndivs$count.probs<-1 - pnbinom(min.counts,mu=fits.allIndivs$mu,size=fits.allIndivs$size)

  #Calculate probability that > prob indivs have count > n
  fits.allIndivs$indiv.probs<-1 - pbinom(perc.indiv*num.indivs,num.indivs,fits.allIndivs$count.probs)

  #Return number of expressed cells per celltype
  return(fits.allIndivs$indiv.probs)
}

#' Estimate the mean and dispersion parameter for each gene
#'
#' Use the DEseq together with the size factor normalization poscounts of DESeq2
#'
#' @param counts.ct Count matrix as in import for DEseq
#'
#' @return list with three data frame containing the normalized mean values,
#' the dispersion parameters and the parameters for the mean-dispersion function
#' estimated from DESeq
#'
#' @import DESeq
#'
#' @export
#'
nbinom.estimation<-function(counts.ct){

  require(DESeq)

  #Create a fake condition matrix ...
  cds<-newCountDataSet(counts.ct,rep("celltype",ncol(counts.ct)))

  #Use reimplementation of DEseq2 method for size factor estimation (poscounts)
  sizeFactors(cds)<-sizeFactorsPosCounts(counts(cds,normalize=FALSE))

  #Estimate dispersion curves
  cds <- estimateDispersions(cds, sharingMode="gene-est-only")

  #Get normalized mean
  norm.mean.ct<-rowMeans(counts(cds,normalized=TRUE))
  mean.train.runs<-data.frame(gene=rownames(cds),mean=norm.mean.ct,
                              stringsAsFactors = FALSE)

  #Get dispersion parameter
  disp.train.runs<- data.frame(gene=rownames(fData(cds)),disp=fData(cds)[,1],
                               stringsAsFactors = FALSE)

  #Get the dispersion function
  disp.fit<-attr(fitInfo(cds)$dispFunc,"coefficients")
  disp.fun.param<-data.frame(asymptDisp=disp.fit["asymptDisp"],
                             extraPois=disp.fit["extraPois"])

  return(list(mean.train.runs,disp.train.runs,disp.fun.param))
}

#' Reimplementation of DEseq2 function for size factors with option poscounts
#'
#' Alternative to the standard size factor, better suited for scRNAseq data
#' with a lot of zero values
#'
#' @param counts UMI count matrix
#'
#' @return size factors to normalize each cell
#'
sizeFactorsPosCounts<-function(counts){

  geoMeanNZ <- function(x) {
    if (all(x == 0)) { log(0) } else { sum(log(x[x > 0])) / length(x)}
  }
  loggeomeans <- apply(counts, 1, geoMeanNZ)

  sizeFactors.pos<- apply(counts, 2, function(cnts)
    exp( median( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) & cnts > 0] ) ) )
  sizeFactors.pos <- sizeFactors.pos/exp(mean(log(sizeFactors.pos)))

  return(sizeFactors.pos)
}

#' Estimate the mixed gamma distribution of the mean values
#'
#' WARNING: the currently used package mixR has a bug (can not find local variable, but only global ...)
#' The method estimates the distribution of all gene mean values using a gamma mixed model
#' of two components. It restricts the total number of fitted genes to a certain number num.genes.kept
#' (remove additional genes with a mean of 0) to fit the same number of genes
#' across different conditions (e.g. cell types).
#'
#' @param mean.vals Vector with all normalized mean values (estimated using nbinom.estimation)
#' @param num.genes.kept Total number of genes used for the fit (remove additional genes with a mean of 0)
#'
#' @return Data frame with the six parameters describing the gamma mixed distributions
#'
#' @import mixR
#'
#' @export
#'
mixed.gamma.estimation<-function(mean.vals, num.genes.kept=21000){

  require(mixR)

  #Check how many 0 genes are in the sample
  zeroGenes<-mean.vals==0
  #cat(paste("Number genes with a mean of 0:",sum(zeroGenes),"\n"))

  #Remove some of the zero genes but keep enough to get num.genes.kept genes in the end
  num.zeros.keep<-num.genes.kept-sum(!zeroGenes)
  zeroGenes[zeroGenes][1:num.zeros.keep]<-FALSE
  mean.vals<-mean.vals[!zeroGenes]

  #Set zero values to a very small number because zero values can not be fitted ...
  mean.vals[mean.vals==0]<-1e-5

  #Split value to speed up the fitting procedure
  split.value<<-0.95

  #Get quantile value for split
  max.q<<-quantile(mean.vals,split.value)

  #Fit the gamma mixed distribution
  mean.vals<<-mean.vals
  fit.mixed<-mixfit(mean.vals,ncomp=2,family="gamma",
                    pi=c(split.value,1-split.value),
                    mu=c(mean(mean.vals[mean.vals<max.q]),mean(mean.vals[mean.vals>max.q])),
                    sd=c(sd(mean.vals[mean.vals<max.q]),sd(mean.vals[mean.vals>max.q])))

  gamma.mixed.fits<-data.frame(pi.c1=fit.mixed$pi[1],pi.c2=fit.mixed$pi[2],
                               mu.c1=fit.mixed$mu[1],mu.c2=fit.mixed$mu[2],
                               sd.c1=fit.mixed$sd[1],sd.c2=fit.mixed$sd[2])

  return(gamma.mixed.fits)
}

#' Sample the mean values using a gamma mixed distribution
#'
#' WARNING: The used package might be changed later
#'
#' @param gamma.parameters Data frame with gamma parameters
#' (fitted using mixed.gamma.estimation)
#' @param nGenes Number of genes to sample
#'
#' @return Vector with simulated mean values
#'
#' @import mixR
#'
#' @export
#'
sample.mean.values<-function(gamma.parameters, nGenes=21000){

  require(mixR)

  #Sample the mean values
  mean.vals<-rmixgamma(nGenes,
                       pi=c(gamma.parameters$pi.c1,gamma.parameters$pi.c2),
                       mu=c(gamma.parameters$mu.c1,gamma.parameters$mu.c2),
                       sd=c(gamma.parameters$sd.c1,gamma.parameters$sd.c2))
  return(mean.vals)
}

#' Sample the dispersion values dependent on mean values using the DESeq function parametrization
#'
#' WARNING: The used package might be changed later
#'
#' @param mean.vals Vector of gene means
#' @param disp.parameter Data frame with parameter of mean-dispersion function
#'
#' @return Vector with simulated mean values
#'
#' @import mixR
#'
#' @export
#'
sample.disp.values<-function(mean.vals,disp.parameter){
  disp.vals<-disp.parameter[1,"asymptDisp"]+disp.parameter[1,"extraPois"]/mean.vals
  return(disp.vals)
}

#' Estimate the relationship between the gamma mixed parameter and the mean UMI counts
gamma.umi.relation<-function(){

}

#' Estimate the relationship between mean UMI counts and mean read counts
umi.read.relation<-function(){

}
