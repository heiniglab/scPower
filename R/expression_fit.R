################################################
# Functions for fitting the expression model
###############################################

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
  #Convert annotation column to factor as it makes handling of missing cts per person easier
  annotation[,colName]<-as.factor(annotation[,colName])
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

#' Estimate gene expression probability based on experimental parameters
#'
#' This function estimates the expression probability of each gene in pseudobulk
#' with a certain cutoff of more than min.counts UMI counts, based on experimental
#' parameters which lead to certain mean and dispersion values for each gene
#'
#' @param nSamples Sample size
#' @param readDepth Target read depth per cell
#' @param nCellsCt Mean number of cells per person and cell type
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param min.counts  More than is number of UMI counts for each gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv  Percentage of individuals that need to have this gene expressed to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#'
#' @return Vector with expression probabilities for each gene
#'
#' @export
#'
estimate.exp.prob.param<-function(nSamples,readDepth,nCellsCt,
                                  read.umi.fit,gamma.mixed.fits,
                                  gamma.probs,ct,disp.fun.param,
                                  min.counts=10,perc.indiv=0.5,
                                  nGenes=21000,samplingMethod="quantiles"){

  #Get mean umi dependent on read depth
  umiCounts<-read.umi.fit$intercept+read.umi.fit$reads*log(readDepth)

  if(umiCounts<=0){
    stop("Read depth too small! UMI model estimates a mean UMI count per cell smaller than 0!")
  }

  #Check if gamma fit data for the cell type exists
  if(! any(gamma.mixed.fits$ct==ct)){
    stop(paste("No gene curve fitting data in the data frame gamma.mixed.fits fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right gamma.mixed.fits object used."))
  }

  #Get gamma values dependent on mean umi
  gamma.fits.ct<-gamma.mixed.fits[gamma.mixed.fits$ct==ct,]
  gamma.fits.ct$fitted.value<-gamma.fits.ct$intercept+gamma.fits.ct$meanUMI*umiCounts

  #Check if gamma probability data for the cell type exists
  if(! any(gamma.probs$ct==ct)){
    stop(paste("No gene curve fitting data in the data frame gamma.probs fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right gamma.probs object used."))
  }

  gamma.parameters<-data.frame(pi.c1=gamma.probs$pi.c1[gamma.probs$ct==ct],
                               pi.c2=1-gamma.probs$pi.c1[gamma.probs$ct==ct],
                               mu.c1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mu.c1"],
                               mu.c2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mu.c2"],
                               sd.c1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd.c1"],
                               sd.c2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd.c2"])

  #Sample means values
  if(samplingMethod=="random"){
    gene.means<-sample.mean.values.random(gamma.parameters,nGenes)
  } else if (samplingMethod=="quantiles"){
    gene.means<-sample.mean.values.quantiles(gamma.parameters,nGenes)
  } else {
    stop("No known sampling method. Use the options 'random' or 'quantiles'.")
  }

  #Check if dispersion data for the cell type exists
  if(! any(disp.fun.param$ct==ct)){
    stop(paste("No dispersion fitting data in the data frame disp.fun.param fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right disp.fun.param object used."))
  }

  #Fit dispersion parameter dependent on mean parameter
  disp.fun<-disp.fun.param[disp.fun.param$ct==ct,]
  gene.disps<-sample.disp.values(gene.means,disp.fun)

  return(estimate.exp.prob.values(gene.means,1/gene.disps,nCellsCt,
                           nSamples,min.counts,perc.indiv))

}

#' Estimate gene expression probability based on mean and dispersion values
#'
#' This function estimates the expression probability of each gene in pseudobulk
#' with a certain cutoff of more than min.counts UMI counts, based on the mean and
#' the disperion value of each gene
#'
#' @param mu Estimated mean value in sc data per gene (vector)
#' @param size Estimated size parameter in sc data from the negative binomial fit (1/dispersion parameter), also per gene (vector)
#' @param nCellsCt Mean number of cells per person and cell type
#' @param nSamples  Total sample size
#' @param min.counts  More than is number of UMI counts for each gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv  Percentage of individuals that need to have this gene expressed to define it as globally expressed
#'
#' @return Vector with expression probabilities for each gene
#'
#' @export
#'
estimate.exp.prob.values<-function(mu,size,nCellsCt,nSamples,
                                   min.counts=10,perc.indiv=0.5){

  #Generate basis data.frame
  fits.allIndivs<-data.frame(mu=mu,size=size,nCellsCt=nCellsCt)

  #Calculate summed NegBinomial distribution
  fits.allIndivs$size<-fits.allIndivs$size*fits.allIndivs$nCellsCt
  fits.allIndivs$mu<-fits.allIndivs$mu*fits.allIndivs$nCellsCt

  #Calculate probability to observe > n counts in one individual
  fits.allIndivs$count.probs<-1 - pnbinom(min.counts,mu=fits.allIndivs$mu,size=fits.allIndivs$size)

  #Calculate probability that > prob indivs have count > n
  fits.allIndivs$indiv.probs<-1 - pbinom(perc.indiv*nSamples,nSamples,fits.allIndivs$count.probs)

  #Return expression probability of each gene
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
  cds <- suppressWarnings(estimateDispersions(cds, sharingMode="gene-est-only"))

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
#' @param default.zero.val Zero values can not be fitted with gamma curves,
#' therefore all zero values are replaced by a very small value
#'
#' @return Data frame with the six parameters describing the gamma mixed distributions
#'
#' @import mixR
#'
#' @export
#'
mixed.gamma.estimation<-function(mean.vals, num.genes.kept=21000,
                                 default.zero.val=1e-5){

  require(mixR)

  #Check how many 0 genes are in the sample
  zeroGenes<-mean.vals==0
  #cat(paste("Number genes with a mean of 0:",sum(zeroGenes),"\n"))

  #Remove some of the zero genes but keep enough to get num.genes.kept genes in the end
  num.zeros.keep<-num.genes.kept-sum(!zeroGenes)

  if(num.zeros.keep<0){
    stop(paste("There are",sum(!zeroGenes),"genes with positive expression.",
               "Increase the num.genes.kept parameter to a value larger than that!"))
  } else if (num.zeros.keep>0){
    zeroGenes[zeroGenes][1:num.zeros.keep]<-FALSE
  }

  mean.vals<-mean.vals[!zeroGenes]

  #Set zero values to a very small number because zero values can not be fitted ...
  mean.vals[mean.vals==0]<-default.zero.val

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

#' Randomly sample the mean values using a gamma mixed distribution
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
sample.mean.values.random<-function(gamma.parameters, nGenes=21000){

  require(mixR)

  #Sample the mean values
  mean.vals<-rmixgamma(nGenes,
                       pi=c(gamma.parameters$pi.c1,gamma.parameters$pi.c2),
                       mu=c(gamma.parameters$mu.c1,gamma.parameters$mu.c2),
                       sd=c(gamma.parameters$sd.c1,gamma.parameters$sd.c2))
  return(mean.vals)
}

#' Draw the mean values using the two quantile distributions of the gamma mixed distribution
#'
#' @param gamma.parameters Data frame with gamma parameters
#' (fitted using mixed.gamma.estimation)
#' @param nGenes Number of genes to sample
#'
#' @return Vector with simulated mean values
#'
#' @export
#'
sample.mean.values.quantiles<-function(gamma.parameters, nGenes=21000){

  #Distribution of genes between the quantiles
  nGenes.c1<-round(nGenes*gamma.parameters$pi.c1)
  nGenes.c2<-nGenes-nGenes.c1

  #Reformate the gamma parameters
  gamma.parameters$rate.c1<-gamma.parameters$mu.c1/gamma.parameters$sd.c1^2
  gamma.parameters$shape.c1<-gamma.parameters$mu.c1*gamma.parameters$rate.c1
  gamma.parameters$rate.c2<-gamma.parameters$mu.c2/gamma.parameters$sd.c2^2
  gamma.parameters$shape.c2<-gamma.parameters$mu.c2*gamma.parameters$rate.c2

  #Calculate which quantile values should be checked to get the total number of genes
  quantiles.c1<-seq(1/(nGenes.c1+1),1-1/(nGenes.c1+1),by=1/(nGenes.c1+1))
  quantiles.c2<-seq(1/(nGenes.c2+1),1-1/(nGenes.c2+1),by=1/(nGenes.c2+1))

  #Sample from each component the quantiles
  means.c1<-qgamma(quantiles.c1,shape=gamma.parameters$shape.c1,
                   rate=gamma.parameters$rate.c1)
  means.c2<-qgamma(quantiles.c2,shape=gamma.parameters$shape.c2,
                   rate=gamma.parameters$rate.c2)

  mean.vals<-c(means.c1,means.c2)

  return(mean.vals)
}

#' Sample the dispersion values dependent on mean values using the DESeq function parametrization
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

#' Calculate the mean UMI count per cell
#'
#' @param countMatrix Raw UMI count matrix with genes as rows and cells as columns
#'
#' @return Mean UMI count per cell
#'
#' @export
#'
meanUMI.calculation<-function(countMatrix){
  return(mean(colSums(countMatrix)))
}

#' Estimate linear relationship between the gamma mixed parameter and the mean UMI counts
#'
#' @param gamma.fits Parameters from the fitted gamma function (pi.c1,pi.c2,mu.c1,mu.c2,sd.c1,sd.c2),
#' calculated in mixed.gamma.estimation, and a column with mean UMI values, can be calculated in
#' mean.umi.calculation
#'
#' @return List with linear fit for gamma parameters mean and sd and mean values for gamma parameters probability
#'
#' @export
#'
umi.gamma.relation<-function(gamma.fits){

  #Fit mean and standard deviation parameters linear dependent on the mean UMI values
  parameter.fits<-NULL
  for(param in c("mu.c1","mu.c2","sd.c1","sd.c2")){

    #Fit a linear relationship
    fit<-lm(data=gamma.fits,as.formula(paste0(param," ~ mean.umi")))

    parameter.fits<-rbind(parameter.fits,
                          data.frame(parameter=param,
                                     intercept=fit$coefficients["(Intercept)"],
                                     meanUMI=fit$coefficients["mean.umi"]))
  }

  #Set probability parameter independent of the mean UMI values
  probability<-data.frame(pi.c1=mean(gamma.fits$pi.c1))

  return(list(parameter.fits,probability))
}

#' Estimate logarithmic relationship between mean UMI counts and mean mapped read counts per cell
#'
#' @param read.umis Data.frame with column for mean UMI values (mean.umi) and
#' column for mean transcriptome mapped reads (transcriptome.mapped.reads)
#'
#' @return Parameters of logarithmic fit between reads and UMIs
#'
#' @export
#'
umi.read.relation<-function(read.umis){
  fit<-lm(data=read.umis,mean.umi~log(transcriptome.mapped.reads))

  read.umi.fit<-data.frame(intercept=fit$coefficients["(Intercept)"],
                           reads=fit$coefficients["log(transcriptome.mapped.reads)"])

  return(read.umi.fit)
}

#' Get median values of parameters from the mean-dispersion fits
#' (no correlation with UMI counts visible)
#'
#' @param disp.funs Data.frame with parameters of mean-dispersion fit
#' (asymptDisp, extraPois)
#'
#' @return Median parameters
#'
#' @export
#'
dispersion.function.estimation<-function(disp.funs){
  disp.fun.general<-data.frame(asymptDisp=median(disp.funs$asymptDisp),
                               extraPois=median(disp.funs$extraPois))

  return(disp.fun.general)
}
