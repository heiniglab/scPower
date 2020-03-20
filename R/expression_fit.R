################################################
# Functions for fitting the expression model
###############################################

#' Reformat count matrix to 3D pseudocount matrix
#'
#' Sum up the counts per cell type and individual for each gene and create a
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
  #Convert annotation column to factor as it makes handling of missing cts per individual easier
  annotation[,colName]<-as.factor(annotation[,colName])
  celltypes<-levels(annotation[,colName])

  expr.array<-array(0, dim=c(nrow(expr.singlets),length(individuals),
                             length(celltypes)),
                    dimnames=c(list(gene.id=rownames(expr.singlets),indiv.id=individuals,
                                    ct.id=celltypes)))

  for(indiv in individuals){
    associated.barcodes<-(annotation$individual==indiv)

    expr.individual <- expr.singlets[,associated.barcodes,drop=FALSE]
    cell.annotation.individual <- annotation[associated.barcodes,colName]
    expr.array[,as.character(indiv),] <- t(apply(expr.individual, 1, tapply, cell.annotation.individual,
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
#' @param min.counts  More than is number of UMI counts for each gene per individual and cell type is required
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

  ## get the percentage of individuals with counts > min.count per gene and cell type
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
#' @param nCellsCt Mean number of cells per individual and cell type
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param min.counts  More than is number of UMI counts for each gene per individual and cell type is required to defined it as expressed in one individual
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
                                  ct,disp.fun.param,
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

  gamma.parameters<-data.frame(p1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="p1"],
                               p2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="p2"],
                               mean1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mean1"],
                               mean2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mean2"],
                               sd1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd1"],
                               sd2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd2"])

  #Convert them to the rateshape variant
  gamma.parameters<-convert.gamma.parameters(gamma.parameters,type="rateshape")

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
#' @param nCellsCt Mean number of cells per individual and cell type
#' @param nSamples  Total sample size
#' @param min.counts  More than is number of UMI counts for each gene per individual and cell type is required to defined it as expressed in one individual
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
#'
#' @param counts.ct Count matrix as in import for DEseq
#' @param sizeFactors Use to different size factor normalization methods, either the "standard"
#' method from DEseq or "poscounts" of DESeq2
#'
#' @return list with three data frame containing the normalized mean values,
#' the dispersion parameters and the parameters for the mean-dispersion function
#' estimated from DESeq
#'
#' @import DESeq
#'
#' @export
#'
nbinom.estimation<-function(counts.ct, sizeFactors="standard"){

  require(DESeq)

  #Create a fake condition matrix ...
  cds<-newCountDataSet(counts.ct,rep("celltype",ncol(counts.ct)))

  if(sizeFactors=="standard"){
    cds<-estimateSizeFactors(cds)
  } else if (sizeFactors=="poscounts"){
    #Use reimplementation of DEseq2 method for size factor estimation (poscounts)
    sizeFactors(cds)<-sizeFactorsPosCounts(counts(cds,normalize=FALSE))
  } else {
    stop("Method for size factor estimation not known. Please use either standard or poscounts!")
  }

  #Check if size factors were calculated correctly
  if(any(is.na(sizeFactors(cds)))){
    stop(paste("Size factor estimation failed! Probably due to sparsity of the matrix!",
               "Try poscounts instead as variable sizeFactors."))
  }

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
#' The estimated distribution will be a mixture of a zero component and two left zensored gamma distributions,
#' using an em fitting procedure implemented in the package. It is recommended to reduce the total number of
#' fitted genes (often in the count matrix a huge fraction of zero genes) using the parameter num.genes.kept.
#'
#' @param mean.vals Vector with all normalized mean values (estimated using nbinom.estimation)
#' @param censoredPoint Censoring point for left censored gamma distributions (if not set,
#' the minimal positive value will be chosen)
#' @param num.genes.kept Total number of genes used for the fit (remove additional genes with a mean of 0)
#' @param return.df If true return data frame with fitted parameters, otherwise an object of class EMResult
#' (containing also the loglikelihoods)
#'
#' @return Data frame with the six parameters describing the gamma mixed distributions
#'
#' @import mixR
#'
#' @export
mixed.gamma.estimation<-function(mean.vals, censoredPoint=NULL,
                                 num.genes.kept=21000, return.df=TRUE){

  #Check how many 0 genes are in the sample
  zeroGenes<-mean.vals==0

  #Remove some of the zero genes but keep enough to get num.genes.kept genes in the end
  num.zeros.keep<-num.genes.kept-sum(!zeroGenes)

  if(num.zeros.keep<0){
    stop(paste("There are",sum(!zeroGenes),"genes with positive expression.",
              "Increase the num.genes.kept parameter to a value larger than that!"))
  } else if (num.zeros.keep>0){
    zeroGenes[zeroGenes][1:num.zeros.keep]<-FALSE
  }

  mean.vals<-mean.vals[!zeroGenes]

  #Set proportions of the distributions
  zero.prop<-sum(mean.vals==0)/(length(mean.vals)*4) #assume 25% of the zeros from the zero component
  outlier.prop<-0.05

  #Set censored point to the minimal positive value if was not set before
  if(is.null(censoredPoint)){
    censoredPoint<-min(mean.vals[mean.vals>0])
  }

  #Initialize first gamma component
  y<-mean.vals[mean.vals > 0 & mean.vals<quantile(1-outlier.prop)]
  Gamma1<-LeftCensoredGamma(shape=mean(y)^2/var(y),rate=mean(y)/var(y),cutoff=censoredPoint)

  #Initialize second gamma component
  y<-mean.vals[mean.vals>=quantile(1-outlier.prop)]
  Gamma2<-LeftCensoredGamma(shape=mean(y)^2/var(y),rate=mean(y)/var(y),cutoff=censoredPoint)

  log<-capture.output({emfit <-em(mean.vals, ncomp=3,
                                  prop=c(zero.prop,1-(zero.prop+outlier.prop),outlier.prop),
                                  model.constructor=c("Zero","Gamma","Gamma"),
                                  models=c(Zero(),Gamma1,Gamma2))})

  if(return.df){
    df<-data.frame(p1=emfit@proportions[1],p2=emfit@proportions[2],
                   s1=emfit@models[[2]]@shape,s2=emfit@models[[3]]@shape,
                   r1=emfit@models[[2]]@rate,r2=emfit@models[[3]]@rate)
    return(df)
  } else {
    return(emfit)
  }
}

#' Randomly sample the mean values using the gamma mixed distribution
#'
#' @param gamma.parameters Data frame with gamma parameters
#' (fitted using mixed.gamma.estimation)
#' @param nGenes Number of genes to sample
#'
#' @return Vector with simulated mean values
#'
#' @export
#'
sample.mean.values.random<-function(gamma.parameters, nGenes=21000){

  if(gamma.parameters$p1<0){
    stop(paste0("Gamma parameter p1 has a value of",
                   gamma.parameters$p1,"smaller than 0!"))
  }

  if(gamma.parameters$p2>1){
    stop(paste0("Gamma parameter p2 has a value of",
                gamma.parameters$p2,"larger than 1!"))
  }

  #Zero Component
  nZeros<-round(nGenes*gamma.parameters$p1)
  vals<-rep(0,nZeros)

  #First Gamma component
  nGamma1<-round(nGenes*gamma.parameters$p2)
  vals<-c(vals,rgamma(n = nGamma1,shape = gamma.parameters$s1,
                      rate = gamma.parameters$r1))

  #Second Gamma component
  nGamma2<-nGenes-nZeros-nGamma1
  vals<-c(vals,rgamma(n = nGamma2,shape = gamma.parameters$s2,
                      rate = gamma.parameters$r2))

  return(vals)
}

#' Draw the mean values using the quantile distributions of the gamma mixed distribution
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

  if(gamma.parameters$p1<0){
    stop(paste0("Gamma parameter p1 has a value of",
                gamma.parameters$p1,"smaller than 0!"))
  }

  if(gamma.parameters$p2>1){
    stop(paste0("Gamma parameter p2 has a value of",
                gamma.parameters$p2,"larger than 1!"))
  }

  #Zero Component
  nZeros<-round(nGenes*gamma.parameters$p1)
  vals<-rep(0,nZeros)

  #First Gamma component
  nGamma1<-round(nGenes*gamma.parameters$p2)
  quantiles.c1<-seq(1/(nGamma1+1),1-1/(nGamma1+1),by=1/(nGamma1+1))
  vals<-c(vals,qgamma(quantiles.c1,shape = gamma.parameters$s1,
                      rate = gamma.parameters$r1))

  #Second Gamma component
  nGamma2<-nGenes-nZeros-nGamma1
  quantiles.c2<-seq(1/(nGamma2+1),1-1/(nGamma2+1),by=1/(nGamma2+1))
  vals<-c(vals,qgamma(quantiles.c2,shape = gamma.parameters$s2,
                      rate = gamma.parameters$r2))

  return(vals)
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

#' Converts two different gamma parameterizations, one using rate
#' and shape (rateshape) and one using mean and sd (meansd)
#'
#' The distributions of the mean values are fitted using the rateshape
#' version, the parameterization over the UMI counts is done using the
#' meansd version.
#'
#' @param gamma.fits Parameters from the fitted gamma function (mixture of two gamma components),
#' either with mean and sd values or with rate and shape values
#' @param type Way of conversion, either to meansd or to rateshape variant
#'
#' @return Converted parameters of the gamma function
#'
#' @export
#'
convert.gamma.parameters<-function(gamma.fits,type="meansd"){
  #Convert rateshape parameterization to meansd
  if(type=="meansd"){
    gamma.fits$mean1<-gamma.fits$s1/gamma.fits$r1
    gamma.fits$sd1<-sqrt(gamma.fits$s1/gamma.fits$r1^2)
    gamma.fits$mean2<-gamma.fits$s2/gamma.fits$r2
    gamma.fits$sd2<-sqrt(gamma.fits$s2/gamma.fits$r2^2)
  #Convert meansd parameterization to rateshape
  } else if (type=="rateshape"){
    gamma.fits$r1<-gamma.fits$mean1 / gamma.fits$sd1^2
    gamma.fits$s1<-gamma.fits$mean1 * gamma.fits$r1
    gamma.fits$r2<-gamma.fits$mean2 / gamma.fits$sd2^2
    gamma.fits$s2<-gamma.fits$mean2 * gamma.fits$r2
  } else {
    stop("Type of gamma conversion not known! Use either meansd or rateshape")
  }

  return(gamma.fits)
}

#' Estimate linear relationship between the gamma mixed parameter and the mean UMI counts
#' or mean read counts (for smartseq data)
#'
#' @param gamma.fits Parameters from the fitted gamma function (p1,p2,mean1,mean2,sd1,sd2),
#' calculated in mixed.gamma.estimation and converted with convert.gamma.parameters,
#' and a column with mean UMI values, can be calculated in mean.umi.calculation
#' @param variable Name of the variable to fit the parameters against (should contain umi counts or read counts)
#'
#' @return Linear fit for each gamma parameter
#'
#' @export
#'
umi.gamma.relation<-function(gamma.fits, variable="mean.umi"){

  #Fit mean and standard deviation parameters linear dependent on the mean UMI values
  parameter.fits<-NULL
  for(param in c("p1","p2","mean1","mean2","sd1","sd2")){

    #Fit a linear relationship
    fit<-lm(data=gamma.fits,as.formula(paste(param,"~",variable)))

    parameter.fits<-rbind(parameter.fits,
                          data.frame(parameter=param,
                                     intercept=fit$coefficients["(Intercept)"],
                                     meanUMI=fit$coefficients[variable]))
  }

  #Delete row names
  rownames(parameter.fits)<-NULL

  return(parameter.fits)
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
