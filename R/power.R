######################################
# Functions for power calculation
#####################################

#' Power calculation for cell type identification
#'
#' Calculate probability to detect at least min.num.cells of a cell type
#' with class frequency cell.type.frac in a sample of size nCells.
#' for each of nSamples individual
#' @param nCells Number of cells measured for each individuum
#' @param min.num.cells Minimal number of the cells for the cell type of interest that should be detected
#' @param cell.type.frac Frequency of the cell type of interest
#' @param nSamples Sample size
#'
#' @return Power to detect the cell type of interest
#'
#' @export
power.detect.celltype<-function(nCells,min.num.cells,cell.type.frac,nSamples){
  return(pnbinom(nCells-min.num.cells,min.num.cells,cell.type.frac)^nSamples)
}

#' Cell sample size calculation for cell type identification
#'
#' Calculate required number of cells per person to detect at least min.num.cells of a cell type
#' with class frequency cell.type.frac for each of nSamples individual.
#' with a probability of prob.cutoff
#' @param prob.cutoff Target power to detect the cell type
#' @param min.num.cells Minimal number of the cells for the cell type of interest that should be detected
#' @param cell.type.frac Frequency of the cell type of interest
#' @param nSamples Sample size (number of individuals)
#'
#' @return Required number of cells per person to reach the target power
#'
#' @export
number.cells.detect.celltype<-function(prob.cutoff,min.num.cells,cell.type.frac,nSamples){
  return(qnbinom(prob.cutoff^(1/nSamples),min.num.cells,cell.type.frac)+min.num.cells)
}

#' Power calculation for a DE/eQTL study  with 10X design (with a restricted number of persons per lane)
#'
#' This function to calculate the detection power for a DE or eQTL study, given DE/eQTL genes from a reference study,
#' in a single cell 10X RNAseq study. The power depends on the cost determining parameter of sample size, number of cells
#' per person and read depth.
#'
#' @param nSamples Sample size
#' @param nCells Number of cells per person
#' @param readDepth Target read depth per cell
#' @param ct.freq Frequency of the cell type of interest
#' @param type (eqtl/de) study
#' @param ref.study Data frame with reference studies to be used for expression ranks and effect sizes
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param samplesPerLane Maximal number of persons per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.UMI.counts Expression defition parameter: more than is number of UMI counts for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane, "linear" if overloading should be
#' modeled explicitly, otherwise "constant". The default value for the parameter multipletRate is matching the option "linear".
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
#'
#' @examples
#' power.general.withDoublets(83,1000,1000,0.2,"de",de.ref.study, "Blueprint (CLL) iCLL-mCLL",
#' 8,read.umi.fit,gamma.mixed.fits,"CD4 T cells",disp.fun.param)
#'
power.general.withDoublets<-function(nSamples,nCells,readDepth,ct.freq,
                                     type,ref.study,ref.study.name,
                                     samplesPerLane,
                                     read.umi.fit,gamma.mixed.fits,
                                     gamma.probs,ct,
                                     disp.fun.param,
                                     mappingEfficiency=0.8,
                                     multipletRate=7.67e-06,multipletFactor=1.82,
                                     min.UMI.counts=10,perc.indiv.expr=0.5,
                                     nGenes=21000,samplingMethod="quantiles",
                                     multipletRateGrowth="linear"){

  #Estimate multiplet fraction dependent on cells per lane
  if(multipletRateGrowth=="linear"){
    multipletFraction<-multipletRate*nCells*samplesPerLane
  } else if (multipletRateGrowth == "constant") {
    multipletFraction<-multipletRate
  } else {
    stop("No known option for multipletRateGrowth. Use the values 'linear' or 'constant'.")
  }

  #Check that the number of cells entered does not provide a multiplet rate of >100%
  if(multipletFraction>=1){
    stop("Too many cells per person entered! Multiplet rate of more than 100%!")
  }

  usableCells<-round((1-multipletFraction)*nCells)
  #Estimate multiplet rate and "real read depth"
  readDepthSinglet<-readDepth*nCells/(usableCells+multipletFactor*(nCells-usableCells))

  #Estimate fraction of correctly mapped reads
  mappedReadDepth<-readDepthSinglet*mappingEfficiency

  #Get the fraction of cell type cells
  ctCells<-round(usableCells*ct.freq)

  #Get mean umi dependent on read depth
  umiCounts<-read.umi.fit$intercept+read.umi.fit$reads*log(mappedReadDepth)

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

  sim.genes<-data.frame(mean=gene.means, disp=gene.disps)

  #Sort simulated genes by mean expression
  sim.genes<-sim.genes[order(sim.genes$mean, decreasing = TRUE),]

  #Calculate the mean per cell type for each individuum
  sim.genes$mean.sum<-sim.genes$mean*ctCells
  sim.genes$disp.sum<-sim.genes$disp/ctCells

  #Calculate for each gene the expression probability with the definition
  #expressed in 50% of the individuals with count > 10
  sim.genes$exp.probs<-estimate.exp.prob.values(sim.genes$mean,1/sim.genes$disp,ctCells,
                                            nSamples=nSamples,min.counts=min.UMI.counts,
                                            perc.indiv=perc.indiv.expr)

  #Calculate the expected number of expressed genes
  exp.genes<-round(sum(sim.genes$exp.probs))

  #Check if the study reference name exists in the data frame
  if(! any(ref.study$name==ref.study.name)){
    stop(paste("No study name in the data frame ref.study fits to the specified reference study name",
               ref.study.name,". Check that the name is correctly spelled and the right ref.study.name object used."))
  }

  #Set the simulated DE genes as the genes at the same rank position as the original DE genes
  ranks<-ref.study$rank[ref.study$name==ref.study.name]

  #Set all DE with rank > 21000 to 21000 (expression anyway nearly 0)
  ranks[ranks>nGenes]<-nGenes

  #Calculate alpha parameter corrected for multiple testing
  alpha<-0.05/exp.genes

  #Choose the DE genes according to the rank
  foundSignGenes<-sim.genes[ranks,]

  #Calculate power
  if(type=="eqtl"){

    #Set the Rsq respectively
    foundSignGenes$Rsq<-ref.study$Rsq[ref.study$name==ref.study.name]

    #Check that the required column Rsq exists
    if(! any(colnames(ref.study)=="Rsq")){
      stop(paste("Column name Rsq missing in the ref.study data frame.",
                 "Please make sure to provide this column for eQTL power analysis."))
    }

    foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.eqtl(foundSignGenes$Rsq[i],
                                                                               alpha,
                                                                               nSamples))
  } else if (type=="de") {

    #Check that the required column FoldChange exists
    if(! any(colnames(ref.study)=="FoldChange")){
      stop(paste("Column name FoldChange missing in the ref.study data frame.",
                 "Please make sure to provide this column for DE power analysis."))
    }

    #Set the fold change respectively
    foundSignGenes$FoldChange<-ref.study$FoldChange[ref.study$name==ref.study.name]

    foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) powerScPop:::power.de(
      floor(nSamples/2),
      foundSignGenes$mean.sum[i],
      foundSignGenes$FoldChange[i],
      1/foundSignGenes$disp.sum[i],
      alpha,3,ssize.ratio=1))
  } else  {
    print('For type parameter only "eqtl" or "de" possible!')
  }

  #Calculate total probability as the DE power times the expression probability
  foundSignGenes$combined.prob<-foundSignGenes$power*foundSignGenes$exp.probs

  power.study<-data.frame(name=ref.study.name,
                          powerDetect=mean(foundSignGenes$combined.prob),
                          exp.probs=mean(foundSignGenes$exp.probs),
                          power=mean(foundSignGenes$power),
                          sampleSize=nSamples,
                          totalCells=nCells,
                          usableCells=usableCells,
                          multipletFraction=multipletFraction,
                          ctCells=ctCells,
                          readDepth=readDepth,
                          readDepthSinglet=readDepthSinglet,
                          mappedReadDepth=mappedReadDepth,
                          expressedGenes=exp.genes)

  return(power.study)
}

#' Power calculation for a DE/eQTL study with 10X design (with a restricted number of cells per lane)
#'
#' This function is a variant of power.general.withDoublets, where not the number of samplesPerLane is given as an
#' parameter, but instead the individuals are distributed over the lanes in a way that restricts the total number of
#' cells per lane instead. This gives also an upper bound for the doublet rate.
#'
#' @param nSamples Sample size
#' @param nCells Number of cells per person
#' @param readDepth Target read depth per cell
#' @param ct.freq Frequency of the cell type of interest
#' @param type (eqtl/de) study
#' @param ref.study Data frame with reference studies to be used for expression ranks and effect sizes
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param cellsPerLane Maximal number of cells per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.UMI.counts Expression defition parameter: more than is number of UMI counts for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane, "linear" if overloading should be
#' modeled explicitly, otherwise "constant". The default value for the parameter multipletRate is matching the option "linear".
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
power.general.restrictedDoublets<-function(nSamples,nCells,readDepth,ct.freq,
                                           type,ref.study,ref.study.name,
                                           cellsPerLane,
                                           read.umi.fit,gamma.mixed.fits,
                                           gamma.probs,ct,
                                           disp.fun.param,
                                           mappingEfficiency=0.8,
                                           multipletRate=7.67e-06,multipletFactor=1.82,
                                           min.UMI.counts=10,perc.indiv.expr=0.5,
                                           nGenes=21000,samplingMethod="quantiles",
                                           multipletRateGrowth="linear"){

  #Distribute persons most optimal over the lanes
  samplesPerLane<-floor(cellsPerLane/nCells)

  if(samplesPerLane==0){
    stop("Allowed number of cells per lane is too low to fit so many cells per person!")
  }

  return(power.general.withDoublets(nSamples,nCells,readDepth,ct.freq,
                             type,ref.study,ref.study.name,
                             samplesPerLane,
                             read.umi.fit,gamma.mixed.fits,
                             gamma.probs,ct,
                             disp.fun.param,
                             mappingEfficiency,
                             multipletRate,multipletFactor,
                             min.UMI.counts,perc.indiv.expr,
                             nGenes,samplingMethod,multipletRateGrowth))
}

#' Power calculation for a DE/eQTL study  with Smart-seq design
#'
#' This function to calculate the detection power for a DE or eQTL study, given DE/eQTL genes from a reference study,
#' in a single cell Smart-seq RNAseq study. The power depends on the cost determining parameter of sample size, number of cells
#' per person and read depth.
#'
#' @param nSamples Sample size
#' @param nCells Number of cells per person
#' @param readDepth Target read depth per cell
#' @param ct.freq Frequency of the cell type of interest
#' @param type (eqtl/de) study
#' @param ref.study Data frame with reference studies to be used for expression ranks and effect sizes
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanReads (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.linear.fit Function to fit the dispersion parameter dependent on the mean (parameter linear dependent on read depth)
#' (required columns: parameter, ct (cell type), intercept, meanReads (slope))
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletFraction Multiplet fraction in the experiment as a constant factor
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.norm.count Expression defition parameter: more than is number of counts per kilobase transcript for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
power.smartseq<-function(nSamples,nCells,readDepth,ct.freq,
                         type,ref.study,ref.study.name,
                         gamma.mixed.fits,gamma.probs,ct,
                         disp.linear.fit,
                         mappingEfficiency=0.8,
                         multipletFraction=0,multipletFactor=1.82,
                         min.norm.count=10,perc.indiv.expr=0.5,
                         nGenes=21000,samplingMethod="quantiles"){

  usableCells<-round((1-multipletFraction)*nCells)
  #Estimate multiplet rate and "real read depth"
  readDepthSinglet<-readDepth*nCells/(usableCells+multipletFactor*(nCells-usableCells))

  #Estimate fraction of correctly mapped reads
  mappedReadDepth<-readDepth*mappingEfficiency

  ctCells<-usableCells*ct.freq

  #Check if gamma fit data for the cell type exists
  if(! any(gamma.mixed.fits$ct==ct)){
    stop(paste("No gene curve fitting data in the data frame gamma.mixed.fits fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right gamma.mixed.fits object used."))
  }

  #Get gamma values dependent on mean reads
  gamma.fits.ct<-gamma.mixed.fits[gamma.mixed.fits$ct==ct,]
  gamma.fits.ct$fitted.value<-gamma.fits.ct$intercept+gamma.fits.ct$meanReads*mappedReadDepth

  if(any(gamma.fits.ct$fitted.value<=0)){
    stop("At least one of the gamma parameter got negative for this read depth.",
         "Choose a higher read depth or a different gamma - read fit.")
  }

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

  sim.genes<-data.frame(mean=gene.means)

  #Calculate the mean per cell type for each individuum
  sim.genes$mean.sum<-sim.genes$mean*ctCells

  #Sort simulated genes by mean expression
  sim.genes<-sim.genes[order(sim.genes$mean, decreasing = TRUE),]

  #Check if the study reference name exists in the data frame
  if(! any(ref.study$name==ref.study.name)){
    stop(paste("No study name in the data frame ref.study fits to the specified reference study name",
               ref.study.name,". Check that the name is correctly spelled and the right ref.study.name object used."))
  }

  #Assign a gene length for each gene (default 5000), for known DE genes the real value
  sim.genes$geneLength<-5000
  temp.ranks<-ref.study[ref.study$rank<nGenes & ref.study$name==ref.study.name,]
  sim.genes$geneLength[temp.ranks$rank]<-temp.ranks$geneLength

  #Check if dispersion data for the cell type exists
  if(! any(disp.linear.fit$ct==ct)){
    stop(paste("No dispersion fitting data in the data frame disp.linear.fit fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right disp.linear.fit object used."))
  }

  #Get dispersion values dependent on mean reads
  disp.fun.ct<-disp.linear.fit[disp.linear.fit$ct==ct,]
  disp.fun.ct$fitted.value<-disp.fun.ct$intercept+disp.fun.ct$meanReads*mappedReadDepth

  #Fit dispersion parameter dependent on mean parameter
  disp.fun<-data.frame(asymptDisp=disp.fun.ct$fitted.value[disp.fun.ct$parameter=="asymptDisp"],
                       extraPois=disp.fun.ct$fitted.value[disp.fun.ct$parameter=="extraPois"],
                       ct=ct)

  #Fit dispersion parameter dependent on mean parameter
  sim.genes$mean.length.transformed<-sim.genes$mean*sim.genes$geneLength/1000
  sim.genes$disp<-sample.disp.values(sim.genes$mean.length.transformed,disp.fun)
  sim.genes$disp.sum<-sim.genes$disp/ctCells

  #Fit also length transformed sum
  sim.genes$mean.length.sum<-sim.genes$mean.length.transformed*ctCells

  #Calculate for each gene the expression probability with the definition
  #expressed in 50% of the individuals with count > 10
  sim.genes$exp.probs<-estimate.exp.prob.values(sim.genes$mean,1/sim.genes$disp,ctCells,
                                            nSamples=nSamples,min.counts=min.norm.count,
                                            perc.indiv=perc.indiv.expr)

  #Calculate the expected number of expressed genes
  exp.genes<-round(sum(sim.genes$exp.probs))

  #Set the simulated DE genes as the genes at the same rank position as the original DE genes
  ranks<-ref.study$rank[ref.study$name==ref.study.name]

  #Set all DE with rank > 21000 to 21000 (expression anyway nearly 0)
  ranks[ranks>nGenes]<-nGenes

  #Calculate alpha parameter corrected for multiple testing
  alpha<-0.05/exp.genes

  #Choose the DE genes according to the rank
  foundSignGenes<-sim.genes[ranks,]

  #Calculate power
  if(type=="eqtl"){

    #Check that the required column Rsq exists
    if(! any(colnames(ref.study)=="Rsq")){
      stop(paste("Column name Rsq missing in the ref.study data frame.",
                 "Please make sure to provide this column for eQTL power analysis."))
    }

    #Set the Rsq respectively
    foundSignGenes$Rsq<-ref.study$Rsq[ref.study$name==ref.study.name]

    foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.eqtl(foundSignGenes$Rsq[i],
                                                                                            alpha,
                                                                                            nSamples))
  } else if (type=="de") {

    #Check that the required column FoldChange exists
    if(! any(colnames(ref.study)=="FoldChange")){
      stop(paste("Column name FoldChange missing in the ref.study data frame.",
                 "Please make sure to provide this column for DE power analysis."))
    }

    #Set the fold change respectively
    foundSignGenes$FoldChange<-ref.study$FoldChange[ref.study$name==ref.study.name]

    foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) powerScPop:::power.de(
      floor(nSamples/2),
      foundSignGenes$mean.length.sum[i],
      foundSignGenes$FoldChange[i],
      1/foundSignGenes$disp.sum[i],
      alpha,3,ssize.ratio=1))
  } else  {
    print('For type parameter only "eqtl" or "de" possible!')
  }

  #Calculate total probability as the DE power times the expression probability
  foundSignGenes$combined.prob<-foundSignGenes$power*foundSignGenes$exp.probs

  power.study<-data.frame(name=ref.study.name,
                          powerDetect=mean(foundSignGenes$combined.prob),
                          exp.probs=mean(foundSignGenes$exp.probs),
                          power=mean(foundSignGenes$power),
                          sampleSize=nSamples,
                          totalCells=nCells,
                          usableCells=usableCells,
                          multipletFraction=multipletFraction,
                          ctCells=ctCells,
                          readDepth=readDepth,
                          readDepthSinglet=readDepthSinglet,
                          mappedReadDepth=mappedReadDepth,
                          expressedGenes=exp.genes)

  return(power.study)
}

#' Optimizing cost parameters to maximize detection power for a given budget and 10X design
#'
#' This function determines the optimal parameter combination for a given budget.
#' The optimal combination is thereby the one with the highest detection power.
#' Of the three parameters sample size, cells per sample and read depth, two need to be set and
#' the third one is uniquely defined given the other two parameters and the overall budget.
#'
#' @param totalBudget Overall experimental budget
#' @param type (eqtl/de) study
#' @param ct Cell type of interest (name should be found in the gamma and dispersion models)
#' @param ct.freq Frequency of the cell type of interest
#' @param costKit Cost for one 10X kit
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param ref.study Data frame with reference studies to be used for effect sizes and ranks
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param samplesPerLane Maximal number of persons per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per person that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.UMI.counts Expression defition parameter: more than is number of UMI counts for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane, "linear" if overloading should be
#' modeled explicitly, otherwise "constant". The default value for the parameter multipletRate is matching the option "linear".
#'
#' @export
#'
#' @examples
#' optimize.constant.budget(10000,seq(1000,10000,by=1000),seq(1000,10000,by=1000),
#' 5600,14032,4100*10^6,0.2,"de",de.ref.study,"Blueprint (CLL) iCLL-mCLL",8,
#' read.umi.fit,gamma.mixed.fits,"CD4 T cells",disp.fun.param)
#'
optimize.constant.budget<-function(totalBudget,type,
                                   ct,ct.freq,
                                   costKit,costFlowCell,readsPerFlowcell,
                                   ref.study,ref.study.name,
                                   samplesPerLane,
                                   read.umi.fit,gamma.mixed.fits,
                                   gamma.probs, disp.fun.param,
                                   nSamplesRange=NULL,
                                   nCellsRange=NULL, readDepthRange=NULL,
                                   mappingEfficiency=0.8,
                                   multipletRate=7.67e-06,multipletFactor=1.82,
                                   min.UMI.counts=10,perc.indiv.expr=0.5,
                                   nGenes=21000,samplingMethod="quantiles",
                                   multipletRateGrowth="linear"){

  #Check that exactly two of the parameters are set and the third one is not defined
  if(sum(is.null(nSamplesRange),is.null(nCellsRange),is.null(readDepthRange))!=1){
    stop(paste("To optimize the experimental design for a given budget,",
               "always exactly one of the parameters nSamplesRange, nCellsRange",
               "and readDepthRange should be set to NULL."))
  }

  #Case 1: estimate the sample size
  if(is.null(nSamplesRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nCellsRange,readDepthRange)
    colnames(param.combis)<-c("nCells","readDepth")

    #Sample size dependent on the budget
    param.combis$nSamples<-sapply(1:nrow(param.combis),
                                  function(i)floor(sampleSizeBudgetCalculation(param.combis$nCells[i],
                                                                               param.combis$readDepth[i],
                                                                               totalBudget,
                                                                               costKit,samplesPerLane,
                                                                               costFlowCell,readsPerFlowcell)))
  #Case 2: estimate the number of cells per persons
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                   function(i)cellsPersonBudgetCalculation(param.combis$nSamples[i],
                                                                          param.combis$readDepth[i],
                                                                          totalBudget,
                                                                          costKit,samplesPerLane,
                                                                          costFlowCell,readsPerFlowcell))
  # Case 3: estimate the read depth
  } else {
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,nCellsRange)
    colnames(param.combis)<-c("nSamples","nCells")

    param.combis$readDepth<-sapply(1:nrow(param.combis),
                                   function(i)readDepthBudgetCalculation(param.combis$nSamples[i],
                                                                         param.combis$nCells[i],
                                                                         totalBudget,
                                                                         costKit,samplesPerLane,
                                                                         costFlowCell,readsPerFlowcell))
  }

  #Remove all combinations where one of the parameters is <=0
  if(any(param.combis$nSamples==0) | any(param.combis$nCells<=0) | any(param.combis$readDepth<=0)){
    warning("Some of the parameter combinations are too expensive and removed from the parameter grid.")
    param.combis<-param.combis[param.combis$nSamples>0 & param.combis$nCells>0 & param.combis$readDepth>0,]
  }

  power.study<-mapply(power.general.withDoublets,
                      param.combis$nSamples,
                      param.combis$nCells,
                      param.combis$readDepth,
                      MoreArgs=list(ct.freq=ct.freq,
                                    multipletRate=multipletRate,
                                    multipletFactor=multipletFactor,
                                    type=type,
                                    ref.study=ref.study,
                                    ref.study.name=ref.study.name,
                                    samplesPerLane=samplesPerLane,
                                    read.umi.fit=read.umi.fit,
                                    gamma.mixed.fits=gamma.mixed.fits,
                                    gamma.probs=gamma.probs,
                                    ct=ct,
                                    disp.fun.param=disp.fun.param,
                                    mappingEfficiency=mappingEfficiency,
                                    min.UMI.counts=min.UMI.counts,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    multipletRateGrowth=multipletRateGrowth))

  power.study<-data.frame(apply(power.study,1,unlist),stringsAsFactors = FALSE)
  power.study[,2:ncol(power.study)]<-apply(power.study[,2:ncol(power.study)],2,as.numeric)

  return(power.study)
}

#' Optimizing cost parameters to maximize detection power for a given budget with
#' library preparation costs per cell
#'
#' This function determines the optimal parameter combination for a given budget.
#' The optimal combination is thereby the one with the highest detection power.
#' Of the three parameters sample size, cells per sample and read depth, two need to be set and
#' the third one is uniquely defined given the other two parameters and the overall budget.
#'
#' @param totalBudget Overall experimental budget
#' @param type (eqtl/de) study
#' @param ct Cell type of interest (name should be found in the gamma and dispersion models)
#' @param ct.freq Frequency of the cell type of interest
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param ref.study Data frame with reference studies to be used for effect sizes and ranks
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param samplesPerLane Maximal number of persons per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per person that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.UMI.counts Expression defition parameter: more than is number of UMI counts for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane, "linear" if overloading should be
#' modeled explicitly, otherwise "constant". The default value for the parameter multipletRate is matching the option "linear".
#'
#' @export
optimize.constant.budget.libPrepCell<-function(totalBudget, type,
                                               ct,ct.freq,
                                               prepCostCell,costFlowCell,readsPerFlowcell,
                                               ref.study,ref.study.name,
                                               samplesPerLane,
                                               read.umi.fit,gamma.mixed.fits,
                                               gamma.probs, disp.fun.param,
                                               nSamplesRange=NULL,
                                               nCellsRange=NULL, readDepthRange=NULL,
                                               mappingEfficiency=0.8,
                                               multipletRate=7.67e-06,multipletFactor=1.82,
                                               min.UMI.counts=10,perc.indiv.expr=0.5,
                                               nGenes=21000,samplingMethod="quantiles",
                                               multipletRateGrowth="linear"){

  #Check that exactly two of the parameters are set and the third one is not defined
  if(sum(is.null(nSamplesRange),is.null(nCellsRange),is.null(readDepthRange))!=1){
    stop(paste("To optimize the experimental design for a given budget,",
               "always exactly one of the parameters nSamplesRange, nCellsRange",
               "and readDepthRange should be set to NULL."))
  }

  #Case 1: estimate the sample size
  if(is.null(nSamplesRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nCellsRange,readDepthRange)
    colnames(param.combis)<-c("nCells","readDepth")

    #Sample size dependent on the budget
    param.combis$nSamples<-sapply(1:nrow(param.combis),
                                  function(i)floor(sampleSizeBudgetCalculation.libPrepCell(param.combis$nCells[i],
                                                                                           param.combis$readDepth[i],
                                                                                           totalBudget,prepCostCell,
                                                                                           costFlowCell,readsPerFlowcell)))
    #Case 2: estimate the number of cells per persons
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                function(i)cellsPersonBudgetCalculation.libPrepCell(param.combis$nSamples[i],
                                                                                    param.combis$readDepth[i],
                                                                                    totalBudget, prepCostCell,
                                                                                    costFlowCell,readsPerFlowcell))
    # Case 3: estimate the read depth
  } else {
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,nCellsRange)
    colnames(param.combis)<-c("nSamples","nCells")

    param.combis$readDepth<-sapply(1:nrow(param.combis),
                                   function(i)readDepthBudgetCalculation.libPrepCell(param.combis$nSamples[i],
                                                                                     param.combis$nCells[i],
                                                                                     totalBudget,prepCostCell,
                                                                                     costFlowCell,readsPerFlowcell))
  }

  #Remove all combinations where one of the parameters is <=0
  if(any(param.combis$nSamples==0) | any(param.combis$nCells<=0) | any(param.combis$readDepth<=0)){
    warning("Some of the parameter combinations are too expensive and removed from the parameter grid.")
    param.combis<-param.combis[param.combis$nSamples>0 & param.combis$nCells>0 & param.combis$readDepth>0,]
  }

  power.study<-mapply(power.general.withDoublets,
                      param.combis$nSamples,
                      param.combis$nCells,
                      param.combis$readDepth,
                      MoreArgs=list(ct.freq=ct.freq,
                                    multipletRate=multipletRate,
                                    multipletFactor=multipletFactor,
                                    type=type,
                                    ref.study=ref.study,
                                    ref.study.name=ref.study.name,
                                    samplesPerLane=samplesPerLane,
                                    read.umi.fit=read.umi.fit,
                                    gamma.mixed.fits=gamma.mixed.fits,
                                    gamma.probs=gamma.probs,
                                    ct=ct,
                                    disp.fun.param=disp.fun.param,
                                    mappingEfficiency=mappingEfficiency,
                                    min.UMI.counts=min.UMI.counts,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    multipletRateGrowth=multipletRateGrowth))

  power.study<-data.frame(apply(power.study,1,unlist),stringsAsFactors = FALSE)
  power.study[,2:ncol(power.study)]<-apply(power.study[,2:ncol(power.study)],2,as.numeric)

  return(power.study)
}

#' Optimizing cost parameters to maximize detection power for a given budget and 10X design
#'
#' This function determines the optimal parameter combination for a given budget.
#' The optimal combination is thereby the one with the highest detection power.
#' Of the three parameters sample size, cells per sample and read depth, two need to be set and
#' the third one is uniquely defined given the other two parameters and the overall budget.
#'
#' @param totalBudget Overall experimental budget
#' @param type (eqtl/de) study
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param ct.freq Frequency of the cell type of interest
#' @param costKit Cost for one 10X kit
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param ref.study Data frame with reference studies to be used for effect sizes and ranks
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param cellsPerLane Maximal number of cells per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell
#' (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per person that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.UMI.counts Expression defition parameter: more than is number of UMI counts for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane, "linear" if overloading should be
#' modeled explicitly, otherwise "constant". The default value for the parameter multipletRate is matching the option "linear".
#'
#' @export
optimize.constant.budget.restrictedDoublets<-function(totalBudget,type,
                                                     ct,ct.freq,
                                                     costKit,costFlowCell,readsPerFlowcell,
                                                     ref.study,ref.study.name,
                                                     cellsPerLane,
                                                     read.umi.fit,gamma.mixed.fits,
                                                     gamma.probs,disp.fun.param,
                                                     nSamplesRange=NULL,
                                                     nCellsRange=NULL, readDepthRange=NULL,
                                                     mappingEfficiency=0.8,
                                                     multipletRate=7.67e-06,multipletFactor=1.82,
                                                     min.UMI.counts=10,perc.indiv.expr=0.5,
                                                     nGenes=21000,samplingMethod="quantiles",
                                                     multipletRateGrowth="linear"){

  #Check that exactly two of the parameters are set and the third one is not defined
  if(sum(is.null(nSamplesRange),is.null(nCellsRange),is.null(readDepthRange))!=1){
    stop(paste("To optimize the experimental design for a given budget,",
               "always exactly one of the parameters nSamplesRange, nCellsRange",
               "and readDepthRange should be set to NULL."))
  }

  #Case 1: estimate the sample size
  if(is.null(nSamplesRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nCellsRange,readDepthRange)
    colnames(param.combis)<-c("nCells","readDepth")

    #Sample size dependent on the budget
    param.combis$nSamples<-sapply(1:nrow(param.combis),
                                  function(i)floor(sampleSizeBudgetCalculation.restrictedDoublets(param.combis$nCells[i],
                                                                                                 param.combis$readDepth[i],
                                                                                                 totalBudget,
                                                                                                 costKit,cellsPerLane,
                                                                                                 costFlowCell,readsPerFlowcell)))
    #Case 2: estimate the number of cells per persons
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                function(i)cellsPersonBudgetCalculation.restrictedDoublets(param.combis$nSamples[i],
                                                                                          param.combis$readDepth[i],
                                                                                          totalBudget,
                                                                                          costKit,cellsPerLane,
                                                                                          costFlowCell,readsPerFlowcell))
    # Case 3: estimate the read depth
  } else {
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,nCellsRange)
    colnames(param.combis)<-c("nSamples","nCells")

    param.combis$readDepth<-sapply(1:nrow(param.combis),
                                   function(i)readDepthBudgetCalculation.restrictedDoublets(param.combis$nSamples[i],
                                                                                           param.combis$nCells[i],
                                                                                           totalBudget,
                                                                                           costKit,cellsPerLane,
                                                                                           costFlowCell,readsPerFlowcell))
  }


  #Remove all combinations where one of the parameters is <=0
  if(any(param.combis$nSamples==0) | any(param.combis$nCells<=0) | any(param.combis$readDepth<=0)){
    warning("Some of the parameter combinations are too expensive and removed from the parameter grid.")
    param.combis<-param.combis[param.combis$nSamples>0 & param.combis$nCells>0 & param.combis$readDepth>0,]
  }


  power.study<-mapply(power.general.restrictedDoublets,
                      param.combis$nSamples,
                      param.combis$nCells,
                      param.combis$readDepth,
                      MoreArgs=list(ct.freq=ct.freq,
                                    multipletRate=multipletRate,
                                    multipletFactor=multipletFactor,
                                    type=type,
                                    ref.study=ref.study,
                                    ref.study.name=ref.study.name,
                                    cellsPerLane=cellsPerLane,
                                    read.umi.fit=read.umi.fit,
                                    gamma.mixed.fits=gamma.mixed.fits,
                                    gamma.probs=gamma.probs,
                                    ct=ct,
                                    disp.fun.param=disp.fun.param,
                                    mappingEfficiency=mappingEfficiency,
                                    min.UMI.counts=min.UMI.counts,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    multipletRateGrowth=multipletRateGrowth))

  power.study<-data.frame(apply(power.study,1,unlist),stringsAsFactors = FALSE)
  power.study[,2:ncol(power.study)]<-apply(power.study[,2:ncol(power.study)],2,as.numeric)

  return(power.study)
}


#' Optimizing cost parameters to maximize detection power for a given budget and Smart-seq design
#'
#' This function determines the optimal parameter combination for a given budget.
#' The optimal combination is thereby the one with the highest detection power.
#' Of the three parameters sample size, cells per sample and read depth, two need to be set and
#' the third one is uniquely defined given the other two parameters and the overall budget.
#'
#' @param totalBudget Overall experimental budget
#' @param type (eqtl/de) study
#' @param ct Cell type of interest (name should be found in the gamma and dispersion models)
#' @param ct.freq Frequency of the cell type of interest
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param ref.study Data frame with reference studies to be used for effect sizes and ranks
#' (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#' (required columns: parameter, ct (cell type), intercept, meanReads (slope))
#' @param gamma.probs Data frame with probability parameter of the gamma distributions
#' (required columns:ct (cell type), pi.c1 (probability of component 1))
#' @param disp.linear.fit Function to fit the dispersion parameter dependent on the mean (parameter linear dependent on read depth)
#' (required columns: parameter, ct (cell type), intercept, meanReads (slope))
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per person that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletFraction Multiplet fraction in the experiment as a constant factor
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.norm.count Expression defition parameter: more than is number of counts per kilobase transcript for each
#' gene per person and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals that need to have this gene expressed
#' to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles or random sampling)
#'
#' @export
#'
optimize.constant.budget.smartseq<-function(totalBudget, type,
                                           ct,ct.freq,
                                           prepCostCell,costFlowCell,readsPerFlowcell,
                                           ref.study,ref.study.name,
                                           gamma.mixed.fits,gamma.probs,
                                           disp.linear.fit,
                                           nSamplesRange=NULL,
                                           nCellsRange=NULL, readDepthRange=NULL,
                                           mappingEfficiency=0.8,
                                           multipletFraction=0,multipletFactor=1.82,
                                           min.norm.count=10,perc.indiv.expr=0.5,
                                           nGenes=21000,samplingMethod="quantiles"){

  #Check that exactly two of the parameters are set and the third one is not defined
  if(sum(is.null(nSamplesRange),is.null(nCellsRange),is.null(readDepthRange))!=1){
    stop(paste("To optimize the experimental design for a given budget,",
               "always exactly one of the parameters nSamplesRange, nCellsRange",
               "and readDepthRange should be set to NULL."))
  }

  #Case 1: estimate the sample size
  if(is.null(nSamplesRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nCellsRange,readDepthRange)
    colnames(param.combis)<-c("nCells","readDepth")

    #Sample size dependent on the budget
    param.combis$nSamples<-sapply(1:nrow(param.combis),
                                  function(i)floor(sampleSizeBudgetCalculation.libPrepCell(param.combis$nCells[i],
                                                                                           param.combis$readDepth[i],
                                                                                           totalBudget,prepCostCell,
                                                                                           costFlowCell,readsPerFlowcell)))
    #Case 2: estimate the number of cells per persons
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                function(i)cellsPersonBudgetCalculation.libPrepCell(param.combis$nSamples[i],
                                                                                    param.combis$readDepth[i],
                                                                                    totalBudget, prepCostCell,
                                                                                    costFlowCell,readsPerFlowcell))
    # Case 3: estimate the read depth
  } else {
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,nCellsRange)
    colnames(param.combis)<-c("nSamples","nCells")

    param.combis$readDepth<-sapply(1:nrow(param.combis),
                                   function(i)readDepthBudgetCalculation.libPrepCell(param.combis$nSamples[i],
                                                                                     param.combis$nCells[i],
                                                                                     totalBudget,prepCostCell,
                                                                                     costFlowCell,readsPerFlowcell))
  }

  #Remove all combinations where one of the parameters is <=0
  if(any(param.combis$nSamples==0) | any(param.combis$nCells<=0) | any(param.combis$readDepth<=0)){
    warning("Some of the parameter combinations are too expensive and removed from the parameter grid.")
    param.combis<-param.combis[param.combis$nSamples>0 & param.combis$nCells>0 & param.combis$readDepth>0,]
  }

  power.study<-mapply(power.smartseq,
                      param.combis$nSamples,
                      param.combis$nCells,
                      param.combis$readDepth,
                      MoreArgs=list(ct.freq=ct.freq,
                                    multipletFraction=multipletFraction,
                                    multipletFactor=multipletFactor,
                                    type=type,
                                    ref.study=ref.study,
                                    ref.study.name=ref.study.name,
                                    gamma.mixed.fits=gamma.mixed.fits,
                                    gamma.probs=gamma.probs,
                                    ct=ct,
                                    disp.linear.fit=disp.linear.fit,
                                    mappingEfficiency=mappingEfficiency,
                                    min.norm.count=min.norm.count,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod))

  power.study<-data.frame(apply(power.study,1,unlist),stringsAsFactors = FALSE)
  power.study[,2:ncol(power.study)]<-apply(power.study[,2:ncol(power.study)],2,as.numeric)

  return(power.study)
}

#' Power calculation for an eQTL gene using the F-test
#'
#' This function calculates the power to detect an eQTL gene.
#' @param heritability Heritability of the trait
#' @param sig.level Significane threshold
#' @param nSamples Sample size
#'
#' @import pwr
#'
#' @return Power to detect the eQTL gene
power.eqtl<-function(heritability, sig.level, nSamples) {
  require(pwr)

  #A sample size larger than 2 is required for the power analysis
  if(nSamples<3){
    return(NA)
  }

  f2 <- heritability / (1 - heritability)
  df.num <- 1 ## dfs of the full model
  df.denom <- nSamples - df.num - 1 ## error dfs
  power<-pwr.f2.test(u=df.num, v=df.denom, f2=f2, sig.level=sig.level)$power
  return(power)
}

# Alternative calculation of power function
# power.eqtl <- function(heritability, sig.level, nSamples) {
#   ## determine the rejection area under the null model (standard normal)
#   reject <- qnorm(1 - sig.level)
#   ## determine the non-centrality paramter
#   z <- sqrt((nSamples * heritability) / (1 - heritability))
#   ## get the probability to be in the rejection area given that alternative is true
#   ## P(reject H0 | H1 true) = P(Z > reject | H1 true)
#   power <- pnorm(reject, lower.tail=FALSE, mean=z)
#   return(power)
# }

#' Power calculation for a DE gene
#'
#' This function calculates the power to detect an DE gene (comparsion of two groups 0 and 1)
#' by using the function power.nb.test of the package MKmisc.
#' @param nSamples.group0 Sample size for group 0
#' @param mu.grou0 Mean value of group 0
#' @param RR effect size of group 1 vs group 0
#' @param theta 1/dispersion parameter of the negative binomial fit
#' @param sig.level Significance threshold
#' @param approach Choose between three different methods implemented in the package for the power calculation (1,2,3)
#' @param ssize.ratio Sample size ratio between group 0 and 1
#'
#' @return Power to detect the DE gene
#'
#' @import MKmisc
#'
power.de<-function(nSamples.group0,mu.group0,RR,theta,sig.level,approach=3,ssize.ratio=1){

  require(MKmisc)

  calc<-power.nb.test(n=nSamples.group0,mu0=mu.group0,RR=RR, duration=1,theta=theta, ssize.ratio=ssize.ratio,
                      sig.level=sig.level,alternative="two.sided",approach=approach)
  return(calc$power)
}

#' Calculate total cost dependent on parameters for 10X design
#'
#' @param nSamples Number of samples
#' @param nCells Cells per person
#' @param readDepth Read depth per cell
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of persons sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param rounding Rounds up the number of used kits and flow cells
#' (which might give a more realistic estimation of costs)
#'
#' @return Total experimental cost dependent on the parameters
#'
#' @export
budgetCalculation<-function(nSamples,nCells,readDepth,
                            costKit,samplesPerLane,
                            costFlowCell,readsPerFlowcell,
                            rounding=FALSE){

  if(rounding){
    totalBudget<-ceiling(nSamples/(6*samplesPerLane))*costKit+
      ceiling(nSamples*nCells*readDepth/readsPerFlowcell)*costFlowCell
  } else {
    totalBudget<-nSamples/(6*samplesPerLane)*costKit+
      nSamples*nCells*readDepth/readsPerFlowcell*costFlowCell
  }
  return(totalBudget)

}

#' Estimate possible sample size depending on the total cost and the other parameters for 10X design
#'
#' A balanced design with two classes is assumed for the sample size calculation.
#'
#' @param nCells Cells per person
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of persons sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Number of samples that can be measured with this budget and other parameters
#'
#' @export
sampleSizeBudgetCalculation<-function(nCells,readDepth,totalCost,
                                      costKit,samplesPerLane,
                                      costFlowCell,readsPerFlowcell){

  #Estimate the maximal sample size dependent on the cost for the other parameter
  samples <- totalCost / (costKit/(6*samplesPerLane) +
                 nCells * readDepth / readsPerFlowcell * costFlowCell)

  #Return only even sample sizes (due to the balanced design)
  return(floor(samples / 2) * 2)
}

#' Estimate possible number of cells per person depending on the total cost
#' and the other parameters for 10X design
#'
#' @param nSamples Number of samples
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of persons sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Cells per person that can be measured with this budget and other parameters
#'
#' @export
cellsPersonBudgetCalculation<-function(nSamples,readDepth,totalCost,
                                     costKit,samplesPerLane,
                                     costFlowCell,readsPerFlowcell){

  nCells <- (totalCost - nSamples/(6*samplesPerLane)*costKit)/
    (nSamples*readDepth/readsPerFlowcell*costFlowCell)

  return(floor(nCells))
}

#' Estimate possible read depth depending on the total cost
#' and the other parameters for 10X design
#'
#' @param nSamples Number of samples
#' @param nCells Cells per person
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of persons sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Read depth that can be measured with this budget and other parameters
#'
#' @export
readDepthBudgetCalculation<-function(nSamples,nCells,totalCost,
                                      costKit,samplesPerLane,
                                      costFlowCell,readsPerFlowcell){

  readDepth <- (totalCost - nSamples/(6*samplesPerLane)*costKit)/
    (nSamples*nCells/readsPerFlowcell*costFlowCell)

  return(floor(readDepth))
}

#' Calculate total cost dependent on parameters (adaptation with cells per lane instead of samples)
#'
#' @param nSamples Number of samples
#' @param nCells Cells per person
#' @param readDepth Read depth per cell
#' @param costKit Cost for one 10X kit
#' @param cellsPerLane Number of cells sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param rounding Rounds up the number of used kits and flow cells
#' (which might give a more realistic estimation of costs)
#'
#' @return Total experimental cost dependent on the parameters
#'
#' @export
budgetCalculation.restrictedDoublets<-function(nSamples,nCells,readDepth,
                                               costKit,cellsPerLane,
                                               costFlowCell,readsPerFlowcell,
                                               rounding=FALSE){

  #Estimate persons per lane
  samplesPerLane<-floor(cellsPerLane/nCells)

  return(budgetCalculation(nSamples,nCells,readDepth,
                           costKit,samplesPerLane, costFlowCell,
                           readsPerFlowcell, rounding))

}

#' Estimate possible sample size depending on the total cost and the other parameters
#' (with a restricted number of cells per lane)
#'
#' Variant of sampleSizeBudgetCalculation, which restricts the cells per lane instead the persons
#' per lane.
#'
#' @param nCells Cells per person
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param cellsPerLane Number of cells sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Number of samples that can be measured with this budget and other parameters
#'
#' @export
sampleSizeBudgetCalculation.restrictedDoublets<-function(nCells,readDepth,totalCost,
                                                        costKit,cellsPerLane,
                                                        costFlowCell,readsPerFlowcell){

  #Estimate persons per lane
  samplesPerLane<-floor(cellsPerLane/nCells)

  #Calculate sample size dependent on the estimated number of persons per lane
  return(sampleSizeBudgetCalculation(nCells,readDepth,totalCost,costKit,samplesPerLane,
                               costFlowCell,readsPerFlowcell))
}

#' Estimate possible number of cells per person depending on the total cost
#' and the other parameters (with a restricted number of cells per lane)
#'
#' Approximation without rounding
#'
#' @param nSamples Number of samples
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param cellsPerLane Number of cells sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Cells per person that can be measured with this budget and other parameters
#'
#' @export
cellsPersonBudgetCalculation.restrictedDoublets<-function(nSamples,readDepth,totalCost,
                                                          costKit,cellsPerLane,
                                                          costFlowCell,readsPerFlowcell){

  #Estimate the maximal sample size dependent on the cost for the other parameter
  nCells <- totalCost / (nSamples*costKit/(6*cellsPerLane) +
                           nSamples * readDepth / readsPerFlowcell * costFlowCell)

  return(floor(nCells))

}

#' Estimate possible read depth depending on the total cost
#' and the other parameters (with a restricted number of cells per lane)
#'
#' @param nSamples Number of samples
#' @param nCells Cells per person
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param cellsPerLane Number of cells sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Read depth that can be measured with this budget and other parameters
#'
#' @export
readDepthBudgetCalculation.restrictedDoublets<-function(nSamples,nCells,totalCost,
                                                        costKit,cellsPerLane,
                                                        costFlowCell,readsPerFlowcell){

  #Estimate persons per lane
  samplesPerLane<-floor(cellsPerLane/nCells)

  return(readDepthBudgetCalculation(nSamples,nCells,totalCost,costKit,samplesPerLane,
                                     costFlowCell,readsPerFlowcell))
}


#' Calculate total cost dependent on parameters for 10X design
#' using library preparation costs per cell
#'
#' @param nSamples Number of samples
#' @param nCells Cells per person
#' @param readDepth Read depth per cell
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param rounding Rounds up the number of used kits and flow cells
#' (which might give a more realistic estimation of costs)
#'
#' @return Total experimental cost dependent on the parameters
#'
#' @export
budgetCalculation.libPrepCell<-function(nSamples,nCells,readDepth,
                                        prepCostCell,
                                        costFlowCell,readsPerFlowcell,
                                        rounding=FALSE){

  if(rounding){
    totalBudget<-prepCostCell*nSamples*nCells+
      ceiling(nSamples*nCells*readDepth/readsPerFlowcell)*costFlowCell
  } else {
    totalBudget<-prepCostCell*nSamples*nCells+
      nSamples*nCells*readDepth/readsPerFlowcell*costFlowCell
  }
  return(totalBudget)

}

#' Estimate possible sample size depending on the total cost,
#' using library preparation costs per cell
#'
#' A balanced design with two classes is assumed for the sample size calculation.
#'
#' @param nCells Cells per person
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Number of samples that can be measured with this budget and other parameters
#'
#' @export
sampleSizeBudgetCalculation.libPrepCell<-function(nCells,readDepth,totalCost,
                                                  prepCostCell,costFlowCell,readsPerFlowcell){

  #Estimate the maximal sample size dependent on the cost for the other parameter
  samples <- totalCost / (prepCostCell * nCells +
                            nCells * readDepth / readsPerFlowcell * costFlowCell)

  #Return only even sample sizes (due to the balanced design)
  return(floor(samples / 2) * 2)
}


#' Estimate possible number of cells per person depending on the total cost,
#' using library preparation costs per cell
#'
#' @param nSamples Number of samples
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Cells per person that can be measured with this budget and other parameters
#'
#' @export
cellsPersonBudgetCalculation.libPrepCell<-function(nSamples,readDepth,totalCost,
                                                   prepCostCell, costFlowCell,readsPerFlowcell){

  nCells <- totalCost / (prepCostCell * nSamples +
                           nSamples * readDepth / readsPerFlowcell * costFlowCell)

  return(floor(nCells))

}

#' Estimate possible read depth depending on the total cost,
#' using library preparation costs per cell
#'
#' @param nSamples Number of samples
#' @param nCells Cells per person
#' @param totalCost Experimental budget
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Read depth that can be measured with this budget and other parameters
#'
#' @export
readDepthBudgetCalculation.libPrepCell<-function(nSamples,nCells,totalCost,
                                                 prepCostCell,
                                                 costFlowCell,readsPerFlowcell){

  readDepth <- (totalCost - prepCostCell*nSamples*nCells)/
    (nSamples*nCells/readsPerFlowcell*costFlowCell)

  return(floor(readDepth))
}
