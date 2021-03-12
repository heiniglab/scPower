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
#' Calculate required number of cells per individual to detect at least min.num.cells of a cell type
#' with class frequency cell.type.frac for each of nSamples individual.
#' with a probability of prob.cutoff
#' @param prob.cutoff Target power to detect the cell type
#' @param min.num.cells Minimal number of the cells for the cell type of interest that should be detected
#' @param cell.type.frac Frequency of the cell type of interest
#' @param nSamples Sample size (number of individuals)
#'
#' @return Required number of cells per individual to reach the target power
#'
#' @export
number.cells.detect.celltype<-function(prob.cutoff,min.num.cells,cell.type.frac,nSamples){
  return(qnbinom(prob.cutoff^(1/nSamples),min.num.cells,cell.type.frac)+min.num.cells)
}

#' Power calculation for a DE/eQTL study  with 10X design (with a restricted number of individuals per lane)
#'
#' This function to calculate the detection power for a DE or eQTL study, given DE/eQTL genes from a reference study,
#' in a single cell 10X RNAseq study. The power depends on the cost determining parameter of sample size, number of cells
#' per individual and read depth.
#'
#' @param nCells Number of cells per individual
#' @param readDepth Target read depth per cell
#' @param ct.freq Frequency of the cell type of interest
#' @param samplesPerLane Maximal number of individuals per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell
#'        depending on the mean readds per cell (required columns: intercept, reads (slope))
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#'        (required columns: parameter, ct (cell type), intercept, meanUMI (slope))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#'        (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome
#'        in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane,
#'        "linear" if overloading should be modeled explicitly, otherwise "constant".
#'        The default value for the parameter multipletRate is matching the option "linear".
#' @param returnResultsDetailed If true, return not only summary data frame,
#'        but additional list with exact probability vectors
#' @inheritParams calculate.probabilities
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
                                     samplesPerLane,read.umi.fit,
                                     gamma.mixed.fits,ct,
                                     disp.fun.param,
                                     mappingEfficiency=0.8,
                                     multipletRate=7.67e-06,multipletFactor=1.82,
                                     min.UMI.counts=3,perc.indiv.expr=0.5,
                                     nGenes=21000,samplingMethod="quantiles",
                                     multipletRateGrowth="linear",
                                     sign.threshold=0.05,
                                     MTmethod="Bonferroni",
                                     useSimulatedPower=TRUE,
                                     simThreshold=4,
                                     speedPowerCalc=FALSE,
                                     returnResultsDetailed=FALSE){

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
    stop("Too many cells per individual entered! Multiplet rate of more than 100%!")
  }

  usableCells<-round((1-multipletFraction)*nCells)
  #Estimate multiplet rate and "real read depth"
  readDepthSinglet<-readDepth*nCells/(usableCells+multipletFactor*(nCells-usableCells))

  #Estimate fraction of correctly mapped reads
  mappedReadDepth<-readDepthSinglet*mappingEfficiency

  #Get the fraction of cell type cells
  ctCells<-round(usableCells*ct.freq)

  #Check that really only one row is given for read.umi.fit
  if(nrow(read.umi.fit)>1){
    stop("Please only enter data frame with one row for read.umi.fit,
         only one fit can be evaluated in one run!")
  }

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

  gamma.parameters<-data.frame(p1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="p1"],
                               p3=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="p3"],
                               mean1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mean1"],
                               mean2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mean2"],
                               sd1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd1"],
                               sd2=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd2"])

  #Convert them to the rateshape variant
  gamma.parameters<-convert.gamma.parameters(gamma.parameters,type="rateshape")

  #Check if dispersion data for the cell type exists
  if(! any(disp.fun.param$ct==ct)){
    stop(paste("No dispersion fitting data in the data frame disp.fun.param fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right disp.fun.param object used."))
  }

  #Fit dispersion parameter dependent on mean parameter
  disp.fun<-disp.fun.param[disp.fun.param$ct==ct,]

  #Calculate expression power depending on gamma parameters
  probs<-calculate.probabilities(nSamples,ctCells,type,
                                 ref.study,ref.study.name,
                                 gamma.parameters,disp.fun,
                                 min.UMI.counts,perc.indiv.expr,
                                 nGenes,samplingMethod,
                                 sign.threshold,MTmethod,
                                 useSimulatedPower,
                                 simThreshold,
                                 speedPowerCalc,
                                 returnResultsDetailed)

  #Return either detailed probabilities for each DE/eQTL gene or only overview
  if(returnResultsDetailed){
    power.study<-data.frame(name=ref.study.name,
                            powerDetect=probs$overview.df$powerDetect,
                            exp.probs=probs$overview.df$exp.probs,
                            power=probs$overview.df$power,
                            sampleSize=nSamples,
                            totalCells=nCells,
                            usableCells=usableCells,
                            multipletFraction=multipletFraction,
                            ctCells=ctCells,
                            readDepth=readDepth,
                            readDepthSinglet=readDepthSinglet,
                            mappedReadDepth=mappedReadDepth,
                            expressedGenes=probs$overview.df$expressedGenes)

    return(list(overview.df=power.study,probs.df=probs$probs.df))

  } else {
    power.study<-data.frame(name=ref.study.name,
                            powerDetect=probs$powerDetect,
                            exp.probs=probs$exp.probs,
                            power=probs$power,
                            sampleSize=nSamples,
                            totalCells=nCells,
                            usableCells=usableCells,
                            multipletFraction=multipletFraction,
                            ctCells=ctCells,
                            readDepth=readDepth,
                            readDepthSinglet=readDepthSinglet,
                            mappedReadDepth=mappedReadDepth,
                            expressedGenes=probs$expressedGenes)
  }

  return(power.study)
}

#' Power calculation for a DE/eQTL study with 10X design (with a restricted number of cells per lane)
#'
#' This function is a variant of power.general.withDoublets, where not the number of samplesPerLane is given as an
#' parameter, but instead the individuals are distributed over the lanes in a way that restricts the total number of
#' cells per lane instead. This gives also an upper bound for the doublet rate.
#'
#' @param cellsPerLane Maximal number of cells per 10X lane
#' @inheritParams power.general.withDoublets
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
power.general.restrictedDoublets<-function(nSamples,nCells,readDepth,ct.freq,
                                           type,ref.study,ref.study.name,
                                           cellsPerLane,read.umi.fit,
                                           gamma.mixed.fits,ct,
                                           disp.fun.param,
                                           mappingEfficiency=0.8,
                                           multipletRate=7.67e-06,multipletFactor=1.82,
                                           min.UMI.counts=3,perc.indiv.expr=0.5,
                                           nGenes=21000,samplingMethod="quantiles",
                                           multipletRateGrowth="linear",
                                           sign.threshold=0.05,
                                           MTmethod="Bonferroni",
                                           useSimulatedPower=TRUE,
                                           simThreshold=4,
                                           speedPowerCalc=FALSE,
                                           returnResultsDetailed=FALSE){

  #Distribute individuals most optimal over the lanes
  samplesPerLane<-floor(cellsPerLane/nCells)

  if(samplesPerLane==0){
    stop("Allowed number of cells per lane is too low to fit so many cells per individual!")
  }

  return(power.general.withDoublets(nSamples,nCells,readDepth,ct.freq,
                             type,ref.study,ref.study.name,
                             samplesPerLane,
                             read.umi.fit,gamma.mixed.fits,ct,
                             disp.fun.param,
                             mappingEfficiency,
                             multipletRate,multipletFactor,
                             min.UMI.counts,perc.indiv.expr,
                             nGenes,samplingMethod,
                             multipletRateGrowth,
                             sign.threshold,MTmethod,
                             useSimulatedPower,
                             simThreshold,
                             speedPowerCalc,
                             returnResultsDetailed))
}

#' Power calculation for a DE/eQTL study with Smart-seq design
#'
#' This function to calculate the detection power for a DE or eQTL study, given DE/eQTL genes from a reference study,
#' in a single cell Smart-seq RNAseq study. The power depends on the cost determining parameter of sample size, number of cells
#' per individual and read depth.
#'
#' @param nSamples Sample size
#' @param nCells Number of cells per individual
#' @param readDepth Target read depth per cell
#' @param ct.freq Frequency of the cell type of interest
#' @param type (eqtl/de) study
#' @param ref.study Data frame with reference studies to be used for expression ranks and effect sizes
#'        (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type
#'        (required columns: parameter, ct (cell type), intercept, meanReads (slope))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.linear.fit Function to fit the dispersion parameter dependent on the mean (parameter linear dependent on read depth)
#'        (required columns: parameter, ct (cell type), intercept, meanReads (slope))
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome
#'        in the end (need to be between 1-0)
#' @param multipletFraction Multiplet fraction in the experiment as a constant factor
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param min.norm.count Expression defition parameter: more than is number of
#'        counts per kilobase transcript for each gene per individual and cell
#'        type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals
#'        that need to have this gene expressed to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking quantiles
#'        or random sampling)
#' @param sign.threshold Significance threshold
#' @param MTmethod Multiple testing correction method (possible options: "Bonferroni","FDR","none")
#' @param useSimulatedPower Option to simulate eQTL power for small mean values to increase accuracy
#'        (only possible for eQTL analysis)
#' @param simThreshold Threshold until which the simulated power is taken instead of the analytic
#' @param speedPowerCalc Option to speed power calculation by skipping all genes with
#'        an expression probability less than 0.01 (as overall power is anyway close to 0)
#' @param returnResultsDetailed If true, return not only summary data frame, but additional list with exact probability vectors
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
power.smartseq<-function(nSamples,nCells,readDepth,ct.freq,
                         type,ref.study,ref.study.name,
                         gamma.mixed.fits,ct,
                         disp.linear.fit,
                         mappingEfficiency=0.8,
                         multipletFraction=0,multipletFactor=1.82,
                         min.norm.count=3,perc.indiv.expr=0.5,
                         nGenes=21000,samplingMethod="quantiles",
                         sign.threshold=0.05,MTmethod="Bonferroni",
                         useSimulatedPower=TRUE,
                         simThreshold=4,
                         speedPowerCalc=FALSE,
                         returnResultsDetailed=FALSE){

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

  if(any(gamma.fits.ct$fitted.value[gamma.fits.ct$parameter %in%
                                    c("mean1","mean2","sd1","sd2")]<=0)){
    stop("At least one of the gamma parameter got negative for this read depth.",
         "Choose a higher read depth or a different gamma - read fit.")
  }

  gamma.parameters<-data.frame(p1=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="p1"],
                               p3=gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="p3"],
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

  #Calculate for each gene the expression probability
  sim.genes$exp.probs<-estimate.exp.prob.values(sim.genes$mean,1/sim.genes$disp,ctCells,
                                            nSamples=nSamples,min.counts=min.norm.count,
                                            perc.indiv=perc.indiv.expr)

  #Calculate the expected number of expressed genes
  exp.genes<-round(sum(sim.genes$exp.probs))

  #Set the simulated DE genes as the genes at the same rank position as the original DE genes
  ranks<-ref.study$rank[ref.study$name==ref.study.name]

  #Set all DE with rank > nGEnes to nGenes (expression anyway nearly 0)
  ranks[ranks>nGenes]<-nGenes

  #Choose the DE genes according to the rank
  foundSignGenes<-sim.genes[ranks,]

  #Calculate alpha parameter corrected for multiple testing
  if(MTmethod=="Bonferroni"){
    alpha<-sign.threshold/exp.genes
  } else if (MTmethod=="none"){
    alpha<-sign.threshold
    #For FDR correction, optimization is done dependent on eqtl/de power later
    #Only first parameters are calculated here
  } else if(MTmethod=="FDR"){
    lowerBound<-sign.threshold/exp.genes
    m0<-exp.genes-round(sum(foundSignGenes$exp.probs))
  } else {
    stop(paste("MTmethod",MTmethod,"is unknown! Please choose between",
               "Bonferroni, FDR and none!"))
  }

  #Calculate power
  if(type=="eqtl"){

    #Check that the required column Rsq exists
    if(! any(colnames(ref.study)=="Rsq")){
      stop(paste("Column name Rsq missing in the ref.study data frame.",
                 "Please make sure to provide this column for eQTL power analysis."))
    }

    #Set the Rsq respectively
    foundSignGenes$Rsq<-ref.study$Rsq[ref.study$name==ref.study.name]

    if(MTmethod=="FDR"){
      #In the extreme case that also with Bonferroni cut-off less than one TP
      #can be found, the optimization is not working, use here the Bonferroni
      #cutoff instead
      if(fdr.optimization(x=lowerBound,
                          fdr=sign.threshold,m0=m0,type=type,
                          exp.vector=foundSignGenes$exp.probs,
                          es.vector=foundSignGenes$Rsq,
                          nSamples=nSamples,
                          mean.vector=foundSignGenes$mean.length.sum,
                          useSimulatedPower=useSimulatedPower,
                          simThreshold=simThreshold)>0){
        alpha<-lowerBound
      } else {
        root<-uniroot(f=fdr.optimization,
                      interval=c(lowerBound,sign.threshold),
                      fdr=sign.threshold,m0=m0,type=type,
                      exp.vector=foundSignGenes$exp.probs,
                      es.vector=foundSignGenes$Rsq,
                      nSamples=nSamples,
                      mean.vector=foundSignGenes$mean.length.sum,
                      useSimulatedPower=useSimulatedPower,
                      simThreshold=simThreshold)

        alpha<-root$root
      }
    } else if (MTmethod=="Bonferroni"){
      #Restrict the Bonferroni for eQTLs cut-off further, assuming 10 independent SNPs per gene
      alpha<-alpha/10
    }

    #Skip power calculation for not expressed genes (if this option is chosen)
    if(speedPowerCalc){
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes), function(i)
        if(foundSignGenes$exp.probs[i]<0.01) {
          return(0)
        }else{
          power.eqtl(foundSignGenes$mean.length.sum[i],
                     foundSignGenes$Rsq[i],
                     alpha,nSamples,
                     useSimulatedPower=useSimulatedPower,
                     simThreshold=simThreshold)})
    } else {
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes),
         function(i) power.eqtl(foundSignGenes$mean.length.sum[i],
                                foundSignGenes$Rsq[i],
                                alpha,nSamples,
                                useSimulatedPower=useSimulatedPower,
                                simThreshold=simThreshold))
    }

  } else if (type=="de") {

    #Check that the required column FoldChange exists
    if(! any(colnames(ref.study)=="FoldChange")){
      stop(paste("Column name FoldChange missing in the ref.study data frame.",
                 "Please make sure to provide this column for DE power analysis."))
    }

    #Set the fold change respectively
    foundSignGenes$FoldChange<-ref.study$FoldChange[ref.study$name==ref.study.name]

    if(MTmethod=="FDR"){
      #In the extreme case that also with Bonferroni cut-off less than one TP
      #can be found, the optimization is not working, use here the Bonferroni
      #cutoff instead
      if(fdr.optimization(x=lowerBound,
                          fdr=sign.threshold,m0=m0,type=type,
                          exp.vector=foundSignGenes$exp.probs,
                          es.vector=foundSignGenes$FoldChange,
                          nSamples=nSamples,
                          mean.vector=foundSignGenes$mean.length.sum,
                          disp.vector = foundSignGenes$disp.sum)>0){
        alpha<-lowerBound
      } else {
        root<-uniroot(f=fdr.optimization,
                      interval=c(lowerBound,sign.threshold),
                      fdr=sign.threshold,m0=m0,type=type,
                      exp.vector=foundSignGenes$exp.probs,
                      es.vector=foundSignGenes$FoldChange,
                      nSamples=nSamples,
                      mean.vector=foundSignGenes$mean.length.sum,
                      disp.vector = foundSignGenes$disp.sum)

        alpha<-root$root
      }
    }

    #Skip power calculation for not expressed genes (if this option is chosen)
    if(speedPowerCalc){
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i)
        if(foundSignGenes$exp.probs[i]<0.01){
          return(0)
        } else {
          foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.de(
            floor(nSamples/2),
            foundSignGenes$mean.length.sum[i],
            foundSignGenes$FoldChange[i],
            1/foundSignGenes$disp.sum[i],
            alpha,3,ssize.ratio=1))
        })
    } else {
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.de(
        floor(nSamples/2),
        foundSignGenes$mean.length.sum[i],
        foundSignGenes$FoldChange[i],
        1/foundSignGenes$disp.sum[i],
        alpha,3,ssize.ratio=1))
    }

  } else  {
    stop('For type parameter only "eqtl" or "de" possible!')
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

  #Return either detailed probabilities for each DE/eQTL gene or only overview
  if(returnResultsDetailed){
    return(list(overview.df=power.study,probs.df=foundSignGenes))
  } else {
    return(power.study)
  }
}

#' Power calculation for a DE/eQTL study with same read depth as the fitted gamma distribution
#'
#' This is a simplified version of the function power.general.withDoublets to be used on a gamma
#' fit not parameterized for UMI/read counts. It evaluates the effect of different samples sizes
#' and cells per person, keeping the same read depth as in the experiment used for fitting.
#'
#' @param nCells Number of cells per individual
#' @param ct.freq Frequency of the cell type of interest
#' @param samplesPerLane Maximal number of individuals per 10X lane
#' @param gamma.parameters Data frame with gamma parameters for each cell type
#' (required columns: ct (cell type), s1, r1, s2, r2, p1, p2/p3 (gamma parameters for both components))
#' @param ct Cell type of interest (name from the gamma mixed models)
#' @param disp.fun.param Function to fit the dispersion parameter dependent on the mean
#' (required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))
#' @param mappingEfficiency Fraction of reads successfully mapped to the transcriptome in the end (need to be between 1-0)
#' @param multipletRate Expected increase in multiplets for additional cell in the lane
#' @param multipletFactor Expected read proportion of multiplet cells vs singlet cells
#' @param multipletRateGrowth Development of multiplet rate with increasing number of cells per lane, "linear" if overloading should be
#' modeled explicitly, otherwise "constant". The default value for the parameter multipletRate is matching the option "linear".
#' @inheritParams calculate.probabilities
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
#'
power.sameReadDepth.withDoublets<-function(nSamples,nCells,ct.freq,
                              type,ref.study,ref.study.name,
                              samplesPerLane,
                              gamma.parameters,ct,
                              disp.fun.param,
                              mappingEfficiency=0.8,
                              multipletRate=7.67e-06,multipletFactor=1.82,
                              min.UMI.counts=3,perc.indiv.expr=0.5,
                              nGenes=21000,samplingMethod="quantiles",
                              multipletRateGrowth="linear",
                              sign.threshold=0.05, MTmethod="Bonferroni",
                              useSimulatedPower=TRUE,
                              simThreshold=4,
                              speedPowerCalc=FALSE,
                              returnResultsDetailed=FALSE){


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
    stop("Too many cells per individual entered! Multiplet rate of more than 100%!")
  }

  usableCells<-round((1-multipletFraction)*nCells)

  #Get the fraction of cell type cells
  ctCells<-round(usableCells*ct.freq)

  #Check if gamma fit data for the cell type exists
  if(! any(gamma.parameters$ct==ct)){
    stop(paste("No gene curve fitting data in the data frame gamma.mixed.fits fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right gamma.mixed.fits object used."))
  }

  gamma.parameters<-gamma.parameters[gamma.parameters$ct==ct,]

  #Check if dispersion data for the cell type exists
  if(! any(disp.fun.param$ct==ct)){
    stop(paste("No dispersion fitting data in the data frame disp.fun.param fits to the specified cell type",
               ct,". Check that the cell type is correctly spelled and the right disp.fun.param object used."))
  }

  #Fit dispersion parameter dependent on mean parameter
  disp.fun<-disp.fun.param[disp.fun.param$ct==ct,]

  #Calculate expression power depending on gamma parameters
  probs<-calculate.probabilities(nSamples,ctCells,type,
                                 ref.study,ref.study.name,
                                 gamma.parameters,disp.fun,
                                 min.UMI.counts,perc.indiv.expr,
                                 nGenes,samplingMethod,
                                 sign.threshold,MTmethod,
                                 useSimulatedPower,
                                 simThreshold,
                                 speedPowerCalc,
                                 returnResultsDetailed)

  #Return either detailed probabilities for each DE/eQTL gene or only overview
  if(returnResultsDetailed){
    power.study<-data.frame(name=ref.study.name,
                            powerDetect=probs$overview.df$powerDetect,
                            exp.probs=probs$overview.df$exp.probs,
                            power=probs$overview.df$power,
                            sampleSize=nSamples,
                            totalCells=nCells,
                            usableCells=usableCells,
                            multipletFraction=multipletFraction,
                            ctCells=ctCells,
                            expressedGenes=probs$overview.df$expressedGenes)

    return(list(overview.df=power.study,probs.df=probs$probs.df))

  } else {
    power.study<-data.frame(name=ref.study.name,
                            powerDetect=probs$powerDetect,
                            exp.probs=probs$exp.probs,
                            power=probs$power,
                            sampleSize=nSamples,
                            totalCells=nCells,
                            usableCells=usableCells,
                            multipletFraction=multipletFraction,
                            ctCells=ctCells,
                            expressedGenes=probs$expressedGenes)
  }

  return(power.study)
}

#' Power calculation for a DE/eQTL study with same read depth as the fitted gamma distribution
#'
#' This function is a variant of power.sameReadDepth.withDoublets, where not the number of samplesPerLane is given as an
#' parameter, but instead the individuals are distributed over the lanes in a way that restricts the total number of
#' cells per lane instead. This gives also an upper bound for the doublet rate.
#' @param cellsPerLane Maximal number of cells per 10X lane
#' @inheritParams power.sameReadDepth.withDoublets
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
#'
power.sameReadDepth.restrictedDoublets<-function(nSamples,nCells,ct.freq,
                              type,ref.study,ref.study.name,
                              cellsPerLane,
                              gamma.parameters,ct,
                              disp.fun.param,
                              mappingEfficiency=0.8,
                              multipletRate=7.67e-06,multipletFactor=1.82,
                              min.UMI.counts=3,perc.indiv.expr=0.5,
                              nGenes=21000,samplingMethod="quantiles",
                              multipletRateGrowth="linear",
                              sign.threshold=0.05, MTmethod="Bonferroni",
                              useSimulatedPower=TRUE,
                              simThreshold=4,
                              speedPowerCalc=FALSE,
                              returnResultsDetailed=FALSE){

  #Distribute individuals most optimal over the lanes
  samplesPerLane<-floor(cellsPerLane/nCells)

  if(samplesPerLane==0){
    stop("Allowed number of cells per lane is too low to fit so many cells per individual!")
  }

  return(power.sameReadDepth.withDoublets(nSamples,nCells,ct.freq,
                                          type,ref.study,ref.study.name,
                                          samplesPerLane,
                                          gamma.parameters,ct,
                                          disp.fun.param,
                                          mappingEfficiency,
                                          multipletRate,multipletFactor,
                                          min.UMI.counts,perc.indiv.expr,
                                          nGenes,samplingMethod,
                                          multipletRateGrowth,
                                          sign.threshold,
                                          MTmethod,
                                          useSimulatedPower,
                                          simThreshold,
                                          speedPowerCalc,
                                          returnResultsDetailed))
}

#' Help function to calculate expression probability, power and detection power
#' for a given gamma distribution plus additional parameters
#'
#' @param nSamples Sample size
#' @param ctCells Number of cells of the target cell type
#' @param type (eqtl/de) study
#' @param ref.study Data frame with reference studies to be used for expression ranks and effect sizes
#'        (required columns: name (study name), rank (expression rank), FoldChange (DE study) /Rsq (eQTL study))
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param gamma.parameters Data frame with gamma parameters, filtered for the correct cell type
#'        (required columns: ct (cell type), s1, r1, s2, r2, p1, p2/p3
#'        (gamma parameters for both components))
#' @param disp.fun Function to fit the dispersion parameter dependent on the mean,
#'        filtered for the correct cell type (required columns: ct (cell type),
#'        asymptDisp, extraPois (both from taken from DEseq))
#' @param min.UMI.counts Expression defition parameter: more than is number of UMI counts for each
#'        gene per individual and cell type is required to defined it as expressed in one individual
#' @param perc.indiv.expr Expression defition parameter: percentage of individuals
#'        that need to have this gene expressed to define it as globally expressed
#' @param nGenes Number of genes to simulate (should match the number of genes used for the fitting)
#' @param samplingMethod Approach to sample the gene mean values (either taking
#'        quantiles or random sampling)
#' @param sign.threshold Significance threshold
#' @param MTmethod Multiple testing correction method
#'        (possible options: "Bonferroni","FDR","none")
#' @param useSimulatedPower Option to simulate eQTL power for small mean values
#'        to increase accuracy (only possible for eQTL analysis)
#' @param simThreshold Threshold until which the simulated power is taken instead
#'        of the analytic (only for the eQTL analysis)
#' @param speedPowerCalc Option to speed power calculation by skipping all genes with
#'        an expression probability less than 0.01 (as overall power is anyway close to 0)
#' @param returnResultsDetailed If true, return not only summary data frame,
#'        but additional list with exact probability vectors
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
calculate.probabilities<-function(nSamples,ctCells,type,
                                  ref.study,ref.study.name,
                                  gamma.parameters,disp.fun,
                                  min.UMI.counts,perc.indiv.expr,
                                  nGenes,samplingMethod,
                                  sign.threshold,MTmethod,
                                  useSimulatedPower,
                                  simThreshold,
                                  speedPowerCalc,
                                  returnResultsDetailed){

  #Sample means values
  if(samplingMethod=="random"){
    gene.means<-sample.mean.values.random(gamma.parameters,nGenes)
  } else if (samplingMethod=="quantiles"){
    gene.means<-sample.mean.values.quantiles(gamma.parameters,nGenes)
  } else {
    stop("No known sampling method. Use the options 'random' or 'quantiles'.")
  }

  gene.disps<-sample.disp.values(gene.means,disp.fun)

  sim.genes<-data.frame(mean=gene.means, disp=gene.disps)

  #Sort simulated genes by mean expression
  sim.genes<-sim.genes[order(sim.genes$mean, decreasing = TRUE),]

  #Calculate the mean per cell type for each individuum
  sim.genes$mean.sum<-sim.genes$mean*ctCells
  sim.genes$disp.sum<-sim.genes$disp/ctCells

  #Calculate for each gene the expression probability
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

  #Set all DE with rank > nGenes to nGenes (expression anyway nearly 0)
  ranks[ranks>nGenes]<-nGenes

  #Choose the DE genes according to the rank
  foundSignGenes<-sim.genes[ranks,]

  #Calculate alpha parameter corrected for multiple testing
  if(MTmethod=="Bonferroni"){
    alpha<-sign.threshold/exp.genes
  } else if (MTmethod=="none"){
    alpha<-sign.threshold
  #For FDR correction, optimization is done dependent on eqtl/de power later
  #Only first parameters are calculated here
  } else if(MTmethod=="FDR"){
    lowerBound<-sign.threshold/exp.genes
    m0<-exp.genes-round(sum(foundSignGenes$exp.probs))
  } else {
    stop(paste("MTmethod",MTmethod,"is unknown! Please choose between",
               "Bonferroni, FDR and none!"))
  }

  #Calculate power
  if(type=="eqtl"){

    #Check that the required column Rsq exists
    if(! any(colnames(ref.study)=="Rsq")){
      stop(paste("Column name Rsq missing in the ref.study data frame.",
                 "Please make sure to provide this column for eQTL power analysis."))
    }

    #Set the Rsq respectively
    foundSignGenes$Rsq<-ref.study$Rsq[ref.study$name==ref.study.name]

    if(MTmethod=="FDR"){
      #In the extreme case that also with Bonferroni cut-off less than one TP
      #can be found, the optimization is not working, use here the Bonferroni
      #cutoff instead
      if(fdr.optimization(x=lowerBound,
                          fdr=sign.threshold,m0=m0,type=type,
                          exp.vector=foundSignGenes$exp.probs,
                          es.vector=foundSignGenes$Rsq,
                          nSamples=nSamples,
                          mean.vector=foundSignGenes$mean.sum,
                          useSimulatedPower=useSimulatedPower,
                          simThreshold=simThreshold)>0){
        alpha<-lowerBound
      } else {
        root<-uniroot(f=fdr.optimization,
                      interval=c(lowerBound,sign.threshold),
                      fdr=sign.threshold,m0=m0,type=type,
                      exp.vector=foundSignGenes$exp.probs,
                      es.vector=foundSignGenes$Rsq,
                      nSamples=nSamples,
                      mean.vector=foundSignGenes$mean.sum,
                      useSimulatedPower=useSimulatedPower,
                      simThreshold=simThreshold)

        alpha<-root$root
      }
    } else if (MTmethod=="Bonferroni"){
      #Restrict the Bonferroni for eQTLs cut-off further, assuming 10 independent SNPs per gene
      alpha<-alpha/10
    }

    #Skip power calculation for not expressed genes (if this option is chosen)
    if(speedPowerCalc){
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes), function(i)
                                   if(foundSignGenes$exp.probs[i]<0.01) {
                                     return(0)
                                     }else{
                                       power.eqtl(foundSignGenes$mean.sum[i],
                                                   foundSignGenes$Rsq[i],
                                                   alpha,nSamples,
                                                   useSimulatedPower,
                                                   simThreshold)})
    } else {
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes),
                                   function(i) power.eqtl(foundSignGenes$mean.sum[i],
                                                          foundSignGenes$Rsq[i],
                                                          alpha,nSamples,
                                                          useSimulatedPower,
                                                          simThreshold))
    }

  } else if (type=="de") {

    #Check that the required column FoldChange exists
    if(! any(colnames(ref.study)=="FoldChange")){
      stop(paste("Column name FoldChange missing in the ref.study data frame.",
                 "Please make sure to provide this column for DE power analysis."))
    }

    #Set the fold change respectively
    foundSignGenes$FoldChange<-ref.study$FoldChange[ref.study$name==ref.study.name]

    if(MTmethod=="FDR"){
      #In the extreme case that also with Bonferroni cut-off less than one TP
      #can be found, the optimization is not working, use here the Bonferroni
      #cutoff instead
      if(fdr.optimization(x=lowerBound,
                          fdr=sign.threshold,m0=m0,type=type,
                          exp.vector=foundSignGenes$exp.probs,
                          es.vector=foundSignGenes$FoldChange,
                          nSamples=nSamples,
                          mean.vector=foundSignGenes$mean.sum,
                          disp.vector = foundSignGenes$disp.sum)>0){
        alpha<-lowerBound
      } else {
        root<-uniroot(f=fdr.optimization,
                      interval=c(lowerBound,sign.threshold),
                      fdr=sign.threshold,m0=m0,type=type,
                      exp.vector=foundSignGenes$exp.probs,
                      es.vector=foundSignGenes$FoldChange,
                      nSamples=nSamples,
                      mean.vector=foundSignGenes$mean.sum,
                      disp.vector = foundSignGenes$disp.sum)

        alpha<-root$root
      }
    }

    #Skip power calculation for not expressed genes (if this option is chosen)
    if(speedPowerCalc){
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i)
          if(foundSignGenes$exp.probs[i]<0.01){
            return(0)
          } else {
            foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.de(
              floor(nSamples/2),
              foundSignGenes$mean.sum[i],
              foundSignGenes$FoldChange[i],
              1/foundSignGenes$disp.sum[i],
              alpha,3,ssize.ratio=1))
          })
    } else {
      foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.de(
        floor(nSamples/2),
        foundSignGenes$mean.sum[i],
        foundSignGenes$FoldChange[i],
        1/foundSignGenes$disp.sum[i],
        alpha,3,ssize.ratio=1))
    }

  } else  {
    stop('For type parameter only "eqtl" or "de" possible!')
  }

  #Calculate total probability as the DE power times the expression probability
  foundSignGenes$combined.prob<-foundSignGenes$power*foundSignGenes$exp.probs

  #Return probabilities and expected number of expressed genes
  results<-data.frame(powerDetect=mean(foundSignGenes$combined.prob),
                     exp.probs=mean(foundSignGenes$exp.probs),
                     power=mean(foundSignGenes$power),
                     expressedGenes=exp.genes)

  #Return either detailed probabilities for each DE/eQTL gene or only overview
  if(returnResultsDetailed){
    return(list(overview.df=results,probs.df=foundSignGenes))
  } else {
    return(results)
  }
}


#' Optimizing cost parameters to maximize detection power for a given budget and 10X design
#'
#' This function determines the optimal parameter combination for a given budget.
#' The optimal combination is thereby the one with the highest detection power.
#' Of the three parameters sample size, cells per sample and read depth, two need to be set and
#' the third one is uniquely defined given the other two parameters and the overall budget.
#'
#' @param totalBudget Overall experimental budget
#' @param costKit Cost for one 10X kit
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per individual that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @inheritParams power.general.withDoublets
#'
#' @return Data frame with overall detection power, power and expression power for
#' each possible parameter combination given the budget and the parameter ranges
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
                                   disp.fun.param,
                                   nSamplesRange=NULL,
                                   nCellsRange=NULL, readDepthRange=NULL,
                                   mappingEfficiency=0.8,
                                   multipletRate=7.67e-06,multipletFactor=1.82,
                                   min.UMI.counts=3,perc.indiv.expr=0.5,
                                   nGenes=21000,samplingMethod="quantiles",
                                   multipletRateGrowth="linear",
                                   sign.threshold=0.05,MTmethod="Bonferroni",
                                   useSimulatedPower=FALSE,
                                   simThreshold=4,
                                   speedPowerCalc=FALSE){

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
  #Case 2: estimate the number of cells per individuals
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                   function(i)cellsBudgetCalculation(param.combis$nSamples[i],
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

  #Check if at least one parameter combination remains
  if(nrow(param.combis)==0){
    stop("The total budget is too low for parameters in the given range!")
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
                                    ct=ct,
                                    disp.fun.param=disp.fun.param,
                                    mappingEfficiency=mappingEfficiency,
                                    min.UMI.counts=min.UMI.counts,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    multipletRateGrowth=multipletRateGrowth,
                                    sign.threshold=sign.threshold,
                                    MTmethod=MTmethod,
                                    useSimulatedPower=useSimulatedPower,
                                    simThreshold=simThreshold,
                                    speedPowerCalc=speedPowerCalc))

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
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per individual that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @inheritParams power.general.withDoublets
#'
#' @return Data frame with overall detection power, power and expression power for
#' each possible parameter combination given the budget and the parameter ranges
#'
#' @export
#'
optimize.constant.budget.libPrepCell<-function(totalBudget, type,
                                               ct,ct.freq,
                                               prepCostCell,costFlowCell,readsPerFlowcell,
                                               ref.study,ref.study.name,
                                               samplesPerLane,
                                               read.umi.fit,gamma.mixed.fits,
                                               disp.fun.param,
                                               nSamplesRange=NULL,
                                               nCellsRange=NULL, readDepthRange=NULL,
                                               mappingEfficiency=0.8,
                                               multipletRate=7.67e-06,multipletFactor=1.82,
                                               min.UMI.counts=3,perc.indiv.expr=0.5,
                                               nGenes=21000,samplingMethod="quantiles",
                                               multipletRateGrowth="linear",
                                               sign.threshold=0.05,MTmethod="Bonferroni",
                                               useSimulatedPower=FALSE,
                                               simThreshold=4,
                                               speedPowerCalc=FALSE){

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
    #Case 2: estimate the number of cells per individuals
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                function(i)cellsBudgetCalculation.libPrepCell(param.combis$nSamples[i],
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

  #Check if at least one parameter combination remains
  if(nrow(param.combis)==0){
    stop("The total budget is too low for parameters in the given range!")
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
                                    ct=ct,
                                    disp.fun.param=disp.fun.param,
                                    mappingEfficiency=mappingEfficiency,
                                    min.UMI.counts=min.UMI.counts,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    multipletRateGrowth=multipletRateGrowth,
                                    sign.threshold=sign.threshold,
                                    MTmethod=MTmethod,
                                    useSimulatedPower=useSimulatedPower,
                                    simThreshold=simThreshold,
                                    speedPowerCalc=speedPowerCalc))

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
#' @param costKit Cost for one 10X kit
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per individual that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @inheritParams power.general.restrictedDoublets
#'
#' @return Data frame with overall detection power, power and expression power for
#' each possible parameter combination given the budget and the parameter ranges
#'
#' @export
#'
optimize.constant.budget.restrictedDoublets<-function(totalBudget,type,
                                                     ct,ct.freq,
                                                     costKit,costFlowCell,readsPerFlowcell,
                                                     ref.study,ref.study.name,
                                                     cellsPerLane,
                                                     read.umi.fit,gamma.mixed.fits,
                                                     disp.fun.param,
                                                     nSamplesRange=NULL,
                                                     nCellsRange=NULL, readDepthRange=NULL,
                                                     mappingEfficiency=0.8,
                                                     multipletRate=7.67e-06,multipletFactor=1.82,
                                                     min.UMI.counts=3,perc.indiv.expr=0.5,
                                                     nGenes=21000,samplingMethod="quantiles",
                                                     multipletRateGrowth="linear",
                                                     sign.threshold=0.05,MTmethod="Bonferroni",
                                                     useSimulatedPower=FALSE,
                                                     simThreshold=4,
                                                     speedPowerCalc=FALSE){

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
    #Case 2: estimate the number of cells per individuals
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                function(i)cellsBudgetCalculation.restrictedDoublets(param.combis$nSamples[i],
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

  #Check if at least one parameter combination remains
  if(nrow(param.combis)==0){
    stop("The total budget is too low for parameters in the given range!")
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
                                    ct=ct,
                                    disp.fun.param=disp.fun.param,
                                    mappingEfficiency=mappingEfficiency,
                                    min.UMI.counts=min.UMI.counts,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    multipletRateGrowth=multipletRateGrowth,
                                    sign.threshold=sign.threshold,
                                    MTmethod=MTmethod,
                                    useSimulatedPower=useSimulatedPower,
                                    simThreshold=simThreshold,
                                    speedPowerCalc=speedPowerCalc))

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
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#' @param nSamplesRange Range of sample sizes that should be tested (vector)
#' @param nCellsRange Range of cells per individual that should be tested (vector)
#' @param readDepthRange Range of read depth values that should be tested (vector)
#' @inheritParams power.smartseq
#'
#' @return Data frame with overall detection power, power and expression power for
#' each possible parameter combination given the budget and the parameter ranges
#'
#' @export
#'
optimize.constant.budget.smartseq<-function(totalBudget, type,
                                           ct,ct.freq,
                                           prepCostCell,costFlowCell,readsPerFlowcell,
                                           ref.study,ref.study.name,
                                           gamma.mixed.fits,
                                           disp.linear.fit,
                                           nSamplesRange=NULL,
                                           nCellsRange=NULL, readDepthRange=NULL,
                                           mappingEfficiency=0.8,
                                           multipletFraction=0,multipletFactor=1.82,
                                           min.norm.count=3,perc.indiv.expr=0.5,
                                           nGenes=21000,samplingMethod="quantiles",
                                           sign.threshold=0.05,MTmethod="Bonferroni",
                                           useSimulatedPower=FALSE,
                                           simThreshold=4,
                                           speedPowerCalc=FALSE){

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
    #Case 2: estimate the number of cells per individuals
  } else if (is.null(nCellsRange)){
    #Build a frame of all possible combinations
    param.combis<-expand.grid(nSamplesRange,readDepthRange)
    colnames(param.combis)<-c("nSamples","readDepth")

    param.combis$nCells<-sapply(1:nrow(param.combis),
                                function(i)cellsBudgetCalculation.libPrepCell(param.combis$nSamples[i],
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

  #Check if at least one parameter combination remains
  if(nrow(param.combis)==0){
    stop("The total budget is too low for parameters in the given range!")
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
                                    ct=ct,
                                    disp.linear.fit=disp.linear.fit,
                                    mappingEfficiency=mappingEfficiency,
                                    min.norm.count=min.norm.count,
                                    perc.indiv.expr=perc.indiv.expr,
                                    nGenes=nGenes,
                                    samplingMethod=samplingMethod,
                                    sign.threshold=sign.threshold,
                                    MTmethod=MTmethod,
                                    useSimulatedPower=useSimulatedPower,
                                    simThreshold=simThreshold,
                                    speedPowerCalc=speedPowerCalc))

  power.study<-data.frame(apply(power.study,1,unlist),stringsAsFactors = FALSE)
  power.study[,2:ncol(power.study)]<-apply(power.study[,2:ncol(power.study)],2,as.numeric)

  return(power.study)
}

#' Wrapper funtion to use either simulated power or power based on F-test
#' (dependent on pseudobulk mean and used parameters)
#'
#' @param count.mean Expression mean in the pseudobulk
#' @param heritability Heritability of the trait
#' @param sig.level Significane threshold
#' @param nSamples Sample size
#' @param useSimulatedPower Option to simulate eQTL power for small mean values
#'        to increase accuracy
#' @param simThreshold Threshold until which the simulated power is taken
#'        instead of the analytic
#'
#' @return Power to detect the eQTL gene
#'
power.eqtl<-function(count.mean,heritability, sig.level, nSamples,
                     useSimulatedPower=TRUE, simThreshold=4){

  #Use a cut-off of 10 for power simulation!
  if(useSimulatedPower){
    if(round(count.mean) == 0){
      return(0)
    } else if (round(count.mean) > simThreshold){
      return(power.eqtl.ftest(heritability, sig.level, nSamples))
    } else {
      index<-paste(round(count.mean),round(heritability,2),nSamples,
                   sep="_")
      if(index %in% rownames(scPower::sim.eqtl.pvals)){
        return(mean(scPower::sim.eqtl.pvals[index,]<sig.level))
      } else {
        # warning(paste0("Simulated p-values not available for the current parameter combination ",
        #                "(",round(count.mean),",",round(heritability,2),",",nSamples,").",
        #                "Calculation from scratch might take a bit!"))
        return(power.eqtl.simulated(count.mean,heritability, sig.level, nSamples))
      }
    }

  } else {
    return(power.eqtl.ftest(heritability, sig.level, nSamples))
  }
}

#' Power calculation for an eQTL gene using the F-test
#'
#' This function calculates the power to detect an eQTL gene.
#' @param heritability Heritability of the trait
#' @param sig.level Significane threshold
#' @param nSamples Sample size
#'
#' @return Power to detect the eQTL gene
#'
power.eqtl.ftest<-function(heritability, sig.level, nSamples) {
  require(pwr)

  #A sample size larger than 2 is required for the power analysis
  if(nSamples<3){
    return(NA)
  }

  f2 <- heritability / (1 - heritability)
  df.num <- 1 ## dfs of the full model
  df.denom <- nSamples - df.num - 1 ## error dfs
  power<-pwr::pwr.f2.test(u=df.num, v=df.denom, f2=f2, sig.level=sig.level)$power
  return(power)
}

#' Power calculation for an eQTL gene using simulations
#'
#' The eQTL power is independent of the mean except for very small mean values,
#' where small effect sizes might not be detectable due to the discrete counts
#' In these cases, the power can instead be simulated.
#'
#' @param count.mean Expression mean in the pseudobulk
#' @param heritability Heritability of the trait
#' @param sig.level Significane threshold
#' @param nSamples Sample size
#' @param rep.times Number of repetitions used for the
#'
#' @return Power to detect the eQTL gene
#'
power.eqtl.simulated<-function(count.mean, heritability, sig.level, nSamples,
                               rep.times=100){

  #Use precalculated size estimates
  size.estimates<-scPower::size.estimates

  #Power for pseudobulk means close to 0 is set to 0 (no simulation possible)
  if(round(count.mean)==0){return(0)}

  #Simulated power
  p.vals<-sapply(1:rep.times,function(i)
    power.eqtl.simulated.help(round(count.mean), heritability, nSamples,
                              size.estimates))

  #Simulated power
  return(mean(p.vals<sig.level))

}

#' Helper function for eQTL simulation power calculation
#'
#' @param count.mean Expression mean in the pseudobulk
#' @param Rsq Heritability of the trait
#' @param nSamples Sample size
#' @param size.estimates Data frame with precalculated size values to speed calculation
#' @param af Allele frequency (sampled if not explicity given)
#'
power.eqtl.simulated.help<-function(count.mean,Rsq,nSamples,
                                    size.estimates,af=NULL){

  #Randomly sample AF if not given
  if(is.null(af)){
    af<-sample(seq(0.1,0.9,by=0.1),1)
  }

  #Genotype distribution in the population
  bb<-round(nSamples*af^2)
  ab<-round(nSamples*2*af*(1-af))
  aa<-nSamples-ab-bb
  genotypes<-c(rep(0,aa),rep(1,ab),rep(2,bb))

  #Calculate beta value and standard error
  beta<-sqrt(Rsq/(2*af*(1-af)))

  #Get dispersion parameter from look-up table if available
  #Look-up dispersion parameter from numerical optimization
  size.vector<-unlist(size.estimates[size.estimates$Rsq==round(Rsq,2) &
                                size.estimates$af==round(af,1) &
                                size.estimates$mean==round(count.mean),4:6])

  #If the look-up value could not be found, calculate it by hand
  if(length(size.vector)==0){
      # warning(paste("Look-up value not found for Rsq",round(Rsq,2),
      #               "af",round(af,1),"mean",round(count.mean),".",
      #               "Optimization will take time."))

      size.vector<-scPower:::estimate.size.simulation(round(count.mean),
                                                      round(Rsq,2),
                                                      round(af,1))
      #Replace size factors with NA with voom approximation!
      sd.error<-sqrt(1-Rsq)

      for(g in 1:3){
        size.vector[g]<-ifelse(is.na(size.vector[g]),
                       1/(sd.error^2-1/(exp(log(count.mean) + beta * (g-1)))),
                       size.vector[g])
      }

  }

  #Sample from a normal distribution dependent on the genotype
  mean.vector<-exp(log(count.mean) + beta * genotypes)
  size.vector<-size.vector[genotypes+1]

  #Suppress warning messages if NA values are generated
  suppressWarnings(counts<-rnbinom(nSamples,mu=mean.vector,size=size.vector))

  simData<-data.frame(count=counts,
                      log.count=log(counts+1),
                      g=genotypes)

  #Check if the counts are all the same
  if(length(setdiff(unique(simData$count),NA))<=1){
    return(1)
  }

  #Calculate p-values
  model<-lm(log.count~g,data=simData)

  #Check if there were enough data points to fit the model
  if("g" %in% row.names(summary(model)$coefficients)){
    return(summary(model)$coefficients["g","Pr(>|t|)"])
  } else{
    return(1)
  }

}


#' Estimation of size parameter (1/dsp) for eQTL power simulation
#'
#' @param count.mean Mean of negative binomial distribution for the counts
#' @param Rsq Effect size
#' @param af Allele frequency
#'
#' @return Numerical optimized size parameter (1/dispersion) to model the standard error
#' as good as possible in the simulation
#'
estimate.size.simulation<-function(count.mean,Rsq,af){

  #Calculate beta and Rsq
  beta<-sqrt(Rsq/(2*af*(1-af)))
  sd.error<-sqrt(1-Rsq)

  #Estimate for genotype the suitable size factor for each genotype
  disp.estimates<-c()
  for(g in c(0,1,2)){
    #Tayler estimate of size as approximation for optimization range
    beta.mean<-exp(log(count.mean) + beta * g)
    size<-1/(sd.error^2-1/beta.mean)

    root<-try(uniroot(f=optimizeSizeFactor, interval=c(0.01,max(size*2,1)),
                      sd.error=sd.error,
                      beta.mean=beta.mean),silent=TRUE)

    #For some very small size factors the optimization does not work
    if(class(root)=="try-error"){
      disp.estimates<-c(disp.estimates,NA)
    } else {
      disp.estimates<-c(disp.estimates,root$root)
    }
  }
  return(disp.estimates)
}


#' Function for numeric optimization of size parameter
#'
#' @param x Tested size factor
#' @param sd.error The targeted standard error
#' @param beta.mean The mean of the counts
#'
#' @return A value close to 0 shows a good accordance of targeted and simulated
#' standard error
optimizeSizeFactor<-function(x,sd.error,beta.mean){
  set.seed(1)
  #Simulate counts
  log.counts<-log(rnbinom(10000,mu=beta.mean,size=x)+1)
  return(sd.error-sd(log.counts))
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
#' by using the function power.nb.test of the package MKmisc. The power for a gene with mean of 0
#' is defined as 0 (the gene will have an expression probability of 0, so the overall detection power is also 0).
#'
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
power.de<-function(nSamples.group0,mu.group0,RR,theta,sig.level,approach=3,ssize.ratio=1){

  require(MKmisc)

  if(mu.group0 == 0){
    return(0)
  } else {
    calc<-MKmisc::power.nb.test(n=nSamples.group0,mu0=mu.group0,RR=RR, duration=1,theta=theta, ssize.ratio=ssize.ratio,
                        sig.level=sig.level,alternative="two.sided",approach=approach)
    return(calc$power)
  }
}

#' Function for numeric optimization to get an FDR corrected significance thresold
#'
#' @param x Adjusted significance threshold to be found during optimization
#' @param fdr Chosen FDR threshold
#' @param m0 Number of null hypothesis
#' @param type Either DE or eQTL power
#' @param exp.vector Vector of expression probabilities for each DEG/eQTL
#' @param es.vector Effect size vector (Fold Change for DEGs, Rsq for eQTLs)
#' @param nSamples Sample size
#' @param mean.vector Mean value for each DEG / eQTL gene
#' @param disp.vector Dispersion value for each DEG (only required for DE power)
#' @param useSimulatedPower Option to simulate eQTL power for small mean values
#'        to increase accuracy (only required for the eQTL power)
#' @param simThreshold Threshold until which the simulated power is taken
#'        instead of the analytic (only required with eQTL power)
#'
#' @return Optimizing this function to 0 will lead to the correct adjusted
#' significance threshold x
#'
fdr.optimization<-function(x,fdr,m0,type,
                           exp.vector,
                           es.vector,
                           nSamples,
                           mean.vector,
                           disp.vector=NULL,
                           useSimulatedPower=TRUE,
                           simThreshold=4){

  #Calculate DE power (similar for eQTL power)
  if(type=="de"){
    power<-sapply(1:length(es.vector), function(i)
      power.de(floor(nSamples/2),mean.vector[i],
               es.vector[i],1/disp.vector[i],x))
  } else if (type=="eqtl"){
    power<-sapply(1:length(es.vector), function(i)
      power.eqtl(mean.vector[i],es.vector[i],x,nSamples,
                 useSimulatedPower,simThreshold))
  } else {
    stop("Type unknown!")
  }

  r1<-sum(power*exp.vector, na.rm=TRUE)

  return(x-(fdr*r1)/(m0*(1-fdr)))
}


#' Calculate total cost dependent on parameters for 10X design
#'
#' @param nSamples Number of samples
#' @param nCells Cells per individual
#' @param readDepth Read depth per cell
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of individuals sequenced per lane
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
#' @param nCells Cells per individual
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of individuals sequenced per lane
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

#' Estimate possible number of cells per individual depending on the total cost
#' and the other parameters for 10X design
#'
#' @param nSamples Number of samples
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of individuals sequenced per lane
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Cells per individual that can be measured with this budget and other parameters
#'
#' @export
cellsBudgetCalculation<-function(nSamples,readDepth,totalCost,
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
#' @param nCells Cells per individual
#' @param totalCost Experimental budget
#' @param costKit Cost for one 10X kit
#' @param samplesPerLane Number of individuals sequenced per lane
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
#' @param nCells Cells per individual
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

  #Estimate individuals per lane
  samplesPerLane<-floor(cellsPerLane/nCells)

  return(budgetCalculation(nSamples,nCells,readDepth,
                           costKit,samplesPerLane, costFlowCell,
                           readsPerFlowcell, rounding))

}

#' Estimate possible sample size depending on the total cost and the other parameters
#' (with a restricted number of cells per lane)
#'
#' Variant of sampleSizeBudgetCalculation, which restricts the cells per lane instead the individuals
#' per lane.
#'
#' @param nCells Cells per individual
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

  #Estimate individuals per lane
  samplesPerLane<-floor(cellsPerLane/nCells)

  #Calculate sample size dependent on the estimated number of individuals per lane
  return(sampleSizeBudgetCalculation(nCells,readDepth,totalCost,costKit,samplesPerLane,
                               costFlowCell,readsPerFlowcell))
}

#' Estimate possible number of cells per individual depending on the total cost
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
#' @return Cells per individual that can be measured with this budget and other parameters
#'
#' @export
cellsBudgetCalculation.restrictedDoublets<-function(nSamples,readDepth,totalCost,
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
#' @param nCells Cells per individual
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

  #Estimate individuals per lane
  samplesPerLane<-floor(cellsPerLane/nCells)

  return(readDepthBudgetCalculation(nSamples,nCells,totalCost,costKit,samplesPerLane,
                                     costFlowCell,readsPerFlowcell))
}


#' Calculate total cost dependent on parameters for 10X design
#' using library preparation costs per cell
#'
#' @param nSamples Number of samples
#' @param nCells Cells per individual
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
#' @param nCells Cells per individual
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


#' Estimate possible number of cells per individual depending on the total cost,
#' using library preparation costs per cell
#'
#' @param nSamples Number of samples
#' @param readDepth Read depth per cell
#' @param totalCost Experimental budget
#' @param prepCostCell Library preparation costs per cell
#' @param costFlowCell Cost of one flow cells for sequencing
#' @param readsPerFlowcell Number reads that can be sequenced with one flow cell
#'
#' @return Cells per individual that can be measured with this budget and other parameters
#'
#' @export
cellsBudgetCalculation.libPrepCell<-function(nSamples,readDepth,totalCost,
                                                   prepCostCell, costFlowCell,readsPerFlowcell){

  nCells <- totalCost / (prepCostCell * nSamples +
                           nSamples * readDepth / readsPerFlowcell * costFlowCell)

  return(floor(nCells))

}

#' Estimate possible read depth depending on the total cost,
#' using library preparation costs per cell
#'
#' @param nSamples Number of samples
#' @param nCells Cells per individual
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
