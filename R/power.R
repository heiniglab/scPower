######################################
# Functions for power calculation
#####################################

#' Power calculation for cell type identification
#'
#' Calculate probability to detect at least min.num.cells of a cell type
#' with class frequency cell.type.frac in a sample of size nCells
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
#' with class frequency cell.type.frac for each of nSamples individual
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

#' Power calculation for a DE/eQTL study
#'
#' This function to calculate the detection power for a DE or eQTL study, given DE/eQTL genes from a reference study
#' in a single cell RNAseq study. The power depends on the cost determining parameter of sample size, number of cells
#' per person and read depth
#' @param nSamples Sample size
#' @param nCells Number of cells per person
#' @param readDepth Target read depth per cell
#' @param ct.freq Frequency of the cell type of interest
#' @param type (eqtl/de) study
#' @param ref.study Data frame with reference studies to be used for effect sizes and ranks (required columns:)
#' @param ref.study.name Name of the reference study. Will be checked in the ref.study data frame for it (as column name).
#' @param personsPerLane Maximal number of persons per 10X lane
#' @param read.umi.fit Data frame for fitting the mean UMI counts per cell depending on the mean readds per cell (required columns: )
#' @param gamma.mixed.fits Data frame with gamma mixed fit parameters for each cell type (required columns: )
#' @param multipletRate
#' @param multipletFactor
#'
#'
#' @return Power to detect the DE/eQTL genes from the reference study in a single cell experiment with these parameters
#'
#' @export
power.general.withDoublets<-function(nSamples,nCells,readDepth,ct.freq,
                                     type,ref.study,ref.study.name,
                                     personsPerLane,
                                     read.umi.fit,gamma.mixed.fits,
                                     multipletRate,multipletFactor,
                                     min.UMI.counts=10,perc.indiv.expr=0.5){

  #Estimate multiplet rate and "real read depth"
  #Estimate multiplet fraction dependent on cells per lane
  multipletFraction<-multipletRate*nCells*personsPerLane
  usableCells<-round((1-multipletFraction)*nCells)
  readDepthSinglet<-readDepth*nCells/(usableCells+multipletFactor*(nCells-usableCells))

  #Get the fraction of cell type cells
  ctCells<-round(usableCells*ct.freq)

  #Get mean umi dependent on read depth
  umiCounts<-read.umi.fit$intercept+read.umi.fit$reads*log(readDepthSinglet)

  #Get gamma values dependent on mean umi
  gamma.fits.ct<-gamma.mixed.fits[gamma.mixed.fits$ct==ct,]
  gamma.fits.ct$fitted.value<-gamma.fits.ct$intercept+gamma.fits.ct$meanUMI*umiCounts

  #Sample means of 21,000 genes
  gene.means<-rmixgamma(nGenes,
                        pi=c(0.95,0.05),
                        mu=c(gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mu.c1"],
                             gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="mu.c2"]),
                        sd=c(gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd.c1"],
                             gamma.fits.ct$fitted.value[gamma.fits.ct$parameter=="sd.c2"]))

  #Fitted dispersion parameter
  gene.disps<-disp.fun.param$asymptDisp[disp.fun.param$ct==ct]+
    disp.fun.param$extraPois[disp.fun.param$ct==ct]/gene.means

  sim.genes<-data.frame(mean=gene.means, disp=gene.disps)

  #Sort simulated genes by mean expression
  sim.genes<-sim.genes[order(sim.genes$mean, decreasing = TRUE),]

  #Calculate the mean per cell type for each individuum
  sim.genes$mean.sum<-sim.genes$mean*ctCells
  sim.genes$disp.sum<-sim.genes$disp/ctCells

  #Calculate for each gene the expression probability with the definition
  #expressed in 50% of the individuals with count > 10
  sim.genes$exp.probs<-estimate.gene.counts(sim.genes$mean,1/sim.genes$disp,ctCells,
                                            min.counts=min.UMI.counts,num.indivs=nSamples,perc.indiv=perc.indiv.expr)

  #Calculate the expected number of expressed genes
  exp.genes<-round(sum(sim.genes$exp.probs))

  #Set the simulated DE genes as the genes at the same rank position as the original DE genes
  ranks<-ref.study$rank[ref.study$name==ref.study.name]

  #Set all DE with rank > 21000 to 21000 (expression anyway nearly 0)
  ranks[ranks>21000]<-21000

  #Calculate alpha parameter corrected for multiple testing
  alpha<-0.05/exp.genes

  #Choose the DE genes according to the rank
  foundSignGenes<-sim.genes[ranks,]

  #Calculate power
  if(type=="eqtl"){

    #Set the Rsq respectively
    foundSignGenes$Rsq<-ref.study$Rsq[ref.study$name==ref.study.name]

    foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.eqtl(foundSignGenes$Rsq[i],
                                                                               alpha,
                                                                               nSamples))
  } else if (type=="de") {

    #Set the fold change respectively
    foundSignGenes$FoldChange<-ref.study$FoldChange[ref.study$name==ref.study.name]

    foundSignGenes$power<-sapply(1:nrow(foundSignGenes),function(i) power.de(
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

  power.study<-data.table(name=ref.study.name,powerDetect=foundSignGenes$combined.prob,
                          exp.probs=foundSignGenes$exp.probs,power=foundSignGenes$power)
  power.study<-power.study[,.(powerDetect=mean(powerDetect),
                              expProb=mean(exp.probs),
                              power=mean(power)),by=.(name)]
  #Get parameters
  power.study$sampleSize<-nSamples
  power.study$totalCells<-nCells
  power.study$usableCells<-usableCells
  power.study$multipletFraction<-multipletFraction
  power.study$ctCells<-ctCells
  power.study$readDepth<-readDepth
  power.study$readDepthSinglet<-readDepthSinglet
  power.study$expressedGenes<-exp.genes

  return(power.study)
}

#' Optimizing cost parameters to maximize detection power for a given budget
#'
#' This function ...
#'
optimize.constant.budget<-function(totalBudget, readDepthRange, cellPropRange){

}

#' Power calculation for an eQTL gene
#'
#' This function calculates the power to detect an eQTL gene.
#' @param heritability Heritability of the trait
#' @param sig.level Significane threshold
#' @param nSamples Sample size
#'
#' @return Power to detect the eQTL gene
power.eqtl <- function(heritability, sig.level, nSamples) {
  ## determine the rejection area under the null model (standard normal)
  reject <- qnorm(1 - sig.level)
  ## determine the non-centrality paramter
  z <- sqrt((nSamples * heritability) / (1 - heritability))
  ## get the probability to be in the rejection area given that alternative is true
  ## P(reject H0 | H1 true) = P(Z > reject | H1 true)
  power <- pnorm(reject, lower.tail=FALSE, mean=z)
  return(power)
}

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
