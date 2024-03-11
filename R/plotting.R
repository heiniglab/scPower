#############################################################
# Functions for plotting to visualize model and results
##########################################################

#' Plotting function to visualize simulated mean values from gamma distribution
#' compare to mean values used for the fitting
#'
#' @param mean.vals Vector with mean values used for fitting the gamma distribution
#' @param gamma.parameters Parameters of the gamma distribution
#' (already filtered for the cell type of interest)
#' @param nGenes number of genes used for fitting the gamma function
#' (parameter num.genes.kept in the function mixed.gamma.estimation)
#' @param lower.dist Lowest value which should be shown in the histogram
#' @param zero.pos Position where the zero component should be plotted at
#' (needs to be smaller than lower.dist)
#' @param num.bins Number of bins in the histogram
#'
#' @return ggplot object of a histogram with logarithmized x axis
#'
#' @export
#'
visualize.gamma.fits<-function(mean.vals,gamma.parameters,nGenes=21000,
                          lower.dist=1e-5,zero.pos=1e-6, num.bins=30){

  require(ggplot2)

  if(zero.pos >= lower.dist){
    stop(paste0("Please enter a value for zero.pos smaller than for lower.dist ",
               "(",lower.dist,")."))
  }

  #Simulate mean values
  sim.mean.vals<-sample.mean.values.quantiles(gamma.parameters, nGenes)

  #Set p1 to a minimal value of 0.01
  gamma.parameters$p1<-max(gamma.parameters$p1,0.01)

  #Calculate from which component each value in the simulation originates (ZI + LZG1 + LZG2)
  nZeros <-round(nGenes*gamma.parameters$p1)
  nGamma1<-round(nGenes*gamma.parameters$p2)
  nGamma2<-round(nGenes*gamma.parameters$p3)
  nGamma3<-nGenes-nZeros-nGamma1-nGamma2

  #Remove zero values from mean vector to reduce it to a size of nGenes
  zeroGenes<-mean.vals==0
  num.zeros.keep<-nGenes-sum(!zeroGenes)

  if(num.zeros.keep<0){
    stop(paste("There are",sum(!zeroGenes),"genes with positive expression.",
               "Increase the nGenes parameter to a value larger than that!"))
  } else if (num.zeros.keep>0){
    zeroGenes[zeroGenes][1:num.zeros.keep]<-FALSE
  }

  mean.vals<-mean.vals[!zeroGenes]

  plot.values <- data.frame(means = c(mean.vals, sim.mean.vals),
                          type = c(rep("Observed", nGenes), rep("sim.mean.zero", nZeros),
                                   rep("sim.mean.c1", nGamma1), rep("sim.mean.c2", nGamma2),
                                   rep("sim.mean.c3", nGamma3)))

  #Set zero values to a very small value to be able to plot on logarithmic scale
  plot.values<-plot.values[plot.values$means==0 | plot.values$means>=lower.dist,]
  plot.values$means[plot.values$means==0]<-zero.pos

  #Remove very large mean values (not shown in the plot)
  plot.values<-plot.values[plot.values$means<1e2,]

  #Order estimated means
  plot.values$type <- factor(plot.values$type,
                           levels = c("Observed", "sim.mean.zero",
                                      "sim.mean.c1", "sim.mean.c2", "sim.mean.c3"))

  g<-ggplot(data=plot.values,aes(x=means, fill=type))+
    geom_histogram(position="dodge",bins=num.bins)+
    scale_x_log10(limits = c(1e-7,1e2),
                  breaks=c(1e-6,1e-4,1e-2,1e0,1e2),
                  labels=c(0,c(1e-4,1e-2,1e0,1e2)))+
    geom_vline(xintercept=lower.dist)+
    theme_bw()+
    xlab("Deseq fitted mean value (logarithmized)")+ylab("Frequency")+
    labs(fill="Mean value")

  return(g)

}

#' Plotting function to visualize calculated power as heatmap, splitted into the
#' separate components of expression probability, power and detection power
#'
#' @param powerDf Data frame with power calculation for different parameter combinations
#' @param var.axis1 Variable which should be plotted at the x axis
#' (possible values: sampleSize, totalCells and readDepth)
#' @param var.axis2 Variable which should be plotted at the y axis
#' (possible values: sampleSize, totalCells and readDepth)
#'
#' @return ggplot object with a heatmap of the two power components plus overall power
#'
#' @export
#'
visualize.power.grid<-function(powerDf,var.axis1="sampleSize",var.axis2="totalCells"){
  require(ggplot2)
  require(reshape2)

  kept.columns<-setdiff(colnames(powerDf),c("powerDetect","exp.probs","power"))
  plot.power<-reshape2::melt(powerDf,id.vars=kept.columns)

  plot.power$variable<-factor(plot.power$variable,
                              levels=c("exp.probs","power","powerDetect"))
  var.labs<-c("Expression probability", "Power", "Detection power")
  names(var.labs)<-c("exp.probs","power","powerDetect")

  plot.power[,var.axis1]<-as.factor(plot.power[,var.axis1])
  plot.power[,var.axis2]<-as.factor(plot.power[,var.axis2])

  g<-ggplot(plot.power,aes_string(x=var.axis1,y=var.axis2,fill="value"))+
      geom_tile()+
      facet_wrap(~variable,ncol=3,labeller = labeller(variable=var.labs))+
      theme_bw()

  return(g)
}
