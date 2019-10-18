#' Estimate gene expression probability
#'
#' This function estimates the expression probability of each gene in pseudobulk
#' with a certain cutoff.
#' 
#' @param mu Estimated mean value in sc data
#' 
#' @export
#' 
estimate.gene.counts<-function(mu,size,cells.per.person,min.counts=10,
                               num.indivs=14,perc.indiv=0.5){
  
  #Generate basis data.frame
  fits.allIndivs<-data.frame(mu=mu,size=size,cells.per.person=cells.per.person)
  
  #Calculate summed NegBinomial distribution
  fits.allIndivs$size<-fits.allIndivs$size*fits.allIndivs$cells.per.person
  fits.allIndivs$mu<-fits.allIndivs$mu*fits.allIndivs$cells.per.person
  
  #Calculate probability to observe > n counts in one individuum
  fits.allIndivs$count.probs<-1 - pnbinom(min.counts,mu=fits.allIndivs$mu,size=fits.allIndivs$size)
  
  #Calculate probability that > prob indivs have count > n
  fits.allIndivs$indiv.probs<-1 - pbinom(perc.indiv*num.indivs,num.indivs,fits.allIndivs$count.probs)
  
  #Return number of expressed cells per celltype
  return(fits.allIndivs$indiv.probs)
}