library(shiny)
library(plotly)
library(ggplot2)

shinyServer(
  function(input, output, session) {
    output$freqPlot<-renderPlotly({

      #Get parameter from GUI
      numSamples<-input$numSamples
      ctFreq<-input$ctFreq
      prob<-input$probCT

      #Number of cells required
      N.min.cells<-seq(input$cellsCT[1],input$cellsCT[2],by=5)

      #Test each possible parameter combination
      parameter.combinations<-expand.grid(prob,N.min.cells,ctFreq,numSamples)
      colnames(parameter.combinations)<-c("prob.cutoff","min.num.cells","cell.type.frac","num.indivs")

      parameter.combinations$sample.size<-mapply(number.cells.detect.celltype,
                                                  parameter.combinations$prob.cutoff,
                                                  parameter.combinations$min.num.cells,
                                                  parameter.combinations$cell.type.frac,
                                                  parameter.combinations$num.indivs)

      #Create plotly plot
      p<-ggplot(parameter.combinations,aes(x=min.num.cells,y=sample.size))+
        geom_line(color="#003E6E")+geom_point(color="#003E6E")+
        xlab("Minimal number of cells from target cell type per individuum")+ylab("Cells per individuum")+
        theme_bw()

      ggplotly(p)

    })

    output$powerPlot<-renderPlotly({

      #Get the parameters
      totalBudget<-input$budget
      readDepthRange<-round(seq(input$rangeReads[1],input$rangeReads[2],
                          length.out = 10))
      cellPersRange<-round(seq(input$rangeCells[1],input$rangeCells[2],
                                length.out = 10))
      type<-input$study
      ref.study.name<-input$ref.study
      ct.freq<-input$ct.freq
      ct<-input$celltype

      costKit<-input$costKit
      costFlowCell<-input$costFlowCell
      readsPerFlowcell<-input$readsPerFlowcell
      personsPerLane<-input$personsLane

      multipletRate<-input$multipletRate
      multipletFactor<-input$multipletFactor
      min.UMI.counts<-input$minUMI
      perc.indiv.expr<-input$percIndiv

      #TODO: Load the current required data fits (need to be chanced for a different version ...)
      # path<-"/Users/katharina.schmid/Documents/singleCellPilotProject/powerScPop/data/"
      path <- system.file("data/", package = "powerScPop")
      read.umi.fit<-readRDS(paste0(path,"readDepthUmiFit.RDS"))
      gamma.mixed.fits<-readRDS(paste0(path,"gamma_umi_fits.RDS"))
      ref.study<-readRDS(paste0(path,"eQTL_ranks_her.RDS"))
      #TODO: this file should be a RDS?
      # disp.fun.param<-read.csv(paste0(path,"dispFunParams_medianEstimate.csv"))
      disp.fun.param<-readRDS(paste0(path,"dispFunParams_medianEstimate.RDS"))

      #Set a name
      ref.study$name<-paste0("Blueprint (",ref.study$cellType,")")

      power.study.plot<-optimize.constant.budget(totalBudget, readDepthRange, cellPersRange,
                                 costKit,costFlowCell,readsPerFlowcell,
                                 ct.freq,type,ref.study,ref.study.name,
                                 personsPerLane,
                                 read.umi.fit,gamma.mixed.fits,ct,
                                 disp.fun.param,
                                 multipletRate,multipletFactor,
                                 min.UMI.counts,perc.indiv.expr)

      power.study.plot$totalCells<-as.factor(power.study.plot$totalCells)
      power.study.plot$readDepth<-as.factor(power.study.plot$readDepth)
      ggplot(power.study.plot,aes(x=readDepth,y=totalCells,fill=powerDetect))+
        geom_tile()+
        ggtitle("Detection power (text shows detection power)")+
        geom_text(aes(label=round(powerDetect,2),color=powerDetect>0.15))+
        scale_colour_manual(values=c("white", "black"))+
        guides(colour=FALSE)

      p<-ggplot(power.study.plot,aes(x=readDepth,y=totalCells,fill=powerDetect))+
        geom_tile()

      ggplotly(p)

    })
  }
)
