library(shiny)
library(plotly)
library(ggplot2)
library(gridExtra)
library(reshape)

shinyServer(

  function(input, output, session) {

    ###############################
    # Power to detect rare cell type

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

      plot_ly(parameter.combinations, x = ~min.num.cells, y = ~sample.size,
              type = 'scatter', mode = 'lines+markers',
              hoverinfo = 'text',
              text = ~paste('Minimal number cells per celltype: ',min.num.cells,
                            '<br> Cells per individuum: ', sample.size))%>%
        layout(xaxis = list(title="Minimal number of cells from target cell type per individuum"),
               yaxis = list(title="Cells per individuum"))

    })

    ###############################
    # Power to detect DE/eQTL genes

    #Set the cell types correctly
    observe({
      data(gammaUmiFits)
      celltypes<-as.character(unique(gamma.mixed.fits$ct))
      choices<-setNames(celltypes,celltypes)
      updateSelectInput(session, "celltype", label = "Target cell type",
                        choices = choices)
    })

    #Set the accessible prior studies correctly (dependent on eQTL/DE)
    observe({
      if(input$study=="eqtl"){
        data(eQTLRefStudy)
        studies<-as.character(unique(eqtl.ref.study$name))
      } else {
        data(DERefStudy)
        studies<-as.character(unique(de.ref.study$name))
      }
      choices<-setNames(studies,studies)
      updateSelectInput(session,"ref.study", label = "Reference study",
                        choices = choices)
    })

    #Calculate power grid for the current version of
    powerFrame <- eventReactive(input$recalc, {

      message("Detection power for the current parameter combination is calculated, please wait.")

      #Reset data of the click event
      session$userData$plotlyInputStore[["plotly_click-powerMap"]]<-NULL

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
      cellsPerLane<-input$cellsLane

      mappingEfficiency<-input$map.eff
      multipletRate<-input$multipletRate
      multipletFactor<-input$multipletFactor
      min.UMI.counts<-input$minUMI
      perc.indiv.expr<-input$percIndiv

      #Load required data sets
      data(readDepthUmiFit) #Relation between reads and UMI
      data(gammaUmiFits) #Relation between UMI and gamma parameters
      data(dispFunParam) #Parameters of mean-dispersion curve

      #Select priors dependent on study type
      if(type=="eqtl"){
        data(eQTLRefStudy)
        ref.study<-eqtl.ref.study

      } else {
        data(DERefStudy)
        ref.study<-de.ref.study
      }

      power.study.plot<-optimize.constant.budget.restrictedDoublets(totalBudget, readDepthRange,
                                                                    cellPersRange,
                                                 costKit,costFlowCell,readsPerFlowcell,
                                                 ct.freq,type,ref.study,ref.study.name,
                                                 cellsPerLane,
                                                 read.umi.fit,gamma.mixed.fits,
                                                 gamma.probs,ct,
                                                 disp.fun.param,mappingEfficiency,
                                                 multipletRate,multipletFactor,
                                                 min.UMI.counts,perc.indiv.expr,
                                                 samplingMethod="quantiles")

      message("Calculation finished.")

      return(power.study.plot)
    })

    output$powerPlot<-renderPlotly({

      power.study.plot<-powerFrame()

      power.study.plot$totalCells<-as.factor(power.study.plot$totalCells)
      power.study.plot$readDepth<-as.factor(power.study.plot$readDepth)

      #Round value to not display to many digits
      power.study.plot$powerDetect<-round(power.study.plot$powerDetect,3)

      #Highlight one point
      s<-event_data("plotly_click", source = "powerMap")
      if (! is.null(s)) {
        #Select study of interest dependent on the click
        max.study<-power.study.plot[power.study.plot$totalCells==s[["y"]] & power.study.plot$readDepth==s[["x"]],]

      } else {
        max.study<-power.study.plot[which.max(power.study.plot$powerDetect),]
      }

      colnames(power.study.plot)[2]<-"Detection.power"


      plot_ly(power.study.plot, x=~readDepth,y=~totalCells,z=~Detection.power,type = "heatmap",
              source="powerMap", hoverinfo = 'text',
              text = ~paste('Read depth: ',readDepth,
                            '<br> Cells per individuum: ',totalCells,
                            '<br> Sample size: ', sampleSize,
                            '<br> Detection power: ', Detection.power))%>%
        layout(annotations =  list(showarrow=TRUE, x = max.study$readDepth,
                                   y = max.study$totalCells,text = "Selected <br> study"),
               xaxis = list(title="Read depth"), yaxis = list(title="Cells per individuum"),
               legend=list(title="Detection power"))

    })


    output$readPlot<-renderPlotly({
      power.study<-powerFrame()

      #Get the parameters
      totalBudget<-input$budget
      costKit<-input$costKit
      costFlowCell<-input$costFlowCell
      readsPerFlowcell<-input$readsPerFlowcell
      personsPerLane<-input$personsLane

      s<-event_data("plotly_click", source = "powerMap")
      if (length(s)) {

        #Select study of interest dependent on the click
        max.study<-power.study[power.study$totalCells==s[["y"]] & power.study$readDepth==s[["x"]],]

      } else {
        #Select study with the maximal values
        max.study<-power.study[which.max(power.study$powerDetect),]
      }

      #Get study type
      if(input$study=="eqtl"){
        powerName<-"eQTL power"
      } else {
        powerName<-"DE power"
      }

      ##############
      #Plot cells per person
      power.study.plot<-power.study[power.study$readDepth==max.study$readDepth,]

      #Replace column names
      colnames(power.study.plot)[2:4]<-c("Detection power","Expression probability",powerName)

      power.study.plot<-melt(power.study.plot,id.vars=c("name","sampleSize","readDepth","totalCells",
                                                        "usableCells","multipletFraction",
                                                        "ctCells","readDepthSinglet",
                                                        "mappedReadDepth","expressedGenes"))

      #Round value to not display to many digits
      power.study.plot$value<-round(power.study.plot$value,3)

      p.cp <- plot_ly(power.study.plot, x = ~totalCells, y = ~value,
                      color=~variable, type = 'scatter', mode = 'lines+markers',
                      legendgroup = ~variable,
                      hoverinfo = 'text',
                      text = ~paste('Cells per individuum: ',totalCells,
                                    '<br> Sample size: ', sampleSize,
                                    '<br>',variable,': ',value))%>%
        layout(xaxis = list(title="Cells per individuum"), yaxis = list(title="Probability"),
               legend=list(x=-.1, y=1.2,orientation = 'h'),
               shapes=list(type='line', x0= max.study$totalCells, x1= max.study$totalCells,
                           y0=0, y1=1, line=list(dash='dot', width=1)))

      ##########
      #Plot read depth
      power.study.plot<-power.study[power.study$totalCells==max.study$totalCells,]

      #Replace column names
      colnames(power.study.plot)[2:4]<-c("Detection power","Expression probability",powerName)

      power.study.plot<-melt(power.study.plot,id.vars=c("name","sampleSize","totalCells",
                                                        "usableCells","multipletFraction",
                                                        "ctCells","readDepth","readDepthSinglet",
                                                        "mappedReadDepth","expressedGenes"))

      #Round value to not display to many digits
      power.study.plot$value<-round(power.study.plot$value,3)

      p.rd <- plot_ly(power.study.plot, x = ~readDepth, y = ~value,
                      color=~variable, type = 'scatter', mode = 'lines+markers',
                      legendgroup = ~variable, showlegend=FALSE,
                      hoverinfo = 'text',
                      text = ~paste('Read depth: ',readDepth,
                                    '<br> Sample size: ', sampleSize,
                                    '<br>',variable,': ',value))%>%
        layout(xaxis = list(title="Read depth"), yaxis = list(title="Probability"),
               shapes=list(type='line', x0= max.study$readDepth, x1= max.study$readDepth,
                           y0=0, y1=1, line=list(dash='dot', width=1)))

      subplot(p.cp,p.rd,shareY=TRUE,titleX=TRUE)

    })


  }
)
