library(shiny)
library(plotly)
library(reshape2)
library(scPower)

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
      updateSelectInput(session, "celltype", label = "Cell type",
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

    #Set labels for the parameter pair properly
    observe({if(input$grid=="sc"){
              updateNumericInput(session,"rangeX_min", label="Samples (min)",
                           value = 10)
              updateNumericInput(session,"rangeX_max", label="Samples (max)",
                         value = 50)
              updateNumericInput(session,"rangeY_min", label="Cells (min)",
                           value = 2000)
              updateNumericInput(session,"rangeY_max", label="Cells (max)",
                                 value = 10000)
            } else if(input$grid=="sr"){
              updateNumericInput(session,"rangeX_min", label="Samples (min)",
                                 value = 10)
              updateNumericInput(session,"rangeX_max", label="Samples (max)",
                                 value = 50)
              updateNumericInput(session,"rangeY_min", label="Reads (min)",
                                value = 100000)
              updateNumericInput(session,"rangeY_max", label="Reads (max)",
                                 value = 500000)
            } else {
              updateNumericInput(session,"rangeX_min", label="Cells (min)",
                                 value = 2000)
              updateNumericInput(session,"rangeX_max", label="Cells (max)",
                                 value = 10000)
              updateNumericInput(session,"rangeY_min", label="Reads (min)",
                                 value = 100000)
              updateNumericInput(session,"rangeY_max", label="Reads (max)",
                                 value = 500000)
            }
    })

    #Calculate power grid for the current version of
    powerFrame <- eventReactive(input$recalc, {
      message("Detection power for the current parameter combination is calculated, please wait.")

      #Reset data of the click event
      session$userData$plotlyInputStore[["plotly_click-powerMap"]]<-NULL

      #Get the parameters
      totalBudget<-input$budget
      rangeX<-round(seq(input$rangeX_min,input$rangeX_max,
                                length.out = input$steps))
      rangeY<-round(seq(input$rangeY_min,input$rangeY_max,
                               length.out = input$steps))

      #Set parameters dependent on pair combination
      selectedPair<-input$grid
      if(selectedPair=="sc"){
        sampleRange<-rangeX
        cellsRange<-rangeY
        readRange<-NULL
      } else if (selectedPair=="sr"){
        sampleRange<-rangeX
        cellsRange<-NULL
        readRange<-rangeY
      } else {
        sampleRange<-NULL
        cellsRange<-rangeX
        readRange<-rangeY
      }

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

      withProgress(
        expr=power.study.plot<-optimize.constant.budget.restrictedDoublets(totalBudget,type,
                                                                    ct, ct.freq,
                                                                    costKit,costFlowCell,readsPerFlowcell,
                                                                    ref.study,ref.study.name,
                                                                    cellsPerLane,
                                                                    read.umi.fit,
                                                                    gamma.mixed.fits,
                                                                    disp.fun.param,
                                                                    nSamplesRange=sampleRange,
                                                                    nCellsRange=cellsRange,
                                                                    readDepthRange=readRange,
                                                                    mappingEfficiency=mappingEfficiency,
                                                                    multipletRate=multipletRate,
                                                                    multipletFactor=multipletFactor,
                                                                    min.UMI.counts=min.UMI.counts,
                                                                    perc.indiv.expr=perc.indiv.expr,
                                                                    samplingMethod="quantiles"),
        message="Calculating power optimization!", value=0.5
      )

      message("Calculation finished.")

      #Save vector with important parameters for plotting
      parameter.vector<-c(selectedPair,totalBudget,costKit,
                          costFlowCell,readsPerFlowcell,type)

      return(list(power.study.plot, parameter.vector))
    }, ignoreNULL = FALSE)

    #Create main plot with optimization grid
    output$powerPlot<-renderPlotly({

      power.study.plot<-powerFrame()[[1]]
      parameter.vector<-powerFrame()[[2]]
      selectedPair<-parameter.vector[1]

      #Set grid dependent on parameter choice
      if(selectedPair=="sc"){
        xAxis<-"sampleSize"
        xAxisLabel<-"Sample size"
        yAxis<-"totalCells"
        yAxisLabel<-"Cells per sample"
      } else if(selectedPair=="sr"){
        xAxis<-"sampleSize"
        xAxisLabel<-"Sample size"
        yAxis<-"readDepth"
        yAxisLabel<-"Read depth"
      } else {
        xAxis<-"totalCells"
        xAxisLabel<-"Cells per sample"
        yAxis<-"readDepth"
        yAxisLabel<-"Read depth"
      }

      power.study.plot[,c(xAxis)]<-as.factor(power.study.plot[,c(xAxis)])
      power.study.plot[,c(yAxis)]<-as.factor(power.study.plot[,c(yAxis)])

      #Round value to not display to many digits
      power.study.plot$powerDetect<-round(power.study.plot$powerDetect,3)

      #Highlight one point
      s<-event_data("plotly_click", source = "powerMap")
      if (! is.null(s)) {
        #Select study of interest dependent on the click
        max.study<-power.study.plot[power.study.plot[,c(xAxis)]==s[["x"]] &
                                      power.study.plot[,c(yAxis)]==s[["y"]],]

      } else {
        max.study<-power.study.plot[which.max(power.study.plot$powerDetect),]
      }

      colnames(power.study.plot)[2]<-"Detection.power"

      plot_ly(power.study.plot, x=as.formula(paste0("~",xAxis)),y=as.formula(paste0("~",yAxis)),
              z=~Detection.power,type = "heatmap",
              source="powerMap", hoverinfo = 'text',
              text = ~paste('Sample size: ', sampleSize,
                            '<br> Cells per individuum: ',totalCells,
                            '<br> Read depth: ',readDepth,
                            '<br> Detection power: ', Detection.power))%>%
        layout(annotations =  list(showarrow=TRUE, x = max.study[,c(xAxis)],
                                   y = max.study[,c(yAxis)],text = "Selected <br> study",
                                   bgcolor  ="white"),
               xaxis = list(title=xAxisLabel), yaxis = list(title=yAxisLabel),
               legend=list(title="Detection power"))


    })

    output$readPlot<-renderPlotly({

      power.study<-powerFrame()[[1]]
      parameter.vector<-powerFrame()[[2]]

      selectedPair<-parameter.vector[1]
      totalBudget<-parameter.vector[2]
      costKit<-parameter.vector[3]
      costFlowCell<-parameter.vector[4]
      readsPerFlowcell<-parameter.vector[5]
      studyType<-parameter.vector[6]

      #Set grid dependent on parameter choice
      if(selectedPair=="sc"){
        xAxis<-"sampleSize"
        xAxisLabel<-"Sample size"
        yAxis<-"totalCells"
        yAxisLabel<-"Cells per sample"
      } else if(selectedPair=="sr"){
        xAxis<-"sampleSize"
        xAxisLabel<-"Sample size"
        yAxis<-"readDepth"
        yAxisLabel<-"Read depth"
      } else {
        xAxis<-"totalCells"
        xAxisLabel<-"Cells per sample"
        yAxis<-"readDepth"
        yAxisLabel<-"Read depth"
      }

      s<-event_data("plotly_click", source = "powerMap")
      if (length(s)) {

        #Select study of interest dependent on the click
        max.study<-power.study[power.study[,xAxis]==s[["x"]] &
                                 power.study[,yAxis]==s[["y"]],]

      } else {
        #Select study with the maximal values
        max.study<-power.study[which.max(power.study$powerDetect),]
      }

      #Get study type
      if(studyType=="eqtl"){
        powerName<-"eQTL power"
      } else {
        powerName<-"DE power"
      }

      ##############
      #Plot cells per person
      power.study.plot<-power.study[power.study[,yAxis]==max.study[,yAxis],]

      #Replace column names
      colnames(power.study.plot)[2:4]<-c("Detection power","Expression probability",powerName)

      power.study.plot<-reshape2::melt(power.study.plot,id.vars=c("name","sampleSize","readDepth","totalCells",
                                                        "usableCells","multipletFraction",
                                                        "ctCells","readDepthSinglet",
                                                        "mappedReadDepth","expressedGenes"))

      #Round value to not display to many digits
      power.study.plot$value<-round(power.study.plot$value,3)

      p.cp <- plot_ly(power.study.plot, x = as.formula(paste0("~",xAxis)), y = ~value,
                      color=~variable, type = 'scatter', mode = 'lines+markers',
                      legendgroup = ~variable,
                      hoverinfo = 'text',
                      text = ~paste('Sample size: ', sampleSize,
                                    '<br> Cells per individuum: ',totalCells,
                                    '<br> Read depth: ',readDepth,
                                    '<br>',variable,': ',value))%>%
        layout(xaxis = list(title=xAxisLabel), yaxis = list(title="Probability"),
               legend=list(x=-.1, y=1.2,orientation = 'h'),
               shapes=list(type='line', x0= max.study[,xAxis], x1= max.study[,xAxis],
                           y0=0, y1=1, line=list(dash='dot', width=1)))

      ##########
      #Plot read depth
      power.study.plot<-power.study[power.study[,xAxis]==max.study[,xAxis],]

      #Replace column names
      colnames(power.study.plot)[2:4]<-c("Detection power","Expression probability",powerName)

      power.study.plot<-reshape2::melt(power.study.plot,id.vars=c("name","sampleSize","totalCells",
                                                        "usableCells","multipletFraction",
                                                        "ctCells","readDepth","readDepthSinglet",
                                                        "mappedReadDepth","expressedGenes"))

      #Round value to not display to many digits
      power.study.plot$value<-round(power.study.plot$value,3)

      p.rd <- plot_ly(power.study.plot, x = as.formula(paste0("~",yAxis)), y = ~value,
                      color=~variable, type = 'scatter', mode = 'lines+markers',
                      legendgroup = ~variable, showlegend=FALSE,
                      hoverinfo = 'text',
                      text = ~paste('Sample size: ', sampleSize,
                                    '<br> Cells per individuum: ',totalCells,
                                    '<br> Read depth: ',readDepth,
                                    '<br>',variable,': ',value))%>%
        layout(xaxis = list(title=yAxisLabel), yaxis = list(title="Probability"),
               shapes=list(type='line', x0= max.study[,yAxis], x1= max.study[,yAxis],
                           y0=0, y1=1, line=list(dash='dot', width=1)))

      subplot(p.cp,p.rd,shareY=TRUE,titleX=TRUE)

    })

  }
)
