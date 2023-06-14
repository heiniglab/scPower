# Loading libraries
lapply(c("DBI", "dplyr", "plotly", "reshape2", "RPostgres", "RPostgreSQL", "scPower", "shiny"), 
      library, character.only = TRUE)

establishDBConnection <- function() {
  return(dbConnect(
      Postgres(),
      dbname = Sys.getenv("POSTGRES_DBNAME"),
      host = Sys.getenv("POSTGRES_HOST"),
      port = as.integer(Sys.getenv("POSTGRES_PORT")),
      user = Sys.getenv("POSTGRES_USER"),
      password = Sys.getenv("POSTGRES_PASSWORD")))
}

constructGammaLinearFit <- function(conn) {
  query <- "
    SELECT g.*, a.id_to_name 
    FROM gamma_linear_fit_results AS g 
    LEFT JOIN main_table AS a 
    ON encode(g.primary_key, 'hex') = a.primary_key
  "

  # Execute the SQL query and fetch results
  result <- dbGetQuery(conn, query)
  names(result)[5] <- "ct"
  result <- result[, c("parameter", "intercept", "mean_umi", "ct")]
  names(result)[3] <- "meanUMI"
  result$intercept <- as.numeric(result$intercept)
  result$meanUMI <- as.numeric(result$meanUMI)

  return(result)
}

constructDisFunParam <- function(conn) {
  query <- "
    SELECT g.*, a.id_to_name 
    FROM disp_fun_estimation_results AS g 
    LEFT JOIN main_table AS a 
    ON encode(g.primary_key, 'hex') = a.primary_key
  "

  # Execute the SQL query and fetch results
  result <- dbGetQuery(conn, query)
  names(result)[4] <- "ct"
  result <- result[, c("ct", "asympt_disp", "extra_pois")]
  names(result)[2] <- "asymptDisp"
  names(result)[3] <- "extraPois"
  result$asymptDisp <- as.numeric(result$asymptDisp)
  result$extraPois <- as.numeric(result$extraPois)
  
  dbDisconnect(conn)
  return(result)
}

formatLabel <- function(choice, index, cutoff) {
  if(index <= cutoff){
    return(choice)
  }
  parts <- strsplit(choice, "_")[[1]]
  sprintf("%s (Assay: %s, Tissue: %s)", parts[3], parts[1], parts[2])
}

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

    ## whether the advanced input options should be displayed or not
    observeEvent(input$advanced, {
      
      if(input$advanced == "no"){
        shinyjs::hide(id = "cost")
        shinyjs::hide(id = "multiplet")
      }else{
        shinyjs::show(id = "cost")
        shinyjs::show(id = "multiplet")
      }
    })
    
    #Set the cell types correctly
    conn <- establishDBConnection()
    query <- "SELECT id_to_name FROM main_table"
    idToName <- as.list(dbGetQuery(conn, query))[[1]]
    dbDisconnect(conn)
    uniqueAssays <- unique(sapply(idToName, function(x) unlist(strsplit(x, "_"))[1]))
    uniqueTissues <- unique(sapply(idToName, function(x) unlist(strsplit(x, "_"))[2]))
    uniqueCelltypes <- unique(sapply(idToName, function(x) unlist(strsplit(x, "_"))[3]))

    # initializing each of the dropdown menu items
    observe({
      data(gammaUmiFits)
      celltypes<-as.character(unique(gamma.mixed.fits$ct))
      celltypes<-append(celltypes, idToName)
      choices <- celltypes

      labels <- sapply(seq_along(choices), function(i) formatLabel(choices[i], i, 7))
      choices <- setNames(choices, labels)

      updateSelectInput(session, "celltype", label = "Cell type",
                        choices = c("B cells", choices))
      updateSelectInput(session, "tissue", choices = c("", sort(uniqueTissues)))
      updateSelectInput(session, "assay", choices = c("", sort(uniqueAssays)))
    })

    # Update cell types when tissue is selected
    observeEvent(input$tissue, {
      if (!(input$tissue == "")) {
        tissue_input <- input$tissue
        filtered_celltypes <- idToName[sapply(idToName, function(x) {
          parts <- strsplit(x, "_")[[1]]
          tissue_part <- parts[2]
          tissue_part == tissue_input
        })]
        filtered_assays <- unique(sapply(filtered_celltypes, function(x) unlist(strsplit(x, "_"))[1]))
        choices <- sapply(seq_along(filtered_celltypes), function(i) formatLabel(filtered_celltypes[i], i, 0))
        
        updateSelectInput(session, "assay", choices = filtered_assays)
        updateSelectInput(session, "celltype", choices = sort(choices))
      }
    })

    # Update cell types when assay is selected
    observeEvent(input$assay, {
      if (!(input$assay == "")) {
        tissue_input <- input$tissue
        assay_input <- input$assay
        filtered_celltypes <- idToName[sapply(idToName, function(x) {
          parts <- strsplit(x, "_")[[1]]
          assay_part <- parts[1]
          tissue_part <- parts[2]
          assay_part == assay_input & tissue_part == tissue_input
        })]
        choices <- sapply(seq_along(filtered_celltypes), function(i) formatLabel(filtered_celltypes[i], i, 0))
        updateSelectInput(session, "celltype", choices = sort(choices))
      }
    })

    #Set the accessible prior studies correctly (dependent on eQTL/DE)
    observe({
      if(input$study=="eqtl"){
        data(eQTLRefStudy)
        studies<-as.character(unique(eqtl.ref.study$name))
        selected<-"Bonferroni"
        shinyjs::hide(id="ssizeratiode")
        shinyjs::show(id="indepsnps")
      } else {
        data(DERefStudy)
        studies<-as.character(unique(de.ref.study$name))
        selected<-"FDR"
        shinyjs::hide(id="indepsnps")
        shinyjs::show(id="ssizeratiode")
      }
      choices<-setNames(studies,studies)
      updateSelectInput(session,"refstudy", label = "Reference study",
                        choices = choices)

      #Choose the preferred MT method
      updateSelectInput(session,"MTmethod",
                        selected=selected)
    })

    #Set labels for the parameter pair properly
    observe({if(input$grid=="sc"){
              updateNumericInput(session,"rangeX_min", label="Total sample size (min)",
                           value = 10)
              updateNumericInput(session,"rangeX_max", label="Total sample size (max)",
                         value = 50)
              updateNumericInput(session,"rangeY_min", label="Cells (min)",
                           value = 2000)
              updateNumericInput(session,"rangeY_max", label="Cells (max)",
                                 value = 10000)
            } else if(input$grid=="sr"){
              updateNumericInput(session,"rangeX_min", label="Total sample size (min)",
                                 value = 10)
              updateNumericInput(session,"rangeX_max", label="Total sample size (max)",
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
      ref.study.name<-input$refstudy
      ct.freq<-input$ctfreq
      ct<-input$celltype

      costKit<-input$costKit
      costFlowCell<-input$costFlowCell
      readsPerFlowcell<-input$readsPerFlowcell
      cellsPerLane<-input$cellsLane

      sign.threshold<-input$pval
      MTmethod<-input$MTmethod

      mappingEfficiency<-input$map.eff
      multipletRate<-input$multipletRate
      multipletFactor<-input$multipletFactor
      min.UMI.counts<-input$minUMI
      perc.indiv.expr<-input$percIndiv

      useSimulatedPower<-input$simPower
      speedPowerCalc<-input$speedCalc
      
      indepSNPs<-input$indepsnps
      ssize.ratio.de<-input$ssizeratiode
      
      reactionsPerKit <- input$reactionsPerKit

      #Load required data sets
      data(readDepthUmiFit) #Relation between reads and UMI
      data(gammaUmiFits) #Relation between UMI and gamma parameters
      data(dispFunParam) #Parameters of mean-dispersion curve

      conn <- establishDBConnection()
      gamma.mixed.fits <- rbind(gamma.mixed.fits, constructGammaLinearFit(conn))
      disp.fun.param <- rbind(disp.fun.param, constructDisFunParam(conn))

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
                                                                    read.umi.fit[read.umi.fit$type=="10X_PBMC_1",],
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
                                                                    samplingMethod="quantiles",
                                                                    sign.threshold =sign.threshold,
                                                                    MTmethod = MTmethod,
                                                                    useSimulatedPower = useSimulatedPower,
                                                                    speedPowerCalc = speedPowerCalc,
                                                                    indepSNPs=indepSNPs,
                                                                    ssize.ratio.de=ssize.ratio.de,
                                                                    reactionsPerKit = reactionsPerKit),
        message="Calculating power optimization!", value=0.5
      )

      message("Calculation finished.")

      #Save vector with important parameters for plotting
      parameter.vector<-c(selectedPair,totalBudget,costKit,
                          costFlowCell,readsPerFlowcell,type)

      return(list(power.study.plot, parameter.vector))
    }, ignoreNULL = FALSE)

    output$downloadData <- downloadHandler(
      filename = function() {
        if (input$fileType == "CSV") {
          "input_parameters.csv"
        } else {
          "input_parameters.tsv"
        }
      },
      content = function(file) {
        # Convert the input parameters to a list
        input_list <- reactiveValuesToList(input)
        
        # Convert each input to a single string
        input_list <- lapply(input_list, function(x) {
          if (is.vector(x)) {
            paste(x, collapse = ", ")
          } else {
            as.character(x)
          }
        })
        
        # Convert the list to a data frame
        df <- data.frame(t(unlist(input_list)), stringsAsFactors = FALSE)
        
        if (input$fileType == "CSV") {
          write.csv(df, file, row.names = FALSE)
        } else {
          write.table(df, file, sep = "\t", row.names = FALSE)
        }
      }
    )

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
        sizeAxis<-"readDepth"
        sizeAxisLabel<-"Read depth"
      } else if(selectedPair=="sr"){
        xAxis<-"sampleSize"
        xAxisLabel<-"Sample size"
        yAxis<-"readDepth"
        yAxisLabel<-"Read depth"
        sizeAxis<-"totalCells"
        sizeAxisLabel<-"Cells per sample"
      } else {
        xAxis<-"totalCells"
        xAxisLabel<-"Cells per sample"
        yAxis<-"readDepth"
        yAxisLabel<-"Read depth"
        sizeAxis<-"sampleSize"
        sizeAxisLabel<-"Sample size"
      }

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

      plot_ly(data=power.study.plot,
              type = "scatter",
              mode="markers",
              x=as.formula(paste0("~",xAxis)),
              y=as.formula(paste0("~",yAxis)),
              color=~Detection.power,
              size=as.formula(paste0("~",sizeAxis)),
              sizes=c(100,500), #choose a size range for the circles,
              source="powerMap",
              hoverinfo = 'text',
              text = ~paste('Sample size: ', sampleSize,
                            '<br> Cells per individuum: ',totalCells,
                            '<br> Read depth: ',readDepth,
                            '<br> Detection power: ', Detection.power)) %>%
        layout(annotations =  list(showarrow=TRUE,
                                   x = max.study[,c(xAxis)],
                                   y = max.study[,c(yAxis)],
                                   text = "Selected <br> study",
                                   bgcolor  ="white"),
               xaxis = list(title=xAxisLabel),
               yaxis = list(title=yAxisLabel),
               legend=list(title="Detection power")
               )
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
