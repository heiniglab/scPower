# Loading libraries
lapply(c("DBI", "dplyr", "plotly", "reshape2", "RPostgres", "RPostgreSQL", "scPower", "shiny",
         "RSQLite", "readr", "magrittr", "zeallot"), 
      library, character.only = TRUE)

establishDBConnection <- function() {
  return(dbConnect(
      Postgres(),
      dbname = "todos",
      host = "localhost",
      port = as.integer("5432"),
      user = "postgres",
      password = "asdasd12x"))
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

offlineQuery <- function(conn) {
  main_table_sql <- read_file(Sys.getenv("MAIN_TABLE_SQLITE"))
  disp_fun_sql <- read_file(Sys.getenv("DISP_FUN_SQLITE"))
  gamma_fits_sql <- read_file(Sys.getenv("GAMMA_LINEAR_FITS_SQLITE"))

  # since sqlite does not have a keyword such as "public", it needs to be removed
  main_table_sql <- gsub("\"public\"\\.", "", main_table_sql)
  disp_fun_sql <- gsub("\"public\"\\.", "", disp_fun_sql)
  gamma_fits_sql <- gsub("\"public\"\\.", "", gamma_fits_sql)

  # Split the SQL script into individual statements
  main_table_sql_statements <- strsplit(main_table_sql, ";", fixed = TRUE)[[1]]
  disp_fun_sql_statements <- strsplit(disp_fun_sql, ";", fixed = TRUE)[[1]]
  gamma_fits_sql_statements <- strsplit(gamma_fits_sql, ";", fixed = TRUE)[[1]]

  for (statement in list(main_table_sql_statements, disp_fun_sql_statements, gamma_fits_sql_statements)) {
    for (sql in statement) {
      try(dbExecute(conn, sql))
    }
  }

  gamma_query <- "
      SELECT g.*, a.id_to_name 
      FROM gamma_linear_fit_results AS g 
      LEFT JOIN main_table AS a
      ON SUBSTRING(g.primary_key, 3) = a.primary_key"

  disp_fun_query <- "
      SELECT d.*, a.id_to_name 
      FROM disp_fun_estimation_results AS d 
      LEFT JOIN main_table AS a 
      ON SUBSTRING(d.primary_key, 3) = a.primary_key"

  gamma_result <- dbGetQuery(conn, gamma_query)
  disp_fun_result <- dbGetQuery(conn, disp_fun_query)

  # gamma changes
  names(gamma_result)[5] <- "ct"
  gamma_result <- gamma_result[, c("parameter", "intercept", "mean_umi", "ct")]
  names(gamma_result)[3] <- "meanUMI"
  gamma_result$intercept <- as.numeric(gamma_result$intercept)
  gamma_result$meanUMI <- as.numeric(gamma_result$meanUMI)

  # dispfun changes
  names(disp_fun_result)[4] <- "ct"
  disp_fun_result <- disp_fun_result[, c("ct", "asympt_disp", "extra_pois")]
  names(disp_fun_result)[2] <- "asymptDisp"
  names(disp_fun_result)[3] <- "extraPois"
  disp_fun_result$asymptDisp <- as.numeric(disp_fun_result$asymptDisp)
  disp_fun_result$extraPois <- as.numeric(disp_fun_result$extraPois)

  dbDisconnect(conn)

  return(list(gamma_result, disp_fun_result))
}

encodeLabel <- function(choice, index, cutoff) {
  if(index <= cutoff){
    return(choice)
  }
  parts <- strsplit(choice, "_")[[1]]
  sprintf("%s (Assay: %s, Tissue: %s)", parts[3], parts[1], parts[2])
}

decodeLabel <- function(encodedChoice) {
  parts <- strsplit(encodedChoice, " \\(Assay: ")[[1]]
  
  if (length(parts) == 1) {
    return(encodedChoice)
  }
  
  cellType <- parts[1]
  subParts <- strsplit(parts[2], ", Tissue: ")[[1]]
  assay <- subParts[1]
  tissue <- substr(subParts[2], 1, nchar(subParts[2]) - 1)
  
  sprintf("%s_%s_%s", assay, tissue, cellType)
}

generateDEScenario <- function(N, ndiff, among_top_N, mean_fc, sd_fc) {
  set.seed(0)
  
  # Initialize the data frame
  de.scenario <- data.frame(name = "Custom", 
                            FoldChange = 2^rnorm(N, sd = sd_fc), 
                            DE = FALSE, 
                            rank = 1:N)
  
  # Randomly assign DE = TRUE to 'ndiff' genes within the top 'among_top_N'
  de.scenario$DE[sample(1:among_top_N, ndiff)] <- TRUE
  
  # Assign fold changes to DE genes based on specified mean and sd
  de.scenario$FoldChange[de.scenario$DE] <- 2^rnorm(ndiff, log2(mean_fc), sd = sd_fc)
  
  # Calculate cumulative fraction of DE genes
  de.scenario$cumFraction <- cumsum(de.scenario$DE) / ndiff

  # only getting DE genes
  de.scenario <- de.scenario[de.scenario$DE == TRUE, ]

  # Filling NA on not to be used fields
  de.scenario$gene <- NA
  de.scenario$FDR <- NA
  de.scenario$geneLength <- NA
  de.scenario$type <- NA
  
  # Reorder columns
  de.scenario <- de.scenario[, c("name", "gene", "FoldChange", "FDR", "cumFraction", "rank", "geneLength", "type")]
  
  return(de.scenario)
}


shinyServer(

  function(input, output, session) {

    onlineToggle <- TRUE
    
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

    # Toggle switch is turned ON/OFF
    observeEvent(input$online, {
      if (input$online) { onlineToggle <<- TRUE } 
      else { onlineToggle <<- FALSE }
    })

    observeEvent(input$refstudy, {
      if(input$refstudy == "Custom") {
        showModal(modalDialog(
          title = "Enter Custom Values for a Reference Study",
          numericInput("n_sim", "Number of Genes", value = 25000),
          numericInput("ndiff", "Number of Relevant Genes", value = 50),
          numericInput("among_top_N", "Ranking Among Top N", value = 5000),
          numericInput("mean_fc", "Fold Change Mean", value = 1.5),
          numericInput("sd_fc", "Fold Change Standard Deviation", value = 0.5),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }
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
    query_organism <- "
        SELECT id_to_name, organism
        FROM main_table
    "
    idToName <- as.list(dbGetQuery(conn, "SELECT id_to_name FROM main_table"))[[1]]
    organism <- dbGetQuery(conn, query_organism)
    organism_list <- setNames(organism$organism, organism$id_to_name)
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

      labels <- sapply(seq_along(choices), function(i) encodeLabel(choices[i], i, 7))
      choices <- setNames(choices, labels)

      updateSelectInput(session, "celltype", label = "Cell type",
                        choices = c("B cells", choices))
      updateSelectInput(session, "tissue", choices = c("Select a tissue..." = "", sort(uniqueTissues)))
      updateSelectInput(session, "assay", choices = c("Select an assay..." = "", sort(uniqueAssays)))
      updateSelectInput(session, "organism", choices = c("Select an organism..." = "", "Homo sapiens", "Mus musculus"))
    })

    observeEvent(input$organism, {
      if (!is.null(input$organism) && input$organism != "") {
        organism_input <- input$organism
        
        # Filter data based on organism selection
        filtered_data <- Filter(function(x) organism_list[x] == organism_input, names(organism_list))
        
        # Extracting assays, tissues, and cell types from filtered data
        assays_tissues_celltypes <- t(sapply(filtered_data, function(x) strsplit(x, "_")[[1]]))
        filtered_assays <- unique(assays_tissues_celltypes[, 1])
        filtered_tissues <- unique(assays_tissues_celltypes[, 2])
        filtered_celltypes <- sapply(seq_along(filtered_data), function(i) encodeLabel(filtered_data[i], i, 0))
        
        # Update UI elements based on filtered options
        updateSelectInput(session, "assay", choices = c("Select an assay..." = "", sort(unique(filtered_assays))))
        updateSelectInput(session, "tissue", choices = c("Select a tissue..." = "", sort(unique(filtered_tissues))))
        updateSelectInput(session, "celltype", choices = c("Select a cell type..." = "", sort(unique(filtered_celltypes))))
      }
    })


    # Update assays and cell types when tissue is selected
    observeEvent(input$tissue, {
      if (!(input$tissue == "")) {
        tissue_input <- input$tissue
        
        # Filter data based on tissue selection, adjusted for organism
        filtered_data <- Filter(function(x) {
          parts <- strsplit(x, "_")[[1]]
          parts[2] == tissue_input
        }, names(organism_list))
        
        # Extract assays and cell types from filtered data, adjust for organism
        filtered_assays <- unique(sapply(filtered_data, function(x) strsplit(x, "_")[[1]][1]))
        filtered_celltypes <- unique(sapply(filtered_data, function(x) strsplit(x, "_")[[1]][3]))
        filtered_organisms <- unique(sapply(filtered_data, function(x) organism_list[x]))
        
        # Update UI elements
        updateSelectInput(session, "assay", choices = c("Select an assay..." = "", sort(unique(filtered_assays))))
        updateSelectInput(session, "celltype", choices = c("Select a cell type..." = "", sort(unique(filtered_celltypes))))
        updateSelectInput(session, "organism", choices = c("Select an organism..." = "", sort(unique(filtered_organisms))))
      }
    })


    # Update tissues and cell types when assay is selected
    observeEvent(input$assay, {
      if (!(input$assay == "")) {
        assay_input <- input$assay
        
        # Filter data based on assay selection, adjusted for organism
        filtered_data <- Filter(function(x) {
          parts <- strsplit(x, "_")[[1]]
          parts[1] == assay_input
        }, names(organism_list))
        
        # Extract tissues and cell types from filtered data, adjust for organism
        filtered_tissues <- unique(sapply(filtered_data, function(x) strsplit(x, "_")[[1]][2]))
        filtered_celltypes <- unique(sapply(filtered_data, function(x) strsplit(x, "_")[[1]][3]))
        filtered_organisms <- unique(sapply(filtered_data, function(x) organism_list[x]))
        
        # Update UI elements
        updateSelectInput(session, "tissue", choices = c("Select a tissue..." = "", sort(unique(filtered_tissues))))
        updateSelectInput(session, "celltype", choices = c("Select a cell type..." = "", sort(unique(filtered_celltypes))))
        updateSelectInput(session, "organism", choices = c("Select an organism..." = "", sort(unique(filtered_organisms))))
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
        studies <- c(studies, "Custom")
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
      ct<-decodeLabel(input$celltype)

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
      
      if(onlineToggle){
        conn <- establishDBConnection()
        gamma.mixed.fits <- rbind(gamma.mixed.fits, constructGammaLinearFit(conn))
        disp.fun.param <- rbind(disp.fun.param, constructDisFunParam(conn))
      }
      else{
        conn <- dbConnect(RSQLite::SQLite(), ":memory:")
        c(gamma_result, disp_fun_result) %<-% offlineQuery(conn)
        
        gamma.mixed.fits <- rbind(gamma.mixed.fits, gamma_result)
        disp.fun.param <- rbind(disp.fun.param, disp_fun_result)
      }

      #Select priors dependent on study type
      if(type=="eqtl"){
        data(eQTLRefStudy)
        ref.study<-eqtl.ref.study

      } else {
        data(DERefStudy)
        ref.study<-de.ref.study
        if(!is.null(input$n_sim)){
          custom.study <- generateDEScenario(input$n_sim, input$ndiff, input$among_top_N, input$mean_fc, input$sd_fc)
          ref.study <- rbind(ref.study, custom.study)
        }
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
          paste("power_results_", Sys.Date(), ".csv", sep="")
        } else {
          paste("power_results_", Sys.Date(), ".tsv", sep="")
        }
      },
      content = function(file) {
        power_results_df <- selectedData()
        
        if (input$fileType == "CSV") {
          write.csv(power_results_df, file, row.names = FALSE)
        } else {
          write.table(power_results_df, file, sep = "\t", row.names = FALSE)
        }
      }
    )

    selectedData <- reactiveVal(NULL)
    
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
        selectedData(max.study)

      } else {
        max.study<-power.study.plot[which.max(power.study.plot$powerDetect),]
        selectedData(max.study)
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
