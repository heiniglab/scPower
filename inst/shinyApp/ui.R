library(shiny)
library(plotly)

options(width = 250)

shinyUI(
  navbarPage('Single cell design help',
    tabPanel("Power to detect DE/eQTL genes",
     sidebarLayout(
       sidebarPanel(

         h3("Study parameters"),
         actionButton("recalc", "Calculate"),
         hr(),

         h4("Comparsion type"),
         radioButtons("study", label = "Study type",
                      choices = list("DE study" = "de", "eQTL study" = "eqtl"),
                      selected = "eqtl"),
         selectInput("celltype", label = "Target cell type",
                     choices = list()),
         numericInput("ct.freq", label = "Cell type frequency",
                      value = 0.25,step=0.05,min=0,max=1),
         selectInput("ref.study", label = "Reference study",
                     choices = list(),
                     selected = "Blueprint (Monocytes)"),
         numericInput("budget", label = "Total budget", value = 30000, step=1000,min=0),
         sliderInput("rangeReads", "Read depth",
                     min = 0, max = 50000,
                     value = c(2000,50000)),
         sliderInput("rangeCells", "Cells per individual",
                     min = 0, max = 12000,
                     value = c(100,12000)),
         hr(),

         h4("Cost and experimental parameters"),
         numericInput("costKit", label = "Cost for 1 10X kit",
                      value = 5600, step=100,min=0),
         numericInput("costFlowCell", label = "Cost for 1 flow cell",
                      value = 14032, step=100,min=0),
         numericInput("readsPerFlowcell", label = "Number reads per flow cell",
                      value = 4100*10^6, step=10000,min=0),
         numericInput("cellsLane", label = "Cells per lane", value = 20000,
                      step=500,min=0),

         h4("Mapping and Multiplet estimation"),
         numericInput("map.eff", label = "Mapping efficiency",
                      value = 0.8,step=0.05,min=0,max=1),
         numericInput("multipletRate", label = "Multiplet rate",
                      value = 7.67e-06,step=1e-6,min=0),
         numericInput("multipletFactor", label = "Multiplet factor",
                      value = 1.82, step=0.1,min=1),

         h4("Expression cutoffs"),
         numericInput("minUMI", label = "Minimal number of UMI per gene",
                      value = 3, step=1,min=1),
         numericInput("percIndiv", label = "Fraction of individuals",
                      value = 0.5,step=0.05,min=0,max=1),
         width = 3
       ),
       mainPanel(
         h3("Power to detect DE/eQTL genes"),
         plotlyOutput("powerPlot"),
         br(),
         div(paste0("The figure represents the detection power that can be gained with each parameter combination of ",
                    "cells per individual and read depth. ",
                    "The third parameter, the sample size, is defined uniquely by the other two parameters and the ",
                    "overall experimental budget. The calculation for a specific parameter combination can take ",
                    "up to 1-2 minutes.")),
         h4("Click on a specific point in the plot to visualize the exact trace:"),
         br(),
         plotlyOutput("readPlot")
       )
       )
    ),
    tabPanel("Power to detect cell types",
             sidebarLayout(
               sidebarPanel(

                 h3("Study parameters"),
                 numericInput("numSamples", label = "Number of samples", value = 50, step=10,min=0),
                 numericInput("ctFreq", label = "Frequency of cell type of interest",
                              value = 0.1,step=0.05,min=0,max=1),
                 numericInput("probCT", label = "Probability to detect cell type",
                              value = 0.95,step=0.05),
                 sliderInput("cellsCT", label = "Minimal number of cells", min = 1,
                             max = 100, value = c(10, 50)),
                 hr(),
                 width = 3
               ),
               mainPanel(
                 h3("Power to detect rare cell types"),
                 plotlyOutput("freqPlot"),
                 br(),
                 div("The figure shows the required number of cells per individual (y-axis, log scale)
                     to detect the minimal number of cells from a target cell type per individuum (x-axis) with a certain
                     probability. The power depends on the total number of individuals and the frequency of the
                     target cell type. Note that the required number of cells per sample only counts
                     correctly measured cells (no doublets etc), so the number is a lower bound for the required cells to be sequenced.
                     ")
                 )
              )
          )
  )
)
