
library(shiny)
library(plotly)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)

# the popovers for `selectInput` options need to be made updateresistant, or they won't show in the app
# https://stackoverflow.com/questions/36965954/shinybs-bspopover-and-updateselectinput
updateResistantPopover <- function(id, title, content, placement = "bottom", trigger = "hover", options = NULL){
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options, content)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      var target = document.querySelector('#", id, "');
      var observer = new MutationObserver(function(mutations) {
        setTimeout(function() {
          shinyBS.addTooltip('", id, "', 'popover', ", options, ");
        }, 200);
      });
      observer.observe(target, { childList: true });
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

header <- dashboardHeader(title = "scPower",
                          tags$li(class = "dropdown",
                                  a(href = 'https://doi.org/10.1038/s41467-021-26779-7',
                                    img(src = "NatureComms_logo.png", height = "30px"),
                                    style = "padding-top:10px; padding-bottom:10px;",
                                    title = "Our preprint")
                          ),
                          tags$li(class = "dropdown",
                                  a(href = 'https://github.com/heiniglab/scPower',
                                    img(src = "GitHub_logo.png", height = "30px"),
                                    style = "padding-top:10px; padding-bottom:10px;",
                                    title = "Our github page")
                                  ),

                          # Set height of dashboardHeader
                          tags$li(class = "dropdown",
                                  tags$style(".main-header {height: 50px}"),
                                  tags$style(".main-header .logo {height: 50px}")
                                  )
                          )

## Sidebar content
sidebar <-   dashboardSidebar(
  sidebarMenu(
    menuItem("Welcome", tabName = "welcome", icon = icon("home")),
    menuItem("Description", tabName = "description", icon = icon("book")),
    menuItem("Detect DE/eQTL genes", tabName = "genes", icon = icon("users")),
    menuItem("Detect cell types", tabName = "celltypes", icon = icon("search")),
    menuItem("License Statement", tabName = "license", icon = icon("file"))
  )
)

body <- ## Body content
  dashboardBody(
    useShinyjs(),
    # make buttons on welcome screen clickable (and simultaneously update active menu panel)
    tags$script(HTML("
        var openTab = function(tabName){
          $('a', $('.sidebar')).each(function() {
            if(this.getAttribute('data-value') == tabName) {
              this.click()
            };
          });
        }
      ")),

    # include our custom CSS style
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom_style.css")
    ),
    tabItems(

      tabItem(tabName="welcome",
              h2("Welcome to scPower"),
              h3("What would you like to do?"),
              tags$style(HTML(".small-box {height: 150px}")),
              fluidRow(
                infoBox("",
                        a("Learn more", onclick = "openTab('description')", href="#"),
                        icon=icon("book")),
                infoBox("",
                        a("Detect DE/eQTL genes", onclick = "openTab('genes')", href="#"),
                        icon=icon("users")),
                infoBox("",
                        a("Detect cell types", onclick = "openTab('celltypes')", href="#"),
                        icon=icon("search")),
              )
      ),

      tabItem(tabName = "description",
              h2("Welcome to scPower"),
              h3("- a statistical framework for design and power analysis of multi-sample single cell transcriptomics experiments-"),
              div(paste("The tool supports the user to set the experimental parameters of cell type",
                        "specific inter-individual DE and eQTL analysis using single cell RNA-seq data.")),
              br(),
              img(src='Figure1.png', align = "center", height="95%", width="95%"),
              br(),
              includeHTML("introduction.html")
              ),

      tabItem(tabName = "genes",
              fluidRow(
                column(width=4,

                  box(width = 0,
                    title = "General parameters",
                    solidHeader = TRUE,
                    status="orange",

                    div(style = "display: flex; align-items: center; justify-content: space-between;",
                      div(style = "width: 66.66%; display: flex; justify-content: space-between; margin-right: 2%;",
                        actionButton("recalc", "Calculate optimal study", icon("paper-plane"),
                                    style = "color: #fff; background-color: #337ab7; border-color: #2e6da4; flex-grow: 1;"),
                        bsPopover("recalc", title = "Calculate optimal study", placement = "top", options = list(container = "body"),
                                  content = "Computes the optimal study design for the given parameter combinations. Can take 1-2 minutes for big grids."),
                      ),
                      div(style = "width: 33.33%;",  
                        downloadButton("downloadData", "", style = "width: 100%;"),
                        bsPopover("downloadData", title = "Download", placement = "top", options = list(container = "body"),
                                  content = "Download the data of the current plot as a CSV or TSV file."),
                        # div(style = "padding-top: 12px; margin-left: 10px;",  
                        #     prettySwitch(
                        #       inputId = "online",
                        #       label = "",
                        #       fill = TRUE,
                        #       value = TRUE,
                        #       status = "primary"
                        #     )
                        # )
                      )
                    ),

                    
                    br(),
                    br(),

                    radioButtons("fileType", "File type to download:", inline = TRUE, choices = c("CSV", "TSV")),
                    bsPopover("fileType", title="File type to download", placement="top", options = list(container = "body"),
                              content="The file type of the downloaded file. CSV is the default, but TSV is recommended for Excel users."),
                    radioButtons("study", label = "Study type", inline=TRUE,
                                 choices = list("DE study" = "de", "eQTL study" = "eqtl"),
                                 selected = "de"),
                    bsPopover("study", title="Study type:", placement="top", options = list(container = "body"),
                              content="For what type of study do you want to design an experiment"),

                    selectInput("organism", label = "Organism",
                                choices = list()),
                    updateResistantPopover("organism", title="Organism", placement="top",
                              content=paste("The organism of interest. The datasets are available for different cell types and different studies.")),
                    selectInput("assay", label = "Assay",
                                choices = list()),
                    updateResistantPopover("assay", title="Assay", placement="top",
                              content=paste("The assay of interest. The datasets are available for different cell types and different studies.")),
                    selectInput("tissue", label = "Tissue",
                                choices = list()),
                    updateResistantPopover("tissue", title="Tissue", placement="top",
                              content=paste("The tissue of interest. The datasets are available for different cell types and different studies.")),
                    selectInput("celltype", label = "Cell type",
                                choices = list()),
                    updateResistantPopover("celltype", title="Cell-type specific analysis", placement="top",
                              content=paste("The expression distribution is selected for this cell type. Only one cell type at once can be analysed.",
                                           "If multiple cell types are of interest (which is often the case), try different cell types,",
                                           "focusing especially on the ones with small cell type frequencies.")),

                    numericInput("ctfreq", label = "Cell type frequency",
                                 value = 0.25,step=0.05,min=0,max=1),
                    bsPopover("ctfreq", title="Cell type frequency", placement="top",
                              content="Frequency of the cell type of interest"),

                    numericInput("ssizeratiode", label="Sample size ratio", value=1, min=0, step=0.05),
                    bsPopover("ssizeratiode", title="Sample size ratio", placement="top",
                                           content="ratio between sample size of group 0 (control group) and group 1 (Ratio=1 in case of balanced design)"),

                    selectInput("refstudy", label = "Reference study",
                                choices = list(),
                                selected = "Blueprint (Monocytes)"),

                    updateResistantPopover("refstudy", title="Reference study", placement="top",
                              content=paste("effect sizes and expression ranks are taken from a reference study, performed on FACS sorted bulk RNA-seq data.",
                              "Different examples are available for DE as well as eQTL studies.")),

                    numericInput("budget", label = "Total budget", value = 50000, step=500,min=0),
                    bsPopover("budget", title="Total budget", placement="top", options = list(container = "body"),
                              content="The total budget available for the sequencing"),

                    selectInput("grid", label = "Parameter grid",
                                choices = list("samples - cells per sample"="sc",
                                               "samples - reads per cell"="sr",
                                               "cells per sample - reads per cell"="cr")),
                    bsPopover("grid", title="Parameter grid", placement="top", options = list(container = "body"),
                              content=paste("all possible combinations for two of the three experimental parameters (sample size, cells per person and read depth) are tested,",
                              "the third parameter is defined uniquely given the other two and the overall budget and will be displayed as circle size. ",
                              "Which of the two shall be tested, can be selcted here. Depending on the selection, the four parameters below are adapted.")),

                    numericInput("rangeX_min",label="Total sample size (min)",value=10),
                    bsPopover("rangeX_min", title=" ", placement="top", options = list(container = "body"),
                              content="Minimal value of the tested ranges for the parameter on the x-Axis."),
                    numericInput("rangeX_max",label="Total sample size (max)",value=50),
                    bsPopover("rangeX_max", title=" ", placement="top", options = list(container = "body"),
                              content="Maximum value of the tested ranges for the parameter on the x-Axis."),

                    numericInput("rangeY_min",label="Cells (min)",value=2000),
                    bsPopover("rangeY_min", title=" ", placement="top", options = list(container = "body"),
                              content="Minimal value of the tested ranges for the parameter on the y-Axis"),
                    numericInput("rangeY_max",label="Cells (max)",value=10000),
                    bsPopover("rangeY_max", title=" ", placement="top", options = list(container = "body"),
                              content="Maximum value of the tested ranges for the parameter on the y-Axis"),

                    numericInput("steps","Steps",value=5,min=0,step=1),
                    bsPopover("steps", title="Steps", placement="top", options = list(container = "body"),
                              content="number of values in the parameter ranges for the parameter grid"),
                    hr(),
                    radioButtons("advanced", label = "Show advanced options", inline=TRUE,
                                 choices = list("yes"="yes","no"="no"), selected="no")
                  ),
                  box(width = 0,
                    title="Cost and experimental parameters",
                    id = "cost",
                    solidHeader = TRUE,
                    status="orange",
                    numericInput("costKit", label = "Cost 10X kit",
                                 value = 5600, step=100,min=0),
                    bsPopover("costKit", title="Cost 10X kit", placement="top", options = list(container = "body"),
                              content="Cost for one 10X Genomics kit"),

                    numericInput("costFlowCell", label = "Cost flow cell",
                                 value = 14032, step=100,min=0),
                    bsPopover("costFlowCell", title="Cost flow cell", placement="top", options = list(container = "body"),
                              content=" Cost for one flow cell"),

                    numericInput("readsPerFlowcell", label = "Number of reads per flow cell",
                                 value = 4100*10^6, step=10000,min=0),

                    numericInput("cellsLane", label = "Cells per lane", value = 8000,
                                 step=500,min=0),
                    bsPopover("cellsLane", title="Cells per lane", placement="top", options = list(container = "body"),
                              content="Number of cells meassured on one 10X lane (dependent on the parameter \"Reactions Per Kit\")"),

                    numericInput("reactionsPerKit", label = "Reactions Per Kit",
                                 value = 6, step = 1, min = 1),
                    bsPopover("reactionsPerKit", title = "Reactions Per Kit", placement = "top", options = list(container = "body"),
                              content = "Number of reactions/lanes on one 10X kit (different kit versions possible)"),

                    hr(),
                    h5("Multiple testing correction"),
                    numericInput("pval",label="P-value",
                                 value=0.05,step=0.01,min=0,max=1),
                    bsPopover("pval", title="P-value", placement="top", options = list(container = "body"),
                              content="Significance threshold"),
                    selectInput("MTmethod", label = "Multiple testing method",
                                choices = list("FWER"="Bonferroni",
                                               "FDR"="FDR",
                                               "none"="none"),
                                selected="FDR"),
                    bsPopover("MTmethod", title="Multiple testing method", placement="top", options = list(container = "body"),
                              content="Possible is adjustment after the family-wise error rate (FWER), after the false discovery rate (FDR) or no adjustment at all (none)."),
                    numericInput("indepsnps", label="Independent SNPs", value=10, min=1, step=1),
                    bsPopover("indepsnps", title="Independent SNPs", placement="top",
                                           content="Number of independent SNPs assumed for each locus (for eQTL Bonferroni multiple testing correction the number of tests are estimated as number expressed genes * indepSNPs)")

                  ),
                  box(width = 0,
                    title="Mapping and Multiplet estimation",
                    id="multiplet",
                    solidHeader = TRUE,
                    status="orange",
                    numericInput("map.eff", label = "Mapping efficiency",
                                 value = 0.8,step=0.05,min=0,max=1),
                    numericInput("multipletRate", label = "Multiplet rate",
                                 value = 7.67e-06,step=1e-6,min=0),
                    bsPopover("multipletRate", title="Multiplet rate", placement="top", options = list(container = "body"),
                              content=paste("Rate factor to calculate the number of multiplets dependent on the number of cells loaded per lane.",
                              "We assume a linear relationship of multiplet fraction = cells per lane * multiplet rate.")),
                    numericInput("multipletFactor", label = "Multiplet factor",
                                 value = 1.82, step=0.1,min=1),
                    bsPopover("multipletFactor", title="Multiplet factor", placement="top", options = list(container = "body"),
                              content="Multiplets have a higher fraction of reads per cell than singlets, the multiplet factor shows the ratio between the reads."),
                    hr(),
                    h5("Expression cutoffs"),
                    numericInput("minUMI", label = "Minimal number of UMI per gene",
                                 value = 3, step=1,min=1),
                    bsPopover(" ", title=" ", placement="top", options = list(container = "body"),
                              content=" "),
                    numericInput("percIndiv", label = "Fraction of individuals",
                                 value = 0.5,step=0.05,min=0,max=1),
                    bsPopover(" ", title=" ", placement="top", options = list(container = "body"),
                              content=" "),
                    hr(),
                    h5("Special parameters"),
                    checkboxInput("speedCalc", "Skip power for lowly expressed genes", value = FALSE),
                    bsPopover("speedCalc", title="Skip power for lowly expressed genes:", placement="top", options = list(container = "body"),
                              content="This parameter will speed the calculation by setting the power of lowly expressed genes (probability smaller than 0.01) automatically to 0.
                              This will have only a marginal effect on the overall power, but of course change the DE/eQTL power estimates."),
                    checkboxInput("simPower", "Use simulated power for eQTLs",value=FALSE),
                    bsPopover("simPower", title="Use simulated power for eQTLs:", placement="top", options = list(container = "body"),
                              content="For genes with small mean values, the method used for eQTL power calculation can get inaccurate.
                              Instead the eQTL power for these small mean values can be estimated via simulation, which however increases the runtime.")
                    )
                  ),

                column(width=8,
                  box(width=0,
                    title="Detection power depending on design parameters",
                    solidHeader = TRUE,
                    status="primary",
                    HTML("Detection power depending on <strong>cells per individual</strong>, <strong>read depth</strong> and <strong>sample size</strong>.",
                    "Display two of those three parameters as x- and y-axis by selecting from the options in <strong>'Parameter grid'</strong>, the third one will be displayed as <strong>circle size</strong>."),
                    br(),
                    br(),
                    HTML("Click the <strong>Calculate optimal study</strong> button to update the plots with the current set of parameters."),
                    #--------------------------
                    plotlyOutput("powerPlot"),
                    #--------------------------
                    br(),
                    HTML("<strong>Click</strong> on a specific point in the plot to visualize the exact trace in the plots below"),
                    br()
                    ),
                  box(width=0,
                    title="Influence of design parameters on individual power components",
                    solidHeader = TRUE,
                    status="primary",
                    HTML("The overall detection power is the result of <strong>expression probability</strong>
                         (probability that the DE/eQTL genes are detected) and <strong>DE power</strong>
                         (probability that the DE/eQTL genes are found significant)."),
                    br(),
                    HTML("Below a visualization how the design choices influence those power components."),
                    plotlyOutput("readPlot"),
                    br(),
                    div("The plots show the influence of the y axis (left) and x axis (right) parameter of the upper plot onto the power of the selected study, while keeping the second parameter constant."),
                    br(),
                    HTML("The dashed lines shows the location of the <strong>selected study</strong>.")
                    )
                  )
                )
              ),

      tabItem(tabName = "celltypes",
              fluidRow(

                box(
                  title="Study parameters",
                  solidHeader = TRUE,
                  status="orange",
                  numericInput("numSamples", label = "Samples", value = 50, step=10,min=0),
                  bsPopover("numSamples", title="Samples: ", placement="top", options = list(container = "body"),
                            content="Number of samples of the study"),
                  numericInput("ctFreq", label = "Cell type frequency",
                               value = 0.1,step=0.05,min=0,max=1),
                  bsPopover("ctFreq", title="Cell type frequency:", placement="top", options = list(container = "body"),
                            content="Cell type frequency of the cell type of interest"),
                  numericInput("probCT", label = "Detection power",
                               value = 0.95,step=0.05),
                  bsPopover("probCT", title="Detection power:", placement="top", options = list(container = "body"),
                            content="The target power that should be reached with the parameter combination"),
                  sliderInput("cellsCT", label = "Minimal number of cells", min = 1,
                              max = 100, value = c(10, 50)),
                  bsPopover("cellsCT", title="Minimal number of cells:", placement="top", options = list(container = "body"),
                            content="How many cells of the target cell type should at least be detected for each individual"),
                  hr(),
                  width = 4
                  ),
                box(
                  title="Required cells per person to detect rare cell types with a certain power",
                  solidHeader = TRUE,
                  status="primary",
                  plotlyOutput("freqPlot"),
                  br(),
                  div("The figure shows the required number of cells per individual (y-axis, log scale)
                  to detect the minimal number of cells from a target cell type per individuum (x-axis) with a certain
                  probability. The power depends on the total number of individuals and the frequency of the
                  target cell type. Note that the required number of cells per sample only counts
                  correctly measured cells (no doublets etc), so the number is a lower bound for the required cells to be sequenced.")
                  )
                )
              ),

      tabItem(tabName = "license",
        h2("The 3-Clause BSD License"),
        HTML("
          <p>Copyright &copy; 2023 Helmholtz Munich</p>
          <p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:</p>
          <ol>
            <li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.</li>
            <li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.</li>
            <li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.</li>
          </ol>
          <p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
        ")
      )
      )
    )


footer <-  tags$footer(class = "main-footer",
                       # leftaligned part

                       p("This tool was developed by the",
                         a("Heinig lab",href="https://www.helmholtz-muenchen.de/icb/research/labs/genetic-and-epigenetic-gene-regulation/overview/index.html"),
                         "at the",
                         a("Helmholtz Center Munich.",href="https://www.helmholtz-muenchen.de/ueber-uns/service/kontakt/index.html"),
                         tags$br(),
                         "The software is also available as",
                         a("R package", href="https://github.com/heiniglab/scPower"),
                         "through Github.",
                         tags$br(),
                         a("Imprint",href="https://www.helmholtz-muenchen.de/impressum/index.html"), "-" ,
                         a("Privacy statement",href="https://www.helmholtz-muenchen.de/en/privacy-statement/index.html")
                         )
                       )

dashboardPage(
  header = header,
  sidebar = sidebar,
  body = body,
  footer = footer
  )
