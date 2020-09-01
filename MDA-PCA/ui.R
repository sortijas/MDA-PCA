library(shinydashboard)
library(DT)
library(mdatools)
library(data.table)


#sidebar design
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Reference Data", tabName = "reference", icon = icon("table")),
    menuItem("Sample Data", tabName = "sample", icon = icon("chart-line"))
  )
)

#_____________________________________________________________________________________


#body design
body <- dashboardBody(
  tabItems(
    
    tabItem(tabName = "reference",
            tabsetPanel(
              tabPanel("Upload",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Upload Reference Data",
                                    # Input: File type ----
                                    radioButtons("filetype1", "File type",
                                                 choices = c(LabSpec = "labspec",
                                                             Tellspec = "tellspec",
                                                             CalFile = "calfile"
                                                 ),
                                                 selected = "labspec"),
                                    
                                    # Input: Select a file ----
                                    fileInput("file1", "Choose File",
                                              multiple = TRUE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    
                                    # Add background file ----
                                    conditionalPanel(
                                      condition = "input.filetype1 == 'tellspec'",
                                      checkboxInput("bg", tags$b("Background Correct")),
                                      conditionalPanel(
                                        condition = "input.bg == '1'",
                                        fileInput("background", "Background File",
                                                  multiple = FALSE,
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv"))
                                      ),
                                      checkboxInput("logtrans", tags$b("Log Transform"))  
                                    ),
                                    
                                    # Horizontal line ----
                                    tags$hr()
                                    
                                )#end box
                                
                         ), #end column 1
                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),
                                
                                box(width=14,
                                    plotOutput("spectra"),
                                    downloadButton("downloadPlot", "Plot"),
                                    downloadButton("downloadData", "Data")
                                ),
                                
                                DTOutput("contents")
                                
                         ) #end column 2
                       ) #end fluidRow
              ), #end tabPanel
              
              tabPanel("PCA",
                       fluidRow(
                         shinyjs::useShinyjs(),
                         column(3,
                                box(width = 12, title = "PCA Parameters",
                                    checkboxInput("center", tags$b("Center Data"), value = TRUE),
                                    checkboxInput("scale", tags$b("Scale Data"), value = FALSE),
                                    checkboxInput("residuals", tags$b("Use Residuals"), value = TRUE),
                                    numericInput("components", tags$b("Principal Components"), value = 5),
                                    downloadButton("downloadCalFile", "Cal File")
                                )#end box
                                
                         ), #end column 1
                         
                         column(3,
                                verbatimTextOutput("pcaSummary")
                         ),
                         
                         column(5,
                                plotOutput("screeplot"),
                                br(),
                                plotOutput("residuals")
                         )
                         
                       ) #end fluidRow
              ) #end tabPanel
              ,
              
              tabPanel("M-Distance",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "M-Dist Parameters",
                                    numericInput("mdistCI", tags$b("Confidence Interval"), value = 0.95, max = 0.99, step = .01)
                                )#end box
                                ,
                                box(width = 12, title = "Download Analysis",
                                    downloadButton("data1Download", "Download")
                                )#end box
                                
                         ), #end column 1
                         
                         column(6,
                                # plotOutput("screeplot"),
                                # br(),
                                # plotOutput("residuals"),
                                DTOutput("mdist")
                         )
                         
                       ) #end fluidRow
              ) #end tabPanel
            ) #end tabsetPanel
    ), #end tabItem
    
    #_____________________________________________________________________________________
    
    #body design
    
    tabItem(tabName = "sample",
            tabsetPanel(
              tabPanel("Upload",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Upload Sample Data",
                                    
                                    # Input: File type ----
                                    radioButtons("filetype2", "File type",
                                                 choices = c(LabSpec = "labspec",
                                                             Tellspec = "tellspec"
                                                 ),
                                                 selected = "labspec"),
                                    
                                    # Input: Select a file ----
                                    fileInput("file2", "Choose CSV File",
                                              multiple = TRUE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    checkboxInput("bg2", tags$b("Background Correct")),
                                    conditionalPanel(
                                      condition = "input.bg == '1'",
                                      fileInput("background2", "Background File",
                                                multiple = FALSE,
                                                accept = c("text/csv",
                                                           "text/comma-separated-values,text/plain",
                                                           ".csv"))
                                    ),
                                    checkboxInput("logtrans2", tags$b("Log Transform"))  
                                    
                                )#end box
                                
                         ), #end column 1
                         
                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),
                                
                                box(width=14,
                                    plotOutput("spectra2"),
                                    downloadButton("downloadPlot2", "Plot"),
                                    downloadButton("downloadData2", "Data")
                                ),
                                
                                DTOutput("contents2")
                                
                         ) #end column 2
                       ) #end fluidRow
              ), #end tabPanel
              
              tabPanel("M-Distance",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Download Analysis",
                                    downloadButton("data2Download", "Download")
                                )#end box
                         ), #end column 1
                         
                         column(6,
                                DTOutput("mdist2")
                         )
                         
                       ) #end fluidRow
              ) #end tabPanel
            ) #end tabsetPanel
    ) #end tabItem
  ) #end tabItems
) #end dasbhoardBody


#_____________________________________________________________________________________

#App arguments
dashboardPage(
  dashboardHeader(title = "PCA MDR"),
  sidebar,
  body
)
