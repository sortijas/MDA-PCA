library(shinydashboard)
library(DT)
library(mdatools)

#sidebar design
sidebar <- dashboardSidebar(
  sidebarMenu(
    # menuItem("About", tabName = "about", icon = icon("futbol")),
    # menuItem("Explore", tabName = "explore", icon = icon("chart-bar")),
    # menuItem("PCA", tabName = "PCA", icon = icon("chart-pie")),
    # menuItem("Model", tabName = "model", icon = icon("chart-line")),
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
                                    # Input: Select a file ----
                                    fileInput("file1", "Choose File",
                                              multiple = FALSE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    
                                    # Horizontal line ----
                                    tags$hr(),
                                    
                                    # Input: File type ----
                                    radioButtons("filetype1", "File type",
                                                 choices = c(LabSpec = "labspec",
                                                             Tellspec = "tellspec"
                                                 ),
                                                 selected = "labspec")
                                    
                                    # # Input: Checkbox if file has header ----
                                    # checkboxInput("header", "Header", TRUE),
                                    # 
                                    # # Input: Select separator ----
                                    # radioButtons("sep", "Separator",
                                    #              choices = c(Comma = ",",
                                    #                          Semicolon = ";",
                                    #                          Tab = "\t"),
                                    #              selected = ","),
                                    # 
                                    # # Input: Select quotes ----
                                    # radioButtons("quote", "Quote",
                                    #              choices = c(None = "",
                                    #                          "Double Quote" = '"',
                                    #                          "Single Quote" = "'"),
                                    #              selected = '"'),
                                    # 
                                    # # Input: Skip lines in record ----
                                    # numericInput("begin", "Data begins on line number",
                                    #              value = 0),
                                    # 
                                    # # Input: Checkbox if need to transpose dataset ----
                                    # checkboxInput("transpose", "Transpose", TRUE)
                                    
                                    
                                )#end box
                                
                                
                         ), #end column 1
                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),
                                plotOutput("spectra"),
                                br(),
                                DTOutput("contents")
                                
                         ) #end column 2
                         
                       ) #end fluidRow
                       
              ), #end tabPanel
              tabPanel("PCA",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "PCA Parameters",
                                    checkboxInput("center", tags$b("Center Data"), value = TRUE),
                                    checkboxInput("scale", tags$b("Scale Data"), value = FALSE),
                                    #checkboxInput("filter", tags$b("Savitzky-Golay Filter"), value = FALSE),
                                    conditionalPanel(
                                      condition = "input.filter == '1'",
                                      sliderInput("sgolay", "Derivative Order", min = 0, max = 2, value = 0),
                                      numericInput("porder", "Polynomial Order", value = 1),
                                      numericInput("window", tags$b("Window"), value = 15)
                                    ),
                                    checkboxInput("residuals", tags$b("Use Residuals"), value = TRUE),
                                    numericInput("components", tags$b("Principal Components"), value = 10)
                                )#end box
                                
                         ), #end column 1
                         column(3,
                                verbatimTextOutput("pcaSummary")
                         ),
                         column(5,
                                plotOutput("screeplot"),
                                br(),
                                plotOutput("residuals")
                                # DTOutput("mdist"),
                                # downloadButton("data1Download", "Download Analysis")
                                
                         )
                         # , #end column 2
                         # column(5,
                         #        plotOutput("screeplot"),
                         #        br(),
                         #        plotOutput("residuals")
                         # 
                         #        )
                         
                       ) #end fluidRow
                       
              ) #end tabPanel
              ,
              tabPanel("M-Distance",
                       fluidRow(
                         column(3,
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
                         # , #end column 2
                         # column(5,
                         #        plotOutput("screeplot"),
                         #        br(),
                         #        plotOutput("residuals")
                         # 
                         #        )
                         
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
                                    # Input: Select a file ----
                                    fileInput("file2", "Choose CSV File",
                                              multiple = FALSE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv"))
                                    ,
                                    
                                    # Input: File type ----
                                    radioButtons("filetype2", "File type",
                                                 choices = c(LabSpec = "labspec",
                                                             Tellspec = "tellspec"
                                                 ),
                                                 selected = "labspec")
                                    
                                    # # Horizontal line ----
                                    # tags$hr(),
                                    # 
                                    # # Input: Checkbox if file has header ----
                                    # checkboxInput("header2", "Header", TRUE),
                                    # 
                                    # # Input: Select separator ----
                                    # radioButtons("sep2", "Separator",
                                    #              choices = c(Comma = ",",
                                    #                          Semicolon = ";",
                                    #                          Tab = "\t"),
                                    #              selected = ","),
                                    # 
                                    # # Input: Select quotes ----
                                    # radioButtons("quote2", "Quote",
                                    #              choices = c(None = "",
                                    #                          "Double Quote" = '"',
                                    #                          "Single Quote" = "'"),
                                    #              selected = '"'),
                                    # 
                                    # # Input: Skip lines in record ----
                                    # numericInput("begin2", "Data begins on line number",
                                    #              value = 0)
                                    # 
                                    # # Input: Checkbox if need to transpose dataset ----
                                    # # checkboxInput("transpose", "Transpose", TRUE)
                                    
                                    
                                )#end box
                                
                         ), #end column 1
                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),
                                plotOutput("spectra2"),
                                br(),
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
                         # , #end column 2
                         # column(5,
                         #        plotOutput("screeplot2"),
                         #        br(),
                         #        plotOutput("residuals2")
                         #        
                         # )
                         
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
