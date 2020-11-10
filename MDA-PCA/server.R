library(shiny)
library(dplyr)
library(readr)
library(DT)
library(chemometrics)
library(mdatools)
library(factoextra)

server <- function(input, output, session) {
  
  #increase file size uploads
  options(shiny.maxRequestSize=30*1024^2)
  
  #stop app when closing
  session$onSessionEnded(function() {
    stopApp()
  })
  
  observe({
    if(input$filetype1=="calfile"){
      shinyjs::disable("center")
      shinyjs::disable("scale")
      shinyjs::disable("residuals")
      shinyjs::disable("components")
      shinyjs::disable("mdistCI")
    } else {
      shinyjs::enable("center")
      shinyjs::enable("scale")
      shinyjs::enable("residuals")
      shinyjs::enable("components")
      shinyjs::enable("mdistCI")
    }
  })
  
  #_________________________________________________________________________________________________________________________________
  #REFERENCE DATA
  
  getRefData <- reactive({
    
    df <- getData(input$file1$datapath, input$file1$name, input$filetype1)[[1]]
    
    metadata <- getData(input$file1$datapath, input$file1$name, input$filetype1)[[2]]
    
    if(!is.null(metadata)) {
      updateCheckboxInput(session, "center", value = metadata$Center[1])
      updateCheckboxInput(session, "scale", value = metadata$Scale[1])
      updateCheckboxInput(session, "residuals", value = metadata$Residuals[1])
      updateNumericInput(session, "components", value = metadata$Factors[1])
      updateNumericInput(session, "mdistCI", value = metadata$CI[1])
    }
    
    if(input$bg == TRUE) {
      df2 <- getData(input$background$datapath, input$background$name, input$filetype1)[[1]]
      
      bgmeans <- colMeans(df2)
      
      for(i in 1:length(bgmeans)) df[,i] <- df[,i]/bgmeans[i]
      
      df
      
    }
    
    if(input$logtrans == TRUE) {
      df <- log(1/df)
    }
    
    df
    
  })
  
  
  #mdist output
  output$mdist <- renderDT({
    getMdist()[[1]]
  }) #end mdist output
  
  #pca output
  output$pcaSummary <- renderPrint({
    summary(getMdist()[[8]])
  }) #end output
  
  
  #reactive
  getMdist <- reactive({
    
    A <- getRefData()
    
    ctr <- if(input$center == TRUE) {
      "TRUE"
    } else {
      "FALSE"
    }
    
    scl <- if(input$scale == TRUE) {
      "TRUE"
    } else {
      "FALSE"
    }
    
    #scale data
    A <- prep.autoscale(A, center = as.logical(ctr), scale =  as.logical(scl))
    
    # #savitzky-golay filtering
    # if(input$filter == TRUE) {
    #   A <- prep.savgol(A, width = as.numeric(input$window), porder = as.numeric(input$porder), dorder = as.numeric(input$sgolay))
    # }
    
    #PCA
    pca <- pca(A, ncomp = input$components)
    PCs <- pca(A)
    
    #get scores
    S <- pca$calres$scores
    
    #get factors
    F <- pca$loadings
    
    #sum of squared residuals
    R <- pca$calres$Q[,input$components]
    
    #mean center error
    Rc <- scale(R, center = as.logical(ctr), scale = as.logical(scl))
    
    #append error
    Sr <- if(input$residuals == TRUE) {
      cbind(S[,1:input$components], Rc)
    } else {
      S[,1:input$components]
    }
    
    #mahalanobis distance
    M <- (t(Sr)%*%Sr)/(dim(A)[1]-1)
    
    #inverse mahalanobis matrix
    M_1 <- solve(M, tol = 1e-70)
    
    #squared m-distance
    D2 <- Sr %*% M_1 %*% t(Sr)
    
    #m-distance
    D <- diag(D2)^(1/2)
    
    #mahalanobis w/RMSG
    RMSG <- ((sum(D^2)/(dim(A)[1]-1)))^(1/2)
    mdist2 <- D/RMSG
    
    output <- as.data.frame(cbind("M-Dist" = mdist2, "Residual" = R))
    
    #M-distance cutoff calculation
    n <- dim(A)[1]
    k <- input$components
    
    f_crit <- qf(input$mdistCI, df1=k, df2=(n-k))
    Dm <- k*(n-1)
    Dma <- f_crit*Dm
    Dmax <- Dma/(n-k)
    DM <- sqrt(Dmax)
    
    cutoffDiff <- mdist2-DM
    PF <- ifelse(mdist2<DM, "Pass", "Fail")
    
    output <- as.data.frame(cbind("M-Dist" = mdist2, "Residual" = R, "Factors" = k, "Cutoff" = DM, 
                                  "Difference" = cutoffDiff, "Outcome" = PF))
    
    #Objects returned
    return(list(output,A, pca, df, Rc, M_1, RMSG, PCs, DM, k))
    
  }) #end getMdist reactive
  
  #spectra output
  output$spectra <- renderPlot(spectraPlot())
  
  spectraPlot <- function()({
    df <- getRefData()
    A <- df
    mdaplot(A, type = 'l')
  })
  
  #screeplot output
  output$screeplot <- renderPlot({
    plotVariance(getMdist()[[3]], type = 'b', show.labels = TRUE, labels = 'values')
    
  })
  
  #screeplot output
  output$residuals <- renderPlot({
    plotResiduals(getMdist()[[3]], show.labels = TRUE)
  })
  
  #data table output
  output$contents <- renderDT({
    df <- getRefData()
    return(df[,1:4])
    options = list(scrollX = TRUE)
  })
  
  #data download handler (m-dist)
  output$data1Download <- downloadHandler(
    filename = "Reference.csv",
    content = function(file) {
      write.csv(getMdist()[[1]], file, row.names = TRUE)
    }
  ) #end data download handler
  
  #data download handler (spectral data)
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("ref_spectral_data.csv")
    },
    content = function(file) {
      write.csv(getRefData(), file, row.names = TRUE)
    }
  ) #end data download handler
  
  #data download handler (calfile data)
  output$downloadCalFile <- downloadHandler(
    filename = function() {
      paste0("ref_calfile.csv")
    },
    content = function(file) {
      df <- getRefData()
      output <- as.data.frame(cbind(df, "Factors" = input$components, "Center" = input$center, 
                                    "Scale" = input$scale, "Residuals" = input$residuals, "CI" = input$mdistCI))
      write.csv(output, file, row.names = TRUE)
    }
  ) #end data download handler
  
  #plot download handler (spectral plot)
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("ref_spectral_plot.png")
    },
    content = function(file) {
      png(file, width=1000)
      spectraPlot()
      dev.off()
    }
  ) #end plot download handler
  
  #_________________________________________________________________________________________________________________________________
  #SAMPLE DATA
  
  getSampleData <- reactive({
    df <- getData(input$file2$datapath, input$file2$name, input$filetype2)[[1]]
    
    if(input$bg2 == TRUE) {
      df2 <- getData(input$background2$datapath, input$background2$name, input$filetype2)[[1]]
      
      bgmeans <- colMeans(df2)
      
      for(i in 1:length(bgmeans)) df[,i] <- df[,i]/bgmeans[i]
      
      df
    }
    
    if(input$logtrans2 == TRUE) {
      df <- log(1/df)
    }
    
    df
    
  })
  
  #spectra output
  output$spectra2 <- renderPlot(spectraPlot2())
  
  spectraPlot2 <- function()({
    A <- getSampleData()
    #A <- df
    mdaplot(A, type = 'l')
  })
  
  #data table output
  output$contents2 <- renderDT({
    df <- getData(input$file2$datapath, input$file2$name, input$filetype2)[[1]]
    return(df[,1:4])
    options = list(scrollX = TRUE)
  })
  
  #Analyze sample data
  getMdist2 <- reactive({
    
    #get reference data frame
    df <- getMdist()[[2]]
    
    #get existing PCA model
    pca <- getMdist()[[3]]
    
    #get scores
    S <- pca$calres$scores
    
    #get factors
    F <- pca$loadings
    
    #sum of squared residuals
    R <- pca$calres$Q[,input$components]
    Rc <- scale(R, center = TRUE, scale = FALSE)
    
    #get new data
    newdata <- getSampleData()
    
    A <- newdata
    
    
    #scale new data
    center <- attr(df,"prep:center")
    scale <- attr(df,"prep:scale")
    
    A2 <- scale(A, center, scale) %*% pca$loadings 
    
    
    #get error
    E <- scale(A,center, scale) - A2%*%t(pca$loadings)
    
    
    #subtract mean of training group residuals
    E2 <- rowSums(E^2)
    E3 <- E2-attr(Rc,"scale")
    
    #append error to new spectral data
    A3 <- cbind(A2,E3)
    
    #get mahalanobis matrix
    M_1 <- getMdist()[[6]]
    
    #
    D2 <- A3 %*% M_1 %*% t(A3)
    
    #
    D <- diag(D2)^(1/2)
    
    #get RMSG
    RMSG <- getMdist()[[7]]
    
    #mahalanobis w/RMSG
    mdist2 <- D/RMSG
    
    DM <- getMdist()[[9]]
    
    k <- getMdist()[[10]]
    
    cutoffDiff <- mdist2-DM
    PF <- ifelse(mdist2<DM, "Pass", "Fail")
    
    output <- as.data.frame(cbind("M-Dist" = mdist2, "Residual" = E2, "Factors" = k, "Cutoff" = DM, "Difference" = cutoffDiff, "Outcome" = PF))
    
    return(list(output,A, pca))
    
  }) #end getMdist2 reactive
  
  #mdist output
  output$mdist2 <- renderDT({
    getMdist2()[[1]]
  }) #end mdist output  
  
  #data download handler
  output$data2Download <- downloadHandler(
    filename = "Sample.csv",
    content = function(file) {
      write.csv(getMdist2()[[1]], file, row.names = TRUE)
    }
  ) #end data download handler
  
  #data download handler (spectral data)
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0("sample_spectral_data.csv")
    },
    content = function(file) {
      write.csv(getSampleData(), file, row.names = TRUE)
    }
  ) #end data download handler
  
  #plot download handler
  output$downloadPlot2 <- downloadHandler(
    filename = function() {
      paste0("sample_spectral_plot.png")
    },
    content = function(file) {
      png(file, width=1000)
      spectraPlot2()
      dev.off()
    }
  ) #end plot download handler
  
  
  #_____________________________________________________________________________________
  #GENERIC FUNCTIONS
  
  mycsvs<- function(files)({
    rbindlist(lapply(files, fread), idcol="file")
  })
  
  
  getRawData <- function(dataset, files)({
    df <- dataset
    df <- as.data.frame(df)
    
    if('scan_sample' %in% df[1,2]) {
      df <- df[-1,]
    } else {
      df
    }
    
    df <- df[!(df[,2]=="scan_sample"),]
    
    #get sample and file names
    file_names <- basename(as.vector(df[,1],"character"))
    sample_names <- basename(as.vector(df[,2],"character"))
    
    #Identify non-numeric columns and delete
    colspace <- c()
    for(i in 1:length(df[1,])) colspace[i] <- ifelse((is.factor(df[1,i])|is.na(df[1,i])),1,0)
    colnums <- sum(colspace)
    df <- df[,-c(1:colnums)]
    
    #delete non-spectra rows
    even <- seq_len(nrow(df)) %% 2
    colnames(df) <- as.vector(df[1,], "character")
    df <- df[seq(2,nrow(df),2),]
    sample_names <- sample_names[seq(2,length(sample_names),2)]
    file_names <- file_names[seq(1,length(file_names),2)]
    
    #extract lot numbers from file names
    df$file <- file_names
    df$id <- factor(df$file, label=files)
    df$id <- paste0(gsub(".*_","",gsub(".csv","",df$id)),"_",sample_names)
    
    #delete non-numeric columns, except id
    drops <- c("","file")
    df <- df[,!(names(df) %in% drops)]
    
    #get colmeans for spectral data by sample
    df <- df %>% group_by(id) %>% summarise_if(is.numeric, ~mean(.x))
    
    #delete id column and use as rownames
    id <- df$id
    df <- df %>% select(-id)
    df <- as.data.frame(df)
    rownames(df) <- id
    
    #return df object
    df 
  })
  
  getData <- function(file, filenames, filetype)({
    if(filetype=="labspec"){
      df <- as.data.frame(t(read.csv(file,
                                     header = FALSE,
                                     sep = "\t",
                                     skip = 33, row.names = 1, stringsAsFactors = FALSE)))
      
      #complete cases
      df <- df[complete.cases(df),]
      
      row.names(df)<-basename(as.vector(df[,1],"character"))
      
      df<-df[,-1]
      
      indx<-sapply(df, is.factor)
      
      df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
      
      df_metadata <- c()
      
      df
      
    } else if(filetype=="tellspec") {
      
      df <- mycsvs(file)
      df <- getRawData(df, filenames)
      
      df_metadata <- c()
      
    } else {
      
      df <- read.csv(file)
      
      metadata <- c("Factors", "Center", "Scale", "Residuals", "CI")
      
      df_metadata <- select(df, metadata)
      
      df <- select(df, -c(metadata))
      
      rownames(df) <- df[,1]
      
      df <- df[,-1]
      
    }
    
    return(list(df, df_metadata))
    
  })
  
}