library(shiny)
library(dplyr)
library(readr)
library(DT)
library(chemometrics)
library(mdatools)
library(factoextra)

server <- function(input, output, session) {
  
  #stop app when closing
  session$onSessionEnded(function() {
    stopApp()
  })
  
  getData <- reactive({
    if(input$filetype1=="labspec"){
      df <- as.data.frame(t(read.csv(input$file1$datapath,
                                     header = FALSE,
                                     sep = "\t",
                                     skip = 33, row.names = 1, stringsAsFactors = FALSE)))
      
      row.names(df)<-basename(as.vector(df[,1],"character"))
      
      df<-df[,-1]
      
      indx<-sapply(df, is.factor)
      
      df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
      
      df
      
    } else {
      df <- read.csv(input$file1$datapath)
      
      names <- basename(as.vector(df[,1],"character"))
      
      df <- df[,-c(1:6)]
      
      even <- seq_len(nrow(df)) %% 2
      
      colnames(df) <- as.vector(df[1,], "character")
      
      df <- data.frame(x=df[!even,])
      
      names <- names[!even]
      
      row.names(df) <- names
      
      df
    }
    
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
    
    df <- getData()
    
    # if(input$transpose == TRUE) {
    #   df2 <- t(df)
    # } else {
    #   df2 <- df
    # }
    
    #spectral data only
    #A <- as.matrix(df[,-c(1,2)])
    
    A <- df
    
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
    
    #savitzky-golay filtering
    if(input$filter == TRUE) {
      A <- prep.savgol(A, width = as.numeric(input$window), porder = as.numeric(input$porder), dorder = as.numeric(input$sgolay))
    }
    
    #PCA
    #(pca <- pca(A, ncomp = input$components))
    pca <- pca(A, ncomp = input$components)
    PCs <- pca(A)
    
    #get scores
    S <- pca$calres$scores
    
    #get factors
    F <- pca$loadings
    
    #get error
    #E <- pca$calres$residuals
    
    #sum of squared residuals
    #R <- rowSums(E)^2
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
    
    #using Moutlier
    #mdist <- Moutlier(Sr, plot = FALSE)
    
    #manual calculation
    M <- (t(Sr)%*%Sr)/(dim(A)[1]-1)
    
    #inverse mahalanobis matrix
    M_1 <- solve(M, tol = 1e-70)
    
    #
    D2 <- Sr %*% M_1 %*% t(Sr)
    #D2 <- S %*% M_1 %*% t(S)
    
    #
    D <- diag(D2)^(1/2)
    
    #mahalanobis w/RMSG
    # mdist2 <- mdist$md/((sum(mdist$md^2)/(dim(A)[1]-1)))^(1/2)
    
    RMSG <- ((sum(D^2)/(dim(A)[1]-1)))^(1/2)
    
    mdist2 <- D/RMSG
    
    #output <- cbind("Lot #" = df[,1],"M-Dist" = mdist$md, "Residual" = R)
    
    # output <- as.data.frame(cbind("Lot #" = as.character(df[,1]),"M-Dist" = mdist2, "Residual" = R))
    
    output <- as.data.frame(cbind("M-Dist" = mdist2, "Residual" = R))
    
    return(list(output,A, pca, df, Rc, M_1, RMSG, PCs))
    
  }) #end getMdist reactive
  
  #spectra output
  output$spectra <- renderPlot({
    
    df <- getData()
    
    # if(input$transpose == TRUE) {
    #   df2 <- t(df)
    # } else {
    #   df2 <- df
    # }
    
    #A <- df[,-1]
    
    A <- df
    
    #A <- t(A)
    
    mdaplot(A, type = 'l')
    
  })
  
  #screeplot output
  output$screeplot <- renderPlot({
    # A <- getMdist()[[2]]
    # pca <- prcomp(A)
    # fviz_eig(pca, addlabels = TRUE)
    plotVariance(getMdist()[[3]], type = 'b', show.labels = TRUE, labels = 'values')
    
  })
  
  #screeplot output
  output$residuals <- renderPlot({
    
    plotResiduals(getMdist()[[3]], show.labels = TRUE)
    
  })
  
  #data table output
  output$contents <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    # if(input$transpose == TRUE) {
    #   t(req(input$file1))
    # } else {
    #   req(input$file1) 
    # }
    
    req(input$file1)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        if(input$filetype1=="labspec"){
          df <- as.data.frame(t(read.csv(input$file1$datapath,
                                         header = FALSE,
                                         sep = "\t",
                                         skip = 33, row.names = 1, stringsAsFactors = FALSE)))
          
          row.names(df)<-basename(as.vector(df[,1],"character"))
          
          df<-df[,-1]
          
          indx<-sapply(df, is.factor)
          
          df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
          
          df
          
        } else {
          df <- read.csv(input$file1$datapath)
          
          names <- basename(as.vector(df[,1],"character"))
          
          df <- df[,-c(1:6)]
          
          even <- seq_len(nrow(df)) %% 2
          
          colnames(df) <- as.vector(df[1,], "character")
          
          df <- data.frame(x=df[!even,])
          
          names <- names[!even]
          
          row.names(df) <- names
          
          df
        }
        
        
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    return(df[,1:4])
    
    
    options = list(scrollX = TRUE)
    
  })
  
  #data download handler
  output$data1Download <- downloadHandler(
    filename = "Reference.csv",
    content = function(file) {
      write.csv(getMdist()[[1]], file, row.names = TRUE)
    }
  ) #end data download handler
  
  #_________________________________________________________________________________________________________________________________
  
  getData2 <- reactive({
    if(input$filetype2=="labspec"){
      df <- as.data.frame(t(read.csv(input$file2$datapath,
                                     header = FALSE,
                                     sep = "\t",
                                     skip = 33, row.names = 1, stringsAsFactors = FALSE)))
      
      row.names(df)<-basename(as.vector(df[,1],"character"))
      
      df<-df[,-1]
      
      indx<-sapply(df, is.factor)
      
      df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
      
      df
      
    } else {
      df <- read.csv(input$file2$datapath)
      
      names <- basename(as.vector(df[,1],"character"))
      
      df <- df[,-c(1:6)]
      
      even <- seq_len(nrow(df)) %% 2
      
      colnames(df) <- as.vector(df[1,], "character")
      
      df <- data.frame(x=df[!even,])
      
      names <- names[!even]
      
      row.names(df) <- names
      
      df
    }
    
    
  })
  
  #spectra output
  output$spectra2 <- renderPlot({
    
    df <- getData2()
    
    # if(input$transpose == TRUE) {
    #   df2 <- t(df)
    # } else {
    #   df2 <- df
    # }
    
    #A <- df[,-1]
    
    #A <- t(A)
    
    A <- df
    
    mdaplot(A, type = 'l')
    
  })
  
  #data table output
  output$contents2 <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    # if(input$transpose == TRUE) {
    #   t(req(input$file1))
    # } else {
    #   req(input$file1) 
    # }
    
    req(input$file2)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        if(input$filetype2=="labspec"){
          df <- as.data.frame(t(read.csv(input$file2$datapath,
                                         header = FALSE,
                                         sep = "\t",
                                         skip = 33, row.names = 1, stringsAsFactors = FALSE)))
          
          row.names(df)<-basename(as.vector(df[,1],"character"))
          
          df<-df[,-1]
          
          indx<-sapply(df, is.factor)
          
          df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))
          
          df
          
        } else {
          df <- read.csv(input$file2$datapath)
          
          names <- basename(as.vector(df[,1],"character"))
          
          df <- df[,-c(1:6)]
          
          even <- seq_len(nrow(df)) %% 2
          
          colnames(df) <- as.vector(df[1,], "character")
          
          df <- data.frame(x=df[!even,])
          
          names <- names[!even]
          
          row.names(df) <- names
          
          df
        }
        
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    return(df[,1:4])
    
  })
  
  #Analyze sample data
  #reactive
  getMdist2 <- reactive({
    
    #get reference data frame
    df <- getMdist()[[2]]
    
    #get existing PCA model
    pca <- getMdist()[[3]]
    
    #get scores
    S <- pca$calres$scores
    
    #get factors
    F <- pca$loadings
    
    #get error
    #E <- pca$calres$residuals
    
    #sum of squared residuals
    #R <- rowSums(E)^2
    R <- pca$calres$Q[,input$components]
    Rc <- scale(R, center = TRUE, scale = FALSE)
    
    #get new data
    newdata <- getData2()
    
    #spectral data only
    
    #A <- as.matrix(newdata[,-1])
    
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
    
    #mahalanobis using Moutlier
    #mdist <- Moutlier(A3, plot = FALSE)
    
    #get mahalanobis matrix
    M_1 <- getMdist()[[6]]
    
    #
    D2 <- A3 %*% M_1 %*% t(A3)
    
    #
    D <- diag(D2)^(1/2)
    
    #get RMSG
    RMSG <- getMdist()[[7]]
    
    #mahalanobis w/RMSG
    #mdist2 <- mdist$md/((sum(mdist$md^2)/(dim(Sr2)[1]-1)))^(1/2)
    mdist2 <- D/RMSG
    
    #output <- cbind("Lot #" = df[,1],"M-Dist" = mdist$md, "Residual" = R)
    
    lots<-rbind(df[,1],newdata[,1])
    
    #output <- as.data.frame(cbind("Lot #" = as.character(lots),"M-Dist" = mdist2, "Residual" = R2))
    #output <- as.data.frame(cbind("Lot #" = as.character(newdata[,1]),"M-Dist" = mdist2, "Residual" = E2))
    output <- as.data.frame(cbind("M-Dist" = mdist2, "Residual" = E2))
    
    return(list(output,A, pca))
    
  }) #end getMdist reactive
  
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
  
}