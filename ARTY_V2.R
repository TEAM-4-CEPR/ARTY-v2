library(shiny)
library(datasets)
library(plotly)
library(dplyr)
BiocManager::install("DESeq2")
library(DESeq2)
ui <- shinyUI(fluidPage(
  titlePanel("ARTY_V2"),
  tabsetPanel(
    tabPanel("Upload File",
             titlePanel("Uploading Files"),
             sidebarLayout(
               sidebarPanel(
                 fileInput('file1', 'Choisir un csv',
                           accept=c('text/csv', 
                                    'text/comma-separated-values,text/plain', 
                                    '.csv')),
                 
                 # added interface for uploading data from
                 # http://shiny.rstudio.com/gallery/file-upload.html
                 tags$br(),
                 checkboxInput('header', 'header', TRUE),
                 selectInput('xcol', 'X:', ""),
                 selectInput('ycol', 'Y', "", selected = ""),
                 selectInput("plot.type","Plot Type:",
                             list(bar = "bar", boxplot = "boxplot")#, histogram = "histogram", density = "density")
                 ),
                 radioButtons('sep', 'Separateur',
                              c(Semicolon=';',
                                Comma=',',
                                Tab='\t'),
                              ';'),
                 radioButtons('quote', 'Quote',
                              c(None='',
                                'Double Quote'='"',
                                'Single Quote'="'"),
                              '"')
                 
               ),
               mainPanel(
                 plotlyOutput('MyPlot')
                 
               )
             )
    ),
    tabPanel("PCA", titlePanel("This is the pca tab"),
             sidebarPanel(
               fileInput('file_pca', 'Choisir un csv',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv')),
               fileInput('file_group_pca', 'Choisir un csv',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv')),
               textInput( 'group_pca','Name of group for pca annotation', ""),
               radioButtons(inputId = "plot_type" , label = "Select the plot", choices = c("2d", "3d" ))),
             mainPanel(
               #tableOutput("test_pca")
               plotlyOutput('Myplot_pca')
               
             )),
    tabPanel("Heatmap", titlePanel("heatmap tab"),
             sidebarPanel(
               fileInput('file_heatmap', 'Choisir un csv',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv'))),mainPanel(
                                    plotlyOutput('MyPlot_heatmap')
                                    
                                  )),
    tabPanel("Design",titlePanel("Differential Expression Analysis"),
             textInput( 'group','Name of group', ""),
             textInput( 'replicate', 'Number of replicate per group',"",),
            
             "This is the DE tab" ,
             sidebarPanel(
             
              tableOutput('DE'),
              downloadButton("downloadData", "Download"),
              
              )),
    tabPanel("DE",titlePanel("Differential Expression Analysis"),
             sidebarPanel(
             fileInput('file2', 'Choisir un csv',
                       accept=c('text/csv', 
                                'text/comma-separated-values,text/plain', 
                                '.csv')),
             fileInput('file3', 'Choisir un csv',
                       accept=c('text/csv', 
                                'text/comma-separated-values,text/plain', 
                                '.csv')),
             downloadButton("downloadData_de", "Download"),),
             mainPanel(
             tableOutput("tableDE"))
    
    
  ), tabPanel("VolcanoPlot",titlePanel("VolcanoPlot"),
              sidebarPanel(
                fileInput('file_de', 'Choisir un csv',
                          accept=c('text/csv', 
                                   'text/comma-separated-values,text/plain', 
                                   '.csv')),
              ),mainPanel(
                #tableOutput("table_vol"),
                plotlyOutput('Myplot_volcano'))
                
              )
  )
))

server <- shinyServer(function(input, output, session) {
  # added "session" because updateSelectInput requires it
  options(warn=-1)
  # options(encoding="UTF-8")
  output$group <- renderText({ input$group })
  output$replicate <- renderText({input$replicate})
  data <- reactive({ 
    req(input$file1) ## ?req #  require that the input is available
    
    inFile <- input$file1 
    
    # tested with a following dataset: write.csv(mtcars, "mtcars.csv")
    # and                              write.csv(iris, "iris.csv")
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep,
                   quote = input$quote)
    for(i in seq_along(df)){
      df[,i] <- as.factor(df[,i])
      
    }
    df[df=="NA"] <- NA
    
    
   
    
    updateSelectInput(session, inputId = 'xcol', label = 'X',
                      choices = names(df), selected = names(df))
    updateSelectInput(session, inputId = 'ycol', label = 'Y',
                      choices = names(df), selected = names(df)[2])
    
    return(df)
  })
  volcano_table <- reactive({
    inFile_volcano <- input$file_de
    df_volcano <- read.csv(inFile_volcano$datapath)
    df_volcano <- mutate_all(df_volcano, function(x) as.numeric(as.character(x)))
    
    df_volcano <- mutate(df_volcano, condition = case_when(
      log2FoldChange > 0 & padj <= 0.05  ~ "upregulated_sig", 
      log2FoldChange < 0 & padj <= 0.05   ~ "downregulated_sig",
      padj > 0.05 ~ "non_sig"
    ))
    df_volcano$padj <- -log10(df_volcano$padj)
    return(df_volcano)
  })
  
  
  heatmap_table <- reactive({
    inFile_heatmap <- input$file_heatmap
    df_heatmap <- read.csv(inFile_heatmap$datapath , row.names="id")
    tt <- as.matrix(t(df_heatmap))
    return(tt)
  })
  pca_table <- reactive({
    inFile_pca <- input$file_pca
    inFile_group_pca <- input$file_group_pca
    df_pca <- read.csv(inFile_pca$datapath , row.names="id")
    df_group_pca <- read.csv(inFile_group_pca$datapath)
    input_group <- input$group_pca
    prin_comp <- prcomp(t(df_pca), rank. = 2)
    components <- prin_comp[["x"]]
    components <- data.frame(components)
    groups <- df_group_pca[,input_group]
    components <- cbind(components, groups)
    
    components$PC2 <- -components$PC2
    return(components)
    })
  
  pca_table_3d <- reactive({
    inFile_pca <- input$file_pca
    inFile_group_pca <- input$file_group_pca
    df_pca <- read.csv(inFile_pca$datapath , row.names="id")
    df_group_pca <- read.csv(inFile_group_pca$datapath)
    input_group <- input$group_pca
    prin_comp <- prcomp(t(df_pca), rank. = 3)
    components <- prin_comp[["x"]]
    components <- data.frame(components)
    groups <- df_group_pca[,input_group]
    components$PC2 <- -components$PC2
    components$PC3 <- -components$PC3
    components = cbind(components, groups)
    tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
    tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)
    return(components)
  })
    
  DE_table <- reactive({ 
    inFile2 <- input$file2 
    inFile3 <- input$file3
    df <- read.csv(inFile2$datapath, row.names="id")
    df1 <- read.csv(inFile3$datapath)
    dds <- DESeqDataSetFromMatrix(countData = df,
                                  colData = df1,
                                  design =  ~ sample)
    dds <- DESeq(dds)
    res <- results(dds)
    return(res)
    
  })
  output$test_pca <- renderTable({pca_table()} , rownames = TRUE)
  data2 <- reactive({data.frame(samplenames  = colnames(data())[-1],
                               sample = c(rep(strsplit(input$group,split = ',')[[1]][1],input$replicate) ,rep(strsplit(input$group,split = ',')[[1]][2] , input$replicate)),
                               replicate= c(seq(1:input$replicate) , seq(1:input$replicate)))
  })
  output$table_vol <- renderTable({volcano_table()})
  output$DE <- renderTable({data.frame(samplenames  = colnames(data())[-1],
                                                 sample = c(rep(strsplit(input$group,split = ',')[[1]][1],input$replicate) ,rep(strsplit(input$group,split = ',')[[1]][2] , input$replicate)),
                                                 replicate= c(seq(1:input$replicate) , seq(1:input$replicate)))
  
  })
  
  output$tableDE <- renderTable({DE_table()} , rownames = TRUE)
  
  output$downloadData <- downloadHandler(
      filename = function () {
        paste("Design.csv", sep = "")
      },
      content = function(file) {
        write.csv(data2(), file , row.names = FALSE)
      }
    )
  output$downloadData_de <- downloadHandler(
    filename = function () {
      paste("DE.csv", sep = "")
    },
    content = function(file) {
      write.csv(DE_table(), file , row.names = FALSE)
    }
  )
  output$Myplot_volcano <- renderPlotly({
    p <- plot_ly(data = volcano_table(), x = ~log2FoldChange, y = ~padj, text = ~id, mode = "markers", color = ~condition) %>% 
      layout(title ="Volcano Plot") 
    p
  })
  output$MyPlot_heatmap <- renderPlotly({
    p <- plot_ly(x=colnames(heatmap_table()), y=rownames(heatmap_table()), z = as.matrix(heatmap_table()), type = "heatmap") %>%
      layout(margin = list(l=120))
    p
  })
  output$Myplot_pca <- renderPlotly({
    if (input$plot_type == "2d") {
      fig <- plot_ly(pca_table(), x = ~PC1, y = ~PC2,  type = 'scatter', color = ~groups, colors = c('#636EFA','#EF553B'),mode = 'markers')%>%
        layout(
          legend=list(title=list(text='color')),
          plot_bgcolor='#e5ecf6',
          xaxis = list(
            title = "0",
            zerolinecolor = "#ffff",
            zerolinewidth = 2,
            gridcolor='#ffff'),
          yaxis = list(
            title = "1",
            zerolinecolor = "#ffff",
            zerolinewidth = 2,
            gridcolor='#ffff'))
      
      fig
    } else if (input$plot_type == "3d") {
      fig <- plot_ly(pca_table_3d(), x = ~PC1, y = ~PC2, z = ~PC3, color = ~groups, colors = c('#636EFA','#EF553B')) %>% add_markers(size = 12)
      fig <- fig %>%
        layout(
          title = 'Total Explained Variance = 99.48',
          scene = list(bgcolor = "#e5ecf6")
        )
      fig
    }
    
  })
  output$MyPlot <- renderPlotly({
    
    if(input$plot.type == "boxplot"){
      p <- ggplot(na.omit(data()), aes(x=na.omit(data())[,input$ycol],y=na.omit(data())[,input$xcol], fill =na.omit(data())[,input$ycol]))+
        geom_boxplot(na.rm = TRUE)+ ggtitle("Boxplot")+
        xlab(input$ycol)+ylab(input$xcol) +scale_fill_discrete(name=input$ycol)
      ggplotly(p) %>% layout(autosize=TRUE)# %>% style(hoverinfo = "none")
    }else{
      if(input$plot.type == "bar"){
        p <-  ggplot(data(),aes(x=as.factor(data()[,input$xcol]),fill=as.factor(data()[,input$ycol]),group = as.factor(data()[,input$ycol])))+
          geom_bar(aes(y=..prop..*100),position=position_dodge())+
          ggtitle("Pourcentage par Groupe")+
          xlab(input$xcol)+ylab("Pourcentage") +scale_fill_discrete(name=input$ycol)
        ggplotly(p) %>% layout(autosize=TRUE) %>% style(hoverinfo = "none")
        # } else {
        #   if(input$plot.type == "histogram"){
        #     p <-  ggplot(data(),aes(x=data()[,input$xcol],fill=data()[,input$ycol],group = data()[,input$ycol]))+
        #       geom_histogram()+
        #       ggtitle("Histogramm")+
        #       xlab(input$xcol) +scale_fill_discrete(name=input$ycol)
        #     ggplotly(p) %>% layout(autosize=TRUE) %>% style(hoverinfo = "none")
        #   }   Not that usefull ....
      }}
    
  })
  
  
  
  
})

shinyApp(ui, server)


