library(shiny)
library(ggplot2)
library(plotly)
library(shinyWidgets)
library(d3heatmap)

ui <- fluidPage(
  titlePanel("Seq and Destroy"),
  
  tabsetPanel(type = "tabs",
              tabPanel("Summery",
                       textOutput("summery"),
                       plotlyOutput("rn4genes"),
                       plotlyOutput("rn6genes")),
              tabPanel("Table", 
                       sidebarPanel(
                         selectInput("dataset",
                                     label="choose the dataset",
                                     choices=c("deg4","deg6")),
                         numericInput("top",
                                      label = "enter the number of genes you want",
                                      value = 10)),
                       mainPanel(
                       tableOutput("table"))
                       ),
              
              tabPanel("volcano plot", plotOutput("volplot",
                                                  click = "plot_click",
                                                  hover = hoverOpts(
                                                    id = "plot_hover")),
                       verbatimTextOutput("click_info"),
                       verbatimTextOutput("hover_info"),
                       plotlyOutput("volplotly")
                       ),
              
              tabPanel("heatmap",
                       d3heatmapOutput("heatmap")
                       ),
              
              tabPanel("search for a gene",
                       sidebarPanel(
                         searchInput(
                           inputId = "search", label = "Enter your text",
                           placeholder = "A placeholder",
                           btnSearch = icon("search"),
                           btnReset = icon("remove"),
                           width = "450px")),
                       mainPanel(   
                       plotOutput("geneplot"),
                       verbatimTextOutput(outputId = "res"),
                       tableOutput("rest"))
              )
  )
)



server <- function(input, output) {

  ## Open input files 
  deg14 <- read.csv("rn4_tr_genesAll.csv", header = TRUE, sep = ",", row.names = 2)
  deg16 <- read.csv("rn6_tr_genesAll.csv", header = TRUE, sep = ",", row.names = 2)
  
  ## Summary tab
  output$summery <- renderText({
    "This is text test!"
  })
  
  output$rn4genes <- renderPlotly({
    plot_ly(data=data.frame('dep' = rownames(decideTests), decideTests$rn4_genes), 
            labels= ~dep, values= ~decideTests.rn4_genes, type = 'pie')
  })
  
  output$rn6genes <- renderPlotly({
    plot_ly(data=data.frame('dep' = rownames(decideTests), decideTests$rn6_genes), 
            labels= ~dep, values= ~decideTests.rn6_genes, type = 'pie')
  })
  
  ## Tabel tab
  datasetInput <- reactive({
    switch(input$dataset,
           "deg4" = deg14,
           "deg6" = deg16)
  })
  
  dataInput <- reactive({
    data.frame(datasetInput()[1:input$top,])
  })
  
  output$table <- renderTable(dataInput(), 
                              rownames = TRUE)
  
  ## Volcano plot tab
  ggpl2 <- reactive ({
    ggplot(degB) +
    geom_point(aes(x=logFC, y=-log10(adj.P.Val)),colour="red")+
    geom_point(aes(x=logFC.1, y=-log10(adj.P.Val.1)), colour="green")+
    labs(title="Hopfully a volcano plot", x="Fold Change", y="P value")
  })
  
  
  output$volplot <- renderPlot({
    ggply
  })
  
  output$volplotly <- renderPlotly({
    ggplotly(ggply )
  })
  
  ggply <- ggplot(deg16) +
    geom_point(aes(x=logFC, y=-log10(PValue)),colour="red")+
    labs(title="Hopefully a volcano plot", x="Fold Change", y="P value")
  
  
  
  output$click_info <- renderPrint({
    cat("input$plot_click:\n")
    str(input$plot_click)
  })
  
  output$hover_info <- renderPrint({
    cat("input$plot_hover:\n")
    str(input$plot_hover)
  })

  ## Heatmap tab
  output$heatmap <- renderD3heatmap({
    d3heatmap(matdf[1:15,1:5])
  })

  
  ## Gene search tab   
  output$geneplot <- renderPlot({
    ggplot(degB) +
      geom_point(aes(x=logFC, y=-log10(adj.P.Val)),colour="red")+
      geom_point(aes(x=logFC.1, y=-log10(adj.P.Val.1)), colour="green")+
      labs(title="Hopefully a volcano plot", x="Fold Change", y="P value")+
      geom_point(data = subset(degB, rownames(degB)==input$search) 
                 ,aes(x=logFC, y=-log10(adj.P.Val)),colour="black")+
      geom_point(data = subset(degB, rownames(degB)==input$search) 
                 ,aes(x=logFC.1, y=-log10(adj.P.Val.1)),colour="white")
  })
  
  searchInputGene <- reactive({
    data.frame(deg14[input$search,])
    data.frame(deg6[input$search,])
  })
  
  output$rest <- renderTable(searchInputGene(), 
                             rownames = TRUE)
  
  
}


shinyApp(ui = ui, server=server)

