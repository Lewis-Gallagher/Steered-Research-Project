# Load libraries
library(shiny)
library(ggplot2)
library(plotly)
library(shinyWidgets)
library(d3heatmap)


# Load data Bowtie pipeline
deg14u <- read.csv("rn4_genes_Unfiltered_uni.csv", header = TRUE, sep = ",", row.names = 2)
deg14 <- read.csv("rn4_genes_filtered_uni.csv", header = TRUE, sep = ",", row.names = 2)
deg16u <- read.csv("rn6_genes_Unfiltered_uni.csv", header = TRUE, sep = ",", row.names = 2)
deg16 <- read.csv("rn6_genes_filtered_uni.csv", header = TRUE, sep = ",", row.names = 2)
decideTests <- read.csv("Bowtie1decideTestsNoTotal.csv", header = TRUE, sep = ",", row.names = 1)
rowcounts <- read.csv("rn4_isoforms_normCounts.csv", header = TRUE, sep = ",", row.names = 2)
summaryTable <- read.csv("summary.csv", header = TRUE, sep = ",", row.names = 1)
summaryBar <- read.csv("summary_bar.csv", header = TRUE, sep = ",")

# Load data alternative pipeline
npdeg14u <- read.csv("np_rn4_geneid_result_uni.csv", header = TRUE, sep = ",", row.names = 1)
npdeg16u <- read.csv("np_rn6_geneid_result_uni.csv", header = TRUE, sep = ",", row.names = 1)
npdecideTests <- read.csv("tuxedo_table.csv", header = TRUE, sep = ",", row.names = 1)
npsummaryTable <- read.csv("tuxedo_table.csv", header = TRUE, sep = ",", row.names = 1)



# Define UI
ui <- fluidPage(
  
  tabsetPanel(type = "tabs",
              # The overview of the results tab
              tabPanel("Summary",
                       fluidRow(
                         style = "border:3px solid #f5f0f1;",
                         column(6,
                                h3("Overview of analysis results"),
                                style = "border:3px solid #f5f0f1;",
                                helpText(""),
                                tableOutput("rn46sum"),
                                plotlyOutput("sumbar")),
                         column(6,
                                style = "border:3px solid #f5f0f1",
                                fluidRow(h3("Library normalisation"),
                                         img(src = "rn4_genesLibSizes.png", height = 330, width = 460)),
                                fluidRow(h3("Controls test"),
                                         img(src = "control_correlation.png", height = 350, width = 490)))
                         ),
                       
                       #DEG
                       fluidRow(
                         h3("Differential expression figures"),
                         style = "border:3px solid #f5f0f1;",
                         column(6,
                                style = "border:3px solid #f5f0f1;",
                                h3("RGSC3.4.66"),
                                h4("Mean-difference plot"),
                                img(src = "rn4_genesMD.png", height = 400, width = 530),
                                h4("Differentially expressed genes percentage"),
                                plotlyOutput("rn4genes")),
                         column(6,
                                style = "border:3px solid #f5f0f1;",
                                h3("Rnor_6.0.92"),
                                h4("Mean-difference plot"),
                                img(src = "rn6_genesMD.png", height = 400, width = 530),
                                h4("Differentially expressed genes percentage"),
                                plotlyOutput("rn6genes"))
                         )
              ),
              # The gene list tab
              tabPanel("DEG Table", 
                       sidebarPanel(
                         selectInput("dataset",
                                     label="Select the refernce genome build",
                                     choices=c("RGSC3.4.66","Rnor_6.92")),
                         numericInput("top",
                                      label = "Number of genes shown",
                                      value = 30)),
                       mainPanel(
                         tableOutput("table"))
              ),
              # Gene expression plots tab
              tabPanel("DEG Plot", 
                       h3("RGSC3.4.66"),
                       selectInput("datasetrn4",
                                   label="Choose the dataset",
                                   choices=c("All genes","DEG")),
                       plotlyOutput("volplotlyrn4", width = "100%", height = "600px"),
                       h3("Rnor_6.0.92"),
                       selectInput("datasetrn6",
                                   label="Choose the dataset",
                                   choices=c("All genes","DEG")),
                       plotlyOutput("volplotlyrn6", width = "100%", height = "600px")
              ),
              # Heatmap tab
              tabPanel("Heatmap",
                       d3heatmapOutput("heatmap", width = "90%", height = "700px")
              ),
              # Gene search tab
              tabPanel("Gene Search",
                       fluidRow(
                         searchInput(
                           inputId = "search", label = "Search for a gene experssion",
                           placeholder = "Gene name in upper case",
                           btnSearch = icon("search"),
                           btnReset = icon("remove"),
                           width = "450px")),
                       fluidRow(
                         column(6,
                                style = "border:3px solid #f5f0f1;",
                                h3("RGSC3.4.66"),
                                tableOutput("resrn4"),
                                plotOutput("geneplotrn4")),
                         column(6,
                                style = "border:3px solid #f5f0f1;",
                                h3("Rnor_6.92"),
                                tableOutput("resrn6"),
                                plotOutput("geneplotrn6")))
            ),
            # Alternative pipeline menu
            navbarMenu("Alternative pipeline",
                       tabPanel("Summary",
                                fluidRow(
                                  style = "border:3px solid #f5f0f1;",
                                  h3("Library normalisation"),
                                  h4("RGSC3.4.66"),
                                  img(src = "rn4.png", height = 330, width = 460),
                                h4("RGSC6.0.92"),
                                img(src = "rn6.png", height = 330, width = 460)),
                       #DEG
                       fluidRow(
                         titlePanel("Differentially expressed genes percentage"),
                         style = "border:3px solid #f5f0f1;",
                         h4("RGSC3.4.66"),
                         plotlyOutput("nprn4genes"),
                         h4("RGSC6.0.92"), 
                         plotlyOutput("nprn6genes"))
                       ),

                       tabPanel("DEG Table", 
                                sidebarPanel(
                                  selectInput("npdataset",
                                              label="Select the reference genome build",
                                              choices=c("RGSC3.4.66","Rnor_6.92")),
                                  numericInput("nptop",
                                               label = "Number of genes shown",
                                               value = 30)),
                                mainPanel(
                                  tableOutput("nptable"))
                       ),
                       tabPanel("DEG Plot", 
                                h3("RGSC3.4.66"),
                                plotlyOutput("npvolplotlyrn4", width = "100%", height = "600px"),
                                h3("Rnor_6.92"),
                                plotlyOutput("npvolplotlyrn6", width = "100%", height = "600px")
                       ),
                       tabPanel("Gene Search",
                                fluidRow(
                                  searchInput(
                                    inputId = "npsearch", label = "Search for a gene experssion",
                                    placeholder = "Gene name in upper case",
                                    btnSearch = icon("search"),
                                    btnReset = icon("remove"),
                                    width = "450px")),
                                fluidRow(
                                  column(6,
                                         style = "border:3px solid #f5f0f1;",
                                         h3("RGSC3.4.66"),
                                         tableOutput("npresrn4"),
                                         plotOutput("npgeneplotrn4")),
                                  column(6,
                                         style = "border:3px solid #f5f0f1;",
                                         h3("Rnor_6.92"),
                                         tableOutput("npresrn6"),
                                         plotOutput("npgeneplotrn6")))
                       )
                       
     )
  )
)


# Define server function
server <- function(input, output) {
  
  ## Summary tab
  
  # Generate a bar chart
  output$sumbar <- renderPlotly({
    ggplotly(
      ggplot(summaryBar, 
             aes(fill=refgen, y=value, x=condition))+
        geom_bar(position="dodge", stat="identity")+
        labs(fill = "Genome" ,y="Number of genes", x= ""))
  })
  
  # Render summary table
  output$rn46sum <- renderTable(
    summaryTable, rownames = TRUE
  )
  
  # Generate pie charts for DGE %
  output$rn4genes <- renderPlotly({
    plot_ly(data=data.frame('dep' = rownames(decideTests), decideTests$rn4_genes), 
            labels= ~dep, values= ~decideTests$rn4_genes, type = 'pie', main = "Wha")
  })
  
  output$rn6genes <- renderPlotly({
    plot_ly(data=data.frame('dep' = rownames(decideTests), decideTests$rn6_genes), 
            labels= ~dep, values= ~decideTests.rn6_genes, type = 'pie')
  })
  
  
  ## Tabel tab
  
  # Reactive function to input the dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "RGSC3.4.66" = deg14,
           "Rnor_6.92" = deg16)
  })
  
  # Reactive function to select the number ofgenes to be displayed
  dataInput <- reactive({
    data.frame(datasetInput()[1:input$top,])
  })
  
  # Render the table
  output$table <- renderTable(dataInput(), 
                              rownames = TRUE)
  
  ## Volcano plot tab
  
  # Generate a ggplot for the datasets of genes both significant and non-signifiacnt (rn4) 
  ggplotlyrn4u <- ggplot(data = deg14u) +
    geom_point(data= deg14u,
               aes(x=logFC, y=-log10(PValue),label1=rownames(deg14u), label2=FDR),
               size = 0.3)+
    geom_hline(yintercept=2.19, linetype="dashed", color = "red")+
    geom_vline(xintercept=1.5, linetype="dashed", color = "red")+
    geom_vline(xintercept=-1.5, linetype="dashed", color = "red")+
    labs(title="Gene expression plot", x="Log Fold Change", y="Log P value")
  
  ggplotlyrn4 <- ggplot(data = deg14) +
    geom_point(data= deg14,
               aes(x=logFC, y=-log10(PValue),label1=rownames(deg14), label2=FDR),
               size = 0.3,
               colour="blue")+
    labs(title="Gene expression plot", x="Log Fold Change", y="Log P value")
  
  # Reactive function to select the dataset
  dataplottInput4 <- reactive({
    switch(input$datasetrn4,
           "All genes" = ggplotlyrn4u,
           "DEG" = ggplotlyrn4)
  })
  
  # Render the plot
  output$volplotlyrn4 <- renderPlotly({
    ggplotly(dataplottInput4())
  })
  
  # Generate a ggplot for the datasets of genes both significant and non-signifiacnt (rn6) 
  ggplotlyrn6u <- ggplot(data = deg16u) +
    geom_point(data= deg16u,
               aes(x=logFC, y=-log10(PValue),label1=rownames(deg16u), label2=FDR),
               size = 0.3)+
    geom_point(data= deg16,
               aes(x=logFC, y=-log10(PValue), label1=rownames(deg16), label2=FDR),
               size = 0.3)+
    geom_hline(yintercept=2.19, linetype="dashed", color = "red")+
    geom_vline(xintercept=1.5, linetype="dashed", color = "red")+
    geom_vline(xintercept=-1.5, linetype="dashed", color = "red")+
    labs(title="Gene expression plot", x="Log Fold Change", y="Log P value")
  
  ggplotlyrn6 <- ggplot(data = deg16) +
    geom_point(data= deg16,
               aes(x=logFC, y=-log10(PValue),label1=rownames(deg16), label2=FDR),
               size = 0.3)+
    geom_point(data= deg16,
               aes(x=logFC, y=-log10(PValue), label1=rownames(deg16), label2=FDR),
               size = 0.3)+
    labs(title="Gene expression plot", x="Log Fold Change", y="Log P value")
  
  # Reactive function to select the dataset
  dataplottInput6 <- reactive({
    switch(input$datasetrn6,
           "All genes" = ggplotlyrn6u,
           "DEG" = ggplotlyrn6)
  })
  
  # Render the plot
  output$volplotlyrn6 <- renderPlotly({
    ggplotly(dataplottInput6())
  })
  
  ## Heatmap tab
  output$heatmap <- renderD3heatmap({
    d3heatmap(rowcounts[1:30,3:7])
  })
  
  
  ## Search tab
  
  # Generate the DE plots
  output$geneplotrn4 <- renderPlot({
    ggplot(data = deg14u) +
      geom_point(data= deg14u,
                 aes(x=logFC, y=-log10(PValue),label1=rownames(deg14u), label2=FDR),
                 size = 0.3,
                 colour="black")+
      geom_hline(yintercept=2.19, linetype="dashed", color = "red")+
      geom_vline(xintercept=1.5, linetype="dashed", color = "red")+
      geom_vline(xintercept=-1.5, linetype="dashed", color = "red")+
      labs(title="Gene expression plot", x="Log Fold Change", y="Log P value")+
      geom_point(data = subset(deg14u, rownames(deg14u)==input$search) 
                 ,aes(x=logFC, y=-log10(PValue)),
                 colour="yellow")
  })
  
  output$geneplotrn6 <- renderPlot({
    ggplot(data = deg16u) +
      geom_point(data= deg16u,
                 aes(x=logFC, y=-log10(PValue),label1=rownames(deg16u), label2=FDR),
                 size = 0.3,
                 colour="black")+
      geom_hline(yintercept=2.19, linetype="dashed", color = "red")+
      geom_vline(xintercept=1.5, linetype="dashed", color = "red")+
      geom_vline(xintercept=-1.5, linetype="dashed", color = "red")+
      labs(title="Gene expression plot", x="Log Fold Change", y="Log P value")+
      geom_point(data = subset(deg16u, rownames(deg16u)==input$search) 
                 ,aes(x=logFC, y=-log10(PValue)),
                 colour="yellow")
  })
  
  # Render tables
  output$resrn4 <- renderTable(deg14u[input$search,], 
                               rownames = TRUE)
  
  output$resrn6 <- renderTable(deg16u[input$search,], 
                               rownames = TRUE)
  
  
  
### Alternative pipeline
  
  ## Summary tab
  output$nprn4genes <- renderPlotly({
    plot_ly(data=data.frame('dep' = rownames(npdecideTests), npdecideTests$rn4_genes), 
            labels= ~dep, values= ~npdecideTests$rn4_genes, type = 'pie', main = "Wha")
  })
  
  output$nprn6genes <- renderPlotly({
    plot_ly(data=data.frame('dep' = rownames(npdecideTests), npdecideTests$rn6_genes), 
            labels= ~dep, values= ~npdecideTests.rn6_genes, type = 'pie')
  })
  
  
  ## Tabel tab
  npdatasetInput <- reactive({
    switch(input$npdataset,
           "RGSC3.4.66" = npdeg14u,
           "Rnor_6.92" = npdeg16u)
  })
  
  npdataInput <- reactive({
    data.frame(npdatasetInput()[1:input$nptop,])
  })
  
  output$nptable <- renderTable(dataInput(), 
                              rownames = TRUE)
  
  ## Volcano plot tab
  output$npvolplotlyrn4 <- renderPlotly({
    ggplotly(ggplot(data = npdeg14u) +
               geom_point(data= npdeg14u,
                          aes(x=logFC, y=-log10(pval),label1=rownames(npdeg14u)),
                          size = 0.3)+
               labs(title="Hopefully a volcano plot for rn4", x="Fold Change", y="P value"))
  })
  
  output$npvolplotlyrn6 <- renderPlotly({
    ggplotly(ggplot(data = npdeg16u) +
               geom_point(data= npdeg16u,
                          aes(x=logFC, y=-log10(pval),label1=rownames(npdeg16u)),
                          size = 0.3)+
               labs(title="Hopefully a volcano plot for rn6", x="Log Fold Change", y="Log P value"))
  })
  
  ## Search tab
  output$npgeneplotrn4 <- renderPlot({
    ggplot(data = npdeg14u) +
      geom_point(data= npdeg14u,
                 aes(x=logFC, y=-log10(pval),label1=rownames(npdeg14u)),
                 size = 0.3)+
      labs(x="Fold Change", y="P value")+
      geom_point(data = subset(npdeg14u, rownames(npdeg14u)==input$search) 
                 ,aes(x=logFC, y=-log10(pval)),
                 colour="yellow")
  })
  
  output$npgeneplotrn6 <- renderPlot({
    ggplot(data = npdeg16u) +
      geom_point(data= npdeg16u,
                 aes(x=logFC, y=-log10(pval),label1=rownames(npdeg16u)),
                 size = 0.3)+
      labs( x="Log Fold Change", y="Log P value")+
      geom_point(data = subset(npdeg16u, rownames(npdeg16u)==input$search) 
                 ,aes(x=logFC, y=-log10(pval)),
                 colour="yellow")
  })
  
  output$npresrn4 <- renderTable(npdeg14u[input$npsearch,], 
                               rownames = TRUE)
  
  output$npresrn6 <- renderTable(npdeg16u[input$npsearch,], 
                               rownames = TRUE)
  
}


shinyApp(ui = ui, server=server)
