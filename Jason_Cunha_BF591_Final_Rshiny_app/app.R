#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


#data source: huntington's brain tissue samples: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810


library(tidyverse)
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this
library(DT) #used for interactive tables

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  #alter font for the entire webpage
  tags$head(
    tags$style(
      HTML("* { font-family: 'Montserrat', sans-serif; }"),
      HTML(".nav-tabs>li>a {color: #22577A}")
    )
  ),
  
  
  #page headings
  h1("BF591 Final Project"),
  h4("Comparisons of gene expression profiles between post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls (Labadorf et al. 2015)."),
  
  
  #creating separate main panels for plot and table
  
  tabsetPanel(id = "tabset",
              
              #Sample info Tab panel
              tabPanel("Samples",
                       
                       sidebarLayout(
                         sidebarPanel(
                           
                           tags$head(tags$style(".btn-file {background-color:#e86866;border-color: #e86866; 
                       }.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("sample_info_upload", paste0("Load sample information"), accept = c(".csv")),
                           
                          
                           
                           ), 
                         
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("Summary", DTOutput("sample_info_summary_table")),
                           
                           tabPanel("Table", DTOutput("sample_data")),
                           
                           tabPanel("Plots", plotOutput("sample_hist")),
                           
                         )) 
                       ),
                       
                       
                       
              ), 
              
              #Counts info Tab panel
              tabPanel("Counts",
                       
                       sidebarLayout(
                         sidebarPanel(
                           
                           tags$head(tags$style(".btn-file {background-color:#e86866;border-color: #e86866; 
                       }.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("count_info_upload", paste0("Load counts data"), accept = c(".csv")),
                           
                           #Slider to include genes with at least X percentile of variance
                           tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                           sliderInput("var_range", "Select the minimum percentile of variance:",
                                       min = 0, max = 100,
                                       value = 50),
                           
                           #Slider to include genes with at least X samples that are non-zero
                           tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                           sliderInput("count_range", "Select the non-zero sample threshold for genes:",
                                       min = 0, max = 69,
                                       value = 0),
                           
                           #action button 
                           actionButton("do", "Filter", width = "100%", icon = icon("filter"), style = "color: black; background-color: #66cc99; border-color: #66cc99")
                           
                           
                         ), 
                         
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("Summary", DTOutput("count_info_summary_table")),
                           
                           tabPanel("Scatter Plots",
                           ),
                           
                           tabPanel("Heatmap",
                           ),
                           
                           tabPanel("PCA",
                           )
                           
                         )) 
                       ),
                       
              ),
              
              #Differential Expression info Tab panel
              tabPanel("DE",
                       
                       sidebarLayout(
                         sidebarPanel(
                           
                           tags$head(tags$style(".btn-file {background-color:#e86866;border-color: #e86866; 
                       }.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("upload", paste0("Load differential expression results"), accept = c(".csv")),
                           
                         ), 
                         
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("DE Results"),
                           
                           tabPanel("Plot",)
                           
                         )) 
                       ),
                       
              ),
              
              #GSEA info Tab panel
              tabPanel("GSEA",
                       
                       sidebarLayout(
                         sidebarPanel(
                           
                           tags$head(tags$style(".btn-file {background-color:#e86866;border-color: #e86866; 
                       }.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("upload", paste0("Load FGSEA results"), accept = c(".csv")),
                           
                         ), 
                         
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("FGSEA Results Barplot"
                                    
                           ),
                           
                           tabPanel("Table",
                                    
                           ),
                           
                           tabPanel("NES Scatter Plot",
                                    
                           ),
                           
                         )) 
                       ),
                       
              ), 
  )
  
)



# Define server logic
server <- function(input, output) {

  options(shiny.maxRequestSize = 30 * 1024^2) #set file size limit high so you don't exceed limit
  
  ##SAMPLE INFO
  
  #' sample_info_data
  #'@details loads sample information dataframe
  
  sample_info_data <- reactive({
    
    req(input$sample_info_upload)
    df <- read.csv(input$sample_info_upload$datapath, stringsAsFactors = FALSE)
    
    return(df)
    
  })
  
  
  #' sample_info_summary
  #'@details takes sample information data, creates special summary table with means of continuous variables
  
  
  sample_info_summary <- function(sum_data){
    
    summary_df <- sample_info_data()%>%
      
      #manually set factor levels
      mutate(Condition = factor(Condition, levels = c("HD", "Control")))
    
    #create new tibble
    summary_tib <- tibble(
      
      "Column Name" = names(summary_df),
      "Type" = sapply(summary_df, class),
      "Mean (sd) or Distinct Values" = lapply(summary_df, function(x) {
        if (is.numeric(x) && !grepl("^mRNA_Seq_reads_", names(summary_df))) {
          mean_val <- round(mean(x, na.rm = TRUE), 2)
          sd_val <- round(sd(x, na.rm = TRUE), 2)
          paste0(mean_val, " (+/- ", sd_val, ")")
        } else {
          paste(unique(x), collapse = ", ")
        }
      })
    )
    
    
    
    
    datatable(summary_tib, rownames = FALSE, options = list(pageLength = nrow(summary_df)))
  }
  
  
  #sample info summary tab
  output$sample_info_summary_table <- renderDT({
    sample_info_summary()
    
  })
  
  
  
  #' sample_data_table
  #'@details takes sample information data, displays it i a sortable table
  
  
  sample_data_table <- function(sum_data){
    
    sample_data_df <- sample_info_data()
    
    
    datatable(sample_data_df, rownames = FALSE, options = list(pageLength = nrow(sample_data_df)))
  }
  
  #sample data table
  output$sample_data <- renderDT({
    sample_data_table()
    
  })
  
  #' continuous_sample_histograms
  #'@details creates histogram plots of continuous variables in the input sample data frame
  
  
  continuous_sample_histograms <- function(sum_data){
    
    sample_data_df <- sample_info_data()
    
    # Extract continuous data and variables
    continuous_data_df <- sample_data_df %>%
      
      select_if(is.numeric)%>%
      
      gather()
    
    
    ggplot(gather(continuous_data_df), aes(value, fill = key))+
      geom_histogram(bins = 10)+
      facet_wrap(~key, scales = "free_x")
    
    
  }
  
  
  #sample summary histograms
  
  output$sample_hist <- renderPlot({
    
    continuous_sample_histograms()
    
  })
  
  
  ##COUNTS INFORMATION
  
  #' counts_data
  #'@details loads RNA counts data frame
  
  counts_data <- reactive({
    
    req(input$count_info_upload)
    count_df <- read.csv(input$count_info_upload$datapath, stringsAsFactors = FALSE)
    
    return(count_df)
    
  })
  
  #' counts_summary_table
  #'@details creates table that summarizes effects of counts filtering, including:
  #'number of samples, total number of genes, number and % of genes passing current filter,
  #'number and % of genes not passing current filter

  counts_summary_table <- function(count_info_data, var_filter, count_filter){

    count_df <- counts_data() 
    
    
    filtered_counts <- count_df %>%
      
      filter(
        rowSums(count_df != 0) > count_filter &
          apply(count_df, 1, function(row) quantile (row, prob = var_filter/100) > 0)
        
      )

  #create new tibble with all the info we want
  count_summary_tib <- tibble(

    "Number of Samples" = ncol(filtered_counts),
    "Total Number of Genes" = nrow(filtered_counts)#,
    # "Number and % of Genes Passing Filters" = ,
    # "Number and % of Genes Not Passing Filters" =
  )

  datatable(count_summary_tib, rownames = FALSE)
}


#counts info summary tab
output$count_info_summary_table <- renderDT({
  
  input$do #need to do this to link slider input to action button
  
  counts_summary_table(count_info_data = counts_data(), var_filter = isolate(input$var_range), count_filter = isolate(input$count_range))
  
})



}

# Run the application 
shinyApp(ui = ui, server = server)
