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
                           fileInput("upload", paste0("Load counts data"), accept = c(".csv")),
                           
                         ), 
                         
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("Summary"),
                           
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
                           
                           tabPanel("Plot",
                                    
                           )
                           
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
  
  
  
  
  #' sample_data_table
  #'@details takes sample information data, displays it i a sortable table
  
  
  sample_data_table <- function(sum_data){
    
    sample_data_df <- sample_info_data()
    
    
    datatable(sample_data_df, rownames = FALSE, options = list(pageLength = nrow(sample_data_df)))
  }
  
  
  #' continuous_sample_histograms
  #'@details creates histogram plots of continuous variables in the input sample data frame
  
  
  continuous_sample_histograms <- function(sum_data){
    
    sample_data_df <- sample_info_data()
    
    # Extract continuous data and variables
    continuous_data_df <- sample_data_df %>%
      
      select_if(is.numeric)%>%
      
      gather() 
    
    print(continuous_data_df)
    
    
    ggplot(gather(continuous_data_df), aes(value, fill = key))+
      geom_histogram(bins = 10)+
      facet_wrap(~key, scales = "free_x")
    
    
  }
  
  
  
  #sample info summary tab
  output$sample_info_summary_table <- renderDT({
    sample_info_summary()
    
  })
  
  #sample data table
  output$sample_data <- renderDT({
    sample_data_table()
    
  })
  
  
  #sample summary histograms
  
  output$sample_hist <- renderPlot({
    
    continuous_sample_histograms()
    
  })
  
  
}



# Run the application 
shinyApp(ui = ui, server = server)
