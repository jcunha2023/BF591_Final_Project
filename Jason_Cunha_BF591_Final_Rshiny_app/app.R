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
library(patchwork) #used for creating multi-panel figures
library(pheatmap) #used for creating clustered counts heatmap

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
  h4("Comparisons of gene expression profiles between post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls (Labadorf et al., 2015)."),
  
  
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
                           
                           #adding two bars to select the passed and failed count filter colors
                           h5("Genes Passing Filters Color"),
                           colourInput("count_passed_filter_color",NULL, "#22577A"), #base point color
                           h5("Genes Not Passing Filters Color"),
                           colourInput("count_failed_filter_color", NULL, "#FFCF56"), #highlight point color
                           
                           #action button 
                           actionButton("do", "Filter", width = "100%", icon = icon("filter"), style = "color: black; background-color: #66cc99; border-color: #66cc99")
                           
                           
                         ), 
                         
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("Summary", DTOutput("count_info_summary_table")),
                           
                           tabPanel("Scatter Plots", plotOutput("count_summary_scatters")
                           ),
                           
                           tabPanel("Heatmap", plotOutput("count_heat")
                           ),
                           
                           tabPanel("PCA", fluidRow(sidebarPanel(
                                    
                                    
                                    #radio buttons to select which principal components to plot on x and y axes
                                    
                                    radioButtons("pca_x_axis", "Choose which principal component to plot on the x-axis",
                                                 c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                                   "PC7", "PC8", "PC9", "PC10"), selected = "PC1"),
                                    
                                    radioButtons("pca_y_axis", "Choose which principal component to plot on the y-axis",
                                                 c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                                   "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                                    
                                    #action button for pca plot
                                    #action button 
                                    actionButton("pca_plot_do", "Plot", width = "100%", icon = icon("brush"), style = "color: black; background-color: #6FCFD1; border-color: #6FCFD1")
                                    
                                    
                                    
                                    ),
                                    #Plot PCA
                                    column(8, plotOutput("count_pca"))
                           ))
                           
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
    count_df <- read.csv(input$count_info_upload$datapath, stringsAsFactors = FALSE, row.names = 1)
    
    return(count_df)
    
  })
  
  #' counts_summary_table
  #'@details creates table that summarizes effects of counts filtering, including:
  #'number of samples, total number of genes, number and % of genes passing current filter,
  #'number and % of genes not passing current filter

  
  
  counts_summary_table <- function(count_info_data, var_filter, count_filter){

    count_df <- counts_data() 
    
    #filtering step
    
    filtered_counts <- count_df %>%
      
      filter(
        rowSums(count_df != 0) > count_filter &
          apply(count_df, 1, function(row) quantile (row, prob = var_filter/100) > 0)
        
      )

  #create new tibble with all the info we want
  count_summary_tib <- tibble(

    "Number of Samples" = ncol(count_df),
    "Total Number of Genes" = nrow(count_df),
    "Number and % of Genes Passing Filters" = paste0(nrow(filtered_counts), " ( ", round((nrow(filtered_counts)/nrow(count_df))*100, 2), "%)"),
    "Number and % of Genes Not Passing Filters" =  paste0(nrow(count_df) - nrow(filtered_counts), " ( ", round(((nrow(count_df) - nrow(filtered_counts))/nrow(count_df))*100, 2), "%)")
  )

  datatable(count_summary_tib, rownames = FALSE)
}


#counts info summary tab
  output$count_info_summary_table <- renderDT({
    
    input$do #need to do this to link slider input to action button
    
    counts_summary_table(count_info_data = counts_data(), var_filter = isolate(input$var_range), count_filter = isolate(input$count_range))
    
})

  #' counts_diagnostic_scatterplots
  #'@details creates two scatter plots that show genes ranked by median count vs. log10(variance) and
  #' genes ranked by median count vs. number of zeros for each gene that passes the filters. Genes
  #' Passing filters are in darker colors, filtered genes in lighter color
  
  
  counts_diagnostic_scatterplots <- function(count_info_data, var_filter, count_filter, color_1, color_2){
    
    #creating a tibble storing the gene medians, variances, the sum of zero 
    #counts for each gene, and rankings
    
    result_tib <- count_info_data %>%
      
      tibble("gene_count_medians" = apply(select(.,-1), 1, median),
             "gene_variance" = apply(select(.,-1), 1, var),
             "sum_zero_counts" = apply(select(.,-1) == 0, 1, sum),
             "ranked_medians" = rank(gene_count_medians)
             
             
             )
    

    med_var_scatter <- ggplot(result_tib, aes(x = ranked_medians, y = log10(gene_variance)))+
      geom_point(aes(color = ifelse(rowSums(count_info_data != 0) > count_filter &
                                      apply(count_info_data, 1, function(row) quantile(row, prob = var_filter / 100) > 0),
                                    "Passed", "Filtered Out")))+
      scale_color_manual(values = c("Filtered Out" = color_2, "Passed"= color_1))+
      theme_bw()

    med_zeros_scatter <- ggplot(result_tib, aes(x = ranked_medians, y = sum_zero_counts))+
      geom_point(aes(color = ifelse(rowSums(count_info_data != 0) > count_filter &
                                      apply(count_info_data, 1, function(row) quantile(row, prob = var_filter / 100) > 0),
                                    "Passed", "Filtered Out")))+
      scale_color_manual(values = c("Filtered Out" = color_2, "Passed"= color_1))+
      theme_bw()

    med_var_scatter / med_zeros_scatter #create figure with both scatter plots

  }
  
  output$count_summary_scatters <- renderPlot({
    
    input$do #connects slider input to action buttion
    
    counts_diagnostic_scatterplots(count_info_data = counts_data(),
                                   var_filter = isolate(input$var_range), 
                                   count_filter = isolate(input$count_range),
                                   color_1 = isolate(input$count_passed_filter_color),
                                   color_2 = isolate(input$count_failed_filter_color))
    
  })
  
  #' counts_filtered_heatmap
  #'@details creates a clustered heatmap of counts remaining after filtering. Counts
  #'are log10-transformed before plotting


counts_filtered_heatmap <- function(count_info_data, var_filter, count_filter) {
  

#   #remove na's, filtering step, convert to matrix



  filtered_counts <- count_info_data %>%

    filter(
      rowSums(count_info_data != 0) > count_filter &
        apply(count_info_data, 1, function(row) quantile (row, prob = var_filter/100) > 0)

    ) 

  filtered_counts_matrix <- as.matrix(filtered_counts)

  #log10 transform the matrix

  # log_transformed_counts <- log10(filtered_counts)
  # 
  # any(is.na(log_transformed_counts))
  # any(is.infinite(log_transformed_counts))
  
  
  # Assuming filtered_counts is your original matrix
  
  # Subset the matrix to include all columns and the first 30 rows
  subset_counts <- filtered_counts_matrix[1:30, ]
  
  # Make heatmap with the subset
  heatmap(as.matrix(subset_counts), color = colorRampPalette(c("white", "blue"))(100), show_colnames = FALSE)
  


}


output$count_heat <- renderPlot({

  input$do #connects slider input to action button

  counts_filtered_heatmap(count_info_data = counts_data(),
                                 var_filter = isolate(input$var_range),
                                 count_filter = isolate(input$count_range))

})

  
  
  #' counts_pca_plot
  #'@details conducts principal components analysis (PCA) and generates a scatterplot
  #' of PCA results for filtered counts
  
  #'
  #' @param data tibble: a (n x _S_) data set
  #' @param meta tibble: sample-level meta information (_S_ x 3)
  #' @param title string: title for plot
  #'
  #' @return ggplot: scatter plot showing each sample in the first two PCs.
  #'
  #' @examples
  #' `plot_pca(data, meta, "Raw Count PCA")`
  
  counts_pca_plot <- function(count_info_data, x_axis, y_axis, title="") {
    
    #put samples in rows


      #do pca on transposed counts matrix (samples as rows, columns as genes)
    pca <- prcomp(
      t(count_info_data),
      center = TRUE,
      scale = FALSE
    )


    #getting pca output and metadata, merging them together by sample names

    pca_output <- data.frame(pca$x[,1:10])



    #incorporate variance_explained in x and y axis labels
    var_explained <- (pca$sdev)^2/sum((pca$sdev)^2) * 100
    
    #extract selected principal components for the x and y axes
    
    x_index <- as.numeric(sub("PC","",x_axis))
    y_index <- as.numeric(sub("PC","",y_axis))
    
    x_component <- var_explained[x_index]
    y_component <- var_explained[y_index]

    # #make plot

    pca_gg <- ggplot(pca_output, aes(x=!!sym(x_axis),y=!!sym(y_axis))) +
      geom_point(color = "#4FB0F5")+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      
      labs(x = paste0(x_axis,': ', round(x_component,1), '% variance'),
           y = paste0(y_axis, ': ', round(y_component,1), '% variance'),
           title = title)

    return(pca_gg)
  }
  
  
  output$count_pca <- renderPlot({

      input$pca_plot_do #connects slider input to action button

      counts_pca_plot(count_info_data = counts_data(), 
                      x_axis = isolate(input$pca_x_axis),
                      y_axis = isolate(input$pca_y_axis))

    })
  


}

# Run the application 
shinyApp(ui = ui, server = server)
