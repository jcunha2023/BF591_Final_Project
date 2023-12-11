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
library('fgsea')
library('biomaRt')

#heatmap libraries (delete the ones you don't use)
library(pheatmap) #used for creating clustered counts heatmap
library(reshape) #used for "melting" counts data to make heatmap
library(ggplotify) ## to convert pheatmap to ggplot2
library(gplots) # used for heatmap.2()

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
                           fileInput("count_info_upload", paste0("Load counts data"), accept = c(".csv", ".tsv")),
                           
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
                           fileInput("deseq_upload", paste0("Load differential expression results"), accept = c(".csv", ".tsv")),
                           
                           
                           #Slider select the adjusted p-value filter
                           tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                           sliderInput("padj_range", "Select the Adjusted P-Value Filtering Threshold:",
                                       min = 0, max = 1, step = 0.01,
                                       value = 0.1),
                           actionButton("deseq_do", "Filter", width = "100%", icon = icon("filter"), style = "color: black; background-color: #66cc99; border-color: #66cc99")
                         ),  
                         
                         mainPanel(tabsetPanel(
                           
                           tabPanel("DE Results", DTOutput("deseq_tab")),
                           
                           tabPanel("Plots",
                                    
                                    tabsetPanel(
                                    tabPanel("Raw P-Values Histogram", plotOutput("pval_hist")),
                                    
                                    tabPanel("Log2 Fold-Changes Histogram", plotOutput("logfc_hist")),
                                    
                                    tabPanel("Top 10 Normalized Counts",),
                                    
                                    tabPanel("Volcano Plot of DE Results", plotOutput("volc_plot"))
                                             
                              
                                             )
                                    
                                    
                                    
                                    )
                           
                         )) 
                       ),
                       
              ),
              
              #GSEA info Tab panel
              tabPanel("GSEA",
                       
                       sidebarLayout(
                         sidebarPanel(
                           
                           h5(HTML('Please load the differential expression results in the DE tab before uploading the files below.')),
                           
                           tags$head(tags$style(".btn-file {background-color:#e86866;border-color: #e86866; 
                       }.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("id2gene_upload", paste0("Load EnsembID-to-MGI Symbol Mappings"), accept = c(".txt")),
                           
                           tags$head(tags$style(".btn-file {background-color:#e86866;border-color: #e86866; 
                       }.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("gene_set_path", paste0("Input the File Path of the Gene Sets of Interest, wrapped in quotation marks (.gmt format)", accept = c(".gmt")))
                           
                         ), 
                         
                         
                         mainPanel(tabsetPanel(

                           tabPanel("FGSEA Results Barplot", fluidRow(sidebarPanel(
                             
                
                             #Slider to select number of top pathways to plot based on adjusted p-value
                             tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                             sliderInput("fgsea_barplot_range", "Select the number of top pathways to plot based on adjusted p-value:",
                                         min = 1, max = 100,
                                         value = 10),
                             actionButton("fgsea_barplot_do", "Plot", width = "100%", icon = icon("brush"), style = "color: black; background-color: #66cc99; border-color: #66cc99")
                           
                
                           ), column(8, plotOutput("fgsea_top_paths")))
                                    
                           ),
                           
                           tabPanel("Table", fluidRow(sidebarPanel(
                             
                             
                             #Slider to filter table by adjusted pvalue
                             tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                             sliderInput("fgsea_table_range", "Select adjusted p-value to filter FGSEA results:",
                                         min = 0, max = 1,
                                         value = 0.1),
                            
                             
                             #radio buttons to select all, positive, or negative NES pathways
                             radioButtons("select_nes_pathways", "Choose which NES Pathways to display",
                                          c("All", "Positive", "Negative"), selected = "All"),
                             
                             #create numeric input for min and max fgsea values
                             numericInput("minval", "Minimum Gene Set Size to Test:", min = 1, value = 15),
                             numericInput("maxval" ,"Maximum Gene Set Size to Test:", min = 1, value = 500), 
                             
                             #filtering action button
                             actionButton("fgsea_table_do", "Filter", width = "100%", icon = icon("filter"), style = "color: black; background-color: #66cc99; border-color: #66cc99"),
                                          
                            #download button to export current filtered and displayed table results
                            
                            downloadLink('download_fgsea_table', 'Download filtered table as a .csv file')
                             
                             
                           ), column(8, DTOutput("fgsea_table")))
                                 
                                    
                        
                           ),
                           
                           tabPanel("NES Scatter Plot", fluidRow(sidebarPanel(
                             
                             
                             #Slider to filter table by adjusted pvalue
                             tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                             sliderInput("fgsea_scatter_range", "Select -log10(adjusted p-value) to filter FGSEA results:",
                                         min = 0, max = 70,
                                         value = 35),
                             actionButton("fgsea_scatter_do", "Filter", width = "100%", icon = icon("filter"), style = "color: black; background-color: #66cc99; border-color: #66cc99"),
                            
                             
                             
                           ), column(8, plotOutput("fgsea_scatter")))
                           
                           
                                    
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
  #'@details takes sample information data, displays it as a sortable table
  
  
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
      theme_bw()+
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
      
      tibble("gene_count_medians" = apply(dplyr::select(.,-1), 1, median),
             "gene_variance" = apply(dplyr::select(.,-1), 1, var),
             "sum_zero_counts" = apply(dplyr::select(.,-1) == 0, 1, sum),
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
  

#  filtering step, convert to matrix, do log10 transformation

  filtered_counts <- count_info_data %>%

    filter(
      rowSums(count_info_data != 0) > count_filter &
        apply(count_info_data, 1, function(row) quantile (row, prob = var_filter/100) > 0)

    ) 

  filtered_counts_matrix <- as.matrix(filtered_counts)
  
  # Reorder columns based on control and disease samples
  control_columns <- grep("^C_", colnames(filtered_counts_matrix))
  disease_columns <- grep("^H_", colnames(filtered_counts_matrix))
  all_columns <- c(control_columns, disease_columns)
  filtered_counts_matrix <- filtered_counts_matrix[, all_columns]
  
  
  # log10_filtered_matrix <- log10(filtered_counts_matrix)


  
  # Subset the matrix to include all columns and the first 30 rows
  subset_counts <- filtered_counts_matrix[1:30, ]
  

  
  #make heatmap 

  heatmap.2(subset_counts, 
            dendrogram = "column",
            key = TRUE,
            density.info = "none",
            trace = "none", 
            Colv = as.dendrogram(hclust(dist(subset_counts))),
            sepwidth = c(0.2, 0.2),         # Set width for column separations
            margins = c(10, 5),             # Adjust top and bottom margins
            main = "Heatmap of Counts",     # Add a main title
            cex.main = 2 ,                   # Set main title font size
            key.xlab = "Samples",  
            key.ylab = "Genes")
  
  # pheatmap(filtered_counts_matrix, scale = "row",
  #          color=colorRampPalette(c("white", "red"))(50))
  
  # #melt the data and add a factor column with HD and Control labels so we can make heatmap
  # 
  # melted_mat <- melt(subset_counts)
  # colnames(melted_mat) <- c("genes", "samples", "counts")
  # mutate(melted_mat, condition = ifelse(grepl("^H_", samples), "HD", "Control"))
  # 
  # 
  # 
  # # Make heatmap with the subset
  # # heatmap(filtered_counts_matrix, color = colorRampPalette(c("white", "blue"))(100), show_colnames = FALSE)
  # 
  # ggplot(melted_mat, aes(x = samples, y = genes, fill = log10(counts)))+
  #   geom_tile()+
  #   scale_x_discrete()
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #     
  #   

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
  
  
  
  ##DIFFERENTIAL EXPRESSION ANALYSIS RESULTS FUNCTIONS

  
  #' load_deseq_output
  #'@details loads the output statistics from a DESEQ2 differential expression analysis.
  
  load_deseq_output <- reactive({
    
    req(input$deseq_upload)
    deseq_df <- read.csv(input$deseq_upload$datapath, stringsAsFactors = FALSE, row.names = 1)
    
    return(deseq_df)
    
  })
  

  #' deseq_results_table
  #'@details takes DESEQ2 differential expression analysis results, adds an additional column
  #'indicating upregulation/downregulation/NS, and outputs a stortable table with gene search functionality
  
  
  deseq_results_table <- function(deseq_mat, padj_threshold){
    
    #load results
  

    #applying padj filter,adding additional up/down/ns column   
    annotated_deseq_res <- deseq_mat %>%
      
      filter(padj <= padj_threshold) %>% 
      
      mutate(Expression_Status = case_when(padj < padj_threshold & log2FoldChange > 0 ~ 'UP',
                                          padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN',
                                          TRUE ~ 'NS')) 

    
    datatable(annotated_deseq_res, rownames = TRUE, options = list(pageLength = nrow(annotated_deseq_res)))
  }
  
  #sample data table
  output$deseq_tab <- renderDT({
    
    
    input$deseq_do #link action button to adjusted pvalue slider
    
    deseq_results_table(deseq_mat = load_deseq_output(),
      padj_threshold = isolate(input$padj_range))
    
  })
  
  
  

    #' Function to plot the unadjusted p-values as a histogram
    #'
    #' @param labeled_results (tibble): Tibble with DESeq2 results
    #'
    #' @return ggplot: a histogram of the raw p-values from the DESeq2 results
    #' @export
    #'
    #' @examples pval_plot <- plot_pvals(labeled_results)
    deseq_pval_histogram <- function(deseq_results) {

      pval_gg <- ggplot(deseq_results, aes(x = pvalue)) +
        geom_histogram(bins = 60, color='black', fill='lightblue2')+
        theme_bw()+
        labs(x = 'P-Value', y = 'Count', title = 'Histogram of Raw P-Values Obtained from DE Analysis (HD vs. Control)')
      pval_gg


    }
  
    output$pval_hist<- renderPlot({
      
      deseq_pval_histogram(deseq_results = load_deseq_output())
      
    })
    
    

      #' Function to plot the log2foldchange from DESeq2 results in a histogram
      #'
      #' @param labeled_results (tibble): Tibble with DESeq2 results 
      #' @param padj_threshold (float): threshold for considering significance (padj)
      #'
      #' @return ggplot: a histogram of log2FC values from genes significant at padj
      #' threshold of 0.1
      #' @export
      #'
      #' @examples log2fc_plot <- plot_log2fc(deseq2_results, .10)
      plot_log2fc <- function(deseq2_results, padj_threshold) {

        #filter input to only include results within defined padj_threshold

        filtered_results <- subset(deseq2_results, padj <= padj_threshold)

        #plotting filtered results
        gg2 <- ggplot(filtered_results, aes(x = log2FoldChange)) +
          geom_histogram(bins = 100, color='black', fill='lightblue2')+
          theme_bw()+
          labs(x = 'log2FoldChange', y = 'Counts', title = 'Histogram of Log2FoldChanges for DE Genes (HD vs. Control)')
        gg2

      }
      
      output$logfc_hist<- renderPlot({
        
        input$deseq_do #links action button to pvalue slider
        
        plot_log2fc(deseq2_results = load_deseq_output(), 
                             padj_threshold = isolate(input$padj_range))
        
      })
      
      
      
      #' Function to generate volcano plot from DESeq2 results. Before plotting, the function
      #'
      #' @param labeled_results (tibble): Tibble with DESeq2 results
      #'
      #' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
      #'   -log10(padj) and labeled by status
      #' @export
      #'
      #' @examples volcano_plot <- plot_volcano(labeled_results)
      #'
      plot_volcano <- function(deseq2_results, padj_threshold) {

        volc_data <- deseq2_results %>%
          
          #adding expression status annotations based on padj threshold
          
          mutate(volc_plot_status = case_when(padj < padj_threshold & log2FoldChange > 0 ~ 'UP',
                                              padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN',
                                              TRUE ~ 'NS')) %>%
          

          #add -log10(padj) column to input tibble
          mutate(`-log10(adjusted p)`=-log10(padj),)


          gg4 <- ggplot(volc_data, aes(x=log2FoldChange,y=`-log10(adjusted p)`, color = volc_plot_status)) +
            geom_point()+
            theme_bw()+
            geom_hline(yintercept = 0, linetype = "dashed") +
            labs(x = 'log2FoldChange', y = '-log10(padj)', title = 'Volcano plot of DESeq2 differential expression results (HD vs. Control)')

          gg4


      }

      output$volc_plot<- renderPlot({
        
        input$deseq_do #links action button to pvalue slider
        
        plot_volcano(deseq2_results = load_deseq_output(),
                     padj_threshold = isolate(input$padj_range)
                     )
        
      })

      
      
      
    
      ##FGSEA GENE ENRICHMENT FUNCTIONS
      
      #' load_id2gene_file
      #'@details loads a text file with mappings between Ensembl IDs and MGI symbols
      
      load_id2gene_file <- reactive({
        
        req(input$id2gene_upload)
        
        
        id2gene <- read_delim(input$id2gene_upload$datapath, delim='\t', col_names = columns, show_col_types = FALSE)
        
        
        
        return(id2gene)

      })
      
      
      #' load_genes
      #'@details loads a .gmt file with gene sets of interest for fgsea

      load_gene_sets<- reactive({

        req(input$gene_set_path)
        
        
        pathwayLines <- strsplit(readLines(input$gene_set_path$datapath), "\t")
        pathways <- lapply(pathwayLines, tail, -2)
        names(pathways) <- sapply(pathwayLines, head, 1)
        
        return(pathways)

      })
      
      output$file_path_output <- renderText({
        
        load_gene_sets()
      })
      
      

      #' Function to generate fgsea results, using a vector ranked by log2FC descending
      #'
      #' @param labeled_results (tibble): Tibble with DESeq2 results 
      #' @param id2gene_file: Path to the file containing the mapping of
      #' ensembl IDs to MGI symbols

      #' @return Named vector with gene symbols as names, and log2FoldChange as values
      #' ranked in descending order
      #' @export
      #'
      #' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')
      
      
      run_fgsea <- function(deseq2_results, gmt_gene_sets, id2gene_file, min_size, max_size, padj_threshold, nes_filter) {
        
        
        #change column names in id2gene file
        colnames(id2gene_file) <- c('gene_symbols', 'ensemblids')


        #filter out NA rows in labeled results, merging labeled results with id2gene tibble
        lab_rank_res <- deseq2_results %>%
          drop_na() %>%
          rownames_to_column(., "genes") %>% #convert gene rownames to a new column for merging purposes
          
          mutate(genes = sub("\\.\\d+$", "", genes)) %>% #strip the version number from ensemblids (needed later on for merging)

          left_join(.,id2gene_file, by = c('genes' = 'ensemblids')) %>%

          mutate(log2fc_ranks = rank(log2FoldChange, ties.method='random')) %>%
          arrange(desc(log2fc_ranks))
        
        
        #create log2fc named vector, use as input for fgsea function
        
        log2fc_vec <- lab_rank_res$log2FoldChange
        
        #adding the gene symbols
        
        names(log2fc_vec) <- lab_rank_res$symbol
        

        
        #hallmark_pathways_fgsea <- fgsea::gmtPathways(gmt_file_path)
        
        #running fgsea
        
        fgsea_results <- fgsea(pathways = gmt_gene_sets, stats = log2fc_vec , minSize = min_size, maxSize= max_size)


        
        #filter input to only include results within defined padj_threshold
        
        fgsea_results <- subset(fgsea_results, padj <= padj_threshold)          
          
          #if/else statements to check if rows have positive or negative NES values
          
          if (nes_filter == "Positive") {
            fgsea_results <- subset(fgsea_results, NES > 0)
          } else if (nes_filter == "Negative") {
            fgsea_results <- subset(fgsea_results, NES < 0)
          }
          
          #convert to tibble, 
          
         fgsea_results <-  as_tibble(fgsea_results)
         

        return(fgsea_results)
        #datatable(fgsea_results, rownames = TRUE, options = list(pageLength = 50))
      }
        
      
      #create reactive expression to use the run_fgsea information
      
      fgsea_output <- reactive({
        
        input$fgsea_table_do #connect table functionality to action button

            run_fgsea(deseq2_results = load_deseq_output(),
                          gmt_gene_sets = load_gene_sets(),
                          id2gene_file = load_id2gene_file(),
                          min_size = isolate(input$minval),
                          max_size = isolate(input$maxval),
                          padj_threshold = isolate(input$fgsea_table_range),
                          nes_filter = isolate(input$select_nes_pathways))
          



      })
      
      

    #calling reactive fgsea variable to generate table
      
      output$fgsea_table <- renderDT({

        input$fgsea_table_do #connect table functionality to action button
        
        fgsea_output()


      })
      
      
      
    # #create a special fgsea output variable for downloading purposes
    #   
    #   data <- renderPrint({
    #     
    #     df <- apply(fgsea_output(), 2, as.character)
    #     
    #     return(df)
    #     
    #   })
    #   
    #   
    #   #add download functionality to run_fgsea
    # 
    #   output$download_fgsea_table <- downloadHandler(
    # 
    #     filename = function(){
    #       paste("filtered_fgsea_table", Sys.Date(), ".csv", sep = "_")
    #     },
    # 
    #     content = function(file){
    #       write.csv(apply(fgsea_output(), 2, as.character), row.names = FALSE)
    # 
    #     })
      
      
      
      


#       #add download functionality to run_fgsea
#       
#       output$download_fgsea_table <- downloadHandler(
#         
#         filename = function(){
#           paste("filtered_fgsea_table", Sys.Date(), ".csv", sep = "_")
#         },
#         
#         content = function(file){
#           write.csv(run_fgsea(deseq2_results = load_deseq_output(),
#                                     gmt_gene_sets = load_gene_sets(),
#                                     id2gene_file = load_id2gene_file(),
#                                     min_size = isolate(input$minval),
#                                     max_size = isolate(input$maxval),
#                                     padj_threshold = isolate(input$fgsea_table_range),
#                                     nes_filter = isolate(input$select_nes_pathways)), file, row.names = FALSE)
#           
#           
#       
#         })
      
      
      #' Function to plot top ten positive NES and top ten negative NES pathways
      #' in a barchart
      #'
      #' @param fgsea_results (tibble): the fgsea results in tibble format returned by
      #'   the previous function
      #' @param num_paths (int): the number of pathways for each direction (top or
      #'   down) to include in the plot. Set this at 10.
      #'
      #' @return ggplot with a barchart showing the top twenty pathways ranked by positive
      #' and negative NES
      #' @export
      #'
      #' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
      top_pathways <- function(fgsea_results, num_paths){

        fgsea_results %>%
          mutate(pathway = forcats::fct_reorder(pathway, NES))

        #getting min and max values for NES in the results tibble, merging them

        lowest_NES <- slice_min(fgsea_results, NES, n = num_paths)

        highest_NES <- slice_max(fgsea_results, NES, n= num_paths)

        merged_top_res <- bind_rows(lowest_NES, highest_NES) %>%
          arrange(NES) %>%
          mutate(pathway = reorder(pathway, NES))
        

        #make_plot
        ggplot(merged_top_res) +
          geom_bar(aes(x=pathway, y=NES, fill = as.character(sign(NES))), stat='identity', show.legend = FALSE) +
          scale_fill_manual(values = c('1' = 'red', '-1' = 'blue')) +

          theme_minimal() +
          ggtitle('fgsea results for Hallmark MSigDB gene sets') +
          ylab('Normalized Enrichment Score (NES)') +
          xlab('')+
          theme(axis.text.x = element_text(size=6))+
          theme(axis.text.y = element_text(size = 4))+
          coord_flip()


      }
      
      output$fgsea_top_paths<- renderPlot({
        
        input$fgsea_barplot_do #links action button to pvalue slider
        
        top_pathways(fgsea_results = fgsea_output(), 
                    num_paths = isolate(input$fgsea_barplot_range))
        
      })
      
      
      #' fgsea_scatter_plot
      #'@details Generates a scatterplot of -log10(p-adjusted) vs. NES of 
      #'fgsea results, with values below padjusted threshold in gray 
      
      #'
      #' @param data tibble: a (n x _S_) data set
      #' @param meta tibble: sample-level meta information (_S_ x 3)
      #' @param title string: title for plot
      #'
      #' @return ggplot: scatter plot showing each sample in the first two PCs.
      #'
      #' @examples
      #' `plot_pca(data, meta, "Raw Count PCA")`
      
      fgsea_scatter_plot <- function(fgsea_results, padj_threshold) {

        
        fgsea_gg <- ggplot(fgsea_results, aes(x= NES, y=-log10(padj))) +
          geom_point(aes(color = ifelse(-log10(padj) < padj_threshold, "FALSE", "TRUE")))+
          scale_color_manual(values = c("TRUE" = 'lightblue2', "FALSE" = 'gray'))+
          theme_bw()+
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
          
          labs(x = "NES",
               y = "-log10(padj)",
               title = "FGSEA NES vs. -log10(padj)",
               color = paste("-log10(padj) <= ", padj_threshold))
        
        return(fgsea_gg)
      }
      
      
      output$fgsea_scatter <- renderPlot({
        
        input$fgsea_scatter_do #connects slider input to action button
        
        fgsea_scatter_plot(fgsea_results = fgsea_output(),
                           padj_threshold = isolate(input$fgsea_scatter_range))
        
      })
      
      
      
      
      
      
      
      
        
  
}










# Run the application 
shinyApp(ui = ui, server = server)
