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
library(reshape) #used for "melting" counts data to make heatmap
library(ggplotify) ## to convert pheatmap to ggplot2

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
                                       min = 0, max = 1,
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
                                    
                                    tabPanel("Volcano Plot of DE Results", plotOutput("volc_plot")))
                                    
                                    
                                    
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
                           fileInput("upload", paste0("Load FGSEA results"), accept = c(".csv", ".tsv")),
                           
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
  

#  filtering step, convert to matrix, do log10 transformation

  filtered_counts <- count_info_data %>%

    filter(
      rowSums(count_info_data != 0) > count_filter &
        apply(count_info_data, 1, function(row) quantile (row, prob = var_filter/100) > 0)

    ) 

  filtered_counts_matrix <- as.matrix(filtered_counts)
  
  # log10_filtered_matrix <- log10(filtered_counts_matrix)


  
  # Subset the matrix to include all columns and the first 30 rows
  subset_counts <- filtered_counts_matrix[1:30, ]
  
  
  #make heatmap 
  
  pheatmap(filtered_counts_matrix, scale = "row",
           color=colorRampPalette(c("white", "red"))(50))
  
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
  

    #adding additional up/down/ns column    
    annotated_deseq_res <- deseq_mat %>%
      mutate(volc_plot_status = case_when(padj < padj_threshold & log2FoldChange > 0 ~ 'UP',
                                          padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN',
                                          TRUE ~ 'NS')) 

    
    datatable(annotated_deseq_res, rownames = TRUE, options = list(pageLength = nrow(annotated_deseq_res)))
  }
  
  #sample data table
  output$deseq_tab <- renderDT({
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
        
        print(volc_data)

          gg4 <- ggplot(volc_data, aes(x=log2FoldChange,y=`-log10(adjusted p)`, color = volc_plot_status)) +
            geom_point()+
            theme_bw()+
            geom_hline(yintercept = 0, linetype = "dashed") +
            labs(x = 'log2FoldChange', y = '-log10(padj)', title = 'Volcano plot of DESeq2 differential expression results (HD vs. Control)')

          gg4


      }

      output$volc_plot<- renderPlot({
        
        plot_volcano(deseq2_results = load_deseq_output(),
                     padj_threshold = isolate(input$padj_range))
        
      })

      
  
}

#' 
#'   #' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#'   #' adds a column to denote plotting status in volcano plot. Column should denote
#'   #' whether gene is either 1. Significant at padj < .10 and has a positive log
#'   #' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#'   #' 3. Not significant at padj < .10. Have the values for these labels be UP,
#'   #' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'   #'
#'   #' @param deseq2_res (df): results from DESeq2 
#'   #' @param padj_threshold (float): threshold for considering significance (padj)
#'   #'
#'   #' @return Tibble with all columns from DESeq2 results and one additional column
#'   #'   labeling genes by significant and up-regulated, significant and
#'   #'   downregulated, and not significant at padj < .10.
#'   #'   
#'   #' @export
#'   #'
#'   #' @examples labeled_results <- label_res(res, .10)
#'   label_res <- function(deseq2_res, padj_threshold) {
#'     
#'     #add deseq2 rownames to new column called genes so they get preserved in tibble conversion
#'     deseq2_res$genes <- rownames(deseq2_res)
#'     
#'     #create dds results tibble
#'     res <- as_tibble(deseq2_res)
#'     
#'     
#'     #delete pthresh when done
#'     
#'     p0_vs_Ad_volc_plot <- res %>%
#'       mutate(volc_plot_status = case_when(padj < padj_threshold & log2FoldChange > 0 ~ 'UP', 
#'                                           padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN',
#'                                           TRUE ~ 'NS')) %>%
#'       relocate(genes) %>%
#'       return(p0_vs_Ad_volc_plot)
#'     
#'   }
#' 
#'   #' Function to plot the unadjusted p-values as a histogram
#'   #'
#'   #' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#'   #' column denoting status in volcano plot
#'   #'
#'   #' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#'   #' @export
#'   #'
#'   #' @examples pval_plot <- plot_pvals(labeled_results)
#'   plot_pvals <- function(labeled_results) {
#'     
#'     gg <- ggplot(labeled_results, aes(x = pvalue)) +
#'       geom_histogram(bins = 60, color='black', fill='lightblue2')+
#'       theme_bw()+
#'       labs(x = 'pvalue', y = 'count', title = 'Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)')
#'     gg
#'     
#'     
#'   }
#'   
#'   
#'   
#'   
#'   #' Function to plot the log2foldchange from DESeq2 results in a histogram
#'   #'
#'   #' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#'   #' column denoting status in volcano plot
#'   #' @param padj_threshold (float): threshold for considering significance (padj)
#'   #'
#'   #' @return ggplot: a histogram of log2FC values from genes significant at padj 
#'   #' threshold of 0.1
#'   #' @export
#'   #'
#'   #' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
#'   plot_log2fc <- function(labeled_results, padj_threshold) {
#'     
#'     #filter input to only include results within defined padj_threshold
#'     
#'     filt_results <- labeled_results %>%
#'       filter(padj < padj_threshold) 
#'     
#'     #plotting filtered results
#'     gg2 <- ggplot(filt_results, aes(x = log2FoldChange)) +
#'       geom_histogram(bins = 100, color='black', fill='lightblue2')+
#'       theme_bw()+
#'       labs(x = 'log2FoldChange', y = 'count', title = 'Histogram of Log2FoldChanges for DE Genes (vP0 vs. vAd)')
#'     gg2
#'     
#'   }
#'   
#'   
#'   
#'   
#'   #' Function to make scatter plot of normalized counts for top ten genes ranked
#' #' by ascending padj
#' #'
#' #' @param labeled_results (tibble): Tibble with DESeq2 results and one
#' #'   additional column denoting status in volcano plot
#' #' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' #' the updated DESeqDataSet object with test results
#' #' @param num_genes (int): Number of genes to plot
#' #'
#' #' @return ggplot: a scatter plot with the normalized counts for each sample for
#' #' each of the top ten genes ranked by ascending padj
#' #' @export
#' #'
#' #' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
#' scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
#' 
#'   #rank labeled_results based on padj (ascending order)
#'   
#'   res_sort <- labeled_results %>% 
#'     arrange(padj)
#'   
#'   #select top genes, using num_genes param
#'   topgenes <- head(res_sort, num_genes)
#'   
#'   #get the normalized counts for top genes
#'   
#'   genes <- topgenes$genes
#'   
#'   top_norm_counts <- counts(dds_obj)[genes,]
#'   
#'   #make dataframe to be used for plot with melt fcn (reshape2 library)
#'   
#'   melted_df <- melt(top_norm_counts, varnames = c('Gene', 'Sample'))
#'   
#'   #plot the data
#'   
#'   gg3 <- ggplot(melted_df, aes(x=Gene, y=log10(value), color = Sample))+
#'     geom_jitter(width = 0.2)+
#'     theme_bw()+
#'     labs(x= 'Genes', y = 'log10(norm_counts)', title = 'Plot of Log10(normalized counts) for top ten DE genes')+
#'     theme(axis.text.x = element_text(angle =90))
#'   gg3
#'   
#'   
#' }
#' 
#' #' Function to generate volcano plot from DESeq2 results
#' #'
#' #' @param labeled_results (tibble): Tibble with DESeq2 results and one
#' #'   additional column denoting status in volcano plot
#' #'
#' #' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#' #'   -log10(padj) and labeled by status
#' #' @export
#' #'
#' #' @examples volcano_plot <- plot_volcano(labeled_results)
#' #' 
#' plot_volcano <- function(labeled_results) {
#' 
#'   volc_data <- labeled_results %>%
#'     
#'     #add -log10(padj) column to input tibble
#'     mutate(`-log10(adjusted p)`=-log10(padj),)
#'     
#'     gg4 <- ggplot(volc_data, aes(x=log2FoldChange,y=`-log10(adjusted p)`,color=`volc_plot_status`)) +
#'       geom_point()+
#'       theme_bw()+
#'       geom_hline(yintercept = 0, linetype = "dashed") +
#'       labs(x = 'log2FoldChange', y = '-log10(padj)', title = 'Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)')
#'     
#'     gg4
#' 
#'   
#' }

  

# Run the application 
shinyApp(ui = ui, server = server)
