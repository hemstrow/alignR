library(shiny)
library(shinyjs)
library(shinyFiles)
source("app_code.R")

root <- c(wd = normalizePath("."), root = "/")

ui <- fluidPage(
  useShinyjs(),
  
  #===============Demultiplexing and initial inputs===============
  titlePanel("Demultiplexing:"),
  sidebarLayout( 
    sidebarPanel(
      fluidRow(selectizeInput("barcodes", 
                              "Where are your barcodes located?", 
                              choices = list("In Headers",
                                             "Start of Reads",
                                             "Seperate fastq File"), 
                              options = list(maxItems = 2),
                           multiple = TRUE),
               checkboxInput("is.paired", "Is your data paired-end?", value = FALSE),
               actionButton("visualize_barcodes", "Visualize and proceed to file input"))
    ),
    mainPanel(
      uiOutput("barcode_locations_image"))
  ),
  
  titlePanel("File Selection: "),
  inputPanel(
    fluidRow(column(5, h4("Fastqs:"), uiOutput("input_file_selector")),
             column(5, h4("Barcodes:"), uiOutput("input_barcode_selector"))
    ),
    h4("Specify sample names?"),
    uiOutput("demultiplex_options"),
    uiOutput("go_demultiplex")
  ),
)





server <- function(input, output, session) {
  
  #==========read in input files===========
  x <- reactiveValues()
  
  # generate image
  barcode.type <- eventReactive(input$visualize_barcodes,{
    img_table <- data.frame(input = c("In Headers",
                                      "Start of Reads",
                                      "Seperate fastq File"),
                            img = c("barcode_in_header", 
                                    "barcode_in_read", 
                                    "barcode_in_r2"))
    
    img <- img_table$img[which(img_table$input %in% input$barcodes)]
    return(img)
  })
  
  
  output$barcode_locations_image <- renderUI({
    res <- list()
    validate(need(!all(c("Start of Reads", "In Headers") %in% input$barcodes), "Barcodes can be either in fastq headers OR at the beginning of reads, but not both!"))
    photo_size <- 300
    if(length(barcode.type()) > 1){
      photo_size <- photo_size - 50*(length(barcode.type()) - 1)
    }
    
    photo_size <- floor(photo_size/length(barcode.type()))
    
    photo_size <- paste0(photo_size, "px")
    
    if('barcode_in_header' %in% barcode.type()){
      res[[length(res) + 1]] <- img(src='barcode_in_header.JPG', height = photo_size)
    }
    if('barcode_in_read' %in% barcode.type()){
      res[[length(res) + 1]] <- img(src='barcode_in_read.JPG', height = photo_size)
    }
    if('barcode_in_r2' %in% barcode.type()){
      res[[length(res) + 1]] <- img(src='barcode_in_r2.JPG', height = photo_size)
    }
    
    if(length(res) > 1){
      new_res <- vector("list", length(res) - 1)
      prog <- 1
      for(i in 1:length(res)){
        new_res[[prog]] <- res[[i]]
        prog <- prog + 1
        if(i != length(res)){
          new_res[[prog]] <- img(src='plus.JPG', height = "50px")
          prog <- prog + 1
        }
      }
      return(new_res)
    }
    else{return(res)}
    
  })
  
  
  
  # generate selector boxes for both the input fastqs and the barcodes
  ## fastqs
  output$input_file_selector <- renderUI({
    validate(need(!all(c("Start of Reads", "In Headers") %in% input$barcodes), ""))
    
    file_R_lab <- "NA"
    fileinputs <- lapply(1:3, function(i){
      if(i == 1){
        file_R_lab <- "R1: read file"
      }
      else if(i == 2){
        if("Seperate fastq File" %in% input$barcodes){
          file_R_lab <- "R2: barcode file"
        }
        else{
          file_R_lab <- "R2: read file"
        }
      }
      else{
        file_R_lab <- "R3: read file"
      }
      shinyFilesButton(paste0("fastq_R", i), label = file_R_lab, multiple = FALSE, accept = c(".fastq", ".fq"), title = ".fastq or .fq file.")
    })
    
    shinyjs::hidden(fileinputs)
  })
  shinyFileChoose(input, 'fastq_R1', root=c(root=root))
  shinyFileChoose(input, 'fastq_R2', root=c(root=root))
  shinyFileChoose(input, 'fastq_R3', root=c(root=root))
  
  
  ## barcodes
  output$input_barcode_selector <- renderUI({
    validate(need(!all(c("Start of Reads", "In Headers") %in% input$barcodes), ""))
    
    barcode_file_lab_index <- "NA"
    barcodeinputs <- lapply(1:2, function(i){
      if(i == 1){
        if("In Headers" %in% input$barcodes){
          barcode_file_lab_index <- "header barcodes"
        }
        else if("Start of Reads" %in% input$barcodes){
          barcode_file_lab_index <- "in-read barcodes"
        }
        else{
          barcode_file_lab_index <- "R2 read barcodes"
        }
      }
      else{
        barcode_file_lab_index <- "R2 read barcodes"
      }
      shinyFilesButton(paste0("barcode_file_", i), label = barcode_file_lab_index, multiple = FALSE, title = "File with barcodes/indices, one per line.")
    })
    
    shinyjs::hidden(barcodeinputs)
  })
  shinyFileChoose(input, 'barcode_file_1', root=c(root=root))
  shinyFileChoose(input, 'barcode_file_2', root=c(root=root))
  
  
  
  observeEvent(eventExpr = input$visualize_barcodes, handlerExpr = {
    if(length(input$barcodes) != 0){
      shinyjs::show("fastq_R1")
    }
    
    
    if("In Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes){
      shinyjs::show("barcode_file_1")
    }
    
    if(input$is.paired){
      shinyjs::show("fastq_R2")
      
      if("Seperate fastq File" %in% input$barcodes){
        shinyjs::show("fastq_R3")
        shinyjs::show("barcode_file_2")
      }
    }
    else if("Seperate fastq File" %in% input$barcodes){
      shinyjs::show("fastq_R2")
      shinyjs::show("barcode_file_2")
    }
  })
  

  output$demultiplex_options <- renderUI({
    sids <- fileInput("sample_ids", "Optional Sample Identifications", placeholder = "Tab delimited key for demultiplexed samples. Two or three columns for: seperate fastq file barcodes, header/in-read barcodes, and sample name.")
    shinyjs::hidden(sids)
  })
  output$go_demultiplex <- renderUI({
    gdb <- actionButton("go_demultiplex", "Demultiplex reads!")
    shinyjs::hidden(gdb)
  })

  
  # locate appropriate files
  observeEvent(eventExpr = 
                 is.list(input$fastq_R1) | 
                 is.list(input$fastq_R2) | 
                 is.list(input$fastq_R3) | 
                 is.list(input$barcode_file_1) | 
                 is.list(input$barcode_file_2), ignoreInit = TRUE, handlerExpr = {
    display <- FALSE
    files <- list()

    # read files for every case:
    ## with in-line or header barcodes
    if(("In Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes) & is.list(input$barcode_file_1) & is.list(input$fastq_R1)){
      browser()
      
      # not paired, dual indexed
      if(!input$is.paired & "Seperate fastq File" %in% input$barcodes & is.list(input$fastq_R2) & is.list(input$barcode_file_2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$barcode_file_2 <- .parse_shinyFiles_path(root, input$barcode_file_2)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R2 <- .parse_shinyFiles_path(root, input$fastq_R2)
      }
      
      # paired, dual indexed
      else if(input$is.paired & "Seperate fastq File" %in% input$barcodes & is.list(input$fastq_R2) & is.list(input$fastq_R3) & is.list(input$barcode_file_2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$barcode_file_2 <- .parse_shinyFiles_path(root, input$barcode_file_2)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R2 <- .parse_shinyFiles_path(root, input$fastq_R2)
        files$R3 <- .parse_shinyFiles_path(root, input$fastq_R3)
      }
      
      # not paired, single indexed
      else if(!input$is.paired & !"Seperate fastq File" %in% input$barcodes){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
      }
      
      else if(input$is.paired & !"Seperate fastq File" %in% input$barcodes & is.list(input$fastq_R2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R3 <- .parse_shinyFiles_path(root, input$fastq_R2) # named R3 for consistancy later
      }
    }
    
    ## with only R2 barcodes
    else if(!("In Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes) & is.list(input$barcode_file_1) & is.list(input$fastq_R1)){
      # not paired
      browser()
      if(!input$is.paired){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R2 <- .parse_shinyFiles_path(root, input$fastq_R2)
      }
      
      # paired
      if(input$is.paired & is.list(input$fastq_R2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R3 <- .parse_shinyFiles_path(root, input$fastq_R2)
      }
    }
    
  })
} 

shinyApp(ui, server)