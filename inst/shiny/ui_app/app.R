library(shiny)
library(shinyjs)
library(shinyFiles)
source("app_code.R")

root <- c(wd = normalizePath("."), root = "/", home = normalizePath("~"))

# things to add: Show head of .fastq files and barcodes once selected

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
               actionButton("start_file_input", "Proceed to file input"))
    ),
    mainPanel(
      uiOutput("barcode_locations_image"))
  ),
  
  titlePanel("File Selection: "),
    
  sidebarLayout(sidebarPanel(fluidRow(h4("Fastqs:"), 
                                      uiOutput("input_file_selector"))),
                mainPanel(htmlOutput("fastq_files_report"))
  ),
  
  sidebarLayout(sidebarPanel(fluidRow(h4("Barcodes:"),
                                     uiOutput("input_barcode_selector"))),
               mainPanel(htmlOutput("barcode_files_report")),
  ),
  textOutput("file_name_header"),
  uiOutput("demultiplex_options"),
  actionButton("go_demultiplex", "Demultiplex reads!")
  
  #===============Alignment============
  titlePanel("Alignment"),
  
  inputPanel(sidebarLayout(
    fluidRow(h4("If you already have demultiplexed data, select RA and RB (forward and reverse) reads below."),
             shinyFilesButton("pre_existing_RA", label = "RA Files", multiple = TRUE, title = "Select RA Reads"),
             shinyFilesButton("pre_existing_RB", label = "RB Files", multiple = TRUE, title = "Select RB Reads"))
    
  ))
  
)





server <- function(input, output, session) {
  
  #==========read in input files for demultiplexing===========
  x <- reactiveValues(files = list(), 
                      files_report = list(), 
                      files_are_good = 0, 
                      have_sample_ids = 0,
                      demultiplexed_files = list(),
                      demultiplexed = 0)
  
  # generate image
  barcode.type <- eventReactive(input$barcodes,{
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
      shinyFilesButton(paste0("fastq_R", i), label = file_R_lab, multiple = TRUE, accept = c(".fastq", ".fq"), title = ".fastq or .fq file.")
    })
    
    shinyjs::hidden(fileinputs)
  })
  shinyFileChoose(input, 'fastq_R1', root=c(root))
  shinyFileChoose(input, 'fastq_R2', root=c(root))
  shinyFileChoose(input, 'fastq_R3', root=c(root))
  
  
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
  shinyFileChoose(input, 'barcode_file_1', root=c(root))
  shinyFileChoose(input, 'barcode_file_2', root=c(root))
  
  # reset if the input type changes!
  observeEvent(input$barcodes, ignoreInit = TRUE, {
    x$files_are_good <- 0
  })
  observeEvent(input$is.paired, ignoreInit = TRUE, {
    x$files_are_good <- 0
  })
  
  observeEvent(eventExpr = input$start_file_input, handlerExpr = {
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

      # not paired, dual indexed
      if(!input$is.paired & "Seperate fastq File" %in% input$barcodes & is.list(input$fastq_R2) & is.list(input$barcode_file_2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$barcode_file_2 <- .parse_shinyFiles_path(root, input$barcode_file_2)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R2 <- .parse_shinyFiles_path(root, input$fastq_R2)
        x$files_are_good <- 1
      }
      
      # paired, dual indexed
      else if(input$is.paired & "Seperate fastq File" %in% input$barcodes & is.list(input$fastq_R2) & is.list(input$fastq_R3) & is.list(input$barcode_file_2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$barcode_file_2 <- .parse_shinyFiles_path(root, input$barcode_file_2)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R2 <- .parse_shinyFiles_path(root, input$fastq_R2)
        files$R3 <- .parse_shinyFiles_path(root, input$fastq_R3)
        x$files_are_good <- 1
      }
      
      # not paired, single indexed
      else if(!input$is.paired & !"Seperate fastq File" %in% input$barcodes){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        x$files_are_good <- 1
      }
      
      else if(input$is.paired & !"Seperate fastq File" %in% input$barcodes & is.list(input$fastq_R2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R3 <- .parse_shinyFiles_path(root, input$fastq_R2) # named R3 for consistancy later
        x$files_are_good <- 1
      }
    }
    
    ## with only R2 barcodes
    else if(!("In Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes) & is.list(input$barcode_file_1) & is.list(input$fastq_R1)){
      # not paired
      if(!input$is.paired){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R2 <- .parse_shinyFiles_path(root, input$fastq_R2)
        x$files_are_good <- 1
      }
      
      # paired
      if(input$is.paired & is.list(input$fastq_R2)){
        files$barcode_file_1 <- .parse_shinyFiles_path(root, input$barcode_file_1)
        files$R1 <- .parse_shinyFiles_path(root, input$fastq_R1)
        files$R3 <- .parse_shinyFiles_path(root, input$fastq_R2)
        x$files_are_good <- 1
      }
    }
    
    x$files <- files
  })
  
  # accept sample IDs and begin
  output$demultiplex_options <- renderUI({
    sib <- shinyFilesButton("sample_ids", label = "Optional Sample Identifications", multiple = FALSE, 
                            title = "Tab delimited key for demultiplexed samples. Two or three columns for: seperate fastq file barcodes, header/in-read barcodes, and sample name.")
    shinyjs::hidden(sib)
  })
  
  output$file_name_header <- renderText("Demultiplexing Options:")
    output$barcode_files_report <- renderUI({
    lapply(x$files_report$barcode, p)
  })
  
  output$fastq_files_report <- renderUI({
    lapply(x$files_report$fastq, p)
  })
  
  observeEvent(x$files_are_good, {
    if(x$files_are_good){
      x$files_report <- .fastq_file_reporter(x$files)
      
      
      shinyjs::show("barcode_files_report")
      shinyjs::show("fastq_files_report")
      shinyjs::show("sample_ids")
      shinyjs::show("go_demultiplex", asis = TRUE)
      shinyjs::show("file_name_header", asis = TRUE)
    }
    else{
      shinyjs::hide("barcode_files_report")
      shinyjs::hide("fastq_files_report")
      shinyjs::hide("go_demultiplex")
      shinyjs::hide("file_name_header")
      shinyjs::hide("sample_ids")
    }
  })
  shinyFileChoose(input, 'sample_ids', root=c(root))

  observeEvent(input$sample_ids,{
    x$files$sample_ids <- .parse_shinyFiles_path(root, input$sample_ids)
    x$have_sample_ids <- 1
  })
  
  #==========run demultiplexing======================
  observeEvent(input$go_demultiplex,{
    browser()
    if(x$files_are_good){
      if("Seperate fastq File" %in% input$barcodes){
        
        if(length(input$barcodes) > 1){
          indices <- readLines(x$files$barcode_file_2)
        }
        else{
          indices <- readLines(x$barcode_file_1)
        }
        
        if(input$is.paired){
          R2 <- x$files$R2
        }
        else{
          R2 <- NULL
        }
        
        x$files$plate_split_files <- plate_split(R1 = x$files$R1, 
                                                 R2 = R2,
                                                 R3 = x$files$R3, 
                                                 indices = indices,
                                                 outfile_prefix = "alignR_gui_plate_split")
        if(length(input$barcodes) == 1){
          x$files$demultiplexed_files <- x$files$plate_split_files
          x$files$plate_split_files <- NULL
        }
      }
      
      if(any(c("In Headers", "Start of Reads") %in% input$barcodes)){
        # fill R1
        if("Seperate fastq File" %in% input$barcodes){
          R1 <- x$files$plate_split_files$R1
        }
        else{
          R1 <- x$files$R1
        }
        
        # fill R2
        if(input$is.paired){
          if("Seperate fastq File" %in% input$barcodes){
            R2 <- x$files$plate_split_files$R3
          }
          else{
            R2 <- x$files$R3
          }
        }
        else{
          R2 <- NULL
        }
        
        # fill sample_names
        if(x$have_sample_ids){
          sample_names <- readLines(x$files$sample_ids)
        }
        else{
          sample_names <- NULL
        }
        
        
        x$files$demultiplexed_files <- demultiplex(R1 = R1,
                                                   R2 = R2,
                                                   sample_names = sample_names,
                                                   barcodes = readLines(x$files$barcode_file_1),
                                                   outfile_prefix = "alignR_gui_demultiplex")
      }
      
      x$demultiplexed <- 1
    }
  })
  
  #==========figure out parameters for aligning==============
} 

shinyApp(ui, server)