library(shiny)
library(shinyjs)
library(shinyFiles)
source("app_code.R")

root <- c(wd = normalizePath("."), root = "/", home = normalizePath("~"))

# things to add: Show head of .fastq files and barcodes once selected

ui <- fluidPage(
  useShinyjs(),
  
  tags$head(tags$style("
      .Fixed_input_row_75{height:75px;}"
  )),
  
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
               checkboxInput("is.paired", "Is your data paired-end?", value = TRUE),
               actionButton("start_file_input", "Proceed to file input"))
    ),
    mainPanel(
      uiOutput("barcode_locations_image"))
  ),
  
  titlePanel("File Selection: "),
  
  h4("Fastqs:"),
  sidebarLayout(sidebarPanel(fluidRow(class = "Fixed_input_row_75",
                                      column(6, hidden(textOutput("fastq_R1_note"))),
                                      column(6, hidden(shinyFilesButton("fastq_R1", "R1 Fastq(s)", title = "Select R1 Fastq file(s)", multiple = TRUE)))),
                             fluidRow(class = "Fixed_input_row_75",
                                      column(6, hidden(textOutput("fastq_R2_note"))),
                                      column(6, hidden(shinyFilesButton("fastq_R2", "R2 Fastq(s)", title = "Select R2 Fastq file(s)", multiple = TRUE)))),
                             fluidRow(class = "Fixed_input_row_75",
                                      column(6, hidden(textOutput("fastq_R3_note"))),
                                      column(6, hidden(shinyFilesButton("fastq_R3", "R3 Fastq(s)", title = "Select R3 Fastq file(s)", multiple = TRUE))))),
    mainPanel(htmlOutput("fastq_files_report"))
  ),
  
  h4("Barcodes:"),
  sidebarLayout(sidebarPanel(fluidRow(class = "Fixed_input_row_75",
                                      column(6, hidden(textOutput("barcode_file_1_note"))),
                                      column(6, hidden(shinyFilesButton("barcode_file_1", "Barcode 1", title = "Select file containing barcodes", multiple = FALSE)))),
                             fluidRow(class = "Fixed_input_row_75",
                                      column(6, hidden(textOutput("barcode_file_2_note"))),
                                      column(6, hidden(shinyFilesButton("barcode_file_2", "Barcode 2", title = "Select file containing barcodes", multiple = FALSE))))),
                mainPanel(htmlOutput("barcode_files_report"))
  ),
  
  
  textOutput("file_name_header"),
  uiOutput("demultiplex_options"),
  actionButton("go_demultiplex", "Demultiplex reads!"),
  
  
  actionButton("input_pre_demultiplexed", "I have pre-demultiplexed fastq files!"),
  hidden(uiOutput("pre_demulti_input_panel")),
  #===============Alignment============
  titlePanel("Alignment"),
  mainPanel(tabsetPanel(type = "tabs",
                        tabPanel("Alignment Options", 
                                 hidden(radioButtons("alignment_radio", 
                                                     label = "Alignment Method", 
                                                     choices = c("reference", "denovo"), 
                                                     selected = character())),
                                 hidden(uiOutput("alignment_panel"))),
                        tabPanel("Demultiplexed File Report", 
                                 fluidRow(column(6, 
                                                 uiOutput("demultiplexed_files_report_RA")),
                                          
                                          column(6, 
                                                 uiOutput("demultiplexed_files_report_RB")))))
  )
  
)





server <- function(input, output, session) {
  
  #==========read in input files for demultiplexing===========
  x <- reactiveValues(files = list(), 
                      files_report = list(), 
                      files_are_good = 0, 
                      have_sample_ids = 0,
                      demultiplexed = 0)
  
  # figure out which images we need
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
  
  # generate images
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
  shinyFileChoose(input, 'fastq_R1', root=c(root))
  shinyFileChoose(input, 'fastq_R2', root=c(root))
  shinyFileChoose(input, 'fastq_R3', root=c(root))
  
  
  ## barcodes
  shinyFileChoose(input, 'barcode_file_1', root=c(root))
  shinyFileChoose(input, 'barcode_file_2', root=c(root))
  
  # reset if the input type changes!
  observeEvent(input$barcodes, ignoreInit = TRUE, {
    x$files_are_good <- 0
  })
  observeEvent(input$is.paired, ignoreInit = TRUE, {
    x$files_are_good <- 0
  })
  
  # show or hide the appropriate input selectors and edit text labels
  ## define consistant notes
  output$fastq_R1_note <- renderText("Forward read files (R1):")
  output$fastq_R3_note <- renderText("Reverse read files (R3):")
  output$barcode_file_2_note <- renderText("Barcodes matching R2 (seperate fastq file barcodes):")
  
  ## define others and show appropriate
  observeEvent(eventExpr = input$start_file_input, handlerExpr = {
    hide("fastq_R1")
    hide("fastq_R2")
    hide("fastq_R3")
    hide("barcode_file_1")
    hide("barcode_file_2")
    hide("fastq_R1_note")
    hide("fastq_R2_note")
    hide("fastq_R3_note")
    hide("barcode_file_1_note")
    hide("barcode_file_2_note")
    
   if(!all(c("Start of Reads", "In Headers") %in% input$barcodes)){
     if(length(input$barcodes) != 0){
       show("fastq_R1")
       show("fastq_R1_note")
     }
     
     
     if("In Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes){
       if("In Headers" %in% input$barcodes){
         output$barcode_file_1_note <- renderText("Header barcodes:")
       }
       else{
         output$barcode_file_1_note <- renderText("In-read barcodes:")
       }
       show("barcode_file_1_note")
       show("barcode_file_1")
     }
     
     if(input$is.paired){
       if(length(input$barcodes) != 0){
         show("fastq_R2")
         show("fastq_R2_note")
       }
       
       if("Seperate fastq File" %in% input$barcodes){
         output$fastq_R2_note <- renderText("Barcode reads (R2):")
         
         show("fastq_R3_note")
         show("barcode_file_2_note")
         show("fastq_R3")
         show("barcode_file_2")
       }
       else{
         output$fastq_R2_note <- renderText("Reverse read files (R2):")
       }
     }
     else if("Seperate fastq File" %in% input$barcodes){
       output$barcode_file_2_note <- renderText("R2 barcodes:")
       output$fastq_R2_note <- renderText("Barcode reads (R2):")
       
       show("barcode_file_2_note")
       show("fastq_R2_note")
       show("fastq_R2")
       show("barcode_file_2")
     }
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
  
  # take sample IDs
  output$demultiplex_options <- renderUI({
    sib <- shinyFilesButton("sample_ids", label = "Optional Sample Identifications", multiple = FALSE, 
                            title = "Tab delimited key for demultiplexed samples. Two or three columns for: seperate fastq file barcodes, header/in-read barcodes, and sample name.")
    hidden(sib)
  })
  
  output$file_name_header <- renderText("Demultiplexing Options:")
  
  # generate reports on barcode and fastq files
  output$barcode_files_report <- renderUI({
    lapply(x$files_report$barcode, p)
  })
  output$fastq_files_report <- renderUI({
    lapply(x$files_report$fastq, p)
  })
  
  # show or hide reports and 'go' buttons based on inputs
  observeEvent(x$files_are_good, {
    if(x$files_are_good){
      x$files_report <- .fastq_file_reporter(x$files)
      
      
      show("barcode_files_report")
      show("fastq_files_report")
      show("sample_ids")
      show("go_demultiplex", asis = TRUE)
      show("file_name_header", asis = TRUE)
    }
    else{
      hide("barcode_files_report")
      hide("fastq_files_report")
      hide("go_demultiplex")
      hide("file_name_header")
      hide("sample_ids")
    }
  })
  shinyFileChoose(input, 'sample_ids', root=c(root))

  # import and note that we have sample ids if we do
  observeEvent(input$sample_ids,{
    x$files$sample_ids <- .parse_shinyFiles_path(root, input$sample_ids)
    x$have_sample_ids <- 1
  })
  
  #==========run demultiplexing======================
  # demultiplex
  observeEvent(input$go_demultiplex,{
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
  
  #==========import existing demultiplexed reads==============
  # panel layout
  output$pre_demulti_input_panel <- renderUI({
      inputPanel(
        column(12,
               checkboxInput("is.paired.pre.existing", "Is your data paired-end?", value = TRUE),
               shinyFilesButton("pre_existing_RA", label = "RA Files", multiple = TRUE, title = "Select RA Reads"),
               hidden(shinyFilesButton("pre_existing_RB", label = "RB Files", multiple = TRUE, title = "Select RB Reads")),
               hidden(actionButton("parse_pre_existing_demultiplexed", "Read In Demultiplexed Files"))))
  })
  
  # show panel only if pre_demultiplexed file input is selected
  observeEvent(input$input_pre_demultiplexed,{
    toggle("pre_demulti_input_panel")
  })
  
  observeEvent(input$is.paired.pre.existing,{
    toggle("pre_existing_RB")
  })
  
  # show 'go' button if correct files have been selected given paired/unpaired
  observeEvent(is.list(input$pre_existing_RA) | is.list(input$pre_existing_RB), ignoreInit = TRUE, {
    if(input$is.paired.pre.existing){
      if(is.list(input$pre_existing_RA) & is.list(input$pre_existing_RB)){
        show("parse_pre_existing_demultiplexed")
      }
      else{
        hide("parse_pre_existing_demultiplexed")
      }
    }
    else{
      if(is.list(input$pre_existing_RA)){
        show("parse_pre_existing_demultiplexed")
      }
      else{
        hide("parse_pre_existing_demultiplexed")
      }
    }
  })
  
  shinyFileChoose(input, 'pre_existing_RA', root=c(root))
  shinyFileChoose(input, 'pre_existing_RB', root=c(root))
  
  # parse in files
  observeEvent(eventExpr = input$parse_pre_existing_demultiplexed,
               ignoreInit = TRUE, handlerExpr = {
                 x$files$demultiplexed_files <- NULL # clear out any others that have been previously entered.
                 
                 x$demultiplexed <- 0
                 if(input$is.paired.pre.existing){
                   if(is.list(input$pre_existing_RA) & is.list(input$pre_existing_RB)){
                     RA <- .parse_shinyFiles_path(root, input$pre_existing_RA)
                     RB <- .parse_shinyFiles_path(root, input$pre_existing_RB)
                     
                     validate(need(length(RA) == length(RB), message = "An equal number of RA and RB files must be provided."))
                     
                     x$files$demultiplexed_files <- vector("list", 2)
                     names(x$files$demultiplexed_files) <- c("RA", "RB")
                     x$files$demultiplexed_files$RA <- RA
                     x$files$demultiplexed_files$RB <- RB
                     
                     x$demultiplexed <- 1
                   }
                 }
                 else{
                   if(is.list(input$pre_existing_RA)){
                     x$files$demultiplexed_files <- list(RA = .parse_shinyFiles_path(root, input$pre_existing_RA))
                     
                     x$demultiplexed <- 1
                   }
                 }
                 
                 # reset our file inputs for validation purposes
                 reset("pre_existing_RB") 
                 reset("pre_existing_RA")
               })
  
  
  #==========report demultiplex results=======================
  # generate reporters
  output$demultiplexed_files_report_RA <- renderUI({
    lapply(x$files_report$demulti$RA, p)
  })
  output$demultiplexed_files_report_RB <- renderUI({
    lapply(x$files_report$demulti$RB, p)
  })
  
  # show or hide reporters
  observeEvent(x$demultiplexed, ignoreInit = TRUE, {
    if(x$demultiplexed){
      x$files_report$demulti <- .fastq_file_reporter(x$files$demultiplexed_files)$demulti
      
      
      show("demultiplexed_files_report_RA")
      show("demultiplexed_files_report_RB")
      
    }
    else{
      hide("demultiplexed_files_report_RA")
      hide("demultiplexed_files_report_RB")
    }
  })
  
  #==========run alignment============================
  # generate the alignment panel based on denovo or reference, paired or unpaired inputs
  output$alignment_panel <- renderUI({
    
    validate(need(input$alignment_radio %in% c("denovo", "reference"), message = "Please select an alignment method."))
    if(input$alignment_radio == "denovo"){
      if(length(x$demultiplexed_files) == 2){
        fluidRow(inputPanel(
          numericInput("M", "M, the number of mismatches allowed between alignments within individuals", 3, 1, step = 1),
          numericInput("n", "n, the number of mismatches allowed between alignments between individuals", 3, 1, step = 1),
          checkboxInput("check_headers", "Check that fastq headers are STACKS acceptable. If you are sure they are OK, uncheck this.", value = TRUE),
          checkboxInput("stacks_cleanup", "Clean-up accessory files after completion.", value = TRUE),
          numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
          checkboxInput("rmdup", "Remove PCR duplicates?", value = TRUE),
          checkboxInput("rmimp", "Remove improperly paired reads? This should usually be FALSE.", value = FALSE),
          numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
        ))
      }
      else{
        fluidRow(inputPanel(
          numericInput("M", "M, the number of mismatches allowed between alignments within individuals", 3, 1, step = 1),
          numericInput("n", "n, the number of mismatches allowed between alignments between individuals", 3, 1, step = 1),
          checkboxInput("check_headers", "Check that fastq headers are STACKS acceptable. If you are sure they are OK, uncheck this.", value = TRUE),
          checkboxInput("stacks_cleanup", "Clean-up accessory files after completion.", value = TRUE),
          numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
          numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
        ))
      }
    }
    else if(input$alignment_radio == "reference"){
      if(length(x$demultiplexed_files) == 2){
        fluidRow(inputPanel(
          shinyFilesButton("reference_genome", "Select Reference Genome", title = "Select Reference Genome", multiple = FALSE),
          numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
          checkboxInput("rmdup", "Remove PCR duplicates?", value = TRUE),
          checkboxInput("rmimp", "Remove improperly paired reads?", value = TRUE),
          numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
        ))
      }
      else{
        fluidRow(inputPanel(
          shinyFilesButton("reference_genome", "Select Reference Genome", title = "Select Reference Genome", multiple = FALSE),
          numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
          numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
        ))
      }
    }
  })
  
  shinyFileChoose(input, 'reference_genome', root=c(root))
  
  observeEvent(x$demultiplexed, ignoreInit = TRUE, {
    show("alignment_radio")
  })
  
  observeEvent(input$alignment_radio, ignoreInit = TRUE,{
    show("alignment_panel")
  })
} 

shinyApp(ui, server)
