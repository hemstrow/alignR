library(shiny)
library(shinyjs)
library(shinyFiles)
source("app_code.R")

root <- c(wd = normalizePath("."), root = "/", home = normalizePath("~"))

# things to add: Show head of .fastq files and barcodes once selected

ui <- navbarPage(
  useShinyjs(),
  
  tags$head(tags$style("
      .Fixed_input_row_75{height:75px;}"
  )),
  
  #===============Demultiplexing and initial inputs===============
  tabPanel(title = "Demultiplexing",
           titlePanel("Demultiplexing"),
           actionButton("input_pre_demultiplexed", "I have demultiplexed fastq files!"),
           
           # barcode selection and input panels
           tags$div(id = "demultiplexing",
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
                    
                    
                    hidden(shinyFilesButton("sample_ids", label = "Optional Sample Identifications", multiple = FALSE, 
                                            title = "Tab delimited key for demultiplexed samples. Two or three columns for: seperate fastq file barcodes, header/in-read barcodes, and sample name.")
                    ),
                    actionButton("go_demultiplex", "Demultiplex reads!")
           ),
           
           
           # pre-demultiplexed input
           hidden(tags$div(id = "pre_demulti_input_panel",
                           inputPanel(
                             column(12,
                                    checkboxInput("is.paired.pre.existing", "Is your data paired-end?", value = TRUE),
                                    hidden(shinyFilesButton("pre_existing_RA", label = "RA Files", multiple = TRUE, title = "Select RA Reads")),
                                    hidden(shinyFilesButton("pre_existing_RB", label = "RB Files", multiple = TRUE, title = "Select RB Reads")),
                                    hidden(actionButton("parse_pre_existing_demultiplexed", "Read In Demultiplexed Files"))))))
  ),
  
  
  #===============Alignment============
  tabPanel(title = "Alignment",
           titlePanel("Alignment"),
           tabsetPanel(type = "tabs",
                       tabPanel("Alignment Options", 
                                actionButton("input_pre_alignments", "I have pre-existing alignments!"),
                                tags$div(id = "alignment_menu",
                                         hidden(selectizeInput("alignment_method", 
                                                               label = "Alignment Method",
                                                               multiple = FALSE, 
                                                               choices = list("reference", "denovo", ""), 
                                                               selected = "")),
                                         hidden(tags$div(id = "denovo_paired_input_panel",
                                                         fluidRow(inputPanel(
                                                           numericInput("M", "M, the number of mismatches allowed between alignments within individuals", 3, 1, step = 1),
                                                           numericInput("n", "n, the number of mismatches allowed between alignments between individuals", 3, 1, step = 1),
                                                           checkboxInput("check_headers", "Check that fastq headers are STACKS acceptable. If you are sure they are OK, uncheck this.", value = TRUE),
                                                           checkboxInput("stacks_cleanup", "Clean-up accessory files after completion.", value = TRUE),
                                                           numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
                                                           checkboxInput("rmdup", "Remove PCR duplicates?", value = TRUE),
                                                           checkboxInput("rmimp", "Remove improperly paired reads? This should usually be FALSE.", value = FALSE),
                                                           numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
                                                         )))),
                                         hidden(tags$div(id = "denovo_single_input_panel",
                                                         fluidRow(inputPanel(
                                                           numericInput("M", "M, the number of mismatches allowed between alignments within individuals", 3, 1, step = 1),
                                                           numericInput("n", "n, the number of mismatches allowed between alignments between individuals", 3, 1, step = 1),
                                                           checkboxInput("check_headers", "Check that fastq headers are STACKS acceptable. If you are sure they are OK, uncheck this.", value = TRUE),
                                                           checkboxInput("stacks_cleanup", "Clean-up accessory files after completion.", value = TRUE),
                                                           numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
                                                           numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
                                                         )))),
                                         hidden(tags$div(id = "reference_paired_input_panel",
                                                         sidebarLayout(sidebarPanel = sidebarPanel(
                                                           fluidRow(inputPanel(
                                                             hidden(shinyFilesButton("reference_genome", "Select Reference Genome", title = "Select Reference Genome", multiple = FALSE)),
                                                             numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
                                                             checkboxInput("rmdup", "Remove PCR duplicates?", value = TRUE),
                                                             checkboxInput("rmimp", "Remove improperly paired reads?", value = TRUE),
                                                             numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
                                                           ))), 
                                                           mainPanel = mainPanel(hidden(textOutput("reference_files_report")))))),
                                         hidden(tags$div(id = "reference_single_input_panel",
                                                         sidebarLayout(sidebarPanel = sidebarPanel(
                                                           fluidRow(inputPanel(
                                                             hidden(shinyFilesButton("reference_genome1", "Select Reference Genome", title = "Select Reference Genome", multiple = FALSE)),
                                                             numericInput("mapQ", "mapQ, minumum acceptable mapping quality score", 5, min = 0, step = 1),
                                                             numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
                                                           ))),
                                                           mainPanel = mainPanel(hidden(textOutput("reference_files_report1")))))),
                                         
                                         hidden(actionButton("run_alignment", "Run Alignment!")),
                                ),
                                # pre-existing alignments
                                hidden(tags$div(id = "pre_alignments_input_panel",
                                                inputPanel(
                                                  column(12,
                                                         hidden(shinyFilesButton("pre_existing_alignments", label = "Select Alignments", multiple = TRUE, title = "Select Alignments (.bam files)")),
                                                         hidden(actionButton("parse_pre_existing_alignments", "Read In Alignment Files"))))))
                                
                       ),
                       tabPanel("Demultiplexed File Report", 
                                fluidRow(column(6, 
                                                hidden(uiOutput("demultiplexed_files_report_RA"))),
                                         
                                         column(6, 
                                                hidden(uiOutput("demultiplexed_files_report_RB")))))
           )
  ),
  
  #==============Genotyping===============
  tabPanel("Genotyping",
           titlePanel("Genotyping"),
           tabsetPanel(type = "tabs",
                       tabPanel("Genotyping Options",
                                inputPanel(
                                  # standard options
                                  fluidRow(column(6,
                                                  selectInput("genotyper", 
                                                              label = "Which genotyper would you like to use?",
                                                              choices = c("GATK", "SAMtools", "SOAPsnp", "phys", "sample"))),
                                           column(6,
                                                  selectInput("doGeno", label = "What kind of genotype outputs do you want?",
                                                              choices = list(`Major and Minor Alleles Only` = "major_minor",
                                                                             `Numeric (0, 1, or 2; the count of the minor allele)` = "numeric",
                                                                             `Allele identities (AA, AG, or GG; the called genotypes)` = "NN",
                                                                             `Posterior probabilites for each genotype` = "all_posteriors",
                                                                             `Posterior probability of only the called genotype` = "called_posterior",
                                                                             `Posterior probabilites for each genotype in binary` = "binary"),
                                                              multiple = FALSE, selected = "NN"))),
                                  
                                  fluidRow(column(6,
                                                  uiOutput("minInd"),
                                                  numericInput("SNP_pval", "SNP p-value cuttoff (sites with p-values indicating they are less likely to be SNPs will be rejected)", 
                                                               value = 1e-08, min = 1e-12, max = 1)),
                                           column(6,
                                                  numericInput("postCutoff", "Posterior Genotype Probability Cuttoff (sites where the genotype call is less confident than this will be rejected)", 
                                                               value = 0.95, min = 0, max = 1),
                                                  numericInput("minQg", "Minimum base-pair quality phred score", value = 20, min = 0, max = 60),
                                                  numericInput("minMapQ", "Minimum mapping quality phred score", value = 20, min = 0, max = 60))),
                                  
                                  fluidRow(column(12,
                                                  checkboxInput("unzip", "Unzip resulting genotype file?", value = TRUE),
                                                  checkboxInput("doVcf", "Make a VCF file from results?", value = FALSE),
                                                  checkboxInput("filter_paralogs", "Filter out potential paralogs?", value = FALSE)))),
                                
                                
                                hidden(tags$div(id = "paralog_filter_menu",
                                                inputPanel(
                                                  
                                                  # paralog options
                                                  
                                                  numericInput("filter_paralogs_alpha", "p-value cuttoff for a site to be considered a paralog", 
                                                               value = 0.0016, min = 1e-12, max = .5),
                                                  numericInput("filter_paralogs_buffer", "Number of bp to exclude from genotyping on either side of a paralog",
                                                               value = 1000, min = 1, max = 1000000, step = 1),
                                                  hidden(shinyFilesButton("paralog_reference", "Select the reference genome to which reads were aligned (if denovo, select the constructed denovo reference)", 
                                                                          title = "Select the reference genome", multiple = FALSE))))),
                                
                                # general options
                                inputPanel(
                                  column(12,
                                         hidden(shinyFilesButton("rf", "Optional: select a file containing regions of the reference/denovo genome to include during genotyping.", # no reason to hide--I'm going to show immediately, but I get the double entry error if I don't
                                                                 title = "Select a file containing regions of the reference/denovo genome to include", multiple = FALSE)),
                                         numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1))),
                       
                       # run
                       hidden(actionButton("run_genotyping", "Run Genotyping!"))),
  
  
  tabPanel("Alignment File Report", 
           fluidRow(hidden(uiOutput("aligned_files_report")))))),
  
  #============navbar options========
  id = "MainTabs"
  #============end UI================
)





server <- function(input, output, session) {
  show("rf")
  #==========read in input files for demultiplexing===========
  x <- reactiveValues(files = list(), 
                      files_report = list(), 
                      files_are_good = 0, 
                      have_sample_ids = 0,
                      demultiplexed = 0,
                      aligned = 0)
  
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
  
  # remove input buttons and inputs if barcode status changes
  observeEvent(input$barcodes,{
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
    x$files_are_good <- 0
    x$files <- list()
    x$files_report <- list()
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
    }
    else{
      hide("barcode_files_report")
      hide("fastq_files_report")
      hide("go_demultiplex")
      hide("sample_ids")
    }
  })
  shinyFileChoose(input, 'sample_ids', root=c(root))

  # import and note that we have sample ids if we do
  observeEvent(is.list(input$sample_ids), ignoreInit = TRUE,{
    if(is.list(input$sample_ids)){
      x$files$sample_ids <- .parse_shinyFiles_path(root, input$sample_ids)
      x$have_sample_ids <- 1
    }
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
      updateTabsetPanel(session, "MainTabs", "Alignment")
    }
  })
  
  #==========import existing demultiplexed reads==============
  # panel layout
  
  # show panel only if pre_demultiplexed file input is selected
  observeEvent(input$input_pre_demultiplexed,{
    toggle("pre_demulti_input_panel")
    toggle("demultiplexing")
    toggle("pre_existing_RA") # no reason to do it this way, but if I don't I get doubled inputs, and this works...
  })
  
  observeEvent(input$is.paired.pre.existing,{
    toggle("pre_existing_RB")
    hide("denovo_paired_input_panel")
    hide("denovo_single_input_panel")
    hide("reference_paired_input_panel")
    hide("reference_single_input_panel")
    reset("reference_genome")
    reset("reference_genome1")
    hide("alignment_method")
    hide("run_alignment")
  })
  
  # show 'go' button if correct files have been selected given paired/unpaired
  observeEvent(is.list(input$pre_existing_RA) | is.list(input$pre_existing_RB), ignoreInit = TRUE, {
    hide("denovo_paired_input_panel")
    hide("denovo_single_input_panel")
    hide("reference_paired_input_panel")
    hide("reference_single_input_panel")
    reset("reference_genome")
    reset("reference_genome1")
    hide("alignment_method")
    hide("run_alignment")
    
    
    
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
                 
                 updateTabsetPanel(session, "MainTabs", "Alignment")
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
      show("alignment_method")
      
    }
    else{
      hide("demultiplexed_files_report_RA")
      hide("demultiplexed_files_report_RB")
      hide("alignment_method")
      hide("run_alignment")
    }
  })
  
  #==========prep alignment============================
  # show or hide correct input panels based on denovo or reference, paired or unpaired inputs
  observeEvent(input$alignment_method, ignoreInit = TRUE, {
    hide("denovo_paired_input_panel")
    hide("denovo_single_input_panel")
    hide("reference_paired_input_panel")
    hide("reference_single_input_panel")
    reset("reference_genome")
    hide("reference_genome")
    reset("reference_genome1")
    hide("reference_genome1")
    hide("run_alignment")

    if(input$alignment_method == "denovo"){
      if(length(x$files$demultiplexed_files) == 2){
        show("denovo_paired_input_panel")
      }
      else{
        show("denovo_single_input_panel")
      }
    }
    else if(input$alignment_method == "reference"){
      if(length(x$files$demultiplexed_files) == 2){
        show("reference_paired_input_panel")
        show("reference_genome")
      }
      else{
        show("reference_single_input_panel")
        show("reference_genome1")
      }
      
    }
  })
  
  shinyFileChoose(input, 'reference_genome', root=c(root))
  shinyFileChoose(input, 'reference_genome1', root=c(root))
  
  
  #==========run alignment================
  # check that we are good to go before offering the option. Using observe here because need to recheck whenever reference genomes change.
  observe({
    hide("run_alignment")
    if(x$demultiplexed){
      if(input$alignment_method == "reference"){
        if(length(x$files$demultiplexed_files) == 2 & is.list(input$reference_genome)){
          x$files$reference_genome <- .parse_shinyFiles_path(root, input$reference_genome)
          show("run_alignment")
        }
        else if(is.list(input$reference_genome1)){
          x$files$reference_genome <- .parse_shinyFiles_path(root, input$reference_genome1)
          show("run_alignment")
        }
      }
      else if(input$alignment_method == "denovo"){
        show("run_alignment")
      }
    }
  })
  
  # generate reporter for the reference genome
  observeEvent(x$files$reference_genome, ignoreInit = TRUE, {
    hide("reference_files_report")
    hide("reference_files_report1")
    if(!is.null(x$files$reference_genome)){
      x$files_report$reference <- .fastq_file_reporter(x$files)$reference
      
      if(length(x$files$demultiplexed_files) == 2){
        output$reference_files_report <- renderText(x$files_report$reference)
        show("reference_files_report")
      }
      else{
        output$reference_files_report1 <- renderText(x$files_report$reference)
        show("reference_files_report1")
      }
    }
  })
  
  # run the alignment
  observeEvent(input$run_alignment,{
    x$aligned <- 0
    x$files$alignments <- list()
    if(input$alignment_method == "denovo"){
      x$files$alignments <- align_denovo(RA_fastqs = x$files$demultiplexed_files$RA,
                                         RB_fastqs = x$files$demultiplexed_files$RB, 
                                         M = input$M, 
                                         n = input$n, 
                                         par = input$par,
                                         check_headers = input$check_headers,
                                         stacks_cleanup = input$stacks_cleanup, 
                                         re_align = TRUE,
                                         mapQ = input$mapQ, 
                                         remove_duplicates = input$rmdup,
                                         remove_improper_pairs = input$rmimp)
      x$aligned <- 1
    }
    else if(input$alignment_method == "reference"){
      x$files$alignments <- align_reference(RA_fastqs = x$files$demultiplexed_files$RA,
                                            RB_fastqs = x$files$demultiplexed_files$RB,
                                            reference = x$files$reference_genome, 
                                            remove_duplicates = input$rmdup, 
                                            mapQ = input$mapQ, 
                                            remove_improper_pairs = input$rmimp, 
                                            par = input$par)
      x$aligned <- 1
    }
    updateTabsetPanel(session, "MainTabs", "Genotyping")
  })
  
  #==========import existing alignments============
  observeEvent(input$input_pre_alignments,{
    toggle("pre_alignments_input_panel")
    toggle("pre_existing_alignments")
    toggle("alignment_menu")
  })
  
  shinyFileChoose(input, id = "pre_existing_alignments", root = c(root))
  
  # import alignments
  observeEvent(is.list(input$pre_existing_alignments), ignoreInit = TRUE, {
    hide("parse_pre_existing_alignments")
    x$aligned <- 0
    
    if(is.list(input$pre_existing_alignments)){
      x$files$alignments <- .parse_shinyFiles_path(root, input$pre_existing_alignments)
      x$aligned <- 1
      show("parse_pre_existing_alignments")
      updateTabsetPanel(session, "MainTabs", "Genotyping")
    }
  })
  
  #==========report on alignment=========================
  output$aligned_files_report <- renderUI({
    lapply(x$files_report$alignments, p)
  })
  
  observeEvent(x$aligned,{
    hide("aligned_files_report")
    if(x$aligned == 1){
      x$files_report$alignments <- .fastq_file_reporter(x$files)$alignments
      show("aligned_files_report")
    }
  })
  
  
  #==========prep for genotyping===========
  shinyFileChoose(input, "paralog_reference", root = root)
  shinyFileChoose(input, "rf", root = root)
  
  # show or hide the paralog filtering menu
  # toggle doesn't behave like I'd expect
  observeEvent(input$filter_paralogs,{
    if(input$filter_paralogs){
      show("paralog_filter_menu")
      show("paralog_reference")
    }
    else{
      hide("paralog_filter_menu")
      hide("paralog_reference")
    }
    
  })
  
  # save the paralog ref if provided
  observeEvent(is.list(input$paralog_reference), ignoreInit = TRUE, {
    hide("paralog_reference_files_report")
    if(is.list(input$paralog_reference)){
      x$files$paralog_reference <- .parse_shinyFiles_path(root, input$paralog_reference)
      show("paralog_reference_files_report")
    }
  })
  
  # report on paralog reference
  observeEvent(x$files$paralog_reference, ignoreNULL = TRUE, ignoreInit = TRUE,{
    x$files_report$paralog_reference <- .fastq_file_reporter(x$files)$paralog_reference
  })
  output$paralog_reference_files_report <- renderText(x$files_report$paralog_reference)
  
  # show or hide the run button
  # observe({
  #   hide("run_genotyping")
  #   if(x$aligned){
  #     if(input$filter_paralogs){
  #       
  #     }
  # 
  #   }
  # })
  
} 

shinyApp(ui, server)
