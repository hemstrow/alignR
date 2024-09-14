library(shiny)
library(shinyjs)
library(shinyFiles)
library(shinybusy)
library(alignR)
source("app_code.R")

root <- c(wd = normalizePath("."), root = "/", home = normalizePath("~"))

# things to add: Show head of .fastq files and barcodes once selected

ui <- navbarPage(

  tags$head(tags$style(".Fixed_input_row_75{height:75px;}"),
            tags$style(".Fixed_input_row_150{height:150px;}")),

  #===============Demultiplexing and initial inputs===============
  tabPanel(title = "Demultiplexing",
           useShinyjs(), # here so it doesn't throw a warning, still works with every tab!
           titlePanel("Demultiplexing"),
           actionButton("input_pre_demultiplexed", "I have demultiplexed fastq files!"),
           actionButton("use_exdata", "Use example data"),

           # barcode selection and input panels
           tags$div(id = "demultiplexing",
                    sidebarLayout(
                      sidebarPanel(
                        fluidRow(selectizeInput("barcodes",
                                                "Where are your barcodes located?",
                                                choices = list("Headers",
                                                               "Start of Reads",
                                                               "R2 fastq file"),
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

                    numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
                    sidebarLayout(sidebarPanel(hidden(shinyFilesButton("sample_ids", label = "Optional Sample Identifications", multiple = TRUE,
                                                                       title = "Tab delimited key(s) for demultiplexed samples, one per fastq input file. Two or three columns for: sample IDs, R2 fastq file barcodes, header/in-read barcodes, in that order."))),
                                  mainPanel(htmlOutput("sample_ids_files_report"))),
                    uiOutput("go_demultiplex")
           ),


           # pre-demultiplexed input
           hidden(tags$div(id = "pre_demulti_input_panel",
                           sidebarLayout(
                             sidebarPanel = sidebarPanel(
                               class = "Fixed_input_row_150",
                               width = 4,
                               column(12,
                                      checkboxInput("is.paired.pre.existing", "Is your data paired-end?", value = TRUE),
                                      hidden(shinyFilesButton("pre_existing_RA", label = "RA Files", multiple = TRUE, title = "Select RA Reads")),
                                      hidden(shinyFilesButton("pre_existing_RB", label = "RB Files", multiple = TRUE, title = "Select RB Reads")),
                                      hidden(actionButton("parse_pre_existing_demultiplexed", "Continue to alignment!")))),
                             mainPanel = mainPanel(
                               width = 8,
                               fluidRow(column(6,htmlOutput("pre_fastq_files_report_RA")), column(6,hidden(htmlOutput("pre_fastq_files_report_RB"))))))))
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

                                         hidden(uiOutput("run_alignment"))),

                                # pre-existing alignments
                                hidden(tags$div(id = "pre_alignments_input_panel",
                                                inputPanel(
                                                  column(12,
                                                         hidden(shinyFilesButton("pre_existing_alignments", label = "Select Alignments", multiple = TRUE, title = "Select Alignments (.bam files)")),
                                                         hidden(actionButton("continue_to_genotyping", "Continue to genotyping!"))))))),

                       tabPanel("Demultiplexed File Report",
                                fluidRow(column(6,
                                                hidden(uiOutput("demultiplexed_files_report_RA"))),

                                         column(6,
                                                hidden(uiOutput("demultiplexed_files_report_RB"))))))
  ),

  #==============Genotyping===============
  tabPanel("Genotyping",
           titlePanel("Genotyping"),
           tabsetPanel(type = "tabs", id = "genotyper_panel",
                       #===================ANGSD=====================
                       tabPanel("Genotyping Options -- ANGSD", value = "ANGSD",
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
                                                checkboxInput("filter_paralogs", "Filter out potential paralogs?", value = FALSE))),


                                hidden(tags$div(id = "paralog_filter_menu",
                                                sidebarLayout(

                                                  # paralog options
                                                  sidebarPanel = sidebarPanel(
                                                    h3("Paralog filtering options:"),
                                                    numericInput("filter_paralogs_alpha", "p-value cuttoff for a site to be considered a paralog",
                                                                 value = 0.0016, min = 1e-12, max = .5),
                                                    numericInput("filter_paralogs_buffer", "Number of bp to exclude from genotyping on either side of a paralog",
                                                                 value = 1000, min = 1, max = 1000000, step = 1),
                                                    hidden(shinyFilesButton("paralog_reference", title = "Select the reference genome to which reads were aligned (if denovo, select the constructed denovo reference)",
                                                                            label = "Select the reference genome", multiple = FALSE)),
                                                    hidden(shinyFilesButton("filter_paralogs_pops", label = "Optional poulation ID file", title = "Select a file containing population IDs for each alignment.", multiple = FALSE))),

                                                  # reports
                                                  mainPanel = mainPanel(
                                                    hidden(uiOutput("paralog_reference_files_report")),
                                                    hidden(uiOutput("filter_paralogs_pops_files_report")))))),

                                # general options
                                sidebarLayout(
                                  sidebarPanel = sidebarPanel(width = 5,
                                                              hidden(shinyFilesButton("rf", "Optional file of gemonic regions to include.", # no reason to hide--I'm going to show immediately, but I get the double entry error if I don't
                                                                                      title = "Select a file containing regions of the reference/denovo genome to include", multiple = FALSE))),
                                  mainPanel = mainPanel(width = 7,
                                                        uiOutput("rf_files_report"))),

                                numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1)),


                       #================GATK==============
                       tabPanel("Genotyping Options -- GATK", value = "GATK",
                                # standard options
                                fluidRow(column(6,
                                                h3("Genome chunking options:"),
                                                hidden(shinyFilesButton("reference_genome_GATK", "Select Reference Genome", title = "Select Reference Genome", multiple = FALSE)),
                                                numericInput("chunk_size", "Size to ghunk genome into blocks for genotyping (whole chromosomes/scaffold if zero)?", value = 0, min = 0),
                                                numericInput("min_chunk_size", "Minimum chromosome/scaffold size to genotype?", value = 0, min = 0)),
                                         column(6,
                                                h3("Hard filtering options:"),
                                                numericInput("QD", "QD: Quality-by-depth cuttoff", value = 2, min = 0),
                                                numericInput("FS", "FS: Fisher strand-bias cuttoff", value = 60, min = 0),
                                                numericInput("SOR", "SOR: Strand odds-ratio cuttoff", value = 3),
                                                numericInput("MQ", "MQ: Mapping quality cuttoff", value = 40, min = 0),
                                                numericInput("MQRankSum", "MQRankSum: Mapping quality variance between ref and alt alleles cuttoff", value = -12.5),
                                                numericInput("ReadPosRankSum", "ReadPosRankSum: Ref/alt allele read position bias cuttoff ", value = -8),
                                                numericInput("min_genotype_quality", "QUAL: minimum genotype quality to retain", value = 13, min = 0)),
                                ),
                                # RG/batch control
                                fluidRow(column(4,
                                                h3("Read-groups:"),
                                                checkboxInput("add_RGs", "Do you need to add read-group tags to your alignments?", TRUE),
                                                textInput("platform", "Sequencing Platform", "ILLUMINA")),
                                         column(4,
                                                h3("Tool paths:"),
                                                hidden(shinyFilesButton("java_path", "path to java", title = "Locate Java", multiple = FALSE)),
                                                hidden(shinyFilesButton("gatk4_path", "path to GATK jar", title = "Locate GATK jar file", multiple = FALSE)),
                                                hidden(shinyFilesButton("picard_path", "path to picard", title = "Locate picard", multiple = FALSE))
                                         ),
                                         column(4,
                                                h3("Misc:"),
                                                numericInput("batch_size", "GenomicsDBImport number of samples to import at once", min = 1, value = 5),
                                                checkboxInput("concatenate_final_vcfs", "Concatenate final .vcf files?", TRUE)))
                       ),

  #===========Reporter===========
  tabPanel("Alignment File Report",
           fluidRow(hidden(uiOutput("aligned_files_report"))))),

  numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),

  # run
  sidebarLayout(sidebarPanel = sidebarPanel(uiOutput("run_genotyping")),
                mainPanel = mainPanel(uiOutput("genotypes_files_report")))),

  #============navbar options========
  id = "alignR"
  #============end UI================
)





server <- function(input, output, session) {

  #==========constant observers and initializations==========
  # initialize reactives
  x <- reactiveValues(files = list(),
                      files_report = list(),
                      files_are_good = 0,
                      have_sample_ids = 0,
                      demultiplexed = 0,
                      aligned = 0,
                      rename_r2 = 0,
                      paired = TRUE,
                      genotyped = 0)

  # file reporter
  observeEvent(x$files,{
    x$files_report <- .fastq_file_reporter(x$files, x$rename_r2)
  })

  # paired sync
  observe({
    updateCheckboxInput(session, "is.paired.pre.existing", value = input$is.paired)
    isolate(x$paired <- input$is.paired)
  })
  observe({
    updateCheckboxInput(session, "is.paired", value = input$is.paired.pre.existing)
    isolate(x$paired <- input$is.paired)
  })


  show("rf")
  #==========read in input files for demultiplexing===========


  # figure out which images we need
  barcode.type <- eventReactive(input$barcodes,{
    img_table <- data.frame(input = c("Headers",
                                      "Start of Reads",
                                      "R2 fastq file"),
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
  })



  # generate images
  output$barcode_locations_image <- renderUI({
    res <- list()
    validate(need(!all(c("Start of Reads", "Headers") %in% input$barcodes), "Barcodes can be either in fastq headers OR at the beginning of reads, but not both!"))
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
  observeEvent(x$paired, ignoreInit = TRUE, {
    x$files_are_good <- 0
  })

  # show or hide the appropriate input selectors and edit text labels
  ## define consistant notes
  output$fastq_R1_note <- renderText("Forward read files (R1):")
  output$fastq_R3_note <- renderText("Reverse read files (R3):")
  output$barcode_file_2_note <- renderText("Barcodes matching R2 (R2 fastq file barcodes):")

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

   if(!all(c("Start of Reads", "Headers") %in% input$barcodes)){
     if(length(input$barcodes) != 0){
       show("fastq_R1")
       show("fastq_R1_note")
     }


     if("Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes){
       if("Headers" %in% input$barcodes){
         output$barcode_file_1_note <- renderText("Header barcodes:")
       }
       else{
         output$barcode_file_1_note <- renderText("In-read barcodes:")
       }
       show("barcode_file_1_note")
       show("barcode_file_1")
     }

     if(x$paired){
       if(length(input$barcodes) != 0){
         show("fastq_R2")
         show("fastq_R2_note")
       }

       if("R2 fastq file" %in% input$barcodes){
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
     else if("R2 fastq file" %in% input$barcodes){
       output$barcode_file_2_note <- renderText("R2 barcodes:")
       output$fastq_R2_note <- renderText("Barcode reads (R2):")

       show("barcode_file_2_note")
       show("fastq_R2_note")
       show("fastq_R2")
       show("barcode_file_2")
     }
   }
  })

  # locate appropriate files, update reporter, check if we are good to proceed based on files and barcode status.
  observeEvent(eventExpr =
                 is.list(input$fastq_R1) |
                 is.list(input$fastq_R2) |
                 is.list(input$fastq_R3) |
                 is.list(input$barcode_file_1) |
                 is.list(input$barcode_file_2), ignoreInit = TRUE, handlerExpr = {

    files <- lapply(list(R1 = input$fastq_R1,
                         R2 = input$fastq_R2,
                         R3 = input$fastq_R3,
                         barcode_file_1 = input$barcode_file_1,
                         barcode_file_2 = input$barcode_file_2),
                    function(x) .parse_shinyFiles_path(root, x))

    files <- purrr::discard(files, is.null)
    x$rename_r2 <- 0


    # read files for every case:
    ## with in-line or header barcodes
    if(("Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes) & is.list(input$barcode_file_1) & is.list(input$fastq_R1)){

      # not paired, dual indexed
      if(!x$paired & "R2 fastq file" %in% input$barcodes & is.list(input$fastq_R2) & is.list(input$barcode_file_2)){
        x$files_are_good <- 1
      }

      # paired, dual indexed
      else if(x$paired & "R2 fastq file" %in% input$barcodes & is.list(input$fastq_R2) & is.list(input$fastq_R3) & is.list(input$barcode_file_2)){
        x$files_are_good <- 1
      }

      # not paired, single indexed
      else if(!x$paired & !"R2 fastq file" %in% input$barcodes){
        x$files_are_good <- 1
      }

      else if(x$paired & !"R2 fastq file" %in% input$barcodes & is.list(input$fastq_R2)){
        names(files)[which(names(files) == "R2")] <- "R3" # renamed R3 for consistency later, since this is reverse reads not barcodes
        x$rename_r2 <- 1
        x$files_are_good <- 1
      }
    }

    ## with only R2 barcodes
    else if(!("Headers" %in% input$barcodes | "Start of Reads" %in% input$barcodes) & is.list(input$barcode_file_1) & is.list(input$fastq_R1)){
      # not paired
      if(!x$paired){
        x$files_are_good <- 1
      }

      # paired
      if(x$paired & is.list(input$fastq_R2)){
        names(files)[which(names(files) == "R2")] <- "R3" # renamed R3 for consistency later, since this is reverse reads not barcodes
        x$rename_r2 <- 1
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
      show("sample_ids")
      show("go_demultiplex", asis = TRUE)
    }
    else{
      hide("sample_ids")
      hide("go_demultiplex")
    }
  })
  shinyFileChoose(input, 'sample_ids', root=c(root))

  # import and note that we have sample ids if we do
  observeEvent(is.list(input$sample_ids), ignoreInit = TRUE,{
    hide("sample_ids_files_report")
    if(is.list(input$sample_ids)){
      x$files$sample_ids <- .parse_shinyFiles_path(root, input$sample_ids)
      x$have_sample_ids <- 1
      show("sample_ids_files_report")
    }
  })

  output$sample_ids_files_report <- renderUI({
    lapply(c("Warning: Input R1/R2/R3 and sample names must all be in the same order! Check before continuing!", x$files_report$sample_ids), p)
  })


  # render the go_demultiplex button here to do some validation.
  output$go_demultiplex <- renderUI({

    length_prob <- 0
    if("R2" %in% names(x$files)){
      if(length(x$files$R2) != length(x$files$R1)){length_prob <- 1}
    }
    if("R3" %in% names(x$files)){
      if(length(x$files$R3) != length(x$files$R1)){length_prob <- 1}
    }
    if("sample_ids" %in% names(x$files)){
      if(length(x$files$R1) != length(x$files$sample_ids)){length_prob <- 1}
    }

    validate(need(expr = length_prob == 0, label = "The number of R1/R2/R3 files must be equal, as must the number of sample ID files if provided."))

    return(actionButton("go_demultiplex", "Demultiplex reads!"))
  })


  #==========run demultiplexing======================
  # demultiplex
  observeEvent(input$go_demultiplex,{
    have_indices <- FALSE

    if(x$files_are_good){

      if(length(input$barcodes) > 1){
        indices <- readLines(x$files$barcode_file_2)
        barcodes <- readLines(x$files$barcode_file_1)
      }
      else if("R2 fastq file" %in% x$barcode_file_1){
        indices <- readLines(x$barcode_file_1)
        barcodes <- NULL
      }
      else{
        indices <- NULL
        barcodes <- readLines(x$barcode_file_1)
      }

      if(x$paired){
        R3 <- x$files$R3
      }
      else{
        R3 <- NULL
      }

      if(x$have_sample_ids){
        sample_names <- vector("list", length = length(x$files$sample_ids))
        for(i in 1:length(sample_names)){
          sample_names[[i]] <- read.table(x$files$sample_ids[i], headers = FALSE)
        }
      }
      else{
        sample_names <- NULL
      }

      barcode_key <- data.frame(input = c("R2 fastq file", "Start of Reads", "Headers"),
                                option = c("R2", "read_start", "header"))
      barcode_locations <- barcode_key$option[match(input$barcodes, barcode_key$input)]

      show_modal_spinner(spin = "circle", color = "#9e5157", text = "Demultiplexing, please wait...")
      x$files$demultiplexed_files <- demultiplex(R1 = x$files$R1,
                                                 R2 = x$files$R2,
                                                 R3 = R3,
                                                 barcode_locations = barcode_locations,
                                                 indices = indices,
                                                 barcodes = barcodes,
                                                 sample_names = sample_names,
                                                 outfile_prefix = "alignR_gui_demulti",
                                                 par = input$par)
      remove_modal_spinner()

      x$demultiplexed <- 1
      updateTabsetPanel(session, "alignR", "Alignment")
    }
  })

  #==========import existing demultiplexed reads==============
  # show panel only if pre_demultiplexed file input is selected
  observeEvent(input$input_pre_demultiplexed,{
    toggle("pre_demulti_input_panel")
    toggle("demultiplexing")
    toggle("pre_existing_RA") # no reason to do it this way, but if I don't I get doubled inputs, and this works...
  })

  observeEvent(x$paired, {
    toggle("pre_existing_RB")
    toggle("pre_fastq_files_report_RB")
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



    if(x$paired){
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
  observeEvent(is.list(input$pre_existing_RA) |
                 is.list(input$pre_existing_RB), ignoreInit = TRUE, ignoreNULL = TRUE,{
                   x$files$demultiplexed_files <- NULL

                   RA <- .parse_shinyFiles_path(root, input$pre_existing_RA)
                   RB <- .parse_shinyFiles_path(root, input$pre_existing_RB)

                   if(is.null(RB) & !is.null(RA)){
                     x$files$demultiplexed_files <- list(RA = .parse_shinyFiles_path(root, input$pre_existing_RA))
                   }
                   else if(!is.null(RB) & !is.null(RA)){
                     x$files$demultiplexed_files <- vector("list", 2)
                     names(x$files$demultiplexed_files) <- c("RA", "RB")
                     x$files$demultiplexed_files$RA <- RA
                     x$files$demultiplexed_files$RB <- RB
                   }

                 })


  output$pre_fastq_files_report_RA <- renderUI({
    lapply(x$files_report$demulti$RA, p)
  })
  output$pre_fastq_files_report_RB <- renderUI({
    lapply(x$files_report$demulti$RB, p)
  })

  # check if we can proceed
  observeEvent(eventExpr = input$parse_pre_existing_demultiplexed,
               ignoreInit = TRUE, handlerExpr = {

                 x$demultiplexed <- 0
                 if(x$paired){
                   if(is.list(input$pre_existing_RA) & is.list(input$pre_existing_RB)){

                     validate(need(length(x$files$demultiplexed_files$RA) == length(x$files$demultiplexed_files$RB), message = "An equal number of RA and RB files must be provided."))

                     x$demultiplexed <- 1
                   }
                 }
                 else{
                   if(is.list(input$pre_existing_RA)){
                     x$demultiplexed <- 1
                   }
                 }

                 # reset our file inputs for validation purposes
                 reset("pre_existing_RB")
                 reset("pre_existing_RA")

                 updateTabsetPanel(session, "alignR", "Alignment")
               })

  #==========read in example data if requested=========
  observeEvent(input$use_exdata,{
    toggle("demultiplexing")
    toggle("pre_existing_RA")

    R1_files <- list.files(system.file("exdata", package = "alignR"), full.names = TRUE)
    file.copy(R1_files, ".")
    R1_files <- basename(R1_files)
    x$files$demultiplexed_files$RB <- R1_files[grep("_R2\\.fastq$", R1_files)]
    x$files$demultiplexed_files$RA <- R1_files[grep("_R1\\.fastq$", R1_files)]
    x$files$reference_genome <- R1_files[grep("\\.fna\\.gz", R1_files)]

    x$demultiplexed <- 1

    # reset our file inputs for validation purposes
    reset("pre_existing_RB")
    reset("pre_existing_RA")

    updateTabsetPanel(session, "alignR", "Alignment")
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
      show("demultiplexed_files_report_RA")
      show("alignment_method")

      if(x$paired){
        show("demultiplexed_files_report_RB")
      }
      else{
        hide("demultiplexed_files_report_RB")
      }
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
      if(x$paired){
        show("denovo_paired_input_panel")
      }
      else{
        show("denovo_single_input_panel")
      }
    }
    else if(input$alignment_method == "reference"){
      if(x$paired){
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
        if(x$paired & is.list(input$reference_genome)){
          x$files$reference_genome <- .parse_shinyFiles_path(root, input$reference_genome)
          show("run_alignment")
        }
        else if(is.list(input$reference_genome1)){
          x$files$reference_genome <- .parse_shinyFiles_path(root, input$reference_genome1)
          show("run_alignment")
        }
        else if(!is.null(x$files$reference_genome)){
          if(file.exists(x$files$reference_genome)){
            show("run_alignment")
          }
        }
      }
      else if(input$alignment_method == "denovo"){
        show("run_alignment")
      }
    }
  })

  # define run_alignment (allow for validation)

  output$run_alignment <- renderUI({
    validate(need(.check_system_install("bwa") & .check_system_install("samtools"), "No system install of bwa or samtools detected. Alignment disabled."))

    if(input$alignment_method == "denovo"){
      validate(need(!isFALSE(.check_system_install("stacks")), "No system install of stacks detected. Denovo alignment disabled."))
    }

    return(actionButton("run_alignment", "Run Alignment!"))
  })


  # generate reporter for the reference genome
  observeEvent(x$files$reference_genome, ignoreInit = TRUE, {
    hide("reference_files_report")
    hide("reference_files_report1")
    if(!is.null(x$files$reference_genome)){

      if(x$paired){
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
    show_modal_spinner(spin = "circle", color = "#9e5157", text = "Aligning, please wait...")
    if(input$alignment_method == "denovo"){
      if(x$paired){
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
      }
      else{
        x$files$alignments <- align_denovo(RA_fastqs = x$files$demultiplexed_files$RA,
                                           RB_fastqs = NULL,
                                           M = input$M,
                                           n = input$n,
                                           par = input$par,
                                           check_headers = input$check_headers,
                                           stacks_cleanup = input$stacks_cleanup,
                                           re_align = TRUE,
                                           mapQ = input$mapQ,
                                           remove_duplicates = input$rmdup,
                                           remove_improper_pairs = input$rmimp)
      }

      x$aligned <- 1
    }
    else if(input$alignment_method == "reference"){
      if(x$paired){
        x$files$alignments <- align_reference(RA_fastqs = x$files$demultiplexed_files$RA,
                                              RB_fastqs = x$files$demultiplexed_files$RB,
                                              reference = x$files$reference_genome,
                                              remove_duplicates = input$rmdup,
                                              mapQ = input$mapQ,
                                              remove_improper_pairs = input$rmimp,
                                              par = input$par)
      }
      else{
        x$files$alignments <- align_reference(RA_fastqs = x$files$demultiplexed_files$RA,
                                              RB_fastqs = NULL,
                                              reference = x$files$reference_genome,
                                              remove_duplicates = input$rmdup,
                                              mapQ = input$mapQ,
                                              remove_improper_pairs = input$rmimp,
                                              par = input$par)
      }

      x$aligned <- 1
    }

    remove_modal_spinner()
    updateTabsetPanel(session, "alignR", "Genotyping")
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
    x$aligned <- 0
    hide("continue_to_genotyping")

    if(is.list(input$pre_existing_alignments)){
      x$files$alignments <- .parse_shinyFiles_path(root, input$pre_existing_alignments)
      x$aligned <- 1
      show("continue_to_genotyping")
    }
  })

  observeEvent(input$continue_to_genotyping, ignoreInit = TRUE,{
    updateTabsetPanel(session, "alignR", "Genotyping")
  })

  #==========report on alignment=========================
  output$aligned_files_report <- renderUI({
    lapply(x$files_report$alignments, p)
  })

  observeEvent(x$aligned,{
    hide("aligned_files_report")
    if(x$aligned == 1){
      show("aligned_files_report")
    }
  })


  #==========prep for genotyping--ANGSD===========
  shinyFileChoose(input, "paralog_reference", root = root)
  shinyFileChoose(input, "filter_paralogs_pops", root = root)
  shinyFileChoose(input, "rf", root = root)

  # show or hide the paralog filtering menu
  observeEvent(input$filter_paralogs,{
    if(input$filter_paralogs){
      show("paralog_filter_menu")
      show("paralog_reference")
      show("filter_paralogs_pops")
    }
    else{
      hide("paralog_filter_menu")
      hide("paralog_reference")
      hide("filter_paralogs_pops")
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
  observeEvent(x$files_report$paralog_reference, ignoreNULL = TRUE, ignoreInit = TRUE,{
    output$paralog_reference_files_report <- renderText(paste0("Reference for paralog filtering: ", x$files_report$paralog_reference))
  })

  # save the paralog pops if provided
  observeEvent(is.list(input$filter_paralogs_pops), ignoreInit = TRUE, {
    hide("filter_paralogs_pops_files_report")
    if(is.list(input$filter_paralogs_pops)){
      x$files$filter_paralogs_pops <- .parse_shinyFiles_path(root, input$filter_paralogs_pops)
      show("filter_paralogs_pops_files_report")
    }
  })

  # report on paralog pops
  observeEvent(x$files_report$filter_paralogs_pops, ignoreNULL = TRUE, ignoreInit = TRUE,{
    output$filter_paralogs_pops_files_report <- renderText(paste0("Pop ID file: ", x$files_report$filter_paralogs_pops))
  })

  # generate minInd selector
  output$minInd <- renderUI({
    validate(need(!is.null(x$files$alignments), message = FALSE))
    max_val <- length(x$files$alignments)
    init <- floor(max_val/2)
    return(numericInput("minInd", "Minimum number of genotyped individuals per locus.", value = init, max = max_val, min = 0, step = 1))
  })

  # save the rf if provided
  observeEvent(is.list(input$rf), {
    hide("rf_files_report")
    if(is.list(input$rf)){
      x$files$rf <- .parse_shinyFiles_path(root, input$rf)
      show("rf_files_report")
    }
  })

  # report on rf
  observeEvent(x$files_report$rf, ignoreNULL = TRUE, ignoreInit = TRUE,{
    output$rf_files_report <- renderText(paste0("Region file: ", x$files_report$rf))
  })

  #===============prep for genotyping, GATK============
  shinyFileChoose(input, "reference_genome_GATK", root = c(root))
  shinyFileChoose(input, "java_path", root = root)
  shinyFileChoose(input, "gatk4_path", root = root)
  shinyFileChoose(input, "picard_path", root = root)

  # show or hide genome input
  observe({
    hide("reference_genome_GATK")
    if(is.null(x$files$reference_genome)){
      show("reference_genome_GATK")
    }
    else{
      if(!file.exists(x$files$reference_genome)){
        show("reference_genome_GATK")
      }
    }
  })

  # parse genome if added
  observeEvent(input$reference_genome_GATK, ignoreNULL = TRUE, ignoreInit = TRUE, {
    x$files$reference_genome <- .parse_shinyFiles_path(root, input$reference_genome_GATK)
  })

  # show or hide add_RG platform
  hide("platform")
  observeEvent(input$add_RGs, {
    if(input$add_RGs){
      show("platform")
    }
    else{
      hide("platform")
    }
  })

  # read path selectors
  show("java_path")
  observeEvent(input$java_path, ignoreNULL = TRUE, ignoreInit = TRUE, {
    x$files$java_path <- .parse_shinyFiles_path(root, input$java_path)
  })
  show("gatk4_path")
  observeEvent(input$gatk4_path, ignoreNULL = TRUE, ignoreInit = TRUE, {
    x$files$gatk4_path <- .parse_shinyFiles_path(root, input$gatk4_path)
  })
  show("picard_path")
  observeEvent(input$picard_path, ignoreNULL = TRUE, ignoreInit = TRUE, {
    x$files$picard_path <- .parse_shinyFiles_path(root, input$picard_path)
  })


  #===============run genotyping==================
  # make the run_genotyping button here for validation
  output$run_genotyping <- renderUI({
    if(x$aligned){
      if(input$genotyper_panel == "ANGSD"){
        if(input$filter_paralogs){
          validate(need(.check_system_install("ngsParalog"), "No ngsParalog installation detected. Genotyping with paralog filtering disabled."))
          validate(need(.check_system_install("samtools"), "No samtools installation detected. Genotyping with paralog filtering disabled."))


          if(!is.null(x$files$paralog_reference)){
            if(!file.exists(x$files$paralog_reference)){
              validate(need(FALSE, "Could not locate paralog reference file"))
            }
          }
          else{
            validate(need(FALSE, "Please provide reference genome for paralog filtering."))
          }

          if(!is.null(x$files$filter_paralog_pops)){
            if(!file.exists(x$files$paralog_reference)){
              validate(need(FALSE, "Could not locate paralog pops file"))
            }
            else{
              validate(need(length(readLines(x$files$filter_paralogs_pops)) == length(x$files$alignments),
                            "The filter paralogs populations file must have the same number of lines (samples) as the number of alignments."))
            }
          }
        }


        validate(need(.check_system_install("angsd"), "No angsd installation detected. Genotyping disabled."))
        if(input$doVcf & input$doGeno %in% c("all_posteriors", "called_posterior", "binary")){
          validate(need(.check_system_install("bcftools"), "No bcftools installation detected. doVcf with non-numeric or allele identities output options disabled."))
        }
      }
      else if(input$genotyper_panel == "GATK"){
        if(is.null(x$files$reference_genome)){
          validate(need(FALSE, "Please provide a reference genome!"))
        }
        else if(!file.exists(x$files$reference_genome)){
          validate(need(FALSE, "Reference genome file not found!"))
        }

        if(input$add_RGs){
          if(is.null(x$files$demultiplexed_files$RA)){
            validate(need(FALSE, ".fastq/.fasta files required to add RGs (R1 reads, see demultiplexing panel)."))
          }
          if(any(!file.exists(x$files$demultiplexed_files$RA))){
            validate(need(FALSE, ".fastq/.fasta files not located."))
          }
        }

        if(is.null(x$files$gatk4_path)){
          validate(need(FALSE, "Please provide the path to the gatk4 .jar!"))
        }
        else if(!file.exists(x$files$gatk4_path)){
          validate(need(FALSE, "gatk4 .jar file not found!"))
        }

        if(is.null(x$files$picard_path)){
          validate(need(FALSE, "Please provide the path to the picard .jar!"))
        }
        else if(!file.exists(x$files$picard_path)){
          validate(need(FALSE, "picard .jar file not found!"))
        }

        if(is.null(x$files$java_path)){
          validate(need(FALSE, "Please provide the path to java!"))
        }
        else if(!file.exists(x$files$java_path)){
          validate(need(FALSE, "java not found!"))
        }

        validate(need(.check_system_install("vcftools"), "vcftools must be installed and on the system path!"))
        validate(need(.check_system_install("samtools"), "samtools must be installed and on the system path!"))
        validate(need(.check_system_install("bwa"), "bwa must be installed and on the system path!"))
        validate(need(.check_system_install("bcftools"), "bcftools must be installed and on the system path!"))
      }

    }

    else{
      validate(need(FALSE, "Please provide alignments!"))
    }


    return(actionButton("run_genotyping", "Run Genotyping!"))
  })

  # run genotyping
  observeEvent(input$run_genotyping,{
    show_modal_spinner(spin = "circle", color = "#9e5157", text = "Genotyping, please wait...")

    if(input$genotyper_panel == "ANGSD"){
      if(is.null(x$files$rf)){
        rf <- FALSE
      }
      else{
        rf <- x$files$rf
      }

      if(!is.null(x$files$filter_paralogs_pops)){
        filter_paralogs_pops <- x$files$filter_paralogs_pops
      }
      else{
        filter_paralogs_pops <- rep("only_pop", length(x$files$alignments))
      }

      if(!is.null(x$files$paralog_reference)){
        paralog_reference <- readLines(x$files$paralog_reference)
      }
      else{
        paralog_reference <- FALSE
      }

      x$files$genotypes <- genotype_bams_ANGSD(x$files$alignments,
                                               outfile_prefix = "alignR_gui_genotypes",
                                               minInd = input$minInd,
                                               SNP_pval = input$SNP_pval,
                                               doGeno = input$doGeno,
                                               genotyper = input$genotyper,
                                               postCutoff = input$postCutoff,
                                               minQ = input$minQg,
                                               minMapQ = input$minMapQ,
                                               unzip = input$unzip,
                                               doVcf = input$doVcf,
                                               filter_paralogs = input$filter_paralogs,
                                               filter_paralogs_alpha = input$filter_paralogs_alpha,
                                               filter_paralogs_buffer = input$filter_paralogs_buffer,
                                               filter_paralogs_populations = filter_paralogs_pops,
                                               filter_paralogs_reference = paralog_reference,
                                               rf = rf,
                                               par = input$par)

      remove_modal_spinner()
      x$genotyped <- 1
    }
    else if(input$genotyper_panel == "GATK"){
      browser()
      x$files$genotypes <- genotype_bams_GATK(x$files$alignments,
                                              fastqs = x$files$demultiplexed_files$RA,
                                              reference = x$files$reference_genome,
                                              par = input$par,
                                              min_chr_size = input$min_chunk_size,
                                              chunk_size = ifelse(input$chunk_size == 0, "chr", input$chunk_size),
                                              QD = input$QD,
                                              FS = input$FS,
                                              SOR = input$SOR,
                                              MQ = input$MQ,
                                              MQRankSum = input$MQRankSum,
                                              ReadPosRankSum = input$ReadPosRankSum,
                                              min_genotype_quality = input$min_genotype_quality,
                                              platform = input$platform,
                                              batch_size = input$batch_size,
                                              add_RGs = input$add_RGs,
                                              java_path = x$files$java_path,
                                              gatk4_path = x$files$gatk4_path,
                                              picard_path = x$files$picard_path,
                                              concatenate_final_vcfs = input$concatenate_final_vcfs)

      remove_modal_spinner()
      x$genotyped <- 1
    }

  })

  #===============generate final report===========
  # report on genotypes
  observeEvent(x$files_report$genotypes, ignoreNULL = TRUE, ignoreInit = TRUE,{
    hide("genotypes_files_report")
    output$genotypes_files_report <- renderText(paste0("Final Genotypes File: ", x$files_report$genotypes))
    if(length(x$files_report$genotypes) > 0){
      show("genotypes_files_report")
    }
  })

}

shinyApp(ui, server)
