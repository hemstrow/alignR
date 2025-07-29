.check_system_install <- function(code) as.logical(Sys.getenv(paste0(code, "_install")))


.parse_shinyFiles_path <- function(root, input_obj){
  if(!is.list(input_obj)){
    return(NULL)
  }
  
  out <- character(length(input_obj[[1]]))
  
  for(i in 1:length(input_obj[[1]])){
    path <- as.list(unlist(input_obj[[1]][[i]]))
    if(path[[1]] == ""){
      path[1] <- NULL
    }
    
    base_dir <- input_obj$root
    base_dir <- root[which(names(root) == base_dir)]
    out[i] <- normalizePath(do.call(file.path, c(base_dir, path)))
  }
  
  return(out)
}

.fastq_file_reporter <- function(files, rename_R3 = FALSE){
  report_fastq <- character()
  if("R1" %in% names(files)){
    report_fastq <- c(report_fastq, paste0("R1 fastq file: ", files$R1))
  }
  if("R2" %in% names(files)){
    report_fastq <- c(report_fastq, paste0("R2 fastq file: ", files$R2))
  }
  if("R3" %in% names(files)){
    if(rename_R3){
      report_fastq <- c(report_fastq, paste0("R2 fastq file: ", files$R3))
    }
    else{
      report_fastq <- c(report_fastq, paste0("R3 fastq file: ", files$R3))
    }
  }
  
  report_barcode <- character()
  if("barcode_file_1" %in% names(files)){
    report_barcode <- c(report_barcode, paste0("Barcode/index file 1:", files$barcode_file_1))
  }
  if("barcode_file_2" %in% names(files)){
    report_barcode <- c(report_barcode, paste0("Barcode/index file 2:", files$barcode_file_2))
  }
  
  report_demulti <- vector("list", length = 2)
  names(report_demulti) <- c("RA", "RB")
  if("demultiplexed_files" %in% names(files)){
    if("RA" %in% names(files$demultiplexed_files)){
      report_demulti$RA <- c(report_demulti$RA, paste0("RA fastq file: ", files$demultiplexed_files$RA))
    }
    if("RB" %in% names(files$demultiplexed_files)){
      report_demulti$RB <- c(report_demulti$RB, paste0("RB fastq file: ", files$demultiplexed_files$RB))
    }
  }
  if(length(report_demulti$RB) == 0){report_demulti$RB <- NULL}
  
  
  report_reference <- character()
  if("reference_genome" %in% names(files)){
    report_reference <- c(report_reference, paste0("Reference Genome: ", files$reference_genome))
  }
  
  bam_report <- character()
  if("alignments" %in% names(files)){
    bam_report <- c(bam_report, files$alignments)
  }
  
  report_paralog_reference <- character()
  if("paralog_reference" %in% names(files)){
    report_paralog_reference <- c(report_paralog_reference, files$paralog_reference)
  }
  
  report_paralog_pops <- character()
  if("filter_paralogs_pops" %in% names(files)){
    report_paralog_pops <- c( report_paralog_pops, files$filter_paralogs_pops)
  }
  
  report_sample_ids <- character()
  if("sample_ids" %in% names(files)){
    report_sample_ids <- c(report_sample_ids, files$sample_ids)
  }
  
  report_rf <- character()
  if("rf" %in% names(files)){
    report_rf <- c(report_rf, files$rf)
  }
  
  report_genotypes <- character()
  if("genotypes" %in% names(files)){
    report_genotypes <- c(report_genotypes, files$genotypes)
  }
  
  out <- list(fastq = report_fastq, barcode = report_barcode, 
              demulti = report_demulti, reference = report_reference,
              alignments = bam_report,
              paralog_reference = report_paralog_reference,
              sample_ids = report_sample_ids,
              filter_paralogs_pops = report_paralog_pops,
              rf = report_rf,
              genotypes = report_genotypes)
  
  return(out)
}

