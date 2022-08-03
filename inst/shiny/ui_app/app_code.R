.parse_shinyFiles_path <- function(root, input_obj){
  path <- as.list(unlist(input_obj[[1]]))
  if(path[[1]] == ""){
    path[1] <- NULL
  }
  
  base_dir <- strsplit(input_obj$root, "\\.")[[1]][2]
  base_dir <- root[which(names(root) == base_dir)]
  return(normalizePath(do.call(file.path, c(base_dir, path))))
}