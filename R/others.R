CheckSeuratObject <- function(ser){
  if(!('celltype' %in% colnames(ser@meta.data))){
    stop("Confirm that cell type metadata has been stored in the 'celltype' column of ser@meta.data!")
  }
  
  if(DefaultAssay(ser) != 'RNA'){
    stop("Confirm that the scRNA-seq data is saved in the 'RNA' assay of ser!")
  }
}

CreatCondaEnv <- function(condaName, python_version, packages){
  
  message(paste0('Creating an conda environment called ', condaName, '...'))
  reticulate::conda_create(envname = condaName, python_version = python_version)
  
  for(pkg in packages){
    message('Installing packages: ', packages, '...')
    reticulate::conda_install(envname = condaName, pkg, pip = TRUE)
  }
  message('Done!')
  python_path <- reticulate::conda_python(envname = condaName)
  return(python_path)
}

GetSysInfo <- function(){
  sys_info <- Sys.info()
  if(sys_info['sysname'] %in% c('Linux', 'linux')){
    platform <- 'linux'
    message(paste0('Running on ', platform, '...'))
  }else if(sys_info['sysname'] %in% c('Windows', 'windows')){
    platform <- 'windows'
    message(paste0('Running on ', platform, '...'))
  }else{
    platform <- sys_info['sysname']
    warning("Running on an unknown operating system. May fail to run.")
  }
  return(platform)
}

Install_NATMI <- function(){
  message('Download NATMI from github...')
  command <- 'git clone https://github.com/asrhou/NATMI.git'
  system(command)
  
  NATMI_dir <- file.path(getwd(), 'NATMI')
  
  if(!file.exists(NATMI_dir)){
    stop('Automatic installation failed, please install it by yourself!')
  }else{
    message('NATMI has been installed!')
    message(paste0('The dir of NATMI: '), NATMI_dir)
    
    #copy_command <- paste('cp -r', paste0(NATMI_dir, '/lrdbs/lrc2p.csv'),
    #                      paste0(NATMI_dir, '/lrdbs/lrc2p_origin.csv'), sep = ' ')
    #system(copy_command)
    #message('The original database lrdbs/lrc2p.csv has been copied to lrc2p_origin.csv!')
    
    return(NATMI_dir)
  }
}

CheckscConnect <- function(PyHome=NULL){
  
  if(is.null(PyHome)){
    PyHome <- Sys.which("python")
    PyHome <- ifelse(PyHome=='', Sys.which("python3"), PyHome)
    if(PyHome==''){
      stop("Can't find python2/3 in the system!")
    }
  }
  message(paste0('Current PyHome is ', PyHome))
  
  reticulate::use_python(PyHome)
  library_exists <- reticulate::py_module_available('scConnect')
  
  if(!library_exists){
    message("scConnect hasn't been installed! Now installing...")
    reticulate::py_install('scConnect', pip = TRUE)
  }
  
  return(PyHome)
}

GetDbPathOfscConnect <- function(PyHome){
  platform <- GetSysInfo()
  #Get library of python
  if(platform == 'linux'){
    lib_path <- gsub('bin\\/python.*', '', PyHome)
    lib_path <- paste0(lib_path, 'lib/')
    lib_path <- dir(lib_path, pattern = "python[0-9]\\.[0-9]?$", full.names = TRUE)[1]
    lib_path <- file.path(lib_path, 'site-packages')
  }else if(platform == 'windows'){
    lib_path <- gsub('python[3]?\\.exe$', '', PyHome)
    lib_path <- paste0(lib_path, '/Lib/site-packages')
  }
  
  db_path <- paste0(lib_path, '/scConnect/data/Gene_annotation/2020-5/')
  return(db_path)
}

CopyDefaultDatabase <- function(species, PyHome){
  
  #species: support human, mouse and more
  if(is.null(species)){
    warning("No 'species' provied, 'species' is set to 'hsapiens'")
    species <- 'hsapiens'
  }else if(species=='human'){
    species <- 'hsapiens'
  }else if(species=='mouse'){
    species <- 'mmusculus'
  }else{
    message("Confirm that 'species' refers to its scientific name, such as 'hsapiens' and 'mmusculus'.")
  }
  
  db_path <- GetDbPathOfscConnect(PyHome)
  
  interactions_fpath <- file.path(db_path, 'interactions.csv')
  interactions_copy_fpath <- sub('interactions.csv', 'interactions_origin.csv', interactions_fpath)
  if(!file.exists(interactions_copy_fpath)){
    copy_command <- paste('cp -r', interactions_fpath, interactions_copy_fpath, sep = ' ')
    system(copy_command)
    message(paste0(interactions_fpath, ' has been copied to ', interactions_copy_fpath))
  }
  
  
  dir_path <- file.path(db_path, species)
  dir_copy_path <- sub(species, paste0(species, '_origin'), dir_path)
  if(!file.exists(dir_copy_path)){
    copy_command <- paste('cp -r', dir_path, dir_copy_path, sep = ' ')
    system(copy_command)
    message(paste0(dir_path, ' has been copied to ', dir_copy_path))
  }
  return(db_path)
}

RecoverDefaultDatabase <- function(species, PyHome){
  #species: support human, mouse and more
  if(is.null(species)){
    warning("No 'species' provied, 'species' is set to 'hsapiens'")
    species <- 'hsapiens'
  }else if(species=='human'){
    species <- 'hsapiens'
  }else if(species=='mouse'){
    species <- 'mmusculus'
  }else{
    message("Confirm that 'species' refers to its scientific name, such as 'hsapiens' and 'mmusculus'.")
  }
  
  db_path <- GetDbPathOfscConnect(PyHome)
  
  interactions_fpath <- file.path(db_path, 'interactions.csv')
  interactions_copy_fpath <- sub('interactions.csv', 'interactions_origin.csv', interactions_fpath)
  if(file.exists(interactions_copy_fpath)){
    copy_command <- paste('cp -rf', interactions_copy_fpath, interactions_fpath, sep = ' ')
    system(copy_command)
    #message(paste0(interactions_fpath, ' has been copied to ', interactions_copy_fpath))
  }
  
  dir_path <- file.path(db_path, species)
  dir_copy_path <- sub(species, paste0(species, '_origin'), dir_path)
  if(file.exists(dir_copy_path)){
    system(paste0('rm -rf ', dir_path))
    copy_command <- paste('cp -rf', dir_copy_path, dir_path, sep = ' ')
    system(copy_command)
    #message(paste0(dir_path, ' has been copied to ', dir_copy_path))
  }
  return(NULL)
}
#' @importFrom utils read.table
ProcssCPDBResult <- function(cpdb.complex, file.path){
  
  result <- read.table(file.path, header = TRUE, sep = "\t")
  result <- result[, c(2, 13:dim(result)[2])]
  
  result <- result %>% pivot_longer(cols = -interacting_pair, names_to = "sr", values_to = "LRscore")
  result <- dplyr::filter(result,  !is.na(LRscore))
  result <- tidyr::separate(data = result, col = sr, into = c("Sender", "Receiver"), sep = "\\.")
  
  complexes <- cpdb.complex$complex_name[which(stringr::str_count(cpdb.complex$complex_name, pattern = '_')==1)]
  for (complex in complexes) {
    change.pair <- which(grepl(complex, result$interacting_pair))
    if(length(change.pair)>0){
      change.complex <- gsub('_', '*', complex)
      result$interacting_pair <- gsub(complex, change.complex, result$interacting_pair)
    }
  }
  
  result <- tidyr::separate(data = result, col = interacting_pair, into = c("Ligand", "Receptor"), sep = "_")
  result$Ligand <- gsub('\\*', '_', result$Ligand)
  result$Receptor <- gsub('\\*', '_', result$Receptor)
  result <- merge(result, cpdb.complex, by.x = "Ligand", by.y = "complex_name", all.x = TRUE)
  result$Ligand[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
  result <- result[,-6]
  result <- merge(result, cpdb.complex, by.x = "Receptor", by.y = "complex_name", all.x = TRUE)
  result$Receptor[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
  result <- result[,-6]
  
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @importFrom utils read.csv
RecoverCpdbComplex <- function(gene.path, protein.path, complex.path){
  
  cpdb.complex <-read.csv(complex.path, na.strings = c('nan'))
  
  rec.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor == 'True')]
  lig.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor != 'True')]
  cpdb.complex <- cpdb.complex[,1:5]
  
  cpdb.gene <- read.csv(gene.path)
  cpdb.gene <- cpdb.gene[,1:2]
  cpdb.gene <- dplyr::distinct(cpdb.gene, .keep_all = TRUE)
  
  cpdb.protein <- read.csv(protein.path)
  cpdb.protein <- cpdb.protein[1:2]
  cpdb.protein <- dplyr::distinct(cpdb.protein, .keep_all = TRUE)
  cpdb.protein$protein_name <- gsub('_HUMAN', '', cpdb.protein$protein_name)
  
  tmp.cpdb.complex <- merge(cpdb.gene, cpdb.complex, by.y = "uniprot_1", by.x = "uniprot", all.y = TRUE)
  tmp.cpdb.complex$gene_name[is.na(tmp.cpdb.complex$gene_name)] <- 
    cpdb.protein$protein_name[match(tmp.cpdb.complex$uniprot[is.na(tmp.cpdb.complex$gene_name)], cpdb.protein$uniprot)]
  tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
  colnames(tmp.cpdb.complex)[1] <- "gene_1"
  
  tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_2", by.x = "uniprot", all.y = TRUE)
  tmp.cpdb.complex$gene_name[is.na(tmp.cpdb.complex$gene_name)] <- 
    cpdb.protein$protein_name[match(tmp.cpdb.complex$uniprot[is.na(tmp.cpdb.complex$gene_name)], cpdb.protein$uniprot)]
  tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
  colnames(tmp.cpdb.complex)[1] <- "gene_2"
  
  tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_3", by.x = "uniprot", all.y = TRUE)
  tmp.cpdb.complex$gene_name[is.na(tmp.cpdb.complex$gene_name)] <- 
    cpdb.protein$protein_name[match(tmp.cpdb.complex$uniprot[is.na(tmp.cpdb.complex$gene_name)], cpdb.protein$uniprot)]
  tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
  colnames(tmp.cpdb.complex)[1] <- "gene_3"
  
  tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_4", by.x = "uniprot", all.y = TRUE)
  tmp.cpdb.complex$gene_name[is.na(tmp.cpdb.complex$gene_name)] <- 
    cpdb.protein$protein_name[match(tmp.cpdb.complex$uniprot[is.na(tmp.cpdb.complex$gene_name)], cpdb.protein$uniprot)]
  tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
  colnames(tmp.cpdb.complex)[1] <- "gene_4"
  
  cpdb.complex <- tmp.cpdb.complex[,5:1]
  #rm(cpdb.gene, tmp.cpdb.complex, cpdb.protein);gc()
  
  cpdb.complex$gene <- paste(cpdb.complex$gene_1, cpdb.complex$gene_2, 
                             cpdb.complex$gene_3, cpdb.complex$gene_4,
                             sep = "&")
  cpdb.complex$gene <- gsub("&NA", "", cpdb.complex$gene)
  cpdb.complex <- cpdb.complex[,c("complex_name", "gene")]
  return(cpdb.complex)
}