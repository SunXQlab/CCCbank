#' @title Transform the database into a format suitable for Domino
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#'
#' @return db_path: The folder path used to save new database
#' 
#' @importFrom data.table fread
#' @importFrom dplyr distinct
#' @importFrom tidyr separate
#' @importFrom utils read.csv write.csv
#' 
#' @export
ChangeDominoDB <- function(priorDatabase = NULL){
  
  db_path <- file.path(getwd(), 'db')
  if(!dir.exists(db_path)){
    dir.create(db_path)
  }
  message(paste0('Path for saving new database: ', db_path))
  
  Ensg_path <- file.path(system.file("extdata", package = "CCCbank"), "ENSGtoSYMBOL.txt.gz")
  Uprot_path <- file.path(system.file("extdata", package = "CCCbank"), "UniprotToSymbol.tsv.gz")
  
  EnsgToSymbol <- read.csv(Ensg_path, header = TRUE)
  EnsgToSymbol <- EnsgToSymbol[, c('Gene.stable.ID', 'Gene.name')]
  EnsgToSymbol <- dplyr::distinct(EnsgToSymbol)
  
  UprotToSymbol <- data.table::fread(Uprot_path, sep = '\t')
  UprotToSymbol <- suppressWarnings(tidyr::separate(UprotToSymbol, `Gene Names`, 'Gene Names', '; '))
  UprotToSymbol <- suppressWarnings(tidyr::separate(UprotToSymbol, `Gene Names`, 'Gene Names', ' '))
  UprotToSymbol$Reviewed <- NULL
  UprotToSymbol$Organism <- NULL
  UprotToSymbol$Length <- NULL
  
  priorDatabase$ligand <- as.character(priorDatabase$ligand)
  priorDatabase$receptor <- as.character(priorDatabase$receptor)
  priorDatabase$lr <- as.character(priorDatabase$lr)
  db <- DominoDBFormat(priorDatabase, EnsgToSymbol, UprotToSymbol)
  
  write.csv(db$interaction, file = paste0(db_path, '/interactions.csv'), row.names = FALSE)
  write.csv(db$genes, file = paste0(db_path, '/genes.csv'), row.names = FALSE)
  write.csv(db$protein, file = paste0(db_path, '/proteins.csv'), row.names = FALSE)
  write.csv(db$complex, file = paste0(db_path, '/complexes.csv'), row.names = FALSE)
  
  return(db_path)
}

#' @title Transform the database into a format suitable for scSeqComm
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#'
#' @return dataframe
#' 
#' @export
ChangescSeqCommDB <- function(priorDatabase=NULL){
  db2 <- priorDatabase
  db2$ligand <- gsub('&', ',', db2$ligand)
  db2$receptor <- gsub('&', ',', db2$receptor)
  db2$lr <- paste(db2$ligand, db2$receptor, sep = '_')
  
  db2$lr <- NULL
  db.lr <- db2
  
  return(db.lr)
}

#' @title Transform the database into a format suitable for ICELLNET
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of ICELLNET and other databases.
#'
#' @return dataframe
#' 
#' @importFrom utils read.csv
#' @importFrom tidyr separate
#' 
#' @export
ChangeICELLNETDB <- function(priorDatabase=NULL, extension=TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  
  origin_db_path <- file.path(system.file("extdata", package = "CCCbank"), "ICELLNETdb.rds")
  origin_db <- readRDS(origin_db_path)
  
  if(extension){
    db1 <- origin_db[,1:5]
    db1$ligand <- paste(db1$`Ligand 1`, db1$`Ligand 2`, sep = '&')
    db1$receptor <- paste(db1$`Receptor 1`, db1$`Receptor 2`, db1$`Receptor 3`, sep = '&')
    db1$ligand <- gsub('&NA','', db1$ligand)
    db1$receptor <- gsub('&NA','', db1$receptor)
    db1$ligand <- gsub('& $','', db1$ligand)
    db1$receptor <- gsub('& $','', db1$receptor)
    
    dup.lrs <- intersect(paste(db1$ligand, db1$receptor, sep = '_'), priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- ICELLNETDBFormat(db2)
      db <- rbind(origin_db, db2)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
    }
  }else{
    db <- ICELLNETDBFormat(priorDatabase)
  }
  return(db)
}

#' @title Transform the database into a format suitable for RNAMagnet
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of RNAMagnet and other databases.
#' 
#' @return dataframe
#' 
#' @importFrom devtools install_github
#' @importFrom stringr str_to_title
#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct
#'
#' @export
ChangeRNAMagnetDB <- function(priorDatabase=NULL, extension=TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  
  if(!requireNamespace('RNAMagnet')){
    message("RNAMagnet hasn't been installed. Now installing...")
    devtools::install_github("veltenlab/rnamagnet")
  }
  
  origin_db <- RNAMagnet::ligandsReceptors_2.0.0
  
  db2 <- priorDatabase
  db2$ligand <- stringr::str_to_title(db2$ligand)
  db2$receptor <- stringr::str_to_title(db2$receptor)
  db2$lr <- paste(db2$ligand, db2$receptor, sep = '_')
  
  if(extension){
    db1 <- origin_db[, 2:3]
    db1 <- tidyr::separate_rows(db1, Ligand.Mouse, sep = '\\|')
    db1 <- tidyr::separate_rows(db1, Receptor.Mouse, sep = '\\|')
    db1 <- dplyr::distinct(db1)
    lr1 <- paste(db1$Ligand.Mouse, db1$Receptor.Mouse, sep = '_')
    
    dup.lrs <- intersect(lr1, db2$lr)
    db2 <- db2[!(db2$lr %in% dup.lrs), ]
    
    if(dim(db2)[[1]]!=0){
      
      db2 <- RNAMagnetDBFormat(db2)
      db <- rbind(origin_db, db2)
      
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
    }
  }else{
    db <- RNAMagnetDBFormat(db2)
  }
  return(db)
}

#' @title Transform the database into a format suitable for CellPhoneDB
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of CellPhoneDB and other databases.
#'
#' @return db_path: The folder path used to save new database
#' 
#' @import reticulate
#' @importFrom dplyr distinct
#' @importFrom tidyr separate
#' @importFrom stringr str_split
#' @importFrom data.table fread
#' @importFrom utils read.csv write.csv
#' 
#' @export
ChangeCellPhoneDB <- function(priorDatabase=NULL, extension=TRUE){
  if(!any(grepl('&', priorDatabase$lr))){
    extension <- TRUE
    message("The parameter 'extension' is set to TRUE, beacause the priorDatabases don't contain any multi-subunit of ligands or receptors. Otherwise, it will cause error!")
  }
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  
  db_path <- file.path(getwd(), 'db')
  if(!dir.exists(db_path)){
    dir.create(db_path)
  }
  message(paste0('Path for saving new database: ', db_path))

  origin_complex <- read.csv(file.path(system.file("extdata", package = "CCCbank"), "CellPhoneDB", 'complexes.csv'))
  origin_gene <- read.csv(file.path(system.file("extdata", package = "CCCbank"), "CellPhoneDB", 'genes.csv'))
  origin_interaction <- read.csv(file.path(system.file("extdata", package = "CCCbank"), "CellPhoneDB", 'interactions.csv'))
  origin_protein <- read.csv(file.path(system.file("extdata", package = "CCCbank"), "CellPhoneDB", 'proteins.csv'))
  
  EnsgToSymbol <- read.csv(file.path(system.file("extdata", package = "CCCbank"), 'ENSGtoSYMBOL.txt.gz'), header = TRUE)
  EnsgToSymbol <- EnsgToSymbol[, c('Gene.stable.ID', 'Gene.name')]
  EnsgToSymbol <- dplyr::distinct(EnsgToSymbol)
  
  UprotToSymbol <- data.table::fread(file.path(system.file("extdata", package = "CCCbank"), 'UniprotToSymbol.tsv.gz'), sep = '\t')
  UprotToSymbol <- suppressWarnings(tidyr::separate(UprotToSymbol, `Gene Names`, 'Gene Names', '; '))
  UprotToSymbol <- suppressWarnings(tidyr::separate(UprotToSymbol, `Gene Names`, 'Gene Names', ' '))
  UprotToSymbol$Reviewed <- NULL
  UprotToSymbol$Organism <- NULL
  UprotToSymbol$Length <- NULL
  
  # handle origin prior database
  if(T){
    
    db1_gene <- dplyr::distinct(origin_gene, gene_name, uniprot, hgnc_symbol, .keep_all = TRUE)
    db1_gene <- db1_gene[,1:2]
    
    db1_complex <- origin_complex[,1:5]
    db1_complex <- merge(db1_complex, db1_gene, by.x = "uniprot_1", by.y = "uniprot")
    db1_complex <- merge(db1_complex, db1_gene, by.x = "uniprot_2", by.y = "uniprot")
    db1_complex <- merge(db1_complex, db1_gene, by.x = "uniprot_3", by.y = "uniprot", all.x = TRUE)
    db1_complex <- suppressWarnings(merge(db1_complex, db1_gene, by.x = "uniprot_4", by.y = "uniprot", all.x = TRUE))
    db1_complex <- db1_complex[, -c(1:4)]
    colnames(db1_complex)[2:5] <- c('gene_1', 'gene_2', 'gene_3', 'gene_4')
    db1_complex$gene <- paste(db1_complex$gene_1, db1_complex$gene_2, db1_complex$gene_3, db1_complex$gene_4, sep = "&")
    db1_complex$gene <- gsub("&NA", "", db1_complex$gene)
    db1_complex <- db1_complex[,c("complex_name", "gene")]
    
    colnames(db1_complex) <- c('partner', 'gene')
    colnames(db1_gene) <- c('gene', 'partner')
    db1_gene <- db1_gene[, c('partner', 'gene')]
    db1_gene <- rbind(db1_gene, db1_complex)
    db1_gene <- dplyr::distinct(db1_gene)
    
    db1_interaction <- origin_interaction
    db1_interaction <- merge(db1_interaction, db1_complex, 
                             by.x = "partner_a", by.y = "partner", all.x = TRUE)
    db1_interaction$partner_a[which(!is.na(db1_interaction$gene))] <- db1_interaction$gene[which(!is.na(db1_interaction$gene))]
    db1_interaction$gene <- NULL
    db1_interaction <- merge(db1_interaction, db1_complex, 
                             by.x = "partner_b", by.y = "partner", all.x = TRUE)
    db1_interaction$partner_b[which(!is.na(db1_interaction$gene))] <- db1_interaction$gene[which(!is.na(db1_interaction$gene))]
    db1_interaction$gene <- NULL
    db1_protein <- origin_protein[, 1:2]
    db1_protein$protein_name <- gsub('_HUMAN', '', db1_protein$protein_name)
    db1_interaction <-  merge(db1_interaction, db1_protein, 
                              by.x = "partner_a", by.y = "uniprot", all.x = TRUE)
    db1_interaction$partner_a[which(!is.na(db1_interaction$protein_name))] <- db1_interaction$protein_name[which(!is.na(db1_interaction$protein_name))]
    db1_interaction$protein_name <- NULL
    db1_interaction <-  merge(db1_interaction, db1_protein, 
                              by.x = "partner_b", by.y = "uniprot", all.x = TRUE)
    db1_interaction$partner_b[which(!is.na(db1_interaction$protein_name))] <- db1_interaction$protein_name[which(!is.na(db1_interaction$protein_name))]
    db1_interaction$protein_name <- NULL
    db1_interaction <- db1_interaction[, 1:2]
    
  }
  
  if(extension){
    dup.lrs1 <- intersect(paste(db1_interaction$partner_a, db1_interaction$partner_b, sep = '_'), priorDatabase$lr)
    dup.lrs2 <- intersect(paste(db1_interaction$partner_b, db1_interaction$partner_a, sep = '_'), priorDatabase$lr) 
    dup.lrs <- unique(c(dup.lrs1, dup.lrs2))
    
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
  }else{
    db2 <- priorDatabase
  }
  
  if(dim(db2)[[1]]!=0){
    db2 <- CellPhoneDBFormat(db2, db1_gene, origin_protein, EnsgToSymbol, UprotToSymbol)
    write.csv(db2$interaction, file = paste0(db_path, '/interactions_add.csv'), row.names = FALSE)
    
    if(!is.null(db2$complex)){
      write.csv(db2$complex, file = paste0(db_path, '/complexes_add.csv'), row.names = FALSE)
    }
    
    if(!is.null(db2$genes)){
      write.csv(db2$genes, file = paste0(db_path, '/genes_add.csv'), row.names = FALSE)
    }
    
    if(!is.null(db2$protein)){
      write.csv(db2$protein, file = paste0(db_path, '/proteins_add.csv'), row.names = FALSE)
    }
    
    # generate new LR database in cpdb conda environment 
    command <- paste0('cellphonedb database generate ',
                      ifelse('interactions_add.csv' %in% list.files(db_path), 
                             paste0('--user-interactions ', paste0(db_path, '/interactions_add.csv ')), ' '),
                      ifelse('complexes_add.csv' %in% list.files(db_path), 
                             paste0('--user-complex ', paste0(db_path, '/complexes_add.csv ')), ' '),
                      ifelse('genes_add.csv' %in% list.files(db_path), 
                             paste0('--user-gene ', paste0(db_path, '/genes_add.csv ')), ' '),
                      ifelse('proteins_add.csv' %in% list.files(db_path), 
                             paste0('--user-protein ', paste0(db_path, '/proteins_add.csv ')), ' '), 
                      ifelse(extension, ' ', '--user-interactions-only '), 
                      '--result-path ', db_path)
    
    if(!('cpdb' %in% reticulate::conda_list()$name)){
      PyHome <- CreatCondaEnv(condaName='cpdb', python_version = '3.7', packages = 'cellphonedb')
    }
    condaName <- 'cpdb'
    
    condaHome <- reticulate::conda_binary()
    condaHome <- sub('\\/conda$', '', condaHome)
    
    
    message('Generating new database...')
    
    platform <- GetSysInfo()
    if(platform=='linux'){
      activateCommand <- paste0('source ', condaHome, '/activate ', condaName)
    }else if(platform=='windows'){
      activateCommand <- paste0('conda.bat activate ', condaName)
    }
    
    run_command <- paste(activateCommand, command,sep = ifelse(platform=='windows', ' & ', ' ; '))
    system(run_command)
    
    complex.path <- paste0(db_path, 'complex_input.csv')
    gene.path <- paste0(db_path, 'gene_input.csv')
    protein.path <- paste0(db_path, 'protein_input.csv')
  }else{
    message('The database to be added is included in the default database.\nUsing the Default database...')
    db_path <- file.path(system.file("extdata", package = "CCCbank"), "CellPhoneDB")
  }
  
  return(db_path)
}

#' @title Transform the database into a format suitable for iTALK
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of iTALK and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#'
#' @return dataframe
#' 
#' @importFrom devtools install_github
#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct
#' 
#' @export
ChangeiTALKDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  
  if(extension){
    
    if(!requireNamespace('iTALK')){
      message("iTALK hasn't been installed. Now installing...")
      devtools::install_github("Coolgenome/iTALK")
    }
    
    origin_db <- iTALK::database
    origin_db$Receptor.ApprovedSymbol[which(is.na(origin_db$Receptor.ApprovedSymbol))] <- origin_db$Receptor.Name[which(is.na(origin_db$Receptor.ApprovedSymbol))]
    lr1 <- paste(origin_db$Ligand.ApprovedSymbol, origin_db$Receptor.ApprovedSymbol, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- iTALKDBFormat(db2, keep_complexes)
      db <- rbind(origin_db, db2)
      db <- dplyr::distinct(db, Ligand.ApprovedSymbol, Receptor.ApprovedSymbol, .keep_all = TRUE)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
    }
    
  }else{
    db <- iTALKDBFormat(priorDatabase, keep_complexes)
  }
  return(db)
}

#' @title Transform the database into a format suitable for CellTalker
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns.
#' @param extension Logic. Whether use the combination of default database of CellTalker and other databases?
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out out.
#'
#' @return dataframe
#' 
#' @importFrom dplyr distinct
#' @importFrom tidyr separate_rows
#' 
#' @export
ChangeCellTalkerDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  
  if(extension){
    origin_db <- celltalker::ramilowski_pairs
    lr1 <- paste(origin_db$ligand, origin_db$receptor, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- CellTalkerDBFormat(db2, keep_complexes)
      db <- rbind(origin_db, db2)
      db <- dplyr::distinct(db, ligand, receptor, .keep_all = TRUE)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
    }
    
  }else{
    db <- CellTalkerDBFormat(priorDatabase, keep_complexes)
  }
  return(db)
}

#' @title Transform the database into a format suitable for SingleCellSignalR
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of SingleCellSignalR and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#'
#' @return dataframe
#' 
#' @importFrom dplyr distinct
#' @importFrom tidyr separate_rows
#' 
#' @export
ChangeSCSRDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  
  if(extension){
    if(!requireNamespace('SingleCellSignalR')){
      message("SingleCellSingalR hasn't been installed. Now installing...")
      if (!requireNamespace("BiocManager"))
      BiocManager::install("SingleCellSignalR")
    }
    db1 <- SingleCellSignalR::LRdb
    lr1 <- paste(db1$ligand, db1$receptor, sep = '_')
    dup.lrs <- intersect(lr1, priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- SCSRDBFormat(db2, keep_complexes)
      db <- rbind(db1, db2)
      db <- dplyr::distinct(db, ligand, receptor, .keep_all = TRUE)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- db1
    }
    
  }else{
    db <- SCSRDBFormat(priorDatabase, keep_complexes)
  }
  
  return(db)
}

#' @title Transform the database into a format suitable for NicheNet
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of NicheNet and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#'
#' @return dataframe
#' 
#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct
#' @importFrom utils download.file
#' 
#' @export
ChangeNicheNetDB <- function(priorDatabase = NULL, extension = TRUE, keep_complexes=TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  
  if(!require('nichenetr', quietly = TRUE)){
    message("nichenetr hasn't been installed. Now installing...")
    devtools::install_github("saeyslab/nichenetr")
    require('nichenetr', quietly = TRUE)
  }
  
  fpath <-  paste0(getwd(), '/NicheNetTemp')
  if(!dir.exists(fpath)){
    dir.create(fpath)
  }
  message(paste0('The files generated during the process are stored in ', fpath))
  
  if(!file.exists(paste0(fpath, '/sig_network.rds'))){
    message('Downloading sig_network from zendo...')
    download.file(url = "https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds", 
                  destfile = file.path(fpath, 'sig_network.rds'))
  }
  sig_network = readRDS(file.path(fpath, 'sig_network.rds'))
  
  if(!file.exists(paste0(fpath, '/gr_network.rds'))){
    message('Downloading gr_network from zendo...')
    download.file(url = "https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds", 
                  destfile = file.path(fpath, 'gr_network.rds'))
  }
  gr_network = readRDS(file.path(fpath, 'gr_network.rds'))

  lr_network <- readRDS(file.path(system.file("extdata", package = "CCCbank"), "NicheNet", 'lr_network.rds'))
  
  if(keep_complexes){
    db2 <- tidyr::separate_rows(priorDatabase, ligand, sep='&')
    db2 <- tidyr::separate_rows(db2, receptor, sep = '&')
    db2$lr <- paste(db2$ligand, db2$receptor, sep = '_')
    db2 <- dplyr::distinct(db2)
  }else{
    db2 <- priorDatabase[!grepl('&', priorDatabase$lr), ]
  }
  
  if(extension){
    diff.lr <- setdiff(db2$lr, paste(lr_network$from, lr_network$to, sep = '_'))
    if(length(diff.lr)>=1){
      db2 <- db2[(db2$lr %in% diff.lr), ]
      db2$lr <- NULL
      colnames(db2) <- c('from', 'to')
      db2$source <- 'user-defined'
      db2$database <- 'user-defined'
      db2 <- rbind(lr_network, db2)
      message('Constructing the ligand-target model...')
      db <- NicheNetDBFormat(db2, sig_network, gr_network)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      if(!file.exists(paste0(fpath, '/ligand_target_matrix.rds'))){
        message('Downloading ligand_target_matrix from zendo...')
        download.file(url = "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds", 
                      destfile = file.path(fpath, 'ligand_target_matrix.rds'))
      }
      ligand_target_matrix = readRDS(file.path(fpath, 'ligand_target_matrix.rds'))
      
      if(!file.exists(paste0(fpath, '/weighted_networks.rds'))){
        message('Downloading weighted_networks from zendo...')
        download.file(url = "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds", 
                      destfile = file.path(fpath, 'weighted_networks.rds'))
      }
      weighted_networks = readRDS(file.path(fpath, 'weighted_networks.rds'))
    }
  }else{
    db2$lr <- NULL
    colnames(db2) <- c('from', 'to')
    db2$source <- 'user-defined'
    db2$database <- 'user-defined'
    message('Constructing the ligand-target model...')
    db <- NicheNetDBFormat(db2, sig_network, gr_network)
  }
  
  if(exists('db')){
    ligand_target_matrix <- db$ligand_target_matrix
    weighted_networks <- db$weighted_networks
    lr_network <- db2
  }
  output <- list(ligand_target_matrix = ligand_target_matrix,
                 weighted_networks = weighted_networks,
                 lr_network = lr_network)
  saveRDS(output, file = paste0(fpath, '/new_db.rds'))
  return(output)
}

#' @title Transform the database into a format suitable for scMLnet
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of scMLnet and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#'
#' @return dataframe
#' 
#' @importFrom dplyr distinct
#' @importFrom tidyr separate_rows
#' @importFrom utils read.table
#' 
#' @export
ChangescMLnetDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  
  if(extension){
    origin_db <- read.table(file.path(system.file("extdata", package = "CCCbank"), "scMLnet", 'LigRec.txt'), header = T)
    colnames(origin_db)[2:3] <- c("source", "target")
    
    lr1 <- paste(origin_db$source, origin_db$target, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- scMLnetDBFormat(db2, keep_complexes)
      LigRecLib <- rbind(origin_db, db2)
      LigRecLib <- dplyr::distinct(LigRecLib, source, target, .keep_all = TRUE)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      LigRecLib <- origin_db
    }
    
  }else{
    LigRecLib <- scMLnetDBFormat(priorDatabase, keep_complexes)
  }
  return(LigRecLib)
}

#' @title Transform the database into a format suitable for NATMI
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of NATMI and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#' @param NATMI_dir Folder path of source code of NATMI. If NULL, it will download the source code of NATMI from github in the current directory.
#'
#' @return NATMI_dir Folder path of source code of NATMI
#' 
#' @importFrom dplyr distinct
#' @importFrom tidyr separate_rows
#' @importFrom utils read.csv write.csv
#' 
#' @export
ChangeNATMIDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE, NATMI_dir =NULL){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  
  if(is.null(NATMI_dir)){
    NATMI_dir <- Install_NATMI()
  }
  
  if(extension){
    origin_db <- read.csv(paste0(NATMI_dir, '/lrdbs/lrc2p.csv'))
    lr1 <- paste(origin_db$Ligand.gene.symbol, 
                 origin_db$Receptor.gene.symbol, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    colnames(origin_db) <- c('Ligand gene symbol', 'Receptor gene symbol')
    if(dim(db2)[[1]]!=0){
      db2 <- NATMIDBFormat(db2, keep_complexes)
      db <- rbind(origin_db, db2)
      db <- dplyr::distinct(db)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
    }
    
  }else{
    db <- NATMIDBFormat(priorDatabase, keep_complexes)
  }
  write.csv(db, file = paste0(NATMI_dir, '/lrdbs/lrdb.csv'), row.names = FALSE)
  return(NATMI_dir)
}

#' @title Transform the database into a format suitable for CellChat
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of CellChat and other databases
#' @param species human, mouse and zebrafish
#'
#' @return list
#' 
#' @importFrom devtools install_github
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr separate
#' 
#' @export
ChangeCellChatDB <- function(priorDatabase=NULL, extension=TRUE, species=NULL){
  #species: support human, mouse and zebrafish
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  species: ', species))
  
  if (!requireNamespace('CellChat')) {
    message("CellChat hasn't been installed. Now installing...")
    devtools::install_github("sqjin/CellChat")
  }
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- 'human'
    origin_db <- CellChat::CellChatDB.human
  }else if(species=='human'){
    origin_db <- CellChat::CellChatDB.human
  }else if(species == 'mouse'){
    origin_db <- CellChat::CellChatDB.mouse
  }else if(species == 'zebrafish'){
    origin_db <- CellChat::CellChatDB.zebrafish
  }else{
    stop('The species provided is not supported by CellChat!\nSupporting species: human, mouse, zebrafish.')
  }
  
  if(extension){
    # Get complexes names and its corresponding genes
    complexes <- origin_db$complex
    complexes$all <- paste(complexes$subunit_1, complexes$subunit_2, 
                           complexes$subunit_3, complexes$subunit_4, sep='_')
    complexes$all <- gsub('__$', '', complexes$all)
    complexes$all <- gsub('_$', '', complexes$all)
    complexes <- tibble::rownames_to_column(complexes, 'name')
    complexes <- complexes[, c('name', 'all')]
    
    # Recover genes of LR interactions
    db1 <- origin_db$interaction
    db1 <- db1[, 3:4]
    db1 <- merge(db1, complexes, by.x = 'ligand', by.y = 'name', all.x = TRUE)
    db1$ligand[!is.na(db1$all)] <- db1$all[!is.na(db1$all)]
    db1$all <- NULL
    db1 <- merge(db1, complexes, by.x = 'receptor', by.y = 'name', all.x = TRUE)
    db1$receptor[!is.na(db1$all)] <- db1$all[!is.na(db1$all)]
    db1$all <- NULL
    db1$receptor <- gsub('_', '&', db1$receptor)
    db1$ligand <- gsub('_', '&', db1$ligand)
    
    lr1 <- paste(db1$ligand, db1$receptor, sep = '_')
    dup.lrs <- intersect(lr1, priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    
    if(dim(db2)[[1]]!=0){
      db2 <- CellChatDBFormat(db2, species)
      
      interaction <- rbind(origin_db$interaction, db2$interaction)
      if(!is.null(db2$complex)){
        complex <- rbind(origin_db$complex, db2$complex)
      }else{
        complex <- origin_db$complex
      }
      if(!is.null(db2$geneInfo)){
        geneInfo <- rbind(origin_db$geneInfo, db2$geneInfo)
      }else{
        geneInfo <- origin_db$geneInfo
      }
      
    }else{
      
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      interaction <- origin_db$interaction
      complex <- origin_db$complex
      geneInfo <- origin_db$geneInfo
      
    }
  }else{
    db2 <- CellChatDBFormat(priorDatabase, species)
    interaction <- db2$interaction
    if(!is.null(db2$complex)){
      complex <- db2$complex
    }else{
      complex <- origin_db$complex
    }
    if(!is.null(db2$geneInfo)){
      geneInfo <- rbind(origin_db$geneInfo, db2$geneInfo)
    }else{
      geneInfo <- origin_db$geneInfo
    }
    
  }
  
  db <- list(interaction = interaction,
             complex = complex,
             cofactor = origin_db$cofactor,
             geneInfo = geneInfo)
  return(db)
}

#' @title Transform the database into a format suitable for Connectome
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with three columns: 'ligand', 'receptor' and 'lr'.
#' @param extension Logic. Whether use the combination of default database of Connectome and other databases.
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#' @param species human, mouse, rat and pig. If other species, extension will be set to FALSE
#'
#' @return dataframe
#' 
#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct
#' 
#' @export
ChangeConnectomeDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE, species = NULL){
  #species: support human, mouse, rat and pig (如果species不是这四个，可以用于任何物种，但extension会被设置为FALSE)
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  message(paste0('  species: ', species))
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- 'human'
    origin_db <- Connectome::ncomms8866_human
  }else if(species == 'human'){
    origin_db <- Connectome::ncomms8866_human
  }else if(species == 'mouse'){
    origin_db <- Connectome::ncomms8866_mouse
  }else if(species == 'pig'){
    origin_db <- Connectome::ncomms8866_pig
  }else if(species == 'rat'){
    origin_db <- Connectome::ncomms8866_rat
  }else{
    message('The species is not supported by the Connectome. The parameter extension is set to FALSE')
    origin_db <- NULL
    extension <- FALSE
  }
  
  if(extension){
    
    origin_db <- origin_db[, c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'mode')]
    lr1 <- paste(origin_db$Ligand.ApprovedSymbol, 
                 origin_db$Receptor.ApprovedSymbol, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- ConnectomeDBFormat(db2, keep_complexes)
      db <- rbind(origin_db, db2)
      db <- dplyr::distinct(db, Ligand.ApprovedSymbol, 
                     Receptor.ApprovedSymbol, .keep_all = TRUE)
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
    }
    
  }else{
    db <- ConnectomeDBFormat(priorDatabase, keep_complexes)
  }
  return(db)
}

#' @title Transform the database into a format suitable for Connectome
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of Connectome and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#' @param species human and mouse
#'
#' @return list
#' 
#' @importFrom dplyr distinct
#' @importFrom tidyr separate_rows
#' 
#' @export
ChangeCytoTalkDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE, species = NULL){
  #species: support human and mouse
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  message(paste0('  species: ', species))
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- 'human'
    origin_pcg = CytoTalk::pcg_human
    origin_db = CytoTalk::lrp_human
  }else if(species == 'human'){
    origin_pcg = CytoTalk::pcg_human
    origin_db = CytoTalk::lrp_human
  }else if(species == 'mouse'){
    origin_pcg = CytoTalk::pcg_mouse
    origin_db = CytoTalk::lrp_mouse
  }else{
    stop('The species provided is not supported by CytoTalk!\nSupporting species: human, mouse.')
  }
  
  if(extension){
    lr1 <- paste(origin_db$ligand, origin_db$receptor, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    if(dim(db2)[[1]]!=0){
      db2 <- CytoTalkDBFormat(db2, keep_complexes)
      lrp <- rbind(origin_db, db2)
      lrp <- dplyr::distinct(lrp, .keep_all = TRUE)
      pcg <- union(origin_pcg, unique(unlist(lrp)))
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      lrp <- origin_db
      pcg <- origin_pcg
    }
  }else{
    db2 <- priorDatabase
    lrp <- CytoTalkDBFormat(db2, keep_complexes)
    pcg <- union(origin_pcg, unique(unlist(lrp)))
  }
  
  db <- list(lrp=lrp, pcg=pcg)
  return(db)
}

#' @title Transform the database into a format suitable for Connectome
#' @author Jiaxin Luo
#' 
#' @param priorDatabase Dataframe with 'ligand', 'receptor' and 'lr' three columns
#' @param extension Logic. Whether use the combination of default database of Connectome and other databases
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.
#' @param species human and mouse
#' @param PyHome The execution file path of python, which has installed scConnect. If NULL, python installed in systems will be found.
#'
#' @return PyHome: The execution file path of python, which has installed scConnect
#' 
#' @import reticulate
#' @importFrom dplyr distinct
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate_rows
#' @importFrom utils read.csv write.table
#' 
#' @export
ChangescConnectDB <- function(priorDatabase=NULL, extension=TRUE, keep_complexes = TRUE, 
                              species = NULL, PyHome = NULL){
  
  message('Settings as follow:')
  message(paste0('  extension: ', extension))
  message(paste0('  keep_complexes: ', keep_complexes))
  message(paste0('  species: ', species))
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'hsapiens'")
    species <- 'hsapiens'
  }else if(species=='human'){
    species <- 'hsapiens'
  }else if(species=='mouse'){
    species <- 'mmusculus'
  }else{
    message("Confirm that 'species' refers to its scientific name, such as 'hsapiens' and 'mmusculus'.")
  }
  
  PyHome <- CheckscConnect(PyHome)
  db_path <- CopyDefaultDatabase(species, PyHome)
  
  # Handle origin prior database ####
  origin_db <- read.csv(paste0(db_path, 'interactions_origin.csv'), sep = ';')
  origin_ligands <- read.csv(paste0(db_path, paste0(species, '_origin'), "/ligands.csv"), sep = ',')
  origin_receptors <- read.csv(paste0(db_path, paste0(species, '_origin'), "/receptors.csv"), sep = ',')
  
  # Get genes of ligands ####
  ligands_db1 <- origin_ligands
  ligands_db1$ligand_symbol <- gsub("\\['", '', ligands_db1$preprogene)
  ligands_db1$ligand_symbol <- gsub("\\']", '', ligands_db1$ligand_symbol)
  ligands_db1$ligand_symbol <- gsub("', '", '|', ligands_db1$ligand_symbol)
  ligands_db1$ligand_symbol[which(ligands_db1$ligand == 'CCL19')] <- 'CCL19'
  ligands_db1$ligand_symbol[which(ligands_db1$ligand == 'Ccl21a')] <- 'CCL21'
  ligands_db1$ligand_symbol[which(ligands_db1$ligand == 'Ccl21b')] <- 'CCL21'
  ligands_db1$ligand_symbol[which(ligands_db1$ligand == 'agouti')] <- 'ASIP'
  ligands_db1$ligand_symbol[which(ligands_db1$ligand == 'relaxin')] <- 'RLN2|RLN1'
  
  # Get genes of receptors ####
  receptors_db1 <- origin_receptors
  receptors_db1$receptor_symbol <- gsub("\\['", '', receptors_db1$gene)
  receptors_db1$receptor_symbol <- gsub("\\']", '', receptors_db1$receptor_symbol)
  receptors_db1 <- dplyr::distinct(receptors_db1, receptor, receptor_symbol, .keep_all = TRUE)
  
  # Recover interactions with genes ####
  ligands_db1_2 <- ligands_db1[which(ligands_db1$ligand_type!='molecule'),]
  db1 <- origin_db[which(origin_db$ligand %in% ligands_db1_2$ligand & origin_db$target %in% receptors_db1$receptor), ]
  db1 <- merge(db1, ligands_db1_2, by = 'ligand')
  db1 <- merge(db1, receptors_db1, by.x = 'target', by.y = 'receptor')
  db1 <- db1[, c('ligand_symbol', 'receptor_symbol')]
  db1 <- tidyr::separate_rows(db1, ligand_symbol, sep = '\\|')
  db1 <- dplyr::distinct(db1)
  db1$lr <- paste(db1$ligand_symbol, db1$receptor_symbol, sep = '_')
  
  if(extension){
    
    lr1 <- paste(db1$ligand_symbol, db1$receptor_symbol, sep = '_')
    dup.lrs <- intersect(lr1,  priorDatabase$lr)
    db2 <- priorDatabase[!(priorDatabase$lr %in% dup.lrs), ]
    
    if(dim(db2)[[1]]!=0){
      db2 <- scConnectDBFormat(db2, keep_complexes, ligands_db1, receptors_db1)
      db <- rbind(origin_db, db2$db)
      db <- dplyr::distinct(db, .keep_all = TRUE)
      ligands <- db2$ligands
      ligands$ligand_symbol <- NULL
      receptors <- db2$receptors
      receptors$receptor_symbol <- NULL
    }else{
      message('The database to be added is included in the default database. \nUsing the Default database... ')
      db <- origin_db
      ligands <- origin_ligands
      ligands$X <- NULL
      receptors <- origin_receptors
      receptors$X <- NULL
    }
    
  }else{
    db2 <- scConnectDBFormat(priorDatabase, keep_complexes, ligands_db1, receptors_db1)
    db <- db2$db
    ligands <- db2$ligands
    ligands$ligand_symbol <- NULL
    receptors <- db2$receptors
    receptors$receptor_symbol <- NULL
  }
  
  db_species_fpath <- paste0(db_path, species)
  if (!dir.exists(db_species_fpath)) {
    dir.create(db_species_fpath)
  }
  
  # overwrite
  rownames(db) <- NULL
  write.table(db, file = paste0(db_path, 'interactions.csv'), sep = ';', row.names = FALSE)
  rownames(ligands) <- NULL
  ligands <- tibble::rownames_to_column(ligands, 'X')
  write.table(ligands, file = paste0(db_species_fpath, '/ligands.csv'), sep = ',', row.names = FALSE, quote = TRUE)
  rownames(receptors) <- NULL
  receptors <- tibble::rownames_to_column(receptors, 'X')
  write.table(receptors, file = paste0(db_species_fpath, '/receptors.csv'), sep = ',', row.names = FALSE, quote = TRUE)
  # copy
  pep_fpath_origin <- paste0(db_path, paste0(species, '_origin'), "/peptide_ligands.csv")
  pep_fpath_new <- paste0(db_species_fpath, "/peptide_ligands.csv")
  mol_fpath_origin <- paste0(db_path, paste0(species, '_origin'), "/molecule_ligands.csv")
  mol_fpath_new <- paste0(db_species_fpath, "/molecule_ligands.csv")
  
  commands <- paste('cp -r ', pep_fpath_origin, ' ', pep_fpath_new)
  system(commands)
  commands <- paste('cp -r ', mol_fpath_origin, ' ', mol_fpath_new)
  system(commands)
  
  return(PyHome)
}