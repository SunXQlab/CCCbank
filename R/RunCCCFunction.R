#' @title Run iTALK method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run iTALK method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param ... other parameters in \code{\link[iTALK]{rawParse}}, except for the following parameter: data
#' 
#' @seealso \code{\link[iTALK]{rawParse}}
#'
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @importFrom devtools install_github
#' @importFrom dplyr distinct
#' @import Seurat
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RuniTALK(ser)
#' 
#' # Run with other LR prior database
#' db <- ChangeiTALKDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
#' result <- RuniTALK(ser, priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeiTALKDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
#' result <- RuniTALK(ser, priorDatabase = db)
#' }
RuniTALK <- function(ser, priorDatabase = NULL, ...){
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('iTALK')){
    message("iTALK hasn't been installed. Now installing...")
    devtools::install_github("Coolgenome/iTALK")
  }
  
  if(!is.null(priorDatabase)){
    database <- priorDatabase
  }else{
    database <- iTALK::database
  }
  
  # Get normalized data and cell type metadata
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  matrix.sc <- as.data.frame(t(matrix.sc))
  matrix.sc$cell_type <- ser$celltype
  
  # Run
  # find top n percent highly expressed genes
  highly_exprs_genes <- iTALK::rawParse(matrix.sc,...)
  # find the ligand-receptor pairs from highly expressed genes
  comm.list<-c('growth factor','other','cytokine','checkpoint')
  
  result <- NULL
  for(comm.type in comm.list){
    res.tmp <- iTALK::FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm.type, database = database)
    res.tmp <- res.tmp[order(res.tmp$cell_from_mean_exprs*res.tmp$cell_to_mean_exprs,decreasing=T),]
    result <- rbind(result,res.tmp)
  }
  
  result$LRscore <- result$cell_from_mean_exprs*result$cell_to_mean_exprs
  colnames(result)[c(1:2, 4,6,8)] <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run CellTalker method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run CellTalker method to infer CCC.
#' 
#' @param ser Seurat object, contains raw and normalized data stored in 'counts' and 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param ... Other parameters in \code{\link[celltalker]{celltalk}} function, except for the following parameters: input_object, metadata_grouping, ligand_receptor_pairs
#' 
#' @seealso \code{\link[celltalker]{celltalk}}
#'
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import dplyr
#' @importFrom devtools install_github
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunCellTalker(ser)
#' 
#' # Run with other LR prior database
#' db <- ChangeCellTalkerDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
#' result <- RunCellTalker(ser, priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeCellTalkerDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
#' result <- RunCellTalker(ser, priorDatabase = db)
#' }
RunCellTalker <- function(ser, priorDatabase = NULL,...){
  
  set.seed(123)
  
  if(!requireNamespace('celltalker')){
    message("celltalker hasn't been installed. Now installing...")
    devtools::install_github("arc85/celltalker")
  }
  
  if(!is.null(priorDatabase)){
    database <- priorDatabase
  }else{
    database <- celltalker::ramilowski_pairs
  }
  
  # Run
  result <- celltalker::celltalk(input_object=ser,
                                 metadata_grouping='celltype',
                                 ligand_receptor_pairs=database,
                                 ...)
  result <- result %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    filter(fdr < 0.05) %>%
    filter(p_val < 0.05) %>%
    filter(interact_ratio > 0)
  
  result <- tidyr::separate(result, 'interaction', c('Ligand', 'Receptor'), sep = '_')
  result <- tidyr::separate(result, 'interaction_pairs', c('Sender', 'Receiver'), sep = '_')
  colnames(result)[5] <- 'LRscore'
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run RNAMagnet method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run RNAMagnet method to infer CCC.
#' 
#' @param ser Seurat object, contains raw data stored in 'counts' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param .version see details in \code{\link[RNAMagnet]{getLigandsReceptors}}
#' @param .cellularCompartment see details in \code{\link[RNAMagnet]{RNAMagnetBase}}
#' @param ... other parameters in \code{\link[RNAMagnet]{RNAMagnetBase}}, except for the following parameters: seurat, .version, .cellularCompartment
#' 
#' @seealso \code{\link[RNAMagnet]{RNAMagnetBase}}
#' @seealso \code{\link[RNAMagnet]{getLigandsReceptors}}
#'
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom devtools install_github
#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct
#' @importFrom stats aggregate
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunRNAMagnet(ser)
#' 
#' # Run with other LR prior database
#' db <- ChangeRNAMagnetDB(priorDatabase=priorDatabase, extension=FALSE)
#' result <- RunRNAMagnet(ser, priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeRNAMagnetDB(priorDatabase=priorDatabase, extension=TRUE)
#' result <- RunRNAMagnet(ser, priorDatabase = db)
#' }
RunRNAMagnet <- function(ser, priorDatabase = NULL, .version = '2.0.0', .cellularCompartment = c("Membrane", "ECM", "Both", "Secreted"), ...){
  #source('./CCCbank/inst/R_scripts/RNAMagnetSignaling_code.R')
  source(system.file("R_scripts", "RNAMagnetSignaling_code.R", package = "CCCbank"))
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('RNAMagnet')){
    message("RNAMagnet hasn't been installed. Now installing...")
    devtools::install_github("veltenlab/rnamagnet")
  }
  
  if(!is.null(priorDatabase)){
    db <- priorDatabase
    #assign db to 'ligandsReceptors_2.0.0' object
    assign("ligandsReceptors_2.0.0", db, envir = .GlobalEnv)
    message("The parameter '.version' is set to '2.0.0', because other database is used.")
    .version <- '2.0.0'
  }else{
    db <- switch(.version, latest = RNAMagnet::ligandsReceptors_2.2.0, 
                  stable =  RNAMagnet::ligandsReceptors_1.0.0, 
                 `1.0.0` =  RNAMagnet::ligandsReceptors_1.0.0,
                 `2.0.0` =  RNAMagnet::ligandsReceptors_2.0.0, 
                 `2.1.0` =  RNAMagnet::ligandsReceptors_2.1.0, 
                 `2.2.0` =  RNAMagnet::ligandsReceptors_2.2.0, 
                 `3.0.0` =  RNAMagnet::ligandsReceptors_3.0.0)
    if(.version == 'latest'){
      .version <- '2.2.0'
    }else if(.version == 'stable'){
      .version <- '1.0.0'
    }
    assign(paste0("ligandsReceptors_", .version), db, envir = .GlobalEnv)
  }
  
  result <- RNAMagnetSignaling(ser, .version = .version,
                               .cellularCompartment = .cellularCompartment, ...)
  
  temp <- c()
  for (sender in unique(ser$celltype)) {
    for (receiver in unique(ser$celltype)) {
      tmp <- RNAMagnet::getRNAMagnetGenes(result, sender, receiver)
      if(dim(tmp)[[1]]!=0){
        tmp$Sender <- sender
        tmp$Receiver <- receiver
        temp <- rbind(temp,tmp)
      }
    }
  }
  
  db <- db[, 1:3]
  db <- dplyr::distinct(db)
  result <- temp
  result <- merge(result, db, by.x = 'pair', by.y = 'Pair.Name')
  result$pair <- NULL
  colnames(result)[4:5] <- c('Ligand', 'Receptor')
  
  rownames(result) <- NULL
  result <- tidyr::separate_rows(result, Ligand, sep = '\\|')
  result <- tidyr::separate_rows(result, Receptor, sep = '\\|')
  result <- aggregate(score~Ligand+Receptor+Sender+Receiver, mean, data=result)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  colnames(result)[5] <- 'LRscore'
  
  return(result)
}

#' @title Run ICELLNET method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run ICELLNET method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param filter.perc see details in \code{\link[icellnet]{sc.data.cleaning}}
#' 
#' @seealso \code{\link[icellnet]{sc.data.cleaning}}
#'
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import dplyr
#' @importFrom devtools install_github
#' @importFrom tidyr separate pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom utils read.csv
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunICELLNET(ser)
#' 
#' # Run with other LR prior database
#' db <- ChangeICELLNETDB(priorDatabase=priorDatabase, extension=FALSE)
#' result <- RunICELLNET(ser, priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeICELLNETDB(priorDatabase=priorDatabase, extension=TRUE)
#' result <- RunICELLNET(ser, priorDatabase = db)
#' }
RunICELLNET <- function(ser, priorDatabase = NULL, filter.perc = NULL){

  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('icellnet')){
    message("icellnet hasn't been installed. Now installing...")
    devtools::install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")
  }
  
  if(!is.null(priorDatabase)){
    db <- priorDatabase
  }else{
    db_path <- file.path(system.file("extdata", package = "CCCbank"), "ICELLNETdb.rds")
    db <- readRDS(db_path)
  }
  
  ## Retrieve gene expression matrix
  # Taking into account the total nb of cells in each cluster
  # filter.perc=5: filter out the percent of genes less than 5%
  if(is.null(filter.perc)){
    filter.perc <- 5
    message(paste0('filter.perc in icellnet::sc.data.cleaning is set to ', filter.perc))
  }
  average.clean= icellnet::sc.data.cleaning(object = ser, db = db, filter.perc = filter.perc, save_file = F)
  ## Apply icellnet pipeline on cluster of interest
  data.icell=as.data.frame(icellnet::gene.scaling(as.data.frame(average.clean), n=1, db=db))
  
  PC.data=as.data.frame(data.icell[, colnames(data.icell)], row.names = rownames(data.icell))
  PC.target=data.frame(Class = colnames(PC.data)[-dim(data.icell)[2]], 
                       ID = colnames(PC.data)[-dim(data.icell)[2]], 
                       Cell_type = colnames(PC.data)[-dim(data.icell)[2]])
  rownames(PC.target) = colnames(PC.data)[-dim(data.icell)[2]]
  
  PC.ct <- colnames(PC.data)[-dim(data.icell)[2]]
  CC.ct <- colnames(PC.data)[-dim(data.icell)[2]]
  
  result.all <- lapply(CC.ct, function(ct){
    ## Compute intercellular communication scores
    score.computation = suppressWarnings(icellnet::icellnet.score(direction = "out", PC.data = PC.data,
                                                                  CC.data = as.data.frame(data.icell[,ct], row.names = rownames(data.icell)),
                                                                  PC.target = PC.target, PC = PC.ct[which(PC.ct!=ct)], CC.type = "RNAseq",
                                                                  PC.type = "RNAseq",  db = db))
    lr <- as.matrix(score.computation[[2]][apply(score.computation[[2]], 1, function(y) any(!is.na(y))),])
    lr <- as.matrix(lr[which(rowSums(lr) > 0),])
    lr
  })
  
  names(result.all) <- colnames(PC.data)[-dim(data.icell)[2]]
  
  result <- lapply(names(result.all), function(ct){
    lr <- result.all[[ct]]
    lr <- as.data.frame(lr)
    lr <- tibble::rownames_to_column(lr, "LR")
    colnames(lr) <- c("LR", paste(ct, colnames(lr)[2:dim(lr)[2]], sep = "_"))
    result.lr <- lr %>% tidyr::pivot_longer(cols = -LR, names_to = "sr", values_to = "LRscore")
    result.lr <- tidyr::separate(data = result.lr, col = sr, into = c("Sender", "Receiver"), sep = "_")
    result.lr <- tidyr::separate(data = result.lr, col = LR, into = c("Ligand", "Receptor"), sep = " / ")
    result.lr <- result.lr[which(result.lr$LRscore>0), ]
    result.lr$Ligand <- gsub(' \\+ ', '&', result.lr$Ligand)
    result.lr$Receptor <- gsub(' \\+ ', '&', result.lr$Receptor)
    result.lr$all <- paste(result.lr$Sender, result.lr$Ligand, result.lr$Receiver, result.lr$Receptor, sep = '_')
    result.lr <- dplyr::distinct(result.lr, all, .keep_all = TRUE)
  })
  result <- do.call(rbind, result)
  
  return(result)
}

#' @title Run NicheNet method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run NicheNet method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param lr If TRUE, result is the dataframe of L-R links.Otherwise, result is the dataframe of L-Target links. Defaults to TRUE
#' @param n_ligands Number of top-ranked ligands, which are further used to predict L-R links or L-Target links. Defaults to 20
#' @param n_targets Number of top-predicted target genes of ligands. If lr is FALSE, the parameter is required. Defaults to 200
#'
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms, if lr is TRUE. Otherwise, it contains  'Ligand', 'Target', 'Sender', 'Receiver' and 'LRscore' five colunms
#' 
#' @export
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @importFrom devtools install_github
#' @importFrom utils download.file
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunNicheNet(ser)
#' 
#' # Run with other LR prior database
#' db <- ChangeNicheNetDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
#' result <- RunNicheNet(ser, priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeNicheNetDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
#' result <- RunNicheNet(ser, priorDatabase = db)
#' }
RunNicheNet <- function(ser, priorDatabase = NULL, lr = TRUE, n_ligands = 20, n_targets = 200){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
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
  
  if(!is.null(priorDatabase)){
    ligand_target_matrix <- priorDatabase$ligand_target_matrix
    weighted_networks <- priorDatabase$weighted_networks
    lr_network <- priorDatabase$lr_network
  }else{
    lr_network <- readRDS(file.path(system.file("extdata", package = "CCCbank"), "NicheNet", 'lr_network.rds'))
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
  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  
  # find DEGs in all celltypes
  DE_table_receiver_all = FindAllMarkers(object = ser, min.pct = 0.05, 
                                         logfc.threshold = 0.15, test.use = "t", 
                                         return.thresh = 0.05)
  
  sender_ct <- unique(ser$celltype)
  
  receiver_ct <- DE_table_receiver_all[which(DE_table_receiver_all$p_val_adj<=0.05), ]$cluster %>%
    unique(.) %>% as.character(.)
  
  message(paste0("Getting top ", n_ligands, ' ligands as best_upstream_ligands!'))
  if(!lr){
    message(paste0("Getting top ", n_targets, ' targets as target genes of ligands!'))
  }
  
  result <- lapply(sender_ct, function(sender){
    tmp.result <- list()
    for (receiver in receiver_ct) {#[which(receiver_ct != sender)]
      # define the sender and receiver celltypes
      expressed_genes_receiver = nichenetr::get_expressed_genes(receiver, ser, pct = 0.05)
      background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
      
      expressed_genes_sender = nichenetr::get_expressed_genes(sender, ser, pct = 0.05)
      
      # find DEGs in receiver cells
      DE_table_receiver <- DE_table_receiver_all[which(DE_table_receiver_all$cluster == receiver), ]
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05) %>% pull(gene)
      geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
      
      # define the potential ligands
      ligands = lr_network %>% pull(from) %>% unique()
      receptors = lr_network %>% pull(to) %>% unique()
      
      expressed_ligands = intersect(ligands,expressed_genes_sender)
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      
      potential_ligands = lr_network %>% 
        filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
        pull(from) %>% 
        unique()
      if(length(potential_ligands)==0){
        next
      }
      # Perform NicheNet ligand activity analysis
      ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, 
                                                               background_expressed_genes = background_expressed_genes,
                                                               ligand_target_matrix = ligand_target_matrix, 
                                                               potential_ligands = potential_ligands)
      ligand_activities = ligand_activities %>% 
        arrange(-pearson) %>% 
        mutate(rank = rank(desc(pearson)))
      
      best_upstream_ligands = ligand_activities %>% top_n(n_ligands, pearson) %>% 
        arrange(-pearson) %>% pull(test_ligand) %>% unique()
      
      cp <- paste(sender, receiver, sep = "_")
      
      if(lr){
        
        ## Receptors of top-ranked ligands
        lr_network_top = lr_network %>% 
          filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
          distinct(from,to)
        
        best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
        
        lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
        colnames(lr_network_top_df_large) <- c('Ligand', 'Receptor', 'LRscore')
        lr_network_top_df_large$Sender <- sender
        lr_network_top_df_large$Receiver <- receiver
        tmp.result[[cp]] <- lr_network_top_df_large
        
      }else{
        # Infer receptors and top-predicted target genes of top-ranked ligands
        ## Active target gene inference
        active_ligand_target_links_df = best_upstream_ligands %>% 
          lapply(nichenetr::get_weighted_ligand_target_links,geneset = geneset_oi, 
                 ligand_target_matrix = ligand_target_matrix, n = n_targets) %>% 
          bind_rows() %>% drop_na()
        colnames(active_ligand_target_links_df) <- c("Ligand", 'Target', 'Score')
        active_ligand_target_links_df$Sender <- sender
        active_ligand_target_links_df$Receiver <- receiver
        tmp.result[[cp]] <- active_ligand_target_links_df
      }
    }
    tmp.result <- do.call(rbind, tmp.result)
    rownames(tmp.result) <- NULL
    tmp.result
  })
  
  result <- do.call(rbind, result)
  
  if(lr){
    message('Return is result of L-R.')
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  }else{
    message('Return is result of L-Targets.')
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Target, sep = '_')
  }
  
  result <- distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run scMLnet method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run scMLnet method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param sender Character, cell types as sender clusters to infer CCC. If NULL, all cell types are used.
#' @param receiver Character, cell types as receiver clusters to infer CCC. If NULL, all cell types are used.
#' @param LigRecLib Dataframe, LR prior database which contains 'source' and 'target' columns. If NULL, the default LR prior database is used.
#' @param RecTFLib Dataframe, R-TF prior database which contains 'source' and 'target' columns. If NULL, the default R-TF prior database is used.
#' @param TFTarLib Dataframe, TF-TG prior database which contains 'source' and 'target' columns. If NULL, the default TF-TG prior database is used.
#'
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender' and 'Receiver' four colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom dplyr distinct
#' @importFrom tidyr separate
#' @importFrom tibble rownames_to_column
#' @importFrom utils read.table packageVersion
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunscMLnet(ser)
#' 
#' # Run with other LR prior database
#' db <- ChangescMLnetDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
#' result <- RunscMLnet(ser, LigRecLib = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangescMLnetDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
#' result <- RunscMLnet(ser, LigRecLib = db)
#' }
RunscMLnet <- function(ser, sender = NULL, receiver = NULL,
                       LigRecLib=NULL, RecTFLib=NULL, TFTarLib=NULL){
  
  set.seed(123)
  
  if(grepl('^4', packageVersion("Seurat"))){
    stop('The version of Seurat should be 3.X, otherwise no results will be inferred.')
  }
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('scMLnet')){
    stop("scMLnet hasn't been installed. \n Please install from https://github.com/SunXQlab/scMLnet2.0!")
  }
  
  if(is.null(LigRecLib)){
    LigRecLib <- read.table(file.path(system.file("extdata", package = "CCCbank"), "scMLnet", 'LigRec.txt'), header = T)
    colnames(LigRecLib)[2:3] <- c("source", "target")
  }
  
  if(is.null(RecTFLib)){
    RecTFLib <- read.table(file.path(system.file("extdata", package = "CCCbank"), "scMLnet", 'RecTF.txt'), header = T)
    colnames(RecTFLib)[1:2] <- c("source", "target")
  }
  
  remove.rec <- setdiff(unique(LigRecLib$target), RecTFLib$source)
  LigRecLib <- LigRecLib[!(LigRecLib$target %in% remove.rec),]
  
  if(is.null(TFTarLib)){
    TFTarLib <- read.table(file.path(system.file("extdata", package = "CCCbank"), "scMLnet", 'TFTargetGene.txt'), header = T)
    colnames(TFTarLib)[1:2] <- c("source", "target")
  }
  
  
  # scMLnet v0.2.0 + database[scMLnet]
  # pacakge: https://github.com/SunXQlab/scMLnet2.0
  GCMat<- GetAssayData(ser, "data", "RNA")
  BarCluTable <- data.frame(Barcode = colnames(ser), Cluster = ser$celltype)
  types <- unique(BarCluTable$Cluster)
  
  if(is.null(sender)){
    message(paste0('Sender: all cell types are used to infer CCC!'))
    sender_ct <- types
  }else{
    message(paste0('Sender: ', paste0(sender, collapse = ', '), ' are used to infer CCC!'))
    sender_ct <- sender
  }
  
  if(is.null(receiver)){
    message(paste0('Receiver: all cell types are used to infer CCC!'))
    receiver_ct <- types
  }else{
    message(paste0('Receiver: ', paste0(receiver, collapse = ', '), ' are used to infer CCC!'))
    receiver_ct <- receiver
  }
  
  result <- list()
  for (LigClu in sender_ct) {
    for (RecClu in receiver_ct[which(receiver_ct != LigClu)]) {
      netList <- tryCatch(scMLnet::RunMLnet(data = GCMat, BarCluTable = BarCluTable,
                                            RecClu = RecClu, LigClu = LigClu,
                                            LigRec.DB = LigRecLib, 
                                            TFTG.DB = TFTarLib,
                                            RecTF.DB = RecTFLib,
                                            Raw.data = FALSE),
                          error=function(e){NA}
      )
      list.names <- paste(LigClu, RecClu, sep = "_")
      result[[list.names]] <- netList
    }
  }
  
  result[which(is.na(result))] <- NULL
  
  result <- lapply(result, function(res){
    res <- res$LigRec
    res
  })
  result <- do.call(rbind, result)
  result <- tibble::rownames_to_column(result, 'sr')
  result$sr <- gsub('\\.[0-9]+', '', result$sr)
  result <- tidyr::separate(result, sr, c('Sender', 'Receiver'), sep = '_')
  result <- result[,-5]
  colnames(result)[3:4] <- c('Ligand', "Receptor")
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run SingleCellSignalR method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run SingleCellSignalR method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species 'human' or 'mouse'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param ... other parameters in \code{\link[SingleCellSignalR]{cell_signaling}}, except for the following parameters: data, genes, cluster, c.names, species, write
#'
#' @seealso \code{\link[SingleCellSignalR]{cell_signaling}}
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom dplyr distinct
#' @importFrom tidyr separate
#' @importFrom utils install.packages
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunSCSR(ser, species = 'human')
#' 
#' # Run with other LR prior database
#' db <- ChangeSCSRDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
#' result <- RunSCSR(ser, species = 'human', priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeSCSRDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
#' result <- RunSCSR(ser, species = 'human', priorDatabase = db)
#' }
RunSCSR <- function(ser, species = NULL, priorDatabase = NULL, ...){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('SingleCellSignalR')){
    message("SingleCellSingalR hasn't been installed. Now installing...")
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    
    BiocManager::install("SingleCellSignalR")
  }
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- "homo sapiens"
  }else if(species == 'human'){
    species <- "homo sapiens"
  }else if(species == 'mouse'){
    species <- "mus musculus"
  }else{
    stop('The species provided is not supported by SingleCellSignalR!\nSupporting species: human, mouse.')
  }
  
  if(!is.null(priorDatabase)){
    assign("LRdb", priorDatabase, envir = .GlobalEnv)
  }else{
    assign("LRdb", SingleCellSignalR::LRdb, envir = .GlobalEnv)
  }
  
  # Get normalized data and cell type metadata
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  meta.sc <- data.frame(celltype = ser$celltype, row.names = colnames(ser))
  
  # Digitize Labels
  i <- 1
  for(ct in unique(meta.sc$celltype)){
    meta.sc[which(meta.sc$celltype == ct), "ct_num"] <- i
    i <- i+1
  }
  c.names <- unique(meta.sc$celltype)
  
  # Run
  signal <- SingleCellSignalR::cell_signaling(data = matrix.sc, 
                                              genes = rownames(matrix.sc),
                                              species = species, 
                                              cluster = meta.sc$ct_num, 
                                              c.names = c.names, 
                                              write = FALSE, ...)
  inter.net <- SingleCellSignalR::inter_network(data = matrix.sc, signal = signal,
                                                genes = rownames(matrix.sc),
                                                cluster = meta.sc$ct_num,
                                                c.names = c.names,
                                                species = species,
                                                write = FALSE)
  # Handle result
  result <- inter.net$`full-network`
  result <- tidyr::separate(result, 'ligand', c('Sender', 'Ligand'), sep = '\\.')
  result <- tidyr::separate(result, 'receptor', c('Receiver','Receptor'), sep = '\\.')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  rownames(result) <- NULL
  
  return(result)
}

#' @title Run CellChat method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run CellChat method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species 'human', 'mouse' or 'zebrafish'
#' @param priorDatabase If NA, run with default LR prior database
#' @param ... other parameters in \code{\link[CellChat]{computeCommunProb}}, except for the following parameter: object
#'
#' @seealso \code{\link[CellChat]{computeCommunProb}}
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import dplyr
#' @importFrom devtools install_github
#' @importFrom tibble rownames_to_column
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunCellChat(ser, species = 'human')
#' 
#' # Run with other LR prior database
#' db <- ChangeCellChatDB(priorDatabase=priorDatabase, extension=FALSE, species = 'human')
#' result <- RunCellChat(ser, species = 'human', priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeCellChatDB(priorDatabase=priorDatabase, extension=TRUE, species = 'human')
#' result <- RunCellChat(ser, species = 'human', priorDatabase = db)
#' }
RunCellChat <- function(ser, species=NULL, priorDatabase = NULL, ...){
  
  #species: support human, mouse and zebrafish
  #...: see details in computeCommunProb function
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
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
  
  if(!is.null(priorDatabase)){
    database <- priorDatabase
  }else{
    database <- origin_db
  }
  
  # Get normalized data and cell type metadata
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  meta.sc <- data.frame(celltype = ser$celltype, row.names = colnames(ser))
  
  # Run
  cellchat <- CellChat::createCellChat(object = matrix.sc, meta = meta.sc, group.by = "celltype")
  cellchat@DB <- database
  cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database
  
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  result <- CellChat::computeCommunProb(cellchat, ...) %>% CellChat::filterCommunication(.) %>% 
    CellChat::subsetCommunication(.)
  
  colnames(result)[1:5] <- c('Sender','Receiver','Ligand', 'Receptor', 'LRscore')
  
  # Get complexes name and its corresponding genes
  complexes <- database$complex
  complexes$all <- paste(complexes$subunit_1, complexes$subunit_2, 
                         complexes$subunit_3, complexes$subunit_4, sep = '&')
  complexes$all <- gsub('&&', '', complexes$all)
  complexes$all <- gsub('&$', '', complexes$all)
  complexes <- tibble::rownames_to_column(complexes, 'complexes')
  complexes <- complexes[, c('complexes', 'all')]
  
  # Replace complexes name with its corresponding genes
  result <- merge(result, complexes, by.x = 'Ligand', by.y = 'complexes', all.x = TRUE)
  result$Ligand[!is.na(result$all)] <- result$all[!is.na(result$all)]
  result$all <- NULL
  result <- merge(result, complexes, by.x = 'Receptor', by.y = 'complexes', all.x = TRUE)
  result$Receptor[!is.na(result$all)] <- result$all[!is.na(result$all)]
  result$all <- NULL
  
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run Connectome method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run Connectome method to infer CCC.
#' 
#' @param ser Seurat object, contains raw and normalized stored in 'counts' and 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species 'human', 'mouse', 'rat' or 'pig'. If provide priorDatabase, species are not limited
#' @param priorDatabase If NULL, run with default LR prior database
#' @param ... other parameters in \code{\link[Connectome]{FilterConnectome}}, except for the following parameter: connectome
#'
#' @seealso \code{\link[Connectome]{FilterConnectome}}
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom devtools install_github
#' @importFrom dplyr distinct
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunConnectome(ser, species = 'human')
#' 
#' # Run with other LR prior database
#' db <- ChangeConnectomeDB(priorDatabase=priorDatabase, extension=FALSE, 
#'                          keep_complexes = TRUE, species = 'human')
#' result <- RunConnectome(ser, species = 'human', priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeConnectomeDB(priorDatabase=priorDatabase, extension=TRUE, 
#'                          keep_complexes = TRUE, species = 'human')
#' result <- RunConnectome(ser, species = 'human', priorDatabase = db)
#' }
RunConnectome <- function(ser, species = NULL, priorDatabase = NULL, ...){
  
  set.seed(123)
  
  if(!requireNamespace('Connectome')){
    message("Connectome hasn't been installed. Now installing...")
    devtools::install_github('msraredon/Connectome', ref = 'master')
  }
  
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
    if(is.null(priorDatabase)){
      stop('The species provided is not supported by CellChat!\nSupporting species: human, mouse, rat, pig.\nIf provide LR prior database, species are not limited')
    }
  }
  
  if(!is.null(priorDatabase)){
    database <- priorDatabase
  }else{
    database <- origin_db
  }
  database <- database[, c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'mode')]
  
  connectome.genes <- union(database$Ligand.ApprovedSymbol, database$Receptor.ApprovedSymbol)
  genes <- connectome.genes[connectome.genes %in% rownames(ser)]
  ser <- ScaleData(ser,features = genes)
  
  # Run
  sc.con <- Connectome::CreateConnectome(ser, LR.database = 'custom', custom.list = database, calculate.DOR = TRUE)
  result <- Connectome::FilterConnectome(sc.con,...)
  
  #result <- result[,c(1:4, 17)]
  colnames(result)[c(1:4, 17)] <- c('Sender', 'Receiver', 'Ligand', 'Receptor', 'LRscore')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run CytoTalk method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run CytoTalk method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species 'human' or 'mouse'
#' @param priorDatabase If NULL, run with default LR prior database
#' @param celltypes Character, cell types used to infer CCC. If NULL, all cell types are used
#' @param ... other parameters in \code{\link[CytoTalk]{run_cytotalk}}, except for the following parameters: lst_scrna, cell_type_a, cell_type_b, pcg, lrp

#'
#' @seealso \code{\link[CytoTalk]{run_cytotalk}}
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import reticulate
#' @importFrom devtools install_github
#' @importFrom dplyr distinct
#' @importFrom utils write.csv combn
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunCytoTalk(ser, species = 'human')
#' 
#' # Run with other LR prior database
#' db <- ChangeCytoTalkDB(priorDatabase=priorDatabase, extension=FALSE, 
#'                        keep_complexes = TRUE, species = 'human')
#' result <- RunCytoTalk(ser, species = 'human', priorDatabase = db)
#' 
#' # Run with default and other LR prior databases
#' db <- ChangeCytoTalkDB(priorDatabase=priorDatabase, extension=TRUE, 
#'                        keep_complexes = TRUE, species = 'human')
#' result <- RunCytoTalk(ser, species = 'human', priorDatabase = db)
#' }
RunCytoTalk <- function(ser, species = NULL, priorDatabase = NULL, celltypes = NULL, ...){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('CytoTalk')){
    message("CytoTalk hasn't been installed. Now installing...")
    devtools::install_github("tanlabcode/CytoTalk")
  }
  
  # create conda environment
  # To install and call Python modules from R.
  if(!('r_reticulate_CytoTalk' %in% reticulate::conda_list()$name)){
    reticulate::conda_create(envname = "r_reticulate_CytoTalk", python_version = "3.7.3")  # Create a new Conda environment to facilitate the Python module installation.
    reticulate::conda_install(envname = "r_reticulate_CytoTalk", "pybind11")  # Install two necessary Python modules for correctly compiling and using the "pcst_fast" Python module.
    reticulate::conda_install(envname = "r_reticulate_CytoTalk", "numpy")
    #reticulate::conda_install(envname = "r_reticulate_CytoTalk", "git+https://github.com/fraenkel-lab/pcst_fast.git", pip = TRUE) # To install the "pcst_fast" module.
    reticulate::conda_install(envname = "r_reticulate_CytoTalk", "pcst_fast", pip = TRUE)
  }
  
  # input
  if(T){
    fpath <-  paste0(getwd(), '/CytoTalkTemp')
    if(!dir.exists(fpath)){
      dir.create(fpath)
    }
    message(paste0('The files generated during the process are stored in ', fpath))
    
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, "counts.csv")
    if(!file.exists(fpath.mat)){
      norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
      write.csv(norm.matrix, fpath.mat, quote=F)
    }
    
    fpath.meta <- paste0(input_fpath,"metadata.csv")
    if(!file.exists(fpath.meta)){
      meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)  
      meta.data <- as.matrix(meta.data)
      write.csv(meta.data, fpath.meta, quote=F, row.names = FALSE)
    }
    
    #rm(norm.matrix, meta.data, ser);
  }
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- 'human'
  }else if(!(species %in% c('human', 'mouse'))){
    stop('The species provided is not supported by CytoTalk!\nSupporting species: human, mouse.')
  }
  
  if(!is.null(priorDatabase)){
    lrp <- priorDatabase$lrp
    pcg <- priorDatabase$pcg
  }else{
    if(species == 'human'){
      pcg = CytoTalk::pcg_human
      lrp = CytoTalk::lrp_human
    }else if(species == 'mouse'){
      pcg = CytoTalk::pcg_mouse
      lrp = CytoTalk::lrp_mouse
    }
  }
  
  lst.sc <- CytoTalk::read_matrix_with_meta(fpath.mat, fpath.meta)
  celltype <- unique(lst.sc$cell_types)
  
  if(is.null(celltypes)){
    message('All cell types are used to infer CCC!')
    celltypes <- celltype
  }else{
    message(paste0(paste0(celltypes, collapse = ','), ' are used to infer CCC!'))
  }
  celltypes = as.character(celltypes)
  comb <- combn(celltypes, 2)
  comb <- t(comb)
  comb <- as.data.frame(comb)
  
  result <- list()
  for (i in 1:dim(comb)[[1]]) {
    tmp <- tryCatch(CytoTalk::run_cytotalk(lst.sc, comb[i, 1], comb[i, 2], 
                                           pcg = pcg, lrp = lrp, ...),
                    error=function(e){NA}
    )
    cp <- paste(comb[i, 1], comb[i, 2], sep = "_")
    result[[cp]] <- tmp
  }
  
  result[which(is.na(result))] <- NULL
  
  result <- lapply(result, function(res){res$pcst$final_network})
  result <- do.call(rbind, result)
  
  if(!is.null(result){
    if(species=='human'){
    result$node1 <- toupper(result$node1)
    result$node2 <- toupper(result$node2)
    result$node1 <- gsub("ORF", "orf", result$node1)
    result$node2 <- gsub("ORF", "orf", result$node2)
  }
  colnames(result)[c(1:4,10)] <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  rownames(result) <- NULL
  }else{
    result <- NA
    message('No ligand-receptor interactions were inferred!')
  }
  
  return(result)
}

#' @title Run scSeqComm method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run scSeqComm method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species 'human' or 'mouse'
#' @param priorDatabase If NULL, run with default LR prior database (origin_lr)
#' @param origin_lr see details in data('scSeqComm'); If NULL, scSeqComm::LR_pairs_Jin_2020 or scSeqComm::LR_pairs_Jin_2020_mouse is used.
#' @param origin_rectf see details in data('scSeqComm'); If NULL, scSeqComm::TF_PPR_KEGG_human or scSeqComm::TF_PPR_KEGG_mouse is used.
#' @param origin_tftg see details in data('scSeqComm'); If NULL, scSeqComm::TF_TG_RegNetwork_High or scSeqComm::TF_TG_RegNetwork_High_mouse is used.
#' @param S_inter_cutoff Numeric between 0 to 1. The cutoff of S_inter is used to filter the L-R pairs.
#' @param S_intra_cutoff Numeric between 0 to 1. The cutoff of S_intra is used to filter the L-R-Pathway links.
#' @param ... other parameters in \code{\link[scSeqComm]{scSeqComm_analyze}}, except for the following parameters: gene_expr, cell_metadata, LR_pairs_DB, TF_reg_DB, R_TF_association
#'
#' @seealso \code{\link[scSeqComm]{scSeqComm_analyze}}
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom devtools install_gitlab
#' @importFrom dplyr distinct filter
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunscSeqComm(ser, species = 'human')
#' 
#' # Run with other LR prior database
#' db <- ChangescSeqCommDB(priorDatabase=priorDatabase)
#' result <- RunscSeqComm(ser, species = 'human', priorDatabase = db)
#' }
#' 
RunscSeqComm <- function(ser, species = NULL, priorDatabase = NULL,
                         origin_lr=NULL, origin_rectf=NULL, origin_tftg=NULL,
                         S_inter_cutoff = 0.8, S_intra_cutoff = 0.8, ...){
  
  #origin_lr, origin_rectf and origin_tftg: see details in scSeqComm
 
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('scSeqComm')){
    message("scSeqComm hasn't been installed. Now installing...")
    devtools::install_gitlab("sysbiobig/scseqcomm")
  }
  
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- 'human'
  }else if(!(species %in% c('human', 'mouse'))){
    stop('The species provided is not supported by scSeqComm!\nSupporting species: human and mouse')
  }
  
  ## Ligand - receptor pairs
  if(is.null(origin_lr)){
    
    if(species=='human'){
      if(is.null(priorDatabase)){
        message('No LR_pairs_DB provided, LR_pairs_DB is set to scSeqComm::LR_pairs_Jin_2020')
        origin_lr <- scSeqComm::LR_pairs_Jin_2020
      }
    }else{
      if(is.null(priorDatabase)){
        message('No LR_pairs_DB provided, LR_pairs_DB is set to scSeqComm::LR_pairs_Jin_2020_mouse')
        origin_lr <- scSeqComm::LR_pairs_Jin_2020_mouse
      }
    }
    
  }
  
  ## Receptor-Transcription factor a-priori association from gene signaling networks
  if(is.null(origin_rectf)){
    
    if(species=='human'){
      message('No R_TF_association provided, R_TF_association is set to scSeqComm::TF_PPR_KEGG_human')
      origin_rectf <- scSeqComm::TF_PPR_KEGG_human
    }else{
      message('No R_TF_association provided, R_TF_association is set to scSeqComm::TF_PPR_KEGG_mouse')
      origin_rectf <- scSeqComm::TF_PPR_KEGG_mouse
    }
    
  }
  
  ## Transcriptional regulatory network
  if(is.null(origin_tftg)){
    if(species=='human'){
      message('No TF_reg_DB provided, TF_reg_DB is set to scSeqComm::TF_TG_RegNetwork_High')
      origin_tftg <- scSeqComm::TF_TG_RegNetwork_High
    }else{
      message('No TF_reg_DB provided, TF_reg_DB is set to scSeqComm::TF_TG_RegNetwork_High_mouse')
      origin_tftg <- scSeqComm::TF_TG_RegNetwork_High_mouse
    }
  }
  
  if(!is.null(priorDatabase)){
    db.lr <- priorDatabase
  }else{
    db.lr <- origin_lr
  }
  
  ## scRNA-seq dataset
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  cell_metadata <- data.frame(Cell_ID = colnames(ser),
                              Cluster_ID = ser$celltype)
  #cell_cluster <- lapply(unique(ser$celltype), function(ct){
  #  meta <- ser@meta.data
  #  cells <- rownames(meta)[which(meta$celltype == ct)]
  #  cells
  #})
  #names(cell_cluster) <- unique(ser$celltype)
  
  ## Intercellular and intracellular signaling analysis
  scSeqComm_res <- scSeqComm::scSeqComm_analyze(gene_expr = matrix.sc,
                                                cell_metadata = cell_metadata,
                                                LR_pairs_DB = db.lr,
                                                TF_reg_DB = origin_tftg,
                                                R_TF_association = origin_rectf, ...)
  message(paste0('The cutoff of S_inter is ', S_inter_cutoff))
  message(paste0('The cutoff of S_intra is ', S_intra_cutoff))
  
  result <- scSeqComm_res$comm_results
  result <- dplyr::filter(result, S_inter>S_inter_cutoff, S_intra>S_intra_cutoff)
  result <- result[, c(1:2,4:5,9)]
  colnames(result) <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result$Ligand <- gsub(',', '&', result$Ligand)
  result$Receptor <- gsub(',', '&', result$Receptor)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}

#' @title Run CellPhoneDB2 method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run CellPhoneDB2 method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param db_path Folder for storing database. If using default database of CellPhoneDB, the db_path should be set to NULL
#' @param cores Number of cores used to run in parallel. Defaults to 1
#' 
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import reticulate
#' @import dplyr
#' @importFrom tidyr separate pivot_longer
#' @importFrom utils write.table read.csv
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunCellPhoneDB2(ser)
#' 
#' # Run with other LR prior database
#' db_path <- ChangeCellPhoneDB(priorDatabase=priorDatabase, extension=FALSE)
#' result <- RunCellPhoneDB2(ser, db_path = db_path)
#' 
#' # Run with default and other LR prior databases
#' db_path <- ChangeCellPhoneDB(priorDatabase=priorDatabase, extension=TRUE)
#' result <- RunCellPhoneDB2(ser, db_path = db_path)
#' }
#' 
RunCellPhoneDB2 <- function(ser, db_path=NULL, cores =1){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  # input
  if(T){
    fpath <-  paste0(getwd(), '/CellPhoneDB2Temp')
    if(!dir.exists(fpath)){
      dir.create(fpath)
    }
    message(paste0('The files generated during the process are stored in ', fpath))
    
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, "counts.txt")
    if(!file.exists(fpath.mat)){
      norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
      write.table(norm.matrix, fpath.mat, sep='\t', quote=F)
    }
    
    fpath.meta <- paste0(input_fpath,"metadata.txt")
    if(!file.exists(fpath.meta)){
      meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)
      colnames(meta.data) <- NULL
      meta.data <- as.matrix(meta.data)
      write.table(meta.data, fpath.meta, sep='\t', quote=F, row.names=F)
    }
    
    #rm(norm.matrix, meta.data, ser);gc()
    
    output_fpath <- paste0(fpath, '/output/')
    if(!dir.exists(output_fpath)){
      dir.create(output_fpath)
    }
  }
  
  if(is.null(db_path)){
    db_path <- file.path(system.file("extdata", package = "CCCbank"), 'CellPhoneDB')
    complex.path <- file.path(db_path,'complexes.csv')
    gene.path <- file.path(db_path, 'genes.csv')
    protein.path <- file.path(db_path, 'proteins.csv')
  }else if(db_path == file.path(system.file("extdata", package = "CCCbank"), 'CellPhoneDB')){
    db_path <- file.path(system.file("extdata", package = "CCCbank"), 'CellPhoneDB')
    complex.path <- file.path(db_path,  'complexes.csv')
    gene.path <- file.path(db_path,  'genes.csv')
    protein.path <- file.path(db_path,  'proteins.csv')
  }else{
    complex.path <- file.path(db_path, 'complex_input.csv')
    gene.path <- file.path(db_path, 'gene_input.csv')
    protein.path <- file.path(db_path, 'protein_input.csv')
  }
  
  if(!('cpdb' %in% reticulate::conda_list()$name)){
    PyHome <- CreatCondaEnv(condaName='cpdb', python_version = '3.7', packages = 'cellphonedb')
  }
  condaName <- 'cpdb'
  
  condaHome <- reticulate::conda_binary()
  condaHome <- sub('\\/conda$', '', condaHome)
  
  platform <- GetSysInfo()
  if(platform=='linux'){
    activateCommand <- paste0('source ', condaHome, '/activate ', condaName)
  }else if(platform=='windows'){
    activateCommand <- paste0('conda.bat activate ', condaName)
  }
  
  database_fpath <- list.files(db_path, pattern = '\\.db', full.names = TRUE)
  if (length(database_fpath)>0) {
    run_command <- paste('cellphonedb method statistical_analysis', fpath.meta, fpath.mat, 
                         '--threshold 0.05', 
                         '--output-path', output_fpath, 
                         '--database', database_fpath, 
                         '--threads', cores,
                         '--counts-data gene_name', 
                         '--debug-seed 123', sep = ' ')
  }else{
    run_command <- paste('cellphonedb method statistical_analysis', fpath.meta, fpath.mat, 
                         '--threshold 0.05', 
                         '--output-path', output_fpath, 
                         '--threads', cores,
                         '--counts-data gene_name', 
                         '--debug-seed 123', sep = ' ')
  }
  
  command <- paste(activateCommand, run_command,sep = ifelse(platform=='windows', ' & ', ' ; '))
  system(command)
  
  # Get Reuslt
  cpdb.complex <- RecoverCpdbComplex(gene.path, protein.path, complex.path)
  
  file.path <- paste(output_fpath, 'significant_means.txt', sep = '/')
  
  result <- ProcssCPDBResult(cpdb.complex, file.path)
  
  return(result)
}

#' @title Run CellPhoneDB3 method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run CellPhoneDB3 method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param db_path Folder for storing database. If using default database of CellPhoneDB, the db_path should be set to NULL
#' @param cores Number of cores used to run in parallel. Defaults to 1 
#' 
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import reticulate
#' @import dplyr
#' @importFrom tidyr separate pivot_longer
#' @importFrom utils write.table read.csv packageVersion
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunCellPhoneDB3(ser)
#' 
#' # Run with other LR prior database
#' db_path <- ChangeCellPhoneDB(priorDatabase=priorDatabase, extension=FALSE)
#' result <- RunCellPhoneDB3(ser, db_path = db_path, condaName = 'cpdb')
#' 
#' # Run with default and other LR prior databases
#' db_path <- ChangeCellPhoneDB(priorDatabase=priorDatabase, extension=TRUE)
#' result <- RunCellPhoneDB3(ser, db_path = db_path, condaName = 'cpdb')
#' }
#' 
RunCellPhoneDB3 <- function(ser, db_path=NULL, cores =1){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  # input
  if(T){
    fpath <-  paste0(getwd(), '/CellPhoneDB3Temp')
    if(!dir.exists(fpath)){
      dir.create(fpath)
    }
    message(paste0('The files generated during the process are stored in ', fpath))
    
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, "counts.txt")
    fpath.meta <- paste0(input_fpath,"metadata.txt")
    fpath.deg <- paste0(input_fpath, 'degs.txt')
    
    if(!file.exists(fpath.mat)){
      norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
      write.table(norm.matrix, fpath.mat, sep='\t', quote=F)
    }
   
    if(!file.exists(fpath.meta)){
      meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)
      colnames(meta.data) <- NULL
      meta.data <- as.matrix(meta.data)
      write.table(meta.data, fpath.meta, sep='\t', quote=F, row.names=F)
    }   
    
    if(!file.exists(fpath.deg)){
      Idents(ser) <- ser$celltype
      DEGs <- FindAllMarkers(ser, test.use = 't',
                             verbose = F, only.pos = T, 
                             random.seed = 123, logfc.threshold = 0.15, 
                             min.pct = 0.05, return.thresh = 0.05)
      if(grepl('^4', packageVersion('Seurat'))){
        fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.15)
        fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')]
      }else if(grepl('^3', packageVersion('Seurat'))){
        fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_logFC > 0.15)
        fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_logFC', 'pct.1', 'pct.2')]
        colnames(fDEGs)[5] <- 'avg_log2FC'
      }
      
      write.table(fDEGs, file = fpath.deg, sep = '\t', quote = F, row.names = F)
    }
    #rm(norm.matrix, meta.data, DEGs, fDEGs, ser);gc()
    
    output_fpath <- paste0(fpath, '/output/')
    if(!dir.exists(output_fpath)){
      dir.create(output_fpath)
    }
  }
  
  if(is.null(db_path)){
    db_path <- file.path(system.file("extdata", package = "CCCbank"), 'CellPhoneDB')
    complex.path <- file.path(db_path,'complexes.csv')
    gene.path <- file.path(db_path, 'genes.csv')
    protein.path <- file.path(db_path, 'proteins.csv')
  }else if(db_path == file.path(system.file("extdata", package = "CCCbank"), 'CellPhoneDB')){
    db_path <- file.path(system.file("extdata", package = "CCCbank"), 'CellPhoneDB')
    complex.path <- file.path(db_path,  'complexes.csv')
    gene.path <- file.path(db_path,  'genes.csv')
    protein.path <- file.path(db_path,  'proteins.csv')
  }else{
    complex.path <- file.path(db_path, 'complex_input.csv')
    gene.path <- file.path(db_path, 'gene_input.csv')
    protein.path <- file.path(db_path, 'protein_input.csv')
  }
  
  if(!('cpdb' %in% reticulate::conda_list()$name)){
    PyHome <- CreatCondaEnv(condaName='cpdb', python_version = '3.7', packages = 'cellphonedb')
  }
  condaName <- 'cpdb'
  
  platform <- GetSysInfo()
  if(platform=='linux'){
    condaHome <- reticulate::conda_binary()
    condaHome <- sub('\\/conda$', '', condaHome)
    activateCommand <- paste0('source ', condaHome, '/activate ', condaName)
  }else if(platform=='windows'){
    activateCommand <- paste0('conda.bat activate ', condaName)
  }
  
  database_fpath <- list.files(db_path, pattern = '\\.db', full.names = TRUE)
  if (length(database_fpath)>0) {
    run_command <- paste('cellphonedb method degs_analysis', fpath.meta, fpath.mat, fpath.deg, 
                         '--threshold 0.05', 
                         '--output-path', output_fpath, 
                         '--database', database_fpath, 
                         '--threads', cores,
                         '--counts-data gene_name', 
                         '--debug-seed 123', sep = ' ')
  }else{
    run_command <- paste('cellphonedb method degs_analysis', fpath.meta, fpath.mat, fpath.deg, 
                         '--threshold 0.05', 
                         '--output-path', output_fpath, 
                         '--threads', cores,
                         '--counts-data gene_name', 
                         '--debug-seed 123', sep = ' ')
  }
  
  command <- paste(activateCommand, run_command,sep = ifelse(platform=='windows', ' & ', ' ; '))
  system(command)
  
  # Get Reuslt
  cpdb.complex <- RecoverCpdbComplex(gene.path, protein.path, complex.path)
  
  file.path <- paste(output_fpath, 'significant_means.txt', sep = '/')
  
  result <- ProcssCPDBResult(cpdb.complex, file.path)
  
  return(result)
}

#' @title Run scConnect method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run scConnect method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species 'human', 'mouse' and others
#' @param PyHome The execution file path of python, which has installed scConnect. If NULL, python installed in systems will be found and install scConnect
#' @param p_cutoff The cutoff of p.value of ligands and receptors. Defaults to 0.05
#' 
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import reticulate
#' @importFrom dplyr distinct
#' @importFrom utils write.csv read.csv
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunscConnect(ser, species ='human')
#' 
#' # Run with other LR prior database
#' PyHome <- ChangescConnectDB(priorDatabase=priorDatabase, extension=FALSE, 
#'                             keep_complexes = TRUE, species = 'human')
#' result <- RunscConnect(ser, species ='human', PyHome = PyHome)
#' 
#' # Run with default and other LR prior databases
#' PyHome <- ChangescConnectDB(priorDatabase=priorDatabase, extension=TRUE, 
#'                             keep_complexes = TRUE, species = 'human')
#' result <- RunscConnect(ser, species ='human', PyHome = PyHome)
#' }
#' 
RunscConnect <- function(ser, species = NULL, PyHome =NULL, p_cutoff = 0.05){
  
  #species: support human, mouse and more
  #PyHome: install scConnect in the PyHome
  #p_cutoff: filter ligands and receptors with p.value less than 0.05
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
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
  
  fpath <-  paste0(getwd(), '/scConnectTemp')
  if(!dir.exists(fpath)){
    dir.create(fpath)
  }
  message(paste0('The files generated during the process are stored in ', fpath))
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  fpath.mat <- paste0(input_fpath, "counts.csv")
  if(!file.exists(fpath.mat)){
    norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
    write.csv(norm.matrix, fpath.mat, quote=F)
  }

  fpath.meta <- paste0(input_fpath, "metadata.csv")
  if(!file.exists(fpath.meta)){
    cell.meta <- data.frame(Cell = rownames(ser@meta.data), Annotation = ser$celltype)
    write.csv(cell.meta, fpath.meta, quote=F, row.names = FALSE)
  }
 
  #rm(norm.matrix, cell.meta, ser);gc()
  
  output_fpath <- paste0(fpath, '/output/')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  output_fpath <- paste0(output_fpath, 'result.csv')
  
  # Check PyHome
  PyHome <- CheckscConnect(PyHome)
  
  # get python script of scConnect
  inst_path <- system.file("python_scripts", package = "CCCbank")
  python_script_path <- file.path(inst_path, "RunscConnect.py")
  
  # Run command
  command <- paste(PyHome, python_script_path, 
                   fpath.mat, fpath.meta, 
                   species, output_fpath, sep = ' ')
  system(command)
  
  # handle result
  if(T){
    result <- read.csv(output_fpath)
    result <- result[,-1]
    
    db_path <- GetDbPathOfscConnect(PyHome)
    ligands_path <- file.path(db_path, species, 'ligands.csv')
    ligands <-  read.csv(ligands_path, sep = ',')
    ligands <- ligands[,-c(1, 3, 9)]
    ligands$ligand_symbol <- gsub("\\['", '', ligands$preprogene)
    ligands$ligand_symbol <- gsub("\\']", '', ligands$ligand_symbol)
    ligands$ligand_symbol <- gsub("', '", '|', ligands$ligand_symbol)
    #ligands$ligand_symbol[which(ligands$ligand == 'CCL19')] <- 'CCL19'
    #ligands$ligand_symbol[which(ligands$ligand == 'Ccl21a')] <- 'CCL21'
    #ligands$ligand_symbol[which(ligands$ligand == 'Ccl21b')] <- 'CCL21'
    #ligands$ligand_symbol[which(ligands$ligand == 'agouti')] <- 'ASIP'
    #ligands$ligand_symbol[which(ligands$ligand == 'relaxin')] <- 'RLN2|RLN1'
    
    result <- merge(ligands, result,by = "ligand")
    colnames(result)[7] <- "ligand_gene"
    result <- result[,-c(2:6)]
    
    receptors_path <- file.path(db_path, species, 'receptors.csv')
    receptors <-  read.csv(receptors_path, sep = ',')
    receptors <- dplyr::distinct(receptors, receptor, gene, .keep_all = FALSE)
    receptors$receptor_symbol <- gsub("\\['", '', receptors$gene)
    receptors$receptor_symbol <- gsub("\\']", '', receptors$receptor_symbol)
    receptors$gene <- NULL
    
    result <- merge(receptors, result, by = 'receptor')
    colnames(result)[2] <- 'receptor_gene'
    message(paste0("Filtering ligands and receptors with p.value less than ", p_cutoff))
    result <- result[which(result$ligand_pval<p_cutoff), ]
    result <- result[which(result$receptor_pval<p_cutoff),]
    result <- result[,c("sender", "reciever", "ligand_gene", "receptor_gene", "score")]
    colnames(result) <- c("Sender", "Receiver", "Ligand", "Receptor", "LRscore")
    result <- tidyr::separate_rows(result, Ligand, sep = '\\|')
    result <- tidyr::separate_rows(result, Receptor, sep = '\\|')
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  RecoverDefaultDatabase(species, PyHome)
  
  return(result)
}

#' @title Run Domino method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run Domino method to infer CCC.
#' 
#' @param ser Seurat object, contains raw and normalized data stored in 'counts' and 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param db_path Folder for storing database. If NULL, using the default database of CellPhoneDB
#' @param fpath.auc File path of the output when running pyscenic
#' @param fpath.reg File path of the output when running pyscenic
#' @param ... other parameters in \code{\link[domino]{build_domino}}, except for the following parameters: dom
#' 
#' @seealso \code{\link[domino]{build_domino}}
#' 
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @import Matrix
#' @import dplyr
#' @importFrom devtools install_github
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_replace_all
#' @importFrom utils read.table
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunDomino(ser, fpath.auc=fpath.auc, fpath.reg=fpath.auc)
#' 
#' # Run with other LR prior database
#' db_path <- ChangeDominoDB(priorDatabase=priorDatabase)
#' result <- RunDomino(ser, db_path = db_path, fpath.auc=fpath.auc, fpath.reg=fpath.auc)
#' }
#' 
RunDomino <- function(ser, db_path=NULL, fpath.auc=NULL, fpath.reg=NULL, ...){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  if(!requireNamespace('domino')){
    message("domino hasn't been installed. Now installing...")
    devtools::install_github('Chris-Cherry/domino')
  }
  
  if(is.null(db_path)){
    db_path <- file.path(system.file("extdata", package = "CCCbank"), "CellPhoneDB")
  }
  
  if(is.null(fpath.auc) & is.null(fpath.reg)){
    stop("No fpath.auc and fpath.reg provided! Please run scenic first!")
  }
  
  # run domino
  if(T){
    auc = t(read.table(fpath.auc, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = ','))
    
    ser <- ser[, colnames(auc)]
    ser <- ScaleData(ser, features = rownames(ser))
    
    # line 200 of rscript of creat_domino function: change tf_genes = df[row, 10] to tf_genes = df[row, 9]
    source(system.file("R_scripts", "create_domino_code.R", package = "CCCbank")) 
    
    dom = create_domino(signaling_db = db_path, 
                        features = auc, counts = ser@assays$RNA@counts, z_scores = ser@assays$RNA@scale.data, 
                        clusters = ser@active.ident, df = fpath.reg)
    dom = domino::build_domino(dom, ...)
  }
  
  # Bulid ligand-receptor dataframe
  if(T){
    lig_to_rec <- data.frame()
    for (rec in names(dom@linkages$rec_lig)) {
      lig <- dom@linkages$rec_lig[[rec]]
      
      if(length(lig)>0){
        lig.rec <- data.frame(ligand = lig, receptor = rep(rec, length(lig)))
      }else{
        lig.rec <- data.frame(ligand = lig, receptor = rec)
      }
      
      lig_to_rec <- rbind(lig_to_rec, lig.rec)
    }
    lig_to_rec <- lig_to_rec[which(lig_to_rec$ligand!=""),]
  }
  
  # Handle result
  if(T){
    result <- lapply(levels(dom@clusters), function(reciever){
      print(reciever)
      # get ligands of sender cells which are communicated with reciever cells
      sender.ligands <- dom@cl_signaling_matrices[[reciever]] %>%
        as.data.frame(.) %>% tibble::rownames_to_column(var = "ligand") %>%
        tidyr::pivot_longer(cols = -ligand, names_to = "sender", values_to = "expression") %>%
        .[which(.$expression>0), ]
      sender.ligands$sender <- stringr::str_replace_all(sender.ligands$sender, "L_", "")
      sender.ligands <- sender.ligands[, -3]
      
      if(dim(sender.ligands)[1] != 0){
        # get tfs of reciever cells
        reciever.tfs <- dom@linkages$clust_tf[[reciever]]
        rec_tf <- data.frame()
        for(tf in reciever.tfs){
          rec <- dom@linkages$tf_rec[[tf]]
          
          if(length(rec)>0){
            rec.tf <- data.frame(receptor = rec, tf = rep(tf, length(rec)))
          }else{
            rec.tf <- data.frame()
          }
          rec_tf <- rbind(rec_tf, rec.tf)
        }
        
        lig_rec_tf <- merge(lig_to_rec, rec_tf, by = "receptor")
        lig_rec_tf <- merge(sender.ligands, lig_rec_tf, by = "ligand")
        lig_rec_tf$reciever <- reciever
        lig_rec_tf$tf <- stringr::str_replace_all(lig_rec_tf$tf, "\\...", "")
      }else{
        lig_rec_tf <- NA
      }
      lig_rec_tf
    })
    
    result[which(is.na(result))] <- NULL
    
    result <- do.call(rbind, result)
    result <- result[,-4]
    colnames(result) <- c('Ligand', 'Sender', 'Receptor', 'Receiver')
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  
  return(result)
}

#' @title Run NATMI method to infer cell-cell communication (CCC).
#' @author Jiaxin Luo
#' @description This function is to run NATMI method to infer CCC.
#' 
#' @param ser Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'
#' @param species "human", "mouse" and more(21 species)
#' @param cores Number of cores used to run in parallel. Defaults to NULL
#' @param filter.perc Numeric, filtering ligands and receptors expressed in less than filter.perc*100 percentage of cells. Defaults to 0.05
#' @param priorDatabase Logic. If FALSE, using default LR prior database to infer CCC
#' @param NATMI_dir Folder path of source code of NATMI. If NA, it will download the source code of NATMI from github in the current directory.
#' @param PyHome The execution file path of python
#' 
#' @seealso \code{\link[domino]{build_domino}}
#' 
#' @return Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom dplyr distinct
#' @importFrom utils write.csv read.csv
#' 
#' @examples
#' \dontrun{
#' # Run with default LR prior database
#' result <- RunNATMI(ser, species ='human', priorDatabase = FALSE)
#' 
#' # Run with other LR prior database
#' NATMI_dir <- ChangeNATMIDB(priorDatabase=NULL, extension=FALSE, keep_complexes = TRUE)
#' result <- RunNATMI(ser, species ='human', priorDatabase = TRUE, NATMI_dir = NATMI_dir)
#' 
#' # Run with default and other LR prior databases
#' NATMI_dir <- ChangeNATMIDB(priorDatabase=NULL, extension=FALSE, keep_complexes = TRUE)
#' result <- RunNATMI(ser, species ='human', priorDatabase = TRUE, NATMI_dir = NATMI_dir)
#' }
#' 
RunNATMI <- function(ser, species = NULL, cores = NULL, filter.perc = 0.05,
                     priorDatabase=FALSE, NATMI_dir = NA, PyHome = NULL){
  
  set.seed(123)
  
  CheckSeuratObject(ser)
  
  species_support <- c("human", "mouse", "rat", "zebrafish", "fruitfly", "chimpanzee", 
                       "dog", "monkey", "cattle", "chicken", "frog", "mosquito", "nematode", 
                       "thalecress", "rice", "riceblastfungus", "bakeryeast", "neurosporacrassa", 
                       "fissionyeast", "eremotheciumgossypii", "kluyveromyceslactis")
  if(is.null(species)){
    message("No 'species' provied, 'species' is set to 'human'")
    species <- 'human'
  }else if(!(species %in% species_support)){
    stop(paste0('The species provided is not supported by NATMI!\nSupporting species: ', paste0(species_support, collapse = ', ')))
  }
  
  fpath <-  paste0(getwd(), '/NATMITemp')
  if(!dir.exists(fpath)){
    dir.create(fpath)
  }
  message(paste0('The files generated during the process are stored in ', fpath))
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  output_fpath <- paste0(fpath, '/output/')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  
  if(is.na(NATMI_dir)){
    NATMI_dir <- Install_NATMI()
  }
  
  fpath.mat <- paste0(input_fpath, "counts.csv")
  if(!file.exists(fpath.mat)){
    write.csv(100 * (exp(as.matrix(GetAssayData(object = ser, assay = "RNA", slot = "data"))) - 1), 
              fpath.mat, row.names = T)
  }
  
  fpath.meta <- paste0(input_fpath, "metadata.csv")
  if(!file.exists(fpath.meta)){
    meta <- data.frame(Cell = colnames(ser), Annotation = ser$celltype)
    write.csv(meta,fpath.meta, row.names = FALSE)
  }
  
  
  platform <- GetSysInfo()
  
  if(!is.null(cores)){
    message(paste0('Running in parallel with ', cores, ' cores...'))
  }else{
    cores <- 1
  }
  
  if(is.null(PyHome)){
    PyHome <- Sys.which("python")
    PyHome <- ifelse(PyHome=='', Sys.which("python3"), PyHome)
    if(PyHome==''){
      stop("Can't find python2/3 in the system!")
    }
  }
  message(paste0('Current PyHome is ', PyHome))
  
  command1 <- paste0('cd ', NATMI_dir)
  command2 <- paste0(PyHome, ' ExtractEdges.py --species ', species,
                     ' --emFile ', fpath.mat,
                     ' --annFile ', fpath.meta,
                     ' --interDB ', ifelse(priorDatabase, 'lrdb', 'lrc2p'),
                     ' --coreNum ', cores,
                     ' --out ', output_fpath)
  if(platform=='windows'){
    command1 <- paste0('cmd.exe /c ',command1)
  }
  command <- paste(command1, command2, sep = ifelse(platform == 'windows', '&', ';'))
  system(command)
  
  result_path <- paste0(output_fpath, ifelse(priorDatabase, 'Edges_lrdb.csv', 'Edges_lrc2p.csv'))
  result <- read.csv(result_path)
  colnames(result)[c(1:4, 17)] <- c("Sender", "Ligand", "Receptor", "Receiver", "LRscore")
  message(paste0('Filtering ligands and receptors expressed in less than ', filter.perc*100, '% of cells...'))
  result <- filter(result, Ligand.detection.rate > filter.perc, 
                   Receptor.detection.rate > filter.perc)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  return(result)
}
