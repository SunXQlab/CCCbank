#' @title Infer cell-cell communication (CCC) with any combination of methods and database
#' @author Jiaxin Luo
#' 
#' @param ser Seurat object, contains raw and normalized data stored in 'counts' and 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'.
#' @param method Character, select one method to infer CCC, and 16 methods are supported.
#' @param species human, mouse or other species.
#' @param databases Character, select one or more LR prior databases of 13 databases. Also, it can be dataframe with two colunms: ligand, receptor. If NULL, the default LR prior database of method will be used.
#' @param extension Logic. When the parameter 'databases' is not NULL, whether use the combination of default database of method and other databases?
#' @param keep_complexes Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit, when the parameter 'databases' is not NULL and the method doesn't consider multi-subunit of ligand and receptor. Otherwise, the interactions containing multi-subunit will be filtered out.
#' @param ... Other parameters in the corresponding method function.
#' 
#' @return Dataframe, contains 'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
#' 
#' @export
#' 
#' @import Seurat
#' @importFrom dplyr distinct
#'
#' @examples
#' \dontrun{
#' # Run CellChat with its default LR prior database
#' result <- CCCbank(ser, method = 'CellChat', species = 'human')
#' 
#' # Run CellChat with database of CellPhoneDBLR
#' result <- CCCbank(ser, method = 'CellChat', species = 'human', 
#'                   database = 'CellPhoneDBLR', extension =FALSE)
#' 
#' # Run CellChat with the combination of its default LR prior database 
#' # and database of CellPhoneDBLR
#' result <- CCCbank(ser, method = 'CellChat', species = 'human', 
#'                   database = 'CellPhoneDBLR', extension =TRUE)
#' 
#' # Run CellChat with the combination of its default LR prior database 
#' # and databases of CellPhoneDB and FANTOM5
#' result <- CCCbank(ser, method = 'CellChat', species = 'human', 
#'                   database = c('CellPhoneDBLR', 'FAMTOM5'), extension =TRUE)
#' 
#' # Run CellChat with user-defined database
#' database <- data.frame(ligand = c('A', 'B', 'C'), 
#'                        receptor = c('D', 'E', 'F'))
#' result <- CCCbank(ser, method = 'CellChat', species = 'human', 
#'                   database = database, extension =FALSE)
#' }
CCCbank <- function(ser, method, species = 'human', databases=NULL, extension =FALSE, keep_complexes=TRUE, ...){
  
  CheckSeuratObject(ser)
  
  methods_available <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", "Connectome", "ICELLNET", "NATMI", "iTALK",
                         "scConnect", "SingleCellSignalR", "CellChat", "RNAMagnet", 
                         "scSeqComm", "NicheNet", "CytoTalk", "scMLnet", "Domino")
  databases_available <- c("CellPhoneDBLR", "CellTalkerDB", "FANTOM5", "ICELLNETDB", "NATMIDB", "scConnectDB", 
                           "SingleCellSignalRDB", "CellChatDB", "RNAMagnetDB", "NicheNetDB", 
                           "CytoTalkDB", "DominoDB", 'CellCallDB')
  
  if(!(method %in% methods_available)){
    stop(paste0('Sorry, ', method, ' is currently not supported!'))
  }else{
    RCCCFunctionName <- paste0('Run', ifelse(method=='SingleCellSignalR', 'SCSR', method))
    RCCCFunction <- get(RCCCFunctionName)
  }
  
  if (is.character(databases)) {
    databases <- lapply(databases, function(database){
      if(!(database %in% databases_available)){
        stop(paste0('Sorry, the database ', database, ' is currently not supported!'))
      }else{
        database_path <- file.path(system.file("extdata", package = "CCCbank"), 'LRPriorDatabases',paste0(database, '.rds'))
        database <- readRDS(database_path)
        database <- database[, c('ligand', 'receptor', 'lr')]
      }
      database
    })
    databases <- do.call(rbind, databases)
    databases <- as.data.frame(databases)
    databases$ligand <- as.character(databases$ligand)
    databases$receptor <- as.character(databases$receptor)
    databases$lr <- as.character(databases$lr)
    databases <- dplyr::distinct(databases)
  }else if(is.data.frame(databases)){
    if(all(c('ligand', 'receptor') %in% colnames(databases))){
      databases$ligand <- as.character(databases$ligand)
      databases$receptor <- as.character(databases$receptor)
      databases$lr <- paste(databases$ligand, databases$receptor)
      databases <- dplyr::distinct(databases)
    }else{
      stop("Comfirm that the provided database contains 'ligand', and 'receptor' colunms!")
    }
  }
  
  # Change LR prior database of method
  if(!is.null(databases)){
    if(grepl('CellPhoneDB', method)){
      ChangeDBFunctionName <- 'ChangeCellPhoneDB'
    }else if(method == 'SingleCellSignalR'){
      ChangeDBFunctionName <- 'ChangeSCSRDB'
    }else{
      ChangeDBFunctionName <- paste0('Change', method, 'DB')
    }
    
    ChangeDBFunction <- get(ChangeDBFunctionName)
    
    message(paste0(Sys.time(), ' Changing LR prior database...'))
    
    if(method %in% c('scSeqComm', 'Domino')){
      message(paste0("The parameters 'extension' and 'keep_complexes' have no effect in ", method))
      db.info <- ChangeDBFunction(priorDatabase=databases)
    }else if(method %in% c('ICELLNET', 'RNAMagnet', 'CellPhoneDB2', 'CellPhoneDB3')){
      message(paste0("The parameter 'keep_complexes' have no effect in ", method))
      db.info <- ChangeDBFunction(priorDatabase=databases, extension = extension)
    }else if(method %in% c('iTALK', 'CellTalker', 'SingleCellSignalR',  'NicheNet', 'scMLnet', 'NATMI')){
      db.info <- ChangeDBFunction(priorDatabase=databases, extension = extension, keep_complexes = keep_complexes)
    }else if(method == 'CellChat'){
      message(paste0("The parameter 'keep_complexes' have no effect in ", method))
      db.info <- ChangeDBFunction(priorDatabase=databases, extension = extension, species = species)
    }else if(method %in% c('Connectome', 'CytoTalk', 'scConnect')){
      db.info <- ChangeDBFunction(priorDatabase=databases, extension = extension, keep_complexes = keep_complexes,
                             species = species)
    }
    
  }else{
    message(paste0(Sys.time(), ' Using the default LR prior database of ', method, '...'))
  }
  
  message(paste0(Sys.time(), ' Runnnig ', method, '...'))
  message(paste0(' More details on other parameters can be found in the ', RCCCFunctionName, ' function...'))
  
  if(exists('db.info')){
    new_db_info <- get('db.info')
  }else{
    new_db_info <- NULL
  }
  
  if(method %in% c('iTALK', 'CellTalker', 'RNAMagnet', 'ICELLNET', 'NicheNet')){
    result <- RCCCFunction(ser = ser, priorDatabase = new_db_info, ...)
  }else if(method == 'scMLnet'){
    result <- RCCCFunction(ser = ser, LigRecLib = new_db_info, ...)
  }else if(method %in% c('SingleCellSignalR', 'CellChat', 'Connectome', 'CytoTalk','scSeqComm')){
    result <- RCCCFunction(ser = ser, species = species, priorDatabase = new_db_info, ...)
  }else if(method %in% c('Domino', 'CellPhoneDB2', 'CellPhoneDB3') ){
    result <- RCCCFunction(ser = ser, db_path = new_db_info,...)
  }else if(method == 'scConnect'){
    result <- RCCCFunction(ser = ser, species = species, PyHome = new_db_info,...)
  }else if(method == 'NATMI'){
    result <- RCCCFunction(ser = ser, species = species, 
                           priorDatabase = ifelse(exists('db.info'), TRUE, FALSE),
                           NATMI_dir = ifelse(exists('db.info'), new_db_info, NA), ...)
  }
  message(paste0(Sys.time(), ' Finish!\nThanks for using CCCbank. Welcome to use again!'))
  return(result)
}
