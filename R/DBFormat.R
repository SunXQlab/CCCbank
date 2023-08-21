iTALKDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- distinct(db2)
  }
  
  colnames(db2) <- c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'Pair.Name')
  db2$Classification <- 'other'
  db2$Ligand.Name <- db2$Ligand.ApprovedSymbol
  db2$Receptor.Name <- db2$Receptor.ApprovedSymbol
  
  db2 <- db2[, c('Pair.Name', 'Ligand.ApprovedSymbol', 'Ligand.Name', 
                 'Receptor.ApprovedSymbol', 'Receptor.Name', 'Classification')]
  return(db2)
}

CellTalkerDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  colnames(db2) <- c('ligand', 'receptor', 'pair')
  
  db2 <- db2[, c('ligand', 'receptor', 'pair')]
  return(db2)
}

RNAMagnetDBFormat <- function(db2){
  colnames(db2)[1:2] <- c("Ligand.Mouse", "Receptor.Mouse")
  db2$Pair.Name <- paste(db2$Ligand.Mouse, db2$Receptor.Mouse, sep = '-')
  db2$lr <- NULL
  
  db2$Source <- 'user-defined'
  db2$ManualAnnotation <- 'Correct'
  db2$Ligand.CC <- 'Both'
  db2$Ligand.GO <- 'Other'
  db2$Reference <- NA
  
  db2 <- db2[, c("Pair.Name", "Ligand.Mouse", "Receptor.Mouse", 
                 "Source", "ManualAnnotation", "Ligand.CC", 
                 "Ligand.GO", "Reference")]
  
  return(db2)
}

ICELLNETDBFormat <- function(db2){
  db2 <- suppressWarnings(tidyr::separate(db2, receptor, c('Receptor 1', 'Receptor 2', 'Receptor 3'), '&'))
  db2 <- suppressWarnings(tidyr::separate(db2, ligand, c('Ligand 1', 'Ligand 2'), '&'))
  db2$lr <- NULL
  db2$Alias <- NA
  db2$Family <- NA
  db2$Cytokines <- NA
  db2$Checkpoints <- NA
  db2$Chemokines <- NA
  db2$Chemokines <- NA
  db2$Other <- NA
  db2$Subfamily <- NA
  db2$Classifications <- 'unknown'
  db2$`PubMed ID` <- NA
  db2$Comments <- NA
  return(db2)
}

SCSRDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  db2$source <- 'user-defined'
  db2$PMIDs <- ''
  db2$lr <- NULL
  
  db2 <- db2[, c('ligand', 'receptor', 'source', 'PMIDs')]
  return(db2)
}

CellChatDBFormat <- function(db2, species=NULL){
  #species: support human, mouse and zebrafish
  if (!requireNamespace('CellChat')) {
    message("CellChat hasn't been installed. Now installing...")
    devtools::install_github("sqjin/CellChat")
  }
  
  if(is.null(species)){
    warning("No 'species' provied, 'species' is set to 'human'")
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
  
  origin_geneInfo <- origin_db$geneInfo
  
  gene_add <- unique(c(db2$ligand, db2$receptor))
  gene_add <- unique(unlist(stringr::str_split(gene_add, '&')))
  gene_add <- setdiff(gene_add, origin_geneInfo$Symbol)
  if(length(gene_add)>0){
    gene_add <- data.frame(Symbol = gene_add,
                           Name = gene_add,
                           EntrezGene.ID = NA,
                           Ensembl.Gene.ID = '',
                           MGI.ID = '',
                           Gene.group.name = '')
    if(species=='human'){
      label <- seq(from = 60000, by = 1, length=dim(gene_add)[[1]])
    }else{
      label <- seq(from = 7000000, by = 1, length=dim(gene_add)[[1]])
    }
    rownames(gene_add) <- paste0(ifelse(species=='human', 'HGNC:', 'MGI:'), label)
  }else{
    gene_add <- NULL
  }
  
  
  db2$lr <- NULL
  db2$ligand <- gsub('&', '_', db2$ligand)
  db2$receptor <- gsub('&', '_', db2$receptor)
  db2$interaction_name <- paste(db2$ligand, db2$receptor, sep = '_')
  db2$lig <- gsub('_', '+', db2$ligand)
  db2$rec <- gsub('_', '+', db2$receptor)
  db2$lig[which(grepl('\\+', db2$lig))] <- paste0('(', db2$lig[which(grepl('\\+', db2$lig))], ')')
  db2$rec[which(grepl('\\+', db2$rec))] <- paste0('(', db2$rec[which(grepl('\\+', db2$rec))], ')')
  db2$interaction_name_2 <- paste(db2$lig, db2$rec, sep = ' - ')
  db2$lig <- NULL; db2$rec <- NULL
  
  db2$pathway_name <- 'Unknown'
  db2$agonist <- ''
  db2$antagonist <- ''
  db2$co_A_receptor <- ''
  db2$co_I_receptor <- ''
  db2$evidence <- 'user-defined'
  db2$annotation <- 'user-defined'
  db2 <- db2[, c("interaction_name", "pathway_name", "ligand", "receptor",
                 "agonist", "antagonist", "co_A_receptor", "co_I_receptor",
                 "evidence", "annotation", "interaction_name_2")]
  complexes_add <- c(db2$ligand[which(grepl('_', db2$ligand))],
                     db2$receptor[which(grepl('_', db2$receptor))])
  complexes_add <- unique(complexes_add)
  if(length(complexes_add)>0){
    complexes_add <- data.frame(complexes = complexes_add)
    complexes_add <- suppressWarnings(tidyr::separate(complexes_add, 'complexes',
                                                      c("subunit_1", "subunit_2", "subunit_3", "subunit_4"),
                                                      '_', remove = FALSE))
    complexes_add$subunit_3[is.na(complexes_add$subunit_3)] <- ''
    complexes_add$subunit_4[is.na(complexes_add$subunit_4)] <- ''
    complexes_add <- tibble::column_to_rownames(complexes_add, 'complexes')
  }else{
    complexes_add <- NULL
  }
  
  
  return(list(interaction = db2,
              complex = complexes_add,
              geneInfo = gene_add))
}

ConnectomeDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  db2$lr <- NULL
  colnames(db2) <- c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol')
  db2$mode <- 'UNCAT'
  
  db2 <- db2[, c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol", "mode" )]
  return(db2)
}

CytoTalkDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  db2$lr <- NULL
  
  db2 <- db2[, c("ligand", "receptor")]
  return(db2)
}

NATMIDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  db2$lr <- NULL
  colnames(db2) <- c('Ligand gene symbol', 'Receptor gene symbol')
  
  db2 <- db2[,c('Ligand gene symbol', 'Receptor gene symbol')]
  return(db2)
}

DominoDBFormat <- function(db, EnsgToSymbol, UprotToSymbol){
  
  genes_tmp <- unique(c(db$ligand, db$receptor))
  
  genes_temp <- unique(unlist(str_split(genes_tmp, '&')))
  db_genes <- data.frame(gene_name = genes_temp,
                         uniprot = UprotToSymbol$Entry[match(genes_temp, UprotToSymbol$`Gene Names`)],
                         hgnc_symbol = genes_temp,
                         ensembl = EnsgToSymbol$Gene.stable.ID[match(genes_temp, EnsgToSymbol$Gene.name)])
  db_genes$uniprot[is.na(db_genes$uniprot)] <- paste0('Uniprot_Unknown', seq(length(which(is.na(db_genes$uniprot)))))
  db_genes$ensembl[is.na(db_genes$ensembl)] <- paste0('ENSG_Unknown', seq(length(which(is.na(db_genes$ensembl)))))
  
  db_complexes <- genes_tmp[grepl('&',genes_tmp)]
  if(length(db_complexes)>0){
    
    genes_temp <- unique(unlist(str_split(db_complexes, '&')))
    Uprot_temp <- db_genes[which(db_genes$hgnc_symbol %in% genes_temp), ]
    Uprot_temp <- Uprot_temp[, 1:2]
    colnames(Uprot_temp) <- c('genes', 'uniprot')
    
    db_complexes <- suppressWarnings(tidyr::separate(data.frame(complex_name = db_complexes), 'complex_name',
                                              c('gene_1', 'gene_2', 'gene_3',	'gene_4'),
                                              sep = '&', remove = FALSE))
    db_complexes <- merge(db_complexes, Uprot_temp, by.x = 'gene_1', by.y = 'genes', all.x = TRUE)
    db_complexes <- merge(db_complexes, Uprot_temp, by.x = 'gene_2', by.y = 'genes', all.x = TRUE)
    db_complexes <- merge(db_complexes, Uprot_temp, by.x = 'gene_3', by.y = 'genes', all.x = TRUE)
    db_complexes <- suppressWarnings(merge(db_complexes, Uprot_temp, by.x = 'gene_4', by.y = 'genes', all.x = TRUE))
    db_complexes <- db_complexes[, -c(1:4)]
    colnames(db_complexes)[2:5] <- c('uniprot_1', 'uniprot_2', 'uniprot_3',	'uniprot_4')
    db_complexes$uniprot_3[which(is.na(db_complexes$uniprot_3))] <- ''
    db_complexes$uniprot_4[which(is.na(db_complexes$uniprot_4))] <- ''
    db_complexes$receptor <- ifelse(db_complexes$complex_name %in% db$receptor, TRUE, FALSE)
    
  }else{
    db_complexes <- data.frame(complex_name ='', uniprot_1='', uniprot_2='', 
                               uniprot_3='', uniprot_4='', receptor='')
  }
  
  
  db_protein <- data.frame(uniprot = db_genes$uniprot ,
                           protein_name = db_genes$hgnc_symbol)
  db_protein$transmembrane <- ''
  db_protein$peripheral <- ''
  db_protein$secreted <- ''
  db_protein$secreted_desc <- ''
  db_protein$secreted_highlight <- ''
  db_protein$receptor <- ''
  db_protein$receptor <- ifelse(db_protein$protein_name %in% db$receptor, 'True', 'False')
  db_protein$receptor_desc <- ''
  db_protein$integrin <- ''
  db_protein$other <- ''
  db_protein$other_desc <- ''
  db_protein$tags <- ''
  db_protein$tags_description <- ''
  db_protein$tags_reason <- ''
  db_protein$pfam <- ''
  db_protein$protein_name <- paste0(db_protein$protein_name, '_HUMAN')
  
  test_genes <- db_genes[, c('uniprot','gene_name')]
  colnames(test_genes) <- c('partner', 'gene')
  
  db_interaction <- db
  db_interaction$lr <- NULL
  db_interaction <- merge(db_interaction, test_genes, 
                          by.x = 'ligand', by.y = 'gene', all.x = TRUE)
  db_interaction <- merge(db_interaction, test_genes, 
                          by.x = 'receptor', by.y = 'gene', all.x = TRUE)
  colnames(db_interaction)[3:4] <- c('partner_a', 'partner_b')
  db_interaction$protein_name_a[!is.na(db_interaction$partner_a)] <- paste0(db_interaction$ligand[!is.na(db_interaction$partner_a)], '_HUMAN')
  db_interaction$protein_name_b[!is.na(db_interaction$partner_b)] <- paste0(db_interaction$receptor[!is.na(db_interaction$partner_b)], '_HUMAN')
  db_interaction$protein_name_a[is.na(db_interaction$protein_name_a)] <- ''
  db_interaction$protein_name_b[is.na(db_interaction$protein_name_b)] <- ''
  
  
  db_interaction$partner_a[is.na(db_interaction$partner_a)] <- db_interaction$ligand[is.na(db_interaction$partner_a)]
  db_interaction$partner_b[is.na(db_interaction$partner_b)] <- db_interaction$receptor[is.na(db_interaction$partner_b)]
  
  db_interaction$annotation_strategy <- 'user-defined'
  db_interaction$source <- 'user-defined'
  db_interaction$id_cp_interaction <- paste0('CPI-User', seq(dim(db_interaction)[[1]]))
  
  db_interaction <- db_interaction[,-c(1:2)]
  db_interaction <- db_interaction[, c(7, 1:6)]
  
  return(list(interaction = db_interaction,
              complex = db_complexes,
              genes = db_genes,
              protein = db_protein))
}

scConnectDBFormat <- function(db2, keep_complexes, ligands_db1, receptors_db1){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  db2$lr <- NULL
  db1_ligands <- ligands_db1
  db1_ligands$X <- NULL
  db1_receptors <- receptors_db1
  db1_receptors$X <- NULL
  
  ligands.diff <- setdiff(unique(db2$ligand), unique(db1_ligands$ligand_symbol))
  receptor.diff <- setdiff(unique(db2$receptor), unique(db1_receptors$receptor_symbol))
  
  if(length(ligands.diff)>0){
    ligands.add <- data.frame(ligand = ligands.diff, 
                              ligand_type = rep('peptide', length(ligands.diff)),
                              preprogene = paste0("['", ligands.diff, "']"),
                              synthesis = rep('', length(ligands.diff)),
                              transport = rep('', length(ligands.diff)),
                              reuptake = rep('', length(ligands.diff)),
                              excluded = rep('', length(ligands.diff)),
                              comment = rep('user-defined', length(ligands.diff)),
                              ligand_symbol = ligands.diff
    )
    ligands_db <- rbind(db1_ligands, ligands.add)
  }else{
    ligands_db <- db1_ligands
  }
  
  if (length(receptor.diff)>0) {
    receptors.add <- data.frame(receptor = receptor.diff,
                                family = rep('unknown', length(receptor.diff)),
                                type = rep('other_protein', length(receptor.diff)),
                                gene = paste0("['", receptor.diff, "']"),
                                receptor_symbol = receptor.diff
    )
    receptors_db <- rbind(db1_receptors, receptors.add)
  }else{
    receptors_db <- db1_receptors
  }
  
  
  db2 <- merge(db2, ligands_db, by.x = 'ligand', by.y = 'ligand_symbol')
  db2 <- db2[, -c(4:10)]
  db2 <- merge(db2, receptors_db, by.x = 'receptor', by.y = 'receptor_symbol')
  db2 <- db2[, c("ligand.y", "receptor.y")]
  colnames(db2) <- c("ligand", "target")
  db2$endogenous <- 'f'
  db2$ligand_species <- 'Human'
  db2$target_species <- 'Unknown'
  db2$pubmed_id <- 'Unknown'
  db2$action <- 'Unknown'
  
  return(list(db = db2,
              ligands = ligands_db,
              receptors = receptors_db))
}

#' @import dplyr
#' @importFrom utils packageVersion
NicheNetDBFormat <- function(db2, sig_network, gr_network){
  
  if(!requireNamespace('nichenetr')){
    message("nichenetr hasn't been installed. Now installing...")
    devtools::install_github("saeyslab/nichenetr")
  }
  
  source_weights_df_new <- rbind(nichenetr::source_weights_df, 
                                 data.frame(source = 'user-defined', weight=1))
  weighted_networks_new = nichenetr::construct_weighted_networks(lr_network = db2,
                                                                 sig_network = sig_network,
                                                                 gr_network = gr_network,
                                                                 source_weights_df = source_weights_df_new)
  
  if(grepl('^1',packageVersion("nichenetr"))){
    weighted_networks_new = nichenetr::apply_hub_corrections(weighted_networks = weighted_networks_new,
                                                             lr_sig_hub = nichenetr::hyperparameter_list$lr_sig_hub,
                                                             gr_hub = nichenetr::hyperparameter_list$gr_hub)
    ligands = unique(db2$from) %>% as.list()
    ligand_target_matrix_new = nichenetr::construct_ligand_target_matrix(weighted_networks = weighted_networks_new, 
                                                                         ligands = ligands, algorithm = "PPR",
                                                                         damping_factor = nichenetr::hyperparameter_list$damping_factor,
                                                                         ltf_cutoff = nichenetr::hyperparameter_list$ltf_cutoff)
  }else if(grepl('^2',packageVersion("nichenetr"))){
    weighted_networks_new = nichenetr::apply_hub_corrections(weighted_networks = weighted_networks_new,
                                                             lr_sig_hub = nichenetr::hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                                             gr_hub = nichenetr::hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))
    ligands = unique(db2$from) %>% as.list()
    ligand_target_matrix_new = nichenetr::construct_ligand_target_matrix(weighted_networks = weighted_networks_new, 
                                                                         ligands = ligands, algorithm = "PPR",
                                                                         damping_factor = nichenetr::hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                                         ltf_cutoff = nichenetr::hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
  }
  
  
  
  return(list(weighted_networks = weighted_networks_new,
              ligand_target_matrix = ligand_target_matrix_new))
}

scMLnetDBFormat <- function(db2, keep_complexes){
  
  complex_lr <- db2[grepl('&', db2$lr), ]
  db2 <- db2[!grepl('&', db2$lr), ]
  
  if ((dim(complex_lr)[[1]]!=0) & keep_complexes) {
    complex_lr <- tidyr::separate_rows(complex_lr, 'ligand', sep = '&')
    complex_lr <- tidyr::separate_rows(complex_lr, 'receptor', sep = '&')
    complex_lr$lr <- paste(complex_lr$ligand, complex_lr$receptor, sep = '_')
    
    if(dim(db2)[[1]]!=0){
      db2 <- rbind(db2, complex_lr)
    }else{
      db2 <- complex_lr
    }
    db2 <- dplyr::distinct(db2)
  }
  
  colnames(db2) <- c('source', 'target', 'Key')
  
  db2 <- db2[, c('Key', 'source', 'target')]
  return(db2)
}

CellPhoneDBFormat <- function(db2, db1_gene, origin_protein, EnsgToSymbol, UprotToSymbol){
  
  genes_tmp <- unique(c(db2$ligand, db2$receptor))
  
  genes_temp <- genes_tmp[!grepl('&', genes_tmp)]
  genes.diff <- setdiff(genes_temp, db1_gene$gene)
  if(length(genes.diff)>0){
    db2_genes <- data.frame(gene_name = genes.diff,
                            uniprot = UprotToSymbol$Entry[match(genes.diff, UprotToSymbol$`Gene Names`)],
                            hgnc_symbol = genes.diff,
                            ensembl = EnsgToSymbol$Gene.stable.ID[match(genes.diff, EnsgToSymbol$Gene.name)])
    db2_genes$uniprot[is.na(db2_genes$uniprot)] <- paste0('Un', db2_genes$gene_name[which(is.na(db2_genes$uniprot))])
    db2_genes$ensembl[is.na(db2_genes$ensembl)] <- paste0('ENSG', db2_genes$gene_name[which(is.na(db2_genes$ensembl))])
  }else{
    db2_genes <- NULL
  }
  
  db2_complexes <- genes_tmp[grepl('&',genes_tmp)]
  complexes.diff <- setdiff(db2_complexes, db1_gene$gene)
  if(length(complexes.diff)>0){
    
    genes_temp <- unique(unlist(stringr::str_split(db2_complexes, '&')))
    Uprot_temp <- data.frame(genes = genes_temp,
                             uniprot = UprotToSymbol$Entry[match(genes_temp, UprotToSymbol$`Gene Names`)],
                             protein_name = UprotToSymbol$`Entry Name`[match(genes_temp, UprotToSymbol$`Gene Names`)])
    Uprot_temp$uniprot[is.na(Uprot_temp$uniprot)] <- Uprot_temp$genes[is.na(Uprot_temp$uniprot)]
    Uprot_temp$protein_name[is.na(Uprot_temp$protein_name)] <- Uprot_temp$genes[is.na(Uprot_temp$protein_name)]
    
    Uprot_tmp <- Uprot_temp[,1:2]
    db2_complexes <-suppressWarnings( tidyr::separate(data.frame(complex_name = db2_complexes), 'complex_name', 
                                                      c('gene_1', 'gene_2', 'gene_3',	'gene_4'), 
                                                      sep = '&', remove = FALSE))
    db2_complexes <- merge(db2_complexes, Uprot_tmp, by.x = 'gene_1', by.y = 'genes', all.x = TRUE)
    db2_complexes <- merge(db2_complexes, Uprot_tmp, by.x = 'gene_2', by.y = 'genes', all.x = TRUE)
    db2_complexes <- merge(db2_complexes, Uprot_tmp, by.x = 'gene_3', by.y = 'genes', all.x = TRUE)
    db2_complexes <- suppressWarnings(merge(db2_complexes, Uprot_tmp, by.x = 'gene_4', by.y = 'genes', all.x = TRUE))
    db2_complexes <- db2_complexes[, -c(1:4)]
    colnames(db2_complexes)[2:5] <- c('uniprot_1', 'uniprot_2', 'uniprot_3',	'uniprot_4')
  }else{
    db2_complexes <- NULL
  }
  
  if(!is.null(db2_complexes)){
    db2_protein <- unique(unlist(db2_complexes[, 2:5]))
    db2_protein <- db2_protein[which(!is.na(db2_protein))]
  }else{
    db2_protein <- c()
  }
  
  protein.diff <- setdiff(db2_protein, origin_protein$uniprot)
  if(length(protein.diff)>0){
    db2_protein <- Uprot_temp[which(Uprot_temp$uniprot %in% protein.diff), ]
    db2_protein <- db2_protein[, 2:3]
    db2_protein$transmembrane <- ''
    db2_protein$peripheral <- ''
    db2_protein$secreted <- ''
    db2_protein$secreted_desc <- ''
    db2_protein$secreted_highlight <- ''
    db2_protein$receptor <- ''
    db2_protein$receptor_desc <- ''
    db2_protein$integrin <- ''
    db2_protein$other <- ''
    db2_protein$other_desc <- ''
    db2_protein$tags <- ''
    db2_protein$tags_description <- ''
    db2_protein$tags_reason <- ''
    db2_protein$pfam <- ''
  }else{
    db2_protein <- NULL
  }
  
  test_genes <- db2_genes[, c('uniprot','gene_name')]
  colnames(test_genes) <- c('partner', 'gene')
  test_genes <- rbind(test_genes, db1_gene)
  test_genes <- dplyr::distinct(test_genes)
  
  db2_interaction <- db2
  db2_interaction$lr <- NULL
  db2_interaction <- merge(db2_interaction, test_genes, by.x = 'ligand', by.y = 'gene', all.x = TRUE)
  db2_interaction$partner[is.na(db2_interaction$partner)] <- db2_interaction$ligand[is.na(db2_interaction$partner)]
  colnames(db2_interaction)[3] <- 'partner_a'
  db2_interaction <- merge(db2_interaction, test_genes, by.x = 'receptor', by.y = 'gene', all.x = TRUE)
  db2_interaction$partner[is.na(db2_interaction$partner)] <- db2_interaction$receptor[is.na(db2_interaction$partner)]
  colnames(db2_interaction)[4] <- 'partner_b'
  db2_interaction$protein_name_a[!grepl('&', db2_interaction$ligand)] <- paste0(db2_interaction$ligand[!grepl('&', db2_interaction$ligand)], '_HUMAN')
  db2_interaction$protein_name_b[!grepl('&', db2_interaction$receptor)] <- paste0(db2_interaction$receptor[!grepl('&', db2_interaction$receptor)], '_HUMAN')
  db2_interaction$annotation_strategy <- 'user-defined'
  db2_interaction$source <- 'user-defined'
  db2_interaction$id_cp_interaction <- paste0('CPI-User', seq(dim(db2_interaction)[[1]]))
  db2_interaction <- db2_interaction[, -c(1:2)]
  
  db2_interaction <- db2_interaction[, c(7, 1:6)]
  
  return(list(interaction = db2_interaction,
              complex = db2_complexes,
              genes = db2_genes,
              protein = db2_protein))
}