RNAMagnetSignaling <- function (seurat, neighborhood.distance = NULL, neighborhood.gradient = NULL,
                                .k = 10, .x0 = 0.5, .minExpression = 10, .minMolecules = 1,
                                .version = "1.0.0", .cellularCompartment = c("Secreted", "Both"), 
                                .manualAnnotation = "Correct") 
{
  RNAMagnetBase(seurat, anchors = NULL, neighborhood.distance, 
                neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, 
                .version, .cellularCompartment, .manualAnnotation, FALSE)
}

#' @import Seurat
#' @importFrom utils install.packages
#' @importFrom devtools install_github
#' @importFrom methods new
#' @importFrom stats cor
RNAMagnetBase <- function (seurat, anchors = NULL, neighborhood.distance = NULL,
                           neighborhood.gradient = NULL, .k = 10, .x0 = 0.5, .minExpression,
                           .minMolecules = 1, .version = "1.0.0", .cellularCompartment,
                           .manualAnnotation = "Correct", .symmetric = F) 
{
  cat("Setting everything up...\n")
  
  if(!requireNamespace('Biobase')){
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    
    BiocManager::install("Biobase")
  }
  
  if(!requireNamespace('Rmagic')){
    devtools::install_github('cran/Rmagic')
  }
  
  if (!grepl("^2", Biobase::package.version("Seurat"))) {
    seurat.ident <- Idents(seurat)
    seurat.cell.names <- colnames(seurat)
    seurat.pca <- Embeddings(seurat, reduction = "pca")
    seurat.raw.data <- GetAssayData(seurat, slot = "counts")
  }
  else {
    seurat.ident <- seurat@ident
    seurat.cell.names <- seurat@cell.names
    seurat.pca <- seurat@dr$pca@cell.embeddings
    seurat.raw.data <- seurat@raw.data
  }
  if (is.null(anchors)) 
    anchors <- as.character(unique(seurat.ident))
  out <- new("rnamagnet", celltype = seurat.ident, params = list(neighborhood.distance = neighborhood.distance, 
                                                                 neighborhood.gradient = neighborhood.gradient, .k = .k, 
                                                                 .x0 = .x0, .minExpression = .minExpression, .minMolecules = .minMolecules, 
                                                                 .cellularCompartment = .cellularCompartment, .manualAnnotation = .manualAnnotation, 
                                                                 .symmetric = .symmetric))
  similarity <- as.matrix(1 - cor(t(seurat.pca[, 1:15])))
  ligrec <- getLigandsReceptors(.version, .cellularCompartment, .manualAnnotation)
  if (.symmetric) 
    ligrec <- makeSymmetric(ligrec)
  filteredGenes <- rownames(seurat.raw.data)[apply(seurat.raw.data[, seurat.cell.names] >= .minMolecules, 1, sum) > .minExpression]
  genes <- unique(c(ligrec$Receptor.Mouse, ligrec$Ligand.Mouse))
  genes <- genes[sapply(genes, function(x) {
    entries <- strsplit(x, "[&|]")
    if (grepl("&", x)) 
      all(entries[[1]] %in% filteredGenes)
    else any(entries[[1]] %in% filteredGenes)
  })]
  genes_formagic <- unlist(strsplit(genes, "[&|]"))
  formagic <- Matrix::t(seurat.raw.data[, seurat.cell.names])
  formagic <- Rmagic::library.size.normalize(formagic)
  formagic <- sqrt(formagic)
  cat("Now running MAGIC to impute dropout values...\n")
  mymagic <- Rmagic::magic(formagic, genes = genes_formagic, seed = 48879)
  mymagic <- as.data.frame(mymagic)
  resolvedRawData <- resolveData(t(mymagic), genes)
  resolvedRawData <- t(apply(resolvedRawData, 1, function(x) (x - min(x))/(max(x) - min(x))))
  annotated_genes <- rownames(resolvedRawData)
  out@mylr <- subset(ligrec, Receptor.Mouse %in% annotated_genes & Ligand.Mouse %in% annotated_genes)
  stepf <- Vectorize(function(x) if (x < 0)
    0
    else x)
  cat("Now running RNAMagnet...\n")
  out@anchors <- do.call(cbind, lapply(anchors, function(id) {
    apply(resolvedRawData[, seurat.ident == id], 1, mean)
  }))
  colnames(out@anchors) <- anchors
  out@interaction <- sapply(anchors, function(pop_l) {
    out@mylr$expression_ligand <- out@anchors[out@mylr$Ligand.Mouse, pop_l]
    sapply(seurat.cell.names, function(cell_r) {
      out@mylr$expression_receptor <- resolvedRawData[out@mylr$Receptor.Mouse, cell_r]
      sum(kernel(out@mylr$expression_ligand, k = .k, x0 = .x0) * kernel(out@mylr$expression_receptor, k = .k, x0 = .x0))
    })
  })
  out@specificity <- t(sapply(rownames(out@interaction), function(cell) {
    x <- out@interaction[cell, ]
    beta <- x/sum(x)
    if (!is.null(neighborhood.distance)) 
      alpha <- apply(out@interaction * (1 - kernel(similarity[cell, ], neighborhood.gradient, x0 = neighborhood.distance)), 2, sum)
    else alpha <- apply(out@interaction, 2, mean)
    alpha <- alpha/sum(alpha)
    beta - alpha
  }))
  rownames(out@specificity) <- rownames(out@interaction)
  out@adhesiveness <- apply(mymagic, 1, function(x) sum(kernel(x, x0 = .x0, k = .k)))
  return(out)
}

getLigandsReceptors <- function (version = "latest", cellularCompartment = c("Membrane", "ECM", "Both", "Secreted"), 
                                 manualAnnotation = "Correct",
                                 ligandClass = c("Other", "Cytokine", "Chemokine", "GrowthFactor", "Interleukin")) 
{
  if (is.character(version)) 
    out <- switch(version, latest = ligandsReceptors_2.2.0, 
                  stable = ligandsReceptors_1.0.0, `1.0.0` = ligandsReceptors_1.0.0, 
                  `2.0.0` = ligandsReceptors_2.0.0, `2.1.0` = ligandsReceptors_2.1.0, 
                  `2.2.0` = ligandsReceptors_2.2.0, `3.0.0` = ligandsReceptors_3.0.0)
  else if (is.data.frame(version)) 
    out <- version
  else stop("Version needs to be a character string or a data frame.")
  if (is.null(out)) 
    error("Version", version, "not supported. See documentation.")
  out <- subset(out, Ligand.CC %in% cellularCompartment & 
                  ManualAnnotation %in% manualAnnotation & Ligand.GO %in% 
                  ligandClass)
  out
}

makeSymmetric <- function (ligrec) 
{
  toadd <- list()
  for (i in 1:nrow(ligrec)) {
    if (!any(ligrec$Ligand.Mouse[i] == ligrec$Receptor.Mouse & 
             ligrec$Receptor.Mouse[i] == ligrec$Ligand.Mouse)) {
      toadd[[length(toadd) + 1]] <- ligrec[i, ]
      toadd[[length(toadd)]]$Receptor.Mouse <- ligrec[i, 
                                                      "Ligand.Mouse"]
      toadd[[length(toadd)]]$Ligand.Mouse <- ligrec[i, 
                                                    "Receptor.Mouse"]
    }
  }
  rbind(ligrec, do.call(rbind, toadd))
}

kernel <- function (x, k = 10, x0 = 0.5)
  1/(1 + exp(-k * (x - x0)))

resolveData <- function (rawdata, lr) 
{
  singleentries <- lr[!grepl("[&|]", lr)]
  use_normdata <- rawdata[singleentries, ]
  doubleentries <- lr[grepl("[&|]", lr)]
  if(length(doubleentries)>0){
    add_normdata <- do.call(rbind, lapply(doubleentries, function(x) {
      if (grepl("&", x)) {
        entries <- strsplit(x, "&", fixed = T)[[1]]
        apply(rawdata[entries, ], 2, min)
      }
      else {
        entries <- strsplit(x, "|", fixed = T)
        entries <- entries[[1]][entries[[1]] %in% rownames(rawdata)]
        if (length(entries) > 1) 
          apply(rawdata[entries, ], 2, max)
        else rawdata[entries, ]
      }
    }))
    rownames(add_normdata) <- doubleentries
    use_normdata <- rbind(use_normdata, add_normdata)
  }
  
  use_normdata[lr, ]
  
}
