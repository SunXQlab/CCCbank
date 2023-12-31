% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunCCCFunction.R
\name{RunNicheNet}
\alias{RunNicheNet}
\title{Run NicheNet method to infer cell-cell communication (CCC).}
\usage{
RunNicheNet(
  ser,
  priorDatabase = NULL,
  lr = TRUE,
  n_ligands = 20,
  n_targets = 200
)
}
\arguments{
\item{ser}{Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'}

\item{priorDatabase}{If NULL, run with default LR prior database}

\item{lr}{If TRUE, result is the dataframe of L-R links.Otherwise, result is the dataframe of L-Target links. Defaults to TRUE}

\item{n_ligands}{Number of top-ranked ligands, which are further used to predict L-R links or L-Target links. Defaults to 20}

\item{n_targets}{Number of top-predicted target genes of ligands. If lr is FALSE, the parameter is required. Defaults to 200}
}
\value{
Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms, if lr is TRUE. Otherwise, it contains  'Ligand', 'Target', 'Sender', 'Receiver' and 'LRscore' five colunms
}
\description{
This function is to run NicheNet method to infer CCC.
}
\examples{
\dontrun{
# Run with default LR prior database
result <- RunNicheNet(ser)

# Run with other LR prior database
db <- ChangeNicheNetDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
result <- RunNicheNet(ser, priorDatabase = db)

# Run with default and other LR prior databases
db <- ChangeNicheNetDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
result <- RunNicheNet(ser, priorDatabase = db)
}
}
\author{
Jiaxin Luo
}
