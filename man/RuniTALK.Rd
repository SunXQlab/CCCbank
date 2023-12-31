% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunCCCFunction.R
\name{RuniTALK}
\alias{RuniTALK}
\title{Run iTALK method to infer cell-cell communication (CCC).}
\usage{
RuniTALK(ser, priorDatabase = NULL, ...)
}
\arguments{
\item{ser}{Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'}

\item{priorDatabase}{If NULL, run with default LR prior database}

\item{...}{other parameters in \code{\link[iTALK]{rawParse}}, except for the following parameter: data}
}
\value{
Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
}
\description{
This function is to run iTALK method to infer CCC.
}
\examples{
\dontrun{
# Run with default LR prior database
result <- RuniTALK(ser)

# Run with other LR prior database
db <- ChangeiTALKDB(priorDatabase=priorDatabase, extension=FALSE, keep_complexes = TRUE)
result <- RuniTALK(ser, priorDatabase = db)

# Run with default and other LR prior databases
db <- ChangeiTALKDB(priorDatabase=priorDatabase, extension=TRUE, keep_complexes = TRUE)
result <- RuniTALK(ser, priorDatabase = db)
}
}
\seealso{
\code{\link[iTALK]{rawParse}}
}
\author{
Jiaxin Luo
}
