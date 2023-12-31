% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunCCCFunction.R
\name{RunICELLNET}
\alias{RunICELLNET}
\title{Run ICELLNET method to infer cell-cell communication (CCC).}
\usage{
RunICELLNET(ser, priorDatabase = NULL, filter.perc = NULL)
}
\arguments{
\item{ser}{Seurat object, contains normalized data stored in 'data' slot of 'RNA' assay, and cell type metadata stored in meta.data as 'celltype'}

\item{priorDatabase}{If NULL, run with default LR prior database}

\item{filter.perc}{see details in \code{\link[icellnet]{sc.data.cleaning}}}
}
\value{
Dataframe, contains  'Ligand', 'Receptor', 'Sender', 'Receiver' and 'LRscore' five colunms.
}
\description{
This function is to run ICELLNET method to infer CCC.
}
\examples{
\dontrun{
# Run with default LR prior database
result <- RunICELLNET(ser)

# Run with other LR prior database
db <- ChangeICELLNETDB(priorDatabase=priorDatabase, extension=FALSE)
result <- RunICELLNET(ser, priorDatabase = db)

# Run with default and other LR prior databases
db <- ChangeICELLNETDB(priorDatabase=priorDatabase, extension=TRUE)
result <- RunICELLNET(ser, priorDatabase = db)
}
}
\seealso{
\code{\link[icellnet]{sc.data.cleaning}}
}
\author{
Jiaxin Luo
}
