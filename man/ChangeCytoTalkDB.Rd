% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChangeDBFunction.R
\name{ChangeCytoTalkDB}
\alias{ChangeCytoTalkDB}
\title{Transform the database into a format suitable for Connectome}
\usage{
ChangeCytoTalkDB(
  priorDatabase = NULL,
  extension = TRUE,
  keep_complexes = TRUE,
  species = NULL
)
}
\arguments{
\item{priorDatabase}{Dataframe with 'ligand', 'receptor' and 'lr' three columns}

\item{extension}{Logic. Whether use the combination of default database of Connectome and other databases}

\item{keep_complexes}{Logic. If TRUE, multi-subunit of ligand and receptor will be splited into single-subunit. Otherwise, the interactions containing multi-subunit will be filtered out.}

\item{species}{human and mouse}
}
\value{
list
}
\description{
Transform the database into a format suitable for Connectome
}
\author{
Jiaxin Luo
}
