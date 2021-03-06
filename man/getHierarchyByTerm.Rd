% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruning.R
\name{getHierarchyByTerm}
\alias{getHierarchyByTerm}
\title{getHierarchyByTerm}
\usage{
getHierarchyByTerm(rootTerm = NA, ontology = NA, interactionType = "is_a")
}
\arguments{
\item{rootTerm}{the term on de root hierachy, this term is used to pruning GO hierarchy and all hierarchy above it is returned}

\item{ontology}{the ontology of GO. select one of this tree options: BP, MF or CC (Biological Process, Molecular Function or Cellular Component)}

\item{interactionType}{the relations between terms, default is 'is_a'}
}
\value{
a list of the hierarchy pruned
}
\description{
this function returns de hierarchy go based on params
}
\examples{
getHierarchyByTerm(rootTerm = 'GO:0008501', ontology = 'BP')

}
