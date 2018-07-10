#' @title getHierarchyByTerm
#'
#' @description this function returns de hierarchy go based on params
#'
#' @param rootTerm the term on de root hierachy, this term is used to pruning GO hierarchy and all hierarchy above it is returned
#' @param ontology the ontology of GO. select one of this tree options: BP, MF or CC (Biological Process, Molecular Function or Cellular Component)
#' @param interactionType the relations between terms, default is 'is_a'
#'
#' @return GO hierarchy with rootTerm at the highest level
#'
#' @examples getHierarchyByTerm(rootTerm = 'GO:0008501', ontology = 'BP')
#'
#' @export getHierarchyByTerm
getHierarchyByTerm <- function(rootTerm = NA, ontology = NA, interactionType = 'is_a'){

  if(is.na(ontology))
    stop("ontology cannot be NA")

  switch (ontology,
    "BP" = {ontology = GOBPCHILDREN},
    "MF" = {ontology = GOMFCHILDREN},
    "CC" = {ontology = GOCCCHILDREN}
  )

  # if(is.na(ontology))
  #   stop("ontology must be MF, BP or CC")

  allTerms<-as.list(ontology)
  termNames <- names(allTerms)
  initialHierarchy <- list()

  if(!is.na(rootTerm)){
    #positon of initialTerm
    initialTermId <- which(termNames%in%rootTerm)

    if(length(initialTermId) >0){
      #childrenTerms of initialTerm
      childrenTerms <- allTerms[[initialTermId]]
      initialHierarchy[[rootTerm]] <- unname(childrenTerms[which(names(childrenTerms) == interactionType)], force = FALSE)
    }else{
      stop("initialTerm not find in ontology")
    }
  }else{
    stop("initialTerm cannot be NA")
  }
  message("** loading hierarchy \r")

  hierarchy = getHierarchy(hierarchy = initialHierarchy,
                           initialHierarchy = initialHierarchy[[1]],
                           allTerms = allTerms,
                           interactionType = interactionType)

  return(hierarchy)

}


#receive go hierarchy object and the last child on hierchy
#call recursivly until there is no have any children
getHierarchy <- function(hierarchy = NA, initialHierarchy = NA, allTerms = NA, interactionType = "is_a"){
  termNames <- names(allTerms)

  message(paste(rep(".", length(hierarchy))),"\r",appendLF=FALSE)

  for(i in 1:length(initialHierarchy)){

    term <- unname(initialHierarchy[i])
    goid <- which(termNames%in%term)
    #if term have any children
    if(length(goid) >0){
      childrenTerms <- allTerms[[term]]
      childrenTerms <- unname(childrenTerms[which(names(childrenTerms) == interactionType)], force = FALSE)
      if(length(childrenTerms) > 0){
        hierarchy[[term]] <- childrenTerms
        hierarchy <- getHierarchy(hierarchy = hierarchy, initialHierarchy = hierarchy[[term]], allTerms= allTerms, interactionType=interactionType)
      }
    }
  }

  return(hierarchy)

}

#filter terms on hierarchy with genes associations
#remove term with no children
#remove genes with no terms associations
#' @title filterTermByGenes
#'
#' @description function that prunes the hierarchy of the ontology gene according to the annotation of the past genes.  removing terms without annotated genes and genes that are not annotated in any term present in the hierarchy.
#'
#' @param terms the terms hierarchy, the result of getHierarchyByTerm
#' @param genes list of genes with their terms
#'
#' @returna list named by the terms and their respective genes
#'
#' @export getHierarchyByTerm
filterTermByGenes  <- function(terms = NA, genes = NA){

  listGene2GO <- list()
  termsVector <- unname(unlist(terms))
  termsVector <- c(names(terms)[1], termsVector)
  #tst <- stack(genes)

  listGene2GO <- genes[(names(genes) %in% termsVector)]
  return(listGene2GO)
}

