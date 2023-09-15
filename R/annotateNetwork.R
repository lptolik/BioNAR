#TODO: convert all annotate_ functions to use either loopOverFiles, or
#      ontologies, or data.frames
#' Remove vertex property.
#'
#' @param GG igraph object
#' @param NAME name of the vertex property to remove
#'
#' @return igraph object with attribute removed
#' @export
#'
#' @examples
#' data(karate, package='igraphdata')
#' vertex_attr_names(karate)
#' m<-removeVertexTerm(karate, 'color')
#' vertex_attr_names(m)
removeVertexTerm <- function(GG, NAME) {
    if (!is.null(get.vertex.attribute(GG, NAME))) {
        GG <- remove.vertex.attribute(GG, name = NAME)
    }
    if (!is.null(get.vertex.attribute(GG, gsub("_", "", NAME)))) {
        GG <- remove.vertex.attribute(GG, name = gsub("_", "", NAME))
    }
    return(GG)
}
COLLAPSE <- ";"
ESC      <- "|"

#Add GeneNames
#' Annotate Human Gene Names
#'
#' For the protein-protein interaction (PPI) or disease gene interaction (DGN)
#' graphs that have EntrezID as a vertex \code{name} this function extract
#' standard name from \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} and annotate
#' vertices.
#'
#' If vertex \code{name} attrubite stores not EntrezID or network is build
#' not from human genes, other \code{\link[AnnotationDbi]{OrgDb-class}}
#' object could be provided in \code{orgDB} and one of
#' \code{\link[AnnotationDbi]{keytypes}} from that object
#' that correspond to the nature of the vertex \code{name} attrubite could
#' be provided in the \code{keytype} attribute.
#'
#' If for some vertices \code{name} attrubite does not match
#' \code{\link[AnnotationDbi]{keys}} with
#' particular \code{\link[AnnotationDbi]{keytypes}} in the
#' \code{orgDB} object, empty string is added as GeneName.
#'
#' @param gg igraph object to annotate
#' @param orgDB ordDB object, by default human is assumed from
#'         \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#' @param keytype type of IDs stored in the \code{name} vertex attribute,
#'         by default \code{ENTREZID} is assumed.
#'
#' @return igraph object with new vertex attribute \code{GeneName}
#' @export
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi mapIds
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' agg<-annotateGeneNames(gg)
annotateGeneNames <-
    function(gg, orgDB = org.Hs.eg.db, keytype = "ENTREZID") {
        ids <- V(gg)$name
        gn <- suppressMessages(AnnotationDbi::mapIds(orgDB, ids,
                                                     column = "SYMBOL",
                                                     keytype = keytype))
        gg <- removeVertexTerm(gg, "GeneName")
        set.vertex.attribute(gg, "GeneName", V(gg), "")
        V(gg)$GeneName <- gn
        return(gg)
    }
#' Get DiseaseTypes
#'
#' Return vector of disease abbreviations for synaptic PPI analysis.
#'
#' @return vector of disease abbreviations for synaptic PPI analysis.
#'
#' @seealso getDiseases
#' @export
#' @examples
#' getDType()
getDType <- function() {
    #---HDO Disease short names
    dtype  <- vector(length = 12)
    dtype[1]   <- "AD"
    dtype[2]   <- "BD"
    dtype[3]   <- "AUT"
    dtype[4]   <- "SCH"
    dtype[5]   <- "ASD"
    dtype[6]   <- "Epi"
    dtype[7]   <- "ID"
    dtype[8]   <- "HTN"
    dtype[9]   <- "HD"
    dtype[10]  <- "PD"
    dtype[11]  <- "FTD"
    dtype[12]  <- "MS"
    #dtype[12]  <- "DMH";
    #dtype[13]  <- "CNSD";
    return(dtype)
}

#' Get HDO disease IDs
#'
#' Return vector of HDO disease IDs for synaptic PPI analysis.
#'
#' @seealso getDType
#' @return vector of disease IDs of interest
#' @export
#'
#' @examples
#' getDiseases()
getDiseases <- function() {
    #---HDO ID DISEASES of INTEREST
    disn    <- vector(length = 12)
    disn[1]  <- "DOID:10652"#Alzheimer's_disease"
    disn[2]  <- "DOID:3312"#bipolar_disorder"
    disn[3]  <- "DOID:12849"#autistic_disorder"
    disn[4]  <- "DOID:5419"#schizophrenia"
    disn[5]  <- "DOID:0060041"#autism_spectrum_disorder
    disn[6]  <- "DOID:1826"#epilepsy_syndrome
    disn[7]  <- "DOID:1059"
    disn[8]  <- "DOID:10763"
    disn[9]  <- "DOID:12858"
    disn[10] <- "DOID:14330"
    disn[11] <- "DOID:9255"
    disn[12] <- "DOID:2377"
    return(disn)
}

#' Utility function to get vertex ids from vertex attributes
#' The function obtain attribute values and check duplicates in it.
#' It fails if any duplicate found.
#'
#' @param gg graph
#' @param idatt attribute name
#'
#' @return \code{idatt} attribute values
getIDs <- function(gg,idatt){
    ids <- get.vertex.attribute(gg,idatt)
    if(any(duplicated(ids))){
        stop("ID attribute ('",idatt,"') should be unique ",
             "for the graph nodes.\n")
    }
    return(ids)
}
#' Generic annotation function
#'
#' Function to build and fill a vertex attribute given an igraph object. Where
#' parameter 'name' is the new vertex attribute name and values are filled from
#' a two column data.frame supplied to 'value' attribute. The first first
#' containing vertex name IDs, and the second the vertex annotation value.
#'
#'
#' As a first step all attributes with provided names will be removed.
#'
#' @param gg igraph object to annotate
#' @param name name of the attribute
#' @param values annotation data.frame
#' @param idatt optional name of the vertex attribute to map to the
#'        annotation \code{data.frame} first column
#'
#' @return igraph object where vertex attribute \code{name} contains
#' annotation terms separated by semicolon.
#' @export
#'
#' @seealso getAnnotationVertexList
#' @examples
#' g1 <- make_star(10, mode="undirected")
#' V(g1)$name <- letters[1:10]
#' m<-rbind(data.frame(ID=letters[1:10], terms=letters[1:10]),
#' data.frame(ID=letters[1:10], terms=LETTERS[1:10]))
#' g2<-annotateVertex(g1, name='cap', values=m)
#' V(g2)$cap
annotateVertex <- function(gg, name, values, idatt='name') {
    #ids <- V(ggm)$name
    ids <- getIDs(gg,idatt)
    vids <- as.character(values[, 1])
    idx <- match(vids, ids)
    nidx <- which(!is.na(idx))
    vids <- vids[nidx]
    val <- as.character(values[nidx, 2])
    uids <- unique(vids)
    gidx <- match(uids, ids)
    annL <- vapply(uids,
                   function(.x)
                       paste(unique(val[vids == .x]), collapse = COLLAPSE),
                   c(ann = ''))
    ggm <- removeVertexTerm(gg, name)
    ggm <- set.vertex.attribute(graph = ggm,
                                name = name,
                                value = '')
    ggm <- set.vertex.attribute(
        graph = ggm,
        name = name,
        index = gidx,
        value = annL
    )
    return(ggm)
}
#' Escapes elements of list in annotation.
#'
#' In situations when a given list of annotation ID terms may not be well
#' formatted, and therefore not be interoperated as unique. For example, given
#' a list of HDO IDs: HDO:14, HDO:143, HDO:1433, and HDO:14330, a grep for the
#' term HDO:14 could return: HDO:143, HDO:1433, HDO:14330. To avoid this all
#' terms should be enclosed in escape characters, which unlikely to find within
#' annotation itself.
#'
#'
#' NOTE: spaces are treated as regular
#' characters, no trimming is applied before or after escaping.
#'
#'
#' @param annVec vector of annotation strings
#' @param col term list separator character
#' @param esc escape character
#'
#' @return vector of annotation strings with elements escaped
#' @export
#'
#' @seealso unescapeAnnotation
#' @examples
#' annVec<-apply(matrix(letters, ncol=13), 2, paste, collapse=';')
#' cbind(annVec, escapeAnnotation(annVec, ';', '|'))
escapeAnnotation <- function(annVec, col = COLLAPSE, esc = ESC) {
    if (any(grepl(esc, annVec, fixed = TRUE))) {
        stop("Either already escaped or escape charecter found in annotation\n")
    }
    annList <- strsplit(annVec, col, fixed = TRUE)
    escFun <- function(.x) {
        if (length(.x) > 0) {
            return(paste0(esc, .x, esc, collapse = ';'))
        } else{
            return("")
        }
    }
    res <- vapply(annList, escFun, c(ann = ''))
    return(res)
}
#' Unescape annotation strings
#'
#' Function to remove all escape characters from annotation strings
#' (opposite to escapeAnnotation).
#'
#' NOTE: spaces are treated as regular
#' characters, no trimming is applied before or after escaping.
#'
#'
#' @param annVec vector of annotation strings
#' @param col list separator character within annotation string
#' @param esc escape character
#'
#' @return vector of annotation strings with removed escape characters
#' @export
#'
#' @seealso escapeAnnotation
#' @examples
#' annVec<-apply(matrix(letters, ncol=13), 2, paste, collapse=';')
#' escVec<-escapeAnnotation(annVec, ';', '|')
#' cbind(annVec, escVec, unescapeAnnotation(escVec, ';', '|'))
unescapeAnnotation <- function(annVec, col = COLLAPSE, esc = ESC) {
    res <- gsub(esc, '', annVec, fixed = TRUE)
    return(res)
}
#' Return vertex list for each term in annotation attribute
#'
#' For different purposes annotation of graph vertices could be
#' represented in three forms:
#' \describe{
#'   \item{Pairs}{dataframe with vertex ID and annotation terms}
#'   \item{Vertex Annotation}{list named with vertex ID and
#'   containing terms annotating each vertex}
#'   \item{Annotation Vertices}{list named with term and
#'   containing vertex IDs}
#' }
#'
#' This function takes Vertex Annotation from vertex attribute
#' and convert it to Annotation Vertices form.
#'
#' @param g graph to get annotation from
#' @param name annotation attribute name
#' @param col list separation character in attribute, by
#' default is \code{;}
#' @param vid attribute to be used as a vertex ID
#'
#' @return named list with annotation in Annotation Vertices form
#' @export
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' avl<-getAnnotationVertexList(gg, 'TopOntoOVGHDOID')
#' head(avl)
getAnnotationVertexList <-
    function(g, name, vid = 'name', col = COLLAPSE) {
        gda <- prepareGDA(g, name)
        vertices <- get.vertex.attribute(g, vid)
        anNames <- getAnnotationList(gda)
        anL <-
            lapply(anNames, function(.a) {
                vertices[grepl(.a, gda, fixed = TRUE)]
            })
        names(anL) <- unescapeAnnotation(anNames)
        return(anL)
    }

#' Extract unique values from annotations.
#'
#' It is not uncommon to find both duplicated vertex annotation terms, and
#' vertices annotated with multiple terms, in a given annotation list. This
#' function creates a vector of unique annotation terms for each vertex given
#' an input annotation list.
#'
#' @param annVec vector of annotation strings
#' @param col list separator character
#' @param sort how to sort the result list
#' @return vector of unique annotation terms
#' @export
#' @seealso getAnnotationVertexList
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' annVec<-V(gg)$TopOntoOVG
#' al<-getAnnotationList(annVec)
#' al
getAnnotationList <- function(annVec,
                              col = COLLAPSE,
                              sort = c('none', 'string', 'frequency')) {
    sort <- match.arg(sort)
    res <- switch (
        sort,
        none = unique(unlist(strsplit(annVec, col))),
        string = sort(unique(unlist(
            strsplit(annVec, col)
        ))),
        frequency = names(sort(table(
            unlist(strsplit(annVec, col))
        ),
        decreasing = TRUE))
    )
    return(res)
}
#Add topOnto_ovg
#' Annotate graph with disease terms
#'
#' The function loads a human disease annotation matrix called \code{dis}, which
#' contains three columns; the first containing gene Entrez IDs, the second gene
#' Human Disease Ontology (HDO) ID terms, the third gene HDO description terms.
#' For human protein-protein interaction (PPI) or disease-gene networks (DGN)
#' that have human Entrez IDs for the igraph vertex name attribute. The function
#' then performs a many-to-one mapping of each matrix row to a network vertex
#' using matching Entrez IDs, filling the vertices attributes
#' \code{TopOnto_OVG_HDO_ID} and \code{TopOnto_OVG}.
#'
#'
#' @param gg igraph object to annotate
#' @param dis annotation matrix in Pairs form
#' @param idatt optional name of the vertex attributes that contains Entrez IDs
#'
#' @return annotated igraph object
#' @export
#' @seealso getAnnotationVertexList
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.csv", package = "BioNAR")
#' tbl <- read.csv(file, sep="\t")
#' gg <- buildNetwork(tbl)
#' # read HDO data extracted from hxin/topOnto.HDO.db for synaptic network
#' afile<-system.file("extdata", "flatfile_human_gene2HDO.csv",
#' package = "BioNAR")
#' dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' agg<-annotateTopOntoOVG(gg, dis)
annotateTopOntoOVG <- function(gg, dis, idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "TopOnto_OVG")
    gg <- removeVertexTerm(gg, "TopOnto_OVG_HDO_ID")
    #--- Set Disease (geneRIF db) attributes in .gml graph
    set.vertex.attribute(gg, "TopOnto_OVG", V(gg), "")
    set.vertex.attribute(gg, "TopOnto_OVG_HDO_ID", V(gg), "")
    disIDS <- dis[, 3]
    disn <- getDiseases()
    dtype <- getDType()
    for (i in seq_along(ids)) {
        ind1 <- which(disIDS == ids[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            #TDOD: refactor this code to work without disn
            disv <- as.vector(dis[ind1, 1])
            indx <- match(disv, disn)
            for (j in seq_along(disv)) {
                if (!is.na(indx[j])) {
                    if (Str1 == "") {
                        Str1 <- as.character(dtype[indx[j]])
                    }
                    else {
                        Str1 <- paste(c(Str1, as.character(dtype[indx[j]])),
                                      collapse = COLLAPSE)
                    }
                    if (Str2 == "") {
                        Str2 <- as.character(disn[indx[j]])
                    }
                    else {
                        Str2 <- paste(c(Str2, as.character(disn[indx[j]])),
                                      collapse = COLLAPSE)
                    }
                }
            }
        }
        Str1 <-
            paste(unique(strsplit(Str1, COLLAPSE)[[1]]), collapse = COLLAPSE)
        Str2 <-
            paste(unique(strsplit(Str2, COLLAPSE)[[1]]), collapse = COLLAPSE)
        V(gg)[i]$TopOnto_OVG <- as.character(Str1)
        V(gg)[i]$TopOnto_OVG_HDO_ID <- as.character(Str2)
    }
    return(gg)
}

#' Add SCHanno synaptic functional groups
#'
#' The function loads an annotation data matrix of functional groups for
#' schizopherina risk genes (1) called anno, which contains three columns; the
#' first containing gene Entrez IDs, the second gene functional group ID terms,
#' the third gene functional group description terms. The function then performs
#' a many-to-one mapping of each matrix row to a network vertex using matching
#' Entrez IDs, filling the \code{SCHanno} vertices attribute.
#'
#' References:
#' 1. Lips E, Cornelisse L, Toonen R, Min J, Hultman C, the International
#' Schizophernia Consortium, Holmans P, Donovan M, Purcell S, Smit A, Verhage M,
#' Sullivan P, Visscher P, D P: Functional gene group analysis identifies
#' synaptic gene groups as risk factor for schizophrenia.
#' Molecular Psychiatry 2012,17:996–1006.
#'
#' @param gg igraph object to annotate
#' @param anno annotation matrix in Pairs form
#' @param idatt optional name of the vertex attributes that contains Entrez IDs
#'
#' @return annotated igraph object
#' @export
#'
#' @seealso getAnnotationVertexList
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.csv", package = "BioNAR")
#' tbl <- read.csv(file, sep="\t")
#' gg <- buildNetwork(tbl)
#' afile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
#' dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' agg<-annotateSCHanno(gg, dis)
annotateSCHanno <- function(gg, anno,idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "SCHanno")
    #--- Set Family attributes in .gml graph
    set.vertex.attribute(gg, "SCHanno", V(gg), "")
    annoIDS <- as.character(anno[, 3])
    type <-
        unique(unlist(strsplit(as.character(unique(
            anno[, 2]
        )), ",")))
    for (i in seq_along(ids)) {
        ind1 <- which(annoIDS == ids[i])
        Str <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str <- as.character(anno[ind1[1], 2])
            }
            else {
                Str <- paste(as.character(anno[ind1, 2]), collapse = COLLAPSE)
            }
        }
        V(gg)[i]$Schanno <- as.character(Str)
    }
    return(gg)
}
#
#' Add presynaptic functional groups
#'
#' Function takes from \code{anno} matrix manually curated presynaptic genes
#' functional annotation derived from
#' Boyken at al. (2013) <doi:10.1016/j.neuron.2013.02.027>
#' and add them to attributes \code{PRESYNAPTIC}.
#'
#' @param gg graph to update
#' @param anno annotation matrix in Pair form
#' @param idatt optional name of the vertex attributes that contains Entrez IDs
#'
#'
#' @return annotated igraph object
#' @export
#' @seealso getAnnotationVertexList
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' sfile<-system.file("extdata", "PresynAn.csv", package = "BioNAR")
#' pres <- read.csv(sfile,skip=1,header=FALSE,strip.white=TRUE,quote="")
#' gg <- annotatePresynaptic(gg, pres)
annotatePresynaptic <- function(gg, anno,idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "PRESYNAPTIC")
    #--- Set Family attributes in .gml graph
    set.vertex.attribute(gg, "PRESYNAPTIC", V(gg), "")
    annoIDS <- as.character(anno[, 3])
    type <-
        unique(unlist(strsplit(as.character(unique(
            anno[, 2]
        )), ",")))
    for (i in seq_along(ids)) {
        ind1 <- which(annoIDS == ids[i])
        Str <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str <- as.character(anno[ind1[1], 2])
            }
            else {
                Str <- paste(as.character(anno[ind1, 2]), collapse = COLLAPSE)
            }
        }
        V(gg)[i]$PRESYNAPTIC <- as.character(Str)
    }
    return(gg)
}
#
#' Add InterPro Family and Domain annotation to the graph vertices
#'
#' Function takes data from \code{annoF} matrix and add them to attributes
#'  \code{InterPro_Family} for term and \code{InterPro_Family_ID} for IDs.
#'
#' Function takes data from \code{annoD} matrix and add them to attributes
#'  \code{InterPro_Domain} for term and \code{InterPro_Domain_ID} for IDs.
#'
#'
#' @param gg graph to update
#' @param annoF family annotation matrix in Pair form
#' @param annoD domain  annotation matrix in Pair form
#' @param idatt optional name of the vertex attributes that contains Entrez IDs
#'
#' @return annotated igraph object
#' @export
#'
#' @seealso getAnnotationVertexList
annotateInterpro <- function(gg, annoF, annoD,idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "InterPro_Family_ID")
    gg <- removeVertexTerm(gg, "InterPro_Family")
    gg <- removeVertexTerm(gg, "InterPro_Domain_ID")
    gg <- removeVertexTerm(gg, "InterPro_Domain")
    #--- Set InterPro_Family attributes in .gml graph
    set.vertex.attribute(gg, "InterPro_Family_ID", V(gg), "")
    set.vertex.attribute(gg, "InterPro_Family", V(gg), "")
    set.vertex.attribute(gg, "InterPro_Domain_ID", V(gg), "")
    set.vertex.attribute(gg, "InterPro_Domain", V(gg), "")
    annoFIDS <- as.character(annoF[, 3])
    annoDIDS <- as.character(annoD[, 3])
    for (i in seq_along(ids)) {
        ind1 <- which(annoFIDS == ids[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str1 <- as.character(annoF[ind1[1], 2])
            }
            else {
                Str1 <- paste(as.character(annoF[ind1, 2]),
                              collapse = COLLAPSE)
            }
            if (length(ind1) == 1) {
                Str2 <- as.character(annoF[ind1[1], 1])
            }
            else {
                Str2 <- paste(as.character(annoF[ind1, 1]),
                              collapse = COLLAPSE)
            }
        }
        V(gg)[i]$InterPro_Family_ID <- as.character(Str2)
        V(gg)[i]$InterPro_Family   <- as.character(Str1)
        ind1 <- which(annoDIDS == ids[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str1 <- as.character(annoD[ind1[1], 2])
            }
            else {
                Str1 <- paste(as.character(annoD[ind1, 2]),
                              collapse = COLLAPSE)
            }
            if (length(ind1) == 1) {
                Str2 <- as.character(annoD[ind1[1], 1])
            }
            else {
                Str2 <- paste(as.character(annoD[ind1, 1]),
                              collapse = COLLAPSE)
            }
        }
        V(gg)[i]$InterPro_Domain_ID <- as.character(Str2)
        V(gg)[i]$InterPro_Domain   <- as.character(Str1)
    }
    return(gg)
}

#' Annotate nodes with GO terms
#'
#' For the protein-protein interaction (PPI) or disease gene interaction (DGN)
#' graphs that have EntrezID as a vertex \code{name} this function extract
#' GeneOntolgy annotation from \code{orgDB}, which should be
#' \code{\link[AnnotationDbi]{OrgDb-class}}, split them into three ontology
#' group (\code{MF},\code{BP},\code{CC}) and annotate vertices with .
#'
#' If vertex \code{name} attrubite stores not EntrezID or network is build
#' not from human genes, other \code{\link[AnnotationDbi]{OrgDb-class}}
#' object could be provided in \code{orgDB} and one of
#' \code{\link[AnnotationDbi]{keytypes}} from that object
#' that correspond to the nature of the vertex \code{name} attrubite could
#' be provided in the \code{keytype} attribute.
#'
#' If for some vertices \code{name} attrubite does not match
#' \code{\link[AnnotationDbi]{keys}} with
#' particular \code{\link[AnnotationDbi]{keytypes}} in the
#' \code{orgDB} object, empty string is added as GeneName.
#'
#' @param gg igraph object to annotate
#' @param orgDB ordDB object, by default human is assumed from
#'         \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#' @param keytype type of IDs stored in the \code{name} vertex attribute,
#'         by default \code{ENTREZID} is assumed.
#' @param idatt optional name of the vertex attributes that contains IDs
#'        matching the \code{keytype}
#'
#' @return igraph object with new vertex attribute \code{GeneName}
#' @export
#'
#' @import GO.db org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @importFrom dplyr filter
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' ggGO <- annotateGOont(gg)
annotateGOont <- function(gg, orgDB = org.Hs.eg.db, keytype = "ENTREZID",
                          idatt = 'name') {
    if (!inherits(orgDB, 'OrgDb')) {
        stop('orgDB suppose to be subclass of OrgDb.')
    }
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    on <- suppressMessages(AnnotationDbi::select(orgDB,
                                ids,
                                columns = c("GO", 'ONTOLOGY'),
                                keytype = keytype))
    ###### MF annotation ######
    mf <- on %>% dplyr::filter(ONTOLOGY == 'MF') %>%
        dplyr::select(!c(EVIDENCE)) %>% unique
    mfid <- mf[, c(keytype, 'GO')]
    gg <- annotateVertex(gg, "GO_MF_ID", mfid,idatt)
    res <- suppressMessages(AnnotationDbi::select(
        GO.db,
        unique(mf$GO),
        column = c('TERM', 'DEFINITION'),
        keytype = 'GOID'
    ))
    mft <- merge(mf, res, by.x = 'GO', by.y = 'GOID')[, c(keytype, 'TERM')]
    gg <- annotateVertex(gg, "GO_MF", mft, idatt)

    ###### BP annotation ######
    bp <- on %>% dplyr::filter(ONTOLOGY == 'BP') %>%
        dplyr::select(!c(EVIDENCE)) %>% unique
    bpid <- bp[, c(keytype, 'GO')]
    gg <- annotateVertex(gg, "GO_BP_ID", bpid,idatt)
    res <- suppressMessages(AnnotationDbi::select(
        GO.db,
        unique(bp$GO),
        column = c('TERM', 'DEFINITION'),
        keytype = 'GOID'
    ))
    bpt <- merge(bp, res, by.x = 'GO', by.y = 'GOID')[, c(keytype, 'TERM')]
    gg <- annotateVertex(gg, "GO_BP", bpt,idatt)

    ###### CC annotation ######
    cc <- on %>% dplyr::filter(ONTOLOGY == 'CC') %>%
        dplyr::select(!c(EVIDENCE)) %>% unique
    ccid <- cc[, c(keytype, 'GO')]
    gg <- annotateVertex(gg, "GO_CC_ID", ccid,idatt)
    res <- suppressMessages(AnnotationDbi::select(
        GO.db,
        unique(cc$GO),
        column = c('TERM', 'DEFINITION'),
        keytype = 'GOID'
    ))
    cct <- merge(cc, res, by.x = 'GO', by.y = 'GOID')[, c(keytype, 'TERM')]
    gg <- annotateVertex(gg, "GO_CC", cct, idatt)

    return(gg)
}


#' Add GO MF annotation to the graph vertices
#'
#' The function loads an annotation data matrix called \code{annoF}, which
#' contains three columns; the first containing gene Entrez IDs, the second gene
#' GO MF ID terms, the third gene GO MF description terms. The function then
#' performs a many-to-one mapping of each matrix row to a network vertex using
#' matching Entrez IDs, filling the vertices attributes \code{GO_MF_ID} and
#' \code{GO_MF}.
#'
#'
#'
#' @param gg graph to update
#' @param annoF annotation matrix in Pair form
#' @param idatt optional name of the vertex attribute to map to the
#'        annotation \code{data.frame} first column
#'
#' @return annotated igraph object
#' @export
#' @seealso getAnnotationVertexList
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' sfile<-system.file("extdata", "flatfile.go.MF.csv", package = "BioNAR")
#' goMF <- read.table(sfile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' sgg <- annotateGoMF(gg, goMF)
annotateGoMF <- function(gg, annoF,idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "GO_MF")
    gg <- removeVertexTerm(gg, "GO_MF_ID")
    #--- Set Disease (geneRIF db) attributes in .gml graph
    set.vertex.attribute(gg, "GO_MF", V(gg), "")
    set.vertex.attribute(gg, "GO_MF_ID", V(gg), "")
    annoFIDS <- as.character(annoF[, 3])
    typeF <-
        unique(unlist(strsplit(as.character(unique(
            annoF[, 2]
        )), ",")))
    for (i in seq_along(ids)) {
        ind1 <- which(annoFIDS == ids[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str1 <- as.character(annoF[ind1[1], 2])
            }
            else {
                Str1 <- paste(as.character(annoF[ind1, 2]),
                              collapse = COLLAPSE)
            }
            if (length(ind1) == 1) {
                Str2 <- as.character(annoF[ind1[1], 1])
            }
            else {
                Str2 <- paste(as.character(annoF[ind1, 1]),
                              collapse = COLLAPSE)
            }
        }
        V(gg)[i]$GO_MF_ID <- as.character(Str2)
        V(gg)[i]$GO_MF    <- as.character(Str1)
    }
    return(gg)
}
#
#' Add GO BP annotation to the graph vertices
#'
#' The function loads an annotation data matrix called \code{annoF}, which
#' contains three columns; the first containing gene Entrez IDs, the second
#' gene GO BP ID terms, the third gene GO BP description terms. The function
#' then performs a many-to-one mapping of each matrix row to a network vertex
#' using matching Entrez IDs, filling the vertices attributes \code{GO_BP_ID}
#' and \code{GO_BP}.
#'
#'
#'
#' @param gg graph to update
#' @param annoF annotation matrix in Pair form
#' @param idatt optional name of the vertex attribute to map to the
#'        annotation \code{data.frame} first column
#'
#' @return annotated igraph object
#' @export
#' @seealso getAnnotationVertexList
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' sfile<-system.file("extdata", "flatfile.go.BP.csv", package = "BioNAR")
#' goBP <- read.table(sfile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' sgg <- annotateGoBP(gg, goBP)
annotateGoBP <- function(gg, annoF,idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "GO_BP")
    gg <- removeVertexTerm(gg, "GO_BP_ID")
    #--- Set Disease (geneRIF db) attributes in .gml graph
    set.vertex.attribute(gg, "GO_BP", V(gg), "")
    set.vertex.attribute(gg, "GO_BP_ID", V(gg), "")
    annoFIDS <- as.character(annoF[, 3])
    typeF <-
        unique(unlist(strsplit(as.character(unique(
            annoF[, 2]
        )), ",")))
    for (i in seq_along(ids)) {
        ind1 <- which(annoFIDS == ids[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str1 <- as.character(annoF[ind1[1], 2])
            }
            else {
                Str1 <- paste(as.character(annoF[ind1, 2]),
                              collapse = COLLAPSE)
            }
            if (length(ind1) == 1) {
                Str2 <- as.character(annoF[ind1[1], 1])
            }
            else {
                Str2 <- paste(as.character(annoF[ind1, 1]),
                              collapse = COLLAPSE)
            }
        }
        V(gg)[i]$GO_BP_ID <- as.character(Str2)
        V(gg)[i]$GO_BP    <- as.character(Str1)
    }
    return(gg)
}

#
#' Add GO CC  annotation to the graph vertices
#'
#' The function loads an annotation data matrix called \code{annoF}, which
#' contains three columns; the first containing gene Entrez IDs, the second
#' gene GO ID terms, the third gene GO CC description terms. The function
#' then performs a many-to-one mapping of each matrix row to a network vertex
#' using matching Entrez IDs, filling the vertices attributes \code{GO_CC_ID}
#' and  \code{GO_CC}.
#'
#'
#' @param gg graph to update
#' @param annoF annotation matrix in Pair form
#' @param idatt optional name of the vertex attribute to map to the
#'        annotation \code{data.frame} first column
#'
#' @return annotated igraph object
#' @export
#' @seealso getAnnotationVertexList
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read.graph(file, format="gml")
#' sfile<-system.file("extdata", "flatfile.go.CC.csv", package = "BioNAR")
#' goCC <- read.table(sfile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' sgg <- annotateGoCC(gg, goCC)
annotateGoCC <- function(gg, annoF,idatt='name') {
    #ids <- V(gg)$name
    ids <- getIDs(gg,idatt)
    gg <- removeVertexTerm(gg, "GO_CC")
    gg <- removeVertexTerm(gg, "GO_CC_ID")

    set.vertex.attribute(gg, "GO_CC", V(gg), "")
    set.vertex.attribute(gg, "GO_CC_ID", V(gg), "")
    annoFIDS <- as.character(annoF[, 3])
    typeF <-
        unique(unlist(strsplit(as.character(unique(
            annoF[, 2]
        )), ",")))
    for (i in seq_along(ids)) {
        ind1 <- which(annoFIDS == ids[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str1 <- as.character(annoF[ind1[1], 2])
            }
            else {
                Str1 <- paste(as.character(annoF[ind1, 2]),
                              collapse = COLLAPSE)
            }
            if (length(ind1) == 1) {
                Str2 <- as.character(annoF[ind1[1], 1])
            }
            else {
                Str2 <- paste(as.character(annoF[ind1, 1]),
                              collapse = COLLAPSE)
            }
        }
        V(gg)[i]$GO_CC_ID <- as.character(Str2)
        V(gg)[i]$GO_CC    <- as.character(Str1)
    }
    return(gg)
}
