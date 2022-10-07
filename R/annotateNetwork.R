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

#' Annotate graph from list of files
#'
#' This function is a syntactic sugar wrapper for the
#' \code{\link{annotateVertex}} function. It could be used to quickly load
#' annotation from the set of files, for example all three branches of GO in
#' one run. Each file suppose to be TSV file (use TAB as a column separator)
#' and contains annotation ID in the first column, annotation term in the
#' second and vertex ID in the third. Names of the columns are ignored.
#'
#' @param GG igraph object
#' @param FILES list of file path strings to read annotation from
#' @param NAME vector of names of the vertex property
#' @param IDS vertex IDs
#' @param addIDS if TRUE NAME_ID property will be added
#'
#' @seealso annotateVertex
#'
#' @return igraph object with vertex attributes from NAME contain annotations
loopOverFiles <- function(GG, FILES, NAME, IDS, addIDS) {
    for (f in seq_along(FILES)) {
        GG <- removeVertexTerm(GG, NAME[f])
        if (addIDS) {
            GG <- removeVertexTerm(GG, sprintf("%s_ID", NAME[f]))
        }
        if (file.exists(sprintf("./%s", FILES[f]))) {
            #--- Set Disease (geneRIF db) attributes in .gml graph
            set.vertex.attribute(GG, NAME[f], V(GG), "")
            annoF    <-
                utils::read.table(
                    FILES[f],
                    sep = "\t",
                    skip = 1,
                    strip.white = TRUE,
                    quote = ""
                )
            oo<-getOO(IDS, annoF)
            GG <-
                set.vertex.attribute(GG, NAME[f], V(GG), as.character(oo[, 2]))
            if (addIDS) {
                GG <- set.vertex.attribute(GG,
                                            sprintf("%s_ID", NAME[f]),
                                            V(GG),
                                            as.character(oo[, 3]))
            }
        }
    }
    return(GG)
}

getOO<-function(IDS, annoF){
    annoFIDS <- as.character(annoF[, 3])
    typeF  <-unique(unlist(strsplit(as.character(unique(annoF[, 2])), ",")))
    oo     <- matrix("", ncol = 3, nrow = length(IDS))
    oo[, 1] <- IDS
    for (i in seq_along(IDS)) {
        ind1 <- which(annoFIDS == IDS[i])
        Str1 <- ""
        Str2 <- ""
        if (length(ind1) != 0) {
            if (length(ind1) == 1) {
                Str1 <- as.character(annoF[ind1[1], 2])
            }else {
                Str1 <- paste(as.character(annoF[ind1, 2]),
                              collapse = COLLAPSE)
            }
            if (length(ind1) == 1) {
                Str2 <- as.character(annoF[ind1[1], 1])
            }else {
                Str2 <- paste(as.character(annoF[ind1, 1]), collapse = COLLAPSE)
            }
            if (grepl(COLLAPSE, Str1)) {
                Str1 <- strsplit(Str1, COLLAPSE)[[1]]
                Str1 <- unique(Str1)
                if (length(Str1) > 1) {
                    Str1 <- paste(as.character(Str1), collapse = COLLAPSE)
                }
            }
            if (grepl(COLLAPSE, Str2)) {
                Str2 <- strsplit(Str2, COLLAPSE)[[1]]
                Str2 <- unique(Str2)
                if (length(Str2) > 1) {
                    Str2 <- paste(as.character(Str2), collapse = COLLAPSE)
                }
            }
            oo[i, 2] <- Str1
            oo[i, 3] <- Str2
        }
    }
    return(oo)
}

#Add GeneNames
#' Annotate Human Gene Names
#'
#' For the protein-protein interaction (PPI) or disease gene interaction (DGN)
#' graphs that have EntrezID as a vertex \code{name} this function extract
#' standard name from \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} and annotate
#' vertices.
#'
#' If vertex \code{name} attrubite stores not EntrezID or network is build
#' not from human genes, other \code{\link[AnnotationDbi]{AnnotationDb}}
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
annotateGeneNames <- function(gg,orgDB=org.Hs.eg.db,keytype = "ENTREZID") {
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
#' Generic annotation function
#'
#' It takes name of the attribute, and two column Pair form annotation
#' data.frame with vertex ID in the first column and annotation term in the
#' second. All terms annotating the same vertex ID will be collapsed with
#' semicolon as term separator.
#'
#' As a first step all attributes with provided names will be removed.
#'
#' @param gg igraph object to annotate
#' @param name name of the attribute
#' @param values annotation data.frame
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
annotateVertex <- function(gg, name, values) {
    ggm <- removeVertexTerm(gg, name)
    ggm <- set.vertex.attribute(graph = ggm,
                                name = name,
                                value = '')
    ids <- V(ggm)$name
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
#' In the case when annotation has not carefully planned, some annotation terms
#' could be substring of other, for example \code{\link{grep}} search fo
#' DOID:14 could return DOID:143, DOID:1433, and DOID:14330. To avoid this a
#' ll terms should be enclosed in escape characters, which unlikely to find
#' within annotation itself.
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
#' Perform opposite to escapeAnnotation operations: remove all escape
#' characters from annotation strings.
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
#' It is not uncommon that some nodes are annotated with list of terms and some
#' terms annotates multiple nodes. This function creates vector of unique terms
#' that were used in annotation.
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
#' For the protein-protein interaction (PPI) or disease gene interaction (DGN)
#' graphs that have EntrezID as a vertex name this function takes annotation
#' from \code{dis} matrix and put it on the vertices with attributes
#' \code{TopOnto_OVG} for terms and \code{TopOnto_OVG_HDO_ID} for IDs.
#'
#' @param gg igraph object to annotate
#' @param dis annotation matrix in Pairs form
#'
#' @return annotated igraph object
#' @export
#' @seealso getAnnotationVertexList
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' # read HDO data extracted from hxin/topOnto.HDO.db for synaptic network
#' afile<-system.file("extdata", "flatfile_human_gene2HDO.csv",
#' package = "BioNAR")
#' dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' agg<-annotateTopOntoOVG(gg, dis)
annotateTopOntoOVG <- function(gg, dis) {
    ids <- V(gg)$name
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
#' Function adds Schizopherina related synaptic gene functional annotation
#' from Lips et al., (2012) <doi:10.1038/mp.2011.117> into \code{SCHanno}
#' vertex attribute.
#'
#' @param gg igraph object to annotate
#' @param anno annotation matrix in Pairs form
#'
#' @return annotated igraph object
#' @export
#'
#' @seealso getAnnotationVertexList
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' afile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
#' dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
#' strip.white=TRUE, quote="")
#' agg<-annotateSCHanno(gg, dis)
annotateSCHanno <- function(gg, anno) {
    ids <- V(gg)$name
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
annotatePresynaptic <- function(gg, anno) {
    ids <- V(gg)$name
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
#'
#' @return annotated igraph object
#' @export
#'
#' @seealso getAnnotationVertexList
annotateInterpro <- function(gg, annoF, annoD) {
    ids <- V(gg)$name
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
#'
#'
#' @param gg
#'
#' @return
#' @export
#'
#' @examples
annotateGOall<-function(gg,orgDB=org.Hs.eg.db,keytype = "ENTREZID"){

}


#' Add GO MF annotation to the graph vertices
#'
#' Function takes data from \code{annoF} matrix and add them to attributes
#'  \code{GO_MF} for term and \code{GO_MF_ID} for IDs.
#'
#'
#' @param gg graph to update
#' @param annoF annotation matrix in Pair form
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
annotateGoMF <- function(gg, annoF) {
    ids <- V(gg)$name
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
#' Function takes data from \code{annoF} matrix and add them to attributes
#'  \code{GO_BP} for term and \code{GO_BP_ID} for IDs.
#'
#'
#' @param gg graph to update
#' @param annoF annotation matrix in Pair form
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
annotateGoBP <- function(gg, annoF) {
    ids <- V(gg)$name
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
#' Function takes data from \code{annoF} matrix and add them to attributes
#'  \code{GO_CC} for term and \code{GO_CC_ID} for IDs.
#'
#'
#' @param gg graph to update
#' @param annoF annotation matrix in Pair form
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
annotateGoCC <- function(gg, annoF) {
    ids <- V(gg)$name
    gg <- removeVertexTerm(gg, "GO_CC")
    gg <- removeVertexTerm(gg, "GO_CC_ID")
    #--- Set Disease (geneRIF db) attributes in .gml graph
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
