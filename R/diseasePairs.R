#' @import WGCNA
qscore <- function(zz, FDR){

    LL <- FDR[FDR[, 1] < as.numeric(zz), 2]

    if( length(LL) != 0 ){ return(LL[stats::end(LL)[1]]); }

    return(1)
}

#' Randomly shuffle annotations
#'
#' This function is a convinience wrapper to \code{\link[base]{sample}} with
#' \code{replace= FALSE}
#'
#' @param GNS annotation list to take data from
#' @param N size of the sample
#'
#' @return random list of \code{GNS} values
#' @export
#'
#' @examples
#' permute(LETTERS, 15)
permute <- function(GNS, N){

    temp <- sample(GNS, N, replace=FALSE)

    return(temp)

}

stars    <- c("*", "**", "***")

#' Auxiliary function to replace NAs with zeros.
#'
#' @param x matrix or vector to process
#'
#' @return matrix or vector with \code{NA}s replaced by zero.
#' @export
#'
#' @examples
#' x<-matrix(NA,nrow = 3,ncol = 3)
#' zeroNA(x)
zeroNA<-function(x){
    x[is.na(x)]<-0
    return(x)
}

#' Utility function to convert string matrices without warnings
#'
#' The function handles special version of string matrix where '.' is
#' used instead of NAN.
#'
#' @param x value to convert
#'
#' @return double version of x, or NA
#' @noRd
dotNumeric<-function(x){
    ifelse('.'==x, NA, as.numeric(x))
}

##
# Calculate each diease-pair overlap/seperation on a selected
# synaptic PPI network models, based on analysis described in:
# Menche, J. et al. Uncovering disease-disease relationships
# through the incomplete interactome.
# Science, 347, (6224):1257601 (2015).
##

#source('../setUp.R')

#--- Find all Gene Disease Associations
#GDA <- V(gg)$TopOntoOVG

#Overlap of Disease A and B in the interactome
# GG   => igraph network
# GDA  => gda data for this graph
# disA => name of disease A
# disA => name of disease B
# OO   => minimum shorest paths for each gda, and each disease
diseaseOverlap <- function(GG, GDA, disA, disB, OO){

    #disease A genes
    IDS1  <- V(GG)$name[grepl(disA, GDA, fixed=TRUE)]
    NIDS1 <- length(IDS1)

    #disease B genes
    IDS2  <- V(GG)$name[grepl(disB, GDA, fixed=TRUE)]
    NIDS2 <- length(IDS2)

    #disease A given B
    paths  <- igraph::distances(GG, IDS1, IDS2, weights=NA)
    dsA    <- as.numeric(as.vector(apply(paths, 1, min)))

    #disease B given A
    paths  <- igraph::distances(GG, IDS2, IDS1, weights=NA)
    dsB    <- as.numeric(as.vector(apply(paths, 1, min)))

    #network-based separation between disease A and B
    dAB <- (sum(dsA) + sum(dsB)) / (NIDS1 + NIDS2)

    #network-based localisation of disease A
    indA <- which(colnames(OO) == disA)
    dA   <- mean(as.numeric(as.vector(OO[OO[, indA[1]] != ".", indA[1]])))

    #network-based localisation of disease B
    indB <- which(colnames(OO) == disB)
    dB   <- mean(as.numeric(as.vector(OO[OO[, indB[1]] != ".", indB[1]])))

    #overlap between disease A and B
    sAB <- as.numeric(dAB) - (as.numeric(dA) + as.numeric(dB)) / 2

    return(sAB)

}

#' Prepare mapping for degree-aware annotation shuffling.
#'
#' Function to randomly shuffle vertex annotation terms, whilst preserving the
#' vertex degree originally found with that annotation term.
#'
#' @param gg graph to analyse
#'
#' @param GDA vertex annotations returned by \code{\link{prepareGDA}}
#' @param dtype list of unique annotation terms to analyze
#'
#' @return mapping matrix between vertices, vertex-degree groups and
#'         annotation terms.
#'
#' @export
#' @seealso prepareGDA
#' @seealso getAnnotationList
#' @seealso sampleDegBinnedGDA
#' @examples
#' options("show.error.messages"=TRUE)
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read_graph(file, format="gml")
#' agg<-annotateGeneNames(gg)
#' gda<-prepareGDA(agg, 'TopOntoOVGHDOID')
#' m<-degreeBinnedGDAs(agg, gda, getAnnotationList(gda))
#' c(dim(m), vcount(agg), length(getAnnotationList(gda)))
#' head(m)
degreeBinnedGDAs <- function(gg, GDA, dtype) {
    deg <- degree(gg)
    bins <- table(deg)
    map <- cbind(names(deg), as.vector(deg))
    map <- cbind(map, match(map[, 2], names(bins)))

    #nGDAs <- length(dtype)

    for (i in seq_along(dtype)) {
        map <- cbind(map, ifelse(grepl(dtype[i], GDA, fixed = TRUE), 1, 0))
    }

    colnames(map) <- c("EntrezID", "Degree", "Bin", dtype)

    return(map)

}

#' Function to randomly shuffle vertex annotation terms, whilst preserving the
#' vertex degree originally found with that annotation term..
#'
#' @param org.map degree-annotation mapping returned by
#'        \code{\link{degreeBinnedGDAs}}
#' @param term annotation term to shuffle
#'
#' @return vertex IDs to assign \code{term} in shuffled annotation
#'
#' @export
#' @seealso degreeBinnedGDAs
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read_graph(file, format="gml")
#' agg<-annotateGeneNames(gg)
#' gda<-prepareGDA(agg, 'TopOntoOVGHDOID')
#' diseases<-getAnnotationList(gda)
#' m<-degreeBinnedGDAs(agg, gda, diseases)
#' sampleDegBinnedGDA(m, diseases[1])
sampleDegBinnedGDA <- function(org.map, term) {
    gda.indx <- match(term, colnames(org.map))

    rnd.gene.set <- NULL

    if (length(gda.indx) > 0) {
        gda.set <- as.vector(org.map[org.map[, gda.indx] == 1, 3])
        Nset <- length(gda.set)

        rnd.gene.set <- rep(NA, Nset)
        map <- org.map
        for (i in seq_len(Nset)) {
            seq.map <- seq(1, dim(map)[1], 1)
            rnd.indx <- seq.map[!is.na(match(as.numeric(map[, 3]),
                                             as.numeric(gda.set[i])))]
            if (length(rnd.indx) > 1) {
                rnd.indx <- as.numeric(sample(rnd.indx))[1]
            }
            if (length(rnd.indx) > 0) {
                rnd.gene.set[i] <- map[rnd.indx, 1]
                map <- map[-rnd.indx,]
            }
        }
    }

    return(rnd.gene.set)

}


#' Function to return vertex annotation from a graph in the Vertex Annotation
#' form and format it for further analysis.
#'
#' @param gg igraph object to take annotation from
#' @param name name of the vertex attribute that contains annotation. If graph
#' has no such vertex attribute an error is thrown..
#'
#' @return escaped annotation in Vertex Annotation form
#' @export
#'
#' @seealso getAnnotationVertexList
#' @seealso escapeAnnotation
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read_graph(file, format="gml")
#' agg<-annotateGeneNames(gg)
#' gda<-prepareGDA(agg, 'TopOntoOVGHDOID')
#' gda<-prepareGDA(agg, 'TopOntoOVGHDOID')
#' head(gda)
prepareGDA <- function(gg, name) {
    if (!name %in% vertex_attr_names(gg)) {
        stop("There is no attribute '", name, "' in the graph.\n")
    }
    gda <- vertex_attr(gg, name)
    gda <- escapeAnnotation(gda)
    return(gda)
}

#' Calculate each disease-disease pair overlap given a list of disease terms.
#'
#' Calculate each disease-disease pair overlap (or separation) on a given PPI
#' network model, based on analysis described in Menche et al. 2015
#'
#' @references Menche, J. et al. Uncovering disease-disease relationships
#'             through the incomplete interactome.
#'             Science, 347, (6224):1257601 (2015).
#'
#'
#' @param gg interactome network as igraph object
#' @param name name of the attribute that stores disease annotation
#' @param diseases list of diseases to match
#' @param permute type of permutations. \code{none} -- no permutation is
#'        applied, \code{random} -- annotation is randomly shuffled,
#'        \code{binned} -- annotation is shuffled in a way to preserve node
#'        degree-annotation relationship by \code{\link{degreeBinnedGDAs}}.
#'
#' @return list with three matrices:
#' * disease_separation -- Ndisease X Ndisease matrix of separations
#' * gene_disease_separation -- Ngenes X Ndisease+2 matrix of gene-disease
#' separation
#' * disease_localisation -- matrix with diseases in rows and number of genes
#' (N), average and standard deviation of gene-disease separation in columns
#' @export
#'
#' @seealso degreeBinnedGDAs
#' @seealso sampleDegBinnedGDA
#' @md
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read_graph(file, format="gml")
#' agg<-annotateGeneNames(gg)
#' p <- calcDiseasePairs(
#' agg,
#' name = "TopOntoOVGHDOID",
#' diseases = c("DOID:10652", "DOID:3312", "DOID:12849"),
#' permute = "n"
#' )
#' p$disease_separation
calcDiseasePairs <- function(gg,
                             name,
                             diseases = NULL,
                             permute = c('none', 'random', 'binned')) {
    permute <- match.arg(permute)
    gda <- prepareGDA(gg, name)
    NN  <- length(which(gda != ""))
    if (is.null(diseases)) {
        diseases <- getAnnotationList(gda, sort = 'freq')
    } else{
        remove <- c()
        diseases <- escapeAnnotation(diseases)
        for (d in seq_along(diseases)) {
            if (!any(grepl(diseases[d], gda))) {
                remove <- c(remove, d)
            }
        }
        if (length(remove) > 0) {
            diseases <- diseases[-remove]
        }
    }
    nDiseases <- length(diseases)
    if (permute != 'none') {
        rgda <- matrix(NA, nrow = vcount(gg), ncol = (nDiseases))
        colnames(rgda) <- c(diseases)
        if (permute == 'binned') {
            map <- degreeBinnedGDAs(gg, gda, diseases)
        }
        for (d in seq_along(diseases)) {
            IDS <- V(gg)$name[grepl(diseases[d], gda, fixed = TRUE)]
            N   <- length(IDS)
            if (permute == 'random') {
                IDS <- BioNAR::permute(V(gg)$name, N) #case
            } else if (permute == 'binned') {
                IDS <- sampleDegBinnedGDA(map, diseases[d])
            }
            rgda[match(IDS, V(gg)$name), d] <- 1
        }
        gda <- apply(rgda, 1, function(.x)
            paste(diseases[!is.na(.x)],
                  collapse = COLLAPSE))
    }
    res           <- matrix(0 , ncol = 4, nrow = nDiseases)
    colnames(res) <- c("Disease", "N", "mean_ds", "SD_ds")
    res[, 1]       <- unescapeAnnotation(diseases)
    oo <- matrix(".", nrow = length(gda), ncol = (nDiseases + 2))
    colnames(oo) <- c("Gene.ID", "Gene.Name", diseases)
    oo[, 1]       <- V(gg)$name#[gda !=""]
    oo[, 2]       <- V(gg)$GeneName#[gda !=""]
    for (d in seq_along(diseases)) {
        IDS <- V(gg)$name[grepl(diseases[d], gda, fixed = TRUE)]
        N   <- length(IDS)
        XX <- igraph::distances(gg, IDS, IDS, weights = NA)
        diag(XX) <- NA
        ds <- apply(XX, 1, min, na.rm = TRUE)
        indX <- match(names(ds), oo[, 1])
        oo[indX, (2 + d)] <- as.vector(ds)
        res[d, 2] <- as.character(N)
        res[d, 3] <- as.character(mean(ds))
        res[d, 4] <- as.character(stats::sd(ds))
    }
    DAB <- matrix(NA, ncol = nDiseases, nrow = nDiseases)
    colnames(DAB) <- diseases
    rownames(DAB) <- diseases
    #--- NOTE ---#
    # DAB is bound by -dmax <= DAB <= dmax
    # where dmax denotes the diameter of the network
    # dmax <- diameter(gg, directed=F)
    #------------#
    ##--- calculate disease-disease overlap
    for (i in seq_along(diseases)) {
        for (j in seq.int(from = i, to = nDiseases)) {
            DAB[i, j] <- 0
            if (i != j) {
                DAB[i, j] <- diseaseOverlap(gg, gda, rownames(DAB)[i],
                                            colnames(DAB)[j], oo)
            }
        }
    }
    colnames(DAB) <- unescapeAnnotation(diseases)
    rownames(DAB) <- unescapeAnnotation(diseases)
    colnames(oo) <-
        c("Gene.ID", "Gene.Name", unescapeAnnotation(diseases))
    if (permute == 'none') {
        oo <- oo[gda != "",]
    }
    oo<-as.data.frame(oo)
    iidx<- c(1,2)
    oo[,-iidx]<-suppressWarnings(lapply(oo[,-iidx],as.numeric))
    res<-as.data.frame(res)
    res[,-1]<-suppressWarnings(lapply(res[,-1],as.numeric))
    return(
        list(
            disease_separation = as.data.frame(DAB),
            gene_disease_separation = oo,
            disease_localisation = res
        )
    )
}

mean0 <- function(x) {
    return(mean(zeroNA(x)))
}
sd0 <- function(x) {
    return(stats::sd(zeroNA(x)))
}
toNum <- function(.x) {
    Nc <- dim(.x$gene_disease_separation)[2]
    apply(.x$gene_disease_separation[, seq(3, Nc)], c(1, 2), dotNumeric)
}

#' Calculate disease-disease pair overlaps on permuted network to estimate
#' its statistical significance
#'
#' Function to calculate the disease-pair overlap characteristics of an inputted
#' network, before applying \code{Nperm} permutations on the disease annotations
#' of #' type "random" or "binned" \code{permute}. From the permuted networks
#' the function estimates the significance of disease overlap: p-value,
#' Bonferoni-adjusted p-value, and q-value in the \code{Disease_overlap_sig}.
#' The function also compares the average disease separation between inputted
#' and permuted networks, and calculates its significance using the Wilcox test
#' and store. Significance of disease-pair overlap and disease separation
#' results are stored in the matrix \code{Disease_location_sig}.
#'
#'
#'
#' Run with care, as large number of permutations could require a lot of
#' memory and be timeconsuming.
#'
#' @param gg interactome network as igraph object
#' @param name name of the attribute that stores disease annotation
#' @param diseases list of diseases to match
#' @param Nperm number of permutations to apply
#' @param permute type of permutations.
#'        \code{random} -- annotation is randomly shuffled,
#'        \code{binned} -- annotation is shuffled in a way to preserve node
#'        degree-annotation relationship by \code{\link{degreeBinnedGDAs}}.
#' @param alpha statistical significance levels
#'
#' @return list of two matrices: \code{Disease_overlap_sig} gives s
#'         tatistics for each pair of disease, and
#'         \code{Disease_location_sig} gives intra-disease statistics
#' @export
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
#' gg <- igraph::read_graph(file, format="gml")
#' agg<-annotateGeneNames(gg)
#' r <- runPermDisease(
#' agg,
#' name = "TopOntoOVGHDOID",
#' diseases = c("DOID:10652", "DOID:3312", "DOID:12849", "DOID:1826"),
#' Nperm = 10,
#' alpha = c(0.05, 0.01, 0.001))
#' r$Disease_location_sig
runPermDisease <- function(gg,
                           name,
                           diseases = NULL,
                           Nperm = 100,
                           permute = c('random', 'binned'),
                           alpha = c(0.05, 0.01, 0.001)) {
    permute <- match.arg(permute)
    resD <-
        calcDiseasePairs(
            gg = gg,
            name = name,
            diseases = diseases,
            permute = 'none'
        )
    ds <- resD$gene_disease_separation
    loc <- resD$disease_localisation
    resL <-
        lapply(seq_len(Nperm), function(.x)
            calcDiseasePairs(
                gg = gg,
                name = name,
                diseases = diseases,
                permute = permute
            ))
    resGDS <- vapply(resL, toNum, FUN.VALUE = toNum(resL[[1]]))
    m <- apply(resGDS, c(1, 2), mean0)
    RANds <-
        cbind(as.data.frame(resL[[1]]$gene_disease_separation[, seq_len(2)]),
              as.data.frame(m))
    disn <- colnames(ds)[3:length(ds[1, ])]
    disease_location_sig           <-
        matrix(0 , ncol = 7, nrow = length(disn))
    colnames(disease_location_sig) <-
        c("HDO.ID",
          "N",
          "mean_ds",
          "SD_ds",
          "Ran_mean_ds",
          "Ran_SD_ds",
          "Utest.pvalue")
    disease_location_sig[, 1] <- disn
    disease_location_sig[, 2] <-
        loc[match(disease_location_sig[, 1], loc[, 1]), 2]
    for (i in seq_along(disn)) {
        indx <- !is.na(ds[, (2 + i)])# != "."
        ids <- ds[indx, 1]
        DS       <- as.numeric(as.vector(ds[indx, (2 + i)]))
        disease_location_sig[i, 3] <- as.numeric(mean(DS))
        disease_location_sig[i, 4] <- as.numeric(stats::sd(DS))
        indy <- match(ids, RANds[, 1])
        RDS <- as.numeric(as.vector(RANds[indy, (2 + i)]))
        disease_location_sig[i, 5] <- as.numeric(mean0(RDS))
        disease_location_sig[i, 6] <- as.numeric(sd0(RDS))
        disease_location_sig[i, 7] <- 1.0
        if (!any(is.infinite(DS)) && !any(is.nan(DS)) &&
            !any(is.na(DS)) &&  !any(is.infinite(RDS)) &&
            !any(is.nan(RDS)) && !any(is.na(RDS))) {
            if (length(DS) != 0 && length(RDS) != 0) {
                wt       <- stats::wilcox.test(DS, RDS)
                disease_location_sig[i, 7] <- as.numeric(wt$p.value)
            }
        }
    }
    disease_location_sig<-as.data.frame(disease_location_sig)
    disease_location_sig[,-1]<-suppressWarnings(lapply(
        disease_location_sig[,-1],as.numeric))
    sAB <- as.matrix(resD$disease_separation)
    RAW_sAB <- vapply(resL, function(.x)
        as.matrix(.x$disease_separation),
        FUN.VALUE = sAB)
    RAN_sAB_mean <- apply(RAW_sAB, c(1, 2), mean0)
    RAN_sAB_sd <- apply(RAW_sAB, c(1, 2), sd0)
    perms <- dim(RAW_sAB)[3]
    Nn    <- length(disn)
    NELE  <- Nn * (Nn + 1) / 2
    Nlevels <- NELE

    CN <-  c(
        "HDO.ID.A",
        "N.A",
        "HDO.ID.B",
        "N.B",
        "sAB",
        "Separated",
        "Overlap",
        "zScore",
        "pvalue",
        "Separation/Overlap.than.chance",
        "Bonferroni",
        "p.adjusted",
        "q-value"
    )
    zs <- matrix(".", nrow = NELE, ncol = length(CN))
    colnames(zs) <- CN
    tests <- matrix(0, nrow = NELE, ncol = perms)
    for (k in 0:(NELE - 1)) {
        i <- floor((2 * Nn + 1 - sqrt((2 * Nn + 1) * (2 * Nn + 1) - 8 * k)) / 2)

        j <- k - Nn * i + i * (i - 1) / 2

        i <- i + 1

        j <- j + i

        zScore <- 0
        if (!is.nan(as.numeric(RAN_sAB_sd[i, j]))) {
            if (as.numeric(RAN_sAB_sd[i, j]) != 0) {
                zScore <- (as.numeric(as.vector(sAB[i, j])) -
                               as.numeric(as.vector(RAN_sAB_mean[i, j]))) /
                    (as.numeric(as.vector(RAN_sAB_sd[i, j])))
            }
            pval <- stats::pnorm(-abs(zScore))
            pval <- 2 * pval
            zs[(k + 1), 1] <- disn[i]
            zs[(k + 1), 2] <-
                as.character(loc[which(loc[, 1] == disn[i]), 2])
            zs[(k + 1), 3] <- disn[j]
            zs[(k + 1), 4] <-
                as.character(loc[which(loc[, 1] == disn[j]), 2])
            zs[(k + 1), 5] <- as.character(sAB[i, j])
            zs[(k + 1), 6] <-
                ifelse((as.numeric(zs[(k + 1), 5]) > 0), "YES", ".")
            zs[(k + 1), 7] <-
                ifelse((as.numeric(zs[(k + 1), 5]) < 0), "YES", ".")
            zs[(k + 1), 8] <- as.character(zScore)
            zs[(k + 1), 9] <- as.character(pval)
            zs[(k + 1), 10] <- ifelse((as.numeric(zs[(k + 1), 8]) < 0),
                                      "Smaller", "larger")
            temp <- "."
            for (x in seq_along(alpha)) {
                if (as.numeric(zs[(k + 1), 9]) <
                    as.numeric(alpha[x] / Nlevels)) {
                    temp <- stars[x]
                }
            }
            zs[(k + 1), 11] <- temp
        } else {
            zs[(k + 1), 1] <- disn[i]
            zs[(k + 1), 2] <-
                as.character(loc[which(loc[, 1] == disn[i]), 2])
            zs[(k + 1), 3] <- disn[j]
            zs[(k + 1), 4] <-
                as.character(loc[which(loc[, 1] == disn[j]), 2])
        }
    }
    zs[, 12] <- stats::p.adjust(as.numeric(zs[, 9]), method = "BY")
    zs[, 13] <- WGCNA::qvalue(as.numeric(zs[, 9]))$qvalue
    zs<-as.data.frame(zs)
    zs[, c("N.A", "N.B", "sAB",
           "zScore", "pvalue",
           "p.adjusted", "q-value")] <- suppressWarnings(
               lapply(zs[,c("N.A","N.B","sAB","zScore",
                            "pvalue","p.adjusted","q-value")],
                      as.numeric))
    return(list(
        Disease_overlap_sig = zs,
        Disease_location_sig = disease_location_sig
    ))
}
