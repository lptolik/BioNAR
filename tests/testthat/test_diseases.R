library(BioNAR)
library(testthat)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read.graph(file, format="gml")
louvainG<-induced_subgraph(gg,V(gg)[V(gg)$louvain%in%c(4,8,10,12)])

test_that('GDA',{
    gda<-prepareGDA(louvainG, 'TopOntoOVGHDOID')
    expect_equal(length(gda),igraph::vcount(louvainG))
    set.seed(100)
    p <- calcDiseasePairs(louvainG,name = "TopOntoOVGHDOID",
                          diseases = c("DOID:10652", "DOID:3312", "DOID:12849"),
                          permute = "n")
    expect_equal(names(p),c("disease_separation",
                            "gene_disease_separation",
                            "disease_localisation"))
    expect_equal(dim(p$disease_localisation),c(3,4))
    expect_equal(dim(p$disease_separation),c(3,3))
    expect_equal(dim(p$gene_disease_separation),c(92,5))
    set.seed(100)
    r <- runPermDisease(louvainG,name = "TopOntoOVGHDOID",
                        diseases = c("DOID:10652", "DOID:3312", "DOID:12849"),
                        Nperm = 10,
                        alpha = c(0.05, 0.01, 0.001))
    expect_equal(names(r),c("Disease_overlap_sig","Disease_location_sig"))
    expect_equal(dim(r$Disease_overlap_sig),c(6,13))
    expect_equal(dim(r$Disease_location_sig),c(3,7))
})
