library(BioNAR)
library(testthat)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "BioNAR")
gg <- igraph::read_graph(file, format="gml")
louvain4<-induced_subgraph(gg,V(gg)[V(gg)$louvain==4])

test_that("GO annotation is added",{
    ggGO <- annotateGOont(gg,orgDB = org.Hs.eg.db, keytype = "ENTREZID")
    expect_equal(length(grep('GO',vertex_attr_names(ggGO))),6)
})

test_that("InterPro annotation is added",{
    ggGO <- annotateGOont(gg,orgDB = org.Hs.eg.db, keytype = "ENTREZID")
    expect_equal(length(grep('GO',vertex_attr_names(ggGO))),6)
})


test_that('Term annotation cleanef from the graph',{
    expect_true('louvain' %in% vertex_attr_names(gg))
    ggT<-removeVertexTerm(gg,'louvain')
    expect_false('louvain' %in% vertex_attr_names(ggT))
    expect_true('TopOntoOVG' %in% vertex_attr_names(gg))
    ggT<-removeVertexTerm(gg,'Top_Onto_OVG')
    expect_false('TopOntoOVG' %in% vertex_attr_names(ggT))
})

test_that('GeneNames annotated',{
    gn<-V(gg)$GeneName
    idx<-grep('DLG',gn)
    gNoN<-removeVertexTerm(gg,'GeneName')
    expect_false('GeneName' %in% vertex_attr_names(gNoN))
    gN<-annotateGeneNames(gNoN)
    expect_true('GeneName' %in% vertex_attr_names(gN))
    expect_equal(V(gN)$GeneName[idx],gn[idx])
})

test_that('DiseaseType',{
    dt<-getDType()
    expect_equal(length(dt),12)
    dd<-getDiseases()
    expect_equal(length(dd),12)
})

test_that('Vertex annotation',{
    adf<-rbind(data.frame(ID=V(gg)$name,val=V(gg)$TopOntoOVG),
               data.frame(ID=V(gg)$name,val=V(gg)$TopOntoOVGHDOID))
    idx<-grep(V(gg)$TopOntoOVG[1],V(gg)$TopOntoOVG)
    ggA<-annotateVertex(gg,'annot',adf)
    expect_true('annot' %in% vertex_attr_names(ggA))
    expect_identical(grep(V(gg)$TopOntoOVG[1],V(ggA)$annot),idx)
    avl<-getAnnotationVertexList(gg, 'TopOntoOVGHDOID')
    al<-getAnnotationList(V(gg)$TopOntoOVGHDOID)
    expect_equal(names(avl),al)
    al1<-getAnnotationList(V(gg)$TopOntoOVGHDOID,sort = 'frequency')
    expect_equal(names(avl)[order(sapply(avl,length),decreasing = TRUE)],al1)
    al2<-getAnnotationList(V(gg)$TopOntoOVGHDOID,sort = 'string')
    expect_equal(sort(names(avl)),al2)
})

test_that('Escape annotation',{
    annVec<-apply(matrix(letters, ncol=13), 2, paste, collapse=';')
    escVec<-escapeAnnotation(annVec, ';', '|')
    expect_equal(length(escVec),length(annVec))
    expect_error(escapeAnnotation(escVec, ';', '|'),'already escaped')
    uescVec<-unescapeAnnotation(escVec, ';', '|')
    expect_equal(annVec,uescVec)
})

test_that('File annotation',{
    # read HDO data extracted from hxin/topOnto.HDO.db for synaptic network
    afile<-system.file("extdata", "flatfile_human_gene2HDO.csv",
    package = "BioNAR")
    dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
    strip.white=TRUE, quote="")
    agg<-annotateTopOntoOVG(louvain4, dis)
    expect_true('TopOnto_OVG' %in% vertex_attr_names(agg))
    expect_true('TopOnto_OVG_HDO_ID' %in% vertex_attr_names(agg))
    expect_equal(getAnnotationList(V(agg)$TopOnto_OVG_HDO_ID,sort = 'string'),
                 getAnnotationList(V(louvain4)$TopOntoOVGHDOID,sort = 'string'))
    afile<-system.file("extdata", "SCH_flatfile.csv", package = "BioNAR")
    dis    <- read.table(afile, sep="\t", skip=1, header=FALSE,
    strip.white=TRUE, quote="")
    agg<-annotateSCHanno(louvain4, dis)
    expect_true('Schanno' %in% vertex_attr_names(agg))
    sfile<-system.file("extdata", "PresynAn.csv", package = "BioNAR")
    pres <- read.csv(sfile,skip=1,header=FALSE,strip.white=TRUE,quote="")
    agg <- annotatePresynaptic(louvain4, pres)
    expect_true('PRESYNAPTIC' %in% vertex_attr_names(agg))
    sfile<-system.file("extdata", "flatfile.go.MF.csv", package = "BioNAR")
    goMF <- read.table(sfile, sep="\t", skip=1, header=FALSE,
    strip.white=TRUE, quote="")
    agg <- annotateGoMF(louvain4, goMF)
    expect_true('GO_MF' %in% vertex_attr_names(agg))
    expect_true('GO_MF_ID' %in% vertex_attr_names(agg))
    sfile<-system.file("extdata", "flatfile.go.BP.csv", package = "BioNAR")
    goBP <- read.table(sfile, sep="\t", skip=1, header=FALSE,
    strip.white=TRUE, quote="")
    agg <- annotateGoBP(louvain4, goBP)
    expect_true('GO_BP' %in% vertex_attr_names(agg))
    expect_true('GO_BP_ID' %in% vertex_attr_names(agg))
    sfile<-system.file("extdata", "flatfile.go.CC.csv", package = "BioNAR")
    goCC <- read.table(sfile, sep="\t", skip=1, header=FALSE,
    strip.white=TRUE, quote="")
    agg <- annotateGoCC(louvain4, goCC)
    expect_true('GO_CC' %in% vertex_attr_names(agg))
    expect_true('GO_CC_ID' %in% vertex_attr_names(agg))
})

test_that('Enrichment',{
    anL<-getAnnotationVertexList(gg, 'TopOntoOVGHDOID')
    res<-clusterORA(gg, alg='louvain', name='TopOntoOVGHDOID', vid='name',
                    alpha=0.1)
    expect_equal(dim(res),c(18,17))
    expect_true(all('louvain'==res$alg))
    expect_equal(unique(res$cl),c(1,3,5,15))
    expect_true(all(198==res$Cn[res$cl==1]))
})
