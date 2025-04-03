RUV_total <- function(raw,pData,fData,k,hkgenes = NULL,exclude = NULL){
  
  ### INPUT: raw - p x n raw expressions with p genes and n samples
  ###        pData - phenotype metadata across samples
  ###        fData - feature metadata across genes
  ###        k - number of dimensions of unwanted variation estimated
  ###        exclude - vector of gene names to exclude
  
  # library(RUVSeq)
  # library(DESeq2)
  # library(limma)
  # library(matrixStats)
	
#   if (!is.null(hkgenes)){
# 	  
#   fData(set)$Class[rownames(set) %in% hkgenes] = 'Housekeeping'
# 	  
# 	  }
  
  fData = fData[rownames(raw),]
  int   = intersect(rownames(raw),rownames(fData))
  fData = fData[int,]
  raw   = raw[int,]

  set <- newSeqExpressionSet(as.matrix(round(raw)),
                             phenoData=pData,
                             featureData=fData)

  cIdx       <- rownames(set)[fData(set)[,1]=="Housekeeping"]
  #browser()
  #cIdx <- rownames(set)[fData(set)$Class == "Housekeeping"]
  cIdx       <- cIdx[!(cIdx %in% exclude)]
  #x   <- as.factor(pData$Group)
  set        <- betweenLaneNormalization(set, which="upper")
  set        <- RUVg(set, cIdx, k=k)
  dds        <- DESeqDataSetFromMatrix(counts(set),colData=pData(set),design=~1)
  rowData(dds) <- fData
  dds        <- estimateSizeFactors(dds)
  dds        <- estimateDispersionsGeneEst(dds)
  cts        <- counts(dds, normalized=TRUE)
  disp       <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)[,"dispGeneEst"] <- disp
  #mcols(dds)$dispGeneEst <- disp
  dds        <- estimateDispersionsFit(dds, fitType="mean")
  vsd        <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat        <- assay(vsd)
  covars     <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat        <- removeBatchEffect(mat, covariates=covars)
  assay(vsd) <- mat
  return(list(set = set, vsd = vsd))


}

imagingQC <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for imaging quality

	fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
	if (!(fovRatio > .75)) {return('Flag')}
	if (fovRatio > .75) {return('No flag')}

}

bindingDensityQC <- function(rcc,low,high){
  
  
  #### INPUT: rcc - input from rcc
  ####         low, high - the lower and upper limits for binding density
  #### OUTPUT: flag for binding density

	bd = as.numeric(rcc$Lane_Attributes[6])
	if(!(bd < high & bd > low)) {return('Flag')}
	if (bd < high & bd > low) {return('No flag')}
	

}

limitOfDetectionQC <- function(rcc,numSD = 0){

  #### INPUT: rcc - input from rcc
  ####         numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: flag for limit of detection
  
	counts = rcc$Code_Summary
	posE = as.numeric(counts$Count[counts$CodeClass == 'Positive'])
	negControls = as.numeric(counts$Count[grepl('Negative',counts$CodeClass)])
	if(  any(!(posE > mean(negControls) + numSD*sd(negControls))) ) {return('Flag')}
	if ( all(posE > mean(negControls) + numSD*sd(negControls)) ) {return('No flag')}

}

positiveLinQC <- function(rcc){

  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for linearity for positive controls
  
  
	counts = rcc$Code_Summary
	posControls = as.numeric(counts$Count[grepl('Positive',counts$CodeClass)])
	known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
	r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
	if(!(r2 > .95) | is.na(r2)) {return('Flag')}
	if(r2 > .95) {return('No flag')}

}

makeRLEplot <- function(data,metadata, title="", int.col="black",
                        no_leg=FALSE, no_title=FALSE,  title_size=24,
                        no_xlab=FALSE, no_ylab=FALSE,
                        add_xlabs=FALSE, x_size=0.5){
  
  #### INPUT: data - matrix of expressions with genes on rows and samples on columns
  ####        metadata - matrix of metadata with a column that corresponds to the colnames of data
  ####        id - colname of sample ids
  #### OUTPUT: ggplot2 RLE plot
  
  data = data - apply(data,1,median)
  # stack = stack(data)
  # colnames(stack)[1] = id
  # stackPlot = merge(stack,metadata,by=id)
  # colnames(stackPlot)[1:2] = c('Sample','values')

  rle_plots = reshape2::melt(data) %>% 
    mutate(group = metadata$group[match(Var2, metadata$id)]) %>%
    ggplot(., aes(x=Var2, y=value, fill=group))+
  #rle_plots = ggplot(data = stackPlot,aes(x = Sample,y = values, color = Group)) +
    geom_boxplot(coef = 6) +
    geom_hline(yintercept=0, linetype="longdash", color=int.col)+
    theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.title = element_text(size = title_size),
          legend.title=element_text(size=20),
          legend.text=element_text(size=20),
          strip.text = element_text(size=24),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey", 
                                      fill = NA, size = .1),
          legend.position = 'bottom',
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    labs(title=title)+
    xlab('Sample')+
    ylab('Median deviation of log expression')+ 
    {if(add_xlabs){
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face="bold", size=rel(x_size)))
    }}+
    {if(no_leg){
      theme(legend.position="none")
      #theme(legend.title=element_blank(), 
      #      legend.text=element_blank())
    }}+
    {if(no_title){
      theme(plot.title = element_blank())
    }}+
    {if(no_ylab){
      theme(#axis.text.y  = element_blank(),
            axis.title.y = element_blank())
    }}+
    {if(no_xlab){
      theme(axis.title.x = element_blank())
    }}+
    ylim(c(-4,4))
  
  return(rle_plots)
  
}
