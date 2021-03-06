# First Script: AD_Comorbidity (Gene Set Enrichment Analysis) 
# by Utpala Nanda Chowdhury 

# 1. Set your work folder
setwd(".../workfolder")

# 2. Load libraries for first script
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

#----------------------------------
# First Example
# A) GSE1297 Alzheimer's Disease (GPL96)

# 3. Download the GEO dataset 
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1297/matrix/"
dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
dataset <- unlist(strsplit(dataset, "\r\n"))
for (ds in dataset) {
  download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
}

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE1297_series_matrix.txt.gz")

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('AD', 12),'CTRL','AD',rep('CTRL', 5),rep('AD', 4),'CTRL',
              rep('AD', 4),'CTRL','CTRL','AD'))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dAD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE1297_table.csv")
nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "Ad study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 10)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE1297.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(AD_GOdata, algorithm = "elim", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultKS.elim),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE1297", useInfo = "all", pdfSW = TRUE)

# 14. Create text file for the correspondence GO terms - genes (this file is mandatory for the script two)
terms <- allRes$GO.ID
genes <- genesInTerm(AD_GOdata,terms)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE1297_correspondence.txt" )
}

#----------------------------------
# Second Example
# A) GSE110226 Alzheimer's Disease (GPL10379)

# 3. Download the GEO dataset 
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110226/matrix/"
dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
dataset <- unlist(strsplit(dataset, "\r\n"))
for (ds in dataset) {
  download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
}

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE110226_series_matrix.txt.gz")
gse <- gse[,sampleNames(gse)!="GSM2982979"]
gse <- gse[,sampleNames(gse)!="GSM2982980"]
gse <- gse[,sampleNames(gse)!="GSM2982981"]
gse <- gse[,sampleNames(gse)!="GSM2982982"]
gse <- gse[,sampleNames(gse)!="GSM2982983"]
gse <- gse[,sampleNames(gse)!="GSM2982984"]
gse <- gse[,sampleNames(gse)!="GSM2982985"]

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('AD', 7),rep('CTRL',6)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dAD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","GeneSymbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$GeneSymbol

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE110226_table.csv")
nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "ad study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 10)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE110226.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(AD_GOdata, algorithm = "elim", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultKS.elim),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE110226", useInfo = "all", pdfSW = TRUE)

# 14. Create text file for the correspondence GO terms - genes (this file is mandatory for the script two)
terms <- allRes$GO.ID
genes <- genesInTerm(AD_GOdata,terms)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE110226_correspondence.txt" )
}

