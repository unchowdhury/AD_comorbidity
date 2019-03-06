# SCRIPT ONE: COELIAC DISEASE EXAMPLES (GENE SET ENRICHMENT ANALYSIS) 
# by Eugenio Del Prete

# 1. Set your work folder
setwd("D:/Bioinformatics/With Ali Moni Aus/AD_Disease_Comorbidity")

# 2. Load libraries for script one
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

#----------------------------------
# A) GSE110226 AD

# 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110226/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

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
# nrow(subtable_result[subtable_result$P.Value<0.05,])

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
                 nodeSize = 5)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE110226.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
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
#----------------------------------

#B) GSE1297

# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1297/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

gse <- getGEO(filename = "GSE1297_series_matrix.txt.gz")
#gse <- gse[,sampleNames(gse)!="GSM2335908"]

d <- factor(c(rep('AD', 12),'CTRL','AD',rep('CTRL', 5),rep('AD', 4),'CTRL',
              rep('AD', 4),'CTRL','CTRL','AD'))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dAD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

write.csv(subtable_result,"GSE1297_table.csv")
 nrow(subtable_result[subtable_result$P.Value<0.05,])

topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

AD_GOdata <- new("topGOdata",
                 description = "Ad study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 5)

n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE1297.csv")

resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE1297", useInfo = "all", pdfSW = TRUE)

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

#------------------------------------------------------------------------

# C) GSE33000

# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33000/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE33000_series_matrix.txt.gz")
for (i in 1:157)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1424",246+i, sep="")]  
}

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('AD', 310),rep('CTRL',157)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dAD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","ORF","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$ORF

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE33000_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>0.5)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
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
write.csv(ug,"GO_GSE33000.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 26)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE33000", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE33000_correspondence.txt" )
}




##############################################################################

#----------------------------------
# D) GSE48350a (hippocampus)


# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE48350_series_matrix.txt.gz")

gse <- gse[,sampleNames(gse)!="GSM300166"]
gse <- gse[,sampleNames(gse)!="GSM300167"]

gse <- gse[,sampleNames(gse)!="GSM300173"]
gse <- gse[,sampleNames(gse)!="GSM300183"]
gse <- gse[,sampleNames(gse)!="GSM300186"]
gse <- gse[,sampleNames(gse)!="GSM300188"]
gse <- gse[,sampleNames(gse)!="GSM300189"]
gse <- gse[,sampleNames(gse)!="GSM300191"]
gse <- gse[,sampleNames(gse)!="GSM300192"]
gse <- gse[,sampleNames(gse)!="GSM300256"]
gse <- gse[,sampleNames(gse)!="GSM300258"]
gse <- gse[,sampleNames(gse)!="GSM300299"]
gse <- gse[,sampleNames(gse)!="GSM300300"]
gse <- gse[,sampleNames(gse)!="GSM300335"]
gse <- gse[,sampleNames(gse)!="GSM300338"]

gse <- gse[,sampleNames(gse)!="GSM300340"]
gse <- gse[,sampleNames(gse)!="GSM300341"]

gse <- gse[,sampleNames(gse)!="GSM318840"]
gse <- gse[,sampleNames(gse)!="GSM350078"]


for (i in 1:7)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",197+i, sep="")]  
}

for (i in 1:5)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",209+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",262+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",280+i, sep="")]
}
for (i in 1:3)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",174+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",178+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",193+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",205+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",215+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",219+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",223+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",227+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",231+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",235+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",239+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",268+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",258+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",272+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",276+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",286+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",290+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",294+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",301+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",305+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",309+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",313+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",317+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",321+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",325+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",329+i, sep="")]
  
}
for (i in 1:11)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",243+i, sep="")]  
}
for (i in 1:15)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1176",195+i, sep="")]  
}

for (i in 1:46)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1176",229+i, sep="")]  
}


# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 43),rep('AD',19)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCTRL-dAD, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE48350a_table.csv")
nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
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
#write.csv(sg,"Gene_GSE48350a.csv")
write.csv(ug,"GO_GSE48350a.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 14)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE48350a", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE48350a_correspondence.txt" )
}


#----------------------------------
# D) GSE48350b (entorhinal cortex)


# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE48350_series_matrix.txt.gz")

gse <- gse[,sampleNames(gse)!="GSM300182"]
gse <- gse[,sampleNames(gse)!="GSM300183"]
gse <- gse[,sampleNames(gse)!="GSM300185"]
gse <- gse[,sampleNames(gse)!="GSM300187"]
gse <- gse[,sampleNames(gse)!="GSM300188"]
gse <- gse[,sampleNames(gse)!="GSM300190"]
gse <- gse[,sampleNames(gse)!="GSM300191"]
gse <- gse[,sampleNames(gse)!="GSM300229"]
gse <- gse[,sampleNames(gse)!="GSM300259"]
gse <- gse[,sampleNames(gse)!="GSM300260"]
gse <- gse[,sampleNames(gse)!="GSM300298"]
gse <- gse[,sampleNames(gse)!="GSM300299"]
gse <- gse[,sampleNames(gse)!="GSM300333"]
gse <- gse[,sampleNames(gse)!="GSM300335"]
gse <- gse[,sampleNames(gse)!="GSM318840"]
gse <- gse[,sampleNames(gse)!="GSM350078"]

for (i in 1:7)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",165+i, sep="")]
  
  gse <- gse[,sampleNames(gse)!=paste("GSM300",196+i, sep="")]
}
for (i in 1:5)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",208+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",222+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",242+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",261+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",279+i, sep="")]
}
for (i in 1:4)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",252+i, sep="")]  
}
for (i in 1:3)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",173+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",177+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",181+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",192+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",204+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",214+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",218+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",230+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",234+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",238+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",248+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",267+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",271+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",275+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",285+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",289+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",293+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",300+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",304+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",308+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",312+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",316+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",320+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",324+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",328+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",338+i, sep="")]
}
for (i in 1:65)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1176",210+i, sep="")]  
}


# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 39),rep('AD',15)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCTRL-dAD, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE48350b_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
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
write.csv(ug,"GO_GSE48350b.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE48350b", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE48350b_correspondence.txt" )
}
#----------------------------------
# D) GSE48350c (superior frontal cortex)


# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE48350_series_matrix.txt.gz")

gse <- gse[,sampleNames(gse)!="GSM300166"]
gse <- gse[,sampleNames(gse)!="GSM300189"]
gse <- gse[,sampleNames(gse)!="GSM300190"]
gse <- gse[,sampleNames(gse)!="GSM300200"]
gse <- gse[,sampleNames(gse)!="GSM300202"]
gse <- gse[,sampleNames(gse)!="GSM300212"]
gse <- gse[,sampleNames(gse)!="GSM300246"]
gse <- gse[,sampleNames(gse)!="GSM300248"]
gse <- gse[,sampleNames(gse)!="GSM300249"]
gse <- gse[,sampleNames(gse)!="GSM300252"]
gse <- gse[,sampleNames(gse)!="GSM300253"]
gse <- gse[,sampleNames(gse)!="GSM300265"]
gse <- gse[,sampleNames(gse)!="GSM300283"]
gse <- gse[,sampleNames(gse)!="GSM300297"]
gse <- gse[,sampleNames(gse)!="GSM300298"]
gse <- gse[,sampleNames(gse)!="GSM300332"]
gse <- gse[,sampleNames(gse)!="GSM300333"]
for (i in 1:8)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",167+i, sep="")]
}
for (i in 1:7)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",180+i, sep="")]
}
for (i in 1:5)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",254+i, sep="")]
}
for (i in 1:3)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",176+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",191+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",195+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",203+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",207+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",213+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",217+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",221+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",225+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",229+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",233+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",237+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",241+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",260+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",266+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",270+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",274+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",278+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",284+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",288+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",292+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",299+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",303+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",307+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",311+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",315+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",319+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",323+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",327+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",337+i, sep="")]
}
for (i in 1:59)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1176",195+i, sep="")]  
}


# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 48),rep('AD',21)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCTRL-dAD, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE48350c_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 5)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE48350c.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE48350c", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE48350c_correspondence.txt" )
}

#-------------------------------------------------
# D) GSE48350d (post-central gyrus )


# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE48350_series_matrix.txt.gz")

gse <- gse[,sampleNames(gse)!="GSM300199"]
gse <- gse[,sampleNames(gse)!="GSM300201"]
gse <- gse[,sampleNames(gse)!="GSM300211"]
gse <- gse[,sampleNames(gse)!="GSM300225"]
gse <- gse[,sampleNames(gse)!="GSM300245"]
gse <- gse[,sampleNames(gse)!="GSM300247"]
gse <- gse[,sampleNames(gse)!="GSM300248"]
gse <- gse[,sampleNames(gse)!="GSM300254"]
gse <- gse[,sampleNames(gse)!="GSM300255"]
gse <- gse[,sampleNames(gse)!="GSM300258"]
gse <- gse[,sampleNames(gse)!="GSM300264"]
gse <- gse[,sampleNames(gse)!="GSM300282"]
gse <- gse[,sampleNames(gse)!="GSM300341"]
gse <- gse[,sampleNames(gse)!="GSM300335"]
gse <- gse[,sampleNames(gse)!="GSM300338"]
gse <- gse[,sampleNames(gse)!="GSM300339"]
gse <- gse[,sampleNames(gse)!="GSM318840"]
gse <- gse[,sampleNames(gse)!="GSM350078"]

for (i in 1:9)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",184+i, sep="")]
}
for (i in 1:8)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",166+i, sep="")]
}
for (i in 1:6)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",295+i, sep="")]
}
for (i in 1:5)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",226+i, sep="")]
}
for (i in 1:3)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM300",175+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",179+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",194+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",202+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",206+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",212+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",216+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",220+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",232+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",236+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",240+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",249+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",259+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",265+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",269+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",273+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",277+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",283+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",287+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",291+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",302+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",306+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",310+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",314+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",318+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",322+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",326+i, sep="")]
  gse <- gse[,sampleNames(gse)!=paste("GSM300",330+i, sep="")]  
}
for (i in 1:34)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1176",195+i, sep="")]  
}
for (i in 1:21)
{
  gse <- gse[,sampleNames(gse)!=paste("GSM1176",254+i, sep="")]  
}

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 43),rep('AD',25)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCTRL-dAD, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE48350d_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 5)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE48350d.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE48350d", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE48350d_correspondence.txt" )
}

#------------------------------------------------
# E. GSE4229


# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4229/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE4229_series_matrix.txt.gz")

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 11),rep('AD',8),rep('CTRL', 7),rep('AD',7),rep('CTRL', 2),rep('AD',3),rep('CTRL', 2)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dAD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","GENE_SYMBOL","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$GENE_SYMBOL

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE4229_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 5)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE4229.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 37)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE4229", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE4229_correspondence.txt" )
}

#------------------------------------------------
# F. GSE12685


# # 3. Download the GEO dataset 
# url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12685/matrix/"
# dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# dataset <- unlist(strsplit(dataset, "\r\n"))
# for (ds in dataset) {
#   download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
# }

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE12685_series_matrix.txt.gz")

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 8),rep('AD',6)))
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
write.csv(subtable_result,"GSE12685_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
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
write.csv(ug,"GO_GSE12685.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 26)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE12685", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE12685_correspondence.txt" )
}


#----------------------------------------------
# G. GSE4226


# 3. Download the GEO dataset 
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/GSE4226/matrix/"
dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
dataset <- unlist(strsplit(dataset, "\r\n"))
for (ds in dataset) {
  download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
}

# 4. Convert download dataset in usable class
gse <- getGEO(filename = "GSE4226_series_matrix.txt.gz")

# 5. Create matrix design for limma and calculate differential expression 
d <- factor(c(rep('CTRL', 7),rep('AD', 7), rep('CTRL', 7),rep('AD', 7)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dAD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

# 6. Create statistics table and subtable
table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","GENE_SYMBOL","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$GENE_SYMBOL

# 7. Save subtable in Excel format
write.csv(subtable_result,"GSE4226_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 8. Choose the LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 9. Create topGO class with annotation
AD_GOdata <- new("topGOdata",
                 description = "AD study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 2)

# 10. Show genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE4226.csv")

# 11. Perform Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

# 12. Compare the tests (not described)
#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

# 13. Create GO terms tree
allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 21)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE4226", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE4226_correspondence.txt" )
}

#--------------------------------------------------------------------------------------
# H) GSE5281

url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/"
dataset <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
dataset <- unlist(strsplit(dataset, "\r\n"))
for (ds in dataset) {
  download.file(paste(url,ds, sep = ""), paste(getwd(), "/",ds,sep = ""))
}

gse <- getGEO(filename = "GSE5281_series_matrix.txt.gz")
#gse <- gse[,sampleNames(gse)!="GSM2335908"]

d <- factor(c(rep('CTRL', 74),rep('AD', 87)))
mod <- model.matrix(~0+d)
fit_1 <- lmFit(gse, mod)
contr <- makeContrasts(dCTRL-dAD, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)

table_result <- topTable(fit_3, n = length(fit_3), adjust="BH", sort.by = "logFC")
subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
geneList <- subtable_result$logFC
names(geneList) <- subtable_result$Gene.Symbol

write.csv(subtable_result,"GSE5281_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

topDiffGenes <- function(allScore){
  return(abs(allScore)>1.0)
}

AD_GOdata <- new("topGOdata",
                 description = "Ad study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 7)

n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(AD_GOdata)
ug <- usedGO(AD_GOdata)
write.csv(ug,"GO_GSE5281.csv")

resultFisher <- runTest(AD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(AD_GOdata, algorithm = "classic", statistic = "ks")

#pvalFis <- score(resultFisher)
#pvalKS <- score(resultKS,whichGO = names(pvalFis))
#cor_pval <- cor(pvalFis,pvalKS)

allRes <- GenTable(AD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes =20)
showSigOfNodes(AD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(AD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "AD1_GSE5281", useInfo = "all", pdfSW = TRUE)

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
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "AD1_GSE5281_correspondence.txt" )
}