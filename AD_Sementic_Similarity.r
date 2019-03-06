# SCRIPT TWO: PATHWAY AND MATRICES (SEMANTIC SIMILIARY) 
# by Eugenio Del Prete

# 1. Set your work folder
#setwd(".../workfolder")
setwd("D:/Bioinformatics/With Ali Moni Aus/AD_Disease_Comorbidity/Output_txt")

# 2. Load libraries for script two
library(GOSemSim)
library(readtext)
library(stringr)
library(factoextra)
library(dendextend)
library(corrplot)
library(RColorBrewer)

# 3. List the available correspondence files 
#path <- c(".../Output_Txt")
path <- c("C:/Users/UtPaL/Dropbox/AD_Disease_Comorbidity/Output_txt")
files <- readtext(path)
files_split <- strsplit(files[,1],"_")
id <- unlist(lapply(files_split, function(x) paste0(x[1],"_",x[2])))

# 4. Create list of genes related to each dataset
go_term <- lapply(files[,2], function(x) str_extract_all(x,"GO:.{7}"))
go_term_sub <- lapply(go_term, function(x) x[[1]][1:5])
names(go_term_sub) <- id
gene_grasp <- function(text_gene){
    aux_1 <- str_replace_all(text_gene,"GO:.{7} |\n|genes:", "")
    aux_2 <- strsplit(aux_1, " ")
    aux_3 <- strsplit(aux_2[[1]][2:6],",")
    aux_4 <- unique(unlist(aux_3))
    return(aux_4)
}
gene_term_sub <- lapply(files[,2],gene_grasp)
names(gene_term_sub) <- id

# saveRDS(gene_term_sub, file = "gene_term_sub_5.Rda")

# 5. Select the gene ontology and prepare the annotation
hsGO <- godata('org.Hs.eg.db', ont="BP")
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)

# 6. Create semantic similarity matrices (GO terms, genes)
len <- length(id)
go_sem_sim_mat = gene_sem_sim_mat <- matrix(data = 0, nrow = len, ncol = len)
rownames(go_sem_sim_mat) = rownames(gene_sem_sim_mat) <- id
colnames(go_sem_sim_mat) = colnames(gene_sem_sim_mat) <- id

for(k in 1:len){
    for(kk in 1:len){
        go_sem_sim_mat[k,kk] <- mgoSim(go_term_sub[[k]], go_term_sub[[kk]], 
                                semData=hsGO, measure="Wang", combine="BMA")
        gene_sem_sim_mat[k,kk] <- clusterSim(gene_term_sub[[k]], gene_term_sub[[kk]], 
                                  semData=hsGO2, measure="Wang", combine="BMA")
        cat("k =",k,"kk =",kk,"per =",
            round(((k-1)*len+kk)/len^2,digits = 3)*100,"%\n")
    }
}

#saveRDS(go_sem_sim_mat,file="go_sem_sim_mat_5.Rda")
#saveRDS(gene_sem_sim_mat,file="gene_sem_sim_mat_5.Rda")


# 7. Plot the semantic similarity matrices
library(DOSE)

pdf(file = "go_mat_5.pdf", width = 10)
simplot(go_sem_sim_mat,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=2,
        font.size=12, xlab="", ylab="")
dev.off()

pdf(file = "gene_mat_5.pdf", width = 10)
simplot(gene_sem_sim_mat,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=2,
        font.size=12, xlab="", ylab="")
dev.off() 

do_id <-  c("DOID:10652","DOID:14330","DOID:12858","DOID:332","DOID:12377",
            "DOID:12217","DOID:9255", 
            "DOID:2377")
do_ac <- c("AD","PD","HD","ALS","SMA","LBD","FTD","MS")

do_sem_sim_mat <- doSim(do_id,do_id,measure="Wang")
rownames(do_sem_sim_mat) = colnames(do_sem_sim_mat) <- do_ac

#saveRDS(do_sem_sim_mat,file="do_sem_sim_mat.Rda")

pdf(file = "do_mat.pdf", width = 10)
simplot(do_sem_sim_mat,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=2,
        font.size=20, xlab="", ylab="")
dev.off()

# 8. Create KEGG Enrichment graph
library(clusterProfiler)

corrisp <- lapply(gene_term_sub,function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
corrisp <- lapply(corrisp,"[","ENTREZID")
corrisp <- lapply(corrisp, unlist)
names(corrisp) = substr(names(corrisp),1,4)

ccluster <- compareCluster(geneClusters = corrisp, fun="enrichKEGG")
pdf(file = "enrich_KEGG_5.pdf", width = 15)
dotplot(ccluster, font.size = 9)
dev.off()
