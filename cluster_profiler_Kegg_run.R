library(clusterProfiler)
library(DOSE)
library(readxl)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project")
transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 1)
epistasis_transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 4)
data1<-read_excel("mart_export.xls")




geneList1 <- as.character(data1[,1]) #logFC
z<-as.character(data1[,3])
names(geneList1) = z #Gene ID
geneList <- sort(geneList1, decreasing = TRUE)
#use head to check gene list 
######################################
geneList1<-data1$`NCBI gene (formerly Entrezgene) ID`
z<-data1$`Gene stable ID`
names(geneList1) = z 
geneList <- sort(geneList1, decreasing = TRUE)
###############

k<-gseKEGG(geneList, organism = "sce")#, #keyType = 'ncbi-geneid', pvalueCutoff = 0.01)

##############################################33
gene_name<-vector()
gene_score<-vector()
for(i in 1:length(data1$`Gene stable ID`)){
  x<-grep(data1$`Gene stable ID`[i], epistasis_transcriptomics$`Systematic Name`)
  if(length(x) != 0 ){
    gene_name[i]<-data1$`Gene stable ID`[i]
    gene_score[i]<-epistasis_transcriptomics$`HLUM=LUM+H`[x]
  }

}

HLUM_LUM_H<-gene_score

#Kegg Enrichment 
shrinkLvV<-cbind(data1,HLUM_LUM_H)

sigGenes <- shrinkLvV$`Gene stable ID`[ shrinkLvV$HLUM_LUM_H  & 
                                      !is.na(shrinkLvV$HLUM_LUM_H) & 
                                      abs(shrinkLvV$HLUM_LUM_H)  ]#> 1 ]
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'sce') #vitis vinifera vvi & arabidopsis ath
head(kk, n=200)

write.table(kk,file = "Yeast_HLUM=LUM+H_pathway.csv", row.names= TRUE, sep=',')

browseKEGG(kk,'sce04111')
browseKEGG(kk,'sce01110')

library(pathview)
Yeast<-read.csv("Yeast_HLUM=LUM+H_pathway.csv")
shrinkLvV<-cbind(data1,HLUM_LUM_H)
logFC <- shrinkLvV$HLUM_LUM_H
names(logFC) <- shrinkLvV$`NCBI gene (formerly Entrezgene) ID`

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HLUM=LUM+H")
for(i in 1:length(Yeast$ID)){
  pathview(gene.data = logFC, 
           pathway.id = as.character(Yeast$ID[i]), 
           species = "sce")#,
  #limit = list(gene=5, cpd=1))
}
###################################################################################################################
############################## Sepcies 2 ##########################################################################
####################################################################################################################
library(clusterProfiler)
library(DOSE)
library(readxl)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project")
transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 1)
epistasis_transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 4)
data1<-read_excel("mart_export.xls")



geneList1 <- as.character(data1[,1]) #logFC
z<-as.character(data1[,3])
names(geneList1) = z #Gene ID
geneList <- sort(geneList1, decreasing = TRUE)
#use head to check gene list 
######################################
geneList1<-data1$`NCBI gene (formerly Entrezgene) ID`
z<-data1$`Gene stable ID`
names(geneList1) = z 
geneList <- sort(geneList1, decreasing = TRUE)
###############

k<-gseKEGG(geneList, organism = "sce")#, #keyType = 'ncbi-geneid', pvalueCutoff = 0.01)

##############################################33
gene_name<-vector()
gene_score<-vector()
for(i in 1:length(data1$`Gene stable ID`)){
  x<-grep(data1$`Gene stable ID`[i], epistasis_transcriptomics$`Systematic Name`)
  if(length(x) != 0 ){
    gene_name[i]<-data1$`Gene stable ID`[i]
    gene_score[i]<-epistasis_transcriptomics$`HLM=HM+L`[x]
  }
  
}

HLUM_LUM_H<-gene_score

#Kegg Enrichment 
shrinkLvV<-cbind(data1,HLUM_LUM_H)

sigGenes <- shrinkLvV$`Gene stable ID`[ shrinkLvV$HLUM_LUM_H  & 
                                          !is.na(shrinkLvV$HLUM_LUM_H) & 
                                          abs(shrinkLvV$HLUM_LUM_H)  ]#> 1 ]
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'sce') #vitis vinifera vvi & arabidopsis ath
head(kk, n=200)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HLM=HM+L")

write.table(kk,file = "Yeast_HLM=HM+L_pathway.csv", row.names= TRUE, sep=',')

#browseKEGG(kk,'sce04111')
#browseKEGG(kk,'sce01110')

library(pathview)
Yeast<-read.csv("Yeast_HLM=HM+L_pathway.csv")
shrinkLvV<-cbind(data1,HLUM_LUM_H)
logFC <- shrinkLvV$HLUM_LUM_H
names(logFC) <- shrinkLvV$`NCBI gene (formerly Entrezgene) ID`

#setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HLUM=LUM+H")
for(i in 1:length(Yeast$ID)){
  pathview(gene.data = logFC, 
           pathway.id = as.character(Yeast$ID[i]), 
           species = "sce")#,
  #limit = list(gene=5, cpd=1))
}
########################################################################################################################
################################################ Species 3 ############################################################
#######################################################################################################################

library(clusterProfiler)
library(DOSE)
library(readxl)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project")
transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 1)
epistasis_transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 4)
data1<-read_excel("mart_export.xls")



geneList1 <- as.character(data1[,1]) #logFC
z<-as.character(data1[,3])
names(geneList1) = z #Gene ID
geneList <- sort(geneList1, decreasing = TRUE)
#use head to check gene list 
######################################
geneList1<-data1$`NCBI gene (formerly Entrezgene) ID`
z<-data1$`Gene stable ID`
names(geneList1) = z 
geneList <- sort(geneList1, decreasing = TRUE)
###############

k<-gseKEGG(geneList, organism = "sce")#, #keyType = 'ncbi-geneid', pvalueCutoff = 0.01)

##############################################33
gene_name<-vector()
gene_score<-vector()
for(i in 1:length(data1$`Gene stable ID`)){
  x<-grep(data1$`Gene stable ID`[i], epistasis_transcriptomics$`Systematic Name`)
  if(length(x) != 0 ){
    gene_name[i]<-data1$`Gene stable ID`[i]
    gene_score[i]<-epistasis_transcriptomics$`HM=M+H`[x]
  }
  
}

HLUM_LUM_H<-gene_score

#Kegg Enrichment 
shrinkLvV<-cbind(data1,HLUM_LUM_H)

sigGenes <- shrinkLvV$`Gene stable ID`[ shrinkLvV$HLUM_LUM_H  & 
                                          !is.na(shrinkLvV$HLUM_LUM_H) & 
                                          abs(shrinkLvV$HLUM_LUM_H)  ]#> 1 ]
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'sce') #vitis vinifera vvi & arabidopsis ath
head(kk, n=200)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HM=M+H")

write.table(kk,file = "Yeast_HM=M+H_pathway.csv", row.names= TRUE, sep=',')

#browseKEGG(kk,'sce04111')
#browseKEGG(kk,'sce01110')

library(pathview)
Yeast<-read.csv("Yeast_HM=M+H_pathway.csv")
shrinkLvV<-cbind(data1,HLUM_LUM_H)
logFC <- shrinkLvV$HLUM_LUM_H
names(logFC) <- shrinkLvV$`NCBI gene (formerly Entrezgene) ID`

#setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HLUM=LUM+H")
for(i in 1:length(Yeast$ID)){
  pathview(gene.data = logFC, 
           pathway.id = as.character(Yeast$ID[i]), 
           species = "sce")#,
  #limit = list(gene=5, cpd=1))
}

########################################################################################################################
################################################ Species 4 ############################################################
#######################################################################################################################

library(clusterProfiler)
library(DOSE)
library(readxl)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project")
transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 1)
epistasis_transcriptomics<-read_excel("nmicrobiol201530-s3.xlsx", sheet = 4)
data1<-read_excel("mart_export.xls")



geneList1 <- as.character(data1[,1]) #logFC
z<-as.character(data1[,3])
names(geneList1) = z #Gene ID
geneList <- sort(geneList1, decreasing = TRUE)
#use head to check gene list 
######################################
geneList1<-data1$`NCBI gene (formerly Entrezgene) ID`
z<-data1$`Gene stable ID`
names(geneList1) = z 
geneList <- sort(geneList1, decreasing = TRUE)
###############

k<-gseKEGG(geneList, organism = "sce")#, #keyType = 'ncbi-geneid', pvalueCutoff = 0.01)

##############################################33
gene_name<-vector()
gene_score<-vector()
for(i in 1:length(data1$`Gene stable ID`)){
  x<-grep(data1$`Gene stable ID`[i], epistasis_transcriptomics$`Systematic Name`)
  if(length(x) != 0 ){
    gene_name[i]<-data1$`Gene stable ID`[i]
    gene_score[i]<-epistasis_transcriptomics$`HU=HUM-M`[x]
  }
  
}

HLUM_LUM_H<-gene_score

#Kegg Enrichment 
shrinkLvV<-cbind(data1,HLUM_LUM_H)

sigGenes <- shrinkLvV$`Gene stable ID`[ shrinkLvV$HLUM_LUM_H  & 
                                          !is.na(shrinkLvV$HLUM_LUM_H) & 
                                          abs(shrinkLvV$HLUM_LUM_H)  ]#> 1 ]
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'sce') #vitis vinifera vvi & arabidopsis ath
head(kk, n=200)

setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HU=HUM-M")

write.table(kk,file = "Yeast_HU=HUM-M_pathway.csv", row.names= TRUE, sep=',')

#browseKEGG(kk,'sce04111')
#browseKEGG(kk,'sce01110')

library(pathview)
Yeast<-read.csv("Yeast_HU=HUM-M_pathway.csv")
shrinkLvV<-cbind(data1,HLUM_LUM_H)
logFC <- shrinkLvV$HLUM_LUM_H
names(logFC) <- shrinkLvV$`NCBI gene (formerly Entrezgene) ID`

#setwd("~/Documents/Iowa_State_University/Courses/BCB 570/Project/BCB-570-Project/HLUM=LUM+H")
for(i in 1:length(Yeast$ID)){
  pathview(gene.data = logFC, 
           pathway.id = as.character(Yeast$ID[i]), 
           species = "sce")#,
  #limit = list(gene=5, cpd=1))
}

