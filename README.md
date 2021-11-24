# C4-yasmeen
##tneurogenomics project 
---
title: "neurogenomics"
author: "yasmien"
date: "11/2/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
```{r}
## Background

This script will use the data from [Engram-specific transcriptome profiling of contextual memory consolidation](https://www.nature.com/articles/s41467-019-09960-x) to find differential gene expression in Engram cells between Fear-conditioned and control experiments.

## Setup for analysis
Packages that need to be installed and loaded are 

+ [tidyverse](https://www.tidyverse.org/packages/)

+ [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

+ [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

+ dendextend

+ reshape2

```{r loading packages, message=FALSE}
# installing packages from the CRAN repo
install.packages('tidyverse')
install.packages('reshape2')
install.packages('dendextend')
# Installing packages from the Bioconductor repo
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
library(tidyverse)
library(DESeq2)
library(edgeR)
library(reshape2)
library(dendextend)
```

```
```{r}
getwd()

```{r import meta data, message=FALSE}
# importing the metadata from the paper
SraRunTable <- read_csv("../Documents/05-PracticalDiffExp/data/SraRunTable.txt")

```
`# cleaning the column names
colnames(SraRunTable) = gsub(" \\(exp\\)", "", colnames(SraRunTable))
```{r}
# swapping out "-" and "+" for "minus" and "plus" becuase it will throw errors otherwise
SraRunTable$Cell_type = gsub("-", "_minus", SraRunTable$Cell_type)
SraRunTable$Cell_type = gsub("\\+", "_plus", SraRunTable$Cell_type)
head(SraRunTable)

```
```{r}
# extracting the GEO accession number for experiment identifier
# getting the list of all count files
file_list <- list.files(path="../Documents/05-PracticalDiffExp/data/counts", full.names = T)

accession = gsub('^.*../Documents/05-PracticalDiffExp/data/counts/\\s*|\\_.*$', '', file_list)
# reading in the gene list from the first count file
genes <- read.table(file_list[1], header=FALSE, sep="\t")[,1] 
# reading in the counts from all the files
counts    <- do.call(cbind,lapply(file_list,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))
colnames(counts) = accession
counts = data.frame(SYMBOL=genes,
                     counts)
head(counts)
```
 

`````{r}
# filter out the htseq stats 
counts = counts[!c(grepl("__no_feature", counts$SYMBOL)| 
                   grepl("__ambiguous", counts$SYMBOL)| 
                   grepl("__too_low_aQual", counts$SYMBOL)|  
                   grepl("__not_aligned", counts$SYMBOL)| 
                   grepl("__alignment_not_unique", counts$SYMBOL)),]
tail(counts)
```
```{r}
# adding read depths to metadata
metadata = data.frame(GEO_Accession = accession,
           depth = colSums(counts[,2:ncol(counts)]) ) %>% 
  left_join(SraRunTable) 
metadata
```

```{r, fig.width=12}
metadata %>% 
  ggplot(aes(x = accession, y = depth, fill = Mouse_ID)) +
  geom_col() +
  coord_flip()+
  facet_wrap(~Cell_type, scales = "free_y")+ 
  ggtitle("Cell Type")
```

```{r fig.width=12}
metadata %>% 
  ggplot(aes(x = accession, y = depth, fill = Mouse_ID)) +
  geom_col() +
  coord_flip()+
  facet_wrap(~Treatment, scales = "free_y")+ 
  ggtitle("Treatment")
```
 
```{r}
##--edgeR--##
# create edgeR object
dgList <- DGEList(counts=counts[,-1], 
                   genes=counts$SYMBOL, 
                   group = metadata$source_name  )
```
```{r}
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)
```

```{r}
countCheck <- countsPerMillion > 1
head(countCheck)
```

```{r}
# filter based on counts per million
keep <- which(rowSums(countCheck) >= 10)
dgList <- dgList[keep,]
genes.filt = genes[keep]
length(genes.filt)
``` 
 if I were to continue using edgeR I would create countmatrix with design and do differential analysis
## Differential Aanlysis
```{r}
##-- switching over to DESeq2--#
counts_filt = dgList$counts
# building the deseq object
dds<-DESeqDataSetFromMatrix(countData = counts_filt,
                            colData = metadata,
                           design = ~ Cell_type + Treatment + Cell_type:Treatment)

# running the deseq model
dds<-DESeq(dds)
```


```{r}
#exporting norm counts
normcounts <- counts(dds, normalized = T)
colnames(normcounts)<- metadata$source_name

write.csv(normcounts, "normalized_counts.csv")

#exporting deseq2 results

res<- results(dds, alpha = 0.05)

summary(res)
resordered= res[order(res$padj),]
write.csv(resordered, 'DESEq_resultsordered.csv')

#countsvslogfold change for verification there is diff genes
plotMA(dds, ylim = c(-6,6))
```

`### visualizing normalized counts
```{r}
vsd <- vst(dds, blind =FALSE)
mat <-assay(vsd)
head(as.data.frame(mat))
```
```{r}
#from here we are clustering data either by 1st hierac, or pca as a sanity check for our data
#add t in order to plot samples by genes and not genes by samples
dend = t(mat) %>% 
  scale %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram 
l = metadata$Treatment[ metadata$GEO_Accession %in% labels(dend)]
dend %>% 
  set("labels", l) %>% 
  plot

```






```{r}
#perform pca

ma.pca = prcomp(t(mat))
summary(ma.pca)
```


```{r} 
#we are making it a data frame in order to plot it nicely
#Extract PCA scrores
scores <-as.data.frame(ma.pca$x)
                   
scores
```

```{r, message=FALSE}
scores %>% 
  mutate(GEO_Accession = rownames(scores)) %>% 
  left_join(metadata) %>% 
  ggplot(aes(PC1, PC2, color = Treatment))+
  geom_point()
scores %>% 
  mutate(GEO_Accession = rownames(scores)) %>% 
  left_join(metadata) %>% 
  ggplot(aes(PC1, PC2, color = Cell_type))+
  geom_point()

```
````{r}
# Getting deseq results 
HomeCagevNonShock= results(dds, contrast = c("Treatment", "HomeCage", "Non Shock"), tidy = TRUE)
FearCondvNonShock= results(dds, contrast = c("Treatment", "Fear Conditioned", "Non Shock"), tidy = TRUE)
FearCondvHomeCage= results(dds, contrast = c("Treatment", "Fear Conditioned", "HomeCage"), tidy = TRUE)

```

```{r}
#reforming results
pval= 0.05
lfc= 1.5
# Home cage v Non Shock treatment
HomeCagevNonShock = HomeCagevNonShock %>%
 mutate(sig = ifelse(log2FoldChange > lfc & padj < pval , "UP", (ifelse(log2FoldChange < lfc & padj < pval, "DOWN", "notsig" )))) %>%
 mutate(SYMBOL= genes.filt)

head(HomeCagevNonShock)

```


```{r}
pval=.05
lfc=1.5
# Home cage v Non Shock treatment
 HomeCagevNonShock = HomeCagevNonShock %>% 
  mutate(sig = ifelse(log2FoldChange > lfc & padj < pval, "UP", (ifelse(log2FoldChange < -lfc & padj < pval, "DOWN", "not sig")) )) %>% 
  mutate(SYMBOL=genes.filt)
head(HomeCagevNonShock)
```

```{r}
#to see no. of genes that are up regulated, downregulated and not sig
table(HomeCagevNonShock$sig)

```
```{r}
#to see genes symbols that are up or down regulated
HomeCagevNonShock$SYMBOL[!(HomeCagevNonShock$sig =="not sig")] #used different results from what I came UP WIth

HomeCagevNonShock$SYMBOL[HomeCagevNonShock$sig == c("UP","DOWN")]



```
```{r}
plotMA(dds,contrast = c("Treatment", "Non Shock", "HomeCage"))
 

```

```{r}
# non shock v Fear conditioned
FearCondvNonShock <- FearCondvNonShock %>% 
  mutate(sig = ifelse(log2FoldChange > lfc & padj < pval, "UP", (ifelse(log2FoldChange < -lfc & padj < pval, "DOWN", "not sig")) ))%>% 
  mutate(SYMBOL=genes.filt)
table(FearCondvNonShock$sig)
```

```{r}
# home cage v fear conditioned
FearCondvHomeCage <- FearCondvHomeCage %>% 
  mutate(sig = ifelse(log2FoldChange > lfc & padj < pval, "UP", (ifelse(log2FoldChange < -lfc & padj < pval, "DOWN", "not sig")) ))%>% 
  mutate(SYMBOL=genes.filt)
table(FearCondvHomeCage$sig)
```
```{r}
#Volcanoplot
FearCondvHomeCage %>%
   filter(!(sig %in% NA)) %>% 
   mutate(l=ifelse(sig == c("UP","DOWN"), SYMBOL," ")) %>%
   ggplot(aes(x= log2FoldChange, y = -log(padj),color=sig, label = l)) +
   geom_point(alpha = 0.3) +
    geom_vline(xintercept = -1.5)+
    geom_vline(xintercept = 1.5)+
    geom_hline(yintercept = -log(0.05)) +
    geom_text()
```


```{r}
geom_text() 
```


```{r}
save.image("Neurotranscriptomics_diffanalysis")
```
```
```{r}
#functional annotation
title: "Practical Functional Annotation"
output: html_notebook
---

```{r loading packages, message=FALSE}
# installing packages from the CRAN repo
install.packages('tidyverse')
install.packages('reshape2')
install.packages('dendextend')
# Installing packages from the Bioconductor repo
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")
BiocManager::install("enrichplot")
# Import Libraries and Functions
library(tidyverse)
library(DESeq2)
library(edgeR)
library(reshape2)
library(dendextend)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
```

```{r}
metadata =
  metadata %>% 
  filter(Treatment %in% "Fear Conditioned")
  counts_filt = counts_filt[,metadata$GEO_Accession]


```
```{r}
#undergoing diff analysis between positive and negative cells

dds = DESeqDataSetFromMatrix(counts_filt, metadata, design =  ~ Cell_type)

dds = DESeq(dds)
```
```{r}
#visualizing normalized counts
 vsd = vst (dds, blind = FALSE)
 mat = assay(vsd)
 head (as.data.frame(mat))
```

```{r}
#CLUSTERING HIERACHIAL

dend = t(mat) %>% 
  scale %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram 
l = metadata$Cell_type[ metadata$GEO_Accession %in% labels(dend)]
dend %>% 
  set("labels", l) %>% 
  plot
  plot
```


```{r}
#clustering
 

pca.ma = prcomp(t(mat.x))
summary(pca.ma)
scores <-as.data.frame(ma.pca$x)
                   
scores


```

```{r}
#visualizing clustering
  
scores %>% 
  mutate(GEO_Accession = rownames(scores)) %>% 
  left_join(metadata) %>% 
  ggplot(aes(PC1, PC2, color = Cell_type))+
  geom_point()
scores

```

```{r}
minusvplus = results(dds,contrast = c("Cell_type","dVenus_plus","dVenus_minus"),tidy=T)
```

```{r}
lfc=1.5
pval = 0.5
minusvplus = minusvplus %>%
  mutate(sig = ifelse(log2FoldChange > lfc & padj< pval, "UP",(ifelse(log2FoldChange < -lfc & padj < pval, "DOWN", "notsig")))) %>%
             mutate(SYMBOL = genes.filt)
minusvplus
table(minusvplus$sig)
```
```{r}
minusvplus%>%
   filter(!(sig %in% NA)) %>% 
   mutate(l=ifelse(sig == c("UP","DOWN"), SYMBOL," ")) %>%
   ggplot(aes(x= log2FoldChange, y = -log(padj),color=sig, label = l)) +
   geom_point(alpha = 0.3) +
    geom_vline(xintercept = -1.5)+
    geom_vline(xintercept = 1.5)+
    geom_hline(yintercept = -log(0.05)) +
    geom_text()
```
```{r}
#annotATION GO, FIRST CREATING GENELIST
GOI= minusvplus$SYMBOL[minusvplus$sig == c("UP")]
genesUP.df = bitr(GOI, fromType = "SYMBOL", toType = c("ENSEMBL","SYMBOL","ENTREZID") , OrgDb = org.Mm.eg.db)
genesUP.df

 
```
## Functional Annotation
```{r}
GOI = minusvplus$SYMBOL[minusvplus$sig %in% c("UP")]
genesUP.df = bitr(GOI, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL", "ENTREZID"), OrgDb = org.Mm.eg.db)
genesUP.df
GOI = minusvplus$SYMBOL[minusvplus$sig %in% c("DOWN")]
genesDOWN.df = bitr(GOI, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL", "ENTREZID"), OrgDb = org.Mm.eg.db)
genesDOWN.df
GOI = minusvplus$SYMBOL[minusvplus$sig %in% c("UP", "DOWN")]
genesDE.df = bitr(GOI, fromType = "SYMBOL", toType = c("ENSEMBL", "SYMBOL", "ENTREZID"), OrgDb = org.Mm.eg.db)
genesDE.df
``````

```{r}
# GO Enrichment analysis 
GO_CC_u =  enrichGO(gene = genesUP.df$ENSEMBL,
         OrgDb = org.Mm.eg.db,
         keyType = "ENSEMBL",
         ont = "CC",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.01,
         qvalueCutoff = 0.05)
GO_CC_u
goplot(GO_CC_u)


```
```{r}
GO_MF_U= enrichGO(gene = genesUP.df$ENSEMBL,
             OrgDb = org.Mm.eg.db,
             keyType = "ENSEMBL",
             ont = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.01,
             qvalueCutoff = 0.05)
goplot(GO_MF_U)
```

```{r}
GO_BP_u =  enrichGO(gene = genesUP.df$ENSEMBL,
         OrgDb = org.Mm.eg.db,
         keyType = "ENSEMBL",
         ont = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.01,
         qvalueCutoff = 0.05)
goplot(GO_BP_u)
```
```{r, fig.width= 12, fig.height=12}
#reactomepathway analysis

s = minusvplus %>% 
  filter(sig %in% "DOWN") %>% 
  right_join(genesDOWN.df) %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneListDOWN = s$log2FoldChange
names(geneListDOWN)= s$ENTREZID
geneListUP <- sort(geneListDOWN, decreasing = TRUE)
head(geneListDOWN)

```

````{r, fig.width= 12, fig.height=12}
pathwayDOWN = enrichPathway(gene =genesDOWN.df$ENTREZID, 
              pvalueCutoff = 0.05, 
              readable=TRUE, 
              organism = "mouse") 

ReactomePA::dotplot(pathwayDOWN)
ReactomePA::cnetplot(pathwayDOWN, Foldchange = geneListDOWN)

```


````{r, fig.width= 12, fig.height=12}
# creating geneList input for pathway analysis
s = minusvplus %>% 
  filter(sig %in% "UP") %>% 
  right_join(genesUP.df) %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneListUP = s$log2FoldChange
names(geneListUP)= s$ENTREZID
geneListUP <- sort(geneListUP, decreasing = TRUE)
head(geneListUP)
pathwayUP = enrichPathway(gene =genesUP.df$ENTREZID, 
              pvalueCutoff = 0.05, 
              readable=TRUE, 
              organism = "mouse")
ReactomePA::dotplot(pathwayUP)
ReactomePA::cnetplot(pathwayUP, Foldchange = geneListDOWN)

```



