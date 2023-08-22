library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()
data <- read.delim ('/Users/sabinamahnesaei/Downloads/data.txt' , header = T)
geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,46:50)]



data[1:10,1:10]

data <- data %>% 
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>% 
  mutate(samples = gsub('\\.', '-', samples)) %>% 
  inner_join(., phenoData, by = c('samples' = 'title')) %>% 
  select(1,3,4) %>% 
  spread(key = 'geo_accession', value = 'counts') %>% 
  column_to_rownames(var = 'ENSEMBLID')

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


data <- data[gsg$goodGenes == TRUE,]
htree <- hclust(dist(t(data)), method = "average")
plot(htree)



pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) # 13284 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()




power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()





grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor




# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




library(dplyr)
df <- data.frame(x = 1:3, y = c("a", "b", "c"))

# Use the mutate function to add a new variable z
df <- df %>% mutate(z = x * 2)

# View the resulting data frame
df


trait <- colData %>%
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>%
  select(8)
  




#binerize cottegry column variab
colData$severity <- factor(colData$severity, levels = c("Convalescent", "Moderate", "Severe", "ICU" ,"Healthy"))

severity.out <-  binarizeCategoricalColumns(colData$severity,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)



trait<- cbind(trait, severity.out)
 #define number







nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)




module.trait.corr <- cor(module_eigengenes, trait, use = "p")
module.trait.corr.pvalues <- corPvalueStudent(module.trait.corr, nSamples)






# visualize module-trait association as a heatmap


heatmap.data <- merge(module_eigengenes, trait, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')


CorLevelPlot( heatmap.data,
              x= names(heatmap.data)[18:22],
              y = names(heatmap.data)[1:17],
              col = c("blue1", "skyblue", "white", "pink" , "red"))




#identifing the color of genes by use this
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>%
filter(`bwnet$colors` == 'turquoise') %>% 
rownames()  










 module.membership.measure <- cor(module_eigengenes, norm.counts, use = "p")
 module.membership.measure.pvalues <- corPvalueStudent(module.membership.measure, nSamples) 

 
 
 
 module.membership.measure.pvalues [1:10, 1:10] 

 
 
 
 
 
 # Calculate the gene significance and associated p-values
 
 
 
 
 
gene.signif.corr <- cor(norm.counts , trait$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvaluse <- corPvalueStudent(gene.signif.corr, nSamples)




gene.signf.corr.pvaluse %>%
as.data.frame() %>%  
arrange(V1) %>%
head(25)
