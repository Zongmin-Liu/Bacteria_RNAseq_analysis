## -----------------------------------------------------------
#Import R packages
library(WGCNA)
library(tidyverse)
library(readr)

##Import gene expression matrix, sample_info
gene_exp <- read.table("genes.TMM.EXPR.matrix", row.names=1)
sample_info <- read.table('sample_info.txt',
                          header = TRUE,row.names=1) %>%
  dplyr::select(-strain,-stage)
##Use tidyverse package to fiter data: select->col，filter->row
## -----------------------------------------------------------

emapper <- read_delim('query_seqs.fa.emapper.annotations', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
   dplyr::select(gene_id = X1, 
                gene_symbol = X9, 
                GO = X10, 
                KO = X12, 
                Pathway = X13, 
                OG = X7, 
                gene_name = X8)

gene_info <- dplyr::select(emapper,  gene_id, gene_symbol, gene_name) %>%
  dplyr::filter(!is.na(gene_name))


write.csv(gene_info, 
          file = 'GeneAnnotion.csv', 
          row.names = F, quote = T)

## -----------------------------------------------------------


## -----------------------------------------------------------
{
  # Transpose
  datExpr0 <- t(gene_exp)
  
  # Filtering of missing data and non-fluctuation data
  gsg <- goodSamplesGenes(
    datExpr0, 
    minFraction = 1/2 #Threshold for the proportion of missing data for genes
  )
  datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  
  # Through clustering, check whether there are obvious abnormal samples, and remove them if necessary
  plot(hclust(dist(datExpr)), 
       cex.lab = 1.5,
       cex.axis = 1.5, 
       cex.main = 2)
  
  #The phenotype must be a number in here, and the correlation coefficient need be calculated
  datExpr[1:4,1:4] 
}


## -----------------------------------------------------------
{
  datTraits <- sample_info
  datTraits[1:4, 1:4]
}

#Looking the best beta value
## ----message=FALSE------------------------------------------
{
  library(WGCNA)
  # Multithreading
  enableWGCNAThreads(nThreads = 6)
  ## disableWGCNAThreads()
  
  # Determine the best power through multiple iterations of power
  sft <- pickSoftThreshold(
    datExpr, 
    powerVector = 1:20, # Trying 1-20
    networkType = "signed"   #"unsigned","singed"
  )
}


## -----------------------------------------------------------
{
  # Drawing
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ggthemes)
  
  fig_power1 <- ggplot(data = sft$fitIndices,
                       aes(x = Power,
                           y = SFT.R.sq)) +
    geom_point(color = 'red') +
    geom_text_repel(aes(label = Power)) +
    geom_hline(aes(yintercept = 0.85), color = 'red') +
    labs(title = 'Scale independence',
         x = 'Soft Threshold (power)',
         y = 'Scale Free Topology Model Fit,signed R^2') +
    theme_few() +
    theme(plot.title = element_text(hjust = 0.5))
  
  fig_power2 <- ggplot(data = sft$fitIndices,
                       aes(x = Power,
                           y = mean.k.)) +
    geom_point(color = 'red') +
    geom_text_repel(aes(label = Power)) +
    labs(title = 'Mean connectivity',
         x = 'Soft Threshold (power)',
         y = 'Mean Connectivity') +
    theme_few()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_grid(fig_power1, fig_power2)
}

# Build the network
## ----------------------------------------------------
{
  net <- blockwiseModules(
    # 0. Input data
    datExpr, 
    
    # 1. Calculate the correlation coefficient
    corType = "pearson", # Correlation coefficient algorithm，pearson|bicor
    
    # 2. Calculate the adjacency matrix
    power = 4, # Obtained beta value
    networkType = "signed", # unsigned | signed | signed hybrid
    
    # 3. Calculate TOM matrix
    TOMType = "signed", # none | unsigned | signed
    saveTOMs = TRUE,
    saveTOMFileBase = "blockwiseTOM",
    
    # 4. Cluster and divide modules
    deepSplit = 4, # 0|1|2|3|4, The larger the value, the more but smaller modules
    minModuleSize = 30,  #30,50,80
    
    # 5. Merge similar modules
    mergeCutHeight = 0.25,  # 0.25.0.35
    numericLabels = FALSE, # Name modules by numbers
    nThreads = 0, # 0 All available threads are used
    maxBlockSize = 100000 # Need to be greater than the number of genes
  )
  # View the number of genes contained in each module
  table(net$colors) 
}

write.csv(table(net$colors), 
          file = 'output/tablecolor.csv', 
          row.names = F, quote = T)

## -----------------------------------------------------------
{
  library(tidyverse)
  wgcna_result <- data.frame(gene_id = names(net$colors),
                             module = net$colors) %>%
    left_join(gene_info, by = c('gene_id' = 'gene_id')) 
  head(wgcna_result)
}

#Save WGCNA analysis results, genes and color patches, annotation relationships
write.csv(wgcna_result, 
          file = 'output/wgcna_result.csv', 
          row.names = F, quote = T)

## -----------------Draw WGCNA color module and gene map------------------------------------------
{
  # Draw cluster diagram and module color at the same time
    plotDendroAndColors(
    dendro = net$dendrograms[[1]], 
    colors = net$colors,
    groupLabels = "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)
}



##Draw advanced heat maps
library(ggcor)
library(dplyr)

MEs <- net$MEs
colnames(MEs) <- str_remove(colnames(MEs), 'ME')
head(MEs)

quickcor(MEs,cor.test=T) +
  geom_square()

#Select three phenotypes: Time, pH, and biofilm ability
datTraits <- sample_info

## ---------------------------------------------------
{
  # Calculate correlation
  moduleTraitCor <- cor(
    net$MEs,
    datTraits,
    use = "p",
    method = 'pearson' # Pay attention to the calculation method of correlation coefficient 'spearman'
  )
  
  # Calculate Pvalue
  moduleTraitPvalue <- corPvalueStudent(
    moduleTraitCor, 
    nrow(datExpr))
}


## -----------------------------------------------------------
{
  # Correlation heatmap
  sizeGrWindow(10,6)
  
  # Connection correlation and pvalue
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  
  
  # heatmap 
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(net$MEs),
                 ySymbols = names(net$MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}


## -----------------------------------------------------------
{
  # Export target module "pink"
  my_modules <- c('pink') 
  
  # Extract the expression matrix of the module
  m_wgcna_result <- filter(wgcna_result, module %in% my_modules)
  write.csv(m_wgcna_result, 
            file = ' m_wgcna_result_pink.csv', 
            row.names = F, quote = T)
  
  m_datExpr <- datExpr[, m_wgcna_result$gene_id]
  
  # Calculate the TOM matrix of the module  
  m_TOM <- TOMsimilarityFromExpr(
    m_datExpr,
    power = 4, #Consistent with the previous value
    networkType = "signed",
    TOMType = "signed")
  
  dimnames(m_TOM) <- list(colnames(m_datExpr), colnames(m_datExpr))
  
  # Export as Cytoscape input file
  cyt <- exportNetworkToCytoscape(
    m_TOM,
    edgeFile = "CytoscapeInput-network-pink.txt",
    weighted = TRUE,
    threshold = 0.2) # threshold 0.2 0.3，0.4，0.5
}

# Input the cyt file to Cytoscape. Then use String database to constucte PPI network.
