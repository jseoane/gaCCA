### Single genotype / single phenotype association

library("doMC")

source("gene-based.cca3.2.R")


# read the gene / snp mapping (from Beedtools)
snps_locusID = read.csv("gene_mapping.csv",header=T)

# remove SNPs non associated with any gene
idnil = which(snps_locusID[,2]=="")
snps_locusID = snps_locusID[-idnil,]

genes = unique(snps_locusID[,2])
l.genes = length(genes)

## generate filtered phenotype matrix


phenD = matrix(data = runif(286524), nrow = 3411, ncol = 84)
## mean 0 sd 1
phenD = scale(phenD)  


## generate genotype data 
## WARNING!!!! TO MUCH DATA, only use with high memory nodes
inputData = matrix(data=sample(3,102275424)-1,nrow=3411,ncol=29984)


# prepare parallel job
registerDoMC()
options(cores=2)
getDoParWorkers()


phenList = as.list(phenD[,3:ncol(phenD)])
pvaluesNew = foreach(i=1:l.genes, .combine = rbind) %dopar%{
    #tic()
    geneName = genes[i]
    #toc()
  
    #cat(paste("Gene ",i,": ",geneName,"\n"))
    #tic()
    snpsid = which(snps_locusID$genesymbol==geneName)
    #toc()
    #tic()
    dna  = as.matrix(inputData[,(1+snpsid)])
    #toc()
    #tic()
    keep.dna=check_mc(dna,F)
    #toc()
    #tic()
    assoc = lapply(phenList,function(x) gene.assoc(dna[,keep.dna],as.matrix(x)))
    #toc()
    
  }

colnames(pvaluesNew)=phenNames
rownames(pvaluesNew)=genes
