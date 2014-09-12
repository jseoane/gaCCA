### genetic algorithm selection of best CCA phenotypes subset

source("gene-based.cca3.2.R")
library("imputation")
library("foreach")
library("CCA")

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




reps = 5 # number of rules to be extracted

# LOAD PARALLEL BACKGROUND
library("doMC")

library(genalg)
source("genealg3.R") # this code modifies the general behaviour of binary genetic algorithm


registerDoMC()
options(cores=2)
getDoParWorkers()

# iterate for each gene
results.genes = foreach(i=1:l.genes) %dopar%{
cat(genes[i])
geneName = genes[i]
cat("\n")

idSnp = which(snps_locusID$genesymbol==geneName)
geneD = as.matrix(inputData[,1+idSnp]) #retrieve genetic data of each gene


geneD.l = dim(geneD)[2]
phenD.l = dim(phenD)[2]

verbose = 0
nind=nrow(geneD)
nsnps=ncol(geneD)

keep.dna=check_mc(geneD,verbose) #remove colinearity
dna.pruned=geneD[,keep.dna]

# define fitness function (if not use forech, can be defined outside)
# calculate the gene based association value for phenotype selection
fitness3 = function(ids){  
  if(all(ids==FALSE))
    return(30)
  else{   
    ix = as.logical(ids)
    return(gene.assoc(dna.pruned,phenD[,ix]))
  }  
}

# calculate all single gene association for all the phenotypes
phen.pvalues = apply(phenD,2,function(f) gene.assoc(dna.pruned,f))
min.pvalue = min(phen.pvalues)

cat(paste(geneName," min uni pvalue: ",min.pvalue,"\n"))
list2 = list()

j=1 # index for repetitions (max rules extracted)
z=1 # index for, in case of non-converge, not infinite loop
it=100 # total number of iterations of the GA
maxIt = 20 # max value for z

while (j<=reps && z<maxIt){
  rbga.results = genalg3(size=phenD.l,  zeroToOneRatio=50,iters=it,elitism= 20,
                        evalFunc=fitness3,popSize =100 , verbose=F)  
  
  cat(paste(geneName," ",rbga.results$best[it]))
  cat("\n")
  if(rbga.results$best[it]<=min.pvalue){ #only use results higher than single max single point association
    list1 = list()
    list1[["pvalue"]]= rbga.results$best[it]
    list1[["phenInd"]] = rbga.results$population[which.min(rbga.results$evaluations),]
    list1[["phenNam"]] = colnames(phenotype[2:dim(phenotype)[2]])[as.logical(rbga.results$population[which.min(rbga.results$evaluations),])]
    #list1[["corY.Yscore"]]=cc(dna.pruned,phenD[,as.logical(rbga.results$population[which.min(rbga.results$evaluations),])])$scores
    list2[[j]]=list1
    j=j+1
  }  
  z = z+1
}
  
list2[["name"]]=geneName
list2[["minpvalue"]]=min.pvalue
list2

}


