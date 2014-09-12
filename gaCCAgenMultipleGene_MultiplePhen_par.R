### genetic algorithm selection of best CCA genes/phen subsets
## All Gen / all phen association


library("doMC")
library(genalg)

source("genealgTwoCromosom.R")
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


reps = 5 # number of rules to be extracted

# define fitness function 
# calculate the gene based association value for gene /phenotype selection
  fitness3 = function(ids){
    dimX = l.genes
    dimY = phentype.l

    verbose = F
  	
      ids = as.logical(ids)
      idsX = ids[1:dimX]
      idsY = ids[(dimX+1):(dimX+dimY)]
    if(all(idsX==FALSE) | all(idsY==FALSE))
	return(30)
    else{
      selGenes = genes[idsX]
      snpsIdx = which(snps_locusID[,2]%in%selGenes)    

        #ix = as.logical(ids)
      return(gene.assocGenSel(as.matrix(inputData[idgen,(1+snpsIdx)]),as.matrix(phenD[,(2+which(idsY))]),verbose))
    }  
  }

 it = 60  # total number of iterations of the GA

# monitor functio to control the population homogeneity using hamming distance
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    #plot(obj, type="hist");
    pop = obj$population
    cat(sum(hamming.distance(pop)))
    cat("\n")
    cat(minEval)
    cat("\n")
    bestpop = pop[which.min(obj$evaluations),]
    cat(genes[as.logical(bestpop[1:l.genes])])
    cat("\n")
    cat(phenotypes[as.logical(bestpop[(l.genes+1):(l.genes+phentype.l)])])
    cat("\n")
  }
  

list2 = list()
  
# population size for the gene/phenotype selection approach
populationSize = 1000
mutationChanceX = 1/l.genes  #number of total genes dependent mutation chance for gene selection

mutationChanceY= 1/phentype.l #number of total phenotype dependent mutation chance for phen selection

for(j in 1:reps){
	rbga.results = genalg3(sizeX=l.genes, sizeY=phentype.l,  zeroToOneRatioX=700,zeroToOneRatioY=50,iters=it, elitism= 100, evalFunc=fitness3,popSize =populationSize , verbose=F,monitorFunc=monitor) 
	
    #cat(paste(phenName," ",rbga.results$best[it]))
    #cat("\n")

    list1 = list()
    list1[["pvalue"]]= rbga.results$best[it]
 
    bestInd = rbga.results$population[which.min(rbga.results$evaluations),]
    
    geneInd = bestInd[1:l.genes]
    list1[["genesInd"]] = geneInd
    list1[["geneNames"]] = genes[as.logical(geneInd)]

    phenInd = bestInd[(l.genes+1):(l.genes+phentype.l)]
    list1[["phenInd"]] = phenInd
    list1[["phenNames"]] = phenotypes[as.logical(phenInd)]
  
    list2[[j]]=list1

 
}

 
