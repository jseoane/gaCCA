### genetic algorithm selection of best CCA genes subset


library("doMC")
library(genalg)


source("gene-based.cca2.R")
source("genealg3.R")
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

results.phen = foreach(i=1:phentype.l) %do%{

  cat(i)
  cat("\n")
  cat(phenotypes[i])
  cat("\n")
  phenName = phenotypes[i]
   
  # define fitness function (if not use forech, can be defined outside)
  # calculate the gene based association value for gene selection
  fitness3 = function(ids){
    verbose = F
  	
    if(all(ids==FALSE))
      return(30)
    else{
      ids = as.logical(ids)
      selGenes = genes[ids]
      snpsIdx = which(snps_locusID[,2]%in%selGenes)    
        #ix = as.logical(ids)
      return(gene.assocGenSel(as.matrix(inputData[,(1+snpsIdx)]),phenD[,(2+i)],verbose))
    }  
  }

  # total number of iterations of the GA
  it = 30

  # monitor functio to control the population homogeneity using hamming distance
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    #plot(obj, type="hist");
    cat(sum(hamming.distance(obj$population)))
    cat("\n")
    cat(minEval)
    cat("\n")
  }
  

  

  list2 = list()
  
   # population size for the gene selection approach
  populationSize = 600

  for(j in 1:reps){
	rbga.results = genalg3(size=l.genes,  zeroToOneRatio=700,iters=it, elitism= 35,
                        evalFunc=fitness3,popSize =populationSize , verbose=F,monitorFunc=monitor) #popsize=300
	
	cat(paste(phenName," ",rbga.results$best[it]))
	cat("\n")

    list1 = list()
    list1[["pvalue"]]= rbga.results$best[it]
    list1[["genesInd"]] = rbga.results$population[which.min(rbga.results$evaluations),]
    list1[["geneNames"]] = genes[as.logical(rbga.results$population[which.min(rbga.results$evaluations),])]
  
    list2[[j]]=list1


   }

  
   list2[["name"]]=phenName

   list2


}







