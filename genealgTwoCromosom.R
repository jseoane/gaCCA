# this function modifies the general behaviour of the original binary genetic algorithm
# of the package genalg with a parallel execution of all the evaluations using foreach (better performance)
# and also, this version evolve simultaneusly two independent populations. Note that each population has its own paramenters 

genalg3 = function (sizeX = 10,sizeY = 10, suggestions = NULL, popSize = 200, iters = 100, 
	  mutationChanceX = NA, elitism = NA, zeroToOneRatioX = 10,
          mutationChanceY = NA,  zeroToOneRatioY = 10, monitorFunc = NULL, 
          evalFunc = NULL, showSettings = FALSE, verbose = FALSE) 
{
  
  
  if (is.null(evalFunc)) {
    stop("A evaluation function must be provided. See the evalFunc parameter.")
  }
  varsX = sizeX
  varsY = sizeY
  if (is.na(mutationChanceX)) {
    mutationChanceX = 1/(varsX + 1)
  }
  if (is.na(mutationChanceY)) {
    mutationChanceY = 1/(varsY + 1)
  }

  if (is.na(elitism)) {
    elitism = floor(popSize/5)
  }

  if (verbose) 
    cat("Testing the sanity of parameters...\n")
  if (popSize < 5) {
    stop("The population size must be at least 5.")
  }
  if (iters < 1) {
    stop("The number of iterations must be at least 1.")
  }
  if (!(elitism < popSize )) {
    stop("The population size must be greater than the elitism.")
  }
  if (showSettings) {
    if (verbose) 
      cat("The start conditions:\n")
    result = list(sizeX = sizeX, sizeY = sizeY, suggestions = suggestions, 
                  popSize = popSize, iters = iters, elitism = elitism, 
                  mutationChanceX = mutationChanceX, mutationChanceY=mutationChanceY)
    class(result) = "rbga"
    cat(summary(result))
  }
  else {
    if (verbose) 
      cat("Not showing GA settings...\n")
  }
  if (varsX > 0 & varsY>0) {
    if (!is.null(suggestions)) {
      if (verbose) 
        cat("Adding suggestions to first population...\n")
      population = matrix(nrow = popSize, ncol = varsX+varsY)
      suggestionCount = dim(suggestions)[1]
      for (i in 1:suggestionCount) {
        population[i, ] = suggestions[i, ]
      }
      if (verbose) 
        cat("Filling others with random values in the given domains...\n")
      for (child in (suggestionCount + 1):popSize) {

        popX = sample(c(rep(0, zeroToOneRatioX), 
                                       1), varsX, rep = TRUE)
	popY = sample(c(rep(0, zeroToOneRatioY), 
                                       1), varsY, rep = TRUE)
        population[child, ] = c(popX,popY)
        while (sum(population[child, 1:varsX]) == 0 | sum(population[child, (varsX+1):(varsX+varsY)]) == 0 ) {
	        popX = sample(c(rep(0, zeroToOneRatioX), 
                                       1), varsX, rep = TRUE)
		popY = sample(c(rep(0, zeroToOneRatioY), 
                                       1), varsY, rep = TRUE)
        	population[child, ] = c(popX,popY)
 
        }
      }
    }
    else {
      if (verbose) 
        cat("Starting with random values in the given domains...\n")
      population = matrix(nrow = popSize, ncol = varsX+varsY)
      for (child in 1:popSize) {
        popX = sample(c(rep(0, zeroToOneRatioX), 
                                       1), varsX, rep = TRUE)
	popY = sample(c(rep(0, zeroToOneRatioY), 
                                       1), varsY, rep = TRUE)
        population[child, ] = c(popX,popY)
        while (sum(population[child, 1:varsX]) == 0 | sum(population[child, (varsX+1):(varsX+varsY)]) == 0 ) {
	        popX = sample(c(rep(0, zeroToOneRatioX), 
                                       1), varsX, rep = TRUE)
		popY = sample(c(rep(0, zeroToOneRatioY), 
                                       1), varsY, rep = TRUE)
        	population[child, ] = c(popX,popY)
 
        }
      }
    }
    bestEvals = rep(NA, iters)
    meanEvals = rep(NA, iters)
    evalVals = rep(NA, popSize)
    for (iter in 1:iters) {
      if (verbose) 
        cat(paste("Starting iteration", iter, "\n"))
      if (verbose) 
        cat("Calculating evaluation values... ")
      
      tmpVal = foreach(ii = 1:popSize,.combine=c)%do%{
        if (is.na(evalVals[ii])) {
            evalFunc(population[ii,]);
        }
        
      }
      idtmp = which(is.na(evalVals))
      evalVals[idtmp]=unlist(tmpVal)
      
          
      bestEvals[iter] = min(evalVals)
      meanEvals[iter] = mean(evalVals)
      if (verbose) 
        cat(" done.\n")
      if (!is.null(monitorFunc)) {
        if (verbose) 
          cat("Sending current state to rgba.monitor()...\n")
        result = list(type = "binary chromosome", sizeX = sizeX, sizeY,sizeY, 
                      popSize = popSize, iter = iter, iters = iters, 
                      population = population, elitism = elitism,
                      mutationChanceX = mutationChanceX,mutationChanceY = mutationChanceY, evaluations = evalVals, 
                      best = bestEvals, mean = meanEvals)
        class(result) = "rbga"
        monitorFunc(result)
      }
      if (iter < iters) {
        if (verbose) 
          cat("Creating next generation...\n")
        newPopulation = matrix(nrow = popSize, ncol = varsX+varsY)
        newEvalVals = rep(NA, popSize)
        if (verbose) 
          cat("  sorting results...\n")
        sortedEvaluations = sort(evalVals, index = TRUE)
        sortedPopulation = matrix(population[sortedEvaluations$ix, 
                                             ], ncol = varsX+varsY)
        if (elitism > 0) {
          if (verbose) 
            cat("  applying elitism...\n")
          newPopulation[1:elitism, ] = sortedPopulation[1:elitism, 
                                                        ]
          newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
        }
        if (varsX > 1) {  #varsX>1 and varsY >1??
          if (verbose) 
            cat("  applying crossover...\n")
          for (child in (elitism + 1):popSize) {
            parentProb = dnorm(1:popSize, mean = 0, sd = (popSize/3))
            parentIDs = sample(1:popSize, 2, prob = parentProb)
            parents = sortedPopulation[parentIDs, ]
            crossOverPointX = sample(0:varsX, 1)
	    crossOverPointY = sample(0:varsY, 1)
            if (crossOverPointX == 0) {
              childX = parents[2, 1:varsX]
              #newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
            }
            else if (crossOverPointX == varsX) {
              childX = parents[1, 1:varsX]              
#              newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
            }
	    else {
              childX = c(parents[1,1:varsX ][1:crossOverPointX], 
                                         parents[2,1:varsX ][(crossOverPointX + 1):varsX])
              
            }
 
            if (crossOverPointY == 0) {
              childY = parents[2, (varsX+1):(varsX+varsY)]
              #newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
            }
            else if (crossOverPointY == varsY) {
              childY = parents[1, (varsX+1):(varsX+varsY)]
              #newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
            }
	    else {
              childY = c(parents[1,(varsX+1):(varsX+varsY) ][1:crossOverPointY], 
                                         parents[2,(varsX+1):(varsX+varsY) ][(crossOverPointY + 1):varsY])
              
            }

	    while(sum(childX)==0){
		childX = sample(c(rep(0, zeroToOneRatioX), 
                                       1), varsX, rep = TRUE)
	    }
	    while(sum(childY)==0){
		childY = sample(c(rep(0, zeroToOneRatioY), 
                                       1), varsY, rep = TRUE)
	    }

	    newPopulation[child,] = c(childX,childY)


          } #close for
        }  #close if vars >1
        else {
          if (verbose) 
            cat("  cannot crossover (#vars=1), using new randoms...\n")
          newPopulation[(elitism + 1):popSize, ] = sortedPopulation[sample(1:popSize, 
                                                                           popSize - elitism), ]
        }
        population = newPopulation
        evalVals = newEvalVals

        if (mutationChanceX > 0) {
          if (verbose) 
            cat("  applying mutations in X... ")
          mutationCount = 0
          for (object in (elitism + 1):popSize) {
            for (var in 1:varsX) {
              if (runif(1) < mutationChanceX) {
                population[object, var] = sample(c(rep(0, 
                                                       zeroToOneRatioX), 1), 1)
                mutationCount = mutationCount + 1
              }
            }
          }
          if (verbose) 
            cat(paste(mutationCount, "mutations applied in X\n"))
        }

        if (mutationChanceY > 0) {
          if (verbose) 
            cat("  applying mutations in Y... ")
          mutationCount = 0
          for (object in (elitism + 1):popSize) {
            for (var in (varsX+1):(varsX+varsY)) {
              if (runif(1) < mutationChanceY) {
                population[object, var] = sample(c(rep(0, 
                                                       zeroToOneRatioY), 1), 1)
                mutationCount = mutationCount + 1
              }
            }
          }
          if (verbose) 
            cat(paste(mutationCount, "mutations applied in Y\n"))
        }

       

      }
    }
  }
  result = list(type = "binary chromosome", sizeX = sizeX, sizeY,sizeY,popSize = popSize, 
                iters = iters, suggestions = suggestions, population = population, 
                elitism = elitism, mutationChanceX = mutationChanceX, mutationChanceY =mutationChanceY,
		evaluations = evalVals,    best = bestEvals, mean = meanEvals)
  class(result) = "rbga"
  return(result)
}









