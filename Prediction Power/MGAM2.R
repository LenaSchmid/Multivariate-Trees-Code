library(mvtnorm)
options(java.parameters = "-Xmx8g")
library(extraTrees)
#### Tree and Forest Functions ####

{
  forestfit <- function(Y,                          ## -> 1 or multiple targetvectors (dim = n X #targets)
                        X,                            ## -> 1 or multiple features (dim = n X #features)
                        ntree = 2,                 ## -> number of trees fitted
                        bos=TRUE,                    ## bootstrap sample or not (default TRUE)
                        min.node.size = 5,            ## -> minimal node size (5 for regression is standard; 
                        ## for other tasks you should alter the default)
                        mtry = NA,                    ## -> number of features used for determining split
                        num.random.cuts = NA,         ## -> number of random cuts per sampled feature
                        cp = 0.01,                     ## -> complexity parameter for minimal increase in 
                        ## splitting criterion to attempt a split
                        split="L2"
  ){
    
    ##
    ## Produce a list with multiple trees based on (different bootstrap samples of)
    ## the original data
    ##
    forest <- lapply(1:ntree,
                     FUN = function(bs){
                       ##
                       ## Sample IDs for bootstrap
                       if(bos){
                         bootstrap.ids <- sample(1:nrow(X), size = nrow(X), replace = TRUE)
                       }
                       else{
                         bootstrap.ids<-1:nrow(X)
                       }
                       
                       ##
                       ## Fit tree on bootstrap sample
                       ##
                       tree <- treefit(Y = data.frame(Y[bootstrap.ids,]),                           
                                       X = data.frame(X[bootstrap.ids,]),                           
                                       min.node.size = min.node.size,           
                                       mtry = mtry,                   
                                       num.random.cuts = num.random.cuts,        
                                       cp = cp,
                                       split=split
                       )
                       
                     })
    
    forest
    
  }
  
  node.split <- function(ID, 
                         tree, 
                         Y, 
                         X, 
                         min.node.size, 
                         mtry, 
                         num.random.cuts,
                         cp,
                         split){
    
    ##
    ## check if current situation is terminal by min.node.size (so if the number of
    ## rows of the data is the min.node.size), if all Y are equal within theirselves 
    ## (so if the variance of all targets is 0)
    if(nrow(X) == min.node.size |
       sum(var(Y)) == 0){
      best.split <- NULL
      
    }else{
      
      ####
      #### Determine possible Splitpoints ####
      ####
      splitpoints <- data.frame()
      
      ##
      ## Determine possible features for the split
      ## If all features should be looked at, then mtry is still NA
      ##
      
      if(is.na(mtry)){
        
        
        splitfeatures <- 1:ncol(X)
        
      }else{
        splitfeatures <- sample(1:ncol(X), mtry)
        ##
        ## [ATTENTION]: features, which can not be splitted anymore, can be drawn
        ## This could result in producing a terminal node, even though more splits
        ## could have been done! To avoid this it should first be checked here, which of 
        ## the features can be splitted and which can not and then sample the numbers
        ## of only the features, which can be further splitted
        ##
        
      }
      ##
      ## This part is only used for the features which are taken as possible 
      ## candidates for a split
      ##
      for(i in splitfeatures){
        
        ##
        ## For this study only numerical values are used
        ##
        ## For categorical values it has to be checked, if X[,i] is a factor and 
        ## then all possible partitions have to be added (or for an ordered factor 
        ## all possible splitpoints have to be added) to the possible splitpoints
        ##
        
        ##
        ## Check if a split is possible in feature i (for categorical feature,
        ## check if feature i has more than 1 level)
        ##
        if(min(X[,i]) != max(X[,i])){
          
          
          ##
          ## If all possible splits of the chosen features should be looked at, then 
          ## num.random.cuts is still NA
          ##
          if(is.na(num.random.cuts)){
            
            ##
            ## Store and order all unique values of feature i
            ##
            unique.values <- unique(X[,i])
            unique.values <- unique.values[order(unique.values)]
            
            ##
            ## Add all possible splitpoints of the chosen feature into the list of 
            ## possible splits
            ##
            splitpoints <- rbind(splitpoints,
                                 data.frame(
                                   feature = i,
                                   Values = (unique.values[1:(length(unique.values)-1)] + unique.values[2:length(unique.values)])/2 ## Mean between two values
                                 )
            )
          }else{
            
            ##
            ## For extraTree method draw num.random.cuts possible cutpoints 
            ## uniformly distributed from min to max from feature i
            ##
            splitpoints <- rbind(splitpoints,
                                 data.frame(
                                   feature = i,
                                   Values = runif(num.random.cuts, 
                                                  min = min(X[, i]),
                                                  max = max(X[, i])) ## Random choice of splitvalue in feature range
                                 )
            )
            
          }## end else from determining the splitting points
          
        }## end the if-check if a split is possible (min(X[,i]) = max(X[,i])?)
        
      }## end for-loop for all candidates for a split
      ####
      #### Evaluate possible splits ####
      ####
      
      ##
      ## Check if there are possible split points and then evaluate them
      ##
      if(nrow(splitpoints) > 0){
        
        ##
        ## Save current L2-loss
        ## It is calculated as a sum of the L2-loss over all targets, because in this
        ## project multiple targets should be looked at, but it also works for a single
        ## target. If another lossfunction should be used (especially in the case of 
        ## not having a regression task)
        ##
        ## Here it is calculated as the sum of all the mean squared errors over the 
        ## different targets in Y in the current node, where the split is looked for
        ## (commented code, for smoothly switching the criterion)
        ##
        # currentloss <- sum(
        #   ##
        #   ## Apply the mean squared errors over all columns of Y
        #   ##
        #   vapply(1:ncol(Y), FUN = function(ycol){
        #     sum(
        #       (Y[,ycol] - mean(Y[,ycol]))^2
        #     )/nrow(Y)
        #   },
        #   FUN.VALUE = 0
        #   )
        # )
        
        #L2 loss
        if(split=="L2"){
          ##
          ## For the L2-loss the calculation can be done more easy
          ##
          ##
          ## Calculate the sum of the trace of the variance/covariance matrix
          ## and multiply by (n-1) to get sum of squared errors
          ##
          currentloss <- sum(diag(var(Y)))*((nrow(Y)-1))/nrow(Y)
          
          ##
          ## Initialize the column for the measured L2-loss after the split
          ##
          splitpoints$l.two.loss <- 0
          ##
          ## Initialize column to measure how many data points will be in left childnode
          ##
          splitpoints$obs.left <- 0
          
          
          templist <- lapply(1:nrow(splitpoints), function(split){
            ##
            ## Check TRUE for all rows which have a value of the current feature 
            ## smaller or equal the value where the split should be tested
            ##
            obsleftnode <- X[,splitpoints$feature[split]] <= splitpoints$Values[split]
            
            ##
            ## Produce list with the number of observations in the possible left node
            ## and the sum of squared errors for each split
            ##
            list(
              sum(obsleftnode),
              
              ##
              ## Apply the mean squared error over all columns of Y and sum up afterwards
              ##
              sum(
                vapply(seq(1,ncol(Y)), FUN = function(ycol){
                  ##
                  ## Add the mean squared error weighted by the number of observations
                  ## in the respective node. This is done by adding up all squared errors
                  ## within a node up and in the end dividing by the total number of
                  ## observations in both nodes together.
                  ##
                  
                  ##
                  ## Calculate sum or squared errors for left childnode divided by the
                  ## total number of observations
                  ##
                  sum(
                    (Y[obsleftnode, ycol] - mean(Y[obsleftnode, ycol]))^2
                  )/nrow(Y) +
                    
                    ##
                    ## Calculate sum or squared errors for right childnode divided by the
                    ## total number of observations
                    ##
                    sum(
                      (Y[!obsleftnode, ycol] - mean(Y[!obsleftnode, ycol]))^2
                    )/nrow(Y)
                },
                FUN.VALUE = 0
                )## close apply for calculating a vector of the L2-loss's of the targets for a specific split
              )## close sum for adding the MSEs for the different targets
            )
          }
          )## close apply for calculating all L2-losses for the different possible splits
        }## close if for split criterion
        
        ## L1 loss with mean
        if(split=="L1"){
          
          currentloss <- sum(
            ##
            ## Apply the L1 errors over all columns of Y
            ##
            vapply(1:ncol(Y), FUN = function(ycol){
              sum(
                abs(Y[,ycol] - mean(Y[,ycol]))
              )/nrow(Y)
            },
            FUN.VALUE = 0
            )
          )
          
          
          ##
          ## Initialize the column for the measured L2-loss after the split
          ##
          splitpoints$l.two.loss <- 0
          ##
          ## Initialize column to measure how many data points will be in left childnode
          ##
          splitpoints$obs.left <- 0
          
          
          templist <- lapply(1:nrow(splitpoints), function(split){
            ##
            ## Check TRUE for all rows which have a value of the current feature 
            ## smaller or equal the value where the split should be tested
            ##
            obsleftnode <- X[,splitpoints$feature[split]] <= splitpoints$Values[split]
            
            ##
            ## Produce list with the number of observations in the possible left node
            ## and the sum of squared errors for each split
            ##
            list(
              sum(obsleftnode),
              
              
              sum(
                vapply(seq(1,ncol(Y)), FUN = function(ycol){
                  
                  sum(
                    abs(Y[obsleftnode, ycol] - mean(Y[obsleftnode, ycol]))
                  )/nrow(Y) +
                    
                    
                    sum(
                      abs(Y[!obsleftnode, ycol] - mean(Y[!obsleftnode, ycol]))
                    )/nrow(Y)
                },
                FUN.VALUE = 0
                )## close apply for calculating a vector of the L2-loss's of the targets for a specific split
              )## close sum for adding the MSEs for the different targets
            )
          }
          )## close apply for calculating all L1-losses for the different possible splits
        }##close if for split criterion
        
        
        ## L1 loss with median
        if(split=="L1m"){
          
          currentloss <- sum(
            ##
            ## Apply the L1 errors over all columns of Y
            ##
            vapply(1:ncol(Y), FUN = function(ycol){
              sum(
                abs(Y[,ycol] - median(Y[,ycol]))
              )/nrow(Y)
            },
            FUN.VALUE = 0
            )
          )
          
          
          ##
          ## Initialize the column for the measured L2-loss after the split
          ##
          splitpoints$l.two.loss <- 0
          ##
          ## Initialize column to measure how many data points will be in left childnode
          ##
          splitpoints$obs.left <- 0
          
          
          templist <- lapply(1:nrow(splitpoints), function(split){
            ##
            ## Check TRUE for all rows which have a value of the current feature 
            ## smaller or equal the value where the split should be tested
            ##
            obsleftnode <- X[,splitpoints$feature[split]] <= splitpoints$Values[split]
            
            ##
            ## Produce list with the number of observations in the possible left node
            ## and the sum of squared errors for each split
            ##
            list(
              sum(obsleftnode),
              
              
              sum(
                vapply(seq(1,ncol(Y)), FUN = function(ycol){
                  
                  sum(
                    abs(Y[obsleftnode, ycol] - median(Y[obsleftnode, ycol]))
                  )/nrow(Y) +
                    
                    
                    sum(
                      abs(Y[!obsleftnode, ycol] - median(Y[!obsleftnode, ycol]))
                    )/nrow(Y)
                },
                FUN.VALUE = 0
                )## close apply for calculating a vector of the L2-loss's of the targets for a specific split
              )## close sum for adding the MSEs for the different targets
            )
          }
          )## close apply for calculating all L1-losses for the different possible splits
        }##close if for split criterion
        
        
        ##
        ## Flatten the resulting list. Every second entry refers to the L2-loss of a split, the others to
        ## the number of observations in the left node
        ##
        splitpoints$l.two.loss <- unlist(templist)[seq(2, nrow(splitpoints)*2, 2)]
        splitpoints$obs.left <- unlist(templist)[seq(1, nrow(splitpoints)*2, 2)]
        
        ####
        #### Check the best splits for possible usage ####
        #### (depending for example on min.node.size)
        ####
        
        ##
        ## Check every split starting with the one with the best increase in the L2-loss
        ## until the best possible split is found
        ##
        bool <- TRUE
        while(bool){
          
          ##
          ## Take best found split
          ##
          split.tested <- splitpoints[which.min(splitpoints$l.two.loss),]
          
          ##
          ## Check if split can be done
          ##
          if(split.tested$obs.left < min.node.size |
             (nrow(Y) - split.tested$obs.left) < min.node.size |
             split.tested$l.two.loss > (1 - cp)*currentloss){
            
            ##
            ## Remove if the split can not be done
            ##
            splitpoints <- splitpoints[-which.min(splitpoints$l.two.loss),]
            
            ##
            ## If no more splits are in the list or if the increase in the criterion for 
            ## even the best found split is to bad, then store no possible split found
            ## (which is setting best.split as NA)
            ##
            if(nrow(splitpoints) == 0 |
               split.tested$l.two.loss > (1 - cp)*currentloss){
              bool <- FALSE
              best.split <- NULL
            }
            
          }else{
            
            bool <- FALSE
            
            best.split <- split.tested
            
          }
          
        }## end while-loop for going over the splits
        
        
      }else{
        
        ##
        ## If no possible splitpoint can be found, then store no possible split found
        ## (which is setting best.split as NA)
        ##
        best.split <- NULL
        
        
      }## end if-else for evaluating the list of possible splitpoints
      
    }## end if-else for checking if current node was already terminal
    
    ##
    ## Return the best found splitpoint and look for follow-up splits
    ##
    
    if(is.null(best.split)){
      
      ##
      ## Make current node terminal
      ##
      tree <- rbind(tree, 
                    data.frame(ID = ID,
                               Feature = NA,
                               Value = NA,
                               numeric = NA,
                               terminal = TRUE,
                               Prediction = I(list(as.vector(colMeans(Y), mode = "numeric")))
                    )
      )
      
    }else{
      
      ##
      ## Store new split in tree
      ## (in our study every feature is numeric, so we set numeric always TRUE;
      ## if in follow up studies you have categorical features, please check before
      ## and add this information right here)
      ##
      tree <- rbind(tree, 
                    data.frame(ID = ID,
                               Feature = names(X)[best.split$feature],
                               Value = best.split$Values,
                               numeric = TRUE,
                               terminal = FALSE,
                               Prediction = I(list(as.vector(colMeans(Y), mode = "numeric")))
                    )
      )
      
      ##
      ## Look for further splits in the left childnode
      ##
      tree <- node.split(ID = ID*2, 
                         tree = tree, 
                         Y = data.frame(Y[X[,best.split$feature] <= best.split$Values,]), 
                         X = data.frame(X[X[,best.split$feature] <= best.split$Values,]), 
                         min.node.size = min.node.size, 
                         mtry = mtry, 
                         num.random.cuts = num.random.cuts,
                         cp = cp,
                         split=split)
      
      ##
      ## Look for further splits in the right childnode
      ##
      tree <- node.split(ID = ID*2+1, 
                         tree = tree, 
                         Y = data.frame(Y[X[,best.split$feature] > best.split$Values,]), 
                         X = data.frame(X[X[,best.split$feature] > best.split$Values,]), 
                         min.node.size = min.node.size, 
                         mtry = mtry, 
                         num.random.cuts = num.random.cuts,
                         cp = cp,
                         split=split)
      
    }
    
    ##
    ## Return the tree
    ##
    tree
    
  }
  
  predict_forest <- function(X, 
                             object){
    
    ##
    ## Apply the predict_tree function on each of the trees and average results
    ## for each observation
    ##
    
    predictions.forest <- lapply(1:length(object), 
                                 FUN = function(curr.tree){
                                   
                                   ##
                                   ## Predict all observations for current tree
                                   ##
                                   
                                   predict_tree(X, object = object[[curr.tree]])
                                   
                                 })
    
    predictions.forest
    
  }
  
  predict_tree <- function(X, 
                           object){
    
    ##
    ## Apply for each observation a search for the terminal node
    ##
    predictions <- lapply(1:nrow(X), FUN = function(observation){
      
      ##
      ## Initialize the "while"-argument and the ID of the current node
      ##
      terminal.pred <- TRUE
      curr.id <- 1
      
      while(terminal.pred){
        
        ##
        ## For each node look if it is terminal, if not look for the ID of the next node
        ##
        if(object[which(object$ID == curr.id), "terminal"]){
          
          ##
          ## If current node is terminal, then stop the while-loop and take the 
          ## "Prediction"-argument as the prediction of current observation
          terminal.pred <- FALSE
          pred <- object[which(object$ID == curr.id), "Prediction"]
          
        }else{
          
          ##
          ## Get the split feature and split value of current node
          ##
          curr.split.feature <- object[which(object$ID == curr.id), "Feature"]
          curr.split.value <- object[which(object$ID == curr.id), "Value"]
          
          ##
          ## If feature value of current observation is smaller of equal the split
          ## value then check out the node with the current ID * 2 and if not 
          ## check out the node with the current ID * 2 + 1
          ##
          if(X[observation, as.character(curr.split.feature)] <= curr.split.value){
            
            curr.id <- curr.id*2
            
          }else{
            
            curr.id <- curr.id*2+1
            
          }## close if-else from checking if Value of observation is bigger oder smaller-equal the split value
          
        }## close if-else from checking if current node is terminal or observation has to go to next node
        
      }## closing while from looking for the terminal node of the observation
      
      ##
      ## After finding the respective terminal node, give out the prediction(s)
      ##
      pred
      
    })## close lapply from going over complete dataset X
    
    ##
    ## Return the predictions
    ##
    matrix(unlist(predictions), byrow = TRUE, ncol = length(predictions[[1]][[1]]))
    
  }
  
  treefit <- function(Y,                            ## -> 1 or multiple targetvectors (dim = n X #targets)
                      X,                            ## -> 1 or multiple features (dim = n X #features)
                      min.node.size = 5,            ## -> minimal node size (5 for regression is standard; 
                      ## for other tasks you should alter the default)
                      mtry = NA,                    ## -> number of features used for determining split
                      num.random.cuts = NA,         ## -> number of random cuts per sampled feature
                      cp = 0.01,                     ## -> complexity parameter for minimal increase in 
                      ## splitting criterion to attempt a split
                      split="L2"
  ){   
    
    
    ##
    ## This code is especially for regression. For other types of trees it is 
    ## mentioned where changes have to be done
    ##
    
    ##
    ## To avoid producing different trees by just ordering X in a different way
    ## first the columns of X are sorted
    ##
    X <- X[, order(names(X))]
    
    ##
    ## Initialize a data.frame for the final tree object
    ##
    final_tree_object <- data.frame()
    final_tree_object <- node.split(ID = 1, 
                                    tree = final_tree_object, 
                                    Y = Y, 
                                    X = X, 
                                    min.node.size = min.node.size,
                                    mtry = mtry,
                                    num.random.cuts = num.random.cuts,
                                    cp = cp,
                                    split=split)
    
    ##
    ## Return final tree object
    ##
    
    final_tree_object
    
  }
  
}


#### First Experiment ####
##
## Number of observations drawn for the experiment
##
n <- 100

##
## Number of covariables
##
p <- 10


k<-5

##
## Number of repetitions for the experiment
##
rep <- 1000

##

set.seed(1234)

for(i in 1:rep){
  ##
  ## Draw random variables
  ##
  rv_T <- rt(n = n,
             df = 2)
  rv_C2 <- runif(n = n,
                 min = 0,
                 max = 2)
  rv_C4 <- runif(n = n,
                 min = 0,
                 max = 4)
  rv_C8 <- runif(n = n,
                 min = 0,
                 max = 8)
  rv_W <- rexp(n = n,
               rate = 1)
  rv_Z <- rnorm(n = n,
                mean = 0,
                sd = 1)
  rv_U <- runif(n = n,
                min = 0,
                max = 1)
  
  ##
  ## Initialize data.frames
  ## (ind = independent, wd = weakly dependent, sd = strongly dependent)
  ##
  ind <- as.data.frame(matrix(data = 0, nrow = n, ncol = p))
  wd <- as.data.frame(matrix(data = 0, nrow = n, ncol = p))
  sd <- as.data.frame(matrix(data = 0, nrow = n, ncol = p))
  
  names(ind) <- paste0("X", 1:p)
  names(wd) <- paste0("X", 1:p)
  names(sd) <- paste0("X", 1:p)
  
  ##
  ## Fill first 5 covariables
  ##
  ind$X1 <- rv_Z
  ind$X2 <- rv_W
  ind$X3 <- rv_T
  ind$X4 <- rv_C4
  ind$X5 <- rv_C8
  
  wd$X1 <- rv_Z + rv_W + rv_T
  wd$X2 <- rv_W
  wd$X3 <- rv_T
  wd$X4 <- rv_C2
  wd$X5 <- rv_C8
  
  sd$X1 <- 0.1*rv_Z + rv_W
  sd$X2 <- rv_W
  sd$X3 <- rv_T
  sd$X4 <- rv_C2
  sd$X5 <- rv_C8
  
  ##
  ## Fill Rest
  ##
  if(p > 5){
    ind[, 6:p] <- runif(n*(p-5),
                        min = 0,
                        max = 1)
    
    wd[, 6:p] <- ind[, 6:p]
    sd[, 6:p] <- ind[, 6:p]
  }
  
  

  ##
  ## Generate Ys for independent data
  ##
  
  Y_ind_0_1 <- cbind(0.1 * sin(ind[, 1]), 0.5*log(ind[,2]), rv_U)+ rmvnorm(n = n,
                                          mean = c(0,0,0),
                                          sigma = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3))
  
  Y_ind_05_1 <-cbind(0.1 * sin(ind[, 1]), 0.5*log(ind[,2]), rv_U) +  rmvnorm(n = n,
                                           mean = c(0,0,0),
                                           sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), nrow = 3, ncol = 3))
  
  Y_ind_05_2 <- cbind(0.1 * sin(ind[, 1]), 0.5*log(ind[,2]), rv_U) + rmvnorm(n = n,
                                           mean = c(0,0,0),
                                           sigma = matrix(c(1,0.5,0.25,0.5,1,0.5,0.25,0.5,1), nrow = 3, ncol = 3))
  
  Y_ind_09_1 <- cbind(0.1 * sin(ind[, 1]), 0.5*log(ind[,2]), rv_U) + rmvnorm(n = n,
                                           mean = c(0,0,0),
                                           sigma = matrix(c(1,0.9,0.9,0.9,1,0.9,0.9,0.9,1), nrow = 3, ncol = 3))
  
  Y_ind_09_2 <- cbind(0.1 * sin(ind[, 1]), 0.5*log(ind[,2]), rv_U)+ rmvnorm(n = n,
                                           mean = c(0,0,0),
                                           sigma = matrix(c(1,0.9,0.81,0.9,1,0.9,0.81,0.9,1), nrow = 3, ncol = 3))
  
  ##
  ## Generate Ys for weakly dependent data
  ##
  Y_wd_0_1 <-  cbind(0.1 * sin(wd[, 1]), 0.5*log(wd[,2]), rv_U) + rmvnorm(n = n,
                                        mean = c(0,0,0),
                                        sigma = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3))
  
  Y_wd_05_1 <- cbind(0.1 * sin(wd[, 1]), 0.5*log(wd[,2]), rv_U) + rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), nrow = 3, ncol = 3))
  
  Y_wd_05_2 <-  cbind(0.1 * sin(wd[, 1]), 0.5*log(wd[,2]), rv_U) + rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.5,0.25,0.5,1,0.5,0.25,0.5,1), nrow = 3, ncol = 3))
  
  Y_wd_09_1 <-  cbind(0.1 * sin(wd[, 1]), 0.5*log(wd[,2]), rv_U) + rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.9,0.9,0.9,1,0.9,0.9,0.9,1), nrow = 3, ncol = 3))
  
  Y_wd_09_2 <-  cbind(0.1 * sin(wd[, 1]), 0.5*log(wd[,2]), rv_U) + rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.9,0.81,0.9,1,0.9,0.81,0.9,1), nrow = 3, ncol = 3))
  
  ##
  ## Generate Ys for strongly dependent data
  ##
  Y_sd_0_1 <- cbind(0.1 * sin(sd[, 1]), 0.5*log(sd[,2]), rv_U) + rmvnorm(n = n,
                                        mean = c(0,0,0),
                                        sigma = matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3))
  
  Y_sd_05_1 <-cbind(0.1 * sin(sd[, 1]), 0.5*log(sd[,2]), rv_U) + rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), nrow = 3, ncol = 3))
  
  Y_sd_05_2 <- cbind(0.1 * sin(sd[, 1]), 0.5*log(sd[,2]), rv_U) + rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.5,0.25,0.5,1,0.5,0.25,0.5,1), nrow = 3, ncol = 3))
  
  Y_sd_09_1 <- cbind(0.1 * sin(sd[, 1]), 0.5*log(sd[,2]), rv_U)+ rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.9,0.9,0.9,1,0.9,0.9,0.9,1), nrow = 3, ncol = 3))
  
  Y_sd_09_2 <- cbind(0.1 * sin(sd[, 1]), 0.5*log(sd[,2]), rv_U)+ rmvnorm(n = n,
                                         mean = c(0,0,0),
                                         sigma = matrix(c(1,0.9,0.81,0.9,1,0.9,0.81,0.9,1), nrow = 3, ncol = 3))
  
  Y_ind <- list(Y_ind_0_1, Y_ind_05_1, Y_ind_05_2, Y_ind_09_1, Y_ind_09_2)
  Y_wd <- list(Y_wd_0_1, Y_wd_05_1, Y_wd_05_2, Y_wd_09_1, Y_wd_09_2)
  Y_sd <- list(Y_sd_0_1, Y_sd_05_1, Y_sd_05_2, Y_sd_09_1, Y_sd_09_2)
  ##
  ## Initialize final data.frame with MSE
  ##
  mse_df <- data.frame()
  
  ##
  ## split data in training and prediction 
  ##
  
  
  indizes_out <- sample((1:(n)), size = n/k, replace = FALSE)
  
  ##
  ## rf Models
  ##
  {
    ind_rf_models_L1 <- lapply(X = Y_ind,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = ind[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           split= "L1"
                                 )
                                 
                               })
    
    wd_rf_models_L1 <- lapply(X = Y_wd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = wd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          split= "L1"
                                )
                                
                              })
    
    sd_rf_models_L1 <- lapply(X = Y_sd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = sd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          split= "L1"
                                )
                                
                              })
    ind_rf_models_L1m <- lapply(X = Y_ind,
                                FUN = function(X){
                                  
                                  forestfit(Y = X[-indizes_out,],
                                            X = ind[-indizes_out,],
                                            ntree = 500,
                                            mtry = max(
                                              floor(p/3),
                                              1
                                            ),
                                            min.node.size = 5,
                                            split= "L1m"
                                  )
                                  
                                })
    
    wd_rf_models_L1m <- lapply(X = Y_wd,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = wd[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           split= "L1m"
                                 )
                                 
                               })
    
    sd_rf_models_L1m <- lapply(X = Y_sd,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = sd[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           split= "L1m"
                                 )
                                 
                               })
    ind_rf_models_L2 <- lapply(X = Y_ind,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = ind[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           split= "L2"
                                 )
                                 
                               })
    
    wd_rf_models_L2 <- lapply(X = Y_wd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = wd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          split= "L2"
                                )
                                
                              })
    
    sd_rf_models_L2 <- lapply(X = Y_sd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = sd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          split= "L2"
                                )
                                
                              })
    
  }
  ##
  ## extraTree Models
  ##
  {
    ind_et_models_L1 <- lapply(X = Y_ind,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = ind[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           num.random.cuts = 1,
                                           bos=FALSE,
                                           split= "L1"
                                 )
                                 
                               })
    
    wd_et_models_L1 <- lapply(X = Y_wd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = wd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          num.random.cuts = 1,
                                          bos=FALSE,
                                          split= "L1"
                                )
                                
                              })
    
    sd_et_models_L1 <- lapply(X = Y_sd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = sd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          num.random.cuts = 1,
                                          bos=FALSE,
                                          split= "L1"
                                )
                                
                              })
    ind_et_models_L1m <- lapply(X = Y_ind,
                                FUN = function(X){
                                  
                                  forestfit(Y = X[-indizes_out,],
                                            X = ind[-indizes_out,],
                                            ntree = 500,
                                            mtry = max(
                                              floor(p/3),
                                              1
                                            ),
                                            min.node.size = 5,
                                            num.random.cuts = 1,
                                            bos=FALSE,
                                            split= "L1m"
                                  )
                                  
                                })
    
    wd_et_models_L1m <- lapply(X = Y_wd,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = wd[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           num.random.cuts = 1,
                                           bos=FALSE,
                                           split= "L1m"
                                 )
                                 
                               })
    
    sd_et_models_L1m <- lapply(X = Y_sd,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = sd[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           num.random.cuts = 1,
                                           bos=FALSE,
                                           split= "L1m"
                                 )
                                 
                               })
    ind_et_models_L2 <- lapply(X = Y_ind,
                               FUN = function(X){
                                 
                                 forestfit(Y = X[-indizes_out,],
                                           X = ind[-indizes_out,],
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           min.node.size = 5,
                                           num.random.cuts = 1,
                                           bos=FALSE,
                                           split= "L2"
                                 )
                                 
                               })
    
    wd_et_models_L2 <- lapply(X = Y_wd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = wd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          num.random.cuts = 1,
                                          bos=FALSE,
                                          split= "L2"
                                )
                                
                              })
    
    sd_et_models_L2 <- lapply(X = Y_sd,
                              FUN = function(X){
                                
                                forestfit(Y = X[-indizes_out,],
                                          X = sd[-indizes_out,],
                                          ntree = 500,
                                          mtry = max(
                                            floor(p/3),
                                            1
                                          ),
                                          min.node.size = 5,
                                          num.random.cuts = 1,
                                          bos=FALSE,
                                          split= "L2"
                                )
                                
                              })
    
  }
  ##
  ## univ. rf Models
  ##
  {
    ind_univ_rf_models_L1 <- lapply(X = Y_ind,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1"
                                        )
                                      )
                                    })
    
    wd_univ_rf_models_L1 <- lapply(X = Y_wd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L1"
                                       )
                                     )
                                   })
    
    sd_univ_rf_models_L1 <- lapply(X = Y_sd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L1"
                                       )
                                     )
                                   })
    ind_univ_rf_models_L1m <- lapply(X = Y_ind,
                                     FUN = function(X){
                                       
                                       list(
                                         forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                   X = ind[-indizes_out,],
                                                   ntree = 500,
                                                   mtry = max(
                                                     floor(p/3),
                                                     1
                                                   ),
                                                   min.node.size = 5,
                                                   split="L1m"
                                         ),
                                         forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                   X = ind[-indizes_out,],
                                                   ntree = 500,
                                                   mtry = max(
                                                     floor(p/3),
                                                     1
                                                   ),
                                                   min.node.size = 5,
                                                   split="L1m"
                                         ),
                                         forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                   X = ind[-indizes_out,],
                                                   ntree = 500,
                                                   mtry = max(
                                                     floor(p/3),
                                                     1
                                                   ),
                                                   min.node.size = 5,
                                                   split="L1m"
                                         )
                                       )
                                     })
    
    wd_univ_rf_models_L1m <- lapply(X = Y_wd,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = wd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = wd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = wd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1m"
                                        )
                                      )
                                    })
    
    sd_univ_rf_models_L1m <- lapply(X = Y_sd,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = sd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = sd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = sd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L1m"
                                        )
                                      )
                                    })
    
    ind_univ_rf_models_L2 <- lapply(X = Y_ind,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L2"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L2"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  split="L2"
                                        )
                                      )
                                    })
    
    wd_univ_rf_models_L2 <- lapply(X = Y_wd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L2"
                                       )
                                     )
                                   })
    
    sd_univ_rf_models_L2 <- lapply(X = Y_sd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 split="L2"
                                       )
                                     )
                                   })
    
    
  }
  ##
  ## univ. extraTree Models
  ##
  {
    ind_univ_et_models_L1 <- lapply(X = Y_ind,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1"
                                        )
                                      )
                                    })
    
    wd_univ_et_models_L1 <- lapply(X = Y_wd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L1"
                                       )
                                     )
                                   })
    
    sd_univ_et_models_L1 <- lapply(X = Y_sd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L1"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L1"
                                       )
                                     )
                                   })
    
    
    ind_univ_et_models_L1m <- lapply(X = Y_ind,
                                     FUN = function(X){
                                       
                                       list(
                                         forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                   X = ind[-indizes_out,],
                                                   ntree = 500,
                                                   mtry = max(
                                                     floor(p/3),
                                                     1
                                                   ),
                                                   min.node.size = 5,
                                                   num.random.cuts = 1,
                                                   bos=FALSE,
                                                   split="L1m"
                                         ),
                                         forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                   X = ind[-indizes_out,],
                                                   ntree = 500,
                                                   mtry = max(
                                                     floor(p/3),
                                                     1
                                                   ),
                                                   min.node.size = 5,
                                                   num.random.cuts = 1,
                                                   bos=FALSE,
                                                   split="L1m"
                                         ),
                                         forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                   X = ind[-indizes_out,],
                                                   ntree = 500,
                                                   mtry = max(
                                                     floor(p/3),
                                                     1
                                                   ),
                                                   min.node.size = 5,
                                                   num.random.cuts = 1,
                                                   bos=FALSE,
                                                   split="L1m"
                                         )
                                       )
                                     })
    
    wd_univ_et_models_L1m <- lapply(X = Y_wd,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = wd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = wd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = wd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1m"
                                        )
                                      )
                                    })
    
    sd_univ_et_models_L1m <- lapply(X = Y_sd,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = sd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = sd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1m"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = sd[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L1m"
                                        )
                                      )
                                    })
    ind_univ_et_models_L2 <- lapply(X = Y_ind,
                                    FUN = function(X){
                                      
                                      list(
                                        forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L2"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L2"
                                        ),
                                        forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                  X = ind[-indizes_out,],
                                                  ntree = 500,
                                                  mtry = max(
                                                    floor(p/3),
                                                    1
                                                  ),
                                                  min.node.size = 5,
                                                  num.random.cuts = 1,
                                                  bos=FALSE,
                                                  split="L2"
                                        )
                                      )
                                    })
    
    wd_univ_et_models_L2 <- lapply(X = Y_wd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = wd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L2"
                                       )
                                     )
                                   })
    
    sd_univ_et_models_L2 <- lapply(X = Y_sd,
                                   FUN = function(X){
                                     
                                     list(
                                       forestfit(Y = data.frame(X[-indizes_out, 1]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 2]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L2"
                                       ),
                                       forestfit(Y = data.frame(X[-indizes_out, 3]),
                                                 X = sd[-indizes_out,],
                                                 ntree = 500,
                                                 mtry = max(
                                                   floor(p/3),
                                                   1
                                                 ),
                                                 min.node.size = 5,
                                                 num.random.cuts = 1,
                                                 bos=FALSE,
                                                 split="L2"
                                       )
                                     )
                                   })
    
    
  }
  
  ##
  ## extraTrees multi-task Models
  ##
  {
    ind_et_mt_models <- lapply(X = Y_ind,
                               FUN = function(X){
                                 
                                 extraTrees(x = rbind(ind[-indizes_out,],
                                                      ind[-indizes_out,],
                                                      ind[-indizes_out,]),
                                            y = c(X[-indizes_out,1],
                                                  X[-indizes_out,2],
                                                  X[-indizes_out,3]),
                                            ntree = 500,
                                            mtry = max(
                                              floor(p/3),
                                              1
                                            ),
                                            nodesize = 5,
                                            numRandomCuts = 1,
                                            tasks = rep(1:3,
                                                        each = nrow(ind[-indizes_out,]))
                                 )
                                 
                               })
    
    wd_et_mt_models <- lapply(X = Y_wd,
                              FUN = function(X){
                                
                                extraTrees(x = rbind(wd[-indizes_out,],
                                                     wd[-indizes_out,],
                                                     wd[-indizes_out,]),
                                           y = c(X[-indizes_out,1],
                                                 X[-indizes_out,2],
                                                 X[-indizes_out,3]),
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           nodesize = 5,
                                           numRandomCuts = 1,
                                           tasks = rep(1:3,
                                                       each = nrow(ind[-indizes_out,]))
                                )
                                
                              })
    
    sd_et_mt_models <- lapply(X = Y_sd,
                              FUN = function(X){
                                
                                extraTrees(x = rbind(sd[-indizes_out,],
                                                     sd[-indizes_out,],
                                                     sd[-indizes_out,]),
                                           y = c(X[-indizes_out,1],
                                                 X[-indizes_out,2],
                                                 X[-indizes_out,3]),
                                           ntree = 500,
                                           mtry = max(
                                             floor(p/3),
                                             1
                                           ),
                                           nodesize = 5,
                                           numRandomCuts = 1,
                                           tasks = rep(1:3,
                                                       each = nrow(ind[-indizes_out,]))
                                )
                                
                              })
    
  }
  
  ##
  ## Calculate goodness of fit #
  ##
  
  ##
  ## rf Models
  
  ##
  {
    mse_ind_rf_models_L1<- lapply(X = 1:length(ind_rf_models_L1),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(ind[indizes_out,], ind_rf_models_L1[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_ind[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
    
    
    mse_wd_rf_models_L1 <- lapply(X = 1:length(wd_rf_models_L1),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(wd[indizes_out,], wd_rf_models_L1[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_wd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
    
    mse_sd_rf_models_L1 <- lapply(X = 1:length(sd_rf_models_L1),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(sd[indizes_out,], sd_rf_models_L1[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_sd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    mse_ind_rf_models_L1m<- lapply(X = 1:length(ind_rf_models_L1m),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(ind[indizes_out,], ind_rf_models_L1m[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_ind[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    
    
    
    mse_wd_rf_models_L1m <- lapply(X = 1:length(wd_rf_models_L1m),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(wd[indizes_out,], wd_rf_models_L1m[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_wd[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    
    mse_sd_rf_models_L1m <- lapply(X = 1:length(sd_rf_models_L1m),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(sd[indizes_out,], sd_rf_models_L1m[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_sd[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    mse_ind_rf_models_L2<- lapply(X = 1:length(ind_rf_models_L2),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(ind[indizes_out,], ind_rf_models_L2[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_ind[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
    
    
    mse_wd_rf_models_L2 <- lapply(X = 1:length(wd_rf_models_L2),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(wd[indizes_out,], wd_rf_models_L2[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_wd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
    mse_sd_rf_models_L2 <- lapply(X = 1:length(sd_rf_models_L2),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(sd[indizes_out,], sd_rf_models_L2[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_sd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
  }
  
  ##
  ## extraTree Models
  ##
  
  {
    mse_ind_et_models_L1 <- lapply(X = 1:length(ind_et_models_L1),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(ind[indizes_out,], ind_et_models_L1[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_ind[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    
    
    
    mse_wd_et_models_L1 <- lapply(X = 1:length(wd_et_models_L1),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(wd[indizes_out,], wd_et_models_L1[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_wd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
    mse_sd_et_models_L1 <- lapply(X = 1:length(sd_et_models_L1),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(sd[indizes_out,], sd_et_models_L1[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_sd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    mse_ind_et_models_L1m <- lapply(X = 1:length(ind_et_models_L1m),
                                    FUN = function(X){
                                      
                                      preds_all <- predict_forest(ind[indizes_out,], ind_et_models_L1m[[X]])
                                      preds <- preds_all[[1]]
                                      for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                      preds <- preds/length(preds_all)
                                      mean(as.vector((preds - Y_ind[[X]][indizes_out,]))^2)
                                      
                                      
                                    })
    
    
    mse_wd_et_models_L1m <- lapply(X = 1:length(wd_et_models_L1m),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(wd[indizes_out,], wd_et_models_L1m[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_wd[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    
    mse_sd_et_models_L1m <- lapply(X = 1:length(sd_et_models_L1m),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(sd[indizes_out,], sd_et_models_L1m[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_sd[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    mse_ind_et_models_L2 <- lapply(X = 1:length(ind_et_models_L2),
                                   FUN = function(X){
                                     
                                     preds_all <- predict_forest(ind[indizes_out,], ind_et_models_L2[[X]])
                                     preds <- preds_all[[1]]
                                     for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                     preds <- preds/length(preds_all)
                                     mean(as.vector((preds - Y_ind[[X]][indizes_out,]))^2)
                                     
                                     
                                   })
    
    
    
    mse_wd_et_models_L2 <- lapply(X = 1:length(wd_et_models_L2),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(wd[indizes_out,], wd_et_models_L2[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_wd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
    mse_sd_et_models_L2 <- lapply(X = 1:length(sd_et_models_L2),
                                  FUN = function(X){
                                    
                                    preds_all <- predict_forest(sd[indizes_out,], sd_et_models_L2[[X]])
                                    preds <- preds_all[[1]]
                                    for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                    preds <- preds/length(preds_all)
                                    mean(as.vector((preds - Y_sd[[X]][indizes_out,]))^2)
                                    
                                    
                                  })
    
  }
  
  ##
  ## univ. rf Models
  ##
  {
    mse_ind_univ_rf_models_L1 <- lapply(X = 1:length(ind_univ_rf_models_L1),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L1[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L1[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L1[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                        })
    
    
    
    mse_wd_univ_rf_models_L1 <- lapply(X = 1:length(wd_univ_rf_models_L1),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L1[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L1[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L1[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
    mse_sd_univ_rf_models_L1 <- lapply(X = 1:length(sd_univ_rf_models_L1),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L1[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L1[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L1[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
    mse_ind_univ_rf_models_L1m <- lapply(X = 1:length(ind_univ_rf_models_L1m),
                                         FUN = function(X){
                                           
                                           ##
                                           ## Predictions for forest for first Y[[X]]
                                           ##
                                           preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L1m[[X]][[1]])
                                           preds <- preds_all[[1]]
                                           for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                           preds <- preds/length(preds_all)
                                           mse1 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,1]))^2)
                                           
                                           ##
                                           ## Predictions for forest for second Y[[X]]
                                           ##
                                           preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L1m[[X]][[2]])
                                           preds <- preds_all[[1]]
                                           for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                           preds <- preds/length(preds_all)
                                           mse2 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,2]))^2)
                                           
                                           ##
                                           ## Predictions for forest for third Y[[X]]
                                           ##
                                           preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L1m[[X]][[3]])
                                           preds <- preds_all[[1]]
                                           for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                           preds <- preds/length(preds_all)
                                           mse3 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,3]))^2)
                                           
                                           mean(mse1, mse2, mse3)
                                           
                                         })
    
    
    
    mse_wd_univ_rf_models_L1m <- lapply(X = 1:length(wd_univ_rf_models_L1m),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L1m[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L1m[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L1m[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                          
                                        })
    
    
    mse_sd_univ_rf_models_L1m <- lapply(X = 1:length(sd_univ_rf_models_L1m),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L1m[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L1m[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L1m[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                          
                                        })
    
    
    mse_ind_univ_rf_models_L2 <- lapply(X = 1:length(ind_univ_rf_models_L2),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L2[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L2[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_rf_models_L2[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                        })
    
    
    
    mse_wd_univ_rf_models_L2 <- lapply(X = 1:length(wd_univ_rf_models_L2),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L2[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L2[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_rf_models_L2[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
    
    mse_sd_univ_rf_models_L2 <- lapply(X = 1:length(sd_univ_rf_models_L2),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L2[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L2[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_rf_models_L2[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
  }
  ##
  ## univ. et Models
  ##
  {
    mse_ind_univ_et_models_L1 <- lapply(X = 1:length(ind_univ_et_models_L1),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L1[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L1[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L1[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                        })
    
    
    
    mse_wd_univ_et_models_L1 <- lapply(X = 1:length(wd_univ_et_models_L1),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L1[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L1[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L1[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
    mse_sd_univ_et_models_L1 <- lapply(X = 1:length(sd_univ_et_models_L1),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L1[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L1[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L1[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
    mse_ind_univ_et_models_L1m <- lapply(X = 1:length(ind_univ_et_models_L1m),
                                         FUN = function(X){
                                           
                                           ##
                                           ## Predictions for forest for first Y[[X]]
                                           ##
                                           preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L1m[[X]][[1]])
                                           preds <- preds_all[[1]]
                                           for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                           preds <- preds/length(preds_all)
                                           mse1 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,1]))^2)
                                           
                                           ##
                                           ## Predictions for forest for second Y[[X]]
                                           ##
                                           preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L1m[[X]][[2]])
                                           preds <- preds_all[[1]]
                                           for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                           preds <- preds/length(preds_all)
                                           mse2 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,2]))^2)
                                           
                                           ##
                                           ## Predictions for forest for third Y[[X]]
                                           ##
                                           preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L1m[[X]][[3]])
                                           preds <- preds_all[[1]]
                                           for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                           preds <- preds/length(preds_all)
                                           mse3 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,3]))^2)
                                           
                                           mean(mse1, mse2, mse3)
                                           
                                         })
    
    
    mse_wd_univ_et_models_L1m <- lapply(X = 1:length(wd_univ_et_models_L1m),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L1m[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L1m[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L1m[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                          
                                        })
    
    mse_sd_univ_et_models_L1m <- lapply(X = 1:length(sd_univ_et_models_L1m),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L1m[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L1m[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L1m[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                          
                                        })
    
    mse_ind_univ_et_models_L2 <- lapply(X = 1:length(ind_univ_et_models_L2),
                                        FUN = function(X){
                                          
                                          ##
                                          ## Predictions for forest for first Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L2[[X]][[1]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse1 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,1]))^2)
                                          
                                          ##
                                          ## Predictions for forest for second Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L2[[X]][[2]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse2 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,2]))^2)
                                          
                                          ##
                                          ## Predictions for forest for third Y[[X]]
                                          ##
                                          preds_all <- predict_forest(ind[indizes_out,], ind_univ_et_models_L2[[X]][[3]])
                                          preds <- preds_all[[1]]
                                          for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                          preds <- preds/length(preds_all)
                                          mse3 <- mean(as.vector((preds - Y_ind[[X]][indizes_out,3]))^2)
                                          
                                          mean(mse1, mse2, mse3)
                                          
                                        })
    
    
    
    mse_wd_univ_et_models_L2 <- lapply(X = 1:length(wd_univ_et_models_L2),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L2[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L2[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(wd[indizes_out,], wd_univ_et_models_L2[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_wd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
    
    mse_sd_univ_et_models_L2 <- lapply(X = 1:length(sd_univ_et_models_L2),
                                       FUN = function(X){
                                         
                                         ##
                                         ## Predictions for forest for first Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L2[[X]][[1]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse1 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,1]))^2)
                                         
                                         ##
                                         ## Predictions for forest for second Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L2[[X]][[2]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse2 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,2]))^2)
                                         
                                         ##
                                         ## Predictions for forest for third Y[[X]]
                                         ##
                                         preds_all <- predict_forest(sd[indizes_out,], sd_univ_et_models_L2[[X]][[3]])
                                         preds <- preds_all[[1]]
                                         for(tree in 2:length(preds_all)){preds <- preds+preds_all[[tree]]}
                                         preds <- preds/length(preds_all)
                                         mse3 <- mean(as.vector((preds - Y_sd[[X]][indizes_out,3]))^2)
                                         
                                         mean(mse1, mse2, mse3)
                                         
                                         
                                       })
  }
  
  ##
  ## extraTrees mt Models
  ##
  {
    mse_ind_et_mt_models <- lapply(X = 1:length(ind_et_mt_models),
                                   FUN = function(X){
                                     
                                     preds_all <- predict(ind_et_mt_models[[1]], 
                                                          newdata = rbind(ind[indizes_out,],
                                                                          ind[indizes_out,],
                                                                          ind[indizes_out,]),
                                                          newtasks = rep(1:3, each = length(indizes_out)))
                                     mean((preds_all - c(Y_ind[[X]][indizes_out,1],
                                                         Y_ind[[X]][indizes_out,2],
                                                         Y_ind[[X]][indizes_out,3]))^2)
                                     
                                     
                                   })
    
    
    
    mse_wd_et_mt_models <- lapply(X = 1:length(wd_et_mt_models),
                                  FUN = function(X){
                                    
                                    preds_all <- predict(wd_et_mt_models[[1]], 
                                                         newdata = rbind(wd[indizes_out,],
                                                                         wd[indizes_out,],
                                                                         wd[indizes_out,]),
                                                         newtasks = rep(1:3, each = length(indizes_out)))
                                    mean((preds_all - c(Y_wd[[X]][indizes_out,1],
                                                        Y_wd[[X]][indizes_out,2],
                                                        Y_wd[[X]][indizes_out,3]))^2)
                                    
                                  })
    
    mse_sd_et_mt_models <- lapply(X = 1:length(sd_et_mt_models),
                                  FUN = function(X){
                                    
                                    preds_all <- predict(ind_et_mt_models[[1]], 
                                                         newdata = rbind(sd[indizes_out,],
                                                                         sd[indizes_out,],
                                                                         sd[indizes_out,]),
                                                         newtasks = rep(1:3, each = length(indizes_out)))
                                    mean((preds_all - c(Y_sd[[X]][indizes_out,1],
                                                        Y_sd[[X]][indizes_out,2],
                                                        Y_sd[[X]][indizes_out,3]))^2)
                                    
                                  })
  }
  
  mse_df <- rbind(mse_df,
                  data.frame(
                    rf_L1 = c(unlist(mse_ind_rf_models_L1),
                              unlist(mse_wd_rf_models_L1),
                              unlist(mse_sd_rf_models_L1)),
                    rf_L1m = c(unlist(mse_ind_rf_models_L1m),
                               unlist(mse_wd_rf_models_L1m),
                               unlist(mse_sd_rf_models_L1m)),
                    rf_L2 = c(unlist(mse_ind_rf_models_L2),
                              unlist(mse_wd_rf_models_L2),
                              unlist(mse_sd_rf_models_L2)),
                    et_L1 = c(unlist(mse_ind_et_models_L1),
                              unlist(mse_wd_et_models_L1),
                              unlist(mse_sd_et_models_L1)),
                    et_L1m = c(unlist(mse_ind_et_models_L1m),
                               unlist(mse_wd_et_models_L1m),
                               unlist(mse_sd_et_models_L1m)),
                    et_L2 = c(unlist(mse_ind_et_models_L2),
                              unlist(mse_wd_et_models_L2),
                              unlist(mse_sd_et_models_L2)),
                    univ_rf_L1 = c(unlist(mse_ind_univ_rf_models_L1),
                                   unlist(mse_wd_univ_rf_models_L1),
                                   unlist(mse_sd_univ_rf_models_L1)),
                    univ_rf_L1m = c(unlist(mse_ind_univ_rf_models_L1m),
                                    unlist(mse_wd_univ_rf_models_L1m),
                                    unlist(mse_sd_univ_rf_models_L1m)),
                    univ_rf_L2 = c(unlist(mse_ind_univ_rf_models_L2),
                                   unlist(mse_wd_univ_rf_models_L2),
                                   unlist(mse_sd_univ_rf_models_L2)),
                    univ_et_L1 = c(unlist(mse_ind_univ_et_models_L1),
                                   unlist(mse_wd_univ_et_models_L1),
                                   unlist(mse_sd_univ_et_models_L1)),
                    univ_et_L1m = c(unlist(mse_ind_univ_et_models_L1m),
                                    unlist(mse_wd_univ_et_models_L1m),
                                    unlist(mse_sd_univ_et_models_L1m)),
                    univ_et_L2 = c(unlist(mse_ind_univ_et_models_L2),
                                   unlist(mse_wd_univ_et_models_L2),
                                   unlist(mse_sd_univ_et_models_L2)),
                    et_mt = c(unlist(mse_ind_et_mt_models),
                              unlist(mse_wd_et_mt_models),
                              unlist(mse_sd_et_mt_models)),
                    Y = rep(c("Y_0_1", "Y_05_1", "Y_05_2", "Y_09_1", "Y_09_2"),
                            times = 3),
                    data = rep(c("ind", "wd", "sd"), each = 5)
                  )
  )
  
  
  
  mse_list <- c(mse_list, list(mse_df))
}

saveRDS(mse_list, file = paste0("mgam2_", n, "_mse.rds"))



