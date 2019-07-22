pkgs <- c("rpart")

sapply(pkgs, require, character.only = TRUE)

data <- read.csv(url("http://bit.ly/flex_iris"),TRUE)

independent   <- 1:4
dependent     <- 7


# divide the dataset into input features and output feature
X             <- data[, independent]
y             <- data[, dependent]

####################################################################################################
# predict function
####################################################################################################
predict.FlexBoost <- function(object, X, type = c("response", "prob"), n_tree = NULL){
  
  # handle args
  type <- match.arg(type)
  
  if(is.null(n_tree)) { tree_seq <- seq_along(object$alphas) } 
  
  else                { tree_seq <- seq(1, n_tree) }
  
  # evaluate score function on sample
  f <- 0
  
  for(i in tree_seq){
    tree       <- object$trees[[i]]
    tree$terms <- object$terms
    pred       <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
    f          <- f + object$alphas[i] * pred
  }
  
  # handle response type
  if(type == "response")  { sign(f) } 
  
  else if(type == "prob") { 1/(1 + exp(-2 * f)) }
}



####################################################################################################
# FlexBoost
####################################################################################################

FlexBoost <- function(X, y, n_rounds, par.k, type, acc.plt, control = rpart.control(cp = -1, maxdepth = 1)){
  
  
  # set 3 ways of parameter k on exp loss function
  k.list      <- c(1/par.k, 1, par.k)
  
  
  # count the number of rows
  n           <- nrow(X)
  
  
  # save parameter k path, globalize path to see after function
  k.path      <- list()
  
  
  # initialize (weight, tree, alpha) on each data
  w           <- list(rep(1/n, n))
  trees       <- list()
  alphas      <- list()
  
  
  # save (weight, tree, alpha) of 3 ways
  temp.w      <- list()
  temp.trees  <- list()
  temp.alphas <- list()
  
  
  # save train accuracy of 3 ways to compare
  temp.result <- list()
  
  # save train accuracy
  acc.result  <- list()
  
  
  # build weak classifiers
  for(i in seq(n_rounds)){
    
    tree <- rpart::rpart(y ~ .,
                         data = data.frame(X), weights = w[[i]],
                         method = "class", control = control,
                         x = FALSE, y = FALSE, model = TRUE)
    
    
    pred      <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
    
    
    # calculate the error of each classifiers
    e         <- sum(w[[i]] * (pred != y))
    
    
    # if error >= 0.5, flip the result
    if(e >= 0.5){e <- 1 - e}
    
    if(e < 1e-5){
      # if first base classifier predicts data perfectly, boosting process ends
      if(i == 1){
        # first base classifier's weight should be 1
        alphas[[i]]     <- 1
        trees[[i]]      <- tree
        terms           <- tree$terms
        acc.result[[i]] <- 1
        break
      }
      break
    }
    
    # count number of 3 ways
    n_count  <- 0
    
    for(i1 in k.list){
      
      # count number of 3 ways
      n_count <- n_count + 1
      
      
      # update weak classifier weight
      alpha   <- 1/(2*i1) * log((1-e)/e)
      
      
      # update and normalize weight of each data
      # save weight of 3ways
      temp.w[[n_count]]     <- w[[i]] * exp(-alpha*pred*y*i1)
      temp.w[[n_count]]     <- temp.w[[n_count]] / sum(temp.w[[n_count]])
      
      
      # Remove formulas since they waste memory
      if(i == 1)  { terms       <- tree$terms }
      
      else        { tree$terms  <- NULL }
      
      
      alphas[[i]] <- alpha
      trees[[i]]  <- tree
      
      # save alpha, tree of 3 ways
      temp.alphas[[n_count]] <- alpha
      temp.trees[[n_count]]  <-tree
      
      
      result                 <- list(terms  = terms,
                                     trees  = trees,
                                     alphas = unlist(alphas))
      
      class(result) <- "FlexBoost"
      
      y_hat                   <- stats::predict(result, X)
      
      
      # save train accuracy
      temp.result[[n_count]]  <- sum(y == y_hat) / length(y)
    }
    
    # clear waste memory
    w[[i]] <- NULL
    
    # if first classifier splits perfect, stop iteration
    if(max(unlist(temp.result)) == 1){
      acc.result[[i]] <- 1
      break
    }
    
    acc.result[[i]] <- max(unlist(temp.result)) 
    
    if (type == 1){
      ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
      
      # (k > 1) => max
      if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k < 1) => max
      else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1) => max
      else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k > 1, k < 1) => max
      else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1, k < 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k = 1, k < 1) => max
      else{
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
      }
      
      
      # initialize each 3 ways value
      temp.w       <- list()
      temp.alphas  <- list()
      temp.trees   <- list()
      temp.result  <- list()
      
    }
    
    else if(type == 2){
      ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
      
      # (k > 1) => max
      if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k < 1) => max
      else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1) => max
      else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k < 1) => max
      else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k = 1, k < 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k = 1, k < 1) => max
      else{
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
      }
      
      
      # initialize each 3 ways value
      temp.w       <- list()
      temp.alphas  <- list()
      temp.trees   <- list()
      temp.result  <- list()
      
    }
    
    else{
      ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
      
      # (k > 1) => max
      if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k < 1) => max
      else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1) => max
      else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k < 1) => max
      else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k = 1, k < 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1, k < 1) => max
      else{
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
      }
      
      
      # initialize each 3 ways value
      temp.w       <- list()
      temp.alphas  <- list()
      temp.trees   <- list()
      temp.result  <- list()
      
    }
    
    
    
    
  }
  
  result        <- list(terms  = terms,
                        trees  = trees,
                        alphas = unlist(alphas),
                        acc = unlist(acc.result))
  
  class(result) <- "FlexBoost"
  
  if(acc.plt == TRUE){
    par(bg = "black")
    par(col = "white")
    
    plot(1:length(alphas), 1 - unlist(acc.result), xlab = "Iteration", ylab = "Train Error", ylim = c(1 - max(unlist(acc.result)), 1 - min(unlist(acc.result))),
         col = "green",
         col.lab = "yellow",
         col.axis = "white")
    lines(1:length(alphas), 1 - unlist(acc.result), col = "green")
  }
  
  
  
  return(result)
}

FlexBoost(X, y,
          n_rounds = 50, # Iteration number
          par.k = 0.3,   # K :: (0 ~ 1)
          type = 3,      # When acc tie occurs,  1 :: k > 1, 2 :: k = 1,  3 :: k < 1 
          acc.plt = TRUE)# Plot Train error
####################################################################################################
# If you have any question about this logic or something,
# Contact :: flywade@skku.edu
####################################################################################################