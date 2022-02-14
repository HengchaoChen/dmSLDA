library(MASS)
library(tidyverse)
library(egg)
library(gridExtra)
library(latex2exp)
library(cowplot)
library(lpSolve)

#---- data.generation ----
data.generation <- function(M, b, mu, Sigma){
  p = nrow(Sigma)
  K = ncol(mu)
  n = K * b
  
  data = list()
  for(m in 1:M){
    data.m = matrix(rep(0, n * p), nrow = n)
    for(k in 1:K){
      data.m[((k-1)*b + 1):(k*b),] = mvrnorm(n = b, mu[,k], Sigma)
    }
    data[[m]] = data.m
  }
  
  mu.list <- list()
  U.list <- list()
  Sigma.list <- list()
  for(m in 1:M){
    mu.m = matrix(rep(0, p * K), nrow = p)
    for(k in 1:K){
      mu.m[,k] = apply(data[[m]][((k-1)*b + 1):(k*b),], 2, mean)
    }
    mu.list[[m]] = mu.m
    
    U.list[[m]] = mu.to.U(mu.m)
    
    data.m.cen = data[[m]]
    for(k in 1:K){
      for(i in ((k-1)*b + 1):(k*b)){
        data.m.cen[i,] = data[[m]][i,] - mu.m[,k]
      }
    }
    Sigma.list[[m]] = t(data.m.cen) %*% data.m.cen / (nrow(data.m.cen) - K)
  }
  return(list(Sigma.list, U.list, data))
}

mu.to.U <- function(mu){
  # balanced
  K = ncol(mu)
  p = nrow(mu)
  mu.bar = apply(mu, 1, mean)
  
  U = matrix(rep(0, (K - 1) * p), nrow = p)
  for(k in 1:(K-1)){
    U[,k] = mu[,k] - mu.bar
  }
  return(U)
}

#---- cen.SLDA ----
cen.SLDA <- function(data, K){
  # balanced setting
  Sigma.list = data[[1]]
  U.list = data[[2]]
  raw.data = data[[3]]
  M = length(Sigma.list)
  n = nrow(raw.data[[1]])
  p = ncol(raw.data[[1]])
  b = n / K
  N = n * M
  # X: data.aggregated
  X = matrix(rep(0, K * b * M * p), ncol = p)
  for(k in 1:K){
    for(m in 1:M){
      X[((k-1)*M*b + (m-1)*b + 1):((k-1)*M*b + m*b),] = raw.data[[m]][((k-1)*b + 1):(k*b),]
    }
  }
  
  # lambda: 3-fold CV
  {
    # three splits
    # split 1
    b.train = round(2/3*M*b)
    b.valid = M*b - b.train
    X.train = matrix(rep(0, b.train * K * p), ncol = p)
    X.valid = matrix(rep(0, b.valid * K * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*M*b + 1):((k-1)*M*b + b.train),]
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*M*b + b.train + 1):(k*M*b),]
    }
    
    mu.train = matrix(rep(0, p * K), nrow = p)
    Sigma.train = matrix(rep(0, p * p), nrow = p)
    X.train.cen = X.train
    for(k in 1:K){
      mu.train[,k] = apply(X.train[((k-1)*b.train + 1):(k*b.train),], 2, mean)
      X.train.cen[((k-1)*b.train + 1):(k*b.train),] = t(t(X.train[((k-1)*b.train + 1):(k*b.train),]) - mu.train[,k])
    }
    Sigma.train = t(X.train.cen) %*% X.train.cen / (nrow(X.train.cen) - K)
    U.train = mu.to.U(mu.train)
    
    mu.valid = matrix(rep(0, p * K), nrow = p)
    Sigma.valid = matrix(rep(0, p * p), nrow = p)
    X.valid.cen = X.valid
    for(k in 1:K){
      mu.valid[,k] = apply(X.valid[((k-1)*b.valid + 1):(k*b.valid),], 2, mean)
      X.valid.cen[((k-1)*b.valid + 1):(k*b.valid),] = t(t(X.valid[((k-1)*b.valid + 1):(k*b.valid),]) - mu.valid[,k])
    }
    Sigma.valid = t(X.valid.cen) %*% X.valid.cen / (nrow(X.valid.cen) - K)
    U.valid = mu.to.U(mu.valid)
    
    X.train.1 = X.train
    X.valid.1 = X.valid
    Sigma.train.1 = Sigma.train
    U.train.1 = U.train
    Sigma.valid.1 = Sigma.valid
    U.valid.1 = U.valid
    
    # split 2
    b.train = M*b - round(1/3*M*b)
    b.valid = M*b - b.train
    X.train = matrix(rep(0, b.train * K * p), ncol = p)
    X.valid = matrix(rep(0, b.valid * K * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*M*b + round(1/3*M*b) + 1):(k*M*b),]
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*M*b + 1):((k-1)*M*b + b.valid),]
    }
    
    mu.train = matrix(rep(0, p * K), nrow = p)
    Sigma.train = matrix(rep(0, p * p), nrow = p)
    X.train.cen = X.train
    for(k in 1:K){
      mu.train[,k] = apply(X.train[((k-1)*b.train + 1):(k*b.train),], 2, mean)
      X.train.cen[((k-1)*b.train + 1):(k*b.train),] = t(t(X.train[((k-1)*b.train + 1):(k*b.train),]) - mu.train[,k])
    }
    Sigma.train = t(X.train.cen) %*% X.train.cen / (nrow(X.train.cen) - K)
    U.train = mu.to.U(mu.train)
    
    mu.valid = matrix(rep(0, p * K), nrow = p)
    Sigma.valid = matrix(rep(0, p * p), nrow = p)
    X.valid.cen = X.valid
    for(k in 1:K){
      mu.valid[,k] = apply(X.valid[((k-1)*b.valid + 1):(k*b.valid),], 2, mean)
      X.valid.cen[((k-1)*b.valid + 1):(k*b.valid),] = t(t(X.valid[((k-1)*b.valid + 1):(k*b.valid),]) - mu.valid[,k])
    }
    Sigma.valid = t(X.valid.cen) %*% X.valid.cen / (nrow(X.valid.cen) - K)
    U.valid = mu.to.U(mu.valid)
    
    X.train.2 = X.train
    X.valid.2 = X.valid
    Sigma.train.2 = Sigma.train
    U.train.2 = U.train
    Sigma.valid.2 = Sigma.valid
    U.valid.2 = U.valid
    
    # split 3
    b.train = round(1/3*M*b) + M*b - round(2/3*M*b)
    b.valid = M*b - b.train
    X.train = matrix(rep(0, b.train * K * p), ncol = p)
    X.valid = matrix(rep(0, b.valid * K * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = rbind(X[((k-1)*M*b + 1):((k-1)*M*b + round(1/3*M*b)),], X[((k-1)*M*b + round(2/3*M*b) + 1):(k*M*b),])
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*M*b + round(1/3*M*b) + 1):((k-1)*M*b + round(2/3*M*b)),]
    }
    
    mu.train = matrix(rep(0, p * K), nrow = p)
    Sigma.train = matrix(rep(0, p * p), nrow = p)
    X.train.cen = X.train
    for(k in 1:K){
      mu.train[,k] = apply(X.train[((k-1)*b.train + 1):(k*b.train),], 2, mean)
      X.train.cen[((k-1)*b.train + 1):(k*b.train),] = t(t(X.train[((k-1)*b.train + 1):(k*b.train),]) - mu.train[,k])
    }
    Sigma.train = t(X.train.cen) %*% X.train.cen / (nrow(X.train.cen) - K)
    U.train = mu.to.U(mu.train)
    
    mu.valid = matrix(rep(0, p * K), nrow = p)
    Sigma.valid = matrix(rep(0, p * p), nrow = p)
    X.valid.cen = X.valid
    for(k in 1:K){
      mu.valid[,k] = apply(X.valid[((k-1)*b.valid + 1):(k*b.valid),], 2, mean)
      X.valid.cen[((k-1)*b.valid + 1):(k*b.valid),] = t(t(X.valid[((k-1)*b.valid + 1):(k*b.valid),]) - mu.valid[,k])
    }
    Sigma.valid = t(X.valid.cen) %*% X.valid.cen / (nrow(X.valid.cen) - K)
    U.valid = mu.to.U(mu.valid)
    
    X.train.3 = X.train
    X.valid.3 = X.valid
    Sigma.train.3 = Sigma.train
    U.train.3 = U.train
    Sigma.valid.3 = Sigma.valid
    U.valid.3 = U.valid
  }  
  lambda.list <- exp(seq(-4, 0, length.out = 81))
  loss.list <- rep(0, 81)
  for(i in 1:81){
    lambda = lambda.list[i]
    W.train.1 = fista(Sigma.train.1, U.train.1, lambda)
    W.train.2 = fista(Sigma.train.2, U.train.2, lambda)
    W.train.3 = fista(Sigma.train.3, U.train.3, lambda)
    loss.list[i] = sum(diag(t(W.train.1) %*% (Sigma.valid.1 %*% W.train.1))) / 2 - sum(diag(t(W.train.1) %*% U.valid.1)) +
      sum(diag(t(W.train.2) %*% (Sigma.valid.2 %*% W.train.2))) / 2 - sum(diag(t(W.train.2) %*% U.valid.2)) +
      sum(diag(t(W.train.2) %*% (Sigma.valid.2 %*% W.train.2))) / 2 - sum(diag(t(W.train.2) %*% U.valid.2))
  }
  
  lambda = lambda.list[which.min(loss.list)]
  
  Sigma = matrix(rep(0, p * p), nrow = p)
  U = matrix(rep(0, p * (K - 1)), nrow = p)
  for(m in 1:M){
    Sigma <- Sigma + Sigma.list[[m]] / M
    U <- U + U.list[[m]] / M
  }
  W <- fista(Sigma, U, lambda)
  
  return(W)
}

#---- dmSLDA ----
dmSLDA <- function(data, K, comm = 3, whole.path = FALSE){
  # balanced setting
  Sigma.list = data[[1]]
  U.list = data[[2]]
  raw.data = data[[3]]
  M = length(Sigma.list)
  n = nrow(raw.data[[1]])
  p = ncol(raw.data[[1]])
  b = n / K
  N = n * M
  
  # not used
  {
    # 3-fold
    X = matrix(rep(0, K * b * M * p), ncol = p)
    for(k in 1:K){
      for(m in 1:M){
        X[((k-1)*M*b + (m-1)*b + 1):((k-1)*M*b + m*b),] = raw.data[[m]][((k-1)*b + 1):(k*b),]
      }
    }
    {
      # three splits
      # split 1
      b.train = round(2/3*M*b)
      b.valid = M*b - b.train
      X.train = matrix(rep(0, b.train * K * p), ncol = p)
      X.valid = matrix(rep(0, b.valid * K * p), ncol = p)
      for(k in 1:K){
        X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*M*b + 1):((k-1)*M*b + b.train),]
        X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*M*b + b.train + 1):(k*M*b),]
      }
      
      X.train.1 = X.train
      X.valid.1 = X.valid
      
      # split 2
      b.train = M*b - round(1/3*M*b)
      b.valid = M*b - b.train
      X.train = matrix(rep(0, b.train * K * p), ncol = p)
      X.valid = matrix(rep(0, b.valid * K * p), ncol = p)
      for(k in 1:K){
        X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*M*b + round(1/3*M*b) + 1):(k*M*b),]
        X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*M*b + 1):((k-1)*M*b + b.valid),]
      }
      
      X.train.2 = X.train
      X.valid.2 = X.valid
      
      # split 3
      b.train = round(1/3*M*b) + M*b - round(2/3*M*b)
      b.valid = M*b - b.train
      X.train = matrix(rep(0, b.train * K * p), ncol = p)
      X.valid = matrix(rep(0, b.valid * K * p), ncol = p)
      for(k in 1:K){
        X.train[((k-1)*b.train + 1):(k*b.train),] = rbind(X[((k-1)*M*b + 1):((k-1)*M*b + round(1/3*M*b)),], X[((k-1)*M*b + round(2/3*M*b) + 1):(k*M*b),])
        X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*M*b + round(1/3*M*b) + 1):((k-1)*M*b + round(2/3*M*b)),]
      }
      
      X.train.3 = X.train
      X.valid.3 = X.valid
    }  
  }
  
  W.list = list()
  lambda.array = rep(0, comm + 1)
  loss.array = rep(0, comm + 1)
  score.array = rep(0, comm + 1)
  
  # initial estimator
  lambda.list = exp(seq(-4.5, -0.5, length.out = 161))
  loss.list = rep(0, 161)
  for(i in 1:161){
    lambda = lambda.list[i]
    W = fista(Sigma.list[[1]], U.list[[1]], lambda)
    loss.list[i] = valid(Sigma.list, U.list, W)
    #loss.list[i] = classifier(W, X.train.1, X.valid.1) + classifier(W, X.train.2, X.valid.2) + classifier(W, X.train.3, X.valid.3)
  }
  lambda.array[1] <-lambda.list[which.min(loss.list)]
  loss.array[1] <- min(loss.list)
  W.list[[1]] <- fista(Sigma.list[[1]], U.list[[1]], lambda.array[1])
  #score.array[1] <- classifier(W.list[[1]], X.train, X.valid)
  #score.array[1] <- valid(Sigma.list, U.list, W.list[[1]])
  score.array[1] <- loss.array[1]
  
  if(comm == 0){
    if(whole.path == TRUE){
      return(list(W.list, lambda.array, loss.array, score.array))
    }
    else if(whole.path == FALSE){
      return(W.list[[1]])
    }
  }
  
  # iterations
  for(t in 1:comm){
    U.shifted <- Sigma.list[[1]] %*% W.list[[t]]
    for(m in 1:M){
      U.shifted <- U.shifted - (Sigma.list[[m]] %*% W.list[[t]] - U.list[[m]]) / M
    }
    
    lambda.list <- exp(seq(log(lambda.array[t]) - 2, log(lambda.array[t]) + 1, length.out = 61))
    loss.list <- rep(0, 61)
    for(i in 1:61){
      lambda <- lambda.list[i]
      W <- fista(Sigma.list[[1]], U.shifted, lambda)
      loss.list[i] <- valid(Sigma.list, U.list, W)
      #loss.list[i] <- classifier(W, X.train.1, X.valid.1) + classifier(W, X.train.2, X.valid.2) + classifier(W, X.train.3, X.valid.3)
    }
    lambda.array[t + 1] <- lambda.list[which.min(loss.list)]
    loss.array[t + 1] <- min(loss.list)
    
    W.list[[t + 1]] <- fista(Sigma.list[[1]], U.shifted, lambda.array[t + 1])
    #score.array[t + 1] <- classifier(W.list[[t + 1]], X.train, X.valid)
    #score.array[t + 1] <- valid(Sigma.list, U.list, W.list[[t + 1]])
    score.array[t + 1] <- loss.array[t + 1]
  }
  
  t.optimal <- which.min(score.array)
  
  if(whole.path == TRUE){
    return(list(W.list, lambda.array, loss.array, score.array))
  }
  else if(whole.path == FALSE){
    return(W.list[[t.optimal]])
  }
}

#---- classifier ----
classifier <- function(W.hat, data.train, data.test){
  p <- nrow(W.hat)
  K <- ncol(W.hat) + 1
  # abandon some W
  if(sum(abs(W.hat)) > 10 * p){
    return(1)
  }
  
  # W.hat -> W.tilde: delete \zero
  W.tilde <- matrix(rep(0, p * (K - 1)), nrow = p)
  count <- 0
  for(i in 1:(K-1)){
    if(is.equal(W.hat[,i], rep(0, p)) == FALSE){
      count <- count + 1
      W.tilde[,count] <- W.hat[,i]
    }
  }
  W.tilde <- W.tilde[,1:count]
  
  if(count == 0){return(1)}
  
  # LDA rule
  n = nrow(data.train)
  b = n / K
  # balanced setting
  data.trans.train = data.train %*% W.tilde
  
  mu.trans = matrix(rep(0, count * K), nrow = count)
  for(k in 1:K){
    if(count == 1){
      mu.trans[,k] <- mean(data.trans.train[((k-1)*b + 1):(k*b),])
    }
    else if(count > 1){
      mu.trans[,k] <- apply(data.trans.train[((k-1)*b + 1):(k*b),], 2, mean)
    }
  }
  
  data.trans.train.cen <- data.trans.train
  for(k in 1:K){
    for(i in ((k-1)*b + 1):(k*b)){
      data.trans.train.cen[i,] = data.trans.train[i,] - mu.trans[,k]
    }
  }
  
  Sigma.trans.train = t(data.trans.train.cen) %*% data.trans.train.cen / n
  if(det(Sigma.trans.train) == 0){
    return(1)
  }
  
  # mcr on data.test
  n.test = nrow(data.test)
  b.test = n.test / K
  # balanced setting
  
  data.trans.test = data.test %*% W.tilde
  class.test = rep(0, n.test)
  for(i in 1:n.test){
    score = rep(0, K)
    for(k in 1:K){
      score[k] = data.trans.test[i,] %*% solve(Sigma.trans.train) %*% mu.trans[,k] - t(mu.trans[,k]) %*% solve(Sigma.trans.train) %*% mu.trans[,k] / 2 + log(1/K)
    }
    class.test[i] = which.max(score)
  }
  
  mcr = 0
  for(k in 1:K){
    for(i in ((k-1)*b.test + 1):(k*b.test)){
      if(class.test[i] != k){
        mcr = mcr + 1/n.test
      }
    }
  }
  
  return(mcr)
}

#---- useful functions ----

fista <- function(Sigma, U, lambda, T = 150){
  if(is.null(ncol(U)) == TRUE){
    K = 2
  }
  else{
    K = ncol(U) + 1
  }
  p = nrow(Sigma)
  L = max(eigen(Sigma)$values)
  
  W.1 <- MASS::ginv(Sigma + diag(rep(0.01, p))) %*% U
  V <- W.1
  a.1 <- 1
  
  for(t in 1:T){
    W.2 <- V - (Sigma %*% V - U) / L
    for(i in 1:p){
      for(j in 1:(K-1)){
        W.2[i,j] <- sign(W.2[i,j]) * max(abs(W.2[i,j]) - lambda / L, 0)
      }
    }
    
    a.2 <- (1 + sqrt(1 + 4 * a.1 ^ 2)) / 2
    
    V <- W.2 + (a.1 - 1) / a.2 * (W.2 - W.1)
    W.1 <- W.2
    a.1 <- a.2
  }
  
  return(W.1)
}

valid <- function(Sigma.list, U.list, W){
  M = length(U.list)
  if(M > 1){
  #p = nrow(W)
  #K = ncol(W) + 1
  #E = matrix(rep(0, p * (K-1)), nrow = p)
  #for(m in 2:M){
  #  E = E + Sigma.list[[m]] %*% W - U.list[[m]]
  #}
  #loss = sum(diag(t(E) %*% E))
  loss = 0
  for(m in 2:M){
    loss = loss + sum(diag(t(W) %*% Sigma.list[[m]] %*% W)) / 2 - sum(diag(t(W) %*% U.list[[m]]))
  }
  }
  else if(M == 1){
    loss = sum(diag(t(W) %*% Sigma.list[[1]] %*% W)) / 2 - sum(diag(t(W) %*% U.list[[1]]))
    warning("There is only one machine.")
  }
  return(loss)
}

is.equal <- function(a, b){
  n = length(a)
  for(i in 1:n){
    if(a[i] != b[i]){
      return(FALSE)
    }
  }
  return(TRUE)
}

#---- oracle.dmSLDA ----
dmSLDA.ora <- function(data, K, W.oracle, comm = 3, whole.path = FALSE){
  # balanced setting
  Sigma.list = data[[1]]
  U.list = data[[2]]
  raw.data = data[[3]]
  M = length(Sigma.list)
  n = nrow(raw.data[[1]])
  p = ncol(raw.data[[1]])
  b = n / K
  N = n * M
  
  W.list = list()
  lambda.array = rep(0, comm + 1)
  loss.array = rep(0, comm + 1)
  score.array = rep(0, comm + 1)
  
  # initial estimator
  lambda.list = exp(seq(-4.5, -0.5, length.out = 81))
  loss.list = rep(0, 81)
  for(i in 1:81){
    lambda = lambda.list[i]
    W = fista(Sigma.list[[1]], U.list[[1]], lambda, T = 150)
    #loss.list[i] = valid(Sigma.list, U.list, W)
    #loss.list[i] = classifier(W, X.train, X.valid)
    loss.list[i] = sum(diag(t(W - W.oracle) %*% (W - W.oracle)))
  }
  lambda.array[1] <-lambda.list[which.min(loss.list)]
  loss.array[1] <- min(loss.list)
  W.list[[1]] <- fista(Sigma.list[[1]], U.list[[1]], lambda.array[1], T = 150)
  #score.array[1] <- classifier(W.list[[1]], X.train, X.valid)
  score.array[1] <- loss.array[1]
  #score.array[1] <- valid(Sigma.list, U.list, W.list[[1]])
  
  if(comm == 0){
    if(whole.path == TRUE){
      return(list(W.list, lambda.array, loss.array, score.array))
    }
    else if(whole.path == FALSE){
      return(W.list[[1]])
    }
  }
  
  # iteration
  for(t in 1:comm){
    U.shifted <- Sigma.list[[1]] %*% W.list[[t]]
    for(m in 1:M){
      U.shifted <- U.shifted - (Sigma.list[[m]] %*% W.list[[t]] - U.list[[m]]) / M
    }
    
    lambda.list <- exp(seq(log(lambda.array[t]) - 2, log(lambda.array[t]) + 1, length.out = 61))
    loss.list <- rep(0, 61)
    for(i in 1:61){
      lambda <- lambda.list[i]
      W <- fista(Sigma.list[[1]], U.shifted, lambda, T = 150)
      #loss.list[i] <- valid(Sigma.list, U.list, W)
      #loss.list[i] <- classifier(W, X.train, X.valid)
      loss.list[i] = sum(diag(t(W - W.oracle) %*% (W - W.oracle)))
    }
    lambda.array[t + 1] <- lambda.list[which.min(loss.list)]
    loss.array[t + 1] <- min(loss.list)
    
    W.list[[t + 1]] <- fista(Sigma.list[[1]], U.shifted, lambda.array[t + 1], T = 150)
    #score.array[t + 1] <- classifier(W.list[[t + 1]], X.train, X.valid)
    score.array[t + 1] <- loss.array[t + 1]
    #score.array[t + 1] <- valid(Sigma.list, U.list, W.list[[t + 1]])
  }
  
  t.optimal <- which.min(score.array)
  
  if(whole.path == TRUE){
    return(list(W.list, lambda.array, loss.array, score.array))
  }
  else if(whole.path == FALSE){
    return(W.list[[t.optimal]])
  }
}

#---- dc.SLDA ----
LPD <- function(Sigma, mu.dif, lambda){
  p = nrow(Sigma)
  f.obj = rep(1, 2 * p)
  f.con.1 = diag(rep(1, 2 * p))
  f.con.2 = cbind(Sigma, - Sigma)
  f.con = rbind(f.con.1, f.con.2, f.con.2)
  f.dir = c(rep(">=", 3 * p), rep("<=", p))
  f.rhs = c(rep(0, 2 * p), mu.dif - lambda, mu.dif + lambda)
  s = lp("min", f.obj, f.con, f.dir, f.rhs)$solution
  beta = s[1:p] - s[(p+1):(2*p)]
  
  return(beta)
}
cen.SLDA.bin <- function(X){
  # balanced
  b = nrow(X) / 2
  p = ncol(X)
  K = 2
  n = K * b
  
  # lambda: 3-fold
  {
    # split 1
    b.valid = round(1/3*b)
    b.train = b - b.valid  
    X.train = matrix(rep(0, K * b.train * p), ncol = p)
    X.valid = matrix(rep(0, K * b.valid * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*b + round(1/3*b) + 1):(k*b),]
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*b + 1):((k-1)*b + round(1/3*b)),]
    }
    
    # binary
    mu.train = matrix(rep(0, p * K), nrow = p)
    X.train.cen = X.train
    for(k in 1:K){
      mu.train[,k] = apply(X.train[((k-1)*b.train + 1):(k*b.train),], 2, mean)
      X.train.cen[((k-1)*b.train + 1):(k*b.train),] = t(t(X.train[((k-1)*b.train + 1):(k*b.train),]) - mu.train[,k])
    }
    Sigma.train = t(X.train.cen) %*% X.train.cen / (nrow(X.train.cen) - K)
    mu.dif.train = mu.train[,1] - mu.train[,2]
    
    X.train.1 = X.train
    X.valid.1 = X.valid
    Sigma.train.1 = Sigma.train
    mu.dif.train.1 = mu.dif.train
    mu.bar.1 = apply(X.train, 2, mean)
    
    # split 2
    b.valid = round(2/3*b) - round(1/3*b)
    b.train = b - b.valid 
    X.train = matrix(rep(0, K * b.train * p), ncol = p)
    X.valid = matrix(rep(0, K * b.valid * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = rbind(X[((k-1)*b + 1):((k-1)*b + round(1/3*b)),], X[((k-1)*b + round(2/3*b) + 1):(k*b),])
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*b + round(1/3*b) + 1):((k-1)*b + round(2/3*b)),]
    }
    
    # binary
    mu.train = matrix(rep(0, p * K), nrow = p)
    X.train.cen = X.train
    for(k in 1:K){
      mu.train[,k] = apply(X.train[((k-1)*b.train + 1):(k*b.train),], 2, mean)
      X.train.cen[((k-1)*b.train + 1):(k*b.train),] = t(t(X.train[((k-1)*b.train + 1):(k*b.train),]) - mu.train[,k])
    }
    Sigma.train = t(X.train.cen) %*% X.train.cen / (nrow(X.train.cen) - K)
    mu.dif.train = mu.train[,1] - mu.train[,2]
    
    X.train.2 = X.train
    X.valid.2 = X.valid
    Sigma.train.2 = Sigma.train
    mu.dif.train.2 = mu.dif.train
    mu.bar.2 = apply(X.train, 2, mean)
    
    # split 3
    b.valid = round(3/3*b) - round(2/3*b)
    b.train = b - b.valid 
    X.train = matrix(rep(0, K * b.train * p), ncol = p)
    X.valid = matrix(rep(0, K * b.valid * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*b + 1):((k-1)*b + round(2/3*b)),]
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*b + round(2/3*b) + 1):((k-1)*b + round(3/3*b)),]
    }
    
    # binary
    mu.train = matrix(rep(0, p * K), nrow = p)
    X.train.cen = X.train
    for(k in 1:K){
      mu.train[,k] = apply(X.train[((k-1)*b.train + 1):(k*b.train),], 2, mean)
      X.train.cen[((k-1)*b.train + 1):(k*b.train),] = t(t(X.train[((k-1)*b.train + 1):(k*b.train),]) - mu.train[,k])
    }
    Sigma.train = t(X.train.cen) %*% X.train.cen / (nrow(X.train.cen) - K)
    mu.dif.train = mu.train[,1] - mu.train[,2]
    
    X.train.3 = X.train
    X.valid.3 = X.valid
    Sigma.train.3 = Sigma.train
    mu.dif.train.3 = mu.dif.train
    mu.bar.3 = apply(X.train, 2, mean)
  }
  
  lambda.list = exp(seq(-4, 0, length.out = 41))
  loss.list = rep(0, 41)
  for(i in 1:41){
    lambda = lambda.list[i]
    beta.1 = LPD(Sigma.train.1, mu.dif.train.1, lambda)
    beta.2 = LPD(Sigma.train.2, mu.dif.train.2, lambda)
    beta.3 = LPD(Sigma.train.3, mu.dif.train.3, lambda)
    loss.list[i] = classifier.bi(beta.1, mu.bar.1, X.valid.1) +
      classifier.bi(beta.2, mu.bar.2, X.valid.2) +
      classifier.bi(beta.3, mu.bar.3, X.valid.3)
  }
  lambda = lambda.list[which.min(loss.list)]
  
  mu = matrix(rep(0, p * K), nrow = p)
  X.cen = X
  for(k in 1:K){
    mu[,k] = apply(X[((k-1)*b + 1):(k*b),], 2, mean)
    X.cen[((k-1)*b + 1):(k*b),] = t(t(X[((k-1)*b + 1):(k*b),]) - mu[,k])
  }
  Sigma = t(X.cen) %*% X.cen / (nrow(X.cen) - K)
  mu.dif = mu[,1] - mu[,2]
  beta = LPD(Sigma, mu.dif, lambda)
  
  mu.bar = apply(X, 2, mean)
  
  return(list(beta, mu.bar, lambda, Sigma, mu.dif))
}
inverse.Sigma <- function(Sigma, lambda){
  p = nrow(Sigma)
  theta = matrix(rep(0, p * p), nrow = p)
  for(i in 1:p){
    e = rep(0, p)
    e[i] = 1
    f.obj = rep(1,  2 * p)
    f.con.1 = diag(rep(-1, 2 * p))
    f.con.2 = cbind(Sigma, - Sigma)
    f.con = rbind(f.con.1, f.con.2, - f.con.2)
    f.dir = rep("<=", 4 * p)
    f.rhs = c(rep(0, 2 * p), e + lambda, lambda - e)
    s = lp("min", f.obj, f.con, f.dir, f.rhs)$solution
    theta.i = s[1:p] - s[(p+1):(2*p)]
    theta[,i] = theta.i
  }
  return(theta)
}
dSLDA.bin <- function(data){
  X.list = data[[3]]
  M = length(X.list)
  p = ncol(X.list[[1]])
  beta.list = list()
  s.vector = rep(0, M)
  for(m in 1:M){
    X = X.list[[m]]
    result = cen.SLDA.bin(X)
    beta.0 = result[[1]]
    lambda = result[[3]]
    Sigma = result[[4]]
    mu.dif = result[[5]]
    theta.0 = inverse.Sigma(Sigma, lambda)
    beta.debiased = beta.0 - theta.0 %*% (Sigma %*% beta.0 - mu.dif)
    beta.list[[m]] = beta.debiased
    s.vector[m] = length(beta.0[beta.0 != 0])
  }
  
  beta = rep(0, p)
  for(m in 1:M){
    beta = beta + beta.list[[m]] / M
  }
  s = max(s.vector)
  t = sort(abs(beta), decreasing = TRUE)[s]
  
  for(i in 1:p){
    if(abs(beta)[i] < t){
      beta[i] = 0
    }
  }

  return(beta)
}
classifier.bi <- function(beta, mu.bar, X){
  # balanced
  # binary
  b = nrow(X) / 2
  p = ncol(X)
  K = 2
  n = K * b
  
  mis.index = rep(0, n)
  for(i in 1:b){
    score.1 = t(X[i,] - mu.bar) %*% beta
    if(score.1 <= 0){
      mis.index[i] = 1
    }
    score.2 = t(X[(b+i),] - mu.bar) %*% beta
    if(score.2 >= 0){
      mis.index[b+i] = 1
    }
  }
  
  mis = sum(mis.index)
  
  return(mis)
}

#---- dc.SVM ----
dSVM <- function(data){
  # balanced
  # centered
  X.list = data[[3]]
  M = length(X.list)
  p = ncol(X.list[[1]])
  n = nrow(X.list[[1]])
  K = 2
  b = n / K
  beta.list = list()
  s.vector = rep(0, M)
  for(m in 1:M){
    X = X.list[[m]]
    result = cen.l1.SVM(X)
    beta.0 = result[[1]]
    H = result[[2]]
    lambda = result[[3]]
    theta.0 = inverse.Sigma(H, lambda)
    
    u = rep(0, p)
    for(i in 1:(n/2)){
      u = u + X[i,] * max(0, 1 - X[i,] %*% beta.0) / n 
    }
    for(i in (n/2+1):n){
      u = u + X[i,] * max(0, 1 + X[i,] %*% beta.0) / n 
    }
    
    beta.debiased = beta.0 - theta.0 %*% u
    beta.list[[m]] = beta.debiased
    s.vector[m] = length(beta.0[beta.0 != 0])
  }
  
  beta = rep(0, p)
  for(m in 1:M){
    beta = beta + beta.list[[m]] / M
  }
  s = max(s.vector)
  t = sort(abs(beta), decreasing = TRUE)[s]
  
  for(i in 1:p){
    if(abs(beta)[i] < t){
      beta[i] = 0
    }
  }
  
  return(beta)
}
LP.svm <- function(X, lambda){
  # balanced
  n = nrow(X)
  K = 2
  b = n / K
  p = ncol(X)
  
  f.obj = c(rep(1/n, n), rep(0, p))
  f.con.1 = cbind(diag(rep(1, n)), matrix(rep(0, n*p), nrow = n))
  f.con.2 = cbind(diag(rep(1, n/2)), diag(rep(0, n/2)), X[1:(n/2),])
  f.con.3 = cbind(diag(rep(0, n/2)), diag(rep(1, n/2)), -X[(n/2+1):n,])
  f.con.4 = cbind(matrix(rep(0, p*n), nrow = p), diag(rep(1, p)))
  f.con = rbind(f.con.1, f.con.2, f.con.3, f.con.4, f.con.4)
  f.dir = c(rep(">=", 2*n+p), rep("<=", p))
  f.rhs = c(rep(0, n), rep(1, n), rep(-lambda, p), rep(lambda, p))
  
  s = lp("min", f.obj, f.con, f.dir, f.rhs)$solution
  beta = s[(n+1):(n+p)]
  
  return(beta)
}
cen.l1.SVM <- function(X){
  # centered
  # balanced
  n = nrow(X)
  K = 2
  b = n / K
  p = ncol(X)
  
  # lambda: 3-fold
  {
    # split 1
    b.valid = round(1/3*b)
    b.train = b - b.valid  
    X.train = matrix(rep(0, K * b.train * p), ncol = p)
    X.valid = matrix(rep(0, K * b.valid * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*b + round(1/3*b) + 1):(k*b),]
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*b + 1):((k-1)*b + round(1/3*b)),]
    }
    
    X.train.1 = X.train
    X.valid.1 = X.valid
    
    # split 2
    b.valid = round(2/3*b) - round(1/3*b)
    b.train = b - b.valid 
    X.train = matrix(rep(0, K * b.train * p), ncol = p)
    X.valid = matrix(rep(0, K * b.valid * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = rbind(X[((k-1)*b + 1):((k-1)*b + round(1/3*b)),], X[((k-1)*b + round(2/3*b) + 1):(k*b),])
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*b + round(1/3*b) + 1):((k-1)*b + round(2/3*b)),]
    }
    
    X.train.2 = X.train
    X.valid.2 = X.valid
    
    # split 3
    b.valid = round(3/3*b) - round(2/3*b)
    b.train = b - b.valid 
    X.train = matrix(rep(0, K * b.train * p), ncol = p)
    X.valid = matrix(rep(0, K * b.valid * p), ncol = p)
    for(k in 1:K){
      X.train[((k-1)*b.train + 1):(k*b.train),] = X[((k-1)*b + 1):((k-1)*b + round(2/3*b)),]
      X.valid[((k-1)*b.valid + 1):(k*b.valid),] = X[((k-1)*b + round(2/3*b) + 1):((k-1)*b + round(3/3*b)),]
    }
    
    X.train.3 = X.train
    X.valid.3 = X.valid
  }
  
  lambda.list = exp(seq(-4, 0, length.out = 41))
  loss.list = rep(0, 41)
  for(i in 1:41){
    lambda = lambda.list[i]
    beta.1 = LP.svm(X.train.1, lambda)
    beta.2 = LP.svm(X.train.2, lambda)
    beta.3 = LP.svm(X.train.3, lambda)
    loss.list[i] = classifier.svm(beta.1, X.valid.1) +
      classifier.svm(beta.2, X.valid.2) +
      classifier.svm(beta.3, X.valid.3)
  }
  lambda = lambda.list[which.min(loss.list)]
  
  beta = LP.svm(X, lambda)
  
  # bandwidth h = 1.06 sigma.hat n^(-1/5)
  e = rep(0, n)
  e[1:(n/2)] = rep(1,n/2) - X[1:(n/2),] %*% beta
  e[(n/2+1):n] = rep(1, n/2) + X[(n/2+1):n,] %*% beta
  sigma.hat = sd(e)
  h = 1.06 * sigma.hat * n ^ (-1/5)
  
  # q ~ dnorm
  H = matrix(rep(0, p * p), nrow = p)
  for(i in 1:(n/2)){
    H = H + (X[i,] %*% t(X[i,])) * diag(rep(dnorm((1 - X[i,] %*% beta) / h) / (n * h), p))
  }
  for(i in (n/2+1):n){
    H = H + X[i,] %*% t(X[i,]) * diag(rep(dnorm((1 + X[i,] %*% beta) / h) / (n * h), p))
  }
  
  return(list(beta, H, lambda))
}
classifier.svm <- function(beta, X){
  # balanced
  b = nrow(X) / 2
  p = ncol(X)
  K = 2
  n = K * b
  
  mis.index = rep(0, n)
  for(i in 1:b){
    score.1 = t(X[i,]) %*% beta
    if(score.1 <= 0){
      mis.index[i] = 1
    }
    score.2 = t(X[(b+i),]) %*% beta
    if(score.2 >= 0){
      mis.index[b+i] = 1
    }
  }
  
  mis = sum(mis.index)
  
  return(mis)
}

