#'Run a logistic regression of vector Y against X using gradient descent or Newton-Raphson
#' @param X the design matrix, scaled appropriately, with first column being 1s for intercept
#' @param Y the response vector of 0s and 1s
#' @param beta the initical beta vector where the algorithm starts
#' @param alpha0 the learning rate for gradient descent
#' @param epsilon the stopping-criteria for the change in beta vector between iterations
#' @param maxit the maximum iterations for the algorithm
#' @param method the optimization method to be used; accepts NR for Newton-Raphson or GD for Grdient Descent
#' @return The function returns the vector of coefficients, their asymptotic variances, and visited log-likelihoods
logistic_regression = function(X, Y, beta, alpha0=0.001, epsilon=10e-7, maxit=400, method="NR"){
  log_likelihood = function(Y, pi){
    sum(Y * log(pi) + (1-Y) * log(1-pi))
  }
  gradient_func = function(X, Y, pi){
    apply(X * rep((pi-Y), k), 2, sum)
  }
  hessian_func = function(X, pi){
    t(X) %*% diag(as.vector(pi * (1-pi))) %*% X
  }

  ll = numeric()
  pi = exp(X %*% beta)/(1+exp(X %*% beta))
  ll[1] = sum(Y * log(pi) + (1-Y) * log(1-pi))
  n = dim(X)[1]
  k = dim(X)[2]

  for(i in 2:maxit){
    pi = exp(X %*% beta)/(1+exp(X %*% beta))
    ll[i] = sum(Y * log(pi) + (1-Y) * log(1-pi))
    gradient = gradient_func(X, Y, pi)


    if(method=="NR"){
      s_hessian = solve(hessian_func(X, pi))
      new_beta = beta - s_hessian %*% gradient
      new_pi = exp(X %*% new_beta)/(1+exp(X %*% new_beta))
      ll[i] = log_likelihood(Y, new_pi)
    }
    else{
      alpha = alpha0
      new_beta = beta - alpha * gradient
      new_pi = exp(X %*% new_beta)/(1+exp(X %*% new_beta))
      ll[i] = log_likelihood(Y, new_pi)
      while(ll[i-1] > ll[i]){
        alpha = alpha/2
        new_beta = beta - alpha * gradient
        new_pi = exp(X %*% new_beta)/(1+exp(X %*% new_beta))
        ll[i] = log_likelihood(Y, new_pi)
      }
    }
    if( sqrt(sum((new_beta - beta)^2))  < epsilon){
      break
    }
    beta = new_beta
  }
  asymp_var = diag(solve(hessian_func(X, new_pi)))

  return(list(beta=new_beta, asymp_var=asymp_var, ll=ll))
}

#'Calculate MLE population distribution of blood types from a sample
#' @param nA number of A phenotypes in sample
#' @param nAB number of AB phenotypes in sample
#' @param nB number of B phenotypes in sample
#' @param nO number of O phenotypes in sample
#' @param pA starting probability of A allele in population
#' @param pB starting probability of B allele in population
#' @param pO starting probability of O allele in population
#' @return The function returns the vector of most likely probabilities of alleles in the population
EM_func = function(nA, nAB, nB, nO, pA, pB, pO, maxit=100){
  n = nA+nAB+nB+nO

  for(i in 1:maxit){
    mAA = nA * pA^2 /(pA^2 + 2*pA*pO)
    mAO = nA * 2*pA*pO/(pA^2 + 2*pA*pO)
    mBB = nB * pB^2/(pB^2 + 2*pB*pO)
    mBO = nB * 2*pB*pO / (pB^2 + 2*pB*pO)


    pA = (2*mAA+mAO+nAB)/(2*n)
    pB = (2*mBB+mBO+nAB)/(2*n)
    pO = (2*nO+mAO+mBO)/(2*n)
  }
  return(matrix(c(pA,pB,pO), nrow=1))
}


#'Simulate from a Hidden State Markov Model
#' @param P the transition probabilities between states
#' @param E the emission probabilities
#' @param pi the starting probabilities
#' @param n the number of samples
#' @return The function will return a sample Y from the specified HMM,
#' the latent states X, and the predicted probabilities gamma
simulate_HMM  = function(P, E, pi, n=100){
  HMM_backward = function(Y, P, E){
    n = length(Y)
    b = matrix(rep(0,2*n), ncol=2)
    b[n,] = c(1,1)
    for(t in (n-1):1){
      b[t,1] = P[1,1]*E[1,Y[t+1]]*b[(t+1),1] + P[1,2]*E[2,Y[t+1]]*b[(t+1),2]
      b[t,2] = P[2,1]*E[1,Y[t+1]]*b[(t+1),1] + P[2,2]*E[2,Y[t+1]]*b[(t+1),2]
    }
    return(b)
  }
  HMM_forward = function(Y, P, E, pi){
    n = length(Y)
    a = matrix(rep(0,2*n), ncol=2)
    a[1,] = pi*E[,Y[1]]
    for(t in 1:(n-1)){
      a[(t+1),1] = E[1,Y[t+1]]*(a[t,1]*P[1,1] + a[t,2]*P[2,1])
      a[(t+1),2] = E[2,Y[t+1]]*(a[t,1]*P[1,2] + a[t,2]*P[2,2])
    }
    return(a)
  }

  X = rep(0,n)
  Y = rep(0,n)

  X[1] = sample(c(1,2), 1, prob=pi)
  Y[1] = sample(1:6, 1, prob=E[X[1],])

  for(i in 2:100){
    X[i] = sample(c(1,2), 1, prob=P[X[i-1],])
    Y[i] = sample(1:6, 1, prob=E[X[i],])
  }

  b = HMM_backward(Y, P, E)
  a = HMM_forward(Y, P, E, pi=c(1/2,1/2))

  b0 = pi[1]*E[1,Y[1]]*b[1,1] + pi[2]*E[2,Y[1]]*b[1,2]
  a0 = sum(a[100,])

  gamma = matrix(0,nrow=100,ncol=2)
  for(i in 1:100){
    gamma[i,] = a[i,] * b[i,] / (a[i,1]*b[i,1] + a[i,2]*b[i,2])
  }

  return(list(Y=Y, X=X, gamma=gamma))
}

#'Implements the Baum-Welsh Algorithm
#' @param Y the vector of observed states
#' @return The function returns the optimized emission probabilities E, the transition probabilities P,
#' the starting probabilities pi, and a matrix of the predicted marginal probabilities
BW = function(Y, maxit=100){
  HMM_backward = function(Y, P, E){
    n = length(Y)
    b = matrix(rep(0,2*n), ncol=2)
    b[n,] = c(1,1)
    for(t in (n-1):1){
      b[t,1] = P[1,1]*E[1,Y[t+1]]*b[(t+1),1] + P[1,2]*E[2,Y[t+1]]*b[(t+1),2]
      b[t,2] = P[2,1]*E[1,Y[t+1]]*b[(t+1),1] + P[2,2]*E[2,Y[t+1]]*b[(t+1),2]
    }
    return(b)
  }
  HMM_forward = function(Y, P, E, pi){
    n = length(Y)
    a = matrix(rep(0,2*n), ncol=2)
    a[1,] = pi*E[,Y[1]]
    for(t in 1:(n-1)){
      a[(t+1),1] = E[1,Y[t+1]]*(a[t,1]*P[1,1] + a[t,2]*P[2,1])
      a[(t+1),2] = E[2,Y[t+1]]*(a[t,1]*P[1,2] + a[t,2]*P[2,2])
    }
    return(a)
  }

  n = length(Y)
  ll = numeric()
  ll[1] = 0
  #initialize P, E and pi
  temp = runif(15)
  P = new_P = matrix(c(temp[1],1-temp[1],temp[2],1-temp[2]), ncol=2, byrow=TRUE)
  E = new_E = matrix(cbind(temp[3:8]/sum(temp[3:8]), temp[9:14]/sum(temp[9:14]) ), nrow=2)
  pi = new_pi = c(temp[15],1-temp[15])

  gamma = matrix(0,nrow=n,ncol=2)

  for(k in 1:maxit){
    b = HMM_backward(Y, P, E)
    a = HMM_forward(Y, P, E, pi)
    #pi[1]*E[1,Y[1]]*b[1,1] + pi[2]*E[2,Y[1]]*b[1,2] == sum(a[n,])

    pY = ll[k] = sum(a[n,])
    for(i in 1:100){
      gamma[i,] = a[i,] * b[i,] / (a[i,1]*b[i,1] + a[i,2]*b[i,2])
    }

    #Update E
    for(j in 1:6){
      new_E[1,j] = sum(gamma[,1]*ifelse(Y==j,1,0))/sum(gamma[,1])
      new_E[2,j] = sum(gamma[,2]*ifelse(Y==j,1,0))/sum(gamma[,2])
    }
    #Update P
    g11 = g12 = g21 = g22 = rep(0,100)
    for(t in 2:n){
      g11[t] = b[t,1]*E[1,Y[t]]*P[1,1]*a[(t-1),1]/pY
      g12[t] = b[t,2]*E[2,Y[t]]*P[1,2]*a[(t-1),1]/pY
      g21[t] = b[t,1]*E[1,Y[t]]*P[2,1]*a[(t-1),2]/pY
      g22[t] = b[t,2]*E[2,Y[t]]*P[2,2]*a[(t-1),2]/pY
    }
    new_P[1,1] = sum(g11[2:n])/sum(g11[2:n]+g12[2:n])
    new_P[1,2] = 1 - new_P[1,1]
    new_P[2,1] = sum(g21[2:n])/sum(g21[2:n]+g22[2:n])
    new_P[2,2] = 1 - new_P[2,1]

    #Update pi
    new_pi = gamma[1,]

    E = new_E
    P = new_P
    pi = new_pi
  }
  return(list(E = E, P = P, pi = as.matrix(pi,ncol=1), gamma=gamma, ll=ll))
}
