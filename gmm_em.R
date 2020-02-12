setwd("C:/ram/s710/final")
library(rattle.data)
library(ClusterR)
library(ggplot2)


source("kmeans.R")

# multivariate normal density
mvn_dens = function(x, mu, Sigma)
{
  D = length(mu)
  Sigma = Sigma + diag(0.0000001, D)
  f_x = (1/sqrt(2 * pi))^D * sqrt(1/(det(Sigma))) * exp(-0.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
  return(f_x)
}

# update parameters estep
update_resp = function(X, matr_means, arr_cov, prop)
{
  N = dim(X)[1]
  K = length(prop)
  resp = matrix(NA, nrow = N, ncol = K)
  for(n in 1:N)
  {
    s = 0
    for(k in 1:K)
    {
      tem = prop[k] * mvn_dens(X[n, ], matr_means[k, ], arr_cov[, , k])
      resp[n, k] = tem
      s = s + tem
    }
    resp[n, ] = as.numeric((1/s)) * resp[n, ]
  }
  return(resp)
}

# update parameters mstep
update_params = function(X, resp, muk, Sigmak, pk)
{
  N = dim(X)[1]
  D = dim(X)[2]
  K = dim(resp)[2]
  for(k in 1:K)
  {
    Nk = sum(resp[, k])
    s_mean = rep(0, D)
    for(n in 1:N)
    {
      s_mean = s_mean + resp[n, k] * X[n, ]
    }
    muk[k, ] = (1/Nk) * s_mean
    s_cov = matrix(0, nrow = D, ncol = D)
    for(n in 1:N)
    {
      s_cov = s_cov + resp[n, k] * (X[n, ] - muk[k, ]) %*% t(X[n, ] - muk[k, ])
    }
    Sigmak[, , k] = (1/Nk) * s_cov
    pk[k] = Nk/N
  }
  return(list("muk" = muk, "Sigmak" = Sigmak, "pk" = pk))
}

# calculates log likelihood
log_likelihood = function(X, muk, Sigmak, pk)
{
  N = dim(X)[1]
  K = dim(muk)[1]
  s_outer = 0
  for(n in 1:N)
  {
    s_inner = 0
    for(k in 1:K)
    {
      s_inner = s_inner + pk[k] * mvn_dens(X[n, ], muk[k, ], Sigmak[, , k])
    }
    s_outer = s_outer + log(s_inner)
  }
  return(s_outer)
}

# gmm function, K >= 2
gmm_em = function(X, K)
{
  init = kmeans_init(X, K)
  if(init$clusters_valid == FALSE)
  {
    return(list("parameters" = NULL, "responsibilities" = NULL, "clusters_valid" = FALSE, "LL" = NULL))
  }
  else
  {
    N = dim(X)[1]
    D = dim(X)[2]
    muk = init$centers # initialize parameters using kmeans
    Sigmak = init$sample_cov
    pk = init$proportions
    llik0 = log_likelihood(X, muk, Sigmak, pk)
    while(TRUE)
    {
      resp = update_resp(X, muk, Sigmak, pk)
      params = update_params(X, resp, muk, Sigmak, pk)
      muk = params$muk
      Sigmak = params$Sigmak
      pk = params$pk
      llik = log_likelihood(X, muk, Sigmak, pk) # calculate log likelihood
      if(abs(llik - llik0) < 0.0001) # convergence criterion based on likelihood
      {
        break
      }
      llik0 = llik
    }
    return(list("parameters" = params, "responsibilities" = resp, "clusters_valid" = TRUE, "LL" = llik0))
  }
}

# function for model selection based on bic
gmm_model_select = function(X, upper_K)
{
  N = dim(X)[1]
  D = dim(X)[2]
  bic = c()
  for(K in 2:upper_K) # start from a minimum of two clusters
  {
    gmm = gmm_em(X, K)
    if(gmm$clusters_valid == TRUE)
    {
      free_param = (K - 1) + K*(D) + K*(D*(D + 1)/2)
      bic = c(bic, -2*gmm$LL + free_param*log(N)) # calculate bic
    }
    else
    {
      break
    }
    print(K)
  }
  return(bic)
}

# iris dataset results
set.seed(100) # set seed
X = as.matrix(iris[, 1:4])
K = 3
gmm = gmm_em(X, K)
gmm_class = apply(gmm$responsibilities, MARGIN = 1, which.max)
gmm_r = GMM(X, 3)
table("actual" = iris$Species, "gmm" = gmm_class)
table("actual" = iris$Species, "kmeans" = kmeans(X, K)$cluster)

gmm$parameters$muk
gmm_r$centroids

# wine dataset results
X = as.matrix(wine[, 2:14])
K = 3
gmm = gmm_em(X, K)
gmm_class = apply(gmm$responsibilities, MARGIN = 1, which.max)
gmm_r = GMM(X, 3)
table("actual" = wine$Type, "gmm" = gmm_class)
table("actual" = wine$Type, "kmeans" = kmeans(X, K)$cluster)

# model selection on artificial gaussian data - run at your own risk very time consuming
D = 5
K = 6
N = 2000
X = matrix(NA, nrow = N, ncol = D)
p = rdirichlet(1, rep(1, K))
Sigmak = array(NA, c(D, D, K))
muk = matrix(NA, nrow = K, ncol = D)
for(k in 1:K)
{
  Sigmak[,,k] = rWishart(1, ceiling(runif(1, 10, 20)), diag(1, D))
  muk[k,] = rnorm(D, 10, 7)
}
for(i in 1:N)
{
  z = sample(1:K, 1, prob = p)
  X[i,] = mvrnorm(1,mu=muk[z,],Sigma=Sigmak[,,z])+mvrnorm(1,mu=rep(0,D),Sigma=diag(0.5,D))
}

# plot bic values against K
vec = gmm_model_select(X, 10)
K = 2:(length(vec)+1)
ggplot(data.frame("bic"=vec,"K"=K),aes(x=K,y=bic))+geom_point()+geom_line()


