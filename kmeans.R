library(rattle.data)

# kmeans algorithm
euc_dist = function(x, y)
{
  return(t(x - y) %*% (x - y))
}

# reassign centers after each iteration
reassign_mean = function(assgn, X, K)
{
  D = dim(X)[2]
  centr = matrix(NA, nrow = K, ncol = D)
  for(k in 1:K)
  {
    ind = which(assgn == k)
    if(length(ind) > 1)
    {
      centr[k, ] = colMeans(X[ind, ])
    }
    else if(length(ind) == 1)
    {
      centr[k, ] = X[ind, ]
    }
  }
  return(centr)
}

# calculate sample covariance matrix
samp_cov = function(X)
{
  N = dim(X)[1]
  D = dim(X)[2]
  col_means = colMeans(X)
  scov = matrix(0, nrow = D, ncol = D)
  for(i in 1:N)
  {
    scov = scov + (1/N) * (X[i, ] - col_means) %*% t(X[i, ] - col_means)
  }
  return(scov)
}

# run an iteration of kmeans
kmeans_rnd = function(X, K)
{
  N = dim(X)[1]
  D = dim(X)[2]
  rnd = sample(1:N, K, replace = FALSE)
  centr = X[rnd, ]
  dist = rep(NA, K)
  assgn = rep(NA, N)
  while(TRUE)
  {
    assgn0 = assgn
    for(i in 1:N)
    {
      for(j in 1:K)
      {
        dist[j] = euc_dist(X[i, ], centr[j, ])
      }
      assgn[i] = which.min(dist)
    }
    centr = reassign_mean(assgn, X, K)
    if(identical(assgn0, assgn))
    {
      break
    }
  }
  if(sum(table(assgn0) == 1) < 1)
  {
    arr_cov = array(NA, c(D, D, K))
    prop = rep(NA, K)
    for(i in 1:K)
    {
      arr_cov[, , i] = samp_cov(X[assgn0 == i, ])
      prop[i] = sum(assgn == i)/N
    }
    return(list("centers" = centr, "sample_cov" = arr_cov, "proportions" = prop, "valid_cluster" = TRUE))
  }
  else
  {
    return(list("centers" = NULL, "sample_cov" = NULL, "proportions" = NULL, "valid_cluster" = FALSE))
  } 
}

# run 10 iterations of kmeans and stop when there each cluster has at least two points
kmeans_init = function(X, K)
{
  flag = 0
  for(i in 1:10)
  {
    km = kmeans_rnd(X, K)
    if(km$valid_cluster == TRUE)
    {
      flag = 1
      break
    }
  }
  if(flag == 1)
  {
    return(list("centers" = km$centers, "sample_cov" = km$sample_cov, "proportions" = km$proportions, "clusters_valid" = TRUE))
  }
  else
  {
    return(list("centers" = NULL, "sample_cov" = NULL, "proportions" = NULL, "clusters_valid" = FALSE))
  }
}
