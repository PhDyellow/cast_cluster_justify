#Takes distance matrix, runs hclust on it, then tests k for a max hubert gamma
hclust_hubert_gamma <- function(x,  method, max_k, sim_mat) {
  #create the tree
  hc <- hclust(x, method = method)
  #cut at 1 to max_k sites
  test_k = seq(1,max_k)
  ctc <- cutree(hc, k = test_k)
  #test hubert gamma at all 50

  gamma <- vapply(seq_along(test_k ), function(k, ct = ctc, s_m = sim_mat) {
    mem_mat <- make_mem_mat(ct[,k])
    return(castcluster::hubert_gamma(s_m,  mem_mat ))
  }, numeric(1))

  return(gamma )

}
library(Matrix )
          #Converts a vector of cluster membership into a matrix, where 1 indicates same cluster, 0 indicates different cluster
make_mem_mat <- function(clust_vec ) {

  n_sites <- length(clust_vec )
  clusts <- unique(clust_vec )
#Create a sparse matrix, then convert to full for hubert gamma
long_df<-purrr::map_dfr(clusts,  \(x,  cl_v = clust_vec ) {
  members <- which(clust_vec ==  x )
  return(  data.frame(tidyr::expand_grid(i =  members, j = members ) , clust =  x))
})
  mem_mat<-as.matrix(sparseMatrix(i =  long_df$i,  j = long_df$j,  x =  long_df$clust,  dims = c(n_sites, n_sites ) ))
return(mem_mat )

}
