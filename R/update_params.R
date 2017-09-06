

################   cov_updates    ###########################

update_params <- function(counts, omega, theta, metadata){
  
    nsamples <- nrow(counts)
    G <- ncol(counts)
    K <- ncol(omega)
  
   if(dim(theta)[2]  != dim(counts)[2]){
    stop("Number of columns in the theta matrix must match the number of columns in counts matrix")
   }
  
   if(dim(theta)[1] != nsamples*K){
     stop("Number of rows in theta must equal number of samples times the number of clusters")
   }
  
    cluster_fac <- rep(1:K, each = nsamples)
    
    library_size <- c(replicate(K, rowSums(counts)))
    tmp1 <- c(omega)
    tmp2 <- sweep(theta, 1, tmp1, "*")
    
    tmp3 <- tmp2
    dim(tmp3) <- c(nsamples, K, G)
    tmp3_normed <- aperm(apply(tmp3, c(1, 3), function(x) {
      if(sum(x) == 0){
        return(rep(1/length(x), length(x)))
      }else{
        return(x/sum(x))
      }}), c(2,1,3))
    
    tmp3_normed_simplified <- tmp3_normed
    dim(tmp3_normed_simplified) <- c(nsamples*K, G)
    
    counts_repped <- do.call(rbind, replicate(K, counts, simplify=FALSE))
    
    latent_counts <- counts_repped * tmp3_normed_simplified
    
    latent_counts_array <- latent_counts
    
    dim(latent_counts_array) <- c(nsamples, K, G)
    
    num_omega <- apply(latent_counts_array, c(1, 2), sum) + (1/K)
    new_omega <- t(apply(num_omega, 1, function(x) return(x/sum(x))))
    
    if(!is.null(metadata)){
      
      num_covars <- ncol(metadata)
      
      if(dim(omega)[1]  != dim(metadata)[1]){
        stop("Number of rows in the metadata must match that of the grades of membership matrix omega")
      }
      
      if(dim(counts)[1]  != dim(metadata)[1]){
        stop("Number of rows in the metadata must match the number of rows in counts matrix")
      }
      
      metadata2 <- do.call(rbind, replicate(K, metadata, simplify=FALSE))
      metadata_cluster <- model.matrix(~factor(cluster_fac)-1)
      metadata3 <- cbind(metadata2, metadata_cluster)
      
      if(!is.null(colnames(metadata))){
        colnames(metadata3) <- c(colnames(metadata), paste0("clus:", 1:K))
      }else{
        colnames(metadata3) <- c(paste0("meta:", 1:dim(metadata)[2]), paste0("clus:", 1:K))
      }
    }else{
      metadata_cluster <- model.matrix(~factor(cluster_fac)-1)
      metadata3 <- metadata_cluster
      colnames(metadata3) <- paste0("clus:", 1:K)
    }
    
    cl <- makeCluster(parallel::detectCores(),type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
    print(cl)
    suppressWarnings(fits <- distrom::dmr(cl, metadata3, latent_counts, verb=1))
    stopCluster(cl)
    
    suppressWarnings(coef_fitted <- as.matrix(coef(fits)))
    metadata4 <- cbind(rep(1, dim(counts_repped)[1]), metadata3)
    
    fitted_val <- exp(as.matrix(metadata4) %*% coef_fitted)
    fitted_val_2 <- t(apply(fitted_val, 1, function(x) return(x/sum(x))))
    new_theta <- fitted_val_2
    
    dist_omega <- mean((omega - new_omega)^2)
    dist_theta <- mean((theta - new_theta)^2)
    
    cat("distance measure: omega - ", dist_omega, ":  theta - ", dist_theta, "\n");
    
    
    ll <- list("theta" =  new_theta,
               "omega" =  new_omega,
               "coefficients" = coef_fitted)
    return(ll)
}




