

################   cov_updates    ###########################

update_theta <- function(counts, omega, theta, metadata){
  
    nsamples <- nrow(counts)
    G <- ncol(counts)
    num_covars <- ncol(metadata)
    K <- ncol(omega)
    
  
   if(dim(omega)[1]  != dim(metadata)[1]){
     stop("Number of rows in the metadata must match that of the grades of membership matrix omega")
   }
  
   if(dim(counts)[1]  != dim(metadata)[1]){
    stop("Number of rows in the metadata must match the number of rows in counts matrix")
   }
  
   if(dim(theta)[2]  != dim(counts)[2]){
    stop("Number of columns in the theta matrix must match the number of columns in counts matrix")
   }
  
   if(dim(theta)[1] != nsamples*K){
     stop("Number of rows in theta must equal number of samples times the number of clusters")
   }
  
    cluster_fac <- rep(1:K, each = nsamples)
    metadata2 <- do.call(rbind, replicate(3, metadata, simplify=FALSE))
    
    library_size <- c(replicate(K, rowSums(counts)))
    tmp1 <- c(omega_2)
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
    
    cl <- makeCluster(parallel::detectCores(),type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
    print(cl)
    
    dim(metadata2)
    metadata_cluster <- model.matrix(~factor(cluster_fac)-1)
    metadata3 <- cbind(metadata2, metadata_cluster)
    if(!is.null(colnames(metadata))){
      colnames(metadata3) <- c(colnames(metadata), paste0("clus:", 1:K))
    }else{
      colnames(metadata3) <- c(paste0("meta:", 1:dim(metadata)[2]), paste0("clus:", 1:K))
    }
    
    fits <- distrom::dmr(cl, metadata3, latent_counts, verb=1)
    stopCluster(cl)
    
    coef_fitted <- coef(fits)
    metadata4 <- cbind(rep(1, dim(counts_repped)[1]), metadata3)
    
    fitted_val <- exp(metadata4 %*% coef_fitted)
    fitted_val_2 <- t(apply(fitted_val, 1, function(x) return(x/sum(x))))
    new_theta <- fitted_val_2
    
    ll <- list("new_theta" =  new_theta,
               "distrom_model" = fits,
               "coefficients" = coef_fitted)
    
    return(ll)
}


# library(ecostructure)
# data <- get(load(system.file("extdata", "HimalayanBirdsData.rda", package = "ecostructure")))
# nsamples <- 38
# K <- 3
# num_covars <- 2
# G <- 304
# omega <- gtools::rdirichlet(nsamples, alpha = rep(1/K, K))
# theta <- gtools::rdirichlet(nsamples*K, alpha = rep(1/(K*G), G))
# metadata <- matrix(sample(1:100, nsamples*num_covars, replace=TRUE), nsamples, num_covars)
# library(distrom)
# out <- update_theta(counts, omega, theta, metadata)
# 

