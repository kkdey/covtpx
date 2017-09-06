

##############  extra functions for the covtpx R package   ##################

nsamples <- 38
K <- 3
num_covars <- 2
G <- 304

omega <- rdirichlet(nsamples, alpha = rep(1/K, K))
theta <- rdirichlet(nsamples*K, alpha = rep(1/(K*G), G))

# theta <- cbind(replicate(nsamples, initopics[,1]), replicate(nsamples, initopics[,2]),
#                replicate(nsamples, initopics[,3]))

cluster_fac <- rep(1:K, each = nsamples)
metadata <- matrix(sample(1:100, nsamples*num_covars, replace=TRUE), nsamples, num_covars)

metadata2 <- do.call(rbind, replicate(3, metadata, simplify=FALSE))

library_size <- c(replicate(3, rowSums(counts)))
omega_2 <- c(omega)

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

counts_repped <- do.call(rbind, replicate(3, counts, simplify=FALSE))

latent_counts <- counts_repped * tmp3_normed_simplified

latent_counts_array <- latent_counts

temp <- latent_counts_array[which(cluster_fac == 1), ] + latent_counts_array[which(cluster_fac == 2), ] + latent_counts_array[which(cluster_fac == 3), ]

dim(latent_counts_array) <- c(nsamples, K, G)

library(distrom)
cl <- makeCluster(4,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)

dim(metadata2)
metadata_cluster <- model.matrix(~factor(cluster_fac)-1)
metadata3 <- cbind(metadata2, metadata_cluster)
colnames(metadata3) <-  c("meta1", "meta2", "clus1", "clus2", "clus3")

fits <- dmr(cl, metadata3, latent_counts, verb=1)
stopCluster(cl)
###  this fit should run  ##########################

########  extracting the coefficients from this fit   #####################

coef_fitted <- coef.dmr(fits)
metadata4 <- cbind(rep(1, dim(counts_repped)[1]), metadata3)

fitted_val <- exp(metadata4 %*% coef_fitted)
fitted_val_2 <- t(apply(fitted_val, 1, function(x) return(x/sum(x))))

new_theta <- fitted_val_2


