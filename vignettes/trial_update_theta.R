

 library(ecostructure)
 library(covtpx)
 data <- get(load(system.file("extdata", "HimalayanBirdsData.rda", package = "ecostructure")))
 taxonomic_counts <- t(exprs(data))
 
 grid_metadata <- pData(phenoData(data))
 head(grid_metadata)
 
 nsamples <- 38
 K <- 5
 num_covars <- 2
 G <- 304
# omega <- gtools::rdirichlet(nsamples, alpha = rep(1/K, K))
# theta <- gtools::rdirichlet(nsamples*K, alpha = rep(1/(K*G), G))
 metadata <- grid_metadata[,1:3]
 library(distrom)
 counts <- taxonomic_counts
 
 iter <- 1
 maxiter <- 50
 
 
 gom_fit <- CountClust::FitGoM(taxonomic_counts, K=5, tol=0.1)
 
 theta_pre <- gom_fit$clust_5$theta
 omega_pre <- gom_fit$clust_5$omega
 
 theta_old <- matrix(0, nsamples*K, dim(theta_pre)[1])
 for(m in 1:dim(theta_pre)[2]){
   theta_old[(nsamples*(m-1)+1):(nsamples*(m)),] <- do.call(rbind, replicate(nsamples, theta_pre[,m], simplify=FALSE))
 }
 
 omega <- omega_pre
 theta <- theta_old
 
 for(iter in 1:maxiter){
   out <- update_params(taxonomic_counts, omega, theta, metadata)
   theta <- out$theta
   omega <- out$omega
   cat("We are at iteration, ", iter)
 }
 
 elevation_metadata=grid_metadata$Elevation;
 east_west_dir = grid_metadata$WorE;
 rownames(omega) <- 1:nrow(omega)
 
 BlockStructure(omega, blocker_metadata = east_west_dir,
                order_metadata = elevation_metadata,
                yaxis_label = "Elevation",
                levels_decreasing = FALSE)
 
 BlockStructure(omega_pre, blocker_metadata = east_west_dir,
                order_metadata = elevation_metadata,
                yaxis_label = "Elevation",
                levels_decreasing = FALSE)
 
 
 
 metadata <- matrix(sample(100, nsamples*num_covars, replace = FALSE), nsamples, num_covars)
 
 omega <- omega_pre
 theta <- theta_old
 
 for(iter in 1:maxiter){
   out <- update_params(taxonomic_counts, omega, theta, metadata)
   theta <- out$theta
   omega <- out$omega
   cat("We are at iteration, ", iter)
 }
 
 elevation_metadata=grid_metadata$Elevation;
 east_west_dir = grid_metadata$WorE;
 rownames(omega) <- 1:nrow(omega)
 
 BlockStructure(omega, blocker_metadata = east_west_dir,
                order_metadata = elevation_metadata,
                yaxis_label = "Elevation",
                levels_decreasing = FALSE)
 
 BlockStructure(omega_pre, blocker_metadata = east_west_dir,
                order_metadata = elevation_metadata,
                yaxis_label = "Elevation",
                levels_decreasing = FALSE)
 
 
 
 omega <- omega_pre
 theta <- theta_old
 
 for(iter in 1:maxiter){
   out <- update_params(taxonomic_counts, omega, theta, metadata=NULL)
   theta <- out$theta
   omega <- out$omega
   cat("We are at iteration, ", iter)
 }
 
 elevation_metadata=grid_metadata$Elevation;
 east_west_dir = grid_metadata$WorE;
 rownames(omega) <- 1:nrow(omega)
 
 BlockStructure(omega, blocker_metadata = east_west_dir,
                order_metadata = elevation_metadata,
                yaxis_label = "Elevation",
                levels_decreasing = FALSE)
 
 BlockStructure(omega_pre, blocker_metadata = east_west_dir,
                order_metadata = elevation_metadata,
                yaxis_label = "Elevation",
                levels_decreasing = FALSE)
 
 
 
 plot(out$coefficients[3,])
 