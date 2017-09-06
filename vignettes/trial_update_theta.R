

 library(ecostructure)
 library(covtpx)
 data <- get(load(system.file("extdata", "HimalayanBirdsData.rda", package = "ecostructure")))
 taxonomic_counts <- t(exprs(data))
 
 grid_metadata <- pData(phenoData(data))
 head(grid_metadata)
 
 nsamples <- 38
 K <- 3
 num_covars <- 2
 G <- 304
 omega <- gtools::rdirichlet(nsamples, alpha = rep(1/K, K))
 theta <- gtools::rdirichlet(nsamples*K, alpha = rep(1/(K*G), G))
 metadata <- grid_metadata[,1:3]
 library(distrom)
 counts <- taxonomic_counts
 out <- update_theta(taxonomic_counts, omega, theta, metadata)

 
 plot(out$coefficients[3,])
 