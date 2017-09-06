

 library(ecostructure)
 library(covtpx)
 data <- get(load(system.file("extdata", "HimalayanBirdsData.rda", package = "ecostructure")))
 taxonomic_counts <- t(exprs(data))
 
 nsamples <- 38
 K <- 3
 num_covars <- 2
 G <- 304
 omega <- gtools::rdirichlet(nsamples, alpha = rep(1/K, K))
 theta <- gtools::rdirichlet(nsamples*K, alpha = rep(1/(K*G), G))
 metadata <- matrix(sample(1:100, nsamples*num_covars, replace=TRUE), nsamples, num_covars)
 library(distrom)
 out <- update_theta(taxonomic_counts, omega, theta, metadata)
