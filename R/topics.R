##### Estimation for Topic Models ######

## intended main function; provides defaults and selects K via marginal lhd
topics <- function(counts,
                   K,
                   shape=NULL,
                   initopics=NULL,
                   tol=0.1,
                   bf=FALSE,
                   kill=2,
                   ord=TRUE,
                   verb=1,
                   admix=TRUE,
                   nbundles=1,
                   use_squarem=FALSE,
                   init.adapt=FALSE,
                   type="full",
                   ind_model_indices = NULL,
                   signatures=NULL,
                   light=1,
                   method_admix=1,
                   sample_init=TRUE,
                   tmax=10000,...)
  ## tpxselect defaults: tmax=10000, wtol=10^(-4), qn=100, grp=NULL,
  ## nonzero=FALSE, dcut=-10, top_genes=100, burn_in=5
{
  X <- CheckCounts(counts)
  p <- ncol(X)
  if(verb>0)
    cat(sprintf("\nEstimating on a %d document collection.\n", nrow(X)))

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }

  if(type == "independent"  &&  is.null(signatures)){
    stop("For an independent model, there has to be a grouping factor data")
  }

  if(type=="independent" && is.null(ind_model_indices)){
    ind_model_indices <- dim(counts)[2]
  }

  ## check the list of candidate K values
  if(prod(K>1)!=1){ stop(cat("use K values > 1\n")) }
  K <- sort(K)

  index_init <- 1:(max(2, min(ceiling(nrow(X)*.05),100)));
  if(sample_init==TRUE){
    samp_length <- length(index_init);
    index_init <- sample(1:nrow(X),samp_length);
  }

  ## initialize
  if(init.adapt==FALSE){

  initopics <- tpxinit(X[index_init,], initopics, K[1],
                       shape, verb, nbundles=1, use_squarem=FALSE, init.adapt)
    #initopics <- t(gtools::rdirichlet(4, rep(1+ 1/K*p, p)))
  }else{
 #   if(change_start_points){
 #      initopics <- tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1]+3,
 #                          shape, verb, nbundles=1, use_squarem=FALSE, init.adapt)
 #      initopics <- initopics[,sort(sample(1:(K[1]+2), K[1], replace=FALSE))];
 #   }else{
      initopics <- tpxinit(X[index_init,], initopics, K[1],
                           shape, verb, nbundles=1, use_squarem=FALSE,
                           init.adapt)
 #    }
  }

  initopics[initopics==1] <- 1 - 1e-14;
  initopics[initopics==0] <- 1e-14;
  initopics <- normalizetpx(initopics, byrow = FALSE)

  if(type=="independent"){
     out <-  tpxThetaGroupInd(initopics, ind_model_indices, signatures)
     initopics <-out$theta;
  }

 # initopics <- initopics[,sort(sample(1:(K[1]+2), K[1], replace=FALSE))];
 # initopics <- initopics[,1:K[1]];
  ## either search for marginal MAP K and return bayes factors, or just fit
  tpx <- tpxSelect(X, K, bf, initopics,
                   alpha=shape, tol, kill, verb, nbundles, use_squarem,
                   type, ind_model_indices, signatures, light, tmax, admix=TRUE,
                   method_admix=1,sample_init=TRUE, grp=NULL, wtol=10^{-4}, qn=100,
                   nonzero=FALSE, dcut=-10,
                   top_genes=150, burn_in=5)

  K <- tpx$K

  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }

  ## topic object
  out <- list(K=K, theta=theta, omega=omega, BF=tpx$BF, D=tpx$D, X=X)
  class(out) <- "topics"
  invisible(out) }
