

##########  test covtpx   #####################


#######  we use an ecological data to test covtpx package   ###############

#install_github("kkdey/ecostructure")

library(ecostructure)

data <- get(load(system.file("extdata", "HimalayanBirdsData.rda",
                             package = "ecostructure")))

taxonomic_counts <- t(exprs(data))
grid_metadata <- pData(phenoData(data))
head(grid_metadata)

shape=NULL
initopics=NULL
tol=0.1
bf=FALSE
kill=2
ord=TRUE
verb=1
admix=TRUE
nbundles=1
use_squarem=FALSE
init.adapt=FALSE
type="full"
ind_model_indices = NULL
signatures=NULL
light=1
method_admix=1
sample_init=TRUE
tmax=10000

wtol=10^(-4)
qn=100
grp=NULL
nonzero=FALSE
dcut=-10
top_genes=100
burn_in=5

counts <- taxonomic_counts
K <- 3

theta <- initopics
