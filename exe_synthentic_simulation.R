# Import some libraries and our source codes ------------------------------
source("R/MCCA.R")
source("R/generate_simulation_data.R")

library(LaplacesDemon)
library(MASS)
library(rTensor)
library(progress)
library(abind)

# set a seed if you need.
set.seed(1)

# Generate a dataset -------------------------------------------------------
# settings for generating the synthentic data
n <- 100 # sample size
p1 <- 50 # dimension of mode-1
p2 <- 50 # dimension of mode-2
num.group <- 3 # the number of groups
# generate dataset for MCCA
data <- generate.sim.data(n, p1, p2, num.group)

# specify the dimensions of latent space
rank1 <- 10
rank2 <- 10
fit.mcca <- mcca(data, ranks=c(rank1,rank2))

# the samples in latent space
data.core <- fit.mcca$core

# the approximated samples
data.est <- fit.mcca$est

# the reconstructed samples : fit.mcca$est + fit.mcca$tnsr.mean
data.reconst <- fit.mcca$tnsr.reconst

# projection matrices
proj <- fit.mcca$V

# latent covariance matrices
cov <- fit.mcca$Lambda

# reconstruction error rate (RER)
rer <- fit.mcca$reconst.rate

# end of file