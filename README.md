# MCCA
Authors: Kohei Yoshikawa and Shuichi Kawano

# Introduction
This is an implementation of MCCA in "Multilinear Common Component Analysis via Kronecker Product Representation".

# Usage example: Demonstration of synthetic dataset
This is a simple example which shows how to use these codes.
Please confirm the example file (exe_synthentic_simulation.R).

## 1. Import some libraries and our source codes.
```R
source("R/MCCA.R")
source("R/generate_simulation_data.R")

library(LaplacesDemon)
library(MASS)
library(rTensor)
library(progress)
library(abind)
```

## 2. Generate a dataset
Set the following parameters, and generate a dataset.
  - n: Sample size 
  - p1: The dimension of mode-1
  - p2: The dimension of mode-2
  - num.group: The number of groups

```R
n <- 100 # sample size
p1 <- 50 # dimension of mode-1
p2 <- 50 # dimension of mode-2

data <- generate.sim.data(n, p1, p2, num.group)
```

## 3. Specify the dimensions of latent space
```R
rank1 <- 10
rank2 <- 10
```

## 4. Perform MCCA
Now, you can perform the MCCA.
```R
fit.mcca <- mcca(data, ranks=c(rank1,rank2))
```

# Licence
These code are free, non-commercial and open source.
