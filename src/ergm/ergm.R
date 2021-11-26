# ERGM Modeling

## Set seed
random_seed <- 42
set.seed(random_seed)

## Set threshold
threshold <- 1.0

## Load data
mep <- read.csv("./data/MEP_data.csv", header = TRUE, sep = ",")
edges <- read.csv("./data/edges_for_against_climate.csv",
                  header = TRUE,
                  sep = ",")

## Add columns
edges$AgreementPercentage <- round(edges$Agreement / 44, digits = 2)

## Create network
##
## Attributes: vertex.names, country, group_abbv
##
edges <- edges[edges[, "AgreementPercentage"] >= threshold, c("MEP_1","MEP_2")]
climate_graph <- igraph::graph_from_data_frame(edges,
                                               vertices = mep,
                                               directed = FALSE)

## Convert network
climate_net <- intergraph::asNetwork(climate_graph)

## Check terms

### Isolates
summary(climate_net ~ isolates)

### Degree
summary(climate_net ~ degree(0:120))

### Triangles
summary(climate_net ~ triangles)

### Triad Census
summary(climate_net ~ triadcensus)

### KStar
summary(climate_net ~ kstar(0:120))

## Create ERGM

### Edges (1sec - CONVERGED)
m0 <- ergm::ergm(climate_net ~ edges)

sink("./models/m0.txt")

#### Model summary
summary(m0)

#### Goodness of Fit
(m0_fit <- ergm::gof(m0))
png("./models/m0_gof.png", width = 1152, height = 585, units = "px")
plot(m0_fit)
dev.off()

sink()

#### Save model
saveRDS(m0, file = "./models/m0.rds")
saveRDS(m0_fit, file = "./models/m0_fit.rds")

### Isolates (5sec - ERROR (number of edges exceeds that in observed))
m1 <- ergm::ergm(climate_net ~ edges + isolates,
                 check.degeneracy = TRUE,
                 control = ergm::control.ergm(MCMC.burnin = 5000,
                                              MCMC.samplesize = 10000,
                                              seed = random_seed,
                                              MCMLE.maxit = 60,
                                              parallel = 4),
                 verbose = TRUE)

### Edges + Isolates + GWDegree (2min - CONVERGED)
m2_1 <- ergm::ergm(climate_net ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE),
                   check.degeneracy = TRUE,
                   control = ergm::control.ergm(MCMC.burnin = 5000,
                                                MCMC.samplesize = 10000,
                                                seed = random_seed,
                                                MCMLE.maxit = 60,
                                                parallel = 4),
                   verbose = TRUE)

sink("./models/m2_1.txt")

#### Model summary
summary(m2_1)

#### Goodness of Fit
(m2_1_fit <- ergm::gof(m2_1))
png("./models/m2_1_gof.png", width = 1152, height = 585, units = "px")
plot(m2_1_fit)
dev.off()

#### MCMC Diagnostics
png("./models/m2_1_mcmc.png", width = 1152, height = 585, units = "px")
ergm::mcmc.diagnostics(m2_1)
dev.off()

sink()

#### Save model
saveRDS(m2_1, file = "./models/m2_1.rds")
saveRDS(m2_1_fit, file = "./models/m2_1_fit.rds")

### Edges + Isolates + AltKStar (10min - NOT CONVERGED (estimation did not converge))
m2_2 <- ergm::ergm(climate_net ~ edges + isolates + altkstar(0.5, fixed = TRUE),
                   check.degeneracy = TRUE,
                   control = ergm::control.ergm(MCMC.burnin = 5000,
                                                MCMC.samplesize = 10000,
                                                seed = random_seed,
                                                MCMLE.maxit = 60,
                                                parallel = 4),
                   verbose = TRUE)

sink("./models/m2_2.txt")

#### Model summary
summary(m2_2)

#### Goodness of Fit
(m2_2_fit <- ergm::gof(m2_2))
png("./models/m2_2_gof.png", width = 1152, height = 585, units = "px")
plot(m2_2_fit)
dev.off()

#### MCMC Diagnostics
png("./models/m2_2_mcmc.png", width = 1152, height = 585, units = "px")
ergm::mcmc.diagnostics(m2_2)
dev.off()

sink()

#### Save model
saveRDS(m2_2, file = "./models/m2_2.rds")
saveRDS(m2_2_fit, file = "./models/m2_2_fit.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") (45min - NO PROGRESS)
m3 <- ergm::ergm(climate_net ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                   nodematch("group_abbv"),
                 check.degeneracy = TRUE,
                 control = ergm::control.ergm(MCMC.burnin = 5000,
                                              MCMC.samplesize = 10000,
                                              seed = random_seed,
                                              MCMLE.maxit = 60,
                                              parallel = 4),
                 verbose = TRUE)

sink("./models/m3.txt")

#### Model summary
summary(m3)

#### Goodness of Fit
(m3_fit <- ergm::gof(m3))
png("./models/m3_gof.png", width = 1152, height = 585, units = "px")
plot(m3_fit)
dev.off()

#### MCMC Diagnostics
png("./models/m3_mcmc.png", width = 1152, height = 585, units = "px")
ergm::mcmc.diagnostics(m3)
dev.off()

sink()

#### Save model
saveRDS(m3, file = "./models/m3.rds")
saveRDS(m3_fit, file = "./models/m3_fit.rds")

### Edges + Isolates + GWDegree + nodematch("country") (2min - CONVERGED)
m4 <- ergm::ergm(climate_net ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                   nodematch("country"),
                 check.degeneracy = TRUE,
                 control = ergm::control.ergm(MCMC.burnin = 5000,
                                              MCMC.samplesize = 10000,
                                              seed = random_seed,
                                              MCMLE.maxit = 60,
                                              parallel = 4),
                 verbose = TRUE)

sink("./models/m4.txt")

#### Model summary
summary(m4)

#### Goodness of Fit
(m4_fit <- ergm::gof(m4))
png("./models/m4_gof.png", width = 1152, height = 585, units = "px")
plot(m4_fit)
dev.off()

#### MCMC Diagnostics
png("./models/m4_mcmc.png", width = 1152, height = 585, units = "px")
ergm::mcmc.diagnostics(m4)
dev.off()

sink()

#### Save model
saveRDS(m4, file = "./models/m4.rds")
saveRDS(m4_fit, file = "./models/m4_fit.rds")

## Compare models
texreg::screenreg(list(m0, m2_1, m2_2, m4))
