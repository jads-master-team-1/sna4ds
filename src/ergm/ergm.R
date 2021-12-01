# ERGM Modeling

## Set seed
random_seed <- 42
set.seed(random_seed)

## Set threshold
threshold <- 1.0

## Load data
mep <- read.csv("./data/MEP_data.csv", header = TRUE, sep = ",")
edges_soc_clim <- read.csv("./data/edges_for_against_social.csv",
                           header = TRUE,
                           sep = ",")
edges_ind_clim <- read.csv("./data/edges_for_against_industry.csv",
                           header = TRUE,
                           sep = ",")

# Drop columns
mep <- mep[c("MEP", "group_abbv", "country")]

# Rename columns
colnames(mep) <- c("mep", "party", "country")
colnames(edges_soc_clim) <- c("mep1", "mep2", "agreement")
colnames(edges_ind_clim) <- c("mep1", "mep2", "agreement")

# Create column
edges_soc_clim$agreement_percent <- round(edges_soc_clim$agreement / 44, digits = 2)
edges_ind_clim$agreement_percent <- round(edges_ind_clim$agreement / 16, digits = 2)

## Create network
uw_edges_soc_clim_100 <- edges_soc_clim[edges_soc_clim[, "agreement_percent"] >= threshold, c("mep1","mep2")]
uw_net_soc_clim_100 <- igraph::graph_from_data_frame(uw_edges_soc_clim_100,
                                                     vertices = mep,
                                                     directed = FALSE)

## Convert network
uw_net_soc_clim_100_nw <- intergraph::asNetwork(uw_net_soc_clim_100)

## Check terms

### Isolates
summary(uw_net_soc_clim_100_nw ~ isolates)

### Degree
summary(uw_net_soc_clim_100_nw ~ degree(0:120))

### Triangles
summary(uw_net_soc_clim_100_nw ~ triangles)

### Triad Census
summary(uw_net_soc_clim_100_nw ~ triadcensus)

### KStar
summary(uw_net_soc_clim_100_nw ~ kstar(0:120))

## Create ERGM

### Edges (1sec - CONVERGED)
m0 <- ergm::ergm(uw_net_soc_clim_100_nw ~ edges)

sink("./models/m0.txt")

#### Model summary
summary(m0)

#### Goodness of Fit
(m0_fit <- ergm::gof(m0))
png("./models/m0_gof.png", width = 1152, height = 585, units = "px")
par(mfrow = c(2, 2))
plot(m0_fit)
par(mfrow = c(1, 1))
dev.off()

sink()

#### Save model
saveRDS(m0, file = "./models/m0.rds")
saveRDS(m0_fit, file = "./models/m0_fit.rds")

### Isolates (5sec - ERROR (number of edges exceeds that in observed))
m1 <- ergm::ergm(uw_net_soc_clim_100_nw ~ edges + isolates,
                 check.degeneracy = TRUE,
                 control = ergm::control.ergm(MCMC.burnin = 5000,
                                              MCMC.samplesize = 10000,
                                              seed = random_seed,
                                              MCMLE.maxit = 60,
                                              parallel = 4),
                 verbose = TRUE)

### Edges + Isolates + GWDegree (2min - CONVERGED)
m2_1 <- ergm::ergm(uw_net_soc_clim_100_nw ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE),
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
par(mfrow = c(2, 2))
plot(m2_1_fit)
par(mfrow = c(1, 1))
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
m2_2 <- ergm::ergm(uw_net_soc_clim_100_nw ~ edges + isolates + altkstar(0.5, fixed = TRUE),
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
par(mfrow = c(2, 2))
plot(m2_2_fit)
par(mfrow = c(1, 1))
dev.off()

#### MCMC Diagnostics
png("./models/m2_2_mcmc.png", width = 1152, height = 585, units = "px")
ergm::mcmc.diagnostics(m2_2)
dev.off()

sink()

#### Save model
saveRDS(m2_2, file = "./models/m2_2.rds")
saveRDS(m2_2_fit, file = "./models/m2_2_fit.rds")

### Edges + Isolates + GWDegree + nodematch("party") (45min - NO PROGRESS)
m3 <- ergm::ergm(uw_net_soc_clim_100_nw ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                   nodematch("party"),
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
par(mfrow = c(2, 2))
plot(m3_fit)
par(mfrow = c(1, 1))
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
m4 <- ergm::ergm(uw_net_soc_clim_100_nw ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
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
par(mfrow = c(2, 2))
plot(m4_fit)
par(mfrow = c(1, 1))
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
