# ERGM

## Set seed
random_seed <- 42
set.seed(random_seed)

## Set threshold
threshold <- 1.0

## Load data
mep <- read.csv("./data/MEP_data.csv", header = TRUE, sep = ",")
edges_climate <- read.csv("./data/edges_for_against_climate.csv",
                          header = TRUE,
                          sep = ",")
edges_industry <- read.csv("./data/edges_for_against_industry.csv",
                           header = TRUE,
                           sep = ",")

## Add columns
edges_climate$AgreementPercentage <- round(edges_climate$Agreement / 44,
                                           digits = 2)
edges_industry$AgreementPercentage <- round(edges_industry$Agreement / 16,
                                            digits = 2)

## Create networks
##
## Attributes: vertex.names, country, group_abbv
##
edges_climate <- edges_climate[
  edges_climate[, "AgreementPercentage"] >= threshold,
  c("MEP_1","MEP_2")
]
network_climate <- igraph::graph_from_data_frame(edges_climate,
                                                 vertices = mep,
                                                 directed = FALSE)

edges_industry <- edges_industry[
  edges_industry[, "AgreementPercentage"] >= threshold,
  c("MEP_1","MEP_2")
]
network_industry <- igraph::graph_from_data_frame(edges_industry,
                                                  vertices = mep,
                                                  directed = FALSE)

## Plot networks
plot(igraph::delete.vertices(network_climate,
                             which(igraph::degree(network_climate) == 0)),
     vertex.size = 4,
     vertex.label = NA,
     edge.size = 0.4,
     layout = igraph::layout_with_graphopt)
plot(igraph::delete.vertices(network_industry,
                             which(igraph::degree(network_industry) == 0)),
     vertex.size = 4,
     vertex.label = NA,
     edge.size = 0.4,
     layout = igraph::layout_with_graphopt)

## Convert networks
network_climate <- intergraph::asNetwork(network_climate)
network_industry <- intergraph::asNetwork(network_industry)

## Check climate terms

### Degree
summary(network_climate ~ degree(0:250))

## Check industry terms

### Degree
summary(network_industry ~ degree(0:250))

## Create climate ERGM

### Edges (5sec - CONVERGED)
m0_climate <- ergm::ergm(network_climate ~ edges)

summary(m0_climate)

(m0_climate_fit <- ergm::gof(m0_climate))
plot(m0_climate_fit)

saveRDS(m0_climate, file = "./models/m0_climate.rds")

### Edges + Isolates (5sec - ERROR (number of edges exceeds that in observed))
m1_climate <- ergm::ergm(network_climate ~ edges + isolates,
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

summary(m1_climate)

(m1_climate_fit <- ergm::gof(m1_climate))
plot(m1_climate_fit)

ergm::mcmc.diagnostics(m1_climate)

saveRDS(m1_climate, file = "./models/m1_climate.rds")

### Edges + Isolates + GWDegree (5min - CONVERGED)
m2_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

summary(m2_climate)

(m2_climate_fit <- ergm::gof(m2_climate))
plot(m2_climate_fit)

ergm::mcmc.diagnostics(m2_climate)

saveRDS(m2_climate, file = "./models/m2_climate.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") (??? - ???)
m3_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                         nodematch("group_abbv", diff = TRUE),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

summary(m3_climate)

(m3_industry_fit <- ergm::gof(m3_climate))
plot(m3_climate_fit)

ergm::mcmc.diagnostics(m3_climate)

saveRDS(m3_climate, file = "./models/m3_climate.rds")

### Edges + Isolates + GWDegree + nodematch("country") (??? - ???)
m4_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                         nodematch("country", diff = FALSE),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

summary(m4_climate)

(m4_industry_fit <- ergm::gof(m4_climate))
plot(m4_climate_fit)

ergm::mcmc.diagnostics(m4_climate)

saveRDS(m4_climate, file = "./models/m4_climate.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") + nodematch("country") (??? - ???)
m5_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                         nodematch("group_abbv", diff = TRUE) + nodematch("country", diff = FALSE),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

summary(m5_climate)

(m5_climate_fit <- ergm::gof(m5_climate))
plot(m5_climate_fit)

ergm::mcmc.diagnostics(m5_climate)

saveRDS(m5_climate, file = "./models/m5_climate.rds")

## Create industry ERGM

### Edges (5sec - CONVERGED)
m0_industry <- ergm::ergm(network_industry ~ edges)

summary(m0_industry)

(m0_industry_fit <- ergm::gof(m0_industry))
plot(m0_industry_fit)

saveRDS(m0_industry, file = "./models/m0_industry.rds")

### Edges + Isolates (10min - CONVERGED)
m1_industry <- ergm::ergm(network_industry ~ edges + isolates,
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

summary(m1_industry)

(m1_industry_fit <- ergm::gof(m1_industry))
plot(m1_industry_fit)

ergm::mcmc.diagnostics(m1_industry)

saveRDS(m1_industry, file = "./models/m1_industry.rds")

### Edges + Isolates + GWDegree (16hr 15min - NO PROGRESS)
m2_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

summary(m2_industry)

(m2_industry_fit <- ergm::gof(m2_industry))
plot(m2_industry_fit)

ergm::mcmc.diagnostics(m2_industry)

saveRDS(m2_industry, file = "./models/m2_industry.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") (??? - ???)
m3_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                          nodematch("group_abbv", diff = TRUE),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

summary(m3_industry)

(m3_industry_fit <- ergm::gof(m3_industry))
plot(m3_industry_fit)

ergm::mcmc.diagnostics(m3_industry)

saveRDS(m3_industry, file = "./models/m3_industry.rds")

### Edges + Isolates + GWDegree + nodematch("country") (??? - ???)
m4_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                          nodematch("country", diff = FALSE),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

summary(m4_industry)

(m4_industry_fit <- ergm::gof(m4_industry))
plot(m4_industry_fit)

ergm::mcmc.diagnostics(m4_industry)

saveRDS(m4_industry, file = "./models/m4_industry.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") + nodematch("country") (??? - ???)
m5_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                          nodematch("group_abbv", diff = TRUE) + nodematch("country", diff = FALSE),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

summary(m5_industry)

(m5_industry_fit <- ergm::gof(m5_industry))
plot(m5_industry_fit)

ergm::mcmc.diagnostics(m5_industry)

saveRDS(m5_industry, file = "./models/m5_industry.rds")
