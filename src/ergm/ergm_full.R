# ERGM (Full)

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

### Triangles
summary(network_climate ~ triangles)

### KStar
summary(network_climate ~ kstar(0:250))

## Check industry terms

### Degree
summary(network_industry ~ degree(0:250))

### Triangles
summary(network_industry ~ triangles)

### KStar
summary(network_industry ~ kstar(0:250))

## Create climate ERGM

### Edges (5sec - CONVERGED)
m0_climate <- ergm::ergm(network_climate ~ edges)

sink("./models/m0_climate.txt")
summary(m0_climate)
(m0_climate_fit <- ergm::gof(m0_climate))
sink()

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

sink("./models/m1_climate.txt")
summary(m1_climate)
(m1_climate_fit <- ergm::gof(m1_climate))
sink()

plot(m1_climate_fit)

sink("./models/m1_climate.txt", append = TRUE)
ergm::mcmc.diagnostics(m1_climate)
sink()

saveRDS(m1_climate, file = "./models/m1_climate.rds")

### Edges + Isolates + GWDegree (5min - CONVERGED)
m2_1_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE),
                           check.degeneracy = TRUE,
                           control = ergm::control.ergm(MCMC.burnin = 5000,
                                                        MCMC.samplesize = 10000,
                                                        seed = random_seed,
                                                        MCMLE.maxit = 60,
                                                        parallel = 4),
                           verbose = TRUE)

sink("./models/m2_1_climate.txt")
summary(m2_1_climate)
(m2_1_climate_fit <- ergm::gof(m2_1_climate))
sink()

plot(m2_1_climate_fit)

sink("./models/m2_1_climate.txt", append = TRUE)
ergm::mcmc.diagnostics(m2_1_climate)
sink()

saveRDS(m2_1_climate, file = "./models/m2_1_climate.rds")

### Edges + Isolates + AltKStar (5min - CONVERGED)
m2_2_climate <- ergm::ergm(network_climate ~ edges + isolates + altkstar(0.5, fixed = TRUE),
                           check.degeneracy = TRUE,
                           control = ergm::control.ergm(MCMC.burnin = 5000,
                                                        MCMC.samplesize = 10000,
                                                        seed = random_seed,
                                                        MCMLE.maxit = 60,
                                                        parallel = 4),
                           verbose = TRUE)

sink("./models/m2_2_climate.txt")
summary(m2_2_climate)
(m2_2_climate_fit <- ergm::gof(m2_2_climate))
sink()

plot(m2_2_climate_fit)

sink("./models/m2_2_climate.txt", append = TRUE)
ergm::mcmc.diagnostics(m2_2_climate)
sink()

saveRDS(m2_2_climate, file = "./models/m2_2_climate.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") (45min - NO PROGRESS)
m3_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                           nodematch("group_abbv"),
                         check.degeneracy = TRUE,
                         # control = ergm::control.ergm(parallel = 4),
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

sink("./models/m3_climate.txt")
summary(m3_climate)
(m3_climate_fit <- ergm::gof(m3_climate))
sink()

plot(m3_climate_fit)

sink("./models/m3_climate.txt", append = TRUE)
ergm::mcmc.diagnostics(m3_climate)
sink()

saveRDS(m3_climate, file = "./models/m3_climate.rds")

### Edges + Isolates + GWDegree + nodematch("country") (5min - CONVERGED)
m4_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                           nodematch("country"),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

sink("./models/m4_climate.txt")
summary(m4_climate)
(m4_climate_fit <- ergm::gof(m4_climate))
sink()

plot(m4_climate_fit)

sink("./models/m4_climate.txt", append = TRUE)
ergm::mcmc.diagnostics(m4_climate)
sink()

saveRDS(m4_climate, file = "./models/m4_climate.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") + nodematch("country") (45min - NO PROGRESS)
m5_climate <- ergm::ergm(network_climate ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                           nodematch("group_abbv") + nodematch("country"),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4),
                         verbose = TRUE)

sink("./models/m5_climate.txt")
summary(m5_climate)
(m5_climate_fit <- ergm::gof(m5_climate))
sink()

plot(m5_climate_fit)

sink("./models/m5_climate.txt", append = TRUE)
ergm::mcmc.diagnostics(m5_climate)
sink()

saveRDS(m5_climate, file = "./models/m5_climate.rds")

## Create industry ERGM

### Edges (5sec - CONVERGED)
m0_industry <- ergm::ergm(network_industry ~ edges)

sink("./models/m0_industry.txt")
summary(m0_industry)
(m0_industry_fit <- ergm::gof(m0_industry))
sink()

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

sink("./models/m1_industry.txt")
summary(m1_industry)
(m1_industry_fit <- ergm::gof(m1_industry))
sink()

plot(m1_industry_fit)

sink("./models/m1_industry.txt", append = TRUE)
ergm::mcmc.diagnostics(m1_industry)
sink()

saveRDS(m1_industry, file = "./models/m1_industry.rds")

### Edges + Isolates + GWDegree (16hr 15min - NO PROGRESS)
m2_1_industry <- ergm::ergm(network_industry ~ edges + gwdegree(decay = 0.5, fixed = FALSE),
                            check.degeneracy = TRUE,
                            control = ergm::control.ergm(MCMC.burnin = 5000,
                                                         MCMC.samplesize = 10000,
                                                         seed = random_seed,
                                                         MCMLE.maxit = 60,
                                                         parallel = 4),
                            verbose = TRUE)

sink("./models/m2_1_industry.txt")
summary(m2_1_industry)
(m2_1_industry_fit <- ergm::gof(m2_1_industry))
sink()

plot(m2_1_industry_fit)

sink("./models/m2_1_industry.txt", append = TRUE)
ergm::mcmc.diagnostics(m2_1_industry)
sink()

saveRDS(m2_1_industry, file = "./models/m2_1_industry.rds")

### Edges + Isolates + AltKStar (10min - NOT CONVERGED (MAX ITER))
m2_2_industry <- ergm::ergm(network_industry ~ edges + altkstar(0.5, fixed = TRUE),
                            check.degeneracy = TRUE,
                            control = ergm::control.ergm(MCMC.burnin = 5000,
                                                         MCMC.samplesize = 10000,
                                                         seed = random_seed,
                                                         MCMLE.maxit = 60,
                                                         parallel = 4),
                            verbose = TRUE)

sink("./models/m2_2_industry.txt")
summary(m2_2_industry)
(m2_2_industry_fit <- ergm::gof(m2_2_industry))
sink()

plot(m2_2_industry_fit)

sink("./models/m2_2_industry.txt", append = TRUE)
ergm::mcmc.diagnostics(m2_2_industry)
sink()

saveRDS(m2_2_industry, file = "./models/m2_2_industry.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") (NA - NOT RUN)
m3_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                            nodematch("group_abbv"),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

sink("./models/m3_industry.txt")
summary(m3_industry)
(m3_industry_fit <- ergm::gof(m3_industry))
sink()

plot(m3_industry_fit)

sink("./models/m3_industry.txt", append = TRUE)
ergm::mcmc.diagnostics(m3_industry)
sink()

saveRDS(m3_industry, file = "./models/m3_industry.rds")

### Edges + Isolates + GWDegree + nodematch("country") (NA - NOT RUN)
m4_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                            nodematch("country"),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

sink("./models/m4_industry.txt")
summary(m4_industry)
(m4_industry_fit <- ergm::gof(m4_industry))
sink()

plot(m4_industry_fit)

sink("./models/m4_industry.txt", append = TRUE)
ergm::mcmc.diagnostics(m4_industry)
sink()

saveRDS(m4_industry, file = "./models/m4_industry.rds")

### Edges + Isolates + GWDegree + nodematch("group_abbv") + nodematch("country") (NA - NOT RUN)
m5_industry <- ergm::ergm(network_industry ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) +
                            nodematch("group_abbv") + nodematch("country"),
                          check.degeneracy = TRUE,
                          control = ergm::control.ergm(MCMC.burnin = 5000,
                                                       MCMC.samplesize = 10000,
                                                       seed = random_seed,
                                                       MCMLE.maxit = 60,
                                                       parallel = 4),
                          verbose = TRUE)

sink("./models/m5_industry.txt")
summary(m5_industry)
(m5_industry_fit <- ergm::gof(m5_industry))
sink()

plot(m5_industry_fit)

sink("./models/m5_industry.txt", append = TRUE)
ergm::mcmc.diagnostics(m5_industry)
sink()

saveRDS(m5_industry, file = "./models/m5_industry.rds")
