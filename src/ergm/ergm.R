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

## Check climate density
(cug_gden_climate <- sna::cug.test(network_climate,
                                   mode = "graph",
                                   FUN = sna::gden,
                                   FUN.args = list(mode = "graph"),
                                   cmode = "size",
                                   reps = 1000))
plot(cug_gden_climate)

## Check industry density
(cug_gden_industry <- sna::cug.test(network_industry,
                                    mode = "graph",
                                    FUN = sna::gden,
                                    FUN.args = list(mode = "graph"),
                                    cmode = "size",
                                    reps = 1000))
plot(cug_gden_industry)

## Create climate ERGM
##
## Attributes: vertex.names, country, group_abbv
##

### Baseline
m0_climate <- ergm::ergm(network_climate ~ edges)

(m0_climate_fit <- ergm::gof(m0_climate))
plot(m0_climate_fit)

### Baseline + Structural Effects
m1_climate <- ergm::ergm(network_climate ~ edges + gwdsp(decay = 1, fixed = TRUE),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4))

(m1_climate_fit <- ergm::gof(m1_climate))
plot(m1_climate_fit)

ergm::mcmc.diagnostics(m1_climate)
