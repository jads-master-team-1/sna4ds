# ERGM

## Set seed
random_seed <- 42
set.seed(random_seed)

## Set threshold
threshold <- 1.0

## Load data
mep <- read.csv("../../data/MEP_data.csv", header = TRUE, sep = ",")
edges_climate <- read.csv("../../data/edges_for_against_climate.csv",
                          header = TRUE,
                          sep = ",")
edges_industry <- read.csv("../../data/edges_for_against_industry.csv",
                           header = TRUE,
                           sep = ",")

## Add columns
edges_climate$AgreementPercentage <- round(edges_climate$Agreement / 44,
                                           digits = 2)
edges_industry$AgreementPercentage <- round(edges_industry$Agreement / 16,
                                            digits = 2)

mepGermany <- mep[mep[, "country"]== 'Germany', ]
names(mepGermany)[1] <- "mep"
edges_germany_climate <- edges_climate[edges_climate$MEP_1 %in% mepGermany$mep, ]
edges_germany_climate <- edges_germany_climate[edges_germany_climate$MEP_2 %in% mepGermany$mep, ]

threshold <- 0.9

## Create networks
edges_germany_climate <- edges_germany_climate[
  edges_germany_climate[, "AgreementPercentage"] >= threshold,
  c("MEP_1","MEP_2")
]
network_germany_climate <- igraph::graph_from_data_frame(edges_germany_climate,
                                                 vertices = mep,
                                                 directed = FALSE)


## Plot networks
plot(igraph::delete.vertices(network_germany_climate,
                             which(igraph::degree(network_germany_climate) == 0)),
     vertex.size = 4,
     vertex.label = NA,
     edge.size = 0.4,
     layout = igraph::layout_with_graphopt)

plot(network_germany_climate,
     vertex.size = 4,
     vertex.label = NA,
     edge.size = 0.4,
     layout = igraph::layout_with_graphopt)


## Convert networks
network_climate_germany_nw <- intergraph::asNetwork(network_germany_climate)

## Check climate density
(cug_gden_climate <- sna::cug.test(network_climate_germany_nw,
                                   mode = "graph",
                                   FUN = sna::gden,
                                   FUN.args = list(mode = "graph"),
                                   cmode = "size",
                                   reps = 1000))
plot(cug_gden_climate)

## Create climate ERGM
##
## Attributes: vertex.names, country, group_abbv
##
summary(network_climate_germany_nw ~ degree)

### Baseline
m0_climate <- ergm::ergm(network_climate_germany_nw ~ edges)

(m0_climate_fit <- ergm::gof(m0_climate))
plot(m0_climate_fit)

### Baseline + Structural Effects
m2_climate <- ergm::ergm(network_climate_germany_nw ~ edges + isolates + gwdegree(decay = 0.5, fixed = FALSE) + gwesp(decay = 0.1, fixed = FALSE),
                         check.degeneracy = TRUE,
                         control = ergm::control.ergm(MCMC.burnin = 5000,
                                                      MCMC.samplesize = 10000,
                                                      seed = random_seed,
                                                      MCMLE.maxit = 60,
                                                      parallel = 4))

#+ gwesp(decay = 1, fixed = TRUE)

(m1_climate_fit <- ergm::gof(m1_climate))
plot(m1_climate_fit)

ergm::mcmc.diagnostics(m1_climate)

# edges + isolates = not terrible
