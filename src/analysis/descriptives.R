edges_climate <- read.csv("../../data/edges_for_against.csv")
mep_climate <- read.csv("../../data/mepdata_climate.csv")
hist(edges_climate$Agreement, col = "lightblue", xlab = "Agreement",main='Histogram of spread of agreement between MEPs')
df <- data.frame(agreement=integer(), density=double())
for (i in 1:44) {
  edges_filter <- edges_climate[edges_climate[, "Agreement"]>i, c("MEP_1","MEP_2")]
  network_filter <- igraph::graph_from_data_frame(edges_filter, vertices=mep_climate, directed=FALSE)
  density_filter <- igraph::edge_density(network_filter)
  de <- list(agreement=i, density=density_filter)
  df = rbind(df,de, stringsAsFactors=FALSE) 
}

plot(df, type="l", col="blue", xlab="Agreement Threshold", ylab="Density", main="Changing density of graphs depending on agreement threshold")
