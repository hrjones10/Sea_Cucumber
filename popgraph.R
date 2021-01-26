### LIBRARIES
library(gstudio)
library(popgraph)
library(igraph)
library(ggplot2)


#### FIRST CONVERT GENOTYPE DATA TO GENEPOP FORMAT

# read in master data and filter
dat <- read_population("filtered_2719neutral_snps_genepop.txt", type = "genepop", header = T)
dat$Population <- substr(dat$ID, 1, 3)
dat$Population

dat.mv <- to_mv(dat)
str(dat.mv)
graph <- popgraph(x = dat.mv, groups = dat$Population)

layout <- layout.fruchterman.reingold(graph)
plot(graph, vertex.size = 8)

## draw the graph: 

V(graph)$names <- c("AK1", "AK2", "AK3", "AK4", "CAL", "CRA", "HOP", "JER", "JUA", "LAS", "LEG", "MAL", "MAZ", "OGD", "PRI", "QUA", "RBY", "REN", "SEL", "SGI", "SHE", "TBL", "TOF", "TOL")

V(graph)$names

nodes <- as.data.frame(V(graph)$names)
nodes
colnames(nodes) <- "names"

nodes$COLOR <- ifelse(nodes$names == "REN", "#d95f02", ifelse(nodes$names == "TBL" | nodes$names == "MAL" | nodes$names == "HOP" | nodes$names == "RBY" | nodes$names == "SHE" | nodes$names == "CRA" | nodes$names == "QUA" | nodes$names == "TOF" | nodes$names == "JER" | nodes$names == "LAS" | nodes$names == "SGI" | nodes$names == "OGD", "#1b9e77", "#7570b3"))

V(graph)$colors <- nodes$COLOR
sites <- read.table("24sites_lat_long.txt", header = T)
graph <- decorate_graph(g, sites, stratum = "SITE")

g <- ggplot()
g <- g + geom_edgeset(aes(x=LONG,y=LAT), graph)
g <- g + labs(x = "Longitude", y = "Latitude")
g <- g + geom_nodeset(aes(x = LONG, y = LAT, color = colors, size = size), graph)
g <- g + scale_color_manual(values = c("#1b9e77", "#7570b3", "#d95f02"))
g <- g + theme_empty()

g
