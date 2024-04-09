## Residue-interaction network based on physical adjacency
## for Nvjp-1
## by HB Guo

library(bio3d)
library(gplots)
library(igraph)

## we write the distance matrices for individual models
## or snapshots from MD trajectories
## this example including 1000 AF2 structure models
## named i.nvjp1.pdb, i in [1, 1000]
for (t in 1:1000) { #main loop

# read the protein
pdb <- read.pdb(paste(t, ".nvjp1.pdb", sep=""))

# protein sequence
seq <- pdbseq(pdb)

# amino acid number
aan <- length(seq)

dist.matrix <- matrix(rep(0, times=(aan)*(aan)), ncol=aan, nrow=aan) #set a matrix

for (i in 1:aan) {
        cbeta <- atom.select(pdb, resno = i, noh=T)
     for (j in 1:aan) {
		     dbeta <- atom.select(pdb, resno = j, noh=T)
     dist.matrix[i,j] <- round(min(dist.xyz(pdb$xyz[cbeta$xyz], pdb$xyz[dbeta$xyz])),3)
     }
}

## make sure diagnal elements are zero (self-interaction)
diag(dist.matrix) <- 0

write.table(dist.matrix, file=paste("csvfull/", t, ".dm.csv", sep=""),
            sep=",",
            row.names=F, col.names=F, quote=F)

} #main loop

## identify differnet categories of res (for coloring, see below)
## (positively charged, negatively charged and aromatic
pos.list <- which(seq[1:aan] %in% c("K", "R"))
neg.list <- which(seq[1:aan] %in% c("D", "E"))
aro.list <- which(seq[1:aan] %in% c("H", "W", "Y", "F"))

## create a template matrix (for a weighted network)
## reading data from the distance matrices
weighted <- matrix(rep(0, times=aan*aan), ncol=aan, nrow=aan)

for (t in 1:1000) {#loop for weighed network

dm <- read.csv(paste("csvfull/", t, ".dm.csv", sep=""), header=F)

dm2 <- as.numeric(unlist(as.matrix(dm)))

## we use 3.5 Angs as the cutoff for interactions
## 1.75/D will round to 1, if the distance is shorter than cutoff
## or round to 0, when the distance is longer than cutoff
## thereby a binary adjencency network can be constructed
nm <- round(1.75/dm2[which(!is.na(dm2))])

dm.mat <- matrix(nm, nrow=aan)
diag(dm.mat) <- 0

## weighed to the weighted matrix
weighted <- weighted + dm.mat
} #end loop for weighted network

## coloring the nodes

colors <- rep("yellow", times=aan)

shapes <- rep("circle", times=aan)

colors[pos.list] <- "blue"
colors[neg.list] <- "red"
colors[aro.list] <- "black"

v.size <- rep(2, times=aan)

wn <- graph.adjacency(matrix(as.numeric(unlist(weighted)), nrow=aan),
                      mode="undirected", weighted=T, diag=F)

## coloring the edges, depending on the persistancy
E(wn)$color[E(wn)$weight >= 750] <- "red" #red for >=75%
E(wn)$color[(E(wn)$weight < 750) & (E(wn)$weight > 250)] <- "gray"
E(wn)$color[E(wn)$weight <= 250] <- "blue" #blue for <=25%


coords <- layout_in_circle(wn)

#a regular (optimized) presentation of the network
png(file="weighted.full.reg.png", width=6, height=6,
    res=600, units="in")
plot(wn, vertex.size=v.size,
     vertex.color=colors, vertex.label=NA,
     vertex.shape=shapes, edge.width=E(wn)$weight/500)
dev.off()

#a circular representation of the network
png(file="weighted.full.cir.png", width=6, height=6,
    res=600, units="in")
plot(wn, vertex.size=v.size, layout=coords,
     vertex.color=colors, vertex.label=NA,
     vertex.shape=shapes,edge.width=E(wn)$weight/500)
dev.off()


#record the weighed adjacency network
write.table(weighted, file="weighted.full.network.csv",
            sep=",",
            row.names=F, col.names=F, quote=F)

length(E(wn)$color[E(wn)$weight >= 750])
length(E(wn)$color[E(wn)$weight <= 250])
length(which((E(wn)$weight > 250) & (E(wn)$weight < 750)))
length(E(wn)$weight)

#histograms
breaks <- hist(E(wn)$weight/10, breaks=20)$breaks
colors <- rep("gray", length(breaks))
colors[breaks >= 75] <- "red"
colors[breaks <= 25] <- "blue"

png("nvjp1.weight.hist.png", height=4, width=4, res=600, units="in")
hist(E(wn)$weight/10, main="", breaks = breaks, col=colors,
     xlab="Pesistency", ylab="Interaction Number")
dev.off()

##this R script was used for the data analysis and figures
##in Guo et al. PLoS One 2024
