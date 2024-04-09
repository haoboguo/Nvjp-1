library(bio3d)

##seven categories of seocndary structure elements
##they can be measured by the program stride
##and are also used by VMD
all <- c("H", "G", "I", "E", "B", "T", "C")
for (i in 1:1000){
  pdb <- read.pdb(paste(i, ".nvjp1.pdb", sep=""))
  seq <- pdbseq(pdb)
  sse <- stride(pdb, exefile="~/stride/stride", resno=TRUE)
  h <- length(which(sse$sse %in% "H"))
  g <- length(which(sse$sse %in% "G"))
  i <- length(which(sse$sse %in% "I"))
  e <- length(which(sse$sse %in% "E"))
  b  <- length(which(sse$sse %in% "B"))
  t  <- length(which(sse$sse %in% "T"))
  c  <- length(which(sse$sse %in% "C"))

  all <- rbind(all, c(h,g,i,e,b,t,c))
}

write.table(all, file="stride.models.csv", sep=",", row.names=F, col.names=F, quote=F)

##now read the data and generate the box plot
##mean+/-sd will be calculated
sse <- read.csv("stride.models.csv", header=T)

png("sse.test.png", width=6, height=4, res=600, units="in")
boxplot(sse/381, lax=2, ylim=c(0,1),
        xlab="Secondary Structure", ylab="Residue Ratio", col=1:7)
dev.off()

mean <- c()
sd <- c()
for (i in 1:7) {
mean <- c(mean, mean(sse[,i]))
sd   <- c(sd, sd(sse[,i]))
}

mean
sd

##this R script was used for the data analysis and figures
##in Guo et al. PLoS One 2024
