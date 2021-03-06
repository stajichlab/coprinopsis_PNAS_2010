orthologsleft <- read.table("orthologs_left.dat",header=T)
orthologsright <- read.table("orthologs_right.dat",header=T)
orthologsctr <- read.table("orthologs_center.dat",header=T)
pdf("orthologs.pdf")
boxplot(orthologsleft$TOTAL,orthologsctr$TOTAL,orthologsright$TOTAL,main="orthologs Density BoxPlot", outline=FALSE, names=c("Left","Center","Right"))
