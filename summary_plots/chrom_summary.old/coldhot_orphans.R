orphanscold <- read.table("orphans_cold.dat",header=T)
orphanshot <- read.table("orphans_hot.dat",header=T)
pdf("coldhot_orphans.pdf")
boxplot(orphanscold$TOTAL,orphanshot$TOTAL,main="orphans Cold-Hot Density BoxPlot", outline=FALSE, names=c("Cold","Hot"))
