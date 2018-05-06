species_specificarms <- read.table("species_specific_arms.dat",header=T)
species_specificctr <- read.table("species_specific_center.dat",header=T)
pdf("species_specific.pdf")
boxplot(species_specificarms$TOTAL,species_specificctr$TOTAL,main="species_specific Density BoxPlot", outline=FALSE, names=c("Arms","Center"))
