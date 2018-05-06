species_specificcold <- read.table("species_specific_cold.dat",header=T)
species_specifichot <- read.table("species_specific_hot.dat",header=T)
pdf("coldhot_species_specific.pdf")
boxplot(species_specificcold$TOTAL,species_specifichot$TOTAL,main="species_specific Cold-Hot Density BoxPlot", outline=FALSE, names=c("Cold","Hot"))
