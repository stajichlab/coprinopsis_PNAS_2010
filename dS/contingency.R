ft <-read.ftable("direction_summary.tab",skip=1,row.var.names="Direction",col.vars=list("Pair" = c("Adj","All")))
fisher.test(ft)
chisq.test(ft)
