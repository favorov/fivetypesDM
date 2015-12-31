library(differential.coverage)
load(file = 'boruta.selected.Rda')
sink('boruta.selected.probes.txt')
print(as.data.frame(closest.gene.start.by.interval(noodles = boruta.selected.probes)))
sink()
sink('boruta.bin.selected.probes.txt')
print(as.data.frame(closest.gene.start.by.interval(noodles = boruta.bin.selected.probes)))
sink()

