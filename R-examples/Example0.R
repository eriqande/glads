rm(list = ls()) #delete later

#Initial population structure
initial.population.size=100
n.loci=1000
n.alleles.per.locus=2

start <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus)

help(initial.struct)

