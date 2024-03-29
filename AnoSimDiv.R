library(MASS)
library(vegan)
library(phyloseq)
library(microbiome)

bac.dist = vegdist(t(as.matrix(otu_table(dada2.new.bac.rare))))
attach(as.data.frame(sample_data(dada2.new.bac.rare)))

bac.ano.tiss = anosim(bac.dist, Tissue)
bac.ano.geno = anosim(bac.dist, Genotype)
bac.ano.loc = anosim(bac.dist, Location)
bac.ano.year = anosim(bac.dist, Year)
bac.ano.txg = anosim(bac.dist, TxG)
bac.ano.txy = anosim(bac.dist, TxY)
bac.ano.txl = anosim(bac.dist, TxL)

bac.OH.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S1", "S5", "S9", "S13", "S17", "S21"), dada2.new.bac.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S1", "S5", "S9", "S13", "S17", "S21"), dada2.new.bac.rare))))
bac.ano.OH.geno = anosim(bac.OH.dist, Genotype)
bac.ano.OH.loc = anosim(bac.OH.dist, Location)
bac.ano.OH.year = anosim(bac.OH.dist, Year)

bac.OG.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S2", "S6", "S10", "S14", "S18", "S22"), dada2.new.bac.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S2", "S6", "S10", "S14", "S18", "S22"), dada2.new.bac.rare))))
bac.ano.OG.geno = anosim(bac.OG.dist, Genotype)
bac.ano.OG.loc = anosim(bac.OG.dist, Location)
bac.ano.OG.year = anosim(bac.OG.dist, Year)

bac.H.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S3", "S7", "S11", "S15", "S19", "S23"), dada2.new.bac.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S3", "S7", "S11", "S15", "S19", "S23"), dada2.new.bac.rare))))
bac.ano.H.geno = anosim(bac.H.dist, Genotype)
bac.ano.H.loc = anosim(bac.H.dist, Location)
bac.ano.H.year = anosim(bac.H.dist, Year)

bac.G.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S4", "S8", "S12", "S16", "S20", "S24"), dada2.new.bac.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S4", "S8", "S12", "S16", "S20", "S24"), dada2.new.bac.rare))))
bac.ano.G.geno = anosim(bac.G.dist, Genotype)
bac.ano.G.loc = anosim(bac.G.dist, Location)
bac.ano.G.year = anosim(bac.G.dist, Year)
  

fun.dist = vegdist(t(as.matrix(otu_table(dada2.new.fun.rare))))
attach(as.data.frame(sample_data(dada2.new.fun.rare)))

fun.ano.tiss = anosim(fun.dist, Tissue)
fun.ano.geno = anosim(fun.dist, Genotype)
fun.ano.loc = anosim(fun.dist, Location)
fun.ano.year = anosim(fun.dist, Year)
fun.ano.txg = anosim(fun.dist, TxG)
fun.ano.txy = anosim(fun.dist, TxY)
fun.ano.txl = anosim(fun.dist, TxL)

fun.OH.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S1", "S5", "S9", "S13", "S17", "S21"), dada2.new.fun.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S1", "S5", "S9", "S13", "S17", "S21"), dada2.new.fun.rare))))
fun.ano.OH.geno = anosim(fun.OH.dist, Genotype)
fun.ano.OH.loc = anosim(fun.OH.dist, Location)
fun.ano.OH.year = anosim(fun.OH.dist, Year)

fun.OG.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S2", "S6", "S10", "S14", "S18", "S22"), dada2.new.fun.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S2", "S6", "S10", "S14", "S18", "S22"), dada2.new.fun.rare))))
fun.ano.OG.geno = anosim(fun.OG.dist, Genotype)
fun.ano.OG.loc = anosim(fun.OG.dist, Location)
fun.ano.OG.year = anosim(fun.OG.dist, Year)

fun.H.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S3", "S7", "S11", "S15", "S19", "S23"), dada2.new.fun.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S3", "S7", "S11", "S15", "S19", "S23"), dada2.new.fun.rare))))
fun.ano.H.geno = anosim(fun.H.dist, Genotype)
fun.ano.H.loc = anosim(fun.H.dist, Location)
fun.ano.H.year = anosim(fun.H.dist, Year)

fun.G.dist = vegdist(t(as.matrix(otu_table(prune_samples(c("S4", "S8", "S12", "S16", "S20", "S24"), dada2.new.fun.rare)))))
attach(as.data.frame(sample_data(prune_samples(c("S4", "S8", "S12", "S16", "S20", "S24"), dada2.new.fun.rare))))
fun.ano.G.geno = anosim(fun.G.dist, Genotype)
fun.ano.G.loc = anosim(fun.G.dist, Location)
fun.ano.G.year = anosim(fun.G.dist, Year)
