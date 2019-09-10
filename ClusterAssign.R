library(vegan)

bac.data = otu_table(dada2.new.bac.001.tiss)
bac.norm = t(bac.data) / colSums(bac.data)
bac.clust = hclust(vegdist(bac.norm), method = "complete")
plot(bac.clust)
rect.hclust(bac.clust, k = 6, border = "red")

fun.data = otu_table(dada2.new.fun.001.tiss)
fun.norm = t(fun.data) / colSums(fun.data)
fun.clust = hclust(vegdist(fun.norm), method = "complete")
plot(fun.clust)
rect.hclust(fun.clust, k = 6, border = "red")

bac.group.data = cbind(abundances(dada2.new.bac.001.tiss, "compositional"), abundances(dada2.new.bac.001.geno, "compositional"), abundances(dada2.new.bac.001.year, "compositional"), abundances(dada2.new.bac.001.loc, "compositional"))
write.csv(bac.group.data, file = "bac_01_data.csv")

fun.group.data = cbind(abundances(dada2.new.fun.001.tiss, "compositional"), abundances(dada2.new.fun.001.geno, "compositional"), abundances(dada2.new.fun.001.year, "compositional"), abundances(dada2.new.fun.001.loc, "compositional"))
write.csv(fun.group.data, file = "fun_01_data.csv")
