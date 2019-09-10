library(MASS)
library(vegan)
library(adaptiveGPCA)

### Need to analyze total data of both the rarefied and un-rarefied data

set.seed(4235421)
dada2.old.bac.pca = ordinate(dada2.old.bac.pseq, "MDS", "bray")
plot_ordination(dada2.old.bac.pseq, dada2.old.bac.pca, color = "Tissue") + geom_point(size = 3)
dada2.old.fun.pca = ordinate(dada2.old.fun.pseq, "MDS", "bray")
plot_ordination(dada2.old.fun.pseq, dada2.old.fun.pca, color = "Tissue") + geom_point(size = 3)

tiff(filename = "BacPCA.tif", units = "px", width = 2000, height = 1188, res = 360)
dada2.new.bac.pca = ordinate(dada2.new.bac.pseq, "MDS", "bray")
plot_ordination(dada2.new.bac.pseq, dada2.new.bac.pca, color = "Tissue") + geom_point(size = 3)
dev.off()

tiff(filename = "FunPCA.tif", units = "px", width = 2000, height = 1188, res = 360)
dada2.new.fun.pca = ordinate(dada2.new.fun.pseq, "MDS", "bray")
plot_ordination(dada2.new.fun.pseq, dada2.new.fun.pca, color = "Tissue") + geom_point(size = 3)
dev.off()

