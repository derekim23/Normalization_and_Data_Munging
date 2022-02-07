#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(x1mo);
nSamples = nrow(x1mo)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(x1mo, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, x2m, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue = apply(moduleTraitPvalue, 2, p.adjust, method = "hochberg")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(x2m),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Cluster-Phenotype Relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(x2m[,'vo2muscleMass']);
names(weight) = "MuscleMass"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(x1mo, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(x1mo, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");



#######
# Define variable weight containing the weight column of datTrait
power = as.data.frame(x2m[,'vo2muscleMass']);
names(power) = "vo2power"

geneTraitSignificance2 = as.data.frame(cor(x1mo, power, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(power), sep="");
names(GSPvalue) = paste("p.GS.", names(power), sep="");
######


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#First identify which color cluster
#Then identify which belong to that cluster p-values
#What's the significance of that membership?
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#######
# p val significant for muscle mass
# positive correlation with muscle mass
# 
abs(geneModuleMembership[moduleGenes, column])
mmpval <- p.adjust(MMPvalue[moduleGenes,column],"hochberg")
names(mmpval) <- rownames(MMPvalue[moduleGenes,])
turq_mem <- mmpval[mmpval<0.05]
names(turq_mem)
turq_pheno_corr <- abs(geneTraitSignificance[rownames(geneTraitSignificance)%in%names(turq_mem),])
turq_mem_corr <- abs(geneModuleMembership[rownames(geneModuleMembership)%in%names(turq_mem),column])
names(turq_mem_corr) <- names(turq_mem)
names(sort(turq_mem_corr,decreasing=T))
write.csv(sort(turq_mem_corr,decreasing=T),"muscleMass_mo2.csv")
#names(turq_mem)[which.max(geneTraitSignificance[rownames(geneTraitSignificance)%in%names(turq_mem),])]

rownames(geneModuleMembership)[which.max(abs(geneModuleMembership[,"MMturquoise"]))]

module2 = "magenta"
column2 = match(module2, modNames)
moduleGenes2 = moduleColors==module2

abs(geneModuleMembership[moduleGenes2, column2])
mmpval2 <- p.adjust(MMPvalue[moduleGenes2,column2],"hochberg")
names(mmpval2) <- rownames(MMPvalue[moduleGenes2,])
magenta_mem <- mmpval2[mmpval2<0.05]
names(magenta_mem)
magenta_pheno_corr <- abs(geneTraitSignificance[rownames(geneTraitSignificance)%in%names(magenta_mem),])
magenta_mem_corr <- abs(geneModuleMembership[rownames(geneModuleMembership)%in%names(magenta_mem),column])
names(magenta_mem_corr) <- names(magenta_mem)
names(sort(magenta_mem_corr,decreasing=T))
write.csv(sort(magenta_mem_corr,decreasing=T),"vo2pwr_mo2.csv")
#geneTraitSignificance[rownames(geneTraitSignificance)%in%names(magenta_mem),]
######


sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


colnames(x1mo)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


colnames(x1mo)[moduleColors=="magenta"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(x1mo)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(geneInfo, file = "geneInfo.csv")

