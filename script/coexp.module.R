#rm(list = ls())

file.ex.var.regressed = snakemake@input[['file_ex_var_regressed']]
file.coexp.module = snakemake@output[['file_coexp_module']]
minModuleSize = snakemake@params[['minModuleSize']]

library(WGCNA)


# load ex_cov_regressed
datExpr = readRDS(file.ex.var.regressed)

# Run WGCNA
### Parameter specification ###
minModuleSize = minModuleSize
if_plot_adjacency_mat_parameter_selection = F
if_plot_only_tree = F
if_plot_color_and_tree = F
if_plot_eigengene_heatmap_tree = F
if_heatmap_of_network = T


### Step1: network construction ###
# determine the paramter in adjacency function: pickSoftThreshold() OR pickHardThreshold()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower = sft$powerEstimate

# network construction
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


### Step2: module detection ###
# tree construction using hierarchical clustering based on TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

# branch cutting using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)

# eigene genes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

moduleLabels = match(mergedColors, c("grey", standardColors(length(unique(mergedColors)))))-1
names(moduleLabels) = colnames(datExpr)
tmp = factor(moduleLabels, c(0, as.numeric(names(sort(table(moduleLabels)[-1], decreasing=T)))), 1:length(unique(moduleLabels))-1 )
moduleLabels = as.numeric(levels(tmp))[tmp]; names(moduleLabels) = names(tmp)

print(table(moduleLabels))
cat("Number of modules:", max(moduleLabels), "\n")

# Save results
result = list(moduleColors = mergedColors,
              moduleLabels = moduleLabels,
              MEs = mergedMEs,
              old_moduleColors = dynamicColors,
              old_moduleLabels = dynamicMods,
              old_MEs = MEs,
              geneTree = geneTree)
saveRDS(result, file = file.coexp.module)
