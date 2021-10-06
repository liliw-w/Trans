rm(list = ls())
require(data.table)

## Input files
file_string_network = "/scratch/midway2/liliw1/MODULES/STRING/9606.protein.links.v11.5.txt.gz"
file_string_info = "/scratch/midway2/liliw1/MODULES/STRING/9606.protein.info.v11.5.txt.gz"
file_meta = "/project2/xuanyao/llw/DGN_PCO.lambda.01_real/result/gene.meta.txt"

## read file
string_network = fread(file_string_network, header = TRUE)
string_info = fread(file_string_info, header = TRUE)
meta = fread(file_meta, header = TRUE)
ind_remain = !(meta$RemovePseu | meta$RemoveAuto)

## df of string protein
df_string_protein = data.frame("9606.protein" = unique(c(string_network$protein1, string_network$protein2)),
                               check.names = FALSE, stringsAsFactors = FALSE) # 19385
df_string_protein$protein = sapply(df_string_protein$"9606.protein",
                                   function(x) strsplit(x, ".", fixed = TRUE)[[1]][2])
df_string_protein$preferred_name = string_info[match(df_string_protein$"9606.protein", string_info$"#string_protein_id"), preferred_name]

## df of genes
df_meta = data.frame("gene" = meta[ind_remain, gene],
                     "GeneNameConv" = meta[ind_remain, GeneNameConv],
                     check.names = FALSE, stringsAsFactors = FALSE)

## df of string interaction
df_interaction = string_network
df_interaction$node1 = string_info[match(df_interaction$protein1, string_info$"#string_protein_id"), preferred_name]
df_interaction$node2 = string_info[match(df_interaction$protein2, string_info$"#string_protein_id"), preferred_name]
df_interaction$node1_accession = sapply(df_interaction$protein1, function(x) strsplit(x, ".", fixed = TRUE)[[1]][2])
df_interaction$node2_accession = sapply(df_interaction$protein2, function(x) strsplit(x, ".", fixed = TRUE)[[1]][2])
df_interaction$node1_annotation = string_info[match(df_interaction$protein1, string_info$"#string_protein_id"), annotation]
df_interaction$node2_annotation = string_info[match(df_interaction$protein2, string_info$"#string_protein_id"), annotation]
df_interaction$score = df_interaction$combined_score/1000
df_interaction$if_included = df_interaction$node1 %in% df_meta$gene & df_interaction$node2 %in% df_meta$gene
df_interaction$one_way = !duplicated(apply(df_interaction[, c("node1", "node2")], 1, function(x) paste(sort(x), collapse = ":")))
df_interaction$one_way2 = !duplicated(apply(df_interaction[, c("protein1", "protein2")], 1, function(x) paste(sort(x), collapse = ":")))


## output
cat("Number of genes in the data included in STRING network: ", sum(df_meta$gene %in% df_string_protein$preferred_name),
    ". Out of ", nrow(df_meta), " genes. \n")
fwrite(df_interaction, "string_interaction_2way.txt.gz",
       quote = FALSE, sep = "\t")
fwrite(df_interaction[df_interaction$if_included, ], "string_interaction_2way_included.txt.gz",
       quote = FALSE, sep = "\t")



#####################################################################
######### PPI network clustering using hierarchical cluster #########
#####################################################################
rm(list = ls())
library(data.table)
library(dynamicTreeCut)

file_mat_interac = "mat_interac.rds"
file_module = "module.rds"


### read interaction matrix data
if(is.null(file_mat_interac)){
  file_df_interaction = "/scratch/midway2/liliw1/MODULES/STRING/string_interaction_2way_included.txt.gz"
  df_interaction=fread(file_df_interaction)
  mat_interac = dcast(df_interaction[, c(4,5,10)], node1 ~ node2, value.var = "score", fun.aggregate=mean, fill = 0)
  mat_interac = as.matrix(mat_interac, rownames = "node1")
  diag(mat_interac) = 1
}else{
  mat_interac = readRDS(file_mat_interac)
}


### distance matrix and construct the tree
diss_interac = 1-mat_interac
geneTree = hclust(as.dist(diss_interac), method = "average")

# branch cutting using dynamic tree cut, as used in WGCNA
minModuleSize = 20
deepSplit = 4
dynamicMods = cutreeDynamic(dendro = geneTree, distM = diss_interac,
                            pamRespectsDendro = FALSE)
#dynamicColors = labels2colors(dynamicMods)

str(dynamicMods)
table(dynamicMods)

### result
cat("As a result,", nrow(mat_interac), "genes consist of",
    max(dynamicMods), "modules, with module size ranging from",
    table(dynamicMods)["1"], "~",
    table(dynamicMods)[as.character(max(dynamicMods))], ". \n")

names(dynamicMods) = colnames(mat_interac)
saveRDS(dynamicMods, file_module)




#####################################################################
######### hierarchical cluster write out to further cluster #########
#####################################################################
rm(list = ls())
library(data.table)

## Input files
file_mat_interac = "mat_interac.rds"
file_module = "module.rds"
module_seq = c(194,193,192)

## read file
mat_interac = readRDS(file_mat_interac)
dynamicMods = readRDS(file_module)

gene_in_mod = names(dynamicMods)[dynamicMods %in% module_seq]
mat_interac_mod = mat_interac[gene_in_mod, gene_in_mod]
mat_interac_mod[lower.tri(mat_interac_mod)] = 0
diag(mat_interac_mod) = 0
ind_interact = which(mat_interac_mod != 0, arr.ind = TRUE)


### result
network_interact = data.table("node1" = rownames(mat_interac_mod)[ind_interact[,1]],
                              "node2" = colnames(mat_interac_mod)[ind_interact[,2]],
                              "score" = (mat_interac_mod[ind_interact]),
                              "node1_Mod" = dynamicMods[rownames(mat_interac_mod)[ind_interact[,1]]],
                              "node2_Mod" = dynamicMods[colnames(mat_interac_mod)[ind_interact[,2]]])
fwrite(network_interact, paste0("network_interact_Mod", paste0(module_seq, collapse = "."), ".txt"),
       quote = FALSE, sep = "\t")
