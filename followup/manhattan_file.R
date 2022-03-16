rm(list = ls())
dir_data = "~/xuanyao_llw/DGN_no_filter_on_mappability/p/"
file_signal = "/project2/xuanyao/llw/DGN_no_filter_on_mappability/FDR/signals.chr.module.perm10.fdr10.txt"
file_indep_signal = "~/xuanyao_llw/DGN_no_filter_on_mappability/postanalysis/indep.signals.chr.module.perm10.fdr10.txt"
file_annotation = "/project2/xuanyao/data/mappability/gencode.v19.annotation.table.txt"
#file_signals_eQTL = "/project2/xuanyao/llw/DGN_PCO.lambda.01/trans_cis/signals_eQTL_FDR10.txt"
#file_signals_sQTL = "/project2/xuanyao/llw/DGN_PCO.lambda.01/trans_cis/signals_sQTL_FDR10.txt"

# prepare gwas data
{

  chr_seq = 1:22


  pvalue_pc = NULL;
  chr_pc = NULL; bp = NULL; snp = NULL;
  for(chr in chr_seq){
    tmp_pc = readRDS(file = paste0(dir_data, "p.module1.chr", chr, ".rds"))
    #names(tmp_pc) = paste0("C", gene_cluster_id, ":", names(tmp_pc))
    pvalue_pc = c(pvalue_pc, tmp_pc)

    chr_pc = c(chr_pc, rep(chr, length(tmp_pc)))

    bp = c(bp, sapply(strsplit(names(tmp_pc), ":"), function(x) as.numeric(x[[2]])))
    snp = c(snp, names(tmp_pc))

    cat("Chr:", chr, "\n")
  }

  gwasResults = data.frame("CHR"=chr_pc, "BP"=bp, "SNP"=snp, "P"=pvalue_pc)
  saveRDS(gwasResults, file = "manhattan.rds")

}



# prepare signals and pvalues
{
  library(data.table)
  dis_cis <- 1e+6

  PCO_trans = fread(file_signal,
                    header = FALSE,
                    col.names = c("module.snp", "pvalue", "qvalue"))
  
  PCO_trans$SNP = sapply(strsplit(PCO_trans$module.snp, ":"), function(x) paste(x[[2]], x[[3]], sep = ":"))

  PCO_trans_uniq = NULL
  for(i in unique(PCO_trans$SNP)){
    tmp = PCO_trans[PCO_trans$SNP == i, ]
    PCO_trans_uniq = rbind(PCO_trans_uniq, tmp[which.min(tmp$pvalue), ])
  }

  # prepare annotation files
  gene_pos = fread(paste0(file_annotation),
                   sep = '\t', header = T)
  gene_pos = as.data.frame(gene_pos[gene_pos$Class %in% c('protein_coding'), ]) #, 'lincRNA'


  PCO_trans_uniq$PCO = PCO_trans_uniq$SNP
  PCO_trans_uniq$elife = "no"
  PCO_trans_uniq$both = "no"
  PCO_trans_uniq$neargene = "NA"
  PCO_trans_uniq$nearestgene = "NA"
  for(i in 1:nrow(PCO_trans_uniq)){

    chr = as.numeric(strsplit(PCO_trans_uniq$SNP[i], ":")[[1]][1])
    snpos = as.numeric(strsplit(PCO_trans_uniq$SNP[i], ":")[[1]][2])

    {
      ind2 = gene_pos$Chromosome == paste0("chr", chr) & abs(gene_pos$Start - snpos) < dis_cis
      dis_snp_gene2 = gene_pos[ind2, "Start"]-snpos
      tmp2 = cbind(gene_pos[ind2, ], 'dis_snp_gene' = dis_snp_gene2)

      if(nrow(tmp2)>0){
        PCO_trans_uniq[i, "neargene"] = paste(tmp2[order(abs(tmp2$dis_snp_gene)), ][, "GeneSymbol"],
                                              collapse = ";")
        PCO_trans_uniq[i, "nearestgene"] = tmp2[order(abs(tmp2$dis_snp_gene)), ][1, "GeneSymbol"]
      }
    }
  }

  fwrite(PCO_trans_uniq, file = "PCO.trans.uniq.txt", quote = FALSE, sep = "\t")
}


# plot manhattan plot
library(qqman)
library(ggrepel)
library(dplyr)


PCO_trans_uniq[PCO_trans_uniq$pvalue==0, "pvalue"] = 1e-20 #pvalues<1e-17 is displayed as 0
gwasResults$P = 1; gwasResults[match(PCO_trans_uniq$SNP, gwasResults$SNP), "P"] = PCO_trans_uniq$pvalue
#gwasResults = readRDS("manhattan.rds")
#PCO_trans_uniq = read.table(file = "PCO.trans.uniq.txt", stringsAsFactors = FALSE, row.names = 1, header = T)
#PCO_trans_indep = read.table("PCO.trans.indep.txt", header = TRUE, stringsAsFactors = FALSE)
#PCO_trans_indep[25, 'neargene']='IKZF1'
{

  #signals_eQTL = as.character(read.table(file_signals_eQTL, stringsAsFactors = FALSE)[[1]])
  #signals_sQTL = as.character(read.table(file_signals_sQTL, stringsAsFactors = FALSE)[[1]])

  anno_dat = PCO_trans_uniq %>% group_by(nearestgene) %>% summarise(SNP_annotate_neargene = SNP[which.min(pvalue)] )
  don <- gwasResults %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( chromosome=BP+tot) %>%

    # Add highlight and annotation information
    mutate(is_draw = ifelse(SNP %in% PCO_trans_uniq$SNP, "yes", "no") ) %>%
    #mutate(is_eQTL = ifelse(SNP %in% signals_eQTL, "yes", "no") ) %>%
    #mutate(is_sQTL = ifelse(SNP %in% signals_sQTL, "yes", "no") ) %>%
    mutate(is_annotate_neargene = ifelse(SNP %in% anno_dat$SNP_annotate_neargene, "yes", "no"))

  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(chromosome) + min(chromosome) ) / 2 )
  label = PCO_trans_uniq[match(don[don$is_annotate_neargene=="yes", "SNP"], PCO_trans_uniq$SNP), ]$nearestgene
  label[label=="NA"] = ""
  
  manhattan_plot <- ggplot(don, aes(x=chromosome, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "grey"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(5, 21) ) +
    geom_point(data=subset(don, is_annotate_neargene=="yes"),
               color="indianred2", shape=19, alpha = 0.7)
  
  manhattan_plot +
    labs(x = "Chromosome", y = "-Log10(P)") +
    geom_text_repel( data=subset(don, is_annotate_neargene=="yes"),
                     aes(label=label),
                     segment.colour="grey",
                     size=3,
                     max.overlaps = 30,
                     fill = NA) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          legend.background = element_rect(color = "black", linetype = "dashed"),
          legend.key.size= unit(0.5, "cm"),
          axis.line = element_line(colour="black"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          axis.text=element_text(colour = "black", size=16),
          axis.title.y = element_text(angle=90,vjust =2, size=16),
          axis.title.x = element_text(vjust = -0.2, size=16) )
  
  manhattan_plot
  
  ggsave("manhattan_coding.png", manhattan_plot)
  system("bash ~/imgcat manhattan_coding.png")

  #ggsave("manhattan.png", manhattan_plot)
  #system("bash ~/imgcat manhattan.png")

}
saveRDS(manhattan_plot, "manhattan_plot.rds")
saveRDS(don, 'don.rds')
saveRDS(PCO_trans_uniq, "PCO_trans_uniq.rds")


########## miscellanous ###########

# prepare annotation files
{
  require(data.table)


  PCO_trans_indep = read.table(file_indep_signal, stringsAsFactors = FALSE, col.names = "SNP")
  PCO_trans_indep$PCO = PCO_trans_indep$SNP; PCO_trans_indep$elife = "no"; PCO_trans_indep$both = "no"
  #PCO_trans_indep[PCO_trans_indep$SNP=="3:56849749", c("elife", "both")] = "3:56849749"
  #PCO_trans_indep[PCO_trans_indep$SNP=="9:126985334", c("elife", "both")] = "9:126985334"

  gene_pos = fread(paste0(file_annotation),
                   sep = '\t', header = T)
  gene_pos = as.data.frame(gene_pos[gene_pos$Class %in% c('protein_coding', 'lincRNA'), ])

  PCO_trans_indep$neargene = "NA"
  for(i in PCO_trans_indep$SNP){
    chr = as.numeric(strsplit(i, ":")[[1]][1])
    snpos = as.numeric(strsplit(i, ":")[[1]][2])

    {
      ind2 = gene_pos$Chromosome == paste0("chr", chr) & abs(gene_pos$Start - snpos)<5*10^5
      dis_snp_gene2 = gene_pos[ind2, "Start"]-snpos
      tmp2 = cbind(gene_pos[ind2, ], 'dis_snp_gene' = dis_snp_gene2)

      if(nrow(tmp2)>0){
        PCO_trans_indep[PCO_trans_indep$SNP==i, "neargene"] = tmp2[order(abs(tmp2$dis_snp_gene)), ][1, "GeneSymbol"]
      }
    }
  }

  #gene_pos[match(PCO_trans_indep$neargene, gene_pos$GeneSymbol), ]

  write.table(PCO_trans_indep, "PCO.trans.indep.txt", quote = FALSE, row.names = FALSE)


}
