###############################################################
### Plot qqplot and histogram for pvalues,#####################
### for each module separately, but with combined chr's########
###############################################################

rm(list = ls())
library(data.table)
library(ggplot2)

dir_p = "/scratch/midway2/liliw1/MODULES/STRING/p/"
Nmodule = 194
module_seq = 1:Nmodule
plot_path = "/scratch/midway2/liliw1/MODULES/STRING/plot/"
if_hist = TRUE
if_qqplot = TRUE


### plot for each module separately, but with combined chr's
for(module in module_seq){
  ### pvalue files names
  file.p = as.character(outer(module, 1:22, FUN = function(x, y) paste0("p.module", x, ".chr", y, ".rds")))
  #file.p = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("p.module", x, ".chr", y, ".rds")))

  ### read p
  p.obs = rbindlist(lapply(file.p, function(x)
  {tmp_y=readRDS(paste0(dir_p, x));
  print(x);
  if(!is.null(tmp_y)){
    as.data.table(setNames(tmp_y, paste0(strsplit(x, '.', fixed = T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)
  }
  }))
  p.obs$V2[p.obs$V2==0] = 1e-20


  ### Draw qqplot
  if(if_qqplot){
    qqplot.hist <- function(input, title = NULL){
      observed <- sort(input)
      lobs <- -(log10(observed))
      expected <- c(1:length(observed))
      lexp <- -(log10(expected / (length(expected+1))))
      df = data.frame(x = lexp, y = lobs, yy = observed)

      res = ggplot(df, aes(x=x, y=y)) + geom_point(size = 0.1) +
        geom_abline(slope = 1, intercept = 0, color="red") +
        coord_cartesian(xlim=c(0,21), ylim=c(0,21)) +
        theme(text=element_text(size=6)) +
        labs(title = title, x = "Expected (-logP)", y = "Observed (-logP)")

      return(res)
    }

    fig_qqplot = qqplot.hist(p.obs$V2, paste0("p_qqplot_M", module))
    ggsave(paste0("p_qqplot_M", module, ".png"), fig_qqplot, path = plot_path)
  }

  ### Draw histogram
  if(if_hist){
    fig_hist = ggplot(p.obs, aes(x = V2)) +
      geom_histogram(breaks = seq(0, 1, 0.01)) +
      coord_cartesian(xlim=c(0, 1), ylim=c(0, 80000)) +
      labs(title = paste0("p_hist_M", module), x = "Observed (P)") +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=5))

    ggsave(paste0("p_hist_M", module, ".png"), fig_hist, path = plot_path)
  }
}


#meta_used_snp = fread('/project2/xuanyao/llw/eQTLGen/eQTLGen.used_snp.meta.txt', header = TRUE)
#fig.all = NULL
#for (dir_p in c("eQTLGen_lessmodule_PCO.lambda1/", "eQTLGen_lessmodule_PCO.lambda.01/", "eQTLGen_moremodule_PCO.lambda1/", "eQTLGen_PCO.lambda.01/")) {
  #file.p = list.files(path = paste0(dir_p, "p"), pattern = "^p.module.*chr.*rds", full.names = TRUE)
  #p.obs = rbindlist(lapply(file.p, function(x)
  #{tmp_y=readRDS(x);
  #print(x);
  #as.data.table(setNames(tmp_y, paste0(strsplit(rev(strsplit(x, '/', fixed = T)[[1]])[1], '.', fixed=T)[[1]][2], ":", names(tmp_y))), keep.rownames=T)
  #}))
  #setnames(p.obs, c('module_snp', 'p'))
  #p.obs$module = sapply(p.obs$module_snp, function(x) strsplit(x, ":")[[1]][1])
  #p.obs$snp_ID = sapply(p.obs$module_snp, function(x) strsplit(x, ":")[[1]][2])
  #p.obs$snp = meta_used_snp[match(p.obs$snp_ID, meta_used_snp$SNP), meta]
  #p.obs$snp_chr = sapply(p.obs$snp, function(x) strsplit(x, ":")[[1]][1])
  #p.obs = as.data.frame(p.obs[order(p.obs$p), ])
  #saveRDS(p.obs, paste0(dir_p, "p.obs.all.rds"))

#  p.obs = readRDS(paste0(dir_p, "p.obs.all.rds"))

#  fig.all = rbind(fig.all, qqplot.hist(p.obs$p, dir_p))
#  cat(dir_p, "done!\n")
#}
#ggsave('qqplot.all.p.png',
#       marrangeGrob(fig.all, nrow=length(fig.all)/2, ncol=2, top = NULL),
#       height = length(fig.all)/2*(6/2), width = 6)


# Nmodule = 19 # 75 # 19
# file.p = as.character(outer(1:Nmodule, 1:22, FUN = function(x, y) paste0("~/xuanyao_llw/eQTLGen/p/p.module", x, ".chr", y, ".rds")))
#str(unique(p.obs[p.obs$V2 == 0, snp]))
#str(unique(p.obs_lessmodule[p.obs_lessmodule$V2 == 0, snp]))
#str(unique(p.obs_01[p.obs_01$V2 == 0, snp]))
