pp4Thre = 0.75

file_coloc_leadSNP <- '/project2/xuanyao/llw/coloc/cis/res_coloc_leadSNP.rds'
file_resColoc_e = file.path("/project2/xuanyao/llw/coloc/cis/cis_e/data/resColoc.txt.gz")
file_resColoc_s = file.path("/project2/xuanyao/llw/coloc/cis/cis_s/data/resColoc.txt.gz")

coloc_leadSNP <- readRDS(file_coloc_leadSNP)
resColoc_e <- fread(file_resColoc_e)
resColoc_s <- fread(file_resColoc_s)

distinct(resColoc_e, Region)
distinct(resColoc_s, Region)


resColoc_e <- resColoc_e %>%
  group_by(Region) %>%
  mutate("gene" = paste(gene, collapse = ";")) %>%
  ungroup() %>%
  distinct(Region, .keep_all = TRUE) %>%
  mutate("if_e" = PP.H4.abf > pp4Thre)

resColoc_s <- resColoc_s %>%
  group_by(Region) %>%
  mutate("gene" = paste(gene, collapse = ";")) %>%
  ungroup() %>%
  distinct(Region, .keep_all = TRUE) %>%
  mutate("if_s" = PP.H4.abf > pp4Thre)

resColoc_se <- full_join(
  resColoc_e,
  resColoc_s,
  by = "Region", suffix = c("_e", "_s")
)


coloc_e <- filter(resColoc_e, if_e)
coloc_s <- filter(resColoc_s, if_s)
coloc_se <- resColoc_se %>%
  filter(if_e | if_s) %>%
  replace_na(
    list(if_e = FALSE, if_s = FALSE)
  )


count(coloc_se, if_e, if_s)


coloc_se %>% filter(!if_e & if_s) %>% View()
coloc_se %>% filter(if_e & !if_s) %>% View()
coloc_se %>% filter(if_e | if_s) %>% View()


saveRDS(
  list(
    "coloc" = distinct(coloc_se, Region) %>% pull(),
    "cis_e_cand" = distinct(coloc_e, Region) %>% pull(),
    "cis_s_cand" = distinct(coloc_s, Region) %>% pull(),
    "trans_all" = distinct(coloc_leadSNP, Region) %>% pull()
  ),
  "/project2/xuanyao/llw/coloc/cis/coloc_upset_set.rds"
)
