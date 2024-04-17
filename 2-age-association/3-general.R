# Identify shared CpGs

# libraries
library(ggplot2)
library(dplyr)

# Remove any chrX/Y CpGs
anno_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
chrX_Y <- unique(c(anno_epic[anno_epic$chr%in% c("chrX", "chrY"),]$Name,
                   anno_450k[anno_450k$chr%in% c("chrX", "chrY"),]$Name))

# load data
load("2-age-association/1-output/age_imm.Rdata")
imm <- corr.values |>
  tibble::rownames_to_column("cg") |> 
  dplyr::filter(!cg %in% chrX_Y)
imm_sig <- imm |> 
  dplyr::filter(p.adjust < 0.05)

load("2-age-association/2-output/age_epi.Rdata")
epi <- corr.values |>
  tibble::rownames_to_column("cg") |> 
  dplyr::filter(!cg %in% chrX_Y)
epi_sig <- epi |> 
  dplyr::filter(p.adjust < 0.05)

# General: significant in both and correlation above 0.2 in either
general <- imm_sig |> 
  # check they are going in the same direction
  dplyr::inner_join(epi_sig, by = "cg", suffix = c("_imm", "_epi")) |> 
  dplyr::filter((corr_epi > 0.2 & corr_imm > 0.2) | (corr_epi < -0.2 & corr_imm < -0.2))

# Epithelial only
epithelial <- epi |> 
  dplyr::inner_join(imm, by = "cg", suffix = c("_epi", "_imm")) |> 
  dplyr::filter(p.adjust_epi < 0.05 & abs(corr_epi) > 0.2 & abs(corr_imm) < 0.02) # 11410 CpGs

immune <- epi |> 
  dplyr::inner_join(imm, by = "cg", suffix = c("_epi", "_imm")) |> 
  dplyr::filter(p.adjust_imm < 0.05 & abs(corr_imm) > 0.2 & abs(corr_epi) < 0.02) # 5336 CpGs

# Non-age associated CpGs 
nonage <- imm |> 
  dplyr::inner_join(epi, by = "cg", suffix = c("_imm", "_epi")) |> 
  dplyr::filter(abs(corr_imm) < 0.02 & abs(corr_epi) < 0.02)

x <- list(general = general$cg,
          epithelial = epithelial$cg,
          immune = immune$cg,
          nonage = nonage$cg)
ggvenn::ggvenn(x)

# Examples
load("0-data/dat.Rdata")
pheno <- dat

# Load beta matrices
load("~/Dropbox/data/3c-buccal/beta_merged.Rdata")
beta_bu <- beta_merged
load("~/Dropbox/data/3c/beta_merged.Rdata")
beta_3c <- beta_merged
load("~/Dropbox/data/3c-blood/beta_merged.Rdata")
beta_bl <- beta_merged

intersect <- intersect(rownames(beta_bu),
                       intersect(rownames(beta_3c), rownames(beta_bl)))
beta <- cbind(beta_bu[intersect,],
              beta_3c[intersect,],
              beta_bl[intersect,])
beta <- beta[,pheno$basename]
identical(pheno$basename, colnames(beta))

epithelial |> arrange(p.adjust_epi) |> dplyr::slice(1) |> pull(cg) # cg11528849
immune |> arrange(p.adjust_imm) |> dplyr::slice(1) |> pull(cg) # cg20786223
g <- general |>
  group_by(cg) |> 
  dplyr::mutate(padj = mean(as.numeric(p.adjust_epi), as.numeric(p.adjust_imm))) |> 
  ungroup() |> 
  arrange(padj) |> 
  dplyr::slice(1) |> pull(cg)
g
g$padj

nonage |> dplyr::slice(1) |> pull(cg)

# cg05024939

pheno$epi_cg11528849 <- as.numeric(beta["cg11528849",])
pheno$imm_cg20786223 <- as.numeric(beta["cg20786223",])
pheno$general_cg05024939 <- as.numeric(beta["cg05024939",])
pheno$nonage_cg24040570 <- as.numeric(beta["cg24040570",])

pheno2 <- pheno |> 
  mutate(icgrp = case_when(ic<0.2 ~ "epithelial",
                           ic>0.2 ~ "immune"))

pheno2 |> 
  ggplot(aes(x = age,
             y = epi_cg11528849)) +
  geom_point() +
  facet_wrap(~icgrp) +
  ggpubr::stat_cor()

pheno2 |> 
  ggplot(aes(x = age,
             y = imm_cg20786223)) +
  geom_point() +
  facet_wrap(~icgrp) +
  ggpubr::stat_cor()

pheno2 |> 
  ggplot(aes(x = age,
             y = general_cg05024939)) +
  geom_point() +
  facet_wrap(~icgrp) +
  ggpubr::stat_cor()

pheno2 |> 
  ggplot(aes(x = age,
             y = nonage_cg24040570)) +
  geom_point() +
  facet_wrap(~icgrp) +
  ggpubr::stat_cor()

save(pheno2, file = "2-age-association/3-output/pheno_examples.Rdata")


# save outputs -----------------
save(epithelial, file = "2-age-association/3-output/epi.Rdata")
save(immune, file = "2-age-association/3-output/imm.Rdata")
general <- general |>
  group_by(cg) |> 
  dplyr::mutate(padj = mean(as.numeric(p.adjust_epi), as.numeric(p.adjust_imm))) |> 
  ungroup() |> 
  arrange(padj)
save(general, file = "2-age-association/3-output/general.Rdata")
save(nonage, file = "2-age-association/3-output/nonage.Rdata")
