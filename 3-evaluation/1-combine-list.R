# Identify overlapping lists of CpGs

# libraries
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)

# 1. Load age related CpGs + double check any chrX genes are removed.
anno_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
chrX_Y <- unique(c(anno_epic[anno_epic$chr%in% c("chrX", "chrY"),]$Name,
                   anno_450k[anno_450k$chr%in% c("chrX", "chrY"),]$Name))

f <- list.files("2-age-association/3-output/")
f <- f[f != "pheno_examples.Rdata"]

for (i in f){
  get(load(paste0("2-age-association/3-output/", i))) |>
           dplyr::filter(!cg %in% chrX_Y)
}

# 2. Load Senescence, PCGT, and Proliferation CpGs

# Senescence
load("1-functional-association/1-senescence/3-output/coef.Rdata")
sen <- df |> 
  dplyr::filter(!cg %in% chrX_Y) |> 
  dplyr::mutate(padj = p.adjust(p, method = "fdr")) |> 
  dplyr::filter(padj < 0.05)

# PCGT
load("1-functional-association/3-pcgt/1-output/pcgt_cpgs_tss200_fetalzero.Rdata")
pcgt <- pcgt_cpgs_tss200_fetalzero |> 
  as.data.frame() |> 
  dplyr::filter(!Name %in% chrX_Y)
nrow(pcgt) # 2577 CpGs
rm(pcgt_cpgs_tss200_fetalzero);gc()

# Proliferation
load("1-functional-association/2-proliferation/3-output/coef.Rdata")
pro <- df |> 
  dplyr::filter(!cg %in% chrX_Y) |> 
  dplyr::mutate(padj = p.adjust(p, method = "fdr")) |> 
  dplyr::filter(padj < 0.05) |> 
  # invert coefficient so direction is more intuitive
  dplyr::mutate(coef = coef*-1)

rm(df, f, i, chrX_Y, anno_450k, anno_epic);gc()

# 3. Overlaps
# 3.1 PCGT
# General
PCGTgen_tmp <- general |> 
  dplyr::filter(corr_epi > 0 & cg %in% pcgt$Name) |> # correlation needs to be positive with age for PCGT; epi + imm go in same direction in this set
  dplyr::mutate(w = 1)

PCGTgen <- PCGTgen_tmp |> 
  pull(w)
names(PCGTgen) <- PCGTgen_tmp$cg
rm(PCGTgen_tmp)

# Epithelial
PCGTepi_tmp <- epithelial |> 
  dplyr::filter(corr_epi > 0 & cg %in% pcgt$Name) |>
  dplyr::mutate(w = 1)
PCGTepi <- PCGTepi_tmp |> 
  pull(w)
names(PCGTepi) <- PCGTepi_tmp$cg
rm(PCGTepi_tmp)

# Epithelial
PCGTimm_tmp <- immune |> 
  dplyr::filter(corr_imm > 0 & cg %in% pcgt$Name) |>
  dplyr::mutate(w = 1)
PCGTimm <- PCGTimm_tmp |> 
  pull(w)
names(PCGTimm) <- PCGTimm_tmp$cg
rm(PCGTimm_tmp)

# Pure: associated with age (general) but NOT in senescence or proliferation
PCGTpure_tmp <- general  |> 
  dplyr::filter(corr_epi > 0 & cg %in% pcgt$Name & !cg %in% c(pro$cg, sen$cg)) |>
  dplyr::mutate(w = 1)
PCGTpure <- PCGTpure_tmp |> 
  pull(w)
names(PCGTpure) <- PCGTpure_tmp$cg
rm(PCGTpure_tmp)

# Non-age
PCGTnonage_tmp <- pcgt |> 
  as.data.frame() |> 
  dplyr::filter(Name %in% nonage$cg) |>
  dplyr::mutate(w = 1)
PCGTnonage <- PCGTnonage_tmp |> 
  pull(w)
names(PCGTnonage) <- PCGTnonage_tmp$Name
rm(PCGTnonage_tmp)

# 3.2. Proliferation
# proliferation values were inverted before: ->
# coefficients indicate higher proliferation
# directionality in age should be same as for DMPs

pro_up <- pro |> 
  dplyr::filter(coef > 0) # up with increasing methylation
pro_down <- pro |> 
  dplyr::filter(coef < 0) # down with increasing methylation

# General - positive (more proliferation) 
PROgen_tmp <- general |> # loss of methylation with age should be associated with an increase in ki67
  dplyr::filter((corr_epi > 0 & cg %in% pro_up$cg) | (corr_epi < 0 & cg %in% pro_down$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi<0, -1, 1))
PROgen_pos <- PROgen_tmp |> 
  pull(w)
names(PROgen_pos) <- PROgen_tmp$cg
rm(PROgen_tmp)

# General - negative (less proliferation)
PROgen_tmp <- general |> # negative association
  dplyr::filter((corr_epi < 0 & cg %in% pro_up$cg) | (corr_epi > 0 & cg %in% pro_down$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi<0, -1, 1))
PROgen_neg <- PROgen_tmp |> 
  pull(w)
names(PROgen_neg) <- PROgen_tmp$cg
rm(PROgen_tmp)

# Epithelial - positive 
PROepi_tmp <- epithelial |> 
  dplyr::filter((corr_epi > 0 & cg %in% pro_up$cg) | (corr_epi < 0 & cg %in% pro_down$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi < 0 , -1, 1))
PROepi <- PROepi_tmp |> 
  pull(w)
names(PROepi) <- PROepi_tmp$cg
PROepi_pos <- PROepi
rm(PROepi_tmp, PROepi)

# Epithelial - negative 
PROepi_tmp <- epithelial |> 
  dplyr::filter((corr_epi < 0 & cg %in% pro_up$cg) | (corr_epi > 0 & cg %in% pro_down$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi < 0 , -1, 1))
PROepi <- PROepi_tmp |> 
  pull(w)
names(PROepi) <- PROepi_tmp$cg
PROepi_neg <- PROepi
rm(PROepi_tmp, PROepi)

# Immune - positive
PROimm_tmp <- immune |> 
  dplyr::filter((corr_imm > 0 & cg %in% pro_up$cg) | (corr_imm < 0 & cg %in% pro_down$cg)) |> 
  dplyr::mutate(w = ifelse(corr_imm < 0 , -1, 1))
PROimm <- PROimm_tmp |> 
  pull(w)
names(PROimm) <- PROimm_tmp$cg
PROimm_pos <- PROimm
rm(PROimm_tmp,PROimm)

# Immune - negative
PROimm_tmp <- immune |> 
  dplyr::filter((corr_imm < 0 & cg %in% pro_up$cg) | (corr_imm > 0 & cg %in% pro_down$cg)) |> 
  dplyr::mutate(w = ifelse(corr_imm < 0 , -1, 1))
PROimm <- PROimm_tmp |> 
  pull(w)
names(PROimm) <- PROimm_tmp$cg
PROimm_neg <- PROimm
rm(PROimm_tmp, PROimm)

# Pure - positive
PROpure_tmp <- general  |> 
  dplyr::filter((corr_epi > 0 & cg %in% pro_up$cg) | (corr_epi < 0 & cg %in% pro_down$cg)) |> 
  dplyr::filter(cg %in% pro$cg & !cg %in% c(pcgt$Name, sen$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi<0, -1, 1))
PROpure <- PROpure_tmp |> 
  pull(w)
names(PROpure) <- PROpure_tmp$cg
rm(PROpure_tmp)
PROpure_pos <- PROpure

# Pure - negative
PROpure_tmp <- general  |> 
  dplyr::filter((corr_epi < 0 & cg %in% pro_up$cg) | (corr_epi > 0 & cg %in% pro_down$cg)) |> 
  dplyr::filter(cg %in% pro$cg & !cg %in% c(pcgt$Name, sen$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi<0, -1, 1))
PROpure <- PROpure_tmp |> 
  pull(w)
names(PROpure) <- PROpure_tmp$cg
PROpure_neg <- PROpure
rm(PROpure_tmp, PROpure)


# Nonage
PROnonage_tmp <- pro |> 
  dplyr::filter(cg %in% nonage$cg) |> 
  dplyr::mutate(w = ifelse(coef<0, -1, 1))
PROnonage <- PROnonage_tmp |> 
  pull(w)
names(PROnonage) <- PROnonage_tmp$cg
rm(PROnonage_tmp)

# 3.3. Senescence
senescence_up <- sen |> 
  dplyr::filter(coef>0)
senescence_down <- sen |> 
  dplyr::filter(coef<0)

# General
SENgen_tmp <- general |>
  dplyr::filter((corr_imm < 0 & cg %in% senescence_down$cg) | (corr_imm > 0 & cg %in% senescence_up$cg)) |> 
  dplyr::mutate(w = ifelse(corr_imm<0, -1, 1))

SENgen <- SENgen_tmp |> 
  pull(w)
names(SENgen) <- SENgen_tmp$cg
rm(SENgen_tmp)

# Epithelial
SENepi_tmp <- epithelial |> 
  dplyr::filter((corr_epi < 0 & cg %in% senescence_down$cg) | (corr_epi > 0 & cg %in% senescence_up$cg)) |> 
  dplyr::mutate(w = ifelse(corr_epi < 0 , -1, 1))
SENepi <- SENepi_tmp |> 
  pull(w)
names(SENepi) <- SENepi_tmp$cg
rm(SENepi_tmp)

# Immune
SENimm_tmp <- immune |> 
  dplyr::filter((corr_imm < 0 & cg %in% senescence_down$cg) | (corr_imm > 0 & cg %in% senescence_up$cg)) |> 
  dplyr::mutate(w = ifelse(corr_imm < 0 , -1, 1))
SENimm <- SENimm_tmp |> 
  pull(w)
names(SENimm) <- SENimm_tmp$cg
rm(SENimm_tmp)

# Pure
SENpure_tmp <- general  |> 
  dplyr::filter((corr_imm < 0 & cg %in% senescence_down$cg) | (corr_imm > 0 & cg %in% senescence_up$cg)) |> 
  dplyr::filter(cg %in% sen$cg & !cg %in% c(pcgt$Name, pro$cg)) |>
  dplyr::mutate(w = ifelse(corr_imm<0, -1, 1))
SENpure <- SENpure_tmp |> 
  pull(w)
names(SENpure) <- SENpure_tmp$cg
rm(SENpure_tmp)

SENnonage_tmp <- sen |> 
  dplyr::filter(cg %in% nonage$cg) |> 
  dplyr::mutate(w = ifelse(coef<0, -1, 1))
SENnonage <- SENnonage_tmp |> 
  pull(w)
names(SENnonage) <- SENnonage_tmp$cg
rm(SENnonage_tmp)

# ALL
AGE_tmp <- general |> 
  filter(!cg %in% c(pcgt$Name, pro$cg, sen$cg))|> 
  dplyr::mutate(w = ifelse(corr_epi<0, -1, 1))
AGE <- AGE_tmp |> 
  pull(w)
names(AGE) <- AGE_tmp$cg
rm(AGE_tmp)

# Epithelial
AGEepi_tmp <- epithelial |> 
  filter(!cg %in% c(pcgt$Name, pro$cg, sen$cg, names(AGE)))|> 
  dplyr::mutate(w = ifelse(corr_epi<0, -1, 1))
AGEepi <- AGEepi_tmp |> 
  pull(w)
names(AGEepi) <- AGEepi_tmp$cg
rm(AGEepi_tmp)

# Immune
AGEimm_tmp <- immune |> 
  filter(!cg %in% c(pcgt$Name, pro$cg, sen$cg, names(AGE)))|> 
  dplyr::mutate(w = ifelse(corr_imm<0, -1, 1))
AGEimm <- AGEimm_tmp |> 
  pull(w)
names(AGEimm) <- AGEimm_tmp$cg
rm(AGEimm_tmp)

saveRDS(AGE, file = "3-evaluation/1-output/AGE.Rds")
saveRDS(AGEepi, file = "3-evaluation/1-output/AGEepi.Rds")
saveRDS(AGEimm, file = "3-evaluation/1-output/AGEimm.Rds")

saveRDS(PCGTgen, file = "3-evaluation/1-output/PCGTgen.Rds")
saveRDS(PCGTpure, file = "3-evaluation/1-output/PCGTpure.Rds")
saveRDS(PCGTnonage, file = "3-evaluation/1-output/PCGTnonage.Rds")
saveRDS(PCGTepi, file = "3-evaluation/1-output/PCGTepi.Rds")
saveRDS(PCGTimm, file = "3-evaluation/1-output/PCGTimm.Rds")

saveRDS(PROnonage, file = "3-evaluation/1-output/PROnonage.Rds")
saveRDS(PROgen_neg, file = "3-evaluation/1-output/PROgen_neg.Rds")
saveRDS(PROpure_neg, file = "3-evaluation/1-output/PROpure_neg.Rds")
saveRDS(PROepi_neg, file = "3-evaluation/1-output/PROepi_neg.Rds")
saveRDS(PROimm_neg, file = "3-evaluation/1-output/PROimm_neg.Rds")
saveRDS(PROgen_pos, file = "3-evaluation/1-output/PROgen_pos.Rds")
saveRDS(PROpure_pos, file = "3-evaluation/1-output/PROpure_pos.Rds")
saveRDS(PROepi_pos, file = "3-evaluation/1-output/PROepi_pos.Rds")
saveRDS(PROimm_pos, file = "3-evaluation/1-output/PROimm_pos.Rds")

saveRDS(SENgen, file = "3-evaluation/1-output/SENgen.Rds")
saveRDS(SENpure, file = "3-evaluation/1-output/SENpure.Rds")
saveRDS(SENnonage, file = "3-evaluation/1-output/SENnonage.Rds")
saveRDS(SENepi, file = "3-evaluation/1-output/SENepi.Rds")
saveRDS(SENimm, file = "3-evaluation/1-output/SENimm.Rds")
