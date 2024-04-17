# Compile list and examples
load("3-mouse/3-output/pcgt_cpgs_tss200_fetalzero.Rdata")
load("3-mouse/2-output/age.Rdata")

# Remove any chrX/chrY cpGs
library(IlluminaMouseMethylationanno.12.v1.mm10)
anno <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10, lociNames = rownames(cors)) |> 
  data.frame() |> 
  dplyr::rename(cg = Name)
excl <- rownames(anno[anno$chr%in%c('chrMT', 'chrX', 'chrY'),])
hist(cors$cor)

age <- cors |> 
  tibble::rownames_to_column("cg") |> 
  dplyr::filter(!cg %in% excl) |> 
  dplyr::filter(pval < 0.05 & cor > 0.000008) |> 
  arrange(pval) |> 
  dplyr::left_join(dplyr::select(anno, cg, chr, pos, GeneName_NCBI, Feature_NCBI)) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(gene = stringr::str_split(GeneName_NCBI, ";", simplify = T)[,1],
                feature = stringr::str_split(Feature_NCBI, ";", simplify = T)[,1]) |> 
  dplyr::select(-c(GeneName_NCBI, Feature_NCBI)) |> 
  dplyr::rename(coef = cor) |>
  ungroup()

nonage <- cors |> 
  tibble::rownames_to_column("cg") |> 
  dplyr::filter(!cg %in% excl) |> 
  dplyr::filter(abs(cor) < 0.000008 & pval > 0.05) |> 
  dplyr::left_join(dplyr::select(anno, cg, chr, pos, GeneName_NCBI, Feature_NCBI)) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(gene = stringr::str_split(GeneName_NCBI, ";", simplify = T)[,1],
                feature = stringr::str_split(Feature_NCBI, ";", simplify = T)[,1]) |> 
  dplyr::select(-c(GeneName_NCBI, Feature_NCBI)) |> 
  dplyr::rename(coef = cor) |>
  ungroup()

mPCGTgen <- age |> 
  dplyr::filter(cg %in% pcgt_cpgs_tss200_fetalzero$Name)

mPCGTnonage <- nonage |> 
  dplyr::filter(cg %in% pcgt_cpgs_tss200_fetalzero$Name)

# Save outputs, append information
save(age, nonage, file = "3-mouse/4-output/age_nonage_annotated.Rdata")
save(mPCGTgen, mPCGTnonage, file = "3-mouse/4-output/mPCGTgen_mPCGTnonage_annotated.Rdata")

# # Compute examples
# load("~/Dropbox/data/mouse-project/beta_final.Rdata")
# load("~/Dropbox/eca/mouse-project/2-uibk/1-preprocessing/1-output/pheno.Rdata")
# 
# ind <- match(colnames(beta_final), pheno$basename)
# pheno <- droplevels(pheno[ind,])
# table(pheno$Experiment)
# 
# # age experiment subset
# pheno_age <- pheno[pheno$Experiment!='exp5',]
# beta_age <- beta_final[,pheno_age$basename]
# identical(pheno_age$basename, colnames(beta_age))
# 
# # Append some age sites
# pheno_age_tmp <- cbind(pheno_age, t(beta_age[age[1:5,]$cg,])) |> 
#   tidyr::pivot_longer(any_of(age[1:5,]$cg),
#                       names_to = "age_cpg",
#                       values_to = "age_cpg_value")
# pheno_age_tmp |> 
#   ggplot(aes(x = age,
#              y = age_cpg_value,
#              colour = Tissue.type)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   facet_wrap(~age_cpg) +
#   ggpubr::stat_cor()
# 
# # Append some age sites
# pheno_nonage <- cbind(pheno_age, t(beta_age[nonage[1:5,]$cg,])) |> 
#   tidyr::pivot_longer(any_of(nonage[1:5,]$cg),
#                       names_to = "age_cpg",
#                       values_to = "age_cpg_value")
# pheno_nonage |> 
#   ggplot(aes(x = age,
#              y = age_cpg_value,
#              colour = Tissue.type)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   facet_wrap(~age_cpg) +
#   ggpubr::stat_cor()
# 
# 
# # Append some age sites
# pheno_pcgt <- cbind(pheno_age, t(beta_age[mPCGTgen[1:5,]$cg,])) |> 
#   tidyr::pivot_longer(any_of(mPCGTgen[1:5,]$cg),
#                       names_to = "age_cpg",
#                       values_to = "age_cpg_value")
# 
# pheno_pcgt |> 
#   ggplot(aes(x = age,
#              y = age_cpg_value,
#              colour = Tissue.type)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   facet_wrap(~age_cpg) +
#   ggpubr::stat_cor()
# 
# pheno_pcgtnonage <- cbind(pheno_age, t(beta_age[PCGTnonage[1:5,]$cg,])) |> 
#   tidyr::pivot_longer(any_of(PCGTnonage[1:5,]$cg),
#                       names_to = "age_cpg",
#                       values_to = "age_cpg_value")
# 
# pheno_pcgtnonage |> 
#   ggplot(aes(x = age,
#              y = age_cpg_value,
#              colour = Tissue.type)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   facet_wrap(~age_cpg) +
#   ggpubr::stat_cor()


# Save coefficients for index
mPCGTgen_w <- rep(1, nrow(mPCGTgen))
names(mPCGTgen_w) <- mPCGTgen$cg
mPCGTnonage_w <- rep(1, nrow(mPCGTnonage))
names(mPCGTnonage_w) <- mPCGTnonage$cg
saveRDS(mPCGTgen_w, file = "3-mouse/4-output/mPCGTgen.Rds")
saveRDS(mPCGTnonage_w, file = "3-mouse/4-output/mPCGTnonage.Rds")
