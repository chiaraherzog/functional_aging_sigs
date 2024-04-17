# Compute AUC for all groups
load("3-evaluation/3-output/data.Rdata")
source("0-data/src/computeAUC.R")


library(dplyr)

x <- data |> 
  filter(group == "current cancer")

data <- data |> 
  dplyr::mutate(celltype = ifelse(group == "current cancer" & celltype == "breast tissue" & type == "Cancer", "breast", celltype))


# Aging --------
age <- data |> dplyr::filter(celltype == "cervical sample" & type == "Control" & group == "current cancer;surrogate tissue") |> 
  dplyr::mutate(type = ifelse(menopause == "Pre", "Control", "normal aging (cervical sample)"),
                group = "normal aging")
dat <- computeAUC(age, adj.ic = T, adj.age = F)
dat_age_cerv <- dat

age <- data |> dplyr::filter(celltype == "buccal sample" & type == "Control" & group == "current cancer;surrogate tissue") |>
  dplyr::mutate(type = ifelse(menopause == "Pre", "Control", "normal aging (buccal sample)"),
                group = "normal aging")
dat <- computeAUC(age, adj.ic = T, adj.age = F)
dat_age_buccal <- dat

age <- data |> dplyr::filter(celltype == "blood sample" & type == "Control" & group == "current cancer;surrogate tissue") |>
  dplyr::mutate(type = ifelse(menopause == "Pre", "Control", "normal aging (blood sample)"),
                group = "normal aging")
dat <- computeAUC(age, adj.ic = T, adj.age = F)
dat_age_blood <- dat
 
# Current cancer --------- 
cur <- data |>
  dplyr::filter(group=="current cancer" & type != "Normal-adjacent") |> 
  dplyr::mutate(celltype = ifelse(group == "current cancer" & celltype == "breast tissue" & type != "Normal-adjacent", "breast", celltype))

for (i in unique(cur$celltype)){
  tmp <- computeAUC(cur[cur$celltype==i,], adj.ic = T, adj.age = T)
  
  if(i == unique(cur$celltype)[1]){
    dat <- tmp
  } else {
    dat <- rbind(dat, tmp)
  }
}
dat_cur <- dat

# Adjacent
cur <- data |>
  dplyr::filter(group=="current cancer" & type %in% c("Control", "Normal-adjacent") & celltype == "breast tissue")
tmp <- computeAUC(cur, adj.ic = T, adj.age = T)
dat_cur_adj <-tmp


# Current cancer, surrogate ---------
cur <- data[data$group=="current cancer;surrogate tissue",]
for (i in unique(cur$celltype)){
  tmp <- computeAUC(cur[cur$celltype==i,], adj.ic = T, adj.age = T)
  
  if(i == unique(cur$celltype)[1]){
    dat <- tmp
  } else {
    dat <- rbind(dat, tmp)
  }
}
dat_cur_sur <- dat

# Future cancer ---------
future <- data[data$group=="future cancer;tissue at risk",]
dat <- computeAUC(future, adj.ic = T, adj.age = T)
dat_future <- dat

# Future cancer; surrogate ---------
future <- data[data$group=="future cancer;surrogate tissue",]
dat <- computeAUC(future, adj.ic = T, adj.age = T)
dat_future_sur <- dat

#  Mortality buccal: all ---------
# Mortality coding
data <- data |> 
  dplyr::mutate(type2 = case_when(group == "mortality" & surv_time %in% c("1-3", "4-6", "7-9") ~ "Dead",
                                  group == "mortality" & !surv_time %in% c("1-3", "4-6", "7-9") ~ "Control",
                                  TRUE ~ type2))

mort <- data |> 
  dplyr::filter(group == "mortality" & celltype=="buccal sample") |> 
  dplyr::filter(type2 == "Control" | (type2 != "Control"))
mort$type <- mort$type2
dat <- computeAUC(mort, refgroup = "type2", adj.ic = T, adj.age = F)
dat_mort <- dat
dat_mort$group <- paste0(dat_mort$group, ";all-cause mortality")

#  Mortality: cancer ---------
mort <- data |> 
  dplyr::filter(group == "mortality" & celltype=="buccal sample") |> 
  dplyr::filter(type2 == "Control" | (type2 != "Control" & !is.na(cause_cancer)))
mort$type <- mort$type2
dat <- computeAUC(mort, refgroup = "type2", adj.ic = T, adj.age = F)
dat_mort_c <- dat
dat_mort_c$group <- paste0(dat_mort_c$group, ";cancer mortality")

#  Mortality: non-cancer ---------
mort <- data |> 
  dplyr::filter(group == "mortality" & celltype=="buccal sample") |> 
  dplyr::filter(type2 == "Control" | (type2 != "Control" & !is.na(cause_other)))
mort$type <- mort$type2
dat <- computeAUC(mort, refgroup = "type2", adj.ic = T, adj.age = F)
dat_mort_o <- dat
dat_mort_o$group <- paste0(dat_mort_o$group, ";other mortality")

#  Mortality blood: all ---------
mort <- data |> 
  dplyr::filter(group == "mortality" & celltype=="blood sample") |> 
  dplyr::filter(type2 == "Control" | (type2 != "Control"))
mort <- data[data$group == "mortality" & data$celltype=="blood sample",]
mort$type <- mort$type2
dat <- computeAUC(mort, refgroup = "type2", adj.ic = F, adj.age = F)
dat_mort_blood <- dat
dat_mort_blood$group <- paste0(dat_mort_blood$group, ";all-cause mortality")

#  Mortality: cancer ---------
mort <- data |> 
  dplyr::filter(group == "mortality" & celltype=="blood sample") |> 
  dplyr::filter(type2 == "Control" | (type2 != "Control" & !is.na(cause_cancer)))
mort$type <- mort$type2
dat <- computeAUC(mort, refgroup = "type2", adj.ic = F, adj.age = F)
dat_mort_c_blood <- dat
dat_mort_c_blood$group <- paste0(dat_mort_c_blood$group, ";cancer mortality")

# risk-increasing exposure  ---------
riskhpv <- data[data$group=="risk-increasing exposure" & data$celltype=="cervical sample",]
table(riskhpv$type, riskhpv$celltype)
dat_hpv <- computeAUC(riskhpv)

risksmk <- data[data$group=="risk-increasing exposure" & data$celltype!="cervical sample",]
table(risksmk$type, risksmk$celltype)
dat_smk <- computeAUC(risksmk, compgroups = c("Former Smoker", "Smoker"))

# risk decreasing exposure ---------
riskred <- data[data$group == "risk-reducing exposure" & grepl("WT", data$type),]
dat <- computeAUC(riskred, paired = T,adj.age = F,
                  refgroup = "Control (WT)")
dat_riskred_wt <- dat

riskred <- data[data$group == "risk-reducing exposure" & grepl("BRCA", data$type),]
dat <- computeAUC(riskred, paired = T,adj.age = F,
                  refgroup = "Control (BRCA)")
dat_riskred_brca <- dat

# OSK reprogramming ---------
osk <- data[data$group=="OSK reprogramming" & data$celltype=="fibroblast" & data$type != "intermediate reprogrammed cells" & data$accession != "NA",]
osk$type <- ifelse(osk$type=="Control", "Control", "iPSC/ESC")
dat_osk1 <- computeAUC(osk, adj.ic = F, adj.age = F)

osk_neuron <- data[data$group=="OSK reprogramming" & data$celltype=="SH-SY5Y neurons",]
dat_osk2 <- computeAUC(osk_neuron, adj.ic = F, adj.age = F)

# Proliferation ---------
prol <- data[data$group=="proliferation",]
for (i in unique(prol$celltype)){
  tmp <- computeAUC(prol[prol$celltype==i,], adj.ic = F, adj.age = F)
  
  if(i == unique(prol$celltype)[1]){
    dat <- tmp
  } else {
    dat <- rbind(dat, tmp)
  }
}
dat_prol <- dat

# Senescence---------
sen <- data[data$group=="senescence",]

for (i in unique(sen$celltype)){
  tmp <- computeAUC(sen[sen$celltype==i,], adj.ic = F, adj.age = F)
  
  if(i == unique(sen$celltype)[1]){
    dat <- tmp
  } else {
    dat <- rbind(dat, tmp)
  }
}
dat_sen <- dat


# Prognosis ---------
prog <- data[data$group=="response",]
dat_prog <- computeAUC(prog,refgroup = "Control", adj.ic = T, adj.age = T, compgroups = "Non-responder")

x <- ls()
x <- x[grepl("^dat_", x)]

save(list = x,
     file = "3-evaluation/4-output/df.Rdata")

table(sen$type, sen$celltype)
osk$type

# Proliferation
table(data$group)


# 
library(ggplot2)
dat |>
  ggplot(aes(x = d,
             y = index,
             colour = interaction(type, celltype))) +
  geom_point()

dat |>
  ggplot(aes(x = auc,
             y = index,
             colour = interaction(type, celltype))) +
  geom_point() +
  geom_errorbar(aes(xmin = cilo,
                    xmax = cihi, width = 0.2))
