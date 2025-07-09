##### Libraries  ---------------------------------
rm(list=ls())

set.seed(2712)


library("caper")
library("patchwork")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library("lattice")
#library(reshape2)
library(adephylo)
library(phylobase)
library(data.table)
library(phytools)
library("corrplot")
library("geiger")
library("ggpubr")
library(janitor)
library(RColorBrewer)
library(ggtreeExtra)
#library(rcartocolor)
library(ggtree)
#library(ggcorrplot)
library(DescTools)

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}




#### My functions   ---------------------------------



PGLS_pvalue <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  pvalue =  formatC(sum_cor$coefficients[8], digits = 3)
  if (pvalue == "   0"){ pvalue = 2.2e-16}
  return(pvalue)
}

PGLS_R2 <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  r2 = formatC(sum_cor$r.squared, digits = 2)
  return(r2)
}

PGLS_lambda <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  PGLS_lambda = sum_cor$param[2]
  return(PGLS_lambda)
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}






##### Load species tree ----------------------

vertebrate_tree <- 
  read.tree("Vertebrate_consensus_tree_nodelabel.nwk")

species_list <- vertebrate_tree$tip.label

##### Load WNT tables  ----------------------

wnt_t1 <- read.table("wnt_copies_perSpeciesclass_perWNTclade_CoelacanthandDipnoi.csv",
                     header=FALSE,
                     sep=",")
wnt_t2 <- read.table("wnt_copies_perSpeciesclass_perWNTclade_bigclasses.csv",
                     header=FALSE,
                     sep=",")

wnt_nb_df <- rbind(wnt_t1, wnt_t2)

colnames(wnt_nb_df) <- c("class", "subfamily", "copy_number", "species")

species_class <- wnt_nb_df %>% dplyr::select(species,class) %>% distinct()

#Complete the dataframe for missing wnt genes

all_combinations <- expand.grid(
  subfamily = unique(wnt_nb_df$subfamily),
  species = unique(wnt_nb_df$species)
)

wnt_nb_df <- wnt_nb_df %>% dplyr::select(-class)

wnt_nb_df_complete <- 
  all_combinations %>%
  left_join(wnt_nb_df, by = c("subfamily", "species")) %>%
  mutate(copy_number = ifelse(is.na(copy_number), 0, copy_number))

wnt_nb_df_complete <- left_join(species_class, wnt_nb_df_complete, by="species")

wnt_nb_df <- wnt_nb_df_complete


length(wnt_nb_df %>% pull(species) %>% unique())
wnt_nb_df %>% pull(copy_number) %>% sum()


wnt_nb_df %>%
  filter(class == "Agnatha") %>%
  pull(species) %>%
  unique()

wnt_nb_df %>%
  filter(class == "Dipnoi") %>%
  pull(species) %>%
  unique()


wnt_nb_df
##### Load BUSCO results ----------------------

BUSCO_df <- 
  read.table("BUSCO_summary.csv", sep=",", header=FALSE)

colnames(BUSCO_df) <- c("species", "Complete", "Complete_single",
                        "Complete_duplicated", "Fragmented", "Missing")

BUSCO_df <- BUSCO_df %>% rowwise() %>% mutate(duplicated_perc = Complete_duplicated/Complete)

BUSCO_df <- left_join(species_class, BUSCO_df, by="species") 
BUSCO_df <- BUSCO_df %>% filter(!is.na(Complete))


BUSCO_df %>%
  ggplot(., aes(x=class, y=duplicated_perc)) +
  geom_boxplot()

BUSCO_df %>% filter(class == "Actinopterygii") %>%
  filter(duplicated_perc > 0.2)

BUSCO_df %>% filter(class == "Amphibia") %>%
  filter(duplicated_perc > 0.2)

non_diploid_sp <- 
  BUSCO_df %>%
  filter(duplicated_perc > 0.2) %>%
  pull(species)

species_class_ploidy <- 
  species_class %>% mutate(ploidy = if_else(
    species %in% non_diploid_sp,
    "non-diploid",
    "diploid"
  ))



ggtree(vertebrate_tree, layout="circular") %<+% species_class_ploidy  +
  geom_tiplab(size=0.3, offset=0.07, fontface=3) +
  aes(color=as.character(ploidy)) 

##### Circle plot subfam x class  ----------------------

#Make a class tree
class_list <- wnt_nb_df %>% pull(class) %>% unique()
sp_sample <- wnt_nb_df %>% group_by(class) %>% slice_head(n=1) %>% dplyr::select(class, species)
class_tree <- keep.tip(vertebrate_tree, sp_sample$species)
class_tree$tip.label <- sp_sample$class[match(class_tree$tip.label, sp_sample$species)]

all_species_list <- wnt_nb_df$species %>% unique()


#Prepare tables of presence absence / min / max 


WNT1_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt1") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))


WNT2_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt2") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))

WNT3_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt3") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))

WNT4_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt4") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))


WNT5_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt5") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))


WNT6_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt6") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))


WNT7_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt7") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))

WNT8_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt8") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))


WNT9_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt9") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))

WNT10_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt10") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))


WNT11_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt11") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))

WNT16_summary <- 
  as.data.frame(
    wnt_nb_df %>%
      filter(subfamily == "Wnt16") %>%
      group_by(class) %>%
      summarise(min_copy_number = min(copy_number),
                max_copy_number = max(copy_number),
                mean_copy_number = mean(copy_number)) %>%
      mutate(prez_abs = if_else(
        mean_copy_number > 0,
        "Presence", 
        "Absence"
      )))




ggtree(class_tree, size=2,
       branch.length="none") +
  #geom_tiplab() + 
  geom_facet(
    panel="Wnt1",
    data=WNT1_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt1",
    data=WNT1_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt1",
    data=WNT1_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt2",
    data=WNT2_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt2",
    data=WNT2_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt2",
    data=WNT2_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt3",
    data=WNT3_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt3",
    data=WNT3_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt3",
    data=WNT3_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt4",
    data=WNT4_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt4",
    data=WNT4_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt4",
    data=WNT4_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt5",
    data=WNT5_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt5",
    data=WNT5_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt5",
    data=WNT5_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt6",
    data=WNT6_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt6",
    data=WNT6_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt6",
    data=WNT6_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt7",
    data=WNT7_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt7",
    data=WNT7_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt7",
    data=WNT7_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt8",
    data=WNT8_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt8",
    data=WNT8_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt8",
    data=WNT8_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt9",
    data=WNT9_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt9",
    data=WNT9_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt9",
    data=WNT9_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt10",
    data=WNT10_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt10",
    data=WNT10_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt10",
    data=WNT10_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt11",
    data=WNT11_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt11",
    data=WNT11_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt11",
    data=WNT11_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt16",
    data=WNT16_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=mean_copy_number),
    color="#abd9e9",
    width = 0.5
  ) +
  geom_facet(
    panel="Wnt16",
    data=WNT16_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=min_copy_number),
    color="#4575b4",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  geom_facet(
    panel="Wnt16",
    data=WNT16_summary,
    geom = geom_point,
    mapping = aes(x=1, alpha=prez_abs, size=max_copy_number),
    color="#d73027",
    shape=1,   # Empty circle
    width = 0.5,
    stroke = 1
  ) +
  scale_alpha_manual(values = c(
    Presence=1,
    Absence =0
  )) +
  scale_size(range = c(1, 12)) +
  theme(legend.position = "none",
        panel.spacing = unit(0.01, "lines")) 


##### Boxplots total WNT number  ----------------------

wnt_nb_df_total <-
  wnt_nb_df %>%
  group_by(species, class) %>%
  summarise(total_wnt = sum(copy_number))

wnt_nb_df_total %>% filter(class == "Dipnoi")

wnt_nb_df_total %>%
  group_by(class) %>%
  summarise(mean_wnt_nb = mean(total_wnt))

#reorder classes to be in the same order as in circle plot

wnt_nb_df_total$class <-
  factor(wnt_nb_df_total$class ,
         levels=c("Agnatha", "Chondrichthyes", "Actinopterygii", "Coelacanth", 
                  "Dipnoi", "Amphibia", "Mammalia", "Lepidosauria", "Aves-Crocodylia",
                  "Testudines"))


#plot

wnt_nb_df_total %>%
  ggplot(., aes(x=class, y=total_wnt)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


##### Compute phylogenetic signal WNT number  ----------------------

tip_order <- vertebrate_tree$tip.label

wnt_nb_df_total_ordered <- 
  enframe(tip_order, name = NULL, value = "species") %>% 
  left_join(wnt_nb_df_total) %>% 
  column_to_rownames("species")

phylosig_WNT <- phylosig(tree = vertebrate_tree, x = wnt_nb_df_total_ordered$total_wnt, method = "lambda", test = T)



##### Import dN/dS data  ----------------------

dNdS_df <- 
  read.table("WNTs.omega.csv", sep=",")
colnames(dNdS_df) <- c("LB", "MLE", "UB", "tree_name", "label", "dN", "dS")
dNdS_df$subfamily <- gsub(".*_", "", dNdS_df$tree_name)
dNdS_df$class <- gsub("_.*", "",(gsub("FINAL_", "", dNdS_df$tree_name)))


#Remove saturation or lack of substitutions

dNdS_df_filt <- 
  dNdS_df %>% 
  filter(dS > 0.01) %>%
  filter(dS < 1) %>%
  filter(dN < 1) %>%
  filter(MLE <= 5)

#Boxplot



dNdS_df_filt %>%
  ggplot(., aes(x=class, y=MLE, fill=subfamily)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 2)


dNdS_df_filt %>%
  ggplot(., aes(x=class, y=MLE, fill=subfamily)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 2)


dNdS_df_filt %>%
  ggplot(., aes(x=subfamily, y=MLE, fill=class)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0, 2)

dNdS_df_filt %>%
  ggplot(., aes(x=subfamily, y=MLE, fill=class)) +
  geom_boxplot() +
  ylim(0, 2)

##### Compare dN/dS of WNT genes between species classes - 1 ----------------------

dNdS_df_filt_mean <- 
  dNdS_df_filt %>%
  group_by(class, subfamily) %>%
  summarise(mean_dNdS = mean(MLE))


#Plot

dNdS_df_filt_mean %>%
  ggplot(., aes(x=class, y=mean_dNdS)) +
  geom_boxplot()

#Are there differences in dN/dS between vertebrate classes ?


all_wilcox_test_class <- 
  pairwise.wilcox.test(
    x=dNdS_df_filt_mean %>% pull(mean_dNdS),
    g=dNdS_df_filt_mean %>% pull(class),
    p.adjust.method="bonferroni"
  )


#No significant differences in the sequence evolutionary rate of WNT genes
# between vertebrate classes

dNdS_df_filt_mean %>%
  ggplot(., aes(x=class, y=mean_dNdS, color=subfamily)) +
  geom_violin() +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2))


dNdS_df_filt_mean %>%
  filter(class == "Agnatha")

##### Compare dN/dS of WNT genes between species classes ----------------------

#Test if there are differences in evolutionary rates between classes

WNT_subfamilies <- dNdS_df_filt %>% pull(subfamily) %>% unique()
WNT_dNdS_vs_class_kw_df <- as.data.frame(NULL)
for(curr_subfam in WNT_subfamilies){
  
  curr_kw_MLE <- 
    kruskal.test(MLE ~ class, 
                 data = dNdS_df_filt %>% filter(subfamily == curr_subfam))
  
  
  curr_df <- as.data.frame(cbind(curr_subfam, curr_kw_MLE$p.value))
  
  colnames(curr_df) <- c("subfamily", "pvalue")
  
  WNT_dNdS_vs_class_kw_df <- rbind(WNT_dNdS_vs_class_kw_df, curr_df)
  
}

WNT_dNdS_vs_class_kw_df$adj.pvalue <- 
  p.adjust(WNT_dNdS_vs_class_kw_df$pvalue, method="bonferroni")

WNT_dNdS_vs_class_kw_df_signif <- 
  WNT_dNdS_vs_class_kw_df %>%
  filter(adj.pvalue < 0.05)


WNT_dNdS_vs_class_kw_df_signif %>%
  arrange(adj.pvalue)

WNT_dNdS_vs_class_dunn_df <- as.data.frame(NULL)
for(curr_line_nb in 1:nrow(WNT_dNdS_vs_class_kw_df_signif)){
  curr_line <- WNT_dNdS_vs_class_kw_df_signif[curr_line_nb,]
  curr_subfam <- curr_line$subfamily
  
  curr_dunn_df <- 
    dunn_test(MLE ~ class, 
              data = dNdS_df_filt %>% filter(subfamily == curr_subfam))
  curr_dunn_df <- 
    as.data.frame(curr_dunn_df) %>%
    mutate(subfamily = curr_subfam) 
  
  WNT_dNdS_vs_class_dunn_df <- 
    rbind(WNT_dNdS_vs_class_dunn_df, curr_dunn_df)
  
}

WNT_dNdS_vs_class_dunn_df_signif <-
  WNT_dNdS_vs_class_dunn_df %>%
  filter(p.adj < 0.05)


WNT_dNdS_vs_class_dunn_df_signif %>%
  filter(subfamily == "Wnt7")
means_to_extract <- 
  WNT_dNdS_vs_class_dunn_df_signif %>%
  dplyr::select(subfamily) %>%
  distinct()

WNT_dNdS_vs_class_dunn_df_signif_mean <- as.data.frame(NULL)
for(row_nb in 1:nrow(means_to_extract)){
  curr_subfam <- means_to_extract[row_nb,]
  
  curr_dNdS_mean <- 
    as.data.frame(dNdS_df_filt %>% 
                    filter(subfamily == curr_subfam) %>%
                    group_by(class) %>%
                    summarise(mean_dNdS = mean(MLE)))
  colnames(curr_dNdS_mean) <- c("class", "mean_dNdS")
  
  curr_df <-
    WNT_dNdS_vs_class_dunn_df_signif %>%
    filter(subfamily == curr_subfam) 
  
  curr_df <- 
    left_join(curr_df, curr_dNdS_mean, 
              by = c("group1" = "class"))
  curr_df <- 
    left_join(curr_df, curr_dNdS_mean, 
              by = c("group2" = "class"))
  colnames(curr_df) <- 
    c("response","group1", "group2", "n1", "n2", "statistic", "p","p.adj","p.adj.signif",
      "subfamily", "mean_dNdS_group1", "mean_dNdS_group2")
  
  
  WNT_dNdS_vs_class_dunn_df_signif_mean <- 
    rbind(WNT_dNdS_vs_class_dunn_df_signif_mean, curr_df)
}


##### Compare dN/dS of WNT genes between subfamilies ----------------------

dNdS_df_filt %>%
  group_by(subfamily) %>%
  summarise(mean_dNdS = mean(MLE)) %>%
  arrange(mean_dNdS)


dNdS_df_filt_mean_per_class <-
  dNdS_df_filt %>%
  group_by(class, subfamily) %>%
  summarise(mean_dNdS = mean(MLE))

dNdS_df_filt_mean_per_class %>%
  ggplot(., aes(x=subfamily, y=mean_dNdS)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(color=class))



all_wilcox_test_subfams <- 
  pairwise.wilcox.test(
    x=dNdS_df_filt_mean_per_class %>% pull(mean_dNdS),
    g=dNdS_df_filt_mean_per_class %>% pull(subfamily),
    p.adjust.method="bonferroni"
  )

my_pvals <- all_wilcox_test_subfams$p.value

my_pvals[my_pvals > 0.05 & !is.na(my_pvals)] <- NA


my_transformed_pvals=-log10(my_pvals)

my_transformed_pvals[is.na(my_transformed_pvals)] <- 0



##### Fate of WNT genes after whole genome duplications in teleostei  ----------------------

non_teleost_sp <- 
  c("Erpetoichthys_calabaricus",
    "Polypterus_senegalus",
    "Polypterus_bichir_lapradei", 
    "Polyodon_spathula",
    "Acipenser_oxyrinchus_oxyrinchus",
    "Acipenser_ruthenus",
    "Amia_calva",
    "Atractosteus_spatula",
    "Lepisosteus_oculatus")

teleost_class_df <- 
  species_class_ploidy %>%
  filter(class == "Actinopterygii") %>%
  filter(! species %in% non_teleost_sp)

teleost_wnt_nb_df <- 
  left_join(teleost_class_df, wnt_nb_df, by=c("class", "species"))


teleost_WGD_sp <- teleost_class_df %>% filter(ploidy != "diploid") %>% pull(species)
teleost_noWGD_sp <- teleost_class_df %>% filter(ploidy == "diploid") %>% pull(species)


##Salmoniformes

WGD_salmo_sp <- c("Oncorhynchus_gorbuscha", "Oncorhynchus_keta", "Oncorhynchus_kisutch", "Oncorhynchus_mykiss",
                  "Oncorhynchus_tshawytscha", "Salvelinus_fontinalis", "Salvelinus_namaycush", "Salvelinus_sp_IW2-2015",
                  "Salmo_salar", "Salmo_trutta", "Hucho_hucho","Brachymystax_lenok_tsinlingensis", "Coregonus_clupeaformis",
                  "Coregonus_sp_balchen", "Coregonus_ussuriensis","Thymallus_thymallus")

outgroup_WGD_salmo_sp <- c("Esox_lucius", "Dallia_pectoralis", "Argentina_silus")


WGD_salmo_sp_wnt <- 
  teleost_wnt_nb_df %>%
  filter(species %in% WGD_salmo_sp) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_WGD = mean(copy_number)) 

outgroup_WGD_salmo_sp_wnt <- 
  teleost_wnt_nb_df %>%
  filter(species %in% outgroup_WGD_salmo_sp) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_noWGD = mean(copy_number)) 


Salmoniformes_WNT <- left_join(WGD_salmo_sp_wnt, outgroup_WGD_salmo_sp_wnt, by="subfamily")
Salmoniformes_WNT <- 
  Salmoniformes_WNT %>%
  mutate(ratio = mean_nb_WGD/mean_nb_noWGD) %>%
  dplyr::select(subfamily, ratio)


##Cyprinidae

WGD_cypr_sp <- c("Sinocyclocheilus_anophthalmus", "Sinocyclocheilus_anshuiensis", "Sinocyclocheilus_grahami",
                 "Sinocyclocheilus_maitianheensis", "Sinocyclocheilus_rhinocerous", "Carassius_auratus",
                 "Carassius_carassius", "Carassius_gibelio", "Cyprinus_carpio", "Cyprinus_carpio_carpio")

WGD_cypr_sp2 <- "Barbus_barbus"
WGD_cypr_sp3 <- "Aspiorhynchus_laticeps"
WGD_cypr_sp4 <- "Oxygymnocypris_stewartii"


outgroup_WGD_cypr_sp <- c("Schizothorax_lantsangensis","Onychostoma_macrolepis", "Acrossocheilus_wenchowensis", 
                          "Enteromius_argenteus", "Enteromius_mattozi", "Enteromius_radiatus",
                          "Enteromius_thamalakanensis", "Enteromius_treurensis", "Enteromius_trimaculatus",
                          "Enteromius_unitaeniatus", "Puntigrus_tetrazona",
                          "Ptychobarbus_kaznakovi", "Schizopygopsis_malacanthus", 
                          "Schizopygopsis_pylzovi", "Gymnocypris_eckloni_scoliostomus")


WGD_cypr_sp_wnt <- 
  teleost_wnt_nb_df %>%
  filter(species %in% WGD_cypr_sp) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_WGD = mean(copy_number)) 
WGD_cypr_sp_wnt2 <- 
  teleost_wnt_nb_df %>%
  filter(species %in% WGD_cypr_sp2) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_WGD = mean(copy_number)) 
WGD_cypr_sp_wnt3 <- 
  teleost_wnt_nb_df %>%
  filter(species %in% WGD_cypr_sp3) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_WGD = mean(copy_number)) 
WGD_cypr_sp_wnt4 <- 
  teleost_wnt_nb_df %>%
  filter(species %in% WGD_cypr_sp4) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_WGD = mean(copy_number)) 


outgroup_WGD_cypr_sp_wnt <- 
  teleost_wnt_nb_df %>%
  filter(species %in% outgroup_WGD_cypr_sp) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_noWGD = mean(copy_number)) 


Cyprinidae_WNT <- left_join(WGD_cypr_sp_wnt, outgroup_WGD_cypr_sp_wnt, by="subfamily")
Cyprinidae_WNT <- Cyprinidae_WNT %>% mutate(ratio = mean_nb_WGD/mean_nb_noWGD) %>% dplyr::select(subfamily, ratio)

Cyprinidae_WNT2 <- left_join(WGD_cypr_sp_wnt2, outgroup_WGD_cypr_sp_wnt, by="subfamily")
Cyprinidae_WNT2 <- Cyprinidae_WNT2 %>% mutate(ratio = mean_nb_WGD/mean_nb_noWGD) %>% dplyr::select(subfamily, ratio)

Cyprinidae_WNT3 <- left_join(WGD_cypr_sp_wnt3, outgroup_WGD_cypr_sp_wnt, by="subfamily")
Cyprinidae_WNT3 <- Cyprinidae_WNT3 %>% mutate(ratio = mean_nb_WGD/mean_nb_noWGD) %>% dplyr::select(subfamily, ratio)

Cyprinidae_WNT4 <- left_join(WGD_cypr_sp_wnt4, outgroup_WGD_cypr_sp_wnt, by="subfamily")
Cyprinidae_WNT4 <- Cyprinidae_WNT4 %>% mutate(ratio = mean_nb_WGD/mean_nb_noWGD) %>% dplyr::select(subfamily, ratio)


##Casto

WGD_casto_sp <- c("Xyrauchen_texanus", "Moxostoma_hubbsi", "Myxocyprinus_asiaticus")

outgroup_WGD_casto_sp <- c("Triplophysa_bombifrons", "Triplophysa_dalaica", "Triplophysa_rosa",
                           "Triplophysa_siluroides", "Triplophysa_tibetana", "Triplophysa_yarkandensis",
                           "Barbatula_barbatula", "Oreonectes_daqikongensis", "Beaufortia_kweichowensis",
                           "Paramisgurnus_dabryanus", "Misgurnus_anguillicaudatus", "Misgurnus_mizolepis")


WGD_casto_sp_wnt <- 
  teleost_wnt_nb_df %>%
  filter(species %in% WGD_casto_sp) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_WGD = mean(copy_number)) 

outgroup_WGD_casto_sp_wnt <- 
  teleost_wnt_nb_df %>%
  filter(species %in% outgroup_WGD_casto_sp) %>%
  group_by(subfamily) %>%
  summarise(mean_nb_noWGD = mean(copy_number)) 


Catostomidae_WNT <- left_join(WGD_casto_sp_wnt, outgroup_WGD_casto_sp_wnt, by="subfamily")
Catostomidae_WNT <- 
  Catostomidae_WNT %>%
  mutate(ratio = mean_nb_WGD/mean_nb_noWGD) %>%
  dplyr::select(subfamily, ratio)


#Combine dataframes

ALL_WNT_ratios_WGD <- 
  rbind(
    Salmoniformes_WNT %>% mutate(clade = "Salmoniformes"),
    Catostomidae_WNT %>% mutate(clade = "Catostomidae"),
    Cyprinidae_WNT %>% mutate(clade = "Cypriniformes_1"),
    Cyprinidae_WNT2 %>% mutate(clade = "Cypriniformes_2"),
    Cyprinidae_WNT3 %>% mutate(clade = "Cypriniformes_3"),
    Cyprinidae_WNT4 %>% mutate(clade = "Cypriniformes_4")
  ) 


ALL_WNT_ratios_WGD %>%
  ggplot(., aes(x=subfamily, y=ratio)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(color = clade), size=2) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  scale_color_manual(values = c(
    Salmoniformes="black",
    Catostomidae = "#A6CEE3", 
    Cypriniformes_1 = "lightpink1" , 
    Cypriniformes_2 = "#FF7F00" , 
    Cypriniformes_3 = "#1F78B4" , 
    Cypriniformes_4 = "#33A02C" ))

#10X6 in portrait


all_wilcox_test_subfams <- 
  pairwise.wilcox.test(
    x=ALL_WNT_ratios_WGD %>% pull(ratio),
    g=ALL_WNT_ratios_WGD %>% pull(subfamily),
    p.adjust.method="bonferroni"
  )

ALL_WNT_ratios_WGD %>%
  group_by(subfamily) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  arrange(mean_ratio)

#==> No differential retentions of WNT genes after WGD !


#Do WNT evolve faster in tetraploid species ?
teleost_sp <- c(teleost_WGD_sp, teleost_noWGD_sp)
dNdS_df_filt_teleost_sp <- as.data.frame(NULL)
for(curr_sp in teleost_sp){
  curr_sp_GC <- paste(curr_sp, "_GC", sep="")
  dNdS_df_filt_teleost_sp <- 
    rbind(dNdS_df_filt_teleost_sp,
          dNdS_df_filt %>% filter(grepl(curr_sp_GC, label)) %>% mutate(species = curr_sp)
    )
  
}


dNdS_df_filt_teleost_sp <- left_join(dNdS_df_filt_teleost_sp, teleost_class_df, by=c("species", "class"))


dNdS_df_filt_teleost_sp %>%
  ggplot(., aes(x=ploidy, y=MLE)) +
  geom_boxplot() +
  ylim(0, 1) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  ylab("dN/dS")

wilcox.test(dNdS_df_filt_teleost_sp %>% filter(ploidy == "diploid") %>% pull(MLE),
            dNdS_df_filt_teleost_sp %>% filter(ploidy == "non-diploid")%>% pull(MLE))


#Are genes in tetraploid under relaxed/intensified selection ? 

RELAX_wgd_df <- 
  read.table("Relax_results.csv",
             header=FALSE,
             sep=",")
colnames(RELAX_wgd_df) <- c("subfamily", "gene", "LRT", "pvalue", "K", "species", "ploidy")

#K < 1 => relaxed ; K > 1 => intensified

RELAX_wgd_df <- 
  RELAX_wgd_df %>%
  mutate(Selection = case_when(
    pvalue < 0.05 & K <= 1 ~ "Relaxed",
    pvalue < 0.05 & K > 1 ~ "Intensified",
    pvalue >= 0.05 ~ "No_selection"
  ))


RELAX_wgd_df_summary <- 
  RELAX_wgd_df %>% 
  group_by(species, ploidy) %>%
  summarise(
    total_genes = n(),
    relaxed_genes = sum(Selection == "Relaxed"),
    intensified_genes = sum(Selection == "Intensified"),
    proportion_relaxed = relaxed_genes / total_genes,
    proportion_intensified = intensified_genes / total_genes
  ) 

RELAX_wgd_df_summary_long <- 
  RELAX_wgd_df_summary %>%
  dplyr::select(species, ploidy, proportion_relaxed, proportion_intensified) %>%
  pivot_longer(cols=c("proportion_relaxed", "proportion_intensified"), names_to = "Selection_type", values_to = "proportion")

RELAX_wgd_df_summary_long %>%
  ggplot(., aes(x=ploidy, y=proportion, fill=Selection_type)) +
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  ylab("Proportion of WNT genes") +
  scale_fill_manual(values = c("proportion_intensified" = "#E69F00", "proportion_relaxed" ="#56B4E9"))



wilcox.test(RELAX_wgd_df_summary %>% filter(ploidy == "diploid") %>% pull(proportion_relaxed),
            RELAX_wgd_df_summary %>% filter(ploidy == "tetraploid") %>% pull(proportion_relaxed))

wilcox.test(RELAX_wgd_df_summary %>% filter(ploidy == "diploid") %>% pull(proportion_intensified),
            RELAX_wgd_df_summary %>% filter(ploidy == "tetraploid") %>% pull(proportion_intensified))

#No difference between the proportion of relaxed genes
#Slight difference in the proportion of intensified genes (maybe diploid hiher)

##=> While there is indeed a higher evolutionary rate of genes in tetraploid species (see dN/dS results in the previous
##=> section), there is probably not enough sufficient time for RELAX to detect significant realaxation or intensification
##=> on these genes. It is thus hard to decipher the contribution of both processes.


##### Compare birth/death rates across classes  ----------------------

wnt_nb_df %>% filter(species == "Danio_rerio")
wnt_nb_df$class %>% unique()

#Import Rates tables
Agnatha_rates <- read.table("Agnatha_TreeRecs_Summary/Rates_Agnatha.csv", sep=",", header=TRUE)
Actinopterygii_rates <- read.table("Rates_Actinopterygii.csv", sep=",", header=TRUE)
Amphibia_rates <- read.table("Rates_Amphibia.csv", sep=",", header=TRUE)
Aves_rates <- read.table("Rates_Aves_Crocodylia.csv", sep=",", header=TRUE)
Chondrichthyes_rates <- read.table("Rates_Chondrichthyes.csv", sep=",", header=TRUE)
Lepidosauria_rates <- read.table("Rates_Lepidosauria.csv", sep=",", header=TRUE)
Mammalia_rates <- read.table("Rates_Mammalia.csv", sep=",", header=TRUE)
Testudines_rates <- read.table("Rates_Testudines.csv", sep=",", header=TRUE)

#Combine rates tables

All_rates <- 
  rbind(
    Agnatha_rates %>% mutate(class = "Agnatha"),
    Actinopterygii_rates %>% mutate(class = "Actinopterygii"),
    Amphibia_rates %>% mutate(class = "Amphibia"),
    Aves_rates %>% mutate(class = "Aves"),
    Chondrichthyes_rates %>% mutate(class = "Chondrichthyes"),
    Lepidosauria_rates %>% mutate(class = "Lepidosauria"),
    Mammalia_rates %>% mutate(class = "Mammalia"),
    Testudines_rates %>% mutate(class = "Testudines")
  )

#write.table(All_rates, "/scicore/home/salzburg/fogg0000/WNT_DataRepository/Notung_results.csv",
#            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

#Remove small branches and terminal branches which lead to overestimations of rates

All_rates_filt <- 
  All_rates %>% 
  filter(branch_length > 2) 


#Compute mean birth and death rates per classes per subfamily and keep
#only classes with at-least 50 observations .. 

All_rates_filt_mean <- 
  All_rates_filt %>%
  group_by(class, subfamily) %>%
  summarise(mean_birth_rate = mean(Birth_rate),
            mean_death_rate = mean(Death_rate),
            count = n()) %>%
  filter(count > 50) %>%
  dplyr::select(-count)



All_rates_filt_mean <-
  All_rates_filt_mean %>%
  pivot_longer(cols = c(mean_birth_rate, mean_death_rate), 
               names_to = "rate_type", 
               values_to = "value")

All_rates_filt_mean <- All_rates_filt_mean %>% filter(subfamily != "Total")

#Alos make an actinopteryigy + amphibian categories without WGD clades

actino_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.goodlength.nwk.nodelabel")
WGD_salmo_sp <- c("Oncorhynchus_gorbuscha", "Oncorhynchus_keta", "Oncorhynchus_kisutch", "Oncorhynchus_mykiss",
                  "Oncorhynchus_tshawytscha", "Salvelinus_fontinalis", "Salvelinus_namaycush", "Salvelinus_sp_IW2-2015",
                  "Salmo_salar", "Salmo_trutta", "Hucho_hucho","Brachymystax_lenok_tsinlingensis", "Coregonus_clupeaformis",
                  "Coregonus_sp_balchen", "Coregonus_ussuriensis","Thymallus_thymallus")
WGD_salmo_nodes <- keep.tip(actino_tree, WGD_salmo_sp)$node.label
WGD_cypr_sp <- c("Sinocyclocheilus_anophthalmus", "Sinocyclocheilus_anshuiensis", "Sinocyclocheilus_grahami",
                 "Sinocyclocheilus_maitianheensis", "Sinocyclocheilus_rhinocerous", "Carassius_auratus",
                 "Carassius_carassius", "Carassius_gibelio", "Cyprinus_carpio", "Cyprinus_carpio_carpio")
WGD_cypr_nodes <- keep.tip(actino_tree, WGD_cypr_sp)$node.label
WGD_cypr_sp2 <- "Barbus_barbus"
WGD_cypr_sp3 <- "Aspiorhynchus_laticeps"
WGD_cypr_sp4 <- "Oxygymnocypris_stewartii"
WGD_casto_sp <- c("Xyrauchen_texanus", "Moxostoma_hubbsi", "Myxocyprinus_asiaticus")
WGD_casto_nodes <- keep.tip(actino_tree, WGD_casto_sp)$node.label

all_wgd_branches <- c(WGD_salmo_sp, WGD_salmo_nodes, WGD_cypr_sp, WGD_cypr_nodes,
                      WGD_cypr_sp2, WGD_cypr_sp3, WGD_cypr_sp4, WGD_casto_sp,
                      WGD_casto_nodes)

actino_wo_wgd <- 
  All_rates_filt %>%
  filter(class == "Actinopterygii") %>%
  filter(! label %in% all_wgd_branches)

actino_wo_wgd$class <- gsub("Actinopterygii", "Actinopterygii_woWGD", actino_wo_wgd$class)

All_rates_filt <- rbind(All_rates_filt, actino_wo_wgd)

amphib_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.goodlength.nwk.nodelabel")
WGD_amphib <- c("Xenopus_laevis", "Xenopus_borealis")
WGD_amphb_nodes <- keep.tip(amphib_tree, WGD_amphib)$node.label
WGD_amphib <- c(WGD_amphib, WGD_amphb_nodes)

amphib_wo_wgd <- 
  All_rates_filt %>%
  filter(class == "Amphibia") %>%
  filter(! label %in% WGD_amphib)

amphib_wo_wgd$class <- gsub("Amphibia", "Amphibia_woWGD", amphib_wo_wgd$class)
All_rates_filt <- rbind(All_rates_filt, amphib_wo_wgd)

All_rates_filt_mean <- 
  All_rates_filt %>%
  group_by(class, subfamily) %>%
  summarise(mean_birth_rate = mean(Birth_rate),
            mean_death_rate = mean(Death_rate),
            count = n()) %>%
  filter(count > 50) %>%
  dplyr::select(-count)



All_rates_filt_mean <-
  All_rates_filt_mean %>%
  pivot_longer(cols = c(mean_birth_rate, mean_death_rate), 
               names_to = "rate_type", 
               values_to = "value")

All_rates_filt_mean <- All_rates_filt_mean %>% filter(subfamily != "Total")


#Plot

All_rates_filt_mean %>%
  filter(rate_type == "mean_birth_rate") %>%
  ggplot(., aes(x=class, y=value, fill=class)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  scale_fill_manual(values = c(
    Agnatha="black",
    Actinopterygii = "#A6CEE3", 
    Actinopterygii_woWGD = "#A6CEE3", 
    Amphibia = "lightpink1" , 
    Amphibia_woWGD = "lightpink1" , 
    Aves = "#FF7F00" , 
    Chondrichthyes = "#1F78B4" , 
    Coelacanth = "#B2DF8A" , 
    Dipnoi = "#33A02C" , 
    Lepidosauria = "#FFFF99" , 
    Mammalia = "#E31A1C" , 
    Testudines = "#6A3D9A")) +
  ylab("Birth rates")



All_rates_filt_mean %>%
  filter(rate_type == "mean_death_rate") %>%
  ggplot(., aes(x=class, y=value, fill=class)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  scale_fill_manual(values = c(
    Agnatha="black",
    Actinopterygii = "#A6CEE3", 
    Actinopterygii_woWGD = "#A6CEE3", 
    Amphibia = "lightpink1" , 
    Amphibia_woWGD = "lightpink1" , 
    Aves = "#FF7F00" , 
    Chondrichthyes = "#1F78B4" , 
    Coelacanth = "#B2DF8A" , 
    Dipnoi = "#33A02C" , 
    Lepidosauria = "#FFFF99" , 
    Mammalia = "#E31A1C" , 
    Testudines = "#6A3D9A")) +
  ylab("Death rates")

#10X3 PORTRAIT
#5X6 PORTRAIT
#Are there differences in birth rates across vertebrate classes ?


all_wilcox_test_class <- 
  pairwise.wilcox.test(
    x=All_rates_filt_mean %>% filter(rate_type == "mean_birth_rate") %>% pull(value),
    g=All_rates_filt_mean %>% filter(rate_type == "mean_birth_rate") %>% pull(class),
    p.adjust.method="bonferroni"
  )

#=> The WNT gene repertoire seem to evolve faster in mammals and ray-finned fishes 
# than in other vertebrates.


##=> for ray finned fishes, remove clades corresponding to whole genome duplications


##=> order by alphabetical + color according to class +
##=> 

All_rates_filt_mean %>%
  filter(class == "Actinopterygii") %>%
  filter(rate_type == "mean_birth_rate") %>%
  arrange(value)


All_rates_filt_mean %>%
  filter(class == "Mammalia") %>%
  filter(rate_type == "mean_birth_rate") %>%
  arrange(value)

##### Figures with rates  ----------------------


All_rates_filt %>%
  group_by(class, subfamily) %>%
  summarise(mean_birth_rate = mean(Birth_rate),
            mean_death_rate = mean(Death_rate)) %>%
  filter(class == "Testudines") %>%
  ggplot(., aes(x=mean_birth_rate, y=mean_death_rate, color=subfamily)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_manual(values = wnt_clades_colors) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position = "none") +
  ylab("Mean death rate") +
  xlab("Mean birth rate")



##### Link birth rates and dN/dS  ----------------------
mean_rates_subfam <- 
  All_rates_filt %>%
  group_by(subfamily) %>%
  summarise(mean_birth_rate = mean(Birth_rate))

mean_dNdS_subfam <- 
  dNdS_df_filt %>%
  group_by(subfamily) %>%
  summarise(mean_dNdS = mean(MLE))

mean_rates_dNdS_subfam <- 
  left_join(mean_rates_subfam, mean_dNdS_subfam, by="subfamily") %>%
  filter(subfamily != "Total")

mean_rates_dNdS_subfam %>%
  arrange(mean_birth_rate)

lm_rates_dNdS <- 
  lm(data = mean_rates_dNdS_subfam, 
     mean_dNdS ~ mean_birth_rate)
summary_lm_rates_dNdS <- summary(lm_rates_dNdS)
rates_dNdS_function <- GLS_function(lm_rates_dNdS)
lm_rates_dNdS_r2 = formatC(summary_lm_rates_dNdS$r.squared, digits = 2)
lm_rates_dNdS_pval = formatC(summary_lm_rates_dNdS$coefficients[8], digits = 3)

wnt_clades_colors <- 
  c("Wnt1"="#a6cee3",
    "Wnt2"="#1f78b4",
    "Wnt3"="#b2df8a",
    "Wnt4"="#33a02c",
    "Wnt5"="#fb9a99",
    "Wnt6"="#e31a1c",
    "Wnt7"="#fdbf6f",
    "Wnt8"="#ff7f00",
    "Wnt9"="#cab2d6",
    "Wnt10"="#6a3d9a",
    "Wnt11"="#ffff99",
    "Wnt16"="#b15928"
  )

mean_rates_dNdS_subfam %>%
  ggplot(., aes(x=mean_birth_rate, y=mean_dNdS, color=subfamily)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_manual(values = wnt_clades_colors) +
  stat_function(fun = rates_dNdS_function, color="black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  ylab("mean dN/dS") +
  xlab("mean birth rate")


##### Correlation dN/dS and copy number  ----------------------

all_sp <- vertebrate_tree$tip.label
dNdS_df_filt_sp <- as.data.frame(NULL)
for(curr_sp in all_sp){
  curr_sp_GC <- paste(curr_sp, "_GC", sep="")
  dNdS_df_filt_sp <- 
    rbind(dNdS_df_filt_sp,
          dNdS_df_filt %>% filter(grepl(curr_sp_GC, label)) %>% mutate(species = curr_sp)
    )
  
}

dNdS_df_mean_sub <- 
  dNdS_df_filt_sp %>%
  group_by(species, subfamily) %>%
  summarise(mean_dNdS = mean(MLE))


dNdS_copynb_df <- 
  left_join(wnt_nb_df, dNdS_df_mean_sub,
            by=c("species", "subfamily")) %>%
  filter(!is.na(mean_dNdS))

dNdS_copynb_df_wide <- 
  as.data.frame(dNdS_copynb_df %>%
                  pivot_wider(names_from = subfamily,
                              values_from = c(copy_number, mean_dNdS)
                  )) 


caper_dNdS_copynb <- 
  comparative.data(
    phy = vertebrate_tree,
    data = dNdS_copynb_df_wide,
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


#Perform pGLS of dN/dS vs copy number
Wnt1_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt1 ~ copy_number_Wnt1, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt2_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt2 ~ copy_number_Wnt2, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt3_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt3 ~ copy_number_Wnt3, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt4_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt4 ~ copy_number_Wnt4, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )


Wnt5_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt5 ~ copy_number_Wnt5, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )


Wnt6_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt6 ~ copy_number_Wnt6, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt7_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt7 ~ copy_number_Wnt7, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt8_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt8 ~ copy_number_Wnt8, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt9_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt9 ~ copy_number_Wnt9, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt10_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt10 ~ copy_number_Wnt10, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Wnt11_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt11 ~ copy_number_Wnt11, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )


Wnt16_dNdS_vs_Nb <-
  pgls(mean_dNdS_Wnt16 ~ copy_number_Wnt16, 
       data = caper_dNdS_copynb , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )



##### Correlation between copy numbers  ----------------------

wnt_nb_df_wide <-
  as.data.frame(wnt_nb_df %>%
                  pivot_wider(names_from = subfamily,
                              values_from = copy_number
                  )) 

#remove species with less than 6 wnt genes

sp_to_remove <- 
  wnt_nb_df %>%
  group_by(species) %>%
  summarise(total_wnt = sum(copy_number)) %>%
  filter(total_wnt < 6) %>% pull(species)

#remove tetraploids species 

sp_to_remove <- c(sp_to_remove, species_class_ploidy %>% filter(ploidy != "diploid") %>% pull(species))


#Winsorize to mitigate the effect of outliers on the pGLS

WNT_subfam_list <- wnt_nb_df$subfamily %>% unique()

for (col in WNT_subfam_list) {
  new_col_name <- paste0(col, "_winsorized")
  #win_col <- Winsorize(x=wnt_nb_df_wide %>% pull(!!sym(col)), probs = c(0.025, 0.975), na.rm = FALSE)
  win_col <- Winsorize(wnt_nb_df_wide %>% pull(!!sym(col)), val = quantile(wnt_nb_df_wide %>% pull(!!sym(col)),  probs = c(0.025, 0.975), na.rm = FALSE))
  wnt_nb_df_wide <- 
    wnt_nb_df_wide %>%
    mutate(!!new_col_name := win_col) 
  
  
}


caper_copynb_all <- 
  comparative.data(
    phy = vertebrate_tree,
    data = wnt_nb_df_wide %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)

caper_copynb_actino <- 
  comparative.data(
    phy = keep.tip(vertebrate_tree, wnt_nb_df_wide %>% filter(class == "Actinopterygii") %>% filter(! species %in% sp_to_remove) %>% pull(species)),
    data = wnt_nb_df_wide %>% filter(class == "Actinopterygii") %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)

caper_copynb_mammals <- 
  comparative.data(
    phy = keep.tip(vertebrate_tree, wnt_nb_df_wide %>% filter(class == "Mammalia") %>% filter(! species %in% sp_to_remove) %>% pull(species)),
    data = wnt_nb_df_wide %>% filter(class == "Mammalia") %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)

caper_copynb_birds <- 
  comparative.data(
    phy = keep.tip(vertebrate_tree, wnt_nb_df_wide %>% filter(class == "Aves-Crocodylia") %>% filter(! species %in% sp_to_remove) %>% pull(species)),
    data = wnt_nb_df_wide %>% filter(class == "Aves-Crocodylia") %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


caper_copynb_amphibia <- 
  comparative.data(
    phy = keep.tip(vertebrate_tree, wnt_nb_df_wide %>% filter(class == "Amphibia") %>% filter(! species %in% sp_to_remove) %>% pull(species)),
    data = wnt_nb_df_wide %>% filter(class == "Amphibia") %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)

caper_copynb_lepidosauria <- 
  comparative.data(
    phy = keep.tip(vertebrate_tree, wnt_nb_df_wide %>% filter(class == "Lepidosauria") %>% filter(! species %in% sp_to_remove) %>% pull(species)),
    data = wnt_nb_df_wide %>% filter(class == "Lepidosauria") %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


caper_copynb_testudines <- 
  comparative.data(
    phy = keep.tip(vertebrate_tree, wnt_nb_df_wide %>% filter(class == "Testudines") %>% filter(! species %in% sp_to_remove) %>% pull(species)),
    data = wnt_nb_df_wide %>% filter(class == "Testudines") %>% filter(! species %in% sp_to_remove),
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


#Lets run pGLS ! 



WNT_subfam_list <- wnt_nb_df$subfamily %>% unique()
WNT_subfam_list <- gsub("$", "_winsorized", WNT_subfam_list)
gene_pairs <- combn(WNT_subfam_list, 2, simplify = FALSE)
pGLS_copynumber_pairwise_df <- as.data.frame(NULL)
for(pair in gene_pairs) {
  first_wnt <- pair[1]
  second_wnt <- pair[2]
  curr_formula <- as.formula(paste(first_wnt, "~", second_wnt))
  
  pGLS_result <- pgls(curr_formula, 
                      data = caper_copynb_all, 
                      lambda = 1)
  
  sum_fit_phy <- summary(pGLS_result)
  PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
  PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
  PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
  if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
  slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
  
  print(PGLS_pvalue)
  
  curr_df <- 
    as.data.frame(
      cbind(
        first_wnt,
        second_wnt,
        PGLS_r2,
        PGLS_pvalue,
        slope))
  
  colnames(curr_df) <- c("subfamily_1", "subfamily_2", "R2", "pvalue", "slope")
  
  pGLS_copynumber_pairwise_df <- rbind(pGLS_copynumber_pairwise_df, curr_df)
  
  
}


pGLS_copynumber_pairwise_df <- 
  read.table("pGLS_copynumber_pairwise_df.csv",
             sep=",",
             header=TRUE)
pGLS_copynumber_pairwise_df$pvalue

pGLS_copynumber_pairwise_df$pvalue <- gsub("< 2.2e-16", "2.2e-16", pGLS_copynumber_pairwise_df$pvalue)
pGLS_copynumber_pairwise_df$pvalue <- as.numeric(pGLS_copynumber_pairwise_df$pvalue)


pGLS_copynumber_pairwise_df %>%
  arrange(R2)


##### Electrona antartica selection test  ----------------------

absrel_electrona_df <- 
  read.table("Electrona_selection_analysis/Electrona_aBSREL_results.csv", 
             sep=",",
             header=FALSE)

colnames(absrel_electrona_df) <- c("branch", "LRT", "adj.pval", "pval")

absrel_electrona_df %>% filter(adj.pval < 0.05)


#=> Electrona WNT8 genes seem to be, for most of them (14/16), under negative selection
#=> only 2 are under diversifying selection



##### Limbloss -  phylogeny  ----------------------


lepido_amphibia_tree <- 
  keep.tip(vertebrate_tree, 
           wnt_nb_df %>% filter(class %in% c("Lepidosauria", "Amphibia")) %>% pull(species) %>% unique())

lepido_amphibia_sp_df <- 
  wnt_nb_df %>% filter(class %in% c("Lepidosauria", "Amphibia")) %>%
  dplyr::select(species, class) %>%
  distinct()

limbloss_sp_list <- 
  scan("limbless_sp_list.txt",
       what="character")


lepido_amphibia_sp_df <- 
  lepido_amphibia_sp_df %>%
  mutate(limb = if_else(
    species %in% limbloss_sp_list,
    "no",
    "yes"
  ))

ggtree(lepido_amphibia_tree, layout="circular") %<+% lepido_amphibia_sp_df + 
  aes(color=limb) +
  geom_tiplab(size = 2, fontface=3) 


##### Limbloss - load some data    ----------------------

#remove species with less than 6 wnt genes

sp_to_remove <- 
  wnt_nb_df %>%
  group_by(species) %>%
  summarise(total_wnt = sum(copy_number)) %>%
  filter(total_wnt < 6) %>% pull(species)

lepido_amphib_sp <- 
  wnt_nb_df %>% filter(class %in% c("Lepidosauria", "Amphibia")) %>%
  dplyr::select(species, class) %>%
  filter(! species %in% sp_to_remove) %>%
  distinct() %>%
  pull(species)


#Now for copy number and dN/dS, also remove tetraploid xenopus species


wnt_nb_df_lepido_amphib <- 
  wnt_nb_df %>% filter(class %in% c("Lepidosauria", "Amphibia")) %>%
  filter(! species %in% sp_to_remove) %>%
  filter(! species %in% c("Xenopus_laevis", "Xenopus_borealis")) %>%
  distinct() 
wnt_nb_df_lepido_amphib_wide <- 
  as.data.frame(wnt_nb_df_lepido_amphib %>%
                  pivot_wider(names_from = subfamily,
                              values_from = copy_number)) 



all_sp <- vertebrate_tree$tip.label
dNdS_df_filt_sp <- as.data.frame(NULL)
for(curr_sp in all_sp){
  curr_sp_GC <- paste(curr_sp, "_GC", sep="")
  dNdS_df_filt_sp <- 
    rbind(dNdS_df_filt_sp,
          dNdS_df_filt %>% filter(grepl(curr_sp_GC, label)) %>% mutate(species = curr_sp)
    )
  
}

dNdS_df_mean_sub <- 
  dNdS_df_filt_sp %>%
  group_by(species, subfamily) %>%
  summarise(mean_dNdS = mean(MLE))

dNdS_df_filt_lepido_amphib <- 
  dNdS_df_mean_sub %>% filter(species %in% lepido_amphib_sp) %>%
  filter(! species %in% sp_to_remove) %>%
  filter(! species %in% c("Xenopus_laevis", "Xenopus_borealis")) %>%
  distinct() 



dNdS_df_filt_lepido_amphib_wide <- 
  as.data.frame(dNdS_df_filt_lepido_amphib %>%
                  pivot_wider(names_from = subfamily,
                              values_from = mean_dNdS)) 


colnames(dNdS_df_filt_lepido_amphib_wide) <- gsub("Wnt", "dNdS_Wnt", colnames(dNdS_df_filt_lepido_amphib_wide))
colnames(wnt_nb_df_lepido_amphib_wide) <- gsub("Wnt", "CN_Wnt", colnames(wnt_nb_df_lepido_amphib_wide))


dNdS_copynb_df_lepido_amphib <- 
  left_join(wnt_nb_df_lepido_amphib_wide , dNdS_df_filt_lepido_amphib_wide, by="species")


dNdS_copynb_df_lepido_amphib <- 
  dNdS_copynb_df_lepido_amphib %>%
  mutate(limb = if_else(
    species %in% limbloss_sp_list,
    "No",
    "Yes"
  ))




##### Limbloss - Subfamilies -  pGLS dN/dS  ----------------------



#Get the copy number per wnt subfamily per species

WNT_copynumber_lepido_amphib <- 
  read.table("count_wnt_species.csv",
             header=FALSE,
             sep=",")
colnames(WNT_copynumber_lepido_amphib) <- c("subfamily", "copy_number", "species")
WNT_copynumber_lepido_amphib_wide <- 
  as.data.frame(WNT_copynumber_lepido_amphib %>%
                  pivot_wider(names_from = subfamily,
                              values_from = copy_number,
                              values_fill = 0)) 



#Get the dN/dS per wnt subfamily per species

all_sp <- vertebrate_tree$tip.label
dNdS_df_filt_sp <- as.data.frame(NULL)
for(curr_sp in all_sp){
  curr_sp_GC <- paste(curr_sp, "_GC", sep="")
  dNdS_df_filt_sp <- 
    rbind(dNdS_df_filt_sp,
          dNdS_df_filt %>% filter(grepl(curr_sp_GC, label)) %>% mutate(species = curr_sp)
    )
  
}

dNdS_df_filt_sp$clade <- gsub(".*_", "", dNdS_df_filt_sp$label)
dNdS_df_filt_sp$clade <- gsub("^1$", "Wnt8a.1", dNdS_df_filt_sp$clade)
dNdS_df_filt_sp$clade <- gsub("^2$", "Wnt8a.2", dNdS_df_filt_sp$clade)



dNdS_df_filt_sp$clade_precise <- 
  gsub(".*_Wnt", "Wnt", dNdS_df_filt_sp$label)




dNdS_df_mean_clade_precise <- 
  dNdS_df_filt_sp %>%
  group_by(species, clade_precise) %>%
  summarise(mean_dNdS = mean(MLE))

dNdS_df_filt_lepido_amphib_precise <- 
  dNdS_df_mean_clade_precise %>% filter(species %in% lepido_amphib_sp) %>%
  filter(! species %in% sp_to_remove) %>%
  filter(! species %in% c("Xenopus_laevis", "Xenopus_borealis")) %>%
  distinct() 



dNdS_df_filt_lepido_amphib_precise_wide <- 
  as.data.frame(dNdS_df_filt_lepido_amphib_precise %>%
                  pivot_wider(names_from = clade_precise,
                              values_from = mean_dNdS)) 


dNdS_copynb_df_lepido_amphib <- 
  dNdS_df_filt_lepido_amphib_precise_wide


dNdS_copynb_df_lepido_amphib <- 
  dNdS_copynb_df_lepido_amphib %>%
  mutate(limb = if_else(
    species %in% limbloss_sp_list,
    "No",
    "Yes"
  ))




Wnt_all_lepido_amphib <- 
  dNdS_copynb_df_lepido_amphib %>%
  mutate(clade = case_when(
    species %in% c("Microcaecilia_unicolor", "Geotrypetes_seraphini", "Rhinatrema_bivittatum") ~ "Caecilian",
    species %in% c("Anniella_stebbinsi", "Rhineura_floridana", "Calyptommatus_sinebrachiatus", "Lerista_edwardsae") ~ "Limbless_lizard",
    limb == "No" & (!species %in% c("Microcaecilia_unicolor", "Geotrypetes_seraphini", "Rhinatrema_bivittatum","Anniella_stebbinsi", "Rhineura_floridana", "Calyptommatus_sinebrachiatus", "Lerista_edwardsae")) ~ "Snakes",
    limb == "Yes"  ~ "Other"
  ))


#Let's launch pGLS ! 


caper_dNdS_copynb <- 
  comparative.data(
    phy = vertebrate_tree,
    data = Wnt_all_lepido_amphib,
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


WNT_subfam_list <- colnames(dNdS_copynb_df_lepido_amphib %>% dplyr::select(-c(limb, species)))
pGLS_limb_df <- as.data.frame(NULL)
for(curr_wnt in WNT_subfam_list) {
  
  
  nrow_caecilian <-
    nrow(Wnt_all_lepido_amphib %>%
           filter(! is.na(.data[[curr_wnt]])) %>%
           filter(clade == "Caecilian"))
  nrow_snakes <- 
    nrow(Wnt_all_lepido_amphib %>%
           filter(! is.na(.data[[curr_wnt]])) %>%
           filter(clade == "Snakes"))
  nrow_limbless_lizard <- 
    nrow(Wnt_all_lepido_amphib %>%
           filter(! is.na(.data[[curr_wnt]])) %>%
           filter(clade == "Limbless_lizard"))  
  
  
  
  curr_formula_dNdS <- as.formula(paste(curr_wnt, " ~ ", "limb", sep=""))
  
  
  if((nrow_caecilian > 0) & (nrow_snakes > 0) & (nrow_limbless_lizard > 0)){
    
    pGLS_dNdS <- pgls(curr_formula_dNdS, data = caper_dNdS_copynb, lambda = "ML")
    
    sum_fit_phy <- summary(pGLS_dNdS)
    PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
    PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
    PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
    if (PGLS_pvalue == "   0"){ PGLS_pvalue = 2.2e-16}
    pGLS_lambda <- sum_fit_phy$param[2]
    print(PGLS_pvalue)
    curr_df_dNdS <- as.data.frame(cbind(curr_wnt, "dNdS", PGLS_r2, PGLS_pvalue, pGLS_lambda))
    colnames(curr_df_dNdS) <- c("subfamily", "response", "R2", "pvalue", "lambda")
    
    
    
    pGLS_limb_df <- rbind(pGLS_limb_df, curr_df_dNdS)
    
  } 
  
}

pGLS_limb_df$pvalue <- as.numeric(pGLS_limb_df$pvalue)
pGLS_limb_df$R2 <- as.numeric(pGLS_limb_df$R2)
pGLS_limb_df$adj.pvalue <- p.adjust(pGLS_limb_df$pvalue, method="fdr")
pGLS_limb_df %>% filter(adj.pvalue < 0.05)


Wnt_all_lepido_amphib_long <- 
  Wnt_all_lepido_amphib %>%
  dplyr::select(limb, Wnt4,Wnt10b) %>%
  pivot_longer(!limb, names_to = "subfamily", values_to = "dNdS")


Wnt_all_lepido_amphib_long %>%
  ggplot(., aes(x=subfamily, y=dNdS, fill=limb)) +
  geom_boxplot() +
  scale_fill_manual(values=c("No" = "#FE6100", "Yes" = "#648FFF")) +
  xlab("Subfamily") +
  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position = "none") 




ggtree(lepido_amphibia_tree, layout="circular", size=1) %<+% lepido_amphibia_sp_df + 
  aes(color=limb) +
  scale_color_manual(values=c("no" = "#FE6100", "yes" = "#648FFF")) +
  geom_tiplab(size = 4, fontface=3) +
  theme(legend.position = "none") +
  xlim(0, 500)


#Seems like there is an accelerated evolution of WNT10b and Wnt4 in species w/o limb ? 


lepido_amphibia_tree_filt <- 
  keep.tip(vertebrate_tree, 
           dNdS_copynb_df_lepido_amphib$species %>% unique())


limb_df <- 
  dNdS_copynb_df_lepido_amphib %>%
  dplyr::select(species, limb)
dNdS_copynb_df_lepido_amphib_g <- 
  dNdS_copynb_df_lepido_amphib %>%
  dplyr::select(-limb)

ggtree(lepido_amphibia_tree_filt, size = 1, layout="circular") %<+% dNdS_copynb_df_lepido_amphib_g + 
  aes(color=Wnt10b) +
  new_scale_fill()  +
  geom_fruit(
    data=limb_df,
    geom=geom_tile,
    mapping=aes(y=species, fill=as.factor(limb)),
    size=5,
    width=5,
    #color="white",
    pwidth=1.5,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  new_scale_fill() + 
  geom_tiplab(size = 2, fontface=3, offset=18)  +
  theme(legend.position = "none") 



ggtree(lepido_amphibia_tree_filt, size = 1, layout="circular") %<+% dNdS_copynb_df_lepido_amphib_g + 
  aes(color=Wnt4) +
  new_scale_fill()  +
  geom_fruit(
    data=limb_df,
    geom=geom_tile,
    mapping=aes(y=species, fill=as.factor(limb)),
    size=5,
    width=5,
    #color="white",
    pwidth=1.5,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  new_scale_fill() + 
  geom_tiplab(size = 2, fontface=3, offset=18)  +
  theme(legend.position = "none") 

#Relaunch pGLS without 

Wnt_all_lepido_amphib_woLerista <- 
  Wnt_all_lepido_amphib %>%
  filter(species != "Lerista_edwardsae")


caper_dNdS_copynb <- 
  comparative.data(
    phy = vertebrate_tree,
    data = Wnt_all_lepido_amphib_woLerista,
    names.col = species, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


pGLS_limb_df_woLerista <- as.data.frame(NULL)
for(curr_wnt in WNT_subfam_list) {
  
  
  nrow_caecilian <-
    nrow(Wnt_all_lepido_amphib %>%
           filter(! is.na(.data[[curr_wnt]])) %>%
           filter(clade == "Caecilian"))
  nrow_snakes <- 
    nrow(Wnt_all_lepido_amphib %>%
           filter(! is.na(.data[[curr_wnt]])) %>%
           filter(clade == "Snakes"))
  nrow_limbless_lizard <- 
    nrow(Wnt_all_lepido_amphib %>%
           filter(! is.na(.data[[curr_wnt]])) %>%
           filter(clade == "Limbless_lizard"))  
  
  
  
  curr_formula_dNdS <- as.formula(paste(curr_wnt, " ~ ", "limb", sep=""))
  
  
  if((nrow_caecilian > 0) & (nrow_snakes > 0) & (nrow_limbless_lizard > 0)){
    
    pGLS_dNdS <- pgls(curr_formula_dNdS, data = caper_dNdS_copynb, lambda = "ML")
    
    sum_fit_phy <- summary(pGLS_dNdS)
    PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
    PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
    PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
    if (PGLS_pvalue == "   0"){ PGLS_pvalue = 2.2e-16}
    pGLS_lambda <- sum_fit_phy$param[2]
    print(PGLS_pvalue)
    curr_df_dNdS <- as.data.frame(cbind(curr_wnt, "dNdS", PGLS_r2, PGLS_pvalue, pGLS_lambda))
    colnames(curr_df_dNdS) <- c("subfamily", "response", "R2", "pvalue", "lambda")
    
    
    
    pGLS_limb_df_woLerista <- rbind(pGLS_limb_df_woLerista, curr_df_dNdS)
    
  } 
  
}

pGLS_limb_df_woLerista$pvalue <- as.numeric(pGLS_limb_df_woLerista$pvalue)
pGLS_limb_df_woLerista$R2 <- as.numeric(pGLS_limb_df_woLerista$R2)
pGLS_limb_df_woLerista$adj.pvalue <- p.adjust(pGLS_limb_df_woLerista$pvalue, method="fdr")
pGLS_limb_df_woLerista %>% filter(adj.pvalue < 0.05)
