######################################################################
# Analysis of standard length and spine data
# Kieran Samuk - Apr 2016
######################################################################

######################################################################
# Libraries
######################################################################

library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("broom")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select<- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#77AB43", "#008FD5", "#BDBDBD")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
morpho_df <- read.csv("data/whtstbk_morphometrics_all.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
morpho_df <- morpho_df %>%
  mutate(id = paste0(population, individual) %>% gsub("CP[a-zA-Z]*", "CP", .))

# join meta data to morpho data

# prep meta data for join

meta_df <- meta_df %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex)

# join in cluster/geno sex and remove landmark data (explored in 01_body...)
morpho_df <- left_join(morpho_df, meta_df, by = "id") %>%
  select(-matches("^[x,y]{1}\\.")) %>%
  select(-matches("^lndmrk"))

######################################################################
# standard length
######################################################################

morpho_df <- morpho_df %>%
  select(-species, -photo.id, -individual, -notes)

morpho_df %>%
  ggplot(aes(y = std.length, x = cluster, color = cluster)) +
  geom_jitter()

morpho_long <- gather(morpho_df, key = trait, value = value, -population, -year, -sex, -id, -cluster, -geno_sex)

morpho_long %>%
  filter(!is.na(cluster)) %>%
  filter(!is.na(geno_sex)) %>%
  filter(value < 10) %>%
  filter(trait != "spine.length.3") %>%
  filter(trait != "std.length.1") %>%
  ggplot(aes(y = value, x = geno_sex, color = cluster)) +
  geom_jitter() +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               geom = "crossbar", width = 0.9, color = "black", fatten = 2) +
  facet_grid(trait~cluster, scales = "free") +
  scale_color_manual(values = whtstbk_palatte)

######################################################################
# size correction
######################################################################

morpho_corr_long <- gather(morpho_df, key = trait, value = value, -population, -year, -sex, -id, -cluster, -geno_sex, -std.length)

morpho_corr_long <- morpho_corr_long %>%
  filter(!is.na(cluster)) %>%
  filter(!is.na(geno_sex)) %>%
  filter(!is.na(std.length)) %>%
  filter(value < 10) %>%
  filter(trait != "spine.length.3") %>%
  filter(trait != "std.length.1") 

resid <- morpho_corr_long %>%
  group_by(cluster, trait) %>%
  do(augment(lm(value ~ std.length, data=.))) %>%
  ungroup %>%
  select(.resid) %>%
  unlist %>% as.numeric

morpho_corr_long$resid_score <- resid

morpho_corr_long %>%
  ggplot(aes(y = resid, x = geno_sex, color = cluster)) +
  geom_jitter() +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               geom = "crossbar", width = 0.9, color = "black", fatten = 2) +
  facet_grid(trait~cluster, scales = "free") +
  scale_color_manual(values = whtstbk_palatte)



