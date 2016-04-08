######################################################################
# Analysis of gill rakers and plate counts (2012)
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("ggplot2")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#BDBDBD")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2012 (rakers, armor)
raker_df <- read.csv("data/2012_rakers_armor.csv")

# standard length for 2012 data
sl_df <- read.csv("data/2012_std_length.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
meta_df <- meta_df %>%
  filter(year == 2012) %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex) %>%
  mutate(id = gsub("whtstbk_gbs_2012_brds_", "", id))

raker_df  <- full_join(raker_df, meta_df, by = "id") %>%
  filter(!is.na(long.dorsal)) %>%
  left_join(sl_df)

wht_pop <- c("SR", "PP", "SF")

raker_df <- raker_df %>%
  mutate(pop = gsub("[0-9]*", "", id)) %>%
  mutate(species = ifelse(pop %in% wht_pop, "white", "common"))

######################################################################
# write collated file
######################################################################

raker_df %>%
  select(id, species, sum.long, short.sum, plate.count, dorsal.spine1, dorsal.spine2, std.length) %>%
  rename(raker_long_count = sum.long, raker_short_count = short.sum, 
         plate_count = plate.count, spine_dorsal1 = dorsal.spine1, 
         spine_dorsal2 = dorsal.spine2, std_length = std.length) %>%
  write.table(file = "data/collated/raw_raker_data.txt", quote = FALSE, row.names = FALSE)

######################################################################
# exploration
######################################################################

raker_lng <- raker_df %>%
  select(id, cluster, geno_sex, pop, species, std.length, plate.count, short.sum, sum.long) %>%
  gather(key = trait, value = value, -id, -cluster, -geno_sex, -std.length, -pop, -species)

######################################################################
# plots
######################################################################
  
raker_lng %>%
  ggplot(aes(x = species, y = value, color = species, label = id)) +
  geom_jitter(size = 2) +
  stat_summary(fun.data = mean_cl_normal, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  facet_wrap(~trait, scale = "free")+
  theme_base()+
  scale_color_manual(values = c("#008FD5", "#BDBDBD"))

raker_lng %>%
  ggplot(aes(x = std.length, y = value, color = species, label = id)) +
  geom_point(size = 2) +
  facet_wrap(~trait, scale = "free")+
  theme_base()+
  scale_color_manual(values = c("#008FD5", "#BDBDBD"))

######################################################################
# linear models
######################################################################

raker_df %>%
  glm(data = ., plate.count~species, family = poisson) %>%
  anova(., test = "Chisq")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)
#NULL                      128     4.7733         
#species  1  0.16189       127     4.6114   0.6874

raker_df %>%
  glm(data = ., short.sum~species, family = poisson) %>%
  anova(., test = "Chisq")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)
#NULL                      128     8.9833         
#species  1  0.18576       127     8.7975   0.6665

raker_df %>%
  glm(data = ., sum.long~species, family = poisson) %>%
  anova(., test = "Chisq")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)
#NULL                      128     13.857         
#species  1  0.51254       127     13.344    0.474

