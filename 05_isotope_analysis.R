######################################################################
# Analysis of isotopic abundance data
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
library("proto")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#77AB43", "#008FD5", "#BDBDBD")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
iso_df <- read.csv("data/isotope_data.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
meta_df <- meta_df %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex)

iso_df  <- left_join(iso_df, meta_df, by = "id")

iso_df <- iso_df %>%
  mutate(group = ifelse(cluster == "wht", "white", "common")) %>%
  mutate(pop = gsub("[^A-Z]*", "", as.character(id))) %>%
  mutate(cn.ratio = C.amount / N.amount)


iso_df$d15N.resid <- lm(data = iso_df, d15N ~ pop*geno_sex, na.action = "na.exclude") %>% residuals
iso_df$d13C.resid <- lm(data = iso_df, d13C ~ pop*geno_sex, na.action = "na.exclude") %>% residuals

iso_df %>%
  select(id, d13C, d15N, d13C.resid, d15N.resid, cn.ratio) %>%
  rename(isotope_d13C = d13C, isotope_d15N = d15N, 
         isotope_d13C_resid = d13C.resid, isotope_d15N_resid = d15N.resid, isotope_CN_ratio = cn.ratio) %>%
  write.table(file = "data/collated/raw_isotope_data.txt", quote = FALSE, row.names = FALSE)

#raw
iso_df %>%
  ggplot(aes(x = d13C, y = d15N, label = id, color = group, shape = geno_sex)) +
  geom_point(size = 4) +
  theme_base()+
  scale_color_manual(values = c("#008FD5", "#BDBDBD"))

# controlled for pop and cn.ratio

iso_df %>%
  filter(!is.na(group)) %>%
  ggplot(aes(x = d13C.resid, y = d15N.resid, label = id, color = group)) +
  geom_vline(xintercept = 0, color = "grey")+
  geom_hline(yintercept = 0, color = "grey")+
  geom_point(size = 4) +
  theme_base()+
  scale_color_manual(values = c("#008FD5", "#BDBDBD"))

iso_df %>%
  filter(!is.na(group)) %>%
  ggplot(aes(x = group, y = d13C.resid, color = group)) +
  geom_jitter(width = 0.5)+
  stat_summary(fun.data = mean_cl_normal, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  theme_base()+
  scale_color_manual(values = c("#008FD5", "#BDBDBD"))

iso_df %>%
  filter(!is.na(group)) %>%
  ggplot(aes(x = group, y = d15N.resid, color = group)) +
  geom_jitter(width = 0.5)+
  stat_summary(fun.data = mean_cl_normal, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  theme_base()+
  scale_color_manual(values = c("#008FD5", "#BDBDBD"))
  

######################################################################
# plot of isotopic data
######################################################################

lm(data = iso_df, d13C ~ pop + group) %>% plot
lm(data = iso_df, d13C ~ pop + group + 0, na.action = "na.omit") %>% summary
lm(data = iso_df, d13C ~  pop + group + 0, na.action = "na.omit") %>% visreg::visreg()

group <- as.numeric(as.factor(iso_df$group))

# kw test
iso_df %>%
  filter(!is.na(group)) %>%
  mutate(group = as.numeric(as.factor(group))) %>%
  kruskal.test(d15N.resid ~ group, data = .)

# kw test
iso_df %>%
  filter(!is.na(group)) %>%
  mutate(group = as.numeric(as.factor(group))) %>%
  kruskal.test(d13C.resid ~ group, data = .)

lm(data = iso_df, d15N ~ group) %>% anova
tmp <- lm(data = iso_df, d15N ~  pop + group, na.action = "na.omit") 
  
  
visreg2d(tmp, xvar = "pop", yvar = "group")

sqrt(iso_df$d15N) %>% hist

#Analysis of Variance Table

# Response: d15N
# Df  Sum Sq Mean Sq F value  Pr(>F)    
# pop         4 146.209  36.552 71.2681 < 2e-16 ***
#   cn.ratio    1   2.181   2.181  4.2526 0.04163 *  
#   cluster     1   1.941   1.941  3.7852 0.05436 .  
# Residuals 106  54.366   0.513  

lm(data = iso_df, d15N ~ pop*group) %>% visreg::visreg()

# Analysis of Variance Table
# 
# Response: d13C
# Df Sum Sq Mean Sq F value  Pr(>F)    
# pop         4 77.712 19.4280 39.1065 < 2e-16 ***
#   cn.ratio    1  2.661  2.6615  5.3573 0.02256 *  
#   cluster     1  0.056  0.0564  0.1134 0.73693    
# Residuals 106 52.660  0.4968    

