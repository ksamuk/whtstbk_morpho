######################################################################
# Analysis of egg and testes data
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
# input data + formatting
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
repro_df <- read.csv("data/eggs_testes_data.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
repro_df <- repro_df %>%
  mutate(population = gsub("2014", "", population)) %>%
  mutate(id = paste0(population, individual) %>% gsub("CP[a-zA-Z]*", "CP", .))

# join meta data to morpho data

# prep meta data for join

meta_df <- meta_df %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex)

# join in cluster/geno sex and remove landmark data (explored in 01_body...)
repro_df  <- left_join(repro_df , meta_df, by = "id") %>%
  select(-membership, -note, -species, -individual, -population)

egg_df <- repro_df %>%
  filter(!is.na(egg.diameter.1)) %>%
  select(-matches("test"))

test_df <- repro_df %>%
  filter(!is.na(teste.length.1)) %>%
  select(-matches("egg"))

######################################################################
# egg data
######################################################################

egg_df <- gather(egg_df, key = egg_number, value = diameter, -egg.number, -weight..mg., -std.length, - depth.1, -id, -cluster, -sex, -sequenced., -geno_sex)

# calculate the mean egg diameter for each individual (there are ~10 each)
mean_dia <- egg_df %>%
  group_by(id) %>%
  summarise(mean_diameter = mean(diameter, na.rm = TRUE), sd_diameter = sd(diameter, na.rm = TRUE))

# join means back to egg data
egg_df <- left_join(mean_dia, egg_df) %>%
  select(id, sex, geno_sex, sequenced., cluster, mean_diameter, sd_diameter, egg.number, weight..mg., std.length, depth.1) %>%
  distinct()

# fix names
names(egg_df) <- c("id", "sex", "geno_sex", "sequenced", "cluster", "mean_diameter", "sd_diameter", 
                   "egg_number", "egg_weight", "std_length", "body_depth")

# group cbr and and mainland commons
egg_df$group <- ifelse(egg_df$cluster == "wht", "white", "common")

# filter outliers
egg_df <- egg_df %>%
  filter(egg_number < 150) %>%
  filter(id != "SR37") %>% 
  filter(!is.na(egg_number))%>% 
  filter(!is.na(std_length))

# residual egg number
resid <- egg_df %>%
  do(augment(lm(egg_number ~ std_length, data=.))) %>%
  ungroup %>%
  select(.resid) %>%
  unlist %>% as.numeric

egg_df$resid_egg_number <- resid 

######################################################################
# testis data
######################################################################

test_df <- repro_df %>%
  filter(!is.na(teste.length.1)) %>%
  filter(!is.na(teste.length.2)) %>%
  select(-matches("egg"))

# combine cbr and mainland commons
test_df$group <- ifelse(test_df$cluster == "wht", "white", "common")

# filter individual who was genotypically determined to be female, but is somehow in this dataset
test_df <- test_df %>%
  filter(!is.na(group)) %>%
  filter(id != "SR128")

# residual weight of testes
resid <- test_df %>%
  do(augment(lm(weight..mg. ~ std.length, data=.))) %>%
  select(.resid) %>%
  unlist %>% as.numeric

test_df$weight_resid <- resid

test_df %>%
  ggplot(aes(x = group, y = weight_resid, color = group, label = id, group = group))+
  #stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2)+
  geom_jitter(width = 0.3)+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black", position = position_dodge()) + 
  stat_summary(fun.y = mean, geom = "point", color = "black", position = position_dodge(), size = 2)
  
######################################################################
# plots
###################################################################### 
  
theme_all <- theme_base()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(color = "white"), 
        plot.background = element_rect(color = "white"),
        axis.ticks.x = element_blank())

# residual egg number
egg_df %>%
  ggplot(aes(x = group, y = resid_egg_number, color = group, label = id, group = group))+
  geom_jitter(width = 0.3)+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size =2)+
  theme_all

# egg weight
egg_df %>%
  ggplot(aes(x = group, y = egg_weight, color = group, label = id, group = group))+
  geom_jitter(width = 0.3)+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size =2)+
  theme_all

# testes mass
test_df %>%
  ggplot(aes(x = group, y = weight_resid, color = group, label = id, group = group))+
  #stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2)+
  geom_jitter(width = 0.3)+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  theme_all
  
