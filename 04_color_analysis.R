######################################################################
# Analysis of photographic color data
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("grDevices")
library("ggplot2")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("agrmt")
library("zoo")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#77AB43", "#008FD5", "#FFFFFF")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
color_df <- read.csv("data/whtstbk_color_data.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
color_df <- color_df %>%
  mutate(id = paste0(population, individual) %>% 
           gsub("CP[a-zA-Z]*", "CP", .) %>% 
           gsub("-2014|2014", "", .))

# join meta data to morpho data

# prep meta data for join

meta_df <- meta_df %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex)

color_df <- left_join(color_df, meta_df, by = "id")

######################################################################
# summarize color histograms
######################################################################

# genius function from stackoverflow for finding local maxima
# source: http://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset

argmax <- function(x, y, w=2, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

#red
tmp_r <- color_df %>%
  select(matches("R_"), id) %>% 
  gather(key = r_bin, value = r_freq, matches("R_"))

tmp_r <- tmp_r %>%
  mutate(r_value = as.numeric(gsub("R_", "", r_bin))) %>%
  group_by(id) %>%
  summarise(r_mean = sum(r_freq*r_value)/sum(r_freq), r_mode_1 = argmax(y = r_freq, x = r_value, w = 10)$x[1], 
            r_mode_2 = ifelse(!is.na(argmax(y = r_freq, x = r_value, w = 20)$x[2]), argmax(y = r_freq, x = r_value, w = 10)$x[2], NA))

#blue
tmp_b <- color_df %>%
  select(matches("B_"), id) %>% 
  gather(key = b_bin, value = b_freq, matches("B_"))

tmp_b <- tmp_b %>%
  mutate(b_value = as.numeric(gsub("B_", "", b_bin))) %>%
  group_by(id) %>%
  summarise(b_mean = sum(b_freq*b_value)/sum(b_freq), b_mode_1 = argmax(y = b_freq, x = b_value, w = 10)$x[1], 
            b_mode_2 = ifelse(!is.na(argmax(y = b_freq, x = b_value, w = 20)$x[2]), argmax(y = b_freq, x = b_value, w = 10)$x[2], NA))

#green
tmp_g <- color_df %>%
  select(matches("G_"), id) %>% 
  gather(key = g_bin, value = g_freq, matches("G_"))

tmp_g <- tmp_g %>%
  mutate(g_value = as.numeric(gsub("G_", "", g_bin))) %>%
  group_by(id) %>%
  summarise(g_mean = sum(g_freq*g_value)/sum(g_freq), g_mode_1 = argmax(y = g_freq, x = g_value, w = 10)$x[1], 
            g_mode_2 = ifelse(!is.na(argmax(y = g_freq, x = g_value, w = 20)$x[2]), argmax(y = g_freq, x = g_value, w = 10)$x[2], NA))

# join all color data

col_df <- left_join(tmp_r, tmp_g)
col_df <- left_join(col_df, tmp_b)
col_df <- left_join(col_df, meta_df)

col_df %>%
  group_by(id) %>%
  mutate(luminance_mean = (r_mean + g_mean + b_mean)) %>%
  mutate(luminance_mode = (r_mean + g_mean + b_mean)) %>%
  rename(rgb_r_mean = r_mean, rgb_g_mean = g_mean, rgb_b_mean = b_mean) %>%
  rename(rgb_r_mode = r_mode_1, rgb_g_mode = g_mode_1, rgb_b_mode = b_mode_1) %>%
  select(id, matches("rgb|lum")) %>%
  write.table(file = "data/collated/raw_color_data.txt", quote = FALSE, row.names = FALSE)

######################################################################
# plotting all color data
######################################################################

#color channels (R–G)/(R+G) and (G–B)/(G+B) (Endler 2012)

col_lng <- col_df %>%
  mutate(luminance = (r_mean + g_mean + b_mean)) %>%
  mutate(r_g_diff = (r_mean-g_mean)/(r_mean+g_mean)) %>%
  mutate(g_b_diff = (g_mean-b_mean)/(g_mean+b_mean)) %>%
  gather(key = trait, value = value, -id, -cluster, -geno_sex) %>%
  mutate(group = ifelse(cluster == "wht", "white", "common"))

col_lng %>%
  filter(!is.na(cluster)) %>%
  filter(trait == "luminance") %>%
  ggplot(aes(x = group, y = value))+
  geom_jitter(width = 0.5, color = "grey")+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  theme_base()

col_lng %>%
  filter(!is.na(cluster)) %>%
  ggplot(aes(x = group, y = value))+
  geom_jitter(width = 0.5, color = "grey")+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  facet_wrap(~trait, scale = "free")+
  theme_base()

######################################################################
# tests of luminance difference sex*species
######################################################################

col_lm <- col_df %>%
  mutate(group = ifelse(cluster == "wht", "white", "common")) %>%
  mutate(luminance = (r_mean + g_mean + b_mean)) %>%
  filter(!is.na(cluster)) %>%
  filter(!is.na(geno_sex)) %>%
  lm(luminance~geno_sex + group + group:geno_sex, data = .)

anova(col_lm)

visreg(col_lm)

col_df %>%
  mutate(group = ifelse(cluster == "wht", "white", "common")) %>%
  mutate(luminance = (r_mean + g_mean + b_mean)) %>%
  ggplot(aes (x = group, y = luminance))+
  geom_jitter(color = "grey", width = 0.5)+
  stat_summary(fun.data = mean_cl_boot, size = 1, geom = "errorbar", width = 0.2, color = "black") + 
  stat_summary(fun.y = mean, geom = "point", color = "black", size = 2)+
  theme_base()
