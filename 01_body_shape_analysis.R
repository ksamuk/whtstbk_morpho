######################################################################
# Analysis of morphometric landmark data
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("nlme")
library("car")
library("broom")
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select<- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#77AB43", "#008FD5", "#FFFFFF")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
morpho_df <- read.csv("data/corrected_landmark_data.csv")

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

morpho_df <- left_join(morpho_df, meta_df, by = "id")

######################################################################
# transform data + generalized procrustes
######################################################################

morpho_df <- data.frame(morpho_df, csize = landmark_gpa$Csize)

# Create a 3D array with just the x/y coordinates 

outlier_individuals <- c("GC15", "MH9")

morpho_sub <- morpho_df %>%
  filter(!is.na(x.1)) %>%
  filter(!is.na(cluster)) %>%
  filter(!(id %in% outlier_individuals)) %>%
  filter(!grepl("GC", id))
  

# select landmarks
# ignoring eye/pelvis (17-19) for now
landmark_df <- morpho_sub %>%
  select(matches("^[x,y]{1}\\.")) %>%
  select(-x.19, -y.19, -x.18, -y.18, -x.17, -y.17, -x.10, -y.10, -x.8, -y.8)

# create landmark array (16 landmarks in 2 columns each)
landmark_array <- arrayspecs(landmark_df, 14, 2)

# generalized grocrustes analysis (scale/align landmarks)
landmark_gpa <- gpagen(landmark_array) 

# confirm gpa is not insane
plot(landmark_gpa)

#landmark_gpa_2d <- two.d.array(landmark_gpa$coords) #2D Data frame of procrustes coordinates

# add in scaling factor into morpho_sub
a <- landmark_gpa$Csize

######################################################################
# pca of landmark data
######################################################################

morpho_pca <- plotTangentSpace(landmark_gpa$coords, groups = factor(morpho_sub$cluster), verbose = TRUE, label = morpho_sub$id)

sherrat_pca_plot(pc1 = 2, pc2 = 3, landmark_gpa = landmark_gpa, groups = factor(morpho_sub$cluster), label = morpho_sub$id)

# screen plot of pcs
pca_scree <- morpho_pca$pc.summary$importance[2,]
names(pca_scree) <- 1:length(pca_scree)
barplot(pca_scree, cex.lab = 0.1)

#parcoord(morpho_pca$pc.scores[,1:16], col=morpho_sub$cluster)

######################################################################
# allometry
######################################################################

# create a data frame for procD operations
gdf <- geomorph.data.frame(landmark_gpa, cluster = morpho_sub$cluster, sex = morpho_sub$geno_sex) # geomorph data frame

# simple test of allometry
whtstbk_allometry <- procD.allometry(coords~Csize, f2 = NULL, f3=NULL, 
                                  logsz = TRUE, data=gdf, iter=499) 

# allow variation between groups
whtstbk_allometry <- procD.allometry(coords~Csize, ~sex*cluster, 
                                  logsz = TRUE, data=gdf, iter=499, RRPP=TRUE)

# anova table
summary(whtstbk_allometry)

# a more sane plot
pca_df <- data.frame(id = morpho_sub$id, sex = morpho_sub$sex, csize = landmark_gpa$Csize,
                     cluster = morpho_sub$cluster, morpho_pca$pc.scores[,1:6])

pca_df_long <- gather(pca_df, key = pc, value = score, -id, -sex, -csize, -cluster)

# how do all the pcs scale with body size?
pca_df_long %>%
  ggplot(aes(y = score, x = csize, color = cluster))+
  geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "black", size = 1, se = FALSE)+
  facet_wrap(~pc)+
  theme_base()+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))

# how do all the pcs compare between species?
pca_df_long %>%
  ggplot(aes(y = score, x = cluster))+
  geom_jitter(color = "grey", width = 0.5)+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width = 0.9, color = "black", fatten = 2)+
  facet_wrap(~pc)+
  theme_base()


# create residual pcs (regress on csize)

resid <- pca_df_long %>%
  group_by(cluster, pc) %>%
  do(augment(lm(score ~ csize, data=.))) %>%
  ungroup %>%
  select(.resid) %>%
  unlist %>% as.numeric

pca_df_long$resid_score <- resid

# plot residuals vs. size (sanity check)
pca_df_long %>%
  ggplot(aes(y = resid_score, x = csize, color = cluster))+
  geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "black", size = 1, se = FALSE)+
  facet_wrap(~pc)+
  theme_base()+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))

# how do the residuals pcs compare between species?
pca_df_long %>%
  ggplot(aes(y = resid_score, x = cluster, color = cluster))+
  geom_jitter(width = 0.5, size = 2)+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width = 0.9, color = "black", fatten = 4)+
  facet_wrap(~pc)+
  theme_base()+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))

# write collated output
pca_df_long %>% 
  select(-score) %>% 
  spread(key = pc, value = resid_score) %>%
  write.table(., file = "data/collated/corrected_pc_scores.txt", quote = FALSE, row.names = FALSE)

pca_df_long %>% 
  select(-resid_score) %>% 
  spread(key = pc, value = score) %>%
  select(-sex, -cluster) %>%
  write.table(., file = "data/collated/raw_pc_scores.txt", quote = FALSE, row.names = FALSE)

######################################################################
# Anova: do pcs differ between groups?
######################################################################

pca_df_resid <- pca_df_long %>%
  select(-score) %>%
  spread(key = pc, value = resid_score)

lm(PC2 ~ sex + cluster, data = pca_df_resid) %>% anova
lm(PC3 ~ sex + cluster, data = pca_df_resid) %>% anova
lm(PC4 ~ sex + cluster, data = pca_df_resid) %>% anova
lm(PC5 ~ sex + cluster, data = pca_df_resid) %>% anova
lm(PC6 ~ sex + cluster, data = pca_df_resid) %>% anova


#mano <- manova(cbind(PC1, PC2, PC3, PC4, PC5, PC6) ~ sex + cluster, data = pca_df_resid)
#summary(mano)
#Anova(mano)

# nope


