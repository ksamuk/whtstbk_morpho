######################################################################
# Analysis of morphometric landmark data
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("nlme")
library("car")
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")
library("ggplot2")
library("ggthemes")

select<- dplyr::select #allows dplyr to use select when MASS package is loaded

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

## GEOMORPH Analysis ## 

# Create a 3D array with just the x/y coordinates 
morpho_sub <- morpho_df %>%
  filter(!is.na(x.1)) %>%
  filter(!is.na(cluster))

# select landmarks
# ignoring eye/pelvis (17-19) for now
landmark_df <- morpho_sub %>%
  select(matches("^[x,y]{1}\\.")) %>%
  select(-x.19, -y.19, -x.18, -y.18, -x.17, -y.17)

# create landmark array (16 landmarks in 2 columns each)
landmark_array <- arrayspecs(landmark_df, 16, 2)

# generalized grocrustes analysis (scale/align landmarks)
landmark_gpa <- gpagen(landmark_array) 

# confirm gpa is not insane
plot(landmark_gpa)

#landmark_gpa_2d <- two.d.array(landmark_gpa$coords) #2D Data frame of procrustes coordinates

######################################################################
# pca of landmark data
######################################################################

morpho_pca <- plotTangentSpace(landmark_gpa$coords, groups = factor(morpho_sub$cluster), verbose = TRUE)

sherrat_pca_plot(pc1 = 1, pc2 = 2, landmark_gpa = landmark_gpa, groups = factor(morpho_sub$cluster))

######################################################################
# direct comparisons of groups
######################################################################

gpa_array <- two.d.array(landmark_gpa$coords)

PC1 <- as.vector(morpho_pca$pc.scores[,1]) #get a vector of PC1 to use as a covariate

#without PC1 bending accounted for

prod_lm <- advanced.procD.lm(gpa_array~PC1, gpa_array~PC1 + morpho_sub$cluster,iter = 1000)

## Linear discriminant function analysis ##

lda.species <- lda(gpa_array, morpho_sub$cluster, CV = TRUE) 

# the percent of variance explained by the LD funcitons
prop.lda <- lda.species$svd^2/sum(lda.species$svd^2) 
prop.lda <- round(prop.lda*100)

# project the original values in lda space
plda <- predict(object = lda.species,
                newdata = gpa_array)

# data frame of projected data with species names
lda.project <- data.frame(cluster = morpho_sub$cluster, 
                          id = morpho_sub$id,
                          ld1 = plda$x[,1], 
                          ld2 = plda$x[,2])

# Plot LDA without "Both" locations
lda.project %>%
  ggplot(aes(color = cluster, x = ld1,y = ld2))+
  geom_point(size = 3) +
  labs(x = paste0("LD1 (", prop.lda[1], "%)"),
       y = paste0("LD2 (", prop.lda[2], "%)")) +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


## Cross Validation of LDA ##

# Method 1: take half of the values and create the LDA, 
#   then plot the rest to check model 

# LDA with CV function (automatic cross validation)
GPA.CW <- GPA
GPA.CW$species <- landmarks$species
GPA.CW$centroid = GPA.landmarks$Csize

GPA.CW<- GPA.CW %>%
  filter(species != "B") #remove "Both" category

landmarks.CW<- landmarks %>%
  filter(species != "B")

GPA.CW <- GPA.CW[-c(39,40)]

lda.species.CV <- lda(GPA.CW, landmarks.CW$species, CV=T)

# percent correct (accuracy of species prediction)
ct <- table(landmarks.CW$species, lda.species.CV$class)
diag(prop.table(ct, 1))

# total percent correct
sum(diag(prop.table(ct)))


# LDA validation Method 2: sample a random half of individuals
#   with equal representation of W and C 

# build a data frame with ID info and procrustes x/y
GPA.50<- data.frame(species = landmarks$species, 
                    population = landmarks$population, 
                    individual = landmarks$individual,
                    ID = landmarks$ID,
                    centroid = GPA.landmarks$Csize)

GPA.50<-merge(GPA.50, GPA, by.x="ID")

GPA.parse<- GPA.50 %>% 
  filter(species != "B") %>%
  group_by(species) 

# sample a random half of the data with equal W and C numbers
GPA.parse <- sample_n(GPA.parse, 403)
GPA.parse1<- GPA.parse[c(1:201, 403:604),]
GPA.parse2 <- GPA.parse[c(202:402,605:806),]

GPA.1 <- GPA.parse1[-c(1:4)]
GPA.2 <- GPA.parse2[-c(1:4)]


# run LDA for each 1/2 sample
lda.CV1 <- lda(GPA.1, as.character(GPA.parse1$species))
lda.CV2 <- lda(GPA.2, as.character(GPA.parse2$species))


#subtract scaling between the two halves and look at differences 
lda.CV <- as.vector(lda.CV1$scaling - lda.CV2$scaling)
plot(lda.CV) #weird LDA units 

# look at the plda and overplotting for each half

# project in lda space
plda.CV1 <- predict(object = lda.CV1,
                    newdata = GPA.1)

plda.CV2 <- predict(object = lda.CV2,
                    newdata = GPA.2)

# data frames of projected data with species names
lda.project.CV1 <- data.frame(species = GPA.parse1$species, 
                              population = GPA.parse1$population, 
                              individual = GPA.parse1$individual,
                              ID = GPA.parse1$ID,
                              ld1 = plda.CV1$x)

lda.project.CV2 <- data.frame(species = GPA.parse2$species, 
                              population = GPA.parse2$population, 
                              individual = GPA.parse2$individual,
                              ID = GPA.parse2$ID,
                              ld1 = plda.CV2$x)

# plot both halves on one plot... looks like they match 
ggplot(lda.project.CV1,aes(fill = species, x = LD1, alpha=0.2)) +
  geom_histogram() +
  geom_histogram(aes(x = lda.project.CV2$LD1, fill = species)) +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


## END cross validation ##


### Genetics VS Morphology ###

# extract centroids
landmark.centroids <- data.frame(population = landmarks$population, 
                                 id = landmarks$ID, 
                                 species = landmarks$species, 
                                 centroid = GPA.landmarks$Csize,
                                 sex= landmarks$sex)

# Plot centroid sizes by population 
landmark.centroids %>%
  ggplot(aes(x = population, y = centroid, color=species)) +
  geom_point()
  facet_grid(~ species, scale = "free", space = "free")


#read in genetic cluster data

cluster.dat <- read.table(file=list.files(pattern="structure",full.names = TRUE), stringsAsFactors = FALSE)
cluster.dat.reduced <- cluster.dat %>%
  filter(k.value.run == 2) %>%
  select(id,X1,X2,membership)

landmark.centroids.rename <- landmark.centroids
landmark.centroids.rename$id <- landmark.centroids$id %>% 
  gsub("2014","",.) %>% gsub("_","",.)

cluster.morph <- left_join(cluster.dat.reduced, landmark.centroids.rename)

cluster.morph <- cluster.morph  %>%
  filter(!is.na(population)) %>%
  filter(!is.na(sex))

cluster.morph$population <- recode(cluster.morph$population, "c('CPSE2014','CPN') = 'CP'")
cluster.morph$population <- gsub("2014", "", cluster.morph$population)
cluster.morph$membership<- recode(cluster.morph$membership, "c(1) = c('common'); else = 'white'")
cluster.morph$sex<- recode(cluster.morph$sex, "c('M') = c('Male'); else = 'Female'")



# Plot centroid size by population again
cluster.morph  %>%
  #filter(population %in% c("CL", "MH", "SF", "SR2014", "RT", "SH")) %>%
  ggplot(aes(x = population, y = centroid, color=factor(membership))) +
  geom_jitter(position = position_jitter(width = .1), size = 4) +
  theme_solarized(base_size = 24, light = FALSE)+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 24, face = 'bold'),
        panel.background = element_rect(fill = "black"),
        legend.title = element_blank(),
        panel.background = element_blank())+
  ylab("Body size")+
  xlab("Population")+
  scale_colour_solarized("blue")+
  facet_grid(sex~.)

ggsave(filename = "body_plot.pdf")
