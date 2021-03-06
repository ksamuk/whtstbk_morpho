#### F2 white stickleback body size and melanocyte analysis #### 


library(ggplot2)
library(dplyr)

### Initial loading 

data <- read.csv("F2 white stickleback morphometrics.csv")

data$family <- as.factor(data$family)

data$melanocytes <- as.integer(data$melanophores)



### Body Size: Length and Depth

morphodat <- data[,c(2:4, 8:11)]

## Graphs 

# Boxplot body depth by family
ggplot(morphodat, aes(x=family, y=body.depth)) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "F2 Family",
    y = "Body Depth from first dorsal spine (cm)") +
  theme_classic() 

# Boxplot standard length by family
morphodat %>%
  filter(family != 7) %>%
  ggplot(aes(x=family, y=standard.length, fill = family, color = "black")) +
    geom_boxplot(na.rm=T) +
    labs(
      x = "F2 Family",
      y = "Standard Length (cm)") +
    theme_solarized(light = FALSE)+
    scale_color_solarized()
    



## Hypothesis testing  

# One-way anova for average population std. length by species 
std.length.lm<-lm(morphodat$standard.length ~ morphodat$family)
summary(std.length.lm)
anova(std.length.lm)

# One-way anova for average population body depth by species 
depth.lm <- lm(morphodat$body.depth ~ morphodat$family)
summary(depth.lm)
anova(depth.lm)


### melanocyte counts analysis ### 

## Graphs 

# melanocyte counts by family
data %>%
  filter(family != 7) %>%
  mutate(family = factor(family)) %>%
  ggplot(aes(x = family, y = melanocytes, fill = factor(family), color = 1)) +
  geom_boxplot()+
    labs(
      x = "F2 Family",
      y = "Number of melanocytes") +
  theme_solarized(base_size = 24, light = FALSE)



## Hypothesis testing: melanocyte counts by family 
melan.lm <- lm(data$melanocytes ~ data$family)
summary(melan.lm)
anova(melan.lm)

ggplot(data, aes(x = melanocytes, y = body.depth, color = family))+
  geom_point(size = 4)

