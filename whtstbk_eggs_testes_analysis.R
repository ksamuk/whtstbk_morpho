###Eggs/testes size analysis###


##load raw data##

######Google sheet method not working#######
# google sheets API
library("devtools")
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")

#devtools::install_github("jennybc/googlesheets")

library("googlesheets")

# list sheets and authenticate with server (will open browswer)
gs_ls()

# load in sheet from google sheets
gonad.dat.gs <- gs_title("East Coast Morphometrics - 2015- eggs-testes data")
gonad.dat <- data.frame(get_via_csv(gonad.dat.gs))

#########

gonad.dat<-read.csv("eggs_testes_data.csv")

##separate eggs and testes data##

gonad.dat<-subset(gonad.dat, egg.diameter.1!='NA' | teste.length.1 !='NA')

eggdat<-subset(gonad.dat, egg.diameter.1!='NA' & membership !='NA' & std.length !='NA') 
eggdat<-subset(eggdat, select= -c(teste.length.1, teste.length.2,
                                 teste.width.1, teste.width.2))

testedat<-subset(gonad.dat, teste.length.1!='NA' & membership !='NA' & std.length !='NA')
testedat<-subset(testedat, select= -c(egg.diameter.1, egg.diameter.2, egg.diameter.3,
                egg.diameter.4, egg.diameter.5, egg.diameter.6, egg.diameter.7,
                egg.diameter.8, egg.diameter.9, egg.diameter.10, egg.number))


##get average weight and size for eggs##

eggdat$ID <- paste(eggdat$population, eggdat$individual, sep="_")

eggdat <- transform(eggdat, avg.diameter = rowMeans(eggdat[,6:15], na.rm=TRUE))

eggdat$avg.weight<-eggdat$weight/10


#graph avg egg size and weight by body size (standard length) 

library("ggplot2")
library("ggthemes")

ggplot(eggdat, aes(x=std.length, y=avg.diameter, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length",
    y = "Avg. egg diameter",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

ggplot(eggdat, aes(x=std.length, y=avg.weight, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) + 
  labs(
    x = "Standard Length",
    y = "Avg. egg weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


#boxplots; standard length controlled through residuals of linear regression

  #diameter
egg.diam.lm <- lm(avg.diameter~std.length, data=eggdat)
eggdat$egg.diam.resid <- residuals(egg.diam.lm)

  #Diameter by Membership
ggplot(eggdat, aes(x=membership, egg.diam.resid, colour=factor(membership))) +
  geom_boxplot() +
  labs(
    x = "Membership",
    y = "Residuals for avg. egg diameter",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


  #weight
egg.weight.lm <- lm(avg.weight~std.length, data=eggdat)
eggdat$egg.weight.resid <- residuals(egg.weight.lm)

  #weight by membership
ggplot(eggdat, aes(x=membership, egg.weight.resid, colour=factor(membership))) +
  geom_boxplot() + 
  labs(
    x = "Membership",
    y = "Residuals for avg. egg weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


##hypothesis testing: difference in egg size & weight by membership

library("nlme")

#size
egg.test<-lme(avg.diameter ~ std.length * membership, data=eggdat, random= ~1|population, na.action=na.omit) 
summary(egg.test) 

anova(egg.test) # type I anova                         

  #linear model fit anova 
egg.diam.anova<-aov(eggdat$avg.diameter ~ eggdat$membership)
summary(egg.diam.anova)


#weight
egg.weight.test<-lme(avg.weight ~ std.length * membership, data=eggdat, random= ~1|population, na.action=na.omit) 
summary(egg.weight.test) 

anova(egg.weight.test) # type I anova                         

  #linear model fit anova 
egg.weight.anova<-aov(eggdat$avg.weight ~ eggdat$membership)
summary(egg.weight.anova)


##get average weight, size, and variance for testes##

testedat$ID <- paste(testedat$population, testedat$individual, sep="_")

testedat <- transform(testedat, avg.size = (rowMeans(testedat[,6:7], na.rm=TRUE)) *
                             rowMeans(testedat[,8:9], na.rm=TRUE))

testedat$avg.weight<-testedat$weight/2

testedat$variance <- abs((testedat[,6]*testedat[,8]) - (testedat[,7]*testedat[,9]))


#graph avg teste size, weight, and variance by body size (standard length) 

library("ggplot2")

#size
ggplot(testedat, aes(x=std.length, y=avg.size, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length",
    y = "Avg. teste size (length*width)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

#weight
testedat %>%
  filter(avg.weight<1) %>%
  filter(!(membership == 1 & std.length < 4)) %>%
ggplot(aes(x=std.length, y=avg.weight, color=factor(membership))) +
  geom_point(size = 4, na.rm = T) +
  stat_smooth(method=lm, na.rm=T, size = 2, se = FALSE) +
  labs(
    x = "Standard Length",
    y = "Avg. teste weight",
    color = "Membership") +
  theme_hc(base_size = 16) +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


##checking other patterns and relationships between variables 

#variance in size
ggplot(testedat, aes(x=std.length, y=variance, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T)

#size vs weight relationship
ggplot(testedat, aes(x=avg.weight, y=avg.size, color=factor(membership))) +
  geom_point(na.rm=T) +
   stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Avg. teste weight",
    y = "Avg. teste size (length*width)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

#different way of calculating average size...  

testedat$avg.size.2 = (((testedat[,6]) *(testedat[,8])) +
                        ((testedat[,7])*(testedat[,9])))/2

  #size 2
ggplot(testedat, aes(x=std.length, y=avg.size.2, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length",
    y = "Avg. teste size (length*width)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))



#boxplots: size and weight differences

  #size lm to get residuals
teste.size.lm <- lm(avg.size ~ std.length, data=testedat, na.action=na.exclude)
testedat$teste.size.resid <- residuals(teste.size.lm)

# size by membership
testedat %>%
  filter(teste.size.resid > 0.3) %>%
ggplot(aes(x=membership, teste.size.resid, colour=factor(membership))) +
  geom_boxplot(na.rm=T) + 
  labs(
    x = "Membership",
    y = "Residuals for avg. teste size",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# weight lm for residuals
teste.weight.lm <- lm(avg.weight~std.length, data=testedat, na.action=na.exclude)
testedat$teste.weight.resid <- residuals(teste.weight.lm)

  #weight by membership
ggplot(testedat, aes(x = membership, y = avg.weight, colour=factor(membership))) +
  xlab("Membership") + ylab("Residuals for  avg. teste weight") +
  geom_boxplot(na.rm=T) +
  labs(
    x = "Membership",
    y = "Residuals for avg. teste weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


##hypothesis testing: difference in teste size & weight by membership

library("nlme")

#size differences 
teste.test <- testedat %>%
filter(!(membership == 1 & std.length < 4)) %>%
filter(avg.weight < 1.0) %>%
lm(avg.weight ~ std.length * membership, data=., na.action=na.omit) %>%
summary() 

visreg2d(teste.test)

anova(teste.test)  #type I anova                        

  #fit linear model anova
teste.size.anova<-aov(testedat$avg.size ~ testedat$membership)
summary(teste.size.anova)


#weight differences

teste.weight.test<-lme(avg.weight ~ std.length * membership, data=testedat, random= ~1|population, na.action=na.omit) 
summary(teste.weight.test) 

anova(teste.weight.test)  # type I anova                        

  #fit linear model anova
teste.weight.anova<-aov(testedat$avg.weight ~ testedat$membership)
summary(teste.weight.anova)




##get other teste data: filter by populations which have been sequenced & SL 4-5cm

data<-read.csv("east_coast_morphometrics_2015.csv")
filtered.data<-filter(data, sequenced.==1)
filtered.data<-filter(filtered.data, species!='B')
filtered.data<-filter(filtered.data, sex=='M')
filtered.data<-filter(filtered.data, std.length>4 & std.length<5)

theme_black=function(base_size=12,base_family="") {
  theme_grey(base_size=base_size,base_family=base_family) %+replace%
    theme(
      # Specify axis options
      axis.line=element_blank(), 
      axis.text.x=element_text(size=base_size*0.8,color="white",
                               lineheight=0.9,vjust=1), 
      axis.text.y=element_text(size=base_size*0.8,color="white",
                               lineheight=0.9,hjust=1), 
      axis.ticks=element_line(color="white",size = 0.2), 
      axis.title.x=element_text(size=base_size,color="white",vjust=1), 
      axis.title.y=element_text(size=base_size,color="white",angle=90,
                                vjust=0.5), 
      axis.ticks.length=unit(0.3,"lines"), 
      axis.ticks.margin=unit(0.5,"lines"),
      # Specify legend options
      legend.background=element_rect(color=NA,fill="black"), 
      legend.key=element_rect(color="white", fill="black"), 
      legend.key.size=unit(1.2,"lines"), 
      legend.key.height=NULL, 
      legend.key.width=NULL,     
      legend.text=element_text(size=base_size*0.8,color="white"), 
      legend.title=element_text(size=base_size*0.8,face="bold",hjust=0,
                                color="white"), 
      legend.position="right", 
      legend.text.align=NULL, 
      legend.title.align=NULL, 
      legend.direction="vertical", 
      legend.box=NULL,
      # Specify panel options
      panel.background=element_rect(fill="black",color = NA), 
      panel.border=element_rect(fill=NA,color="white"), 
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank(), 
      panel.margin=unit(0.25,"lines"),  
      # Specify facetting options
      strip.background=element_rect(fill="grey30",color="grey10"), 
      strip.text.x=element_text(size=base_size*0.8,color="white"), 
      strip.text.y=element_text(size=base_size*0.8,color="white",
                                angle=-90), 
      # Specify plot options
      plot.background=element_rect(color="black",fill="black"), 
      plot.title=element_text(size=base_size*1.2,color="white"), 
      plot.margin=unit(c(1,1,0.5,0.5),"lines")
    )
}

