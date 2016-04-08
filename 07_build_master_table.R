######################################################################
# Combine all phenotypic data into master table :o
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("dplyr")
library("tidyr")
library("broom")
library("ggplot2")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#BDBDBD")

######################################################################
# input data
######################################################################

# get file names for files to be combined
pheno_files <- list.files("data/collated", pattern = "raw", full.names = TRUE)

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

meta_df <- meta_df %>%
  rename(geno_sex = sex) %>%
  filter(!(pop %in% c("DK", "LC")))

meta_df$sequenced <- 1

######################################################################
# joining data frames
######################################################################

# all but raker data are individuals from 2014 (so can join on ids directly)
files_2014 <- grep("raker", pheno_files, invert = TRUE, value = TRUE)

pheno_dfs <- lapply(files_2014, read.table, header = TRUE, stringsAsFactors = FALSE)

# initial joint of meta_df to 2014 data

# initialize pheno_df
pheno_df <- meta_df

for (i in pheno_dfs){
  
  pheno_df <- full_join(pheno_df, i, by = c("id"))
}

# fix missing and year info and remove duplicated columns
pheno_df <- pheno_df %>%
  mutate(geno_sex = as.character(geno_sex)) %>%
  mutate(sex = as.character(sex)) %>%
  mutate(year = ifelse(is.na(year), 2014, year))

# fix sex info
pheno_df$joint_sex <- paste0(pheno_df$geno_sex, pheno_df$sex) %>% 
  gsub("NA*", "", .) %>%
  substr(1,1) %>%
  ifelse(. == "", NA, .)

# scrub 2012 ids of weird slug
pheno_df <- pheno_df %>%
  mutate(id = gsub("whtstbk_gbs_2012_brds_", "", id))

# join in 2012 raker data (carefully)
raker_file <- grep("raker", pheno_files, value = TRUE)
raker_df <- read.table(raker_file, header = TRUE, stringsAsFactors = FALSE) 

raker_df$species <- ifelse(raker_df$species == "common", "cmn", "wht")
raker_df$year <- 2012

pheno_df <- full_join(pheno_df, raker_df)

# add in population codes
pheno_df <- pheno_df %>%
  mutate(pop = gsub("[^A-Z]*", "", id) %>% substr(1,2))
  
write.table(pheno_df, "data/pheno_df_master.txt", quote = FALSE, row.names = FALSE)

######################################################################
# joining data frames
######################################################################



