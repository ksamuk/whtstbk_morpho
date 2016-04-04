

# quick join and format of 2012 data 

install.packages("geomorph")
library("geomorph")
library("dplyr")

tps_file <- "2012_data/whtstbk_2012_shapedata.TPS"


# first, extract scales
tps_con <- file(tps_file, open = "r")
tps_raw <- readLines(tps_con)
close(tps_con)

max_num <- 31
x <- seq_along(tps_raw)
tps_raw <- split(tps_raw, ceiling(x/max_num))

get_id_scale <- function(x){
  
  id <- x[29] %>% gsub("IMAGE=|.jpg", "", .)
  scale <- x[31] %>% gsub("SCALE=", "", .) %>% as.numeric
  
  data.frame(id, scale)
}

scale_df <- lapply(tps_raw, get_id_scale)
scale_df <- bind_rows(scale_df)

# next, get shape data
tps_raw <- split(tps_raw, ceiling(seq_along(tps_raw)/31))
tps_dat <- readland.tps("2012_data/whtstbk_2012_shapedata.TPS",specID = "imageID")

# don't want to do this, scales fish :P
GPA.landmarks <- gpagen(tps_dat) 
GPA.2D <- two.d.array(tps_dat)


# for SL, just want snount and caudal peduncle ('tail')

get_nose_tail <- function(x){
  
  id <- x[29] %>% gsub("IMAGE=|.jpg", "", .)
  scale <- x[31] %>% gsub("SCALE=", "", .) %>% as.numeric
  
  data.frame(id, scale)
}

#nose [17,]
#tail [1,]

# standard length for 2012 data

#noses (x,y)
nose <- data.frame(GPA.2D[,17:18])
names(nose) <- c("nose_x", "nose_y")
caud <- data.frame(GPA.2D[,1:2])
names(caud) <- c("caud_x", "caud_y")

sl_df <- data.frame(id = row.names(caud), nose, caud)

sl_df <- sl_df %>%
  mutate(sl = sqrt(((caud_x - nose_x)^2) + ((caud_y - nose_y))^2))

sl_df <- left_join(sl_df, scale_df)

#sl_df$std.length <- (sl_df$sl*sl_df$scale)
sl_df$std.length <- (sl_df$sl)

sl_df <- sl_df %>% 
  select(id, std.length)

sl_df <- sl_df %>%
  mutate(pop = gsub("[^A-Z]", "", id)) %>%
  mutate(id = gsub("[^0-9]", "", id)) %>%
  mutate(year = 2012) %>%
  select(pop, year, id, std.length)

boxplot(sl_df$std.length~sl_df$pop)
  

write.table(sl_df, "2012_data/whtstbk_2012_standard_length.txt", row.names = FALSE, quote = FALSE )

