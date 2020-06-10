#!/usr/bin/env Rscript
#
## ---------------------------
##
## Script name: all_datasets_plot.R
##
## Purpose of script: plot all datasets together
## 
## Author: Leonidas Liakos
##
## Date Created: 01/06/2020
##
##
## ---------------------------
##
##
##
## ---------------------------

library(readr)
library(ggplot2)
library(config)

source(here::here("R", 'myfunctions.R'))



####====== Read settings ==================== ####
myconfig <- config::get(file = here::here("R", "config_active.yml"),config ="default")
name<-myconfig$ACTIVE
cfg <- config::get(file = here::here("R", "config.yml"),config = name)
####========================================= ####


####====== Read csv ==================== ####
csv_files <- list.files(here::here("output", name), "*dataset*.*csv", full.names = T)

DO_dataframes                 <- readr::read_csv(here::here("output", name,sprintf("%s_dataframe_datasets_DO_year_%s.csv",name, cfg$YEAR)))
median_shift_roman_dataframes <- readr::read_csv(here::here("output", name,sprintf("%s_dataframe_datasets2_roman_%s.csv", name,cfg$YEAR)))
cao_dataframes                <- readr::read_csv(here::here("output", name,sprintf("%s_dataframe_datasets2_cao_%s.csv", name,cfg$YEAR)))


####====== merge csvs in tibble ==================== ####
mytbl <- DO_dataframes %>%  
  dplyr::left_join(median_shift_roman_dataframes, by=c("mydates")) %>% 
  dplyr::left_join(cao_dataframes, by=c("mydates")) %>% 
  dplyr::select(mydates, name,original.x,median_shift_corrected.x, roman_corrected,cao_corrected, meanDNB_corrected ) %>% 
  `colnames<-`(c("mydates", "name", "original","median_shift", "roman","cao", "DO"))


# convert to long format for ggplot
data_long <- tidyr::gather(mytbl, dataset, value, original:DO) %>%
    dplyr::filter(!is.na(value))



centroid <-
  as(ext, "SpatialPolygons") %>% 
  st_as_sf() %>% 
  st_centroid() %>% 
  st_set_crs(2100) %>% 
  st_transform("+init=epsg:4326") %>% 
  st_coordinates() 
colnames(centroid) <- c("lon", "lat")



days365 <- seq(min(mytbl$mydates),max(mytbl$mydates),by="day")
moon_ill <- data.frame(li = sapply(
    days365,
    FUN = function(x) {
      janiczek(
        lon = centroid[, "lon"],
        lat = centroid[, "lat"],
        year = lubridate::year(x),
        month = lubridate::month(x),
        day = lubridate::day(x),
        time_format=0,
        hour_minutes = "0100",
        sky_conditions = 1
      )
    }
)) %>% mutate(state = ifelse(li < 0.0005, NA,"Full" ))

runs<-rle(as.vector(moon_ill$state))
myruns = which(runs$values == "Full")
runs.lengths.cumsum = cumsum(runs$lengths)
ends = runs.lengths.cumsum[myruns]
newindex = ifelse(myruns>1, myruns-1, 0)
starts = runs.lengths.cumsum[newindex] + 1
if (0 %in% newindex) starts = c(1,starts)
starts_moon_dates <- days365[starts]
starts_end_dates  <- days365[ends]

dateRanges_moons <- data.frame(from = starts_moon_dates,
                         to = starts_end_dates) %>% na.omit()


Sys.setlocale("LC_ALL", 'en_US.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")


# full year with original and corrected data and FULL MOONS on vetrical bars
ggplot2::ggplot(data = data_long) +
  scale_fill_manual('',
                    values = 'grey',
                    guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_rect(
      data = dateRanges_moons,
      aes(
          xmin = from,
          xmax = to,
          ymin = -Inf,
          ymax = Inf#,
          #fill = bar_legend #uncomment to show in legend
      ),
      alpha = 0.5
  ) +
  #ggtitle(sprintf("%s, 2018, Mean VIIRS DNB Radiance",name)) +
  geom_line( data = data_long,aes(x = as.Date(mydates , "X%Y.%m.%d"), y =value,linetype = dataset), size=0.3) +
  scale_linetype_manual(name = NULL,#"DNB Datasets:",
                      values = c("original" = "solid", "median_shift" = "dashed","roman"= "longdash","cao"= "dotted", "DO"= "twodash"),
                      labels = c("original" = "Original", "median_shift" = "Median Shift","roman"= "BRDF","cao"= "Reflectance", "DO"= "Dark Oject")) +
  xlab("Date") +
  ylab("DNB") +
  scale_y_continuous(limit = c(0, max(data_long$value))) +theme_gray(base_size = 10)+
  theme(axis.text=element_text(size=12),
      axis.title=element_text(size=12),
      legend.position=c(0.7,0.8))+
 coord_cartesian(xlim = (c(
    as.Date("2018-06-25"), as.Date("2018-07-31")
  )))

ggsave(
    sprintf('all_datasets_%s_%s.tif',name,cfg$YEAR),
    plot = last_plot(),
    device = "tiff",
    path = here::here("output", name),
    scale = 1,
    width = 25,
    height = 10,
    units = c("cm"),
    dpi = 1200,
    limitsize = TRUE
)
