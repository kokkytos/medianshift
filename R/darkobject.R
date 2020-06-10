library(cluster)
library(dplyr)
library(rasterVis)
library(purrr)
library(ggplot2)
library(lunar)
library(raster)
library(lubridate)
library(sf)
library(parallel)
library(furrr)
library(suncalc)
library(stringr)


source(here::here("R", 'myfunctions.R'))


####====== Read settings ==================== ####
myconfig <- config::get(file = here::here("R", "config_active.yml"),config ="default")
name<-myconfig$ACTIVE
cfg <- config::get(file = here::here("R", "config.yml"),config = name)
cfg2 <- config::get(file = here::here("R", "config.yml"),config = myconfig$ACTIVE_DATASET)
####========================================= ####

# set some settings
options(scipen=999)


# generate extent object for csv data for AOI -----------------------------
ext <-getExtentfromCsv(name)
ext_DO <-getExtentfromCsv('DarkObj') #Dark object extent


# Read/crop and save raster stacks for AOI --------------------------------

# Read
DNB_na.original  <- raster::stack(here::here("data","grd", cfg2[[cfg$DNB_geotiffs_DIR]], cfg$original_DNB_grd))


# Crop

DNB_na.DO  <- crop(DNB_na.original, ext_DO) #original data Dark object
DNB_na.original  <- crop(DNB_na.original, ext) #original data
mydates <- as.Date(names( raster::stack(DNB_na.original) ), "X%Y.%m.%d")


# find and remove 100% NA from DO
percentage_DO <-calc_percentage_NA(DNB_na.DO) 
DNB_na.DO <- subset(DNB_na.DO, which(percentage_DO<100))

# find date with min phase
DO_min_phase_date <-
  suncalc::getMoonIllumination(
    date = as.Date(names(raster::stack(DNB_na.DO)), "X%Y.%m.%d"),
    keep = c("fraction", "phase", "angle")
  ) %>% dplyr::filter(phase == min(phase)) %>% pull(date)

# calculate median of Dark Object
DNB_na.DO_sub_Median <-
  cellStats(subset(DNB_na.DO, which(
    as.Date(names(raster::stack(DNB_na.DO)), "X%Y.%m.%d") == DO_min_phase_date
  )), median, na.rm = T)

DO_median <-   cellStats(DNB_na.DO, median, na.rm = T)

days365<-seq(min(mydates), max(mydates), "days")

tibl <-  tibble(mymedian = DO_median, date = as.Date(names(DO_median), "X%Y.%m.%d")) %>%
  dplyr::right_join(tibble(date = days365)) 

s1 <- splinefun(x = seq_along(1:length(tibl$mymedian)),y=tibl$mymedian, method="monoH.FC")
f  <- s1( seq_along(1:length(tibl$mymedian)))
tibl <- tibl %>%  dplyr::mutate(spline_median=as.integer(f))  


tibl <- suncalc::getMoonIllumination(date = days365,
                                     keep = c("fraction", "phase", "angle")) %>% dplyr::select(date, phase) %>% 
  inner_join(tibl) %>% 
  mutate(coef =spline_median - DNB_na.DO_sub_Median) %>%
  mutate(coef = dplyr::if_else(coef < 0., 0., coef)) %>% filter(date %in% as.Date(names(DNB_na.original), "X%Y.%m.%d")) 

coef <- tibl %>%  dplyr::pull(coef)


# apply correction
DNB_na.original_corrected <- DNB_na.original - coef
DO_corrected_df <-
  data.frame(meanDNB_corrected=cellStats(DNB_na.original_corrected, mean, na.rm = T),
             meanDNB=cellStats(DNB_na.original, mean, na.rm = T),
             dates = as.Date(names(DNB_na.original_corrected), "X%Y.%m.%d")) %>% as_tibble() %>% tidyr::drop_na() 



# metrics
mean_diff_do_corrected  <- round(mean_diff(DO_corrected_df$meanDNB_corrected),3)
write.csv2(
  data.frame(mean_diff = mean_diff_do_corrected,
             name=name
             ),
  file = here::here(
    "output",
    name,
    sprintf("Mean_Difference_DO_%s.txt", name)
  ),
  row.names = F
)

fit_corrected <- lm(formula = meanDNB_corrected ~ poly(dates, 2), DO_corrected_df)
summary(fit_corrected)
sink(here::here("output", name, sprintf("fit_DO_%s.txt",name)))
print(summary(fit_corrected))
sink()  # returns output to the console




#export data for compare
df_export <-data.frame(mydates=days365,name=name) %>% as_tibble() %>% 
    dplyr::left_join(DO_corrected_df, by=c("mydates"="dates"))

readr::write_csv(df_export, path = here::here("output", name,sprintf("%s_dataframe_datasets_DO_year_%s.csv",name, cfg$YEAR)))


####### GGPLOT ####################### 

# convert to long format for ggplot
data_long <- tidyr::gather(DO_corrected_df, dataset, value, meanDNB_corrected:meanDNB) %>%
    dplyr::filter(!is.na(value))


days365 <-  seq(min(DO_corrected_df$dates),max(DO_corrected_df$dates),by="day")
mi <- getMoonIllumination(date = days365, 
                                 keep = c("fraction", "phase", "angle")) %>% 
  mutate(state = ifelse(phase >=0.4 & phase<=0.6, "Full", NA))

runs<-rle(as.vector(mi$state))
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
























centroid <-
  as(ext, "SpatialPolygons") %>% 
  st_as_sf() %>% 
  st_centroid() %>% 
  st_set_crs(st_crs(crs(DNB_na.original_corrected) )) %>% 
  st_transform("+init=epsg:4326") %>% 
  st_coordinates() 
colnames(centroid) <- c("lon", "lat")

days365 <-  seq(min(DO_corrected_df$dates),max(DO_corrected_df$dates),by="day")
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



g <- ggplot2::ggplot(data = data_long) +
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
    ggtitle(sprintf("%s, 2018, Mean VIIRS DNB Radiance, Dark Object Method",name)) +
    geom_line( data=data_long, aes(x=dates, y=value, linetype=dataset), size=0.3) +
    scale_linetype_manual(name = NULL,#"DNB Datasets:",
                        values = c("meanDNB_corrected" = "solid", "meanDNB" = "dashed"),
                        labels = c("meanDNB_corrected" ="Corrected", "meanDNB" = "Original")) +
    xlab("Date") +
    ylab("DNB") +
    #scale_y_continuous(limit = c(0, max(data_long$value))) +theme_gray(base_size = 10)+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position=c(0.2,0.8))



ggsave(
      sprintf('fig_ggplot_DO_corrected_%s_%s.tif',name,cfg$YEAR),
      plot = g,
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 25,
      height = 10,
      units = c("cm"),
      dpi = 1200,
      limitsize = TRUE
)
