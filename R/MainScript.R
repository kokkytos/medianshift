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
library(zoo)


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



# Read/crop and save raster stacks for AOI --------------------------------

# Read
DNB_na.original  <- raster::stack(here::here("data",cfg2$grd_dir, cfg2[[cfg$DNB_geotiffs_DIR]], cfg$original_DNB_grd))
DNB_na <- raster::stack(here::here("data",cfg2$grd_dir, cfg2[[cfg$DNB_geotiffs_DIR]], cfg$DNB_grd))

# Crop
DNB_na.original  <- crop(DNB_na.original, ext) #original data
DNB_na  <- crop(DNB_na, ext) # lunar corrected (roman or cao method)
mydates <- as.Date(names( raster::stack(DNB_na.original) ), "X%Y.%m.%d")


# φίλτρο για ένα έτος
filter_dates <-
  mydates >= sprintf("%s-01-01", cfg$YEAR) &
  mydates <= sprintf("%s-12-31", cfg$YEAR)

DNB_na.original   <- subset(DNB_na.original  , which(filter_dates)) #keep only 2018
DNB_na  <- subset(DNB_na  , which(filter_dates)) #keep only 2018
mydates <- mydates[filter_dates]


#  κράτα μόνο τα layers που έχουν πλήθος > από non NA ---------------------

countPixelsNonNa <- purrr::map(as.list(DNB_na.original),
                               ~ {
                                 raster::cellStats(.x, function(x, ...)
                                   sum(!is.na(x))) %>% setNames(names(.x))
                               }) %>%
  unlist()

nonNA_theshold_count  <- 30 # όριο πάνω από το οπόιο κρατάμε layers
DNB_na.original       <-  DNB_na.original[[which(countPixelsNonNa >= nonNA_theshold_count)]]
DNB_na                <-  DNB_na[[which(countPixelsNonNa >= nonNA_theshold_count)]]


### extract dates from layer names of raster stack
mydates <- as.Date(names(raster::stack(DNB_na.original)), "X%Y.%m.%d")
# .............................................................................




# Lunar Data
moon_df <- data.frame(mydates = mydates, month = lubridate::month(mydates), year=lubridate::year(mydates), phase=lunar.illumination(mydates))

newmoons <- moon_df %>%
    dplyr::group_by(year, month) %>%
    dplyr::mutate(phase = min(phase)) %>%
    dplyr::left_join( moon_df, by=(c("month"="month", "year"="year", "phase"="phase"))) %>%
    dplyr::distinct(mydates.y) %>%
    dplyr::pull(mydates.y)

fullmoons <- moon_df %>%
    dplyr::group_by(year, month) %>%
    dplyr::mutate(phase = max(phase)) %>%
    dplyr::left_join( moon_df, by=(c("month"="month", "year"="year", "phase"="phase"))) %>%
    dplyr::distinct(mydates.y) %>%
    dplyr::pull(mydates.y)


# διορθωμένα raster stacks
DNB_na.original_corrected <- correct_DNB_30Days(DNB_na.original, newmoons)

# υπολογισμός mean διορθωμένων raster raster stacks
DNB_na.original_stat <- cellStats(DNB_na.original, mean, na.rm = T)
DNB_na.original_corrected_stat <- cellStats(DNB_na.original_corrected, mean, na.rm = T)
DNB_na_stat <- cellStats(DNB_na, mean, na.rm = T)


# just a plot with all time series
plot(DNB_na.original_stat, type="l")
lines(DNB_na_stat, col="green" )
lines(DNB_na.original_corrected_stat, col="blue" )







#  Metric 1. mean difference of the average value of a scene in on --------


mean_diff_original  <- mean_diff(DNB_na.original_stat)
mean_diff_corrected <-mean_diff(DNB_na.original_corrected_stat)
mean_diff_corrected_brdf <-mean_diff(DNB_na_stat)

(v_mean_diff <-
  c(mean_diff_original,
    mean_diff_corrected_brdf,
    mean_diff_corrected
    ) %>% setNames(c(
      "Mean diff Original",
      "Mean diff Corrected", # Brdf or Reflectance corrected
      "Mean diff Median Shift Corrected"

    )))





# Metric 2. The deviation to the second order regression line, fit --------

# Δοκιμάζουμε fit στα αδιόρθωτα και διορθωμένα δεδομένα, second order regression.
# Στόχος είναι στα διορθωμένα δεδομένα να έχουμε καλυτερο R2 και μικρότερα κατάλοιπα(?)
# Ελέγχουμε για R2 (Multiple R-squared) και για Residuals.
############# second order regression #############
reg_df <- data.frame(dnb=DNB_na.original_stat, dnb_corrected=DNB_na.original_corrected_stat,dnb_corrected_roman=DNB_na_stat, time=c(1:length(mydates)))
fit_original <- lm(formula = dnb ~ poly(time, 2), reg_df)
fit_corrected <- lm(formula = dnb_corrected ~ poly(time, 2), reg_df)
fit_corrected_roman <- lm(formula = dnb_corrected_roman ~ poly(time, 2), reg_df)

par(mfrow=c(1,2))

plot(x = reg_df$time, y=reg_df$dnb, xlab = "Days", ylab = "Original DNB")
predictedcounts <- predict(fit_original,list(time=reg_df$time), time2=reg_df$time^2)
lines(reg_df$time, predictedcounts,  lwd = 3, xlab = "Time (s)", ylab = "Counts",col = "blue")


plot(x = reg_df$time, y=reg_df$dnb_corrected, xlab = "Days", ylab = "Corrected DNB")
predictedcounts <- predict(fit_corrected,list(time=reg_df$time), time2=reg_df$time^2)
lines(reg_df$time, predictedcounts,  lwd = 3, xlab = "Time (s)", ylab = "Counts",col = "blue")

plot(x = reg_df$time, y=reg_df$dnb_corrected_roman, xlab = "Days", ylab = "BRDF Corrected DNB")
predictedcounts <- predict(fit_corrected_roman,list(time=reg_df$time), time2=reg_df$time^2)
lines(reg_df$time, predictedcounts,  lwd = 3, xlab = "Time (s)", ylab = "Counts",col = "blue")

summary(fit_original)
summary(fit_corrected_roman)
summary(fit_corrected)



ylim=c(min(reg_df$dnb), max(reg_df$dnb))
g_secOrdReg_1 <- ggplot_g_second_order_regression(reg_df,"dnb", "Original DNB" )
g_secOrdReg_2 <- ggplot_g_second_order_regression(reg_df, "dnb_corrected", "Corrected DNB" )



#### Batch save ggplots
plots <- list(g_secOrdReg_1, g_secOrdReg_2)
invisible(
  lapply(
    seq_along(plots), 
    function(x) {
    ggsave(
      sprintf('fig.5_%s_ggplot_second_order_regression_%s_%s.tif',x, name, cfg$YEAR),
      plot = plots[[x]],
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 9,
      height = 5,
      units = c("cm"),
      dpi = 600,
      limitsize = TRUE
)}))




###### End of second order regression #######################





# generate dates sequence for one year
days365<-seq(min(mydates), max(mydates), "days")

# combine mean (or SoL), full dates year
my_dataframe <-data.frame(mydates=days365) %>%
    dplyr::left_join(data.frame(mydates=mydates,
                                original=DNB_na.original_stat,
                                median_shift_corrected=DNB_na.original_corrected_stat
                                ))

#export data for compare
df_export <-data.frame(mydates=days365,name=name) %>%
    dplyr::left_join(data.frame(mydates=mydates,
                                original=DNB_na.original_stat,
                                median_shift_corrected=DNB_na.original_corrected_stat,
                                roman_corrected=DNB_na_stat
                                ))

readr::write_csv(df_export, path = here::here("output", name,sprintf("%s_dataframe_datasets2_%s_%s.csv",name, myconfig$ACTIVE_DATASET, cfg$YEAR)))

# convert to long format for ggplot
data_long <- tidyr::gather(my_dataframe, dataset, value, original:median_shift_corrected) %>%
    dplyr::filter(!is.na(value))




# ggplot ------------------------------------------------------------------

# days365 <- seq(min(mydates),max(mydates),by="day")
# mi <- suncalc::getMoonIllumination(date = days365, 
#                                  keep = c("fraction", "phase", "angle")) %>% 
#   mutate(state = ifelse(phase >=0.4 & phase<=0.6, "Full", NA))
# 
# runs<-rle(as.vector(mi$state))
# myruns = which(runs$values == "Full")
# runs.lengths.cumsum = cumsum(runs$lengths)
# ends = runs.lengths.cumsum[myruns]
# newindex = ifelse(myruns>1, myruns-1, 0)
# starts = runs.lengths.cumsum[newindex] + 1
# if (0 %in% newindex) starts = c(1,starts)
# starts_moon_dates <- days365[starts]
# starts_end_dates  <- days365[ends]
# 
# dateRanges_moons <- data.frame(from = starts_moon_dates,
#                          to = starts_end_dates) %>% na.omit()


centroid <-
  as(ext, "SpatialPolygons") %>% 
  st_as_sf() %>% 
  st_centroid() %>% 
  st_set_crs(st_crs(crs(DNB_na.original_corrected) )) %>% 
  st_transform("+init=epsg:4326") %>% 
  st_coordinates() 
colnames(centroid) <- c("lon", "lat")

days365 <- seq(min(mydates),max(mydates),by="day")
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
        sky_conditions = 1,
        unknown = -3
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
(
  g_year_fullmoons  <-
    ggplot_Original_Corrected_Data(
      data_long,
      dateRanges_moons,
      "0.4<=moon phase<= 0.6,\n(0.5=Full Moon)"
    )
) 

# a PERIOD of year with original and corrected data and  FULL MOONS on vetrical bars
(g_month_fullmoons <-
    g_year_fullmoons + coord_cartesian(xlim = (c(
      as.Date(sprintf("%s-06-25",cfg$YEAR)), as.Date(sprintf("%s-07-31",cfg$YEAR))
    )))) #zoom to date range


# a PERIOD of year with original and corrected data and  FULL MOONS on vertical dotted lines


mi <- suncalc::getMoonIllumination(days365) %>% mutate(id=c(1:365))
#find local maxima
a <-rollapply(as.zoo(mi$fraction), 3, function(x) which.max(x)==2) %>% data.frame(id=index(.), FM=.)

fullmoon_dates <-mi %>% left_join(a) %>% filter(FM==T) %>% pull(date)



g_month_fullmoons2 <- ggplot2::ggplot(data = data_long) +
    # scale_fill_manual('',
    #                   values = 'grey',
    #                   guide = guide_legend(override.aes = list(alpha = 1))) +
    # geom_rect(
    #     data = dateRanges_moons,
    #     aes(
    #         xmin = from,
    #         xmax = to,
    #         ymin = -Inf,
    #         ymax = Inf#,
    #         #fill = bar_legend #uncomment to show in legend
    #     ),
    #     alpha = 0.5
    # ) +
    #ggtitle(sprintf("%s, 2018, Mean VIIRS DNB Radiance",name)) +

    geom_line( data = data_long,aes(x = as.Date(mydates , "X%Y.%m.%d"), y =value,linetype = dataset), size=0.3) +
    scale_linetype_manual(name = NULL,#"DNB Datasets:",
                        values = c("median_shift_corrected" = "solid", "original" = "dashed"),
                        labels = c("Median shift", "Original")) +
    xlab("Date") +
    ylab("DNB") +
    scale_y_continuous(limit = c(0, max(data_long$value))) +theme_gray(base_size = 7)+
    theme(axis.text=element_text(size=7),
        axis.title=element_text(size=7),
        legend.position=c(0.65,0.8))+
    coord_cartesian(xlim = (c(
      as.Date("2018-06-25"), as.Date("2018-07-31")
    )))+
    geom_vline(xintercept = fullmoon_dates+1, linetype="dotted", 
                color = "black", size=1)




### save all plots as tiff

plots <- list(g_year_fullmoons,g_month_fullmoons2)
invisible(
  lapply(
    seq_along(plots), 
    function(x) {
    ggsave(
      sprintf('fig.6_ggplot_original_corrected_%s_%s_%s_%s.tif',name,x,cfg$YEAR,myconfig$ACTIVE_DATASET),
      plot = plots[[x]],
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 9,
      height = 5,
      units = c("cm"),
      dpi = 600,
      limitsize = TRUE
) }))

ggsave(
      sprintf('fig.4_ggplot_original_corrected_%s_%s_%s.tif',name,"month",cfg$YEAR,myconfig$ACTIVE_DATASET),
      plot = g_month_fullmoons,
      device = "tiff",
      path = here::here("output", name),
      scale = 1,
      width = 20,
      height = 10,
      units = c("cm"),
      dpi = 600,
      limitsize = TRUE
)

# ..............................................................................





# histograms New/Full moons -----------------------------------------------



newmoon_date <-newmoons[1]
g1 <- ggplot_hist(DNB_na.original, newmoon_date,sprintf("New Moon (original), %s",newmoon_date))
g2 <- ggplot_hist(DNB_na.original_corrected, newmoon_date,sprintf("New Moon (median shift corrected), %s",newmoon_date))


full_date <-fullmoons[1] #fullmoons[5]
g3 <- ggplot_hist(DNB_na.original, full_date,sprintf("Full Moon (original), %s",full_date))
g4 <- ggplot_hist(DNB_na.original_corrected, full_date,sprintf("Full Moon (median shift corrected), %s",full_date))



### Save all plots as tiff
plots <- list(g1, g2, g3, g4)
invisible(
  lapply(
    seq_along(plots), 
    function(x) {
   ggsave(
    sprintf('hist_g%s_%s.tif',x, cfg$YEAR),
    plot = plots[[x]],
    device = "tiff",
    path = here::here("output", name),
    scale = 1,
    width = 9,
    height = 5,
    units = c("cm"),
    dpi = 600,
    limitsize = TRUE
)
       }))





##################  Roman correction ggplot ##########################################
# combine mean (or SoL), full dates year
my_dataframe <-data.frame(mydates=days365) %>%
    dplyr::left_join(data.frame(mydates=mydates,
                                original=DNB_na.original_stat,
                                roman_corrected=DNB_na_stat
                                ))


# convert to long format for ggplot
data_long <- tidyr::gather(my_dataframe, dataset, value, original:roman_corrected) %>%
    dplyr::filter(!is.na(value))


# full year with original and corrected data and FULL MOONS on vertical bars
g <-ggplot2::ggplot(data = data_long) +
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
        ggtitle(sprintf("%s, %s, Mean VIIRS DNB Radiance (%s)",name, cfg$YEAR, myconfig$ACTIVE_DATASET)) +
    geom_line( data = data_long,aes(x = as.Date(mydates , "X%Y.%m.%d"), y =value,linetype = dataset), size=0.3) +
    scale_linetype_manual(name = NULL,#"DNB Datasets:",
                        values = c("roman_corrected" = "solid", "original" = "dashed"),
                        labels = c("Original", "Corrected")) +
    xlab("Date") +
    ylab("DNB") +
    scale_y_continuous(limit = c(0, max(data_long$value))) +theme_gray(base_size = 10)+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position=c(0.2,0.8))


ggsave(
      sprintf('fig_ggplot_%s_corrected_%s_%s.tif',name,myconfig$ACTIVE_DATASET, cfg$YEAR),
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
