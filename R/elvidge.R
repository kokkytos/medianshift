library(stringr)
library(sf)
library(raster)
library(lubridate)
library(dplyr)
library(suncalc)
library(ggplot2)
library(stringr)

####====== Read settings ==================== ####
myconfig <- config::get(file = here::here("R", "config_active.yml"),config ="default")
name<-myconfig$ACTIVE
cfg <- config::get(file = here::here("R", "config.yml"),config = name)
####========================================= ####


source(here::here("R", 'myfunctions.R'))

# set some settings
options(scipen=999)


# Set directories variables -----------------------------------------------
time_geotiffs_DIR <- here::here("data","geotiffs_elvidge", "BlackMarble_VNP46A1_BIG_BRDF_vnp43ma4v001_h20v05")
geotiffs_DIR <- here::here("data","geotiffs_elvidge", "BlackMarble_VNP46A1_BIG_BRDF_vnp43ma4v001_h20v05")


# generate extent object for csv data for AOI -----------------------------
ext <-getExtentfromCsv(name)


# Read/crop and save raster stacks for AOI --------------------------------
time_files <- list.files(time_geotiffs_DIR , pattern = sprintf("^VNP46A1.A%s.*\\.DNB*.*time_2100.tif$",cfg$YEAR), full.names = T) #lunar corrected files (roman method)
original_dnb_files <- list.files(geotiffs_DIR, pattern = sprintf("^VNP46A1.A%s.*original*.*tif$",cfg$YEAR), full.names = T) #original dnb files.no lunar correction



#Extract Dates from filenames 
mydates_time <- extract_original_dates(time_files)
mydates_time_secondary = extract_dates(original_dnb_files)

# Crop
time_stack  <- raster::stack(time_files) %>% 
  crop(ext) %>%  
  setNames(mydates_time) 

# find na percentage
percentage <-calc_percentage_NA(time_stack) 


# Calculate moon illuminance (lux, janiczek)
lunar_illumination <- batch_moon_illuminance(time_stack)

# Plot  moon illuminance
lunar_illumination_sub <-lunar_illumination[!is.na(lunar_illumination)]
li_bck<-lunar_illumination_sub
li_bck[li_bck>0.0005]<-NA
plot(lunar_illumination_sub, type='l')
points(x= li_bck, type = "p" ,col='red', cex=0.1)



# DNB
mydates_dnb <- extract_original_dates(original_dnb_files)

original_dnb_stack <- raster::stack(original_dnb_files) %>%
  crop( ext)

# remove dnb files according to lunar illuminance
original_dnb_stack_elvidge <- original_dnb_stack %>% 
  subset(which(lunar_illumination<0.0005)) %>% setNames(mydates_time_secondary[which(lunar_illumination<0.0005)])

# σε περίπτωση διπλών ημερομηνιών κάνε mean με stackApply. Όρισε ημερομηνίες με κριτηριο την ώρα
original_dnb_stack_elvidge <-
  stackApply(original_dnb_stack_elvidge,
             as.character(mydates_time_secondary)[which(lunar_illumination < 0.0005)],
             mean,
             na.rm = T) %>% 
  setNames(unique(as.character(mydates_time_secondary)[which(lunar_illumination < 0.0005)]))


#  κράτα μόνο τα layers που έχουν πλήθος > από non NA ---------------------
countPixelsNonNa <- purrr::map(as.list(original_dnb_stack_elvidge),
                               ~ {
                                 raster::cellStats(.x, function(x, ...)
                                   sum(!is.na(x))) %>% setNames(names(.x))
                               }) %>%
  unlist()
nonNA_theshold_count  <- 30 # όριο πάνω από το οπόιο κρατάμε layers
original_dnb_stack_elvidge       <-  original_dnb_stack_elvidge[[which(countPixelsNonNa >= nonNA_theshold_count)]]


original_dnb_stack_elvidge_df <-
  data.frame(meanDNB=cellStats(original_dnb_stack_elvidge, mean, na.rm = T),
             dates = as.Date(names(original_dnb_stack_elvidge), "X%Y.%m.%d")) %>% as_tibble() %>% tidyr::drop_na() 

#mydates <- mydates[lubridate::year(mydates) == cfg$YEAR] 


mean_diff_original  <- round(mean_diff(original_dnb_stack_elvidge_df$meanDNB),3)
write.csv2(
  data.frame(mean_diff = mean_diff_original,
             name=name
             ),
  file = here::here(
    "output",
    name,
    sprintf("Mean_Difference_elvidge_%s.txt", name)
  ),
  row.names = F
)

fit_corrected <- lm(formula = meanDNB ~ poly(dates, 2), original_dnb_stack_elvidge_df)
summary(fit_corrected)
sink(here::here("output", name, sprintf("fit_elvidge_%s.txt",name)))
print(summary(fit_corrected))
sink()  # returns output to the console



























    ####### GGPLOT      ################
####### 

mytibble <-
  data.frame(
    mydates_time = mydates_time,
    mydates_dnb = mydates_dnb,
    mydates_time_secondary = mydates_time_secondary,
    lunar_illumination = lunar_illumination
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  dplyr::left_join(original_dnb_stack_elvidge_df, by = c('mydates_time_secondary' = "dates"))



days365 <-  seq(min(mytibble$mydates_time_secondary),max(mytibble$mydates_time_secondary),by="day")
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



Sys.setlocale("LC_ALL", 'en_US.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")
ggplot2::ggplot(data = mytibble) +
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
    ggtitle(sprintf("%s, 2018, Mean VIIRS DNB Radiance",name)) +
    geom_line( data=mytibble, aes(x=mydates_time_secondary, y=meanDNB, group=1)) +
  
    xlab("Date") +
    ylab("DNB") +
    #scale_y_continuous(limit = c(0, max(data_long$value))) +theme_gray(base_size = 10)+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position=c(0.2,0.8))


















t <- time_stack[[13]]
#t[is.na(t)] <- 1
#polygon <- sf::st_centroid(rasterToPolygons(t,dissolve=T) %>% sf::st_as_sf())

NonNAindex <- which(!is.na(t[]))
firstNonNA <- min(NonNAindex)
p1 <-st_transform(xyFromCell(t, firstNonNA, spatial=T) %>% sf::st_as_sf() ,"+init=epsg:4326")

lat <- st_coordinates(p1)[,"Y"]
lon <- st_coordinates(p1)[,"X"]
year <- lubridate::year((as.Date(names( t ), "X%Y.%m.%d")))
month <- lubridate::month((as.Date(names( t ), "X%Y.%m.%d")))
day <- lubridate::day((as.Date(names( t ), "X%Y.%m.%d")))
time_format <-0
sky_conditions <- 1
mytime <-t[][firstNonNA]
hour_minutes <- paste(str_pad(floor(mytime), 2, pad = "0"), str_pad(floor(round((mytime-floor(mytime))*60)), 2, pad = "0"), sep="")
# unknown <- -3
hour_minutes<-"0300"

(m.il <-janiczek(lon,lat, year, month, day, time_format, hour_minutes,sky_conditions))


# example to call janiczek from R
# result <- system("printf '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' -77 39 1986 9 13 1 1955 1 -3 | bwbasic ./R/janiczek.bas", intern=T)

#janiczek  example parameters
lon <- 23.727539
lat <- 37.983810
year <- 1986
month <- 9
day <- 13
time_format <-0
sky_conditions <- 1
hour_minutes <- 1955
unknown <- -3
(m.il <-janiczek(lon,lat, year, month, day, time_format, hour_minutes,sky_conditions))


