
library(here)
library(mapview)
library(tidyterra)
library(tidyverse)
library(rgdal)
library(KrigR)
library(terra)
library(magrittr)
library(rts)
library(spatialEco)
library(Kendall)
library(ggpubr)

# mne = montenegro

date_from <- "1991-01-01"
date_to <- "2020-01-01"

dates <- seq(as.Date(date_from), as.Date(date_to), by = "day")
years <- substr(seq(as.Date(date_from), as.Date(date_to), by = "year"), 1, 4)
months <- unique(substr(dates, 6, 7))

mne_dem <- terra::rast(here::here("data", "dem", "mn_dem.tif"))
mne_dem <- terra::aggregate(mne_dem, 6)
mne_ext <- terra::ext(mne_dem)


pet_rasters_month <- terra::rast(here::here("data", "PET", "pet_per_month.tif"))
pet_rasters_month <- terra::crop(pet_rasters_month, mne_ext)
pet_rasters_month <- terra::resample(pet_rasters_month, mne_dem)
pet_rasters_month <- terra::mask(pet_rasters_month, mne_dem)

# plot(mne_dem)

pre_rasters_month <- terra::rast(here::here("data", "ERA5_per_month", "ERA5_per_month.tif"))
pre_rasters_month <- terra::crop(pre_rasters_month, mne_ext)
pre_rasters_month <- terra::resample(pre_rasters_month, mne_dem)
pre_rasters_month <- terra::mask(pre_rasters_month, mne_dem)

# =========================== Computing Aridity Index ==========================

AI_rasters <- pre_rasters_month/pet_rasters_month

plot(AI_rasters[[2:6]])



#======================== MannKendall functions ================================

mkfun <- function(i) { Kendall::MannKendall(i)$tau}
mkfun_pvalue <- function(i) { Kendall::MannKendall(i)$sl}

#======================== Sequences of dates ===================================
date_from <- "1991-01-01"
date_to <- "2020-12-31"
month_dates <- seq(as.Date(date_from), as.Date(date_to), by = "month")
years_dates <- seq(as.Date(date_from), as.Date(date_to), by = "years")

#===============================================================================



# ============ add time to AI rasters ==========================================
terra::time(AI_rasters, tstep="months") <- rep(1:12, 30)


# ============= Computing mk statistics ========================================
mk_months_rasters <- terra::tapp(AI_rasters, "months", mkfun)
mk_months_rasters_p <- terra::tapp(AI_rasters, "months", mkfun_pvalue)

plot(mk_months_rasters_p[[12]])


names(mk_months_rasters) <- month.name
names(mk_months_rasters_p) <- month.name

# ======================== Crop to mne mask ====================================

mk_months_rasters <- terra::crop(mk_months_rasters, mne_ext)
mk_months_rasters <- terra::mask(mk_months_rasters, mne_dem)
plot(mk_months_rasters[[1]])

mk_months_rasters_p <- terra::crop(mk_months_rasters_p, mne_ext)
mk_months_rasters_p <- terra::mask(mk_months_rasters_p, mne_dem)
plot(mk_months_rasters_p[[1]])

#===============================================================================


#====================== Mean monthly AI ========================================
terra::time(AI_rasters, tstep="months") <- rep(1:12, 30)
AI_months_mean_rasters <- terra::tapp(AI_rasters, "months", mean)
names(AI_months_mean_rasters) <- month.name
plot(AI_months_mean_rasters[[1]]) # Plot mean January AI 

#===============================================================================

#============ CALCULATION FOR SEASONS ==========================================

#=================== WINTER ===================================================
terra::time(pre_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")
terra::time(pet_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")

dte <- time(pet_rasters_month)
m <- as.numeric(format(dte, "%m"))
winter_months <- (m < 3 | m > 11)

winter_pre_rasters_month <- pre_rasters_month[[winter_months]]
winter_pet_rasters_month <- pet_rasters_month[[winter_months]]

d <- terra::time(winter_pre_rasters_month)
y <- as.numeric(format(d, "%Y"))
m <- as.numeric(format(d, "%m"))

# Sum winter than yearly mean
winter_pre_rasters_month <- terra::tapp(winter_pre_rasters_month, y, sum)
winter_pet_rasters_month <- terra::tapp(winter_pet_rasters_month, y, sum)

AI_winter <- winter_pre_rasters_month/winter_pet_rasters_month
AI_winter_mean <- terra::app(AI_winter, mean)

winter_mktest_rasters <- terra::app(AI_winter, fun=mkfun, cores = 1)
winter_mktest_rasters_p <- terra::app(AI_winter, fun=mkfun_pvalue, cores = 1)
# plot(winter_mktest_rasters_p)

winter_mktest_rasters <- terra::crop(winter_mktest_rasters, mne_ext)
winter_mktest_rasters <- terra::mask(winter_mktest_rasters, mne_dem)
plot(winter_mktest_rasters[[1]])

winter_mktest_rasters_p <- terra::crop(winter_mktest_rasters_p, mne_ext)
winter_mktest_rasters_p <- terra::mask(winter_mktest_rasters_p, mne_dem)
plot(winter_mktest_rasters_p[[1]])

#=================== SPRING =========================================
terra::time(pre_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")
terra::time(pet_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")

dte <- time(pet_rasters_month)
m <- as.numeric(format(dte, "%m"))
spring_months <- (m %in% c(3, 4, 5))

spring_pre_rasters_month <- pre_rasters_month[[spring_months]]
spring_pet_rasters_month <- pet_rasters_month[[spring_months]]

d <- terra::time(spring_pre_rasters_month)
y <- as.numeric(format(d, "%Y"))
m <- as.numeric(format(d, "%m"))

# Sum spring than yearly mean
spring_pre_rasters_month <- terra::tapp(spring_pre_rasters_month, y, sum)
spring_pet_rasters_month <- terra::tapp(spring_pet_rasters_month, y, sum)

AI_spring <- spring_pre_rasters_month/spring_pet_rasters_month
AI_spring_mean <- terra::app(AI_spring, mean)

spring_mktest_rasters <- terra::app(AI_spring, fun=mkfun, cores = 1)
spring_mktest_rasters_p <- terra::app(AI_spring, fun=mkfun_pvalue, cores = 1)

spring_mktest_rasters <- terra::crop(spring_mktest_rasters, mne_ext)
spring_mktest_rasters <- terra::mask(spring_mktest_rasters, mne_dem)
plot(spring_mktest_rasters[[1]])

spring_mktest_rasters_p <- terra::crop(spring_mktest_rasters_p, mne_ext)
spring_mktest_rasters_p <- terra::mask(spring_mktest_rasters_p, mne_dem)
plot(spring_mktest_rasters_p[[1]])

#=================== SUMMER =========================================
terra::time(pre_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")
terra::time(pet_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")

dte <- time(pet_rasters_month)
m <- as.numeric(format(dte, "%m"))
summer_months <- (m %in% c(6, 7, 8))
summer_pre_rasters_month <- pre_rasters_month[[summer_months]]
summer_pet_rasters_month <- pet_rasters_month[[summer_months]]

d <- terra::time(summer_pre_rasters_month)
y <- as.numeric(format(d, "%Y"))
m <- as.numeric(format(d, "%m"))

# Sum summer than yearly mean
summer_pre_rasters_month <- terra::tapp(summer_pre_rasters_month, y, sum)
summer_pet_rasters_month <- terra::tapp(summer_pet_rasters_month, y, sum)

AI_summer <- summer_pre_rasters_month/summer_pet_rasters_month
AI_summer_mean <- terra::app(AI_summer, mean)

summer_mktest_rasters <- terra::app(AI_summer, fun=mkfun, cores = 1)
summer_mktest_rasters_p <- terra::app(AI_summer, fun=mkfun_pvalue, cores = 1)

summer_mktest_rasters <- terra::crop(summer_mktest_rasters, mne_ext)
summer_mktest_rasters <- terra::mask(summer_mktest_rasters, mne_dem)
plot(summer_mktest_rasters[[1]])

summer_mktest_rasters_p <- terra::crop(summer_mktest_rasters_p, mne_ext)
summer_mktest_rasters_p <- terra::mask(summer_mktest_rasters_p, mne_dem)
plot(summer_mktest_rasters_p[[1]])

#=================== AUTUMN =========================================
terra::time(pre_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")
terra::time(pet_rasters_month) <- seq(as.Date("1991-1-1"), as.Date("2020-12-31"), by="month")

dte <- time(pet_rasters_month)
m <- as.numeric(format(dte, "%m"))
autumn_months <- (m %in% c(9, 10, 11))
autumn_pre_rasters_month <- pre_rasters_month[[autumn_months]]
autumn_pet_rasters_month <- pet_rasters_month[[autumn_months]]

d <- terra::time(autumn_pre_rasters_month)
y <- as.numeric(format(d, "%Y"))
m <- as.numeric(format(d, "%m"))

# Sum autumn than yearly mean
autumn_pre_rasters_month <- terra::tapp(autumn_pre_rasters_month, y, sum)
autumn_pet_rasters_month <- terra::tapp(autumn_pet_rasters_month, y, sum)

AI_autumn <- autumn_pre_rasters_month/autumn_pet_rasters_month
AI_autumn_mean <- terra::app(AI_autumn, mean)

autumn_mktest_rasters <- terra::app(AI_autumn, fun=mkfun, cores = 1)
autumn_mktest_rasters_p <- terra::app(AI_autumn, fun=mkfun_pvalue, cores = 1)

autumn_mktest_rasters <- terra::crop(autumn_mktest_rasters, mne_ext)
autumn_mktest_rasters <- terra::mask(autumn_mktest_rasters, mne_dem)
plot(autumn_mktest_rasters[[1]])

autumn_mktest_rasters_p <- terra::crop(autumn_mktest_rasters_p, mne_ext)
autumn_mktest_rasters_p <- terra::mask(autumn_mktest_rasters_p, mne_dem)
plot(autumn_mktest_rasters_p[[1]])

#============================== Seasons rasters ================================

AI_seasons_rasters <- c(AI_winter_mean, AI_spring_mean, AI_summer_mean, AI_autumn_mean)
names(AI_seasons_rasters) <- c("Winter", "Spring", "Summer", "Autumn")

seasons_mktest_rasters <- c(winter_mktest_rasters, spring_mktest_rasters, summer_mktest_rasters, autumn_mktest_rasters)
names(seasons_mktest_rasters) <- c("Winter", "Spring", "Summer", "Autumn")

seasons_mktest_rasters_p <- c(winter_mktest_rasters_p, spring_mktest_rasters_p, summer_mktest_rasters_p, autumn_mktest_rasters_p)
names(seasons_mktest_rasters_p) <- c("Winter", "Spring", "Summer", "Autumn")

#============================== Yearly =========================================
pre_rasters_year <- terra::tapp(pre_rasters_month, "years", sum)
pet_rasters_year <- terra::tapp(pet_rasters_month, "years", sum)

AI_year <- pre_rasters_year/pet_rasters_year
AI_year_mean <- terra::app(AI_year, mean)

year_mktest_rasters <- terra::app(AI_year, fun=mkfun, cores = 1)
year_mktest_rasters_p <- terra::app(AI_year, fun=mkfun_pvalue, cores = 1)

year_mktest_rasters <- terra::crop(year_mktest_rasters, mne_ext)
year_mktest_rasters <- terra::mask(year_mktest_rasters, mne_dem)
plot(year_mktest_rasters[[1]])

year_mktest_rasters_p <- terra::crop(year_mktest_rasters_p, mne_ext)
year_mktest_rasters_p <- terra::mask(year_mktest_rasters_p, mne_dem)
plot(year_mktest_rasters_p[[1]])

#===============================================================================
#============================== PREPARING GGPLOTS ==============================
#===============================================================================

#========================= Monthly AI plot =====================================

AI_months_mean_plots <- as.list(rep(NA, times = 12))

for(i in 1:12){
  AI_months_mean_plots[[i]] <- ggplot() +
    geom_spatraster(data = AI_months_mean_rasters[[i]]) +
    scale_fill_whitebox_c(palette = "gn_yl", direction = -1,  na.value = 'transparent') + 
    geom_spatraster_contour(data = AI_months_mean_rasters[[i]]) + 
    labs(title = month.name[[i]], fill="AI") +
    theme_bw()
}

AI_months_mean_plots[[6]]

AI_months_mean_all_plots <- ggarrange(AI_months_mean_plots[[1]], AI_months_mean_plots[[2]], AI_months_mean_plots[[3]], AI_months_mean_plots[[4]], AI_months_mean_plots[[5]], AI_months_mean_plots[[6]],
                                      AI_months_mean_plots[[7]], AI_months_mean_plots[[8]], AI_months_mean_plots[[9]], AI_months_mean_plots[[10]], AI_months_mean_plots[[11]], AI_months_mean_plots[[12]], ncol = 3, nrow = 4)

ggsave(plot = AI_months_mean_all_plots, filename = here::here("plots", "AI_months_mean_all_plots.jpg"), dpi = 600, width = 10, height = 10)


#================ Seasonal AI PLOT =============================================

AI_seasons_rasters_plots <- as.list(rep(NA, times = 4))
for(i in 1:4){
  AI_seasons_rasters_plots[[i]] <- ggplot() +
    geom_spatraster(data = AI_seasons_rasters[[i]]) +
    scale_fill_whitebox_c(palette = "gn_yl", direction = -1,  na.value = 'transparent') + 
    geom_spatraster_contour(data = AI_seasons_rasters[[i]]) + 
    labs(title = names(AI_seasons_rasters[[i]]), fill="AI") +
    theme_bw()
}

AI_seasons_rasters_plots[[4]]

AI_seasons_plots <- ggarrange(AI_seasons_rasters_plots[[1]], AI_seasons_rasters_plots[[2]], AI_seasons_rasters_plots[[3]], AI_seasons_rasters_plots[[4]], ncol = 2, nrow = 2)

ggsave(plot = AI_seasons_plots, filename = here::here("plots", "AI_seasons_plots.jpg"), dpi = 600, width = 6, height = 6)

#================ Yearly AI plot ===================================================

AI_year_rasters_plot <- ggplot() +
  geom_spatraster(data = AI_year_mean) +
  scale_fill_whitebox_c(palette = "gn_yl", direction = -1, na.value = 'transparent') + 
  geom_spatraster_contour(data = AI_year_mean) + 
  labs(fill="AI") +
  theme_bw()

AI_year_rasters_plot

ggsave(plot = AI_year_rasters_plot, filename = here::here("plots", "AI_year_rasters_plot.jpg"), dpi = 600, width = 4, height = 4)

#===============================================================================

#================ mktest plots  ====================================================

#================ mktest monthly ===================================================

mktest_months_rasters_plot <- ggplot() +
  geom_spatraster(data = mk_months_rasters) +
  facet_wrap(~lyr, ncol = 4) + 
  scale_fill_whitebox_c(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="MK tau") +
  theme_bw()

mktest_months_rasters_plot

ggsave(plot = mktest_months_rasters_plot, filename = here::here("plots",  "mktest_months_rasters_plot.jpg"), dpi = 600, width = 10, height = 10)

#=========== Significant months ================================================

sig_mktest_months_rasters_plot <- ggplot() +
  geom_spatraster(data = mk_months_rasters[[c(1, 4, 6)]]) +
  facet_wrap(~lyr, ncol = 3) + 
  scale_fill_whitebox_c(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="MK tau") +
  theme_bw()

sig_mktest_months_rasters_plot

ggsave(plot = sig_mktest_months_rasters_plot, filename = here::here("plots", "sig_mktest_months_rasters_plot.jpg"), dpi = 600, width = 8, height = 5)

#================ mktest seasonal ==============================================

mktest_seasons_rasters_plot <- ggplot() +
  geom_spatraster(data = seasons_mktest_rasters) +
  facet_wrap(~lyr, ncol = 2) + 
  scale_fill_whitebox_c(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="MK tau") +
  theme_bw()

mktest_seasons_rasters_plot

ggsave(plot = mktest_seasons_rasters_plot, filename = here::here("plots", "mktest_seasons_rasters_plot.jpg"), dpi = 600, width = 6, height = 6)

#============= Significan seasons (winter and autumn) ==========================
winter_autumn_mktest_seasons_rasters_plot <- ggplot() +
  geom_spatraster(data = seasons_mktest_rasters[[c(1, 4)]]) +
  facet_wrap(~lyr, ncol = 2) + 
  scale_fill_whitebox_c(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="MK tau") +
  theme_bw()

winter_autumn_mktest_seasons_rasters_plot

ggsave(plot = winter_autumn_mktest_seasons_rasters_plot, filename = here::here("plots", "winter_autumn_mktest_seasons_rasters_plot.jpg"), dpi = 600, width = 8, height = 5)

#================ mktest yearly plot ===================================================

mktest_year_rasters_plot <- ggplot() +
  geom_spatraster(data = year_mktest_rasters) +
  scale_fill_whitebox_c(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="MK tau") +
  theme_bw()

mktest_year_rasters_plot

ggsave(plot = mktest_year_rasters_plot, filename = here::here("plots", "mktest_year_rasters_plot.jpg"), dpi = 600, width = 4, height = 4)


#===============================================================================

#================ mktest p-value ====================================================

#================ mktest p-value monthly ===================================================

mk_significant_rasters_negative <- (mk_months_rasters < 0) * (mk_months_rasters_p < 0.1)*-1
mk_significant_rasters_positive <- (mk_months_rasters > 0) * (mk_months_rasters_p < 0.1)

mk_significant_rasters <- mk_significant_rasters_negative + mk_significant_rasters_positive

# TRUE + TRUE = 2 - Negative
# TRUE + FALSE = 1 - None
# FALSE + TRUE = 1 - Positive
# FALSE + FALSE = 0 - None

#mk_significant_rasters_none <- mk_significant_rasters == 3


cls <- data.frame(id=-1:1, Significance = c("Negative", "None", "Positive"))


for(i in 1:12){
  levels(mk_significant_rasters[[i]]) <- cls
}

names(mk_significant_rasters) <- month.name

mk_significant_rasters_plot <- ggplot() +
  geom_spatraster(data = mk_significant_rasters) +
  facet_wrap(~lyr, ncol = 4) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = -1, na.value = 'transparent', labels = c("Negative", "None", "Positive")) + 
  labs(fill="Significance") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"))

mk_significant_rasters_plot

ggsave(plot = mk_significant_rasters_plot, filename = here::here("plots", "mk_significant_rasters_plot.jpg"), dpi = 600, width = 10, height = 10)



months_mk_significant_rasters_plot <- ggplot() +
  geom_spatraster(data = mk_significant_rasters[[c(1, 4, 6)]]) +
  facet_wrap(~lyr, ncol = 3) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = -1, na.value = 'transparent', labels = c("Negative", "None", "Positive")) + 
  labs(fill="Significance") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"))

months_mk_significant_rasters_plot

ggsave(plot = months_mk_significant_rasters_plot, filename = here::here("plots", "months_mk_significant_rasters_plot.jpg"), dpi = 600, width = 8, height = 5)


other_mk_significant_rasters_plot <- ggplot() +
  geom_spatraster(data = mk_significant_rasters[[-c(1, 4, 6)]]) +
  facet_wrap(~lyr, ncol = 3) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = -1, na.value = 'transparent', labels = c("Negative", "None", "Positive")) + 
  labs(fill="Significance") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"))

other_mk_significant_rasters_plot

ggsave(plot = other_mk_significant_rasters_plot, filename = here::here("plots", "other_mk_significant_rasters_plot.jpg"), dpi = 600, width = 8, height = 5)

#===============================================================================


mktest_months_rasters_p_plot <- ggplot() +
  geom_spatraster(data = mk_months_rasters_p < 0.1) +
  facet_wrap(~lyr, ncol = 4) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="p-value < 0.1") +
  theme_bw()

mktest_months_rasters_p_plot

ggsave(plot = mktest_months_rasters_p_plot, filename = here::here("plots", "mktest_months_rasters_p_plot.jpg"), dpi = 600, width = 10, height = 10)

sig_mktest_months_rasters_p_plot <- ggplot() +
  geom_spatraster(data = mk_months_rasters_p[[c(1, 4, 6)]] < 0.1) +
  facet_wrap(~lyr, ncol = 3) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="p-value < 0.1") +
  theme(strip.text = element_text(size = 60)) + 
  theme_bw()

sig_mktest_months_rasters_p_plot

ggsave(plot = sig_mktest_months_rasters_p_plot, filename = here::here("plots", "sig_mktest_months_rasters_p_plot.jpg"), dpi = 600, width = 8, height = 5)


#================ mktest p-value seasonal ===================================================

mktest_seasons_rasters_p_plot <- ggplot() +
  geom_spatraster(data = seasons_mktest_rasters_p < 0.1) +
  facet_wrap(~lyr, ncol = 2) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="p-value < 0.1") +
  theme_bw()

mktest_seasons_rasters_p_plot

ggsave(plot = mktest_seasons_rasters_p_plot, filename = here::here("plots", "mktest_seasons_rasters_p_plot.jpg"), dpi = 600, width = 6, height = 6)

#===================== Significant seasons p - plot ============================

winter_autumn_mktest_seasons_rasters_p_plot <- ggplot() +
  geom_spatraster(data = seasons_mktest_rasters_p[[c(1, 4)]] < 0.1) +
  facet_wrap(~lyr, ncol = 2) + 
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="p-value < 0.1") +
  theme_bw()

winter_autumn_mktest_seasons_rasters_p_plot

ggsave(plot = winter_autumn_mktest_seasons_rasters_p_plot, filename = here::here("plots",  "winter_autumn_mktest_seasons_rasters_p_plot.jpg"), dpi = 600, width = 8, height = 5)



#================ mktest p-value yearly ===================================================

mktest_year_rasters_p_plot <- ggplot() +
  geom_spatraster(data = year_mktest_rasters_p < 0.1) +
  scale_fill_whitebox_d(palette = "bl_yl_rd", direction = 1, na.value = 'transparent') + 
  labs(fill="p-value < 0.1") +
  theme_bw()

mktest_year_rasters_p_plot

ggsave(plot = mktest_year_rasters_p_plot, filename = here::here("plots",  "mktest_year_rasters_p_plot.jpg"), dpi = 600, width = 4, height = 4)

