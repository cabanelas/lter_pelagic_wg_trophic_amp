#############################################################
#############          Pelagic Synthesis WG        ##########
#############                NES                   ##########
#############         Trophic Amplification        ##########
#############################################################
# updated MAR-2025
# by: Alexandra Cabanelas
# R version 4.3.3
## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##

library(tidyverse) #v2.0.0
library(runner) #for running mean calc; v0.4.3
#library(zoo) has other runmean functions

# ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
# can download EcoMon data from 
#https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513 
# NCEI Accession 0187513 v3.3 

abu <- read.csv(file.path("raw",
                          "EcoMon_v3_8_wDateStrata.csv"))

## ------------------------------------------ ##
#            Tidy Data -----
## ------------------------------------------ ##

# add season column  
abu <- abu %>%
  mutate(season = case_when(between(month, 3, 5) ~ "Spring",
                            between(month, 6, 8) ~ "Summer",
                            between(month, 9, 11) ~ "Fall",
                            TRUE ~ "Winter"))

# remove the "outside" region 
abu <- abu %>% 
  filter(region != 0)

# assign region names
abu <- abu %>%
  mutate(region_name = case_when(region == 1 ~ "MAB",
                                 region == 2 ~ "MAB", #SNE
                            #naming SNE as MAB to align w phyto ts
                                 region == 3 ~ "GB",
                                 region == 4 ~ "GOM",
                                 TRUE ~ "Outside"))

## ------------------------------------------ ##
#            FFish Data -----
## ------------------------------------------ ##
## Larval Forage Fish

# atlantic herring == Clupea harengus
# atlantic mackerel == scomber scombrus
# atlantic butterfish == peprilus spp
# sand lance == ammodytes spp

# alewife == not in data = river h = alosa
# blueback herring == not in data = river h = alosa

# select columns  
ffish <- abu %>%
  select(
    cruise_name, station, zoo_gear, ich_gear, lat, lon, region, region_name, 
    strata_47, strata_46, strata_26, date, month, day, year, doy, time, depth, 
    season, sfc_temp, sfc_salt, btm_temp, btm_salt, 
    cluhar_10m2, # atlantic herring == Clupea harengus
    scosco_10m2, # atlantic mackerel == scomber scombrus
    pepspp_10m2, # atlantic butterfish == peprilus spp
    ammspp_10m2 # sand lance == ammodytes spp
    #cluhar_100m3, scosco_100m3, pepspp_100m3, ammspp_100m3
  )
# nofish_abnd; bretyr_abnd::lopame_abnd

# when the ich_gear cell is empty it seems like they did not sample for larval f
abu %>%
  summarise(across(bretyr_10m2:lopame_10m2, ~ sum(is.na(.))))
sum(abu$ich_gear == "")
# 3425 rows that dont have icht data

ffish <- ffish %>% 
  filter(ich_gear != "")

# there are 2 entries that are the same for:
ffish %>%
  count(cruise_name, station, region_name, date) %>%
  filter(n > 1)
#cruise == PC2104   station ==  122 
#cruise == PC2106   station ==   88
# the fish values are the same for both entries; the zp are diff 
# seems like 2 different zooplankton samples (different volumes) but same fish vals
ffish <- ffish %>%
  arrange(cruise_name, station, region_name, date, 
          desc(zoo_gear == "6B3Z")) %>% 
  distinct(cruise_name, station, region_name, date, .keep_all = TRUE)
#28995

## ------------------------------------------ ##
#            Total FFish Abundance -----
## ------------------------------------------ ##
# USING 10m2

# sum abundances of diff fishes to get a total of forage fish 
ffish <- ffish %>%
  group_by(cruise_name, station, region_name, date) %>% #grouping by sample
  mutate(fishsum = sum(cluhar_10m2, 
                       scosco_10m2, 
                       pepspp_10m2, 
                       ammspp_10m2, na.rm = T))

####
ggplot(ffish, aes(x = year, y = fishsum)) +
  geom_point() +
  facet_wrap(~region_name, scales = "free")

ggplot(ffish, aes(x = year, y = fishsum)) +
  geom_point() +
  facet_grid(season~region_name, scales = "free")
 

#############################################################
######   ----     Trophic Amplification:
# 1) log10(x+ (min/2))    for each station ; x = sum trophic level 
# 2) average across stations for a cruise or year/season
# 3) take a running mean with timespan ~ longest lived taxon 
# 4) compute st.dev. of time-series 
#############################################################

################################################
## ------------------------------------------ ##
#           1) Data Transformation -----
## ------------------------------------------ ##
# find the minimum non-zero value for each region and season

ffish <- ffish %>%
  group_by(region_name, season) %>%
  mutate(min_nonzero = min(fishsum[fishsum > 0], 
                           na.rm = TRUE)) %>%
  ungroup()

# Log10 transformation
ffish$log10_FF <- log10(ffish$fishsum + ffish$min_nonzero/2)

ggplot(ffish, aes(x = year, y = log10_FF)) +
  geom_point() +
  facet_grid(season~region_name, scales = "free")


## 1978-1988
## split into two time series intervals to match phyto ts available
ts1 <- ffish %>%
  filter(between(year, 1978, 1987))

ts2 <- ffish %>%
  filter(between(year, 1998, max(year)))


## ------------------------------------------ ##
#           2) Average -----
## ------------------------------------------ ##
# average across stations for a cruise or year/season

## ----
# Entire time series
## ----
ffish <- ffish %>%
  group_by(region_name, year, season) %>% #by cruise too?no
  mutate(log_mean_FF = mean(log10_FF, na.rm = T)) %>%
  ungroup()

ffish_a <- ffish %>%
  group_by(region_name, year, season) %>%
  summarize(log_mean_FF = mean(log10_FF, na.rm = TRUE), .groups = "drop")

## ----
# Time series 1 = 1978-1987
## ----
ts1 <- ts1 %>%
  group_by(region_name, year, season) %>% 
  mutate(log_mean_FF = mean(log10_FF, na.rm = T)) %>%
  ungroup()

## ----
# Time series 2 = 1998-now
## ----
ts2 <- ts2 %>%
  group_by(region_name, year, season) %>% 
  mutate(log_mean_FF = mean(log10_FF, na.rm = T)) %>%
  ungroup()

ts2_a <- ts2 %>%
  group_by(region_name, year, season) %>%
  summarize(log_mean_FF = mean(log10_FF, na.rm = TRUE), .groups = "drop")

##
ggplot(ts1, aes(x = year, y = log_mean_FF)) +
  geom_point() +
  facet_grid(season~region_name) 


## ------------------------------------------ ##
#    3) Running mean = 5 years  -----
## ------------------------------------------ ##
# take a running mean with timespan ~ longest lived taxon
# **im not entirely sure about this; have to read on runner mean run

## ----
# Entire time series 
## ----
rm_fish_tsFull <- ffish %>%
  distinct(year, region_name, season, .keep_all = T) %>% #i need this since run mean does it by order
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  mutate(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k=5))
# mutate(running_mean = rollmean(aver, k = 5, fill = NA, align = "right"))

# same as above but less columns 
rm_fish_tsFull_a <- ffish_a %>%
  #distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  #summarize(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k=5))
  mutate(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k = 5))

## ----
# Time series 1 = 1978-1987 
## ----
rm_fish_ts1 <- ts1 %>%
  distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  mutate(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k=5))

## ----
# Time series 2 = 1998-now
## ----
rm_fish_ts2 <- ts2 %>%
  distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>%
  arrange(year) %>%
  mutate(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k=5))

rm_fish_ts2_a <- ts2_a %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  mutate(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k = 5))


## ------------------------------------------ ##
#    4) Compute SD   -----
## ------------------------------------------ ##
# compute st.dev. of time-series 

## ----
# Entire time series
## ----
ffish <- ffish %>%
  #distinct(region_name, season, .keep_all = TRUE) %>%
  group_by(region_name, season) %>% 
  mutate(sd_FF = sd(log_mean_FF, na.rm = T)) 
#write.csv(ffish, "output/trophamp_FFISH_NES.csv")

#SD only
sd_foragefish <- ffish %>%
  select(region_name, season, sd_FF) %>%
  distinct(region_name, season, .keep_all = T) %>%
  arrange(desc(region_name))
#write.csv(sd_foragefish, "output/SD_trophamp_FFISH_NES.csv") 

distinct_tsFull_ff <- ffish %>%
  distinct(region_name, season, year, .keep_all = T)

rm_fish_tsFull_FINAL <- rm_fish_tsFull_a %>%
  left_join(distinct_tsFull_ff %>% 
              select(-log_mean_FF), 
              by = c("region_name", "season", "year")) %>%
  select(all_of(setdiff(names(.), names(rm_fish_tsFull_a))), 
         everything()) %>%
  #reorder col
  relocate(region_name, .after = region) %>%
  relocate(year, .after = day) %>%
  relocate(season, .after = depth) %>%
  relocate(log_mean_FF, .after = log10_FF)
#write.csv(rm_fish_tsFull_FINAL, "output/trophamp_FFISH_runmean_NES.csv")

## ----
# Time series 1 = 1978-1987
## ----
ts1_sd <- ts1 %>%
  #distinct(region_name, season, .keep_all = TRUE) %>%
  group_by(region_name, season) %>%
  mutate(sd_FF = sd(log_mean_FF, na.rm = T)) 

## ----
# Time series 2 = 1998-now
## ----
ts2 <- ts2 %>%
  group_by(region_name, season) %>%
  mutate(sd_FF = sd(log_mean_FF, na.rm = T)) 
#write.csv(ts2, "output/trophamp_FFISH_NES_1998_2021.csv")

distinct_ts2_ff <- ts2 %>%
  distinct(region_name, season, year, .keep_all = T)

rm_fish_ts2_FINAL <- rm_fish_ts2_a %>%
  left_join(distinct_ts2_ff %>% 
            select(-log_mean_FF), 
            by = c("region_name", "season", "year")) %>%
  select(all_of(setdiff(names(.), names(rm_fish_ts2_a))), 
         everything()) %>%
  #reorder col
  relocate(region_name, .after = region) %>%
  relocate(year, .after = day) %>%
  relocate(season, .after = depth) %>%
  relocate(log_mean_FF, .after = log10_FF)
#write.csv(rm_fish_ts2_FINAL, "output/trophamp_FFISH_runmean_NES_1998_2021.csv")

## ------------------------------------------ ##
#    Plots   -----
## ------------------------------------------ ##
ggplot(ts1, aes(x = year, y=sd_FF)) +
  geom_point(color = "red") +
  facet_grid(season~region_name)

ggplot() +
  geom_point(data = ts1, aes(x = year, y = log_mean_FF)) +
  geom_point(data = ts1_sd, aes(x = year, y = sd_FF,color= "red")) +
  facet_grid(season~region_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----
# plot - FULL TS
## ----

ggplot(ffish, aes(x = year, y = log10_FF)) +
  geom_point() +
  facet_grid(season~region_name, scales = "free")

ggplot() +
  geom_point(data = ffish, aes(x = year, y = log10_FF)) +
  geom_point(data = ffish, aes(x = year, y = log_mean_FF),
             color = "red", shape = 1, size = 3) +
  geom_line(data = rm_fish_tsFull, aes(x = year, y = runMean_FF), 
            color = "red", 
            size = 2) +
  ggtitle("FISH") + 
  facet_grid(season~region_name, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())


ggplot() +
  #geom_point(data = subF2, aes(x = year, y = log10Sum)) +
  geom_line(data = ffish, aes(x = year, y = log_mean_FF),
             color = "black") +
  geom_line(data = rm_fish_tsFull, aes(x = year, y = runMean_FF), 
            color = "red", 
            size = 2) +
  ggtitle("LARVAL FISH") + 
  facet_grid(season~region_name, 
             #scales = "free"
             ) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, 
                                    #colour = "orange", 
                                    #angle = 90
                                    ))


rm(list = setdiff(ls(), "abu"))

#############################################################
#############################################################
## ------------------------------------------ ##
#      Zooplankton DISPLACEMENT VOLUME -----
## ------------------------------------------ ##
#############################################################

# subset columns for zooplankton
ZP <- abu %>%
  select(
    cruise_name, station, zoo_gear, ich_gear, lat, lon, region, region_name, 
    strata_47, strata_46, strata_26, date, month, day, year, doy, time, depth, 
    season, sfc_temp, sfc_salt, btm_temp, btm_salt, volume_1m2)
# volume_1m2 = Zooplankton Displacement Volume (ml) per 1m2 of surface area

sum(is.na(abu$volume_1m2))
sum(is.na(abu$volume_100m3))
#3538 NAs
sum(abu$zoo_gear == "")
#3441

abu %>%
  summarise(across(ctyp_10m2:pnepau_10m2, ~ sum(is.na(.))))
# for some reason there are observations that have zp abundance but the 
# volume is still NA - so maybe they didnt measure vol sometimes

ZP <- ZP %>% 
  drop_na(volume_1m2)

## ------------------------------------------ ##
#           1) Data Transformation -----
## ------------------------------------------ ##
# find the minimum non-zero value for each region and season

ZP <- ZP %>%
  group_by(region_name, season) %>%
  mutate(min_nonzero = min(volume_1m2[volume_1m2 > 0],
                           na.rm = TRUE)) %>%
  ungroup()

# Log transformation 
ZP$log10_ZP <- log10(ZP$volume_1m2 + ZP$min_nonzero/2)

ggplot(ZP, aes(x = year, y = log10_ZP)) +
  geom_point() +
  facet_grid(season~region_name, scales = "free")

## 1978-1988
## split into two 
ts1_zp <- ZP %>%
  filter(between(year, 1978, 1987))

ts2_zp <- ZP %>%
  filter(between(year, 1998, max(year)))

## ------------------------------------------ ##
#           2) Average -----
## ------------------------------------------ ##
# average across stations for a cruise or year/season

## ----
# Entire time series
## ----
ZP <- ZP %>%
  group_by(region_name, year, season) %>% 
  mutate(log_mean_ZP = mean(log10_ZP, na.rm = T)) %>%
  ungroup()

ZP_a <- ZP %>%
  group_by(region_name, year, season) %>%
  summarize(log_mean_ZP = mean(log10_ZP, na.rm = TRUE), .groups = "drop")

## ----
# Time series 1 = 1978-1987
## ----
ts1_zp <- ts1_zp %>%
  group_by(region_name, year, season) %>% 
  mutate(log_mean_ZP = mean(log10_ZP, na.rm = T)) %>%
  ungroup()

## ----
# Time series 2 = 1998-now
## ----
ts2_zp <- ts2_zp %>%
  group_by(region_name, year, season) %>% 
  mutate(log_mean_ZP = mean(log10_ZP, na.rm = T)) %>%
  ungroup()

ts2_zp_a <- ts2_zp %>%
  group_by(region_name, year, season) %>%
  summarize(log_mean_ZP = mean(log10_ZP, na.rm = TRUE), .groups = "drop")

##
ggplot(ts1_zp, aes(x = year, y = log_mean_ZP)) +
  geom_point() +
  facet_grid(season~region_name) 


## ------------------------------------------ ##
#    3) Running mean = 5 years  -----
## ------------------------------------------ ##
# take a running mean with timespan ~ longest lived taxon
# smoothing over life span of organisms 

## ----
# Entire time series 
## ----
rm_zp_tsFULL <- ZP %>%
  distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  mutate(runMean_ZP = runner::mean_run(log_mean_ZP, na_rm = TRUE, k=5))
  #mutate(running_mean = rollmean(log_mean_ZP, k = 5, fill = NA, align = "right"))

# same as above but less columns 
rm_zp_tsFULL_a <- ZP_a %>%
  #distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  #summarize(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k=5))
  mutate(runMean_ZP = runner::mean_run(log_mean_ZP, na_rm = TRUE, k = 5))

## ----
# Time series 1 = 1978-1987 
## ----
rm_zp_ts1 <- ts1_zp %>%
  distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>%
  mutate(runMean_ZP = runner::mean_run(log_mean_ZP, na_rm = TRUE, k=5))

## ----
# Time series 2 = 1998-now
## ----
rm_zp_ts2 <- ts2_zp %>%
  distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region_name, season, year) %>%
  group_by(region_name, season) %>%
  mutate(runMean_ZP = runner::mean_run(log_mean_ZP, na_rm = TRUE, k=5))

rm_zp_ts2_a <- ts2_zp_a %>%
  arrange(region_name, season, year) %>% 
  group_by(region_name, season) %>% 
  mutate(runMean_ZP = runner::mean_run(log_mean_ZP, na_rm = TRUE, k = 5))

## ------------------------------------------ ##
#    4) Compute SD   -----
## ------------------------------------------ ##
# compute st.dev. of time-series 

## ----
# Entire time series
## ----
ZP <- ZP %>%
  group_by(region_name, season) %>%
  mutate(sd_ZP = sd(log_mean_ZP, na.rm = T))
#write.csv(ZP, "output/trophamp_ZP_NES.csv")

#SD only
sd_ZP <- ZP %>% 
  select(region_name, season, sd_ZP) %>%
  distinct(region_name, season, .keep_all = T) %>%
  arrange(desc(region_name))
#write.csv(sd_ZP, "output/SD_trophamp_ZP_NES.csv")

distinct_tsFull_ZP <- ZP %>%
  distinct(region_name, season, year, .keep_all = T)

rm_zp_tsFull_FINAL <- rm_zp_tsFULL_a %>%
  left_join(distinct_tsFull_ZP %>% 
              select(-log_mean_ZP), 
            by = c("region_name", "season", "year")) %>%
  select(all_of(setdiff(names(.), names(rm_zp_tsFULL_a))), 
         everything()) %>%
  #reorder col
  relocate(region_name, .after = region) %>%
  relocate(year, .after = day) %>%
  relocate(season, .after = depth) %>%
  relocate(log_mean_ZP, .after = log10_ZP)
#write.csv(rm_zp_tsFull_FINAL, "output/trophamp_ZP_runmean_NES.csv")

## ----
# Time series 1 = 1978-1987
## ----
ts1_zp <- ts1_zp %>%
  #distinct(region_name, season, .keep_all = TRUE) %>%
  group_by(region_name, season) %>%
  mutate(SD_ZP = sd(log_mean_ZP, na.rm = T)) 

## ----
# Time series 2 = 1998-now
## ----
ts2_zp <- ts2_zp %>%
  group_by(region_name, season) %>%
  mutate(sd_ZP = sd(log_mean_ZP, na.rm = T)) 
#write.csv(ts2_zp, "output/trophamp_ZP_NES_1998_2021.csv")

distinct_ts2_zp <- ts2_zp %>%
  distinct(region_name, season, year, .keep_all = T)

rm_zp_ts2_FINAL <- rm_zp_ts2_a %>%
  left_join(distinct_ts2_zp %>% 
              select(-log_mean_ZP), 
            by = c("region_name", "season", "year")) %>%
  select(all_of(setdiff(names(.), names(rm_zp_ts2_a))), 
         everything()) %>%
  #reorder col
  relocate(region_name, .after = region) %>%
  relocate(year, .after = day) %>%
  relocate(season, .after = depth) %>%
  relocate(log_mean_ZP, .after = log10_ZP)
#write.csv(rm_zp_ts2_FINAL, "output/trophamp_ZP_runmean_NES_1998_2021.csv")


## ------------------------------------------ ##
#    Plots   -----
## ------------------------------------------ ##
ggplot() +
  geom_point(data = ZP, aes(x = year, y = log10_ZP)) +
  geom_point(data = ZP, aes(x = year, y = log_mean_ZP),
             color = "red", shape = 1, size = 2) +
  geom_line(data = rm_zp_tsFULL, aes(x = year, y = runMean_ZP), 
            color = "red", 
            size = 2) +
  ggtitle("ZP") + 
  facet_grid(season~region_name, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())

ggplot() +
  #geom_point(data = sub_ZP_v2, aes(x = year, y = log10Sum)) +
  #geom_point(data = ZP, aes(x = year, y = aver),
  #           color = "red", shape = 1, size = 2) +
  geom_line(data = ZP, aes(x = year, y = log_mean_ZP),
            color = "black") +
  geom_line(data = rm_zp_tsFULL, aes(x = year, y = runMean_ZP), 
            color = "red", 
            size = 2) +
  ggtitle("ZP") + 
  facet_grid(season~region_name, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())
