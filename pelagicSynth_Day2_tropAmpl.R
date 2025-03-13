#############################################################
#############      Pelagic Synthesis Meeting       ##########
#############                DAY 2                 ##########
#############             02-NOV-2023              ##########
#############         Trophic Amplification        ##########
#############################################################
## Packages
library(tidyverse)
library(runner) #for running mean calc

## Data
abu <- read.csv(file.path("raw","EcoMon_v3_8_wDateStrata.csv"))
#############################################################

## Larval Fish

# atlantic herring == Clupea harengus
# alewife == not in data = river h = alosa
# blueback herring == not in data = river h = alosa
# atlantic mackerel == scomber scombrus
# atlantic butterfish == peprilus spp
# sand lance == ammodytes spp

# right now looking at 4 larval fishes, potentially look at 
# other options. NOAA trawl data?? 
#ecodata::aggregate_biomass
#ecodata::seabird_ne

## Zooplankton calcs start on line # 

# pick columns needed 
subF <- abu[, c(1:21, 208, 229, 230, 235,
                254, 275, 276, 281)]

# long format - probably dont need 
#subF_l <- subF %>%
#  pivot_longer(
#    cols = cluhar_10m2:ammspp_100m3,
#    names_to = "taxa",
#    values_to = "ind")

# sum abundances of diff fishes to get a total of forage fish 
subF1 <- subF %>%
  group_by(cruise_name, station, region, date) %>% #grouping by sample
  mutate(fishsum = sum(cluhar_10m2, scosco_10m2,pepspp_10m2,ammspp_10m2, na.rm = T))

# add season to data 
subF1 <- subF1 %>%
  mutate(season = case_when(between(month, 3, 5) ~ "Spring",
                            between(month, 6, 8) ~ "Summer",
                            between(month, 9, 11) ~ "Fall",
                            TRUE ~ "Winter"))

# this takes a long long looooong time; check for alternatives 
subF1 <- subF1 %>%
  mutate(Region = case_when(region == 1 ~ "MAB",
                            region == 2 ~ "SNE",
                            region == 3 ~ "GB",
                            region == 4 ~ "GOM",
                            TRUE ~ "Outside"))

# remove the "outside" region 
subF2 <- subF1 %>% filter(region != 0)

####
ggplot(subF2, aes(x = year, y = fishsum)) +
  geom_point() +
  facet_wrap(~Region, scales = "free")

ggplot(subF2, aes(x = year, y = fishsum)) +
  geom_point() +
  facet_grid(season~Region, scales = "free")
 
#############################################################
######   ----     Trophic Amplification:
# 1) log10(x+ (min/2))    for each station ; x = sum trophic level 
# 2) average across stations for a cruise or year/season
# 3) take a running mean with timespan ~ longest lived taxon 
# 4) compute st.dev. of time-series 
#############################################################

#############################################################
## 1) Data Transformation 
# Find the minimum non-zero value for data transform for each region and season
subF2 <- subF2 %>%
  group_by(Region, season) %>%
  filter(fishsum != 0) %>% #exclude 0 
  mutate(min_nonzero = min(fishsum)) %>% 
  ungroup()

# Log transformation
subF2$LogofSum <- log(subF2$fishsum + subF2$min_nonzero/2)

ggplot(subF2, aes(x = year, y = LogofSum)) +
  geom_point() +
  facet_grid(season~Region, scales = "free")

###################################
## 1978-1988
## split into two time series/time chunks 
ts1 <- subF2 %>%
  filter(between(year, 1978, 1987))

ts2 <- subF2 %>%
  filter(between(year, 1998, max(year)))

#############################################################
## 2) Average

# Entire time series
tsFull <- subF2 %>%
  group_by(Region, year, season) %>% #by cruise too?no
  mutate(aver = mean(LogofSum, na.rm = T)) %>%
  ungroup()

# Time series 1 = 1978-1987
ts1 <- ts1 %>%
  group_by(Region, year, season) %>% #by cruise too?no
  mutate(aver = mean(LogofSum, na.rm = T)) %>%
  ungroup()

# Time series 2 = 1988-now
ts2 <- ts2 %>%
  group_by(Region, year, season) %>% #by cruise too?no
  mutate(aver = mean(LogofSum, na.rm = T)) %>%
  ungroup()

##
ggplot(ts1, aes(x = year, y = aver)) +
  geom_point() +
  facet_grid(season~Region) 

#############################################################
## 3) Running mean = 5 years 

#library(zoo)

# entire time series 
#rm_fish_tsFull <- tsFull %>%
#  distinct(year, Region, season, .keep_all = T) %>%
#  group_by(Region, season) %>%
#  arrange(year) %>%
#  mutate(running_mean = rollmean(aver, k = 5, fill = NA, align = "right"))

# first 4 yrs missing so trying to fix 
#result <- tsFull %>%
#  distinct(year, Region, season, .keep_all = T) %>%
#  group_by(Region, season) %>%
#  arrange(year) %>%
#  mutate(
#    running_mean = ifelse(year <= 1980, rollmean(aver, k = 4, fill = NA, align = "right"), 
#                          rollmean(aver, k = 5, fill = NA, align = "right")))

# Entire time series - NEED TO RENAME
result1 <- tsFull %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>% 
  mutate(runMean = runner::mean_run(aver, na_rm = TRUE, k=5))

# Time series 1 = 1978-1987 
#rm_fish_ts1 <- ts1 %>%
#  distinct(year, Region, season, .keep_all = T) %>%
#  group_by(Region, season) %>%
#  arrange(year) %>%
#  mutate(running_mean = rollmean(aver, k = 5, fill = NA, align = "right"))

rm_fish_ts1 <- ts1 %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>% 
  mutate(runMean = runner::mean_run(aver, na_rm = TRUE, k=5))

# Time series 2 = 1988-now
rm_fish_ts2 <- ts2a %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>%
  arrange(year) %>%
  mutate(running_mean = rollmean(aver, k = 5, fill = NA, align = "right"))

#############################################################
## 4) Compute SD 
## i need to do for timeseries

# Entire time series
tsFull_sd <- tsFull %>%
  #distinct(region, season, .keep_all = TRUE) %>%
  group_by(Region, season) %>%
  mutate(SD_ts = sd(aver, na.rm = T)) 

#organize print neatly 
tsFull_sd_1 <- tsFull_sd %>%
  select(Region, season, SD_ts) %>%
  distinct(Region, season, .keep_all = T) %>%
  arrange(desc(Region))

#write.csv(tsFull_sd_1, "output/SD_trophamp_FFISH.csv")
#tsFull_sd_2 <- tsFull_sd_1 %>% filter(Region != "Outside")
#write.csv(tsFull_sd_2, "output/SD_trophamp_FFISHv2.csv")

# Time series 1 = 1978-1987
ts1_sd <- ts1 %>%
  #distinct(region, season, .keep_all = TRUE) %>%
  group_by(Region, season) %>%
  mutate(SD_ts = sd(aver, na.rm = T)) 

#organize print neatly 
ts1_sd_1 <- ts1_sd %>% 
  select(Region, season, SD_ts) %>%
  distinct(Region, season, .keep_all = T) %>%
  arrange(desc(Region))

# Time series 2 = 1988-now --- CHECK!!!!
ts2_sd <- ts2 %>%
  distinct(region, year, season, .keep_all = TRUE) %>% # NO ?
  group_by(region, season) %>%
  mutate(SD_ts = sd(aver, na.rm = T)) 

############################################
ggplot(ts1sd, aes(x = year, y=SD_ts)) +
  geom_point(color = "red") +
  facet_grid(season~Region)

ggplot() +
  geom_point(data = ts1, aes(x = year, y = aver)) +
  geom_point(data = ts1sd, aes(x = year, y = SD_ts,color= "red")) +
  facet_grid(season~Region) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#####################################################
#####################################################
############## full plot - FULL TS
ggplot(subF2, aes(x = year, y = LogofSum)) +
  geom_point() +
  facet_grid(season~Region, scales = "free")

subF2 <- subF2 %>% filter(region != 0)
tsFull<- tsFull %>% filter(region != 0)
rm_fish_tsFull <- rm_fish_tsFull %>% filter(region != 0)

ggplot() +
  geom_point(data = subF2, aes(x = year, y = LogofSum)) +
  geom_point(data = tsFull, aes(x = year, y = aver),
             color = "red", shape = 1, size = 3) +
  geom_line(data = result1, aes(x = year, y = runMean), 
            color = "red", 
            size = 2) +
  ggtitle("FISH") + 
  facet_grid(season~Region, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())


ggplot() +
  #geom_point(data = subF2, aes(x = year, y = LogofSum)) +
  geom_line(data = tsFull, aes(x = year, y = aver),
             color = "black") +
  geom_line(data = result1, aes(x = year, y = runMean), 
            color = "red", 
            size = 2) +
  ggtitle("LARVAL FISH") + 
  facet_grid(season~Region, 
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


#############################################################
#############################################################
##             Zooplankton DISPLACEMENT VOLUME             ##
#############################################################

sub_ZP_v <- abu[, c(1:22)] #subset for columns 

sub_ZP_v <- sub_ZP_v %>%
  mutate(season = case_when(between(month, 3, 5) ~ "Spring",
                            between(month, 6, 8) ~ "Summer",
                            between(month, 9, 11) ~ "Fall",
                            TRUE ~ "Winter"))

sub_ZP_v <- sub_ZP_v %>%
  mutate(Region = case_when(region == 1 ~ "MAB",
                            region == 2 ~ "SNE",
                            region == 3 ~ "GB",
                            region == 4 ~ "GOM",
                            TRUE ~ "Outside"))

# remove the "outside" region 
#sub_ZP_v <- sub_ZP_v %>% filter(region != 0)

#############################################################
## 1) Data Transformation 
# Find the minimum non-zero value for data transform for each region and season
sub_ZP_v2 <- sub_ZP_v %>%
  group_by(Region, season) %>%
  filter(volm2sum != 0) %>%
  mutate(min_nonzero = min(volm2sum)) %>%
  ungroup()

# Log transformation 
sub_ZP_v2$LogofSum <- log(sub_ZP_v2$volm2sum + sub_ZP_v2$min_nonzero/2)

ggplot(sub_ZP_v2, aes(x = year, y = LogofSum)) +
  geom_point() +
  facet_grid(season~Region, scales = "free")

## 1978-1988
## split into two 
ts1v <- sub_ZP_v2 %>%
  filter(between(year, 1978, 1987))

ts2v <- sub_ZP_v2 %>%
  filter(between(year, 1998, max(year)))

#############################################################
## 2) Average 

# Entire time series
tsFULLv <- sub_ZP_v2 %>%
  group_by(Region, year, season) %>% #by cruise too?no
  mutate(aver = mean(LogofSum, na.rm = T)) %>%
  ungroup()

# Time series 1 = 1978-1987
ts1v <- ts1v %>%
  group_by(Region, year, season) %>% #by cruise too?no
  mutate(aver = mean(LogofSum, na.rm = T)) %>%
  ungroup()

# Time series 2 = 1988-now
ts2av <- ts2v %>%
  group_by(Region, year, season) %>% #by cruise too?no
  mutate(aver = mean(LogofSum, na.rm = T)) %>%
  ungroup()

##
ggplot(ts1v, aes(x = year, y = aver)) +
  geom_point() +
  facet_grid(season~Region) 

#############################################################
## 3) Running mean = 5 years
# smoothing over life span of organisms 

# Entire time series
rm_result_tsFULLv <- tsFULLv %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>%
  arrange(year) %>%
  mutate(running_mean = rollmean(aver, k = 5, fill = NA, align = "right"))

#actual running mean without loosing first 4 yrs 
result2 <- tsFULLv %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>% 
  mutate(runMean = runner::mean_run(aver, na_rm = TRUE, k=5))

# Time series 1 = 1978-1987
rm_result_ts1 <- ts1v %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>%
  arrange(year) %>%
  mutate(running_mean = rollmean(aver, 
                                 k = 5,
                                 fill = NA, 
                                 align = "right"))

# Time series 2 = 1988-now
rm_result_ts2 <- ts2av %>%
  distinct(year, Region, season, .keep_all = T) %>%
  group_by(Region, season) %>%
  arrange(year) %>%
  mutate(running_mean = rollmean(aver, k = 5, fill = NA, align = "right"))

#############################################################
## 4) Compute SD

# Entire time series
tsFULLv_SD <- tsFULLv %>%
  #distinct(region, season, .keep_all = TRUE) %>%
  group_by(Region, season) %>%
  mutate(SD_ts = sd(aver, na.rm = T)) 

tsFULLvSD1 <- tsFULLv_SD %>% 
  select(Region, season, SD_ts) %>%
  distinct(Region, season, .keep_all = T) %>%
  arrange(desc(Region))
#write.csv(tsFULLvSD1, "output/SD_trophamp_ZP.csv")

# Time series 1 = 1978-1987
ts1sdavol <- ts1v %>%
  #distinct(region, season, .keep_all = TRUE) %>%
  group_by(Region, season) %>%
  mutate(SD_ts = sd(aver, na.rm = T)) 

ts1sda1vol <- ts1sdavol %>% 
  select(Region, season, SD_ts) %>%
  distinct(Region, season, .keep_all = T) %>%
  arrange(desc(Region))

# Time series 2 = 1988-now
ts2sdavol <- ts2av %>%
  #distinct(region, season, .keep_all = TRUE) %>%
  group_by(Region, season) %>%
  mutate(SD_ts = sd(aver, na.rm = T)) 

ts2sdavol <- ts2sdavol %>% 
  select(Region, season, SD_ts) %>%
  distinct(Region, season, .keep_all = T) %>%
  arrange(desc(Region))



#####
#sub_ZP_v2 <- sub_ZP_v2 %>% filter(region != 0)
#tsFULLv<- tsFULLv %>% filter(region != 0)
#rm_result_tsFULLv <- rm_result_tsFULLv %>% filter(region != 0)

ggplot() +
  geom_point(data = sub_ZP_v2, aes(x = year, y = LogofSum)) +
  geom_point(data = tsFULLv, aes(x = year, y = aver),
             color = "red", shape = 1, size = 2) +
  geom_line(data = rm_result_tsFULLv, aes(x = year, y = running_mean), 
            color = "red", 
            size = 2) +
  ggtitle("ZP") + 
  facet_grid(season~Region, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())

ggplot() +
  #geom_point(data = sub_ZP_v2, aes(x = year, y = LogofSum)) +
  #geom_point(data = tsFULLv, aes(x = year, y = aver),
  #           color = "red", shape = 1, size = 2) +
  geom_line(data = tsFullv, aes(x = year, y = aver),
            color = "black") +
  geom_line(data = rm_result_tsFULLv, aes(x = year, y = running_mean), 
            color = "red", 
            size = 2) +
  ggtitle("ZP") + 
  facet_grid(season~Region, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())








#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
cal <- ecodata::calanus_stage %>% 
  dplyr::filter(EPU == "MAB") %>% 
  tidyr::separate(Var, into = c("Var", "season"), sep = "-") %>% 
  filter(Var %in% c("c3", "c4", "c5", "adt"))

#############################################################cal$Var <- factor(cal$Var, levels = c("c3", "c4", "c5", "adt"))
cal$season <- factor(cal$season, levels = c("Spring", "Summer", "Fall"))

cal %>% 
  ggplot2::ggplot(aes(x = Time, y = Value, color = Var, fill = Var)) +
  ggplot2::annotate("rect", fill = shade.fill, alpha = shade.alpha,
                    xmin = x.shade.min , xmax = x.shade.max, #issue here ! 
                    ymin = -Inf, ymax = Inf)+
  ggplot2::geom_bar(stat = "identity")+
  ggplot2::facet_wrap(~season)+
  ggplot2::ylab("Calanus Stage (N/100m^3)") +
  ggplot2::xlab(element_blank())+
  ggplot2::ggtitle("MAB Calanus Stage Abundance") +
  ggplot2::theme(legend.position = "bottom", 
                 legend.title = element_blank())+
  ecodata::theme_facet()+
  scale_fill_manual(values = c("steelblue1","steelblue3", "coral1", "coral3"))+
  scale_color_manual(values = c("steelblue1","steelblue3", "coral1", "coral3"))+
  ecodata::theme_title()

bird <- ecodata::seabird_ne
help("seabird_ne")
ggplot(bird, aes(x = Time, Value)) +
  geom_point()

ch <- ecodata::chl_pp


## turn these into loops - the plots 
#abuF <- abu[, c(1:21, 205:251)]

fishL <- abuF %>%
  pivot_longer(
    cols = pnepau_100m3:lopame_10m2,
    names_to = "taxa",
    values_to = "ind"
  )

ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("bretyr_10m2"),],
             aes(x = year, y=ind))

cluhar_10m2 <- fishL %>% filter(taxa == "cluhar_10m2")

cluhar_10m2$r4r <- cluhar_10m2$ind^(1/4)

ggplot() + 
  geom_point(data=cluhar_10m2,
             aes(x = year, y=r4r))

pepspp_10m2 <- fishL %>% filter(taxa == "pepspp_10m2")

pepspp_10m2$r4r <- pepspp_10m2$ind^(1/4)

ggplot() + 
  geom_point(data=pepspp_10m2,
             aes(x = year, y=ind))


ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("cluhar_10m2"),],
             aes(x = year, y=ind))

#looks good 
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("cycspp_10m2"),],
             aes(x = year, y=ind))

#maybe
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("diaspp_10m2"),],
             aes(x = year, y=ind))


ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("cermad_10m2"),],
             aes(x = year, y=ind))

ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("benspp_10m2"),],
             aes(x = year, y=ind))

#pretty abundance missing a couple of yrs 
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("urospp_10m2"),],
             aes(x = year, y=ind))

#abundant after 2000
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("enccim_10m2"),],
             aes(x = year, y=ind))

#abundant throughout
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("gadmor_10m2"),],
             aes(x = year, y=ind))


ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("melaeg_10m2"),],
             aes(x = year, y=ind))

#no
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("polvir_10m2"),],
             aes(x = year, y=ind))

ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("meralb_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("merbil_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("centstr_10m2"),],
             aes(x = year, y=ind))

ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("pomsal_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("cynreg_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("leixan_10m2"),],
             aes(x = year, y=ind))

ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("menspp_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("micund_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("tauads_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("tauoni_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("auxspp_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("scosco_10m2"),],
             aes(x = year, y=ind))
#maybe
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("pepspp_10m2"),],
             aes(x = year, y=ind))

#yes- rockfishes
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("sebspp_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("prispp_10m2"),],
             aes(x = year, y=ind))
#maybe
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("myoaen_10m2"),],
             aes(x = year, y=ind))
#sculpin a lot 
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("myooct_10m2"),],
             aes(x = year, y=ind))
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("prispp_10m2"),],
             aes(x = year, y=ind))

#sand lance 
ggplot() + 
  geom_point(data=fishL[fishL$taxa %in% c("ammspp_10m2"),],
             aes(x = year, y=ind))

#"phogun_10m2" "ulvsub_10m2"  "anaspp_10m2"  "citarc_10m2"  "etrspp_10m2" 
#"syaspp_10m2"  "botspp_10m2"  "hipobl_10m2"  "parden_10m2" 
#"pseame_10m2"  "hippla_10m2"  "limfer_10m2"  "glycyn_10m2" 
#"scoaqu_10m2"  "sypspp_10m2"  "lopame_10m2" 