##allometry of AM vs EMC trees
###using the BAAD dataset (Falster et al 2015 Ecology)
##Written by Ashley Lang and Fiona Jevon

#only need to do this once:
#install.packages("devtools")
#devtools::install_github("richfitz/datastorr")
#devtools::install_github("traitecoevo/baad.data")

#load libraries
library(baad.data)
library(tidyverse)
library(cowplot)
library(sp)
library(raster)
library(corrplot)
library(lme4)
library(lmerTest)
library(sjPlot)
library(ggeffects)
library(arm)
library(MuMIn)
theme_set(theme_cowplot())

###get FUNGALROOT data----
FR_occurrences <- read_csv("FR_occurrences.csv")
FR_measurements <- read_csv("FR_measurements.csv") %>%
  pivot_wider(names_from = measurementType, values_from = measurementValue) %>% 
  rename(Myc_type='Mycorrhiza type')

#summarise Fungalroot database  to get the number of observations of different myc types for each species
fungal_root <- left_join(FR_occurrences, FR_measurements, by = "CoreID") %>%
  dplyr::select(order, family, scientificName, Myc_type) %>%
  separate(scientificName, c("genus", "species"), extra = "drop", fill = "right") %>%
  group_by(order, family, genus, species, Myc_type) %>%
  summarise(n = n())

#make Myc type groupings: 
unique(fungal_root$Myc_type)

ECMS = c("EcM, AM undetermined", "EcM, no AM colonization")
ERCS = c("ErM, AM", "ErM, EcM", "ErM")

fungal_root$myc_group <- ifelse(fungal_root$Myc_type == "AM", "AM", ifelse(fungal_root$Myc_type == "EcM,AM", "ECM/AM", ifelse(fungal_root$Myc_type %in% ECMS, "ECM", ifelse(fungal_root$Myc_type %in% ERCS, "ERC", "Other"))))

#now summarise to the speceis level: choose the most common mycorrhizal type associated with each species of plant
fungal_root_sp <- fungal_root %>%
  group_by(order, family, genus, species, myc_group) %>%
  summarise(number= sum(n)) %>%
  slice_max(order_by = number, n = 1) %>% # need to fix ties
  arrange(genus, species) %>% 
  filter(myc_group != "ECM/AM" & myc_group != "Other" & myc_group != "ERC") %>% 
  group_by(order, family, genus, species) %>% 
  mutate(dupe = n()>1) %>% 
  filter(dupe==F)

duplicates= fungal_root_sp %>% 
  filter(dupe==TRUE)#Now we have 83 species with equal #s of AM and ECM observations
#need to come up with a way of dealing with speceis that have ties: using genus level associations?
#Problem with using Genus level associations: many are from genera with mixtures of AM and ECM species or dual symbionts
#My thought is...treat them how we're treating dual symbionts (remove from data set). I have done this above with the last filter in that pipe for fungal_root_sp

#create df based on FungalRoot that has the most common mycorrhizal type associated with each genus of plant 
# fungal_root_genus <- fungal_root %>%
#   group_by(order, family, genus, myc_group) %>%
#   summarise(number= n()) %>%
#   slice_max(order_by = number, n = 1)

#get allometry database into a dataframe & add mycorrhizal associations
baad <- baad.data::baad_data()
dict <- as.data.frame(baad$dictionary)
baad_df <- as.data.frame(baad$data) %>% #Q from Ashley: What is the difference between species and speciesMatched? - see dict file: it looks like it is just a checking column for comparing this to other databases. I switched to using the regular species column here (I don't think it matters)  
  dplyr::select(studyName, latitude, longitude, species, vegetation, map, mat, pft, a.lf, h.t, d.bh, m.lf, m.st, m.so, m.rt, m.to, 	ma.ilf) %>%
  separate(species, c("genus", "species"), extra = "drop", fill = "right") %>%
  left_join(fungal_root_sp, by = c("genus", "species")) %>% 
  mutate(LmTm = m.lf/ m.to,
         RmTm = m.rt/m.to,
         LMRM = m.lf/m.rt,
         pft = as.factor(pft)) %>% 
   filter( myc_group != "ERC" & myc_group != "Other" & pft != "DG")# %>%
#   dplyr::select(studyName, pft, Temp, mat, Prec, vegetation, myc_group, h.t, m.so, m.to, m.rt, m.lf, LmSo, LmTm, LmSm, LaSm, LmLa, RmTm) 

sum <- baad_df %>%
  group_by(genus, species, myc_group) %>%
  summarise(n = n()) %>%
  arrange(genus, species)


###MAT/MAP----
###Need worldclim data because most of the BAAD database don't have MAT/MAP
r <- raster::getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")

coords <- data.frame(x=baad_df$longitude, y=baad_df$latitude) %>%
  na.omit() %>%
  distinct()

points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points)
df <- cbind.data.frame(coordinates(points),values)%>%
  rename(latitude = y, longitude = x) %>% 
  mutate(Temp=Temp/10)#units: MAT is now in deg C. Precip is in mm

#Also some missing values for temp/precip that maybe we should back-fill.

###filtering----
#Now the full version of the data with baad_df, myc types, and climate:
full_df = baad_df %>% 
  left_join(df, by=c("latitude", "longitude")) %>%
  filter(RmTm>0, h.t >.5, myc_group != "ECM/AM") %>% 
  unite(study_species, studyName, genus, species, sep="_", remove=F) %>% 
  mutate(log_ht = log(h.t),
         log_LMRM = log(LMRM),
         leaf_habit= case_when(pft=="EA" | pft== "EG" ~ "evergreen",
                               pft=="DA" ~ "deciduous") )

#check distributions
hist(full_df$LmTm)
hist(full_df$RmTm)
hist(full_df$LMRM)#this needs to be log transformed
hist(full_df$h.t) #this needs to be log transformed
plot(full_df$LmTm~full_df$RmTm)


#summarise group numbers for pfts, ecosystems, myc types
sub <- full_df %>%
  group_by(pft, myc_group, vegetation) %>%
  summarise(n = n())


###Making some plots:----
ggplot(full_df, aes(x = log(h.t), y = RmTm)) +
  geom_point(aes(color = myc_group)) +
  geom_smooth(aes(color = myc_group), method = "gam") +
  facet_grid(.~ pft, scales = "free") 

ggplot(full_df, aes(x = myc_group, y = RmTm)) +
geom_point(aes(color = myc_group), size=3, alpha=0.5)+
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.8, position="dodge")+
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.3, position = position_dodge(width = 0.8))+
  facet_grid(.~ pft, scales = "free") 

#Any strong correlations between continuous variables?
full_df_cor=full_df %>% 
  dplyr::select(Temp, Prec, LmTm, RmTm,  log_LMRM, log_ht) %>% 
  drop_na()

corrplot::corrplot(cor(full_df_cor), method="number", type="upper")
#Negative correlation with height and LmTm, height and RmTm (smaller trees have more leaves *and* roots relative to trunk)
#Positive corr. between Precip and Temp and LmTm and Temp (leaf mass is a larger component of total mass in hotter places-- & here hotter also seems to mean wetter but less strong corr. with precip)

####models-----
#make LMMS
#make separate models for each pft
#temperature not interactive
#standardize the continuous variables
#use ggeffects to pull out effect of myc type

##root mass/total mass models
#Evergreen Angiosperms
R_EA_m1 <-lmer(RmTm~log_ht + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "EA"))
summary(R_EA_m1)
r.squaredGLMM(R_EA_m1)
tab_model(R_EA_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#use ggpredict to plot just the effects of myc group (basically removing height)
plot(ggpredict(R_EA_m1, terms = c("Temp", "myc_group")))
ggplot(data = subset(full_df, pft == "EA"), aes(x = Temp, y = RmTm)) + geom_point(aes(color = myc_group))

#use ggpredict to plot just the effects of myc group (basically removing temp)
plot(ggpredict(R_EA_m1, terms = c("log_ht", "myc_group")))
ggplot(data = subset(full_df, pft == "EA"), aes(x = log_ht, y = RmTm)) + geom_point(aes(color = myc_group))

#use ggpredict to plot just the effects of myc group (remove Temp and height)
plot(ggpredict(R_EA_m1, terms = c("myc_group")))
ggplot(data = subset(full_df, pft == "EA"), aes(x = myc_group, y = RmTm)) + geom_boxplot(aes(color = myc_group))


#Deciduous Angiosperms
R_DA_m1 <-lmer(RmTm~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "DA"))
summary(R_DA_m1)
r.squaredGLMM(R_DA_m1)
tab_model(R_DA_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#use ggpredict to plot just the effects of myc group (basically removing height)
plot(ggpredict(R_DA_m1, terms = c("Temp", "myc_group")))
ggplot(data = subset(full_df, pft == "DA"), aes(x = Temp, y = RmTm)) + geom_point(aes(color = myc_group))

#use ggpredict to plot just the effects of myc group (remove Temp and height)
plot(ggpredict(R_DA_m1, terms = c("myc_group")))

#Evergreen Gymnosperms
R_EG_m1 <-lmer(RmTm~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "EG"))
summary(R_EG_m1)
r.squaredGLMM(R_EG_m1)
tab_model(R_EG_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#use ggpredict to plot just the effects of myc group (basically removing height)
plot(ggpredict(R_EG_m1, terms = c("Temp", "myc_group")))
ggplot(data = subset(full_df, pft == "EG"), aes(x = Temp, y = RmTm)) + geom_point(aes(color = myc_group))

#use ggpredict to plot just the effects of myc group (remove Temp and height)
plot(ggpredict(R_EG_m1, terms = c("myc_group")))


##leaf mass/total mass models
L_EA_m1 <-lmer(LmTm~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "EA"))
tab_model(L_EA_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")
plot(ggpredict(L_EA_m1, terms = c("Temp", "myc_group")))

L_DA_m1 <-lmer(LmTm~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "DA"))
tab_model(L_DA_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")
plot(ggpredict(L_DA_m1, terms = c("Temp", "myc_group")))

L_EG_m1 <-lmer(LmTm~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "EG"))
tab_model(L_EG_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")
plot(ggpredict(L_EG_m1, terms = c("Temp", "myc_group")))

#Leaf:root models
B_EA_m1 <-lmer(log(LMRM)~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "EA"))
tab_model(B_EA_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")
plot(ggpredict(B_EA_m1, terms = c("Temp", "myc_group")))


B_DA_m1 <-lmer(log(LMRM)~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "DA"))
tab_model(B_DA_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")
plot(ggpredict(B_DA_m1, terms = c("Temp", "myc_group")))

B_EG_m1 <-lmer(log(LMRM)~log(h.t) + myc_group + Temp + (1|study_species), data = subset(full_df, pft == "EG"))
tab_model(B_EG_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")
plot(ggpredict(B_EG_m1, terms = c("Temp", "myc_group")))
