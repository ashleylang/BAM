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
library(car)
library(ggpubr)
library(ggmap)
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
                               pft=="DA" ~ "deciduous"),
         evo_group= case_when(pft=="EA" ~ "angiosperm",
                              pft== "EG" ~ "gymnosperm",
                              pft=="DA" ~ "angiosperm") )

#check distributions
hist(full_df$LmTm)
hist(full_df$RmTm)
hist(full_df$LMRM)#this needs to be log transformed
hist(full_df$h.t) #this needs to be log transformed
plot(full_df$LmTm~full_df$RmTm)


#summarise group numbers for pfts, ecosystems, myc types
sub <- full_df %>%
  group_by(Temp, Prec, myc_group, leaf_habit, latitude, longitude) %>%
  summarise(n=n())


###Making some plots:----
AM_ECM=c( '#8597FE', '#7CAE31')
#Any strong correlations between continuous variables?
full_df_cor=full_df %>% 
  dplyr::select(Temp, Prec, LmTm, RmTm,  log_LMRM, log_ht) %>% 
  drop_na()

corrplot::corrplot(cor(full_df_cor), method="number", type="upper")
#Negative correlation with height and LmTm, height and RmTm (smaller trees have more leaves *and* roots relative to trunk)
#Positive corr. between Precip and Temp and LmTm and Temp (leaf mass is a larger component of total mass in hotter places-- & here hotter also seems to mean wetter but less strong corr. with precip)

# plot myc group on climate space- check for confounding? 
# Color of myc type, shape for leaf habit, plotted on MAP/MAT

clim_space=ggplot(sub,aes(x= Temp, y=Prec))+
  geom_point(aes(colour=myc_group,shape=leaf_habit, size=n), alpha=0.65)+
  scale_colour_manual(labels = c("AM", "ECM"),values=AM_ECM)+
  scale_shape_manual(labels = c("Deciduous", "Evergreen"), values=c(16,17))+
  labs(x=expression("Mean Annual Temperature ("*degree*C*")"), y="Mean Annual Precipitation (mm)")+
  theme(legend.title=element_blank(),legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(0, "mm"), 
        legend.text=element_text(size=8), legend.position = c(.04, .85), 
        axis.title=element_text(size=11), axis.text=element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size=3), order=1), shape = guide_legend(override.aes = list(size=4), order=2))+
  scale_size(range = c(3,8), guide="none")
clim_space
#ggsave("Figure_1_climate.pdf", path="~/BAM", width= 88, height= 180, units="mm")
#^ Was trying to figure out how to save this to our GitHub repository; failed

#map
world <- map_data("world")

map=ggplot(data=world)+
  geom_polygon(aes(x=long, y=lat, group=group), fill="white", color="black", size=0.2)+
  coord_fixed(1.3)+
  theme(axis.text=element_blank(), axis.line = element_blank(), 
        axis.ticks=element_blank(), axis.title=element_blank(), 
        panel.background = element_rect(fill = "white"), 
        #panel.border = element_rect(linetype="solid", fill=NA), legend.position="none")+
        panel.border = element_blank(), legend.position="none")+
  geom_point(aes(x = longitude, y = latitude,color=myc_group), data = sub, size = 1)+
  scale_colour_manual(values=AM_ECM)
map

ggarrange(map, clim_space, nrow=1, ncol=2, labels=c("a", "b")) #may look odd with diff. screen widths; sized for saving as pdf using ggsave below.

ggsave("Figure_1.pdf", height=80,width=180, units="mm")


####models-----
#make LMMS
#run models with standardized and unstandardized coefficients
#use ggeffects to pull out effect of myc type

##Root mass/total mass model
R_full_model <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df)
vif(R_full_model)
summary(R_full_model)
tab_model(R_full_model, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#remove insignificant interaction terms: log_ht*myc_group, Temp*myc_group, log_ht*Temp, Temp*leaf_habit (marginally significant)
R_reduced_m1 <-lmer(RmTm~log_ht*leaf_habit + myc_group + Temp + Prec  + (1|study_species), data = full_df)
vif(R_reduced_m1)
summary(R_reduced_m1)
tab_model(R_reduced_m1, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)

#use ggeffect to calculate the marginal effects of myc group and height 
R_E1 <- ggeffect(R_reduced_m1, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")

##Leaf mass/total mass models
#full model (all terms and interactions)
L_full_model <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df)
vif(L_full_model)
summary(L_full_model)

#reduced model: remove Temp*leaf_habit, log_ht*Temp, myc_group*Temp
L_reduced_m1 <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + (1|study_species), data = full_df)
vif(L_reduced_m1)
summary(L_reduced_m1)
tab_model(L_reduced_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#use ggeffect to calculate the marginal effects of myc group and height 
L_E1 <- ggeffect(L_reduced_m1, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")

#Leaf:root full model
B_full_model <-lmer(log_LMRM ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df)
vif(B_full_model)
summary(B_full_model)

#reduced model: remove myc_group*temp and leafhabit*Temp
B_reduced_m1 <-lmer(log_LMRM ~ log_ht*leaf_habit + log_ht*myc_group +  log_ht*Temp  + Prec + (1|study_species), data = full_df)
vif(B_reduced_m1)
summary(B_reduced_m1)
tab_model(B_reduced_m1, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#use ggeffect to calculate the marginal effects of myc group and height 
B_E1 <- ggeffect(B_reduced_m1, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")

#Figure 2: run through ggarrange
#plot the marginal effects and the raw data for Rm/Tm
a <- ggplot() +
  geom_point(data = full_df, aes(x = log_ht, y = RmTm, color = myc_group, shape=leaf_habit), alpha = .2, size=3) +
  geom_line(data = R_E1, aes(x = x, y = predicted, color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal Type", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17))+
  #  theme(legend.position = c(.65, .75), legend.title= element_text(hjust=0.5))+
  theme(legend.position = c(.45, .55), legend.title= element_text(hjust=0.5))+
  labs(x= "ln(Height)", y="Root mass / Total mass" )+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))

#plot the marginal effects and the raw data for Lm/Tm
b <- ggplot() +
  geom_point(data = full_df, aes(x = log_ht, y = LmTm, color = myc_group, shape=leaf_habit), alpha = .2, size=3) +
  geom_line(data = L_E1, aes(x = x, y = predicted, color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal Type", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17)) +
  theme(legend.position = "none")+
  labs(x= "ln(Height)", y="Leaf mass / Total mass" )+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))
  

#plot the marginal effects and the raw data for Lm/Rm
c <- ggplot() +
  geom_point(data = full_df, aes(x = log_ht, y = log_LMRM, color = myc_group, shape=leaf_habit), alpha = .2, size=3) +
  geom_line(data = B_E1, aes(x = x, y = predicted, color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal Type", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17)) +
  theme(legend.position = "none", axis.title.y=element_text(size=12.5))+
  labs(x= "ln(Height)", y="ln(Leaf mass / Root mass)" )+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))
  

#cowplot::plot_grid(a, b, c, nrow = 1, labels="auto")

#here's a way to get the three figs and legend as their own quadrants of a square:
leg=get_legend(a)
a=a+theme(legend.position="none")
ggarrange(a, b, c, leg, labels=c("a", "b", "c", " "), nrow=2, ncol=2)
#ggsave("Figure_2.pdf", path="/Users/ashleylang/Documents/GitHub/BAM/", height=160,width=180, units="mm")
