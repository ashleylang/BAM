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
library(raster)
library(sp)
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

#now summarise to the speceis level: chose the myc 
fungal_root_sp <- fungal_root %>%
  group_by(order, family, genus, species, myc_group) %>%
  summarise(number= n()) %>%
  slice_max(order_by = number, n = 1)

#create df based on FungalRoot that has the most common mycorrhizal type associated with each genus of plant (to use in cases where there is no obs in Fugalroot for a speceis in BAAD)
fungal_root_genus <- fungal_root %>%
  group_by(order, family, genus, myc_group) %>%
  summarise(number= n()) %>%
  slice_max(order_by = number, n = 1)


#get allometry database into a dataframe & add mycorrhizal associations
baad <- baad.data::baad_data()
dict <- as.data.frame(baad$dictionary)
baad_df <- as.data.frame(baad$data) %>% #Q from Ashley: What is the difference between species and speciesMatched?
  dplyr::select(studyName, latitude, longitude, species, vegetation, map, mat, pft, a.lf, h.t, d.bh, m.lf, m.st, m.so, m.rt, m.to, 	ma.ilf) %>%
  separate(species, c("genus", "species"), extra = "drop", fill = "right") %>%
  left_join(fungal_root_sp, by = c("genus", "species")) %>% 
  mutate(LmSm = m.lf/ m.st,
         LmSo = m.lf/ m.so,
         LmTm = m.lf/ m.to,
         LaSm = a.lf/ m.st,
         LmLa = m.lf/ a.lf,
         RmTm = m.rt/m.to,
         pft = as.factor(pft)) %>% 
   filter( myc_group != "ERC" & myc_group != "Other" & pft != "DG")# %>%
#   dplyr::select(studyName, pft, Temp, mat, Prec, vegetation, myc_group, h.t, m.so, m.to, m.rt, m.lf, LmSo, LmTm, LmSm, LaSm, LmLa, RmTm) 

# sub <- bigdf %>%
#   group_by(pft) %>%
#   summarise(n_distinct(species)) 

###MAT/MAP----
###Need worldclim data because most of the BAAD database don't have MAT/MAP
library(raster)
library(sp)

r <- getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")

coords <- data.frame(x=baad_df$longitude, y=baad_df$latitude) %>% 
  na.omit() %>%
  distinct()

points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points)
df <- cbind.data.frame(coordinates(points),values)%>%
  rename(latitude = y, longitude = x)

#sooo....what are the units of temp??
#Also some missing values for temp/precip that maybe we should back-fill.

#Now the full version of the data with baad_df, myc types, and climate:
full_df = baad_df %>% 
  left_join(df, by=c("latitude", "longitude"))

###Making some plots:----
  
a <- ggplot(sub, aes(x = log(h.t), y = LmTm)) +
  geom_point(aes(color = myc_group)) +
  geom_smooth(aes(color = myc_group)) +
  facet_grid(.~ pft, scales = "free") +
  labs(x = "MAT (c)", y = "Root biomass/Total biomass") +
  theme(legend.position = "none")

b <- ggplot(sub, aes(x = Prec, y = RmTm)) +
  geom_point(aes(color = myc_group)) +
  geom_smooth(aes(color = myc_group), method = "lm") +
  #facet_grid(.~ pft, scales = "free") +
  labs(x = "MAP (mm)", y = "Root biomass/Total biomass") +
  theme(legend.position = "none")

plot_grid(a, b, nrow = 1)

summary(lm(RmTm ~ log(h.t) + myc_group*Prec +Temp*Prec, data = sub))

ggplot(sub, aes(x = Temp, y = log(h.t))) +
  geom_point(aes(color = myc_group)) +
  geom_smooth(aes(color = myc_group), method = "lm") 

c <- ggplot(sub, aes(x = log(h.t), y = LmTm)) +
  geom_point(aes(color = myc_group)) +
  geom_smooth(aes(color = myc_group), method = "lm") +
  labs(x = "Log(tree height)", y = "Leaf mass/Total mass") +
  theme(legend.position = "none")

d <- ggplot(sub, aes(x = log(h.t), y = RmTm)) +
  geom_point(aes(color = myc_group)) +
  geom_smooth(aes(color = myc_group), method = "lm") +
  labs(x = "Log(tree height)", y = "Root mass/Total mass") +
  theme(legend.position = "none")

plot_grid(c,d, nrow = 1)


