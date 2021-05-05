##allometry of AM vs EMC trees
###using the BAAD dataset (Falster et al 2015 Ecology)

#only need to do this once:
install.packages("devtools")
devtools::install_github("richfitz/datastorr")
devtools::install_github("traitecoevo/baad.data")

library(baad.data)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

baad <- baad.data::baad_data()
dict <- as.data.frame(baad$dictionary)
baad_df <- as.data.frame(baad$data)

baad_df <- baad_df %>%
  dplyr::select(latitude, longitude, species, vegetation, map, mat, speciesMatched, pft, a.lf, h.t, d.bh, m.lf, m.st, m.so, m.rt, m.to, 	ma.ilf) %>%
  separate(col = speciesMatched, into = c("genus", "sp"), sep = " ")


###get FUNGALROOT data----
FR_occurrences <- read_csv("FR_occurrences.csv")
FR_measurements <- read_csv("FR_measurements.csv") %>%
  pivot_wider(names_from = measurementType, values_from = measurementValue) 

colnames(FR_measurements)[3] <- "Myc_type"

fungal_root <- left_join(FR_occurrences, FR_measurements, by = "CoreID") %>%
  select(order, family, genus, scientificName, Myc_type) %>%
  group_by(order, family, genus, scientificName, Myc_type) %>%
  summarise(n = n())

ECMS = c("EcM, AM undetermined", "EcM, no AM colonization")
ERCS = c("ErM, AM", "ErM, EcM", "ErM")

fungal_root$myc_group <- ifelse(fungal_root$Myc_type == "AM", "AM", ifelse(fungal_root$Myc_type == "ECM,AM", "ECM,AM", ifelse(fungal_root$Myc_type %in% ECMS, "ECM", ifelse(fungal_root$Myc_type %in% ERCS, "ERC", "Other"))))

fungal_root_family <- fungal_root %>%
  group_by(order, family, myc_group) %>%
  summarise(number= n()) %>%
  slice_max(order_by = number, n = 1)

fungal_root_genus <- fungal_root %>%
  group_by(order, family, genus, myc_group) %>%
  summarise(number= n()) %>%
  slice_max(order_by = number, n = 1)

###MAT/MAP----
###also need worldclim data because most of the BAAD database don't have MAT/MAP
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
df <- cbind.data.frame(coordinates(points),values)
df <- df %>%
  rename(latitude = y, longitude = x)

baad_df2 <- left_join(baad_df, df)

##add myc typ from FR genus to BAAD----
bigdf <- left_join(baad_df2, fungal_root_genus, by = "genus")
bigdf$LmSm = bigdf$m.lf/ bigdf$m.st
bigdf$LmSo = bigdf$m.lf/ bigdf$m.so
bigdf$LmTm = bigdf$m.lf/ bigdf$m.to
bigdf$LaSm = bigdf$a.lf/ bigdf$m.st
bigdf$LmLa = bigdf$m.lf/bigdf$a.lf
bigdf$RmTm = bigdf$m.rt/bigdf$m.to
bigdf$pft = as.factor(bigdf$pft)


forests <- c("TropRF", "TropSF", "TempRF", "TempF", "BorF")

bigdf$FT <- ifelse(bigdf$vegetation == "BorF", "Boreal", ifelse(bigdf$vegetation %in% c("TempF", "TempRF"), "Temperate", ifelse(bigdf$vegetation %in% c("TropRF", "TropSF"), "Tropical", "nonF")))

sub <- bigdf %>%
  filter(vegetation %in% forests & myc_group != "ERC" & myc_group != "Other" & pft != "DG") %>%
  dplyr::select(pft, Temp, mat, Prec, vegetation, FT, myc_group, h.t, m.so, m.to, m.rt, m.lf, LmSo, LmTm, LmSm, LaSm, LmLa, RmTm) 

sub_s <- sub %>%
  group_by(pft) %>%
  summarise(n_distinct())
  
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


