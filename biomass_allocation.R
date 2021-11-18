##Tree biomass allocation differs by mycorrhizal association
#For submission to Ecology
#Ashley Lang and Fiona Jevon

#databases used in this analysis: 
### BAAD database 
#(Falster, Daniel S., Remko A. Duursma, Masae I. Ishihara, Diego R. Barneche, Richard G. FitzJohn, Angelica Vårhammar, Masahiro Aiba et al. "BAAD: a Biomass And Allometry Database for woody plants." Ecology 96, no. 5 (2015): 1445.)

### FungalRoot database 
#Soudzilovskaia, Nadejda A., Stijn Vaessen, Milagros Barcelo, Jinhong He, Saleh Rahimlou, Kessy Abarenkov, Mark C. Brundrett, Sofia IF Gomes, Vincent Merckx, and Leho Tedersoo. "FungalRoot: global online database of plant mycorrhizal associations." New Phytologist 227, no. 3 (2020): 955-966.


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

###Read in FungalRoot data----
FR_occurrences <- read_csv("FR_occurrences.csv")
FR_measurements <- read_csv("FR_measurements.csv") %>%
  pivot_wider(names_from = measurementType, values_from = measurementValue) %>% 
  rename(Myc_type='Mycorrhiza type')

#summarise FungalRoot database  
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

#Choosing the most common mycorrhizal type associated with each species of plant
fungal_root_sp <- fungal_root %>%
  group_by(order, family, genus, species, myc_group) %>%
  summarise(number= sum(n)) %>%
  slice_max(order_by = number, n = 1) %>% 
  arrange(genus, species) %>% 
  filter(myc_group != "ECM/AM" & myc_group != "Other" & myc_group != "ERC") %>% 
  group_by(order, family, genus, species) %>% 
  mutate(dupe = n()>1) %>% 
  filter(dupe==F)

## Read in BAAD database, remove groups with poor coverage: myc types other than AM and ECM, and decidious gymnosperms
baad <- baad.data::baad_data()
dict <- as.data.frame(baad$dictionary)
baad_df <- as.data.frame(baad$data) %>% 
  dplyr::select(studyName, latitude, longitude, species, vegetation, map, mat, pft, a.lf, h.t, d.bh, m.lf, m.st, m.so, m.rt, m.to, 	ma.ilf) %>%
  separate(species, c("genus", "species"), extra = "drop", fill = "right") %>%
  left_join(fungal_root_sp, by = c("genus", "species")) %>% 
  mutate(LmTm = m.lf/ m.to,
         RmTm = m.rt/m.to,
         SmTm= m.st/m.to,
         pft = as.factor(pft)) %>% 
   filter( myc_group != "ERC" & myc_group != "Other" & pft != "DG")

###Extract MAT/MAP from WorldClim----
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


###filtering----
#Now the full version of the data with baad_df, myc types, and climate:
full_df = baad_df %>% 
  left_join(df, by=c("latitude", "longitude")) %>%
  filter(RmTm>0, h.t >.5, myc_group != "ECM/AM", family!= "Ericaceae") %>% 
  unite(study_species, studyName, genus, species, sep="_", remove=F) %>% 
  mutate(log_ht = log(h.t),
         leaf_habit= case_when(pft=="EA" | pft== "EG" ~ "evergreen",
                               pft=="DA" ~ "deciduous"),
         evo_group= case_when(pft=="EA" ~ "angiosperm",
                              pft== "EG" ~ "gymnosperm",
                              pft=="DA" ~ "angiosperm") ) 

###Figure 1:----

AM_ECM=c("#C49F50", "#91BBA8")

sub <- full_df %>%
  group_by(Temp, Prec, leaf_habit, myc_group, longitude, latitude) %>%
  summarise(n = n())

clim_space <- ggplot(sub, aes(x= Temp, y=Prec))+
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

 #map
world <- map_data("world")

map=ggplot(data=world)+
  geom_polygon(aes(x=long, y=lat, group=group), fill="white", color="black", size=0.2)+
  coord_fixed(1.3)+
  theme(axis.text=element_blank(), axis.line = element_blank(), 
        axis.ticks=element_blank(), axis.title=element_blank(), 
        panel.background = element_rect(fill = "white"), 
        panel.border = element_blank(), legend.position="none")+
  geom_point(aes(x = longitude, y = latitude,color=myc_group), data = sub, size = 1)+
  scale_colour_manual(values=AM_ECM)

#Figure 1:
ggarrange(map, clim_space, nrow=1, ncol=2, labels=c("a", "b")) 


####models-----
#make clean dataset for models 
full_df_mod <- full_df %>%
  dplyr::select(RmTm, LmTm, SmTm, log_ht, leaf_habit, myc_group, Temp, Prec, study_species, family) %>%
  drop_na() %>% 
  separate(study_species, into=c("Study", "Genus", "Species"), sep="_", remove=F) %>% 
  unite(SppName, c(Genus, Species), sep="_")
#1429 observations


full_df_mod %>% group_by(myc_group) %>% summarise(fam_num = n_distinct(family))
#families by myc type : 28 AM, 6 ECM

full_df_mod %>% group_by(leaf_habit) %>% summarise(fam_num = n_distinct(family))
#families by leaf habit: 18 decidous, 16 evergreen

full_df_mod %>% group_by(leaf_habit, myc_group) %>% summarise(n = n())

#check to see if variables need to be transformed
hist(full_df_mod[full_df_mod$myc_group=="AM",]$RmTm)
hist(full_df_mod$LmTm)
hist(full_df_mod$SmTm)

##Root mass/total mass model----

#test the simplest model to get an idea of the coefficient for myc
R_mini_model <-lmer(RmTm ~ myc_group  + (1|study_species), data = full_df_mod)
summary(R_mini_model)

#next, test full model with all interaction terms
R_full_model <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(R_full_model)
summary(R_full_model)

#same but with plant family added
R_full_model_f <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + family + (1|study_species), data = full_df_mod)
vif(R_full_model_f)
summary(R_full_model_f)
AIC(R_full_model_f)
#tab_model(R_full_model_f, show.se = TRUE, show.ci = FALSE, digits = 3, digits.re = 3, show.std = "std2")

#remove insignificant interaction terms, testing for AIC improvement: log_ht*myc_group, Temp*myc_group, log_ht*Temp, Temp*leaf_habit (marginally significant)
R_reduced1 <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced1)
AIC(R_reduced1)
R_reduced2 <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced2)
AIC(R_reduced2)
R_reduced3 <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + log_ht*Temp + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced3)
AIC(R_reduced3)
R_reduced4 <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + Temp + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced4)
AIC(R_reduced4)
vif(R_reduced4)

AIC(R_reduced1, R_reduced2, R_reduced3, R_reduced4)
#use ggeffect to calculate the marginal effects of myc group on RM/TM across tree heights
R_E1 <- ggeffect(R_reduced4, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
R_E1$height=exp(R_E1$x)

##Leaf mass/total mass models

#first, test the simplest model to get an idea of the coefficient for myc
L_mini_model <-lmer(LmTm ~ myc_group  + (1|study_species), data = full_df_mod)
summary(L_mini_model)

#full model (all terms and interactions, including plant family)
L_full_model <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(L_full_model)
summary(L_full_model)
AIC(L_full_model)

#remove insignificant interaction terms, testing for AIC improvement
L_reduced1 <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + Temp*leaf_habit + Prec  + (1|study_species), data = full_df_mod)
summary(L_reduced1)
AIC(L_reduced1)

L_reduced2 <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + Prec  + (1|study_species), data = full_df_mod)
summary(L_reduced2)
AIC(L_reduced2)

L_reduced3 <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec  + (1|study_species), data = full_df_mod)
summary(L_reduced3)
AIC(L_reduced3)
vif(L_reduced3)

AIC(L_reduced1, L_reduced2, L_reduced3)
#use ggeffect to calculate the marginal effects of myc group acorss tree heights 
L_E1 <- ggeffect(L_reduced3, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
L_E1$height=exp(L_E1$x)

#Stem mass/Total mass full model
#first, test the simplest model to get an idea of the coefficient for myc
S_mini_model <-lmer(SmTm ~ myc_group  + (1|study_species), data = full_df_mod)
summary(S_mini_model)

S_full_model <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(S_full_model)
summary(S_full_model)
AIC(S_full_model)

#AKL 11/1 reduced model: remove myc_group*temp, logheight*temp, and leafhabit*Temp
#with plant family
#S_reduced_m2 <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + family +(1|study_species), data = full_df_mod)
#without plant family
S_reduced1 <- lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced1)
AIC(S_reduced1)

S_reduced2 <- lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced2)
AIC(S_reduced2)

S_reduced3 <- lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced3)
AIC(S_reduced3)

S_reduced4 <- lmer(SmTm ~ log_ht*leaf_habit + myc_group + Temp + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced4)
AIC(S_reduced4)

S_reduced5 <- lmer(SmTm ~ log_ht + leaf_habit + myc_group + Temp + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced5)
AIC(S_reduced5)

#confirm final reduced model has lowest AIC:
AIC(S_reduced1,S_reduced2, S_reduced3, S_reduced4, S_reduced5, S_full_model)


#use ggeffect to calculate the marginal effects of myc group and height 
S_E1 <- ggeffect(S_reduced5, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
S_E1$height=exp(S_E1$x)

# Best, reduced model output acorss all three tissue types
tab_model(L_reduced3, S_reduced5, R_reduced4, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)

#Figure 2:
#plot the marginal effects and the raw data for Rm/Tm

a <- ggplot() +
  geom_point(data = full_df, aes(x = h.t, y = RmTm, color = myc_group, shape=leaf_habit), alpha = .15, size=3) +
  geom_line(data = R_E1, aes(x = height, y = predicted, color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal Type", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17))+
    theme(legend.position = "none")+
  scale_x_continuous(trans='log2')+
  labs(x= "Tree height (m)", y="Root mass fraction" )+
  ylim(0,1)+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))

#plot the marginal effects and the raw data for Lm/Tm
b <- ggplot() +
  geom_point(data = full_df, aes(x = h.t, y = LmTm, color = myc_group, shape=leaf_habit), alpha = .15, size=3) +
  geom_line(data = L_E1, aes(x = height, y = predicted, color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal\nType", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17)) +
  scale_x_continuous(trans='log2')+
  theme(legend.position = c(.67, .75), 
        legend.title= element_text(hjust=0.5, size=10), 
        legend.text=element_text(size=9),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"), 
        legend.spacing.y = unit(1, "mm"))+
  labs(x= "Tree height (m)", y="Leaf mass fraction" )+
  ylim(0,1)+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))
  

#plot the marginal effects and the raw data for Sm/Tm
c <- ggplot() +
  geom_point(data = full_df, aes(x = h.t, y = SmTm, color = myc_group, shape=leaf_habit), alpha = .15, size=3) +
  geom_line(data = S_E1, aes(x = height, y = predicted, color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal Type", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17)) +
  theme(legend.position = "none")+
  scale_x_continuous(trans='log2', )+
  labs(x= "Tree height (m)", y="Stem mass fraction" )+
  ylim(0,1)+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))


R_E2 <- ggeffect(R_reduced4, terms = c("log_ht[0:1]", "myc_group"), type = "random")
L_E2 <- ggeffect(L_reduced3, terms = c("log_ht[0:1]", "myc_group"), type = "random")
S_E2 <- ggeffect(S_reduced5, terms = c("log_ht[0:1]", "myc_group"), type = "random")

L_E2=as.data.frame(L_E2) %>% 
  mutate(model="Leaf")
S_E2=as.data.frame(S_E2) %>% 
  mutate(model="Stem")

d_data=as.data.frame(R_E2) %>% 
  mutate(model="Root") %>% 
  rbind(L_E2) %>% 
  rbind(S_E2) %>% 
  filter(x=="0") %>% 
  group_by(group) %>% 
  mutate(ypos=case_when(model=="Leaf" ~ sum(predicted),
                        model=="Stem" ~ sum(predicted[model != "Leaf"]),
                        model=="Root" ~ predicted)) %>% 
  mutate(scaled_proportion=case_when(group=="AM" ~ ypos/0.9907910,
                                     group=="ECM" ~ ypos/0.9871502)) %>% 
  mutate(scaled_predicted=case_when(group=="AM" ~ predicted/0.9907910,
                                    group=="ECM" ~ predicted/0.9871502))

d_data$model=factor(d_data$model, levels=c("Leaf", "Stem", "Root"))
d=ggplot(data=d_data, aes(x=group, y=scaled_predicted, fill=model))+
  geom_bar(position="stack",
           stat="identity",
           colour="black")+
   scale_fill_manual(values=c("gray92", "gray60", "gray25"))+
  geom_errorbar(aes(ymin = scaled_proportion-std.error, ymax = scaled_proportion+std.error), width = 0.2, position="identity")+
  labs(x=" ", y= "Predicted proportion\nof total biomass")+
  theme(legend.title=element_blank())


ggarrange( b, c, a, d, labels=c("a", "b", "c" , "d"), nrow=2, ncol=2)
#ggsave("Figure_2.tiff")