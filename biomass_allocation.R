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

###climate and geography of dataset----

AM_ECM=c("#C49F50", "#91BBA8")

sub <- full_df %>%
  group_by(Temp, Prec, leaf_habit, myc_group, longitude, latitude) %>%
  summarise(n = n())

clim_space <- ggplot(sub, aes(x= Temp, y=Prec))+
  geom_point(aes(colour=myc_group,shape=leaf_habit, size=n), alpha=0.65)+
  scale_colour_manual(labels = c("AM", "ECM"),values=AM_ECM)+
  scale_shape_manual(labels = c("Deciduous", "Evergreen"), values=c(16,17))+
  labs(x=expression("Mean Annual Temperature ("*degree*C*")"), y="Mean Annual\nPrecipitation (mm)")+
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
        panel.border = element_blank(), legend.position="none", 
        plot.margin = unit(c(.01,.01,.01,.01), "lines"))+
  geom_point(aes(x = longitude, y = latitude,color=myc_group), data = sub, size = 1)+
  scale_colour_manual(values=AM_ECM)


####model prep-----
#make clean dataset for models 
full_df_mod <- full_df %>%
  dplyr::select(RmTm, LmTm, SmTm, m.to, log_ht, leaf_habit, myc_group, Temp, Prec, study_species, family, order) %>%
  drop_na() %>% 
  separate(study_species, into=c("Study", "Genus", "Species"), sep="_", remove=F) %>% 
  unite(SppName, c(Genus, Species), sep="_")
#1429 observations

full_df_mod %>% group_by(myc_group) %>% summarise(fam_num = n_distinct(SppName))
#species by myc type : 34 AM, 21 ECM

unique(full_df_mod$order)
#16 orders

full_df_mod %>% group_by(myc_group) %>% summarise(order_num = n_distinct(order))
#orders by myc type : 15 AM, 4 ECM

full_df_mod %>% group_by(leaf_habit) %>% summarise(order_num = n_distinct(order))
#orders by leaf habit: 18 decidous, 16 evergreen

full_df_mod %>% group_by(leaf_habit, myc_group) %>% summarise(n = n())
#group numbers

#check variable distributions
hist(full_df_mod$RmTm)
hist(full_df_mod[full_df_mod$myc_group=="AM",]$RmTm)
hist(full_df_mod$LmTm)
hist(full_df_mod$SmTm)
hist(full_df_mod$m.to) #needs transformation
full_df_mod$log_totbio = log(full_df_mod$m.to)
hist(full_df_mod$log_totbio)

#Total biomass model----
T_mini_model <- lmer(log_totbio ~ log_ht*myc_group + (1|study_species), data = full_df_mod)
summary(T_mini_model)

#test full model with all interaction terms
T_full_model <-lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(T_full_model)
summary(T_full_model)

#remove insignificant interaction term (myc_group*Temp), testing for AIC improvement.
T_reduced1 <- lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(T_reduced1)
#all remaining terms are significant

#check that the removal of interaction term also decreased AIC)
AIC(T_full_model, T_reduced1)

#get effects from final model 
T_E1 <- ggeffect(T_reduced1, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
T_E1$height=exp(T_E1$x)

#for supplement: 
#repeat process with plant order added to model
T_mini_model_o <- lmer(log_totbio ~ log_ht*myc_group + order + (1|study_species), data = full_df_mod)
summary(T_mini_model_o)
#no significant plant orders

#full model with plant order added
T_full_model_o <-lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
vif(T_full_model_o)
summary(T_full_model_o)
#note: no significant plant orders

#remove insignificant interaction (myc_group*Temp)
T_reduced1_o <- lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
summary(T_reduced1_o)

#remove insignificant interaction (log_ht:Temp)
T_reduced2_o <- lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + leaf_habit*Temp + Prec + order+ (1|study_species), data = full_df_mod)
summary(T_reduced2_o)

#remove insignificant interaction (leaf_habit:Temp)
T_reduced3_o <- lmer(log_totbio ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(T_reduced3_o)

#check AIC of various models with plant order
AIC(T_full_model_o, T_reduced1_o, T_reduced2_o, T_reduced3_o)
#T_reduced3_o is the best

##Root mass/total mass model----

#test the simplest model to get an idea of the coefficient for myc
R_mini_model <-lmer(RmTm ~ myc_group  + (1|study_species), data = full_df_mod)
summary(R_mini_model)

#next, test full model with all interaction terms
R_full_model <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(R_full_model)
summary(R_full_model)

#remove insignificant interaction terms:  Temp*myc_group
R_reduced1 <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced1)

#remove insignificant interaction term:  Temp*log_ht
R_reduced2 <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced2)

#remove insignificant interaction term: Temp*log_ht
R_reduced3 <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced3)

#remove insignificant interaction terms:  Temp*leaf_habit
R_reduced4 <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + Temp + Prec + (1|study_species), data = full_df_mod)
summary(R_reduced4)
vif(R_reduced4)

#Test various reduuded models for AIC improvement
AIC(R_full_model, R_reduced1, R_reduced2, R_reduced3, R_reduced4)
#R_reduced4 is the best model

R_E1 <- ggeffect(R_reduced4, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
R_E1$height=exp(R_E1$x)


#for supplement: repeat process with plant order added to model
R_full_model_o <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
vif(R_full_model_o)
summary(R_full_model_o)

#remove insignificant interaction terms:  Temp*leaf_habit
R_reduced1_o <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(R_reduced1_o)

#remove insignificant interaction term:  Temp*myc_group
R_reduced2_o <-lmer(RmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(R_reduced2_o)

#remove insignificant interaction term: Myc_group*log_ht
R_reduced3_o <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + log_ht*Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(R_reduced3_o)

#remove insignificant interaction terms:  Temp*leaf_habit
R_reduced4_o <-lmer(RmTm ~ log_ht*leaf_habit + myc_group + Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(R_reduced4_o)
vif(R_reduced4_o)

#Test various reducded models for AIC improvement
AIC(R_full_model_o, R_reduced1_o, R_reduced2_o, R_reduced3_o, R_reduced4_o)
#R_reduced4_o is the best model

##Leaf mass/total mass model----
#first, test the simplest model to get an idea of the coefficient for myc
L_mini_model <-lmer(LmTm ~ myc_group  + (1|study_species), data = full_df_mod)
summary(L_mini_model)

#full model (all terms and interactions)
L_full_model <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(L_full_model)
summary(L_full_model)

#remove insignificant interaction terms: Myc_group*Temp
L_reduced1  <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(L_reduced1)

#remove insignificant interaction terms: Leaf_habit*Temp
L_reduced2  <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Prec + (1|study_species), data = full_df_mod)
summary(L_reduced2)

#remove insignificant interaction terms: log_ht*Temp
L_reduced3  <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + (1|study_species), data = full_df_mod)
summary(L_reduced3)
vif(L_reduced3)

AIC(L_full_model, L_reduced1, L_reduced2, L_reduced3)
#L_reduced3 is the best model

#use ggeffect to calculate the marginal effects of myc group across tree heights 
L_E1 <- ggeffect(L_reduced3, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
L_E1$height=exp(L_E1$x)

#repeat process but including plant order
#full model (all terms and interactions, including plant order)
L_full_model_o <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
vif(L_full_model_o)
summary(L_full_model_o)

#remove insignificant interactions: height*temp
L_reduced1_o <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
summary(L_reduced1_o)

#remove insignificant interactions: height*temp
L_reduced2_o <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + Prec + order + (1|study_species), data = full_df_mod)
summary(L_reduced2_o)

#remove insignificant interactions: myc_group*temp
L_reduced3_o <-lmer(LmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(L_reduced3_o)
vif(L_reduced3_o)

AIC(L_full_model_o, L_reduced1_o, L_reduced2_o, L_reduced3_o)
#L_reduced3_o is best model

#Stem mass/total mass full model----
#first, test the simplest model to get an idea of the coefficient for myc
S_mini_model <-lmer(SmTm ~ myc_group  + (1|study_species), data = full_df_mod)
summary(S_mini_model)

#test full model with all variables and reasonable interactions
S_full_model <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
vif(S_full_model)
summary(S_full_model)

#remove insignificant interactions: height*myc_group
S_reduced1 <-lmer(SmTm ~ log_ht*leaf_habit +  Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced1)

#remove insignificant interactions: temp*myc_group
S_reduced2 <-lmer(SmTm ~ log_ht*leaf_habit + myc_group + log_ht*Temp + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced2)

#remove insignificant interactions: temp*height
S_reduced3 <-lmer(SmTm ~ log_ht*leaf_habit + myc_group + Temp*leaf_habit + Prec + (1|study_species), data = full_df_mod)
summary(S_reduced3)

#confirm final reduced model has lowest AIC:
AIC(S_full_model, S_reduced1, S_reduced2, S_reduced3)
#S_reduced3 is best

#use ggeffect to calculate the marginal effects of myc group and height 
S_E1 <- ggeffect(S_reduced3, terms = c("log_ht[-.7:3.5]", "myc_group"), type = "random")
S_E1$height=exp(S_E1$x)

#test full model with all variables and reasonable interactions, including plant order
S_full_model_o <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
summary(S_full_model_o)

#remove insignifncant interactions: myc_group*temp
S_reduced1_o <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Temp*leaf_habit + Prec + order + (1|study_species), data = full_df_mod)
summary(S_reduced1_o)

#remove insignifncant interactions: leaf_habit*temp
S_reduced2_o <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + log_ht*Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(S_reduced2_o)

#remove insignifncant interactions: log_ht*temp
S_reduced3_o <-lmer(SmTm ~ log_ht*leaf_habit + log_ht*myc_group + Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(S_reduced3_o)

#remove insignifncant interactions: log_ht*myc_group
S_reduced4_o <-lmer(SmTm ~ log_ht*leaf_habit + myc_group + Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(S_reduced4_o)

#remove insignifncant interactions: log_ht*leaf_habit
S_reduced5_o <-lmer(SmTm ~ log_ht + leaf_habit + myc_group + Temp + Prec + order + (1|study_species), data = full_df_mod)
summary(S_reduced5_o)

AIC(S_full_model_o, S_reduced1_o, S_reduced2_o, S_reduced3_o, S_reduced4_o, S_reduced5_o)
#S_reduced5_o is the best

#combining models into a table and creating Fig 2----
#Table 1: three tissue types
tab_model(L_reduced3, S_reduced3, R_reduced4,  show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)

#Table S1: total biomass model with and without plant order
tab_model(T_reduced1, T_reduced3_o, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)

#Table S2: best model output across all three tissue types with plant order 
tab_model(L_reduced3_o, S_reduced5_o, R_reduced4_o, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)

#check for changes in models that include order
tab_model(L_reduced3, L_reduced3_o, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)
tab_model(S_reduced3, S_reduced5_o, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)
tab_model(R_reduced4, R_reduced4_o, show.se = TRUE, show.ci = FALSE, show.std = "std2", digits = 3, digits.re = 3)

#Figure 2:
#plot the marginal effects and the raw data
full_df$log_totbio = log(full_df$m.to)

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

#plot the marginal effects and the raw data for total biomass: not for use in paper
e <- ggplot() +
  geom_point(data = full_df, aes(x = h.t, y = m.to, color = myc_group, shape=leaf_habit), alpha = .15, size=3) +
  geom_line(data = T_E1, aes(x = height, y = exp(predicted), color = group), size = 1.5) +
  scale_colour_manual(name="Mycorrhizal Type", values=AM_ECM)+
  scale_shape_manual(name="Leaf Habit", labels = c("Deciduous", "Evergreen"), values=c(16,17)) +
  theme(legend.position = "none")+
  scale_x_continuous(trans='log2' )+
  scale_y_continuous(trans='log10', breaks = c(0.01, 0.1, 1, 10, 100, 1000), labels=function(n){format(n, scientific = FALSE)})+
  labs(x= "Tree height (m)", y="Total biomass (kg)" )+
  guides(color = guide_legend(override.aes = list(alpha=1), order=1), shape = guide_legend(override.aes = list(alpha=1), order=2))

#best models: L_reduced3, S_reduced3, R_reduced4, T_reduced1
R_E2 <- ggeffect(R_reduced4, terms = c("log_ht[0:1]", "myc_group"), type = "random")
L_E2 <- ggeffect(L_reduced3, terms = c("log_ht[0:1]", "myc_group"), type = "random")
S_E2 <- ggeffect(S_reduced3, terms = c("log_ht[0:1]", "myc_group"), type = "random")

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
  mutate(scaled_proportion=case_when(group=="AM" ~ ypos/0.9778042,
                                     group=="ECM" ~ ypos/1.0061589)) %>% 
  mutate(scaled_predicted=case_when(group=="AM" ~ predicted/0.9778042,
                                    group=="ECM" ~ predicted/1.0061589))

d_data$model=factor(d_data$model, levels=c("Leaf", "Stem", "Root"))
d=ggplot(data=d_data, aes(x=group, y=scaled_predicted, fill=model))+
  geom_bar(position="stack",
           stat="identity",
           colour="black")+
   scale_fill_manual(values=c("gray92", "gray60", "gray25"))+
  geom_errorbar(aes(ymin = scaled_proportion-std.error, ymax = scaled_proportion+std.error), width = 0.2, position="identity")+
  labs(x=" ", y= "Predicted proportion\nof total biomass")+
  theme(legend.title=element_blank())
d
fig2 <- ggarrange( b, c, a, d, labels=c("a", "b", "c" , "d"), nrow=2, ncol=2)
ggsave("Figure_2.pdf", plot = fig2, width = 8 , height = 7, units = c("in"))

##phylogenetic tree for Fig 1----
AM_ECM=c("#C49F50", "#91BBA8")

spp_sum <- full_df_mod %>%
  group_by(SppName, myc_group, leaf_habit, order, family) %>%
  summarise(n = n_distinct())

#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
library(ape)

#create df with species of interest in appropriate format
sp_phy <- full_df_mod %>%
  group_by(SppName, family) %>%
  summarise(n = n_distinct()) %>%
  mutate(species = str_replace(SppName, "_", " ")) %>%
  separate(SppName, into = c("genus", "spp"), sep = "_") %>%
  dplyr::select(species,genus, family) %>%
  mutate(family = replace(family, family == "Aceraceae", "Sapindaceae"),
         family = replace(family, family == "Taxodiaceae", "Cupressaceae")) #fix old familiy names

#generate tree
tree.a <- phylo.maker(sp.list = sp_phy, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#reordering to color by myc type and shape by leaf habit
reorder_idx <- match(spp_sum$SppName, tree.a$scenario.3$tip.label) 
colors <- spp_sum %>%
  arrange(reorder_idx) %>%
  dplyr::select(myc_group)
col <- as.factor(colors$myc_group)
mycol <- c("#C49F50", "#91BBA8")[col]

shps <- spp_sum %>%
  arrange(reorder_idx) %>%
  dplyr::select(leaf_habit)
shps2 <- as.factor(shps$leaf_habit)
myshps <- c(16,17)[shps2]

#make figure with labelled phylogenetic tree of the tree species used in this analysis
plot(tree.a$scenario.3, tip.color = mycol, cex = .6, label.offset = 12, no.margin = TRUE)
tiplabels(pch = myshps, col = mycol, cex = .7, adj = 7)
x <- recordPlot()
phylogeny <- as_grob(x)

#Figure 1:
fig1ab <- ggarrange(map, clim_space, nrow = 2, labels = c("b", "c"))
fig1 = ggarrange(phylogeny, fig1ab, widths = c(.5, 1), ncol = 2,  labels = c("a", ""))
fig1
#ggsave("Figure_1.pdf", plot = fig1, width = 9, height = 7, units = c("in"))
