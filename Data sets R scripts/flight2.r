# Flight assay data analysis for flight column
if(!require(psych)){install.packages("psych")}
if(!require(multcompView)){install.packages("multcompView")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(FSA)){install.packages("FSA")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(phia)){install.packages("phia")}
if(!require(car)){install.packages("car")} ## allows Type I, II, III sum of squares
install.packages("gmodels") # a family of alternative funbction for summarizing data
install.packages("tidyverse")
install.packages("magrittr")
install.packages("dplyr") 
install.packages("psych")

library(psych)
library(magrittr)
library(dplyr)
library(tidyverse)
library(gmodels)
library(ggplot2)
library(car)

setwd('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/')
## Two folders with data: batch_21-07-15 and batch_21_07_30. Need to adjust read.csv and write.csv directories accordingly.
## Changed 'deficiency' variable to 'nuc' in all_flies.csv output from flight analysis.

flight<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_30/all_flies.csv', header = TRUE)
head(flight)
dim(flight)

## Remove all w1118 mtDNA samples 21-07_30 to make this set fully factorial
flight <- flight[which(flight$mito != "w1118"), ]
dim(flight)

flight$mito = as.factor(flight$mito)
flight$nuc = as.factor(flight$nuc)
flight$mitonuc = as.factor(paste(flight$mito, flight$nuc, sep="_")) 
flight$sex= as.factor(flight$sex)
flight$fly= as.factor(flight$fly)
#flight$day= as.factor(flight$day)

flightm=subset(flight, sex == 'm')
#flightm10=subset(flightm, day =='10')
#flightm24=subset(flightm, day =='24')
flightf=subset(flight, sex == 'f')
#flightf10=subset(flightf, day =='10')
#flightf24=subset(flightf, day =='24')


dim(flight)
head(flight)
dim(flightf)
dim(flightm)

## set margin inches "mai" to accommodate vertical x label names
par("mai" = c(1.5, 1.25, 0.5, 0.5))
par(mgp=c(3,1,0))
boxplot(Y~mito*nuc*sex,data = flight, las=2, xlab="", main="Flight Performance by mtDNA, Nuclear Genotype and Sex", ylab="Flight Performance (Landing Height)", col=c("red","pink", "yellow", "green", "blue"))

fitm <- aov(Y~mito*nuc,data = flightm)
summary(fitm)
Tukeym <- TukeyHSD(fitm)
Tukeym
#write to a table.  Not working
#Tukeym_table <- as.data.frame(Tukeym)
#write.csv(Tukeym_table ,file="Tukeym.csv")

## reverse order of mito & nuc
fitma <- aov(Y~nuc*mito,data = flightm)
summary(fitma)

fitf <- aov(Y~mito*nuc,data = flightf)
summary(fitf)
Tukeyf <- TukeyHSD(fitf)
Tukeyf
write.csv(Tukeyf)

## reverse order of mito & nuc
fitfa <- aov(Y~nuc*mito,data = flightf)
summary(fitma)

## linear models for Type III SSQ
fit2fmod3 <- lm(Y~mito * nuc, data = flightf, contrasts = list(mito = "contr.sum", nuc = "contr.sum"))
fit2fmod3_result <- Anova(fit2fmod3, type="III")
print(fit2fmod3_result)
## ERROR: Error in Anova_III_lm(mod, error, singular.ok = singular.ok, ...) : 
## there are aliased coefficients in the model

fit2mmod3 <- lm(Y~mito * nuc, data = flightm, contrasts = list(mito = "contr.sum", nuc = "contr.sum"))
fit2mmod3_result <- Anova(fit2mmod3, type="III")
print(fit2mmod3_result)

#model2 = lm(Y ~ mito*nuc*sex,data = subset(flight, mito != 'si1' & nuc != 'w1118'))
#Defs is not the same - check for what genotypes are missing is missing
#Defs=subset(flight, mito != 'w1118') #& nuc != 'si1'
#Defs=subset(flight, mito != 'si1' & nuc != 'w1118')
#summary(Defs)
#model2 = lm(Y ~ mito*nuc*sex,data = subset(flight, mito != 'si1' & nuc != 'w1118'))
#model3 = lm(Y ~ mito*nuc*sex,data = subset(flight, mito != 'si1' & mito != 'ore'))
#modelm = lm(Y ~ mito*nuc,data = subset(flightm, mito != 'si1' & mito != 'ore'))
#modelf = lm(Y ~ mito*nuc,data = subset(flightf, mito != 'si1' & mito != 'ore'))
modelm = lm(Y ~ mito*nuc,data = flightm)
modelf = lm(Y ~ mito*nuc,data = flightf)

# Anova(modelm, type = "II")
summary(modelm)

# Anova(modelf, type = "II")
summary(modelf)


# Male interaction plot
interaction.plot(x.factor     = flightm$nuc,
                 trace.factor = flightm$mito,
                 response     = flightm$Y,
                 fun = mean,
                 type="b",
                 col=c("black","red","green", "blue", "pink"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15, 16, 18),     ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o",
                 las=0, xlab="Nuclear Genotype", main="Male Flight Performance by Mitonuclear Genotype", ylab="Flight Performance (Landing Height)")


interaction.plot(x.factor     = flightf$nuc,
                 trace.factor = flightf$mito,
                 response     = flightf$Y,
                 fun = mean,
                 type="b",
                 col=c("black","red","green", "blue", "pink"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15, 16, 18),     ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o",
                 las=0, xlab="Nuclear Genotype", main="Female Flight Performance by Mitonuclear Genotype", ylab="Flight Performance (Landing Height)")

#glimpse(flightf) not a function

## ggplot box plots from https://sebastiansauer.github.io/vis_interaction_effects/
ggplot(flight) +
  aes(x = mitonuc, y = Y) +
  geom_boxplot() +
  facet_wrap(~sex)

## ggplot interactions from https://ggplot2tutor.com/tutorials/interaction_plot
## Male Interaction ggplot
flightm %>% 
  group_by(mito, nuc) %>% 
  #summarise(mnY = mean(Y)) -> flightm2  script from web site above
  summarise(mnY = ci(Y)[1],  # Joaquin's BETTER version
            mnY_low = ci(Y)[2],
            mnY_high = ci(Y)[3]) -> flightm2

flightm2

flightm2 %>% 
  ggplot() +
  aes(x = nuc, y = mnY, ymin = mnY_low, ymax = mnY_high, color = mito) +
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  #theme_classic() +
  labs(
    title = "Male Flight Performance by MitoNuclear Genotype",
    #subtitle = paste0("Nuclear genotype has largest effect\n",
    #                  "Mito x Nuclear Interactions Present"),
    x = "Nuclear Genotype",
    y = "Flight Performance (Landing Height)"
  ) 

# Female interaction ggplot
## ggplot interactions from https://ggplot2tutor.com/tutorials/interaction_plot

flightf %>% 
  group_by(mito, nuc) %>% 
  #summarise(mnY = mean(Y)) -> flightf2  script from web site above
  summarise(mnY = ci(Y)[1],  # Joaquin's BETTER version
            mnY_low = ci(Y)[2],
            mnY_high = ci(Y)[3]) -> flightf2

flightf2

flightf2 %>% 
  ggplot() +
  aes(x = nuc, y = mnY, ymin = mnY_low, ymax = mnY_high, color = mito) +
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  #theme_classic() +
  labs(
    title = "Female Flight Performance by MitoNuclear Genotype",
    #subtitle = paste0("Nuclear genotype has largest effect\n",
    #                  "Mito x Nuclear Interactions Present"),
    x = "Nuclear Genotype",
    y = "Flight Performance (Landing Height)"
  ) 

## Leah's CV bootstrap code email 2025-06-20
## Edited from Climbing analyses - no 'day' variable, 'geno' in Climbing = 'nuc' in Flight

cv = function(x) sd(x) / mean(x)

boot_cv = function(x, n=1000){
  replicate(n,cv(sample(x, replace = TRUE)))
}

mean_flight = flight %>%
  group_by(mito,nuc,sex) %>%
  summarise(Y = mean(Y))

mean_flight_1491set.df <- data.frame(mean_flight)
write.csv(mean_flight_1491set.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_30/mean_flight1491set.csv')

boot_flight = mean_flight %>% 
  group_by(nuc,sex) %>%
  summarise(boot = list(boot_cv(Y))) %>%
  tidyr::unnest(boot)

cv_flight = mean_flight %>%
  group_by(nuc,sex) %>%
  summarise(cv = cv(Y)) 

cv_flightm = cv_flight %>% filter(sex=="m") %>% ungroup()
cv_flightf = cv_flight %>% filter(sex=="f")%>% ungroup()

cv_by_sex = cv_flightm %>% rename(cv_m = cv) %>% select(-sex) %>%
  left_join(cv_flightf %>% rename(cv_f = cv) %>% select(-sex), join_by(nuc))
write.csv(cv_by_sex, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_30/cv_by_sex_1491set.csv')


ggplot(cv_by_sex, aes(x = cv_f, y = cv_m, color = nuc))+
  geom_point(size=5) +
  geom_abline(intercept = 0, Y = 1) +
  xlim(0,.3) +
  ylim(0,.3) +
  theme_bw() +
  labs(color="Nuclear Genotype")

pbootflight = ggplot(boot_flight, aes(y = boot, x = nuc, fill = sex, color=sex)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.3),size=0.3) +
  geom_boxplot(color="black", alpha=0.5, outliers = FALSE)  + 
#  facet_wrap(~day) +
  labs(title = "Flight Performance (Landing Height)",
       x = "Genotype, Sex",
       y = "Bootstrapped CV") +
  theme_minimal()
pbootflight

dim(flight)

## OBSOLETE CODE FOR CALCULATING Coefficient of Variation 
library(psych)
## MALE Code
flightmmeans <- describeBy(Y ~ mitonuc, data = flightm) ## Actually need means separated by mito to get mito-CV in each condition
flightmmeans
flightmmeans2 <- do.call("rbind", flightmmeans)
flightmmeans2
flightmmeans.df <- data.frame(flightmmeans2)
flightmmeans.df$cv = flightmmeans.df$sd / flightmmeans.df$mean

## write to a file that allows taking means of mtDNAs across nuc and day
write.csv(flightmmeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_15/flightmmeans.csv')
## Edit this by hand in Excel = malemeans to have a nuc column added for CV calculation: CV across mtDNAs within a nuc = a deficiency
## READ this hand-edited file to allow mito means for mito-CV calculation 
flymmeans<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_15/flightmmeans.csv')
flymmitomeans <- describeBy(mean ~ nuc, data = flymmeans)
flymmitomeans2 <- do.call("rbind", flymmitomeans)
flymmitomeans.df <- data.frame(flymmitomeans2)
flymmitomeans.df$cv = flymmitomeans2$sd / flymmitomeans2$mean
flymmitomeans.df

## FEMALE Code
flightfmeans <- describeBy(Y ~ mitonuc, data = flightf) ## Actually need means separated by mito to get mito-CV in each condition
flightfmeans
flightfmeans2 <- do.call("rbind", flightfmeans)
flightfmeans2
flightfmeans.df <- data.frame(flightfmeans2)
flightfmeans.df$cv = flightfmeans.df$sd / flightfmeans.df$mean

## write to a file that allows taking means of mtDNAs across nuc and day
write.csv(flightfmeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_15/flightfmeans.csv')
## Edit this by hand in Excel: add a column that puts the nuc term from mitonuc into its own column so means and s.d. by Df can be calculated then cv column calculated and added
## READ this hand-edited file (use the same name bnut with new column added) to allow mito means for mito-CV calculation 
flyfmeans<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_15/flightfmeans.csv')
flyfmitomeans <- describeBy(mean ~ nuc, data = flyfmeans)
flyfmitomeans2 <- do.call("rbind", flyfmitomeans)
flyfmitomeans.df <- data.frame(flyfmitomeans2)
flyfmitomeans.df$cv = flyfmitomeans2$sd / flyfmitomeans2$mean
flyfmitomeans.df

write.csv(flymmitomeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_15/flymaleCVs.csv' )
write.csv(flyfmitomeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Flight/batch_21_07_15/flyfemaleCVs.csv' )

## PLOT Female x Male CVs 
ggplot() + geom_point(size=5, aes(x=flyfmitomeans.df$cv, y=flymmitomeans.df$cv)) +
  labs(
    title = "Female vs. Male mtDNA Coefficient of Varitation for Flight",
    #subtitle = paste0("4 mtDNAs in Nuclear Dfs 1491, 2959, 7837, w1118"),
    subtitle = paste0("3 mtDNAs in Nuclear Dfs 7744, 8469, w1118"),
    x = "Female Mitonuclear-Age CV",
    y = "Male Mitonuclear-Age CV)") +
  xlim(0,0.2) + ylim(0,0.2)

## trying to facet wrap by variable 'nuc'
flightf %>% 
  group_by(mito, nuc) %>% 
  #summarise(mnY = mean(Y)) -> flightf2 Old code
  summarise(mnY = ci(Y)[1],  # Joaquin's BETTER version
            mnY_low = ci(Y)[2],
            mnY_high = ci(Y)[3]) -> flightf2Df


flightf2Df %>% 
  ggplot() +
  aes(x = nuc, y = mnY, ymin = mnY_low, ymax = mnY_high, color = mito) +  ## Joaquin's BETTER version
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  facet_wrap(~nuc) +
  labs(
    title = "Male Flight Performance MitoNuclear Genotype",
    subtitle = paste0("Nuclear genotype has largest effect\n",
                      "Mito x Nuclear Interactions Present"),
    x = "Nuclear  Genotype",
    y = "Flight Performance (Landing Height)"
  ) 

## trying to facet wrap by variable 'sex'
flightf %>% 
  group_by(mito, nuc) %>% 
  #summarise(mnY = mean(Y)) -> flightf2 Old code
  summarise(mnY = ci(Y)[1],  # Joaquin's BETTER version
            mnY_low = ci(Y)[2],
            mnY_high = ci(Y)[3]) -> flightf2Df


flightf2Df %>% 
  ggplot() +
  aes(x = nuc, y = mnY, ymin = mnY_low, ymax = mnY_high, color = mito) +  ## Joaquin's BETTER version
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  facet_wrap(~nuc) +
  labs(
    title = "Female Flight Performance MitoNuclear Genotype",
    subtitle = paste0("Nuclear genotype has largest effect\n",
                      "Mito x Nuclear Interactions Present"),
    x = "Nuclear  Genotype",
    y = "Flight Performance (Landing Height)"
  ) 



## Joaquin code for distributions in interaction plot, by facets 8/4/2021
install.packages("hrbrthemes")
install.packages("ggdist")
library(hrbrthemes)
library(ggdist)

#flightm %>%
flight %>%
  ggplot(aes(
    x=nuc,
    y=Y,
    fill = mito # pold code had 'as.factor(mito)'
  )) +
  #  ggdist::stat_halfeye(
  #    adjust = 5, # large number smooths things, small number allows modality
  #    width = .6, 
  #    .width = 0, 
  #    justification = -.3, 
  #    point_colour = NA) + 
  geom_point(
    aes(group = mito),
    size = 0.9,
    alpha = .1,
    #    position = position_jitter( - this jitter puts all outliers in central position
    position = position_dodge(0.7)  # position_dodge puts outliers on each box
    #      seed = 1, width = .1)
  )  +
  geom_boxplot(
    #   geom_violin(
    size = 0.2,
    width = .7, 
    outlier.shape = NA
  ) + 
  theme_bw() + 
  #theme(legend.position = "none") +
  #coord_flip() +
  ylab("Flight Performance (Landing Height)") +
  xlab("Nuclear Genotype") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~sex) +  # joint sex x time effect in facests  see top of code block
  labs(
    title = "Female and Male Flight Performance by MitoNuclear Genotype",
    subtitle = paste0("Females show more mtDNA effects\n",
                      "Males show less mtDNA effects or Mito x Nuclear Interaction"),
    x = "Nuclear  Genotype",
    y = "Flight Performance (Landing Height)"
  )


## lines below error : data must be in data frame, not and S3 object.  Did this pass 'aes() to the 'data' argument?
ggplot(aes(flightm, Y)) +
  geom_line(size = 1.2, aes(group = nuc, color = nuc)) +
  geom_point(size = 2.6, aes(color = nuc), shape = 15)


