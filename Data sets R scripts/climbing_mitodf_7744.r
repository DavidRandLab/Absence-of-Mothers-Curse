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

setwd('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing')

## NOTE: analyses for two different results folders: 21_07_30 and batch8906_9061.  NOTE different version in 21_07_30 for ROIs
## climb1 = results from analysis with 1 ROI in freeclimber,that covered bothe Day_10 and Day_22 images
## climb2 = results from analysis with 2 ROIs, a separate  ROIfor Day_10 vs Day_22
## climb1<-read.csv('/Users/davidrand/Dropbox (Brown)/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/2021_07_30_Day_10+22one_results.csv', header = TRUE)

## USE 2 ROI data set - more  reliable - see comparisons of using one ROI vs. 2 separate ROIs
## NOTE: change working directories, and write.csv paths to match different raw csv, sample batch, Deficiency names and dates.
## Read data for batch8906_9061 batch of Dfs - NOTE: Days = 7 and 24
## climb2<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/results_8906_9061_Day_7&24.csv', header = TRUE)

## USE 2 ROI data set - more  reliable - see comparisons of using one ROI vs. 2 separate ROIs
## NOTE: change working directories, and write.csv paths to match different raw csv, sample batch, Deficiency names and dates.
## Read data for 2021_07_30 batch of Dfs
## climb2<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/2021_07_30_Day_10+22sep_results.csv', header = TRUE)

## NOTE: This data sets adds 2021_07_30_1vial and _3vial data stets for addtional climbs of low yield genotypes missing from 2021_07_30_results/2021_07_30_Day_10+22sep_results.csv
## climb2<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/2021_07_30_Day_10+22sep_w1vial_3vial_results.csv', header = TRUE)

## Read data for 2021_07_15 batch of Dfs 7744 8469 w1118
## NOTE MISSING: Day 24 genotypes siI_7744_f and _m, ore_w1118_f and _m, siI_w1118_f and _m
climb2<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_15_results/21_07_15_results.csv', header = TRUE)

head(climb2)
dim(climb2)

## Remove vial_IDs with "_all"
climb2 <- climb2[ grep("_all", climb2$vial_ID, invert = TRUE) , ]
dim(climb2)

## Remove mtDNAs ore and si1 to make a 3 mito X 3 nuc factorial 21-07_15
climb2 <- climb2[which(climb2$mito != "ore"), ]
climb2 <- climb2[which(climb2$mito != "si1"), ]
dim(climb2)


## Analysis with 2 separate ROIs
climb2$mito = as.factor(climb2$mito)
climb2$geno = as.factor(climb2$geno)
climb2$sex= as.factor(climb2$sex)
climb2$day= as.factor(climb2$day)
climb2$rep= as.factor(climb2$rep)
climb2$mitogeno = as.factor(paste(climb2$mito, climb2$geno, sep="_")) 
climb2$mitonucsex = as.factor(paste(climb2$mito, climb2$geno, climb2$sex, sep="_"))
climb2$mitonucsexday = as.factor(paste(climb2$mito, climb2$geno, climb2$sex, climb2$day, sep="_"))


climb2m=subset(climb2, sex == 'm')
climb2m10=subset(climb2m, day =='10')
climb2m24=subset(climb2m, day =='24')
## climb2m22=subset(climb2m, day =='22')
climb2f=subset(climb2, sex == 'f')
climb2f10=subset(climb2f, day =='10')
climb2f24=subset(climb2f, day =='24')
## climb2f22=subset(climb2f, day =='22')

dim(climb2)
dim(climb2m)
dim(climb2f)
dim(climb2m10)
dim(climb2m24)
dim(climb2f10)
dim(climb2f24)

dim(climb2)

## set margin inches "mai" to accommodate vertical x label names
par("mai" = c(1.5, 1.25, 0.5, 0.5))
par(mgp=c(3,1,0))
boxplot(slope~sex*day*mito*geno, data = climb2, las=2, xlab="", main="Climbing Speed by Sex (Purple, Green), Age, Mito Nuclear genotype", ylab="Climbing Speed (slope)", col=c("purple","green"))
## boxplot(slope~mito*geno*day*sex, data = climb2, las=2, xlab="", main="Climbing Speed by Sex (Purple, Green), Age, Mito Nuclear genotype", ylab="Climbing Speed (slope)", col=c("purple","green"))

## Fit ANOVA models
## fit1m <- aov(slope~mito*geno*day,data = climb1m)
## summary(fit1m)

fit2f <- aov(slope~mito*geno*day,data = climb2f)
summary(fit2f)

fit2fa <- aov(slope~geno*mito*day,data = climb2f)
summary(fit2fa)

## fit2flm <- lm(slope~mito*geno*day,data = climb2f)
## summary(fit2flm)fit2m <- aov(slope~mito*geno*day,data = climb2m) ## climb2m uses csv file based on two (2) ROI in film analyses for males (m).
fit2m <- aov(slope~mito*geno*day,data = climb2m)
summary(fit2m)

fit2ma <- aov(slope~geno*mito*day,data = climb2m) ## reverse order of mito & geno to compare order of aov type 1 SSQ
summary(fit2ma)

table(climb2f$mito, climb2f$geno, climb2f$day)
table(climb2m$mito, climb2m$geno, climb2m$day)

## fit2mlm <- lm(slope~mito*geno*day,data = climb2m)
## summary(fit2mlm)

## Type III sum of squares models (see https://rcompanion.org/rcompanion/d_04.html)
## options(contrasts = c("contr.sum", "contr.poly"))
## model.3 = lm(Y ~ A + B + A:B,
##              contrasts=list(A="contr.sum", B="contr.sum"))
## Anova(model.3, type="III")

fit2fmod3 <- lm(slope~mito * geno * day, data = climb2f, contrasts = list(mito = "contr.sum", geno = "contr.sum", day = "contr.sum"))
fit2fmod3_result <- Anova(fit2fmod3, type="III")
print(fit2fmod3_result)

fit2mmod3 <- lm(slope~mito * geno * day, data = climb2m, contrasts = list(mito = "contr.sum", geno = "contr.sum", day = "contr.sum"))
fit2mmod3_result <- Anova(fit2mmod3, type="III")
print(fit2mmod3_result)

## Interaction plots by sex and day (age)
## ggplot interactions from https://ggplot2tutor.com/tutorials/interaction_plot

## MALE Df background Facet wrap by 'day'
## climb2mtrim<-subset(climb2m, slope <5) ## TRIMMING outlier slope values  - does NOT fix w1118;7837 outlier
climb2m %>% 
  group_by(mito,geno,day) %>% 
  #summarise(mnslope = mean(slope)) -> climbm2 Old code
  summarise(mnslope = ci(slope)[1],  # Joaquin's BETTER version
            mnslope_low = ci(slope)[2],
            mnslope_high = ci(slope)[3]) -> climb2m2day


climb2m2day %>% 
  ggplot() +
  aes(x = geno, y = mnslope, ymin = mnslope_low, ymax = mnslope_high, color = mito) +  ## Joaquin's BETTER version
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  facet_wrap(~day) +
  #scale_fill_discrete(breaks=c('w1118', 'Zim', 'siI', 'sm21',  'yak')) +
  #theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  ylim(0.5,6.0) +
  labs(
    title = "Male Climbing of Mito-Deficiency Genotypes at Ages 10 or 24 Days",
#    subtitle = paste0("Significant Mito x Nuclear Interaction"),
    x = "Deficiency Genotype",
    y = "Climbing Speed (cm/second)"
  ) 

## FEMALE Cy OR Df background Facet wrap by 'day'
## climb2ftrim<-subset(climb2f, slope <5)
climb2f %>% 
  group_by(mito, geno, day) %>% 
  #summarise(mnslope = mean(slope)) -> climbm2 Old code
  summarise(mnslope = ci(slope)[1],  # Joaquin's BETTER version
            mnslope_low = ci(slope)[2],
            mnslope_high = ci(slope)[3]) -> climb2f2day


climb2f2day %>% 
  ggplot() +
  aes(x = geno, y = mnslope, ymin = mnslope_low, ymax = mnslope_high, color = mito) +  ## Joaquin's BETTER version
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  facet_wrap(~day) +
  #scale_fill_discrete(breaks=c('w1118', 'Zim', 'siI', 'sm21',  'yak')) +
  #theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  ylim(0.5,6.0) +
  #theme_classic() +
  labs(
    title = "Female Climbing of Mito-Deficiency Genotypes at Ages 10 or 24 Days",
#    subtitle = paste0("Significant Mito x Nuclear Interactions"),
    x = "Deficiency Genotype",
    y = "Climbing Speed (cm/sec)"
  ) 


## Leah's CV bootstrap code email 2025-06-20

cv = function(x) sd(x) / mean(x)

boot_cv = function(x, n=1000){
  replicate(n,cv(sample(x, replace = TRUE)))
}

mean_climb = climb2 %>%
  group_by(mito,geno,sex,day) %>%
  summarise(slope = mean(slope))

mean_climb_7744set_balanced.df <- data.frame(mean_climb)
write.csv(mean_climb_7744set_balanced.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_15_results/mean_climb_7744set_balanced.csv')

boot_climb = mean_climb %>% 
  group_by(geno,sex,day) %>%
  summarise(boot = list(boot_cv(slope))) %>%
  tidyr::unnest(boot)

cv_climb = mean_climb %>%
  group_by(geno,sex,day) %>%
  summarise(cv = cv(slope)) 

cv_climbm = cv_climb %>% filter(sex=="m") %>% ungroup()
cv_climbf = cv_climb %>% filter(sex=="f")%>% ungroup()

cv_by_sex = cv_climbm %>% rename(cv_m = cv) %>% select(-sex) %>%
  left_join(cv_climbf %>% rename(cv_f = cv) %>% select(-sex), join_by(geno,day))
write.csv(cv_by_sex, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_15_results/cv_by_sex_7744_balanced.csv')


ggplot(cv_by_sex, aes(x = cv_f, y = cv_m, color = geno))+
  geom_point(size=5) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(0,.3) +
  ylim(0,.3) +
  theme_bw() +
  labs(color="Nuclear Genotype")

pbootclimb = ggplot(boot_climb, aes(y = boot, x = geno, fill = sex, color=sex)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.3),size=0.3) +
  geom_boxplot(color="black", alpha=0.5, outliers = FALSE)  + 
  facet_wrap(~day) +
  labs(title = "Climbing",
       x = "Genotype,Environment",
       y = "Bootstrapped CV") +
  theme_minimal()
pbootclimb

dim(climb2)


climb2m=subset(climb2, sex == 'm')
climb2m10=subset(climb2m, day =='10')
climb2m22=subset(climb2m, day =='22')
## climb2m7=subset(climb2m, day =='7')
## climb2m24=subset(climb2m, day =='24')

climb2f=subset(climb2, sex == 'f')
climb2f10=subset(climb2f, day =='10')
climb2f22=subset(climb2f, day =='22')
## climb2f7=subset(climb2f, day =='7')
## climb2f24=subset(climb2f, day =='24')

dim(climb2)
dim(climb2m)
dim(climb2f)
dim(climb2m10)
dim(climb2m22)
dim(climb2f10)
dim(climb2f22)
## dim(climb2m7)
## dim(climb2m24)
## dim(climb2f7)
## dim(climb2f24)


##
## OBSOLETE CODE FOR CV calculations - requires custom file editing
## Means by  mito-geno-day
## describeBy(mydata, group,...)
## estatistica <- describeBy(pag,list(pag$Jogo)) cannot use write.table orwrite.csv
## put all list elements into one data frame with do.call() and rbind() and then write to file. This will make data frame where group names will be added before original variable names.
## estatistica2<-do.call("rbind",estatistica)
## estatistica2

male2means <- describeBy(slope ~ mitonucsexday, data = climb2m) ## Actually need means separated by mito to get mito-CV in each condition
male2means
male2means2 <- do.call("rbind", male2means)
male2mean.df <- data.frame(male2means2) ## NOTE 'mean' vs. 'means': male2mean.df is local to RAM; male2means.csv gets written to hard drive
male2mean.df$cv = male2mean.df$sd / male2mean.df$mean

##  descriptive statistics for the clin.trial data, broken down separately by therapy type. 
##  The command I would use here is:
##    describeBy( x=clin.trial, group=clin.trial$therapy )

## write to a file that allows taking means of mtDNAs across nuc and day
write.csv(male2mean.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/male2means.csv')
## write.csv(male2mean.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/male2means.csv')
## Edit this by hand in Excel = male2means to have a NU_DAY column added for CV calculation.
## READ this hand-edited file to allow mito means for mito-CV calculation 
male2means<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/male2means.csv')
## male2means<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/male2means.csv')
male2mitomeans <- describeBy(mean ~ nucday, data = male2means)
male2mitomeans2 <- do.call("rbind", male2mitomeans)
male2mitomeans.df <- data.frame(male2mitomeans2)
male2mitomeans.df$cv = male2mitomeans2$sd / male2mitomeans2$mean
male2mitomeans.df

## This calculates CV across mitonucsexday
female2means <- describeBy(slope ~ mitonucsexday, data = climb2f)
female2means
female2means2 <- do.call("rbind", female2means)
female2mean.df <- data.frame(female2means2)
female2mean.df$cv = female2mean.df$sd / female2mean.df$mean

## write to a file that allows taking means of mtDNAs across nuc and day
## write.csv(female2mean.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/female2means.csv')
write.csv(female2mean.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/female2means.csv')

## Edit this by hand in Excel = femalemeans to have a NUC_DAY column added for CV calculation.
## READ this hand-edited file to allow mito means for mito-CV calculation 
## female2means<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/female2means.csv')
female2means<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/female2means.csv')
female2mitomeans <- describeBy(mean ~ nucday, data = female2means)
female2mitomeans2 <- do.call("rbind", female2mitomeans)
female2mitomeans.df <- data.frame(female2mitomeans2)
female2mitomeans.df$cv = female2mitomeans2$sd / female2mitomeans2$mean
female2mitomeans.df

## write.csv(male2mitomeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/male2CVs.csv')
## write.csv(female2mitomeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/2021_07_30_results/female2CVs.csv')
write.csv(male2mitomeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/male2CVs .csv' )
write.csv(female2mitomeans.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mito Deficiency Screen/Mito Deficiency Climbing/batch8906_9061/female2CVs .csv' )

ggplot() + geom_point(size=5, aes(x=female2mitomeans.df$cv, y=male2mitomeans.df$cv)) +
  labs(
    title = "Female vs. Male mtDNA Coefficient of Varitation for Climbing Speed",
    subtitle = paste0("4 mtDNAs in Nuclear Dfs 1491,4959,7837,w1118 at Ages 10 and 22 Days"),
    x = "Female Mitonuclear-Age CV",
    y = "Male Mitonuclear-Age CV)") +
  xlim(0,0.5) + ylim(0,0.5)


#NOTE: 4/22/22025 different quantitative results with input data as ...one.csv vs. ... sep.csv: comparing how different ROIs interpret climbing movies
# Good news is that points are all below the 1:1 line rejecting Mothers Curese, but pointare in different locations.

## Compare with batch_8906_9061 ggplot:
ggplot() + geom_point(aes(x=femalemitomeans.df$cv, y=malemitomeans.df$cv)) +
  ggtitle('Female vs. Male Genetic CV for Climbing Speed Across 15 Mitonuclear Genotypes at Ages 7 and 24 Days') +
  xlab('Female Mitonuclear-Age CV') +
  ylab('Male Mitonuclear-Age CV') +
  xlim(0,0.5) + ylim(0,0.5)

# Check the correlation between the two version of ROI definition

## climb1all <- climb1[ grep("_all", climb1$vial_ID) , ]
climb2all <- climb2[ grep("_all", climb2$vial_ID) , ]
## dim(climb1all)
dim(climb2all)

## climb1 <- climb1[ grep("_all", climb1$vial_ID, invert = TRUE) , ]
climb2 <- climb2[ grep("_all", climb2$vial_ID, invert = TRUE) , ]
## dim(climb1)
dim(climb2)

head(climb1$slope)

reps1v2 <- data.frame(slope1 = climb1$slope,  
                      slope2 = climb2$slope)

## ggplot version
#load necessary libraries
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)

#create plot using individual 'reps' points, with regression line, regression equation, and R-squared
ggplot(data=reps1v2, aes(x=climb1$slope, y=climb2$slope)) +
  geom_smooth(method="lm") +
  geom_point() +
  xlim(0,9) +ylim(0,9)
stat_regline_equation(label.x=3, label.y=9) +
  stat_cor(aes(label=..rr.label..), label.x=3, label.y=8)


plot(slope1 ~ slope2,   # jitter offsets points so you can see them all
     data=reps1v2,  
     pch = 16,                 # shape of points
     cex = 1.0,                # size of points
     xlab="Climbing Slope rep values with One ROI",
     ylab="Climbing Slope rep values with Separate ROI")


## create a data frame with just the 'all' values for each genotype
all1v2 <- data.frame(slope1 = climb1all$slope,  
                     slope2 = climb2all$slope)

#create plot using only the 'all' rows  with regression line, regression equation, and R-squared
ggplot(data=all1v2, aes(x=slope1, y=slope2)) +
  geom_smooth(method="lm") +
  geom_point() +
  xlim(0,3) +ylim(0,3) +
  stat_regline_equation(label.x=0.5, label.y=2.5) +
  stat_cor(aes(label=..rr.label..), label.x=0.5, label.y=2)


plot(slope1 ~ slope2,   # jitter offsets points so you can see them all
     data=all1v2,  
     pch = 16,                 # shape of points
     cex = 1.0,                # size of points
     xlab="Climbing Slope All with One ROI",
     ylab="Climbing Slope All with Separate ROI")


## Analysis using the 1 ROI approach
climb1$mito = as.factor(climb1$mito)
climb1$geno = as.factor(climb1$geno)
climb1$sex= as.factor(climb1$sex)
climb1$day= as.factor(climb1$day)
climb1$rep= as.factor(climb1$rep)
climb1$mitogeno = as.factor(paste(climb1$mito, climb1$geno, sep="_")) 
climb1$mitonucsex = as.factor(paste(climb1$mito, climb1$geno, climb1$sex, sep="_"))
climb1$mitonucsexday = as.factor(paste(climb1$mito, climb1$geno, climb1$sex, climb1$day, sep="_"))

climb1m=subset(climb1, sex == 'm')
climb1m10=subset(climb1m, day =='10')
climb1m22=subset(climb1m, day =='22')
climb1f=subset(climb1, sex == 'f')
climb1f10=subset(climb1f, day =='10')
climb1f22=subset(climb1f, day =='22')

dim(climb1)
dim(climb1m)
dim(climb1f)
dim(climb1m10)
dim(climb1m22)
dim(climb1f10)
dim(climb1f22)

## set margin inches "mai" to accommodate vertical x label names
par("mai" = c(3, 1.25, 0.5, 0.5))
par(mgp=c(3,1,0))
boxplot(slope~sex*day*mito*geno, data = climb1, las=2, xlab="", main="Climbing Speed by Sex (Purple, Green), Age, Mito Nuclear genotype", ylab="Climbing Speed (slope)", col=c("purple","green"))


## Analysis with 2 separate ROIs
climb2$mito = as.factor(climb2$mito)
climb2$geno = as.factor(climb2$geno)
climb2$sex= as.factor(climb2$sex)
climb2$day= as.factor(climb2$day)
climb2$rep= as.factor(climb2$rep)
climb2$mitogeno = as.factor(paste(climb2$mito, climb2$geno, sep="_")) 
climb2$mitonucsex = as.factor(paste(climb2$mito, climb2$geno, climb2$sex, sep="_"))
climb2$mitonucsexday = as.factor(paste(climb2$mito, climb2$geno, climb2$sex, climb2$day, sep="_"))

climb2m=subset(climb2, sex == 'm')
climb2m10=subset(climb2m, day =='10')
climb2m22=subset(climb2m, day =='22')
climb2f=subset(climb2, sex == 'f')
climb2f10=subset(climb2f, day =='10')
climb2f22=subset(climb2f, day =='22')

dim(climb2)
dim(climb2m)
dim(climb2f)
dim(climb2m10)
dim(climb2m22)
dim(climb2f10)
dim(climb2f22)

## set margin inches "mai" to accommodate vertical x label names
par("mai" = c(3, 1.25, 0.5, 0.5))
par(mgp=c(3,1,0))
boxplot(slope~sex*day*mito*geno, data = climb2, las=2, xlab="", main="Climbing Speed by Sex (Purple, Green), Age, Mito Nuclear genotype", ylab="Climbing Speed (slope)", col=c("purple","green"))


## Fit ANOVA models
## fit1m <- aov(slope~mito*geno*day,data = climb1m)
## summary(fit1m)

fit2m <- aov(slope~mito*geno*day,data = climb2m)
summary(fit2m)

## fit1f <- aov(slope~mito*geno*day,data = climb1f)
## summary(fit1f)

fit2f <- aov(slope~mito*geno*day,data = climb2f)
summary(fit2f)

#Tukeym <- TukeyHSD(fit1m)
#Tukeym

# No females in tis data set - fitf <- aov(slope~mito*geno,data = climbf)
# No females in tis data set - summary(fitf)
# No females in tis data set - Tukeyf <- TukeyHSD(fitf)
# No females in tis data set - Tukeyf

#model2 = lm(slope ~ mito*geno*sex,data = subset(climb, mito != 'si1' & geno != 'w1118'))
#Defs=subset(climb, mito != 'si1' & geno != 'w1118')
#summary(Defs)
#model2 = lm(slope ~ mito*geno*sex,data = subset(climb, mito != 'si1' & geno != 'w1118'))
#model3 = lm(slope ~ mito*geno*sex,data = subset(climb, mito != 'si1' & mito != 'ore'))
#modelm = lm(slope ~ mito*geno*day,data = subset(climbm, mito != 'si1' & mito != 'ore'))
#modelf = lm(slope ~ mito*geno*day,data = subset(climbf, mito != 'si1' & mito != 'ore'))

#Anova(modelf, type = "II")
#summary(modelf)


interaction.plot(x.factor     = climbm12$geno,
                 trace.factor = climbm12$diet,
                 response     = climbm12$slope,
                 fun = mean,
                 type="b",
                 col=c("black","red","green", "blue"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15, 16),     ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 lty = 1,  # line type
                 lwd = 2,  # line width
                 #leg.bt = "o",
                 las=0, xlab="Nuclear Genotype", main="Male Climbing Speed for Day27 by Diet", ylab="Climbing Speed (slope)")

#glimpse(climbf) not a function

## ggplot box plots from https://sebastiansauer.github.io/vis_interaction_effects/
ggplot(climbm) +
  aes(x = genodiet, y = slope) +
  geom_boxplot() +
  facet_wrap(~day)

## ggplot interactions from https://ggplot2tutor.com/tutorials/interaction_plot

## MALE Df background Facet wrap by 'day'
## climb2mtrim<-subset(climb2m, slope <5) ## TRIMMING outlier slope values  - does NOT fix w1118;7837 outlier
climb2m %>% 
  group_by(mito,geno,day) %>% 
  #summarise(mnslope = mean(slope)) -> climbm2 Old code
  summarise(mnslope = ci(slope)[1],  # Joaquin's BETTER version
            mnslope_low = ci(slope)[2],
            mnslope_high = ci(slope)[3]) -> climb2m2day


climb2m2day %>% 
  ggplot() +
  aes(x = geno, y = mnslope, ymin = mnslope_low, ymax = mnslope_high, color = mito) +  ## Joaquin's BETTER version
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  facet_wrap(~day) +
  #scale_fill_discrete(breaks=c('w1118', 'Zim', 'siI', 'sm21',  'yak')) +
  #theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  ylim(0.5,4.0) +
  labs(
    title = "Male Climbing of Mito-Deficiency Genotypes at Ages 10 or 22 Days",
    subtitle = paste0("Significant Mito x Nuclear Interaction"),
    x = "Deficiency Genotype",
    y = "Climbing Speed (cm/second)"
  ) 

## FEMALE Cy OR Df background Facet wrap by 'day'
## climb2ftrim<-subset(climb2f, slope <5)
climb2f %>% 
  group_by(mito, geno, day) %>% 
  #summarise(mnslope = mean(slope)) -> climbm2 Old code
  summarise(mnslope = ci(slope)[1],  # Joaquin's BETTER version
            mnslope_low = ci(slope)[2],
            mnslope_high = ci(slope)[3]) -> climb2f2day


climb2f2day %>% 
  ggplot() +
  aes(x = geno, y = mnslope, ymin = mnslope_low, ymax = mnslope_high, color = mito) +  ## Joaquin's BETTER version
  geom_line(aes(group = mito)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  facet_wrap(~day) +
  #scale_fill_discrete(breaks=c('w1118', 'Zim', 'siI', 'sm21',  'yak')) +
  #theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  ylim(0.5,4.0) +
  #theme_classic() +
  labs(
    title = "Female Climbing of Mito-Deficiency Genotypes at Ages 10 or 22 Days",
    subtitle = paste0("Significant Mito x Nuclear Interactions"),
    x = "Deficiency Genotype",
    y = "Climbing Speed (cm/sec)"
  ) 


## Subset specific pairs of genotypes to assess interaction with diet.

climbno315=subset(climb, geno != '315')
dim(climbno315)
climbno31580133=subset(climbno315, geno != '80133')
climb765x80133=subset(climbno31580133, geno != '80133x315')

dim(climbno31580133)
dim(climb765x80133)

fit765x80133 <- aov(slope~day*geno*diet,data = climb765x80133)
summary(fit765x80133)
Tukey765 <- TukeyHSD(fit765x80133)
Tukey765

## ggplot box plots from https://sebastiansauer.github.io/vis_interaction_effects/
ggplot(climb765x80133) +
  aes(x = genodiet, y = slope) +
  geom_boxplot() +
  facet_wrap(~day)

interaction.plot(x.factor     = climb765x80133$geno,
                 trace.factor = climb765x80133$diet,
                 response     = climb765x80133$slope,
                 fun = mean,
                 type="b",
                 col=c("black","red","green", "blue"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15, 16),     ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 lty = 1,  # line type
                 lwd = 2,  # line width
                 #leg.bt = "o",
                 las=0, xlab="Nuclear Genotype", main="Male Climbing Speed for Day27 by Diet", ylab="Climbing Speed (slope)")

## Joaquin code for distributions in interaction plot, by facets 8/4/2021
install.packages("hrbrthemes")
install.packages("ggdist")
library(hrbrthemes)
library(ggdist)

#climbm %>%
climb %>%
  ggplot(aes(
    x=geno,
    y=slope,
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
  ylab("Climbing Speed (Slope)") +
  xlab("Nuclear Genotype") +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(sex~day) +  # joint sex x time effect in facests  see top of code block
  labs(
    title = "Female and Male Climbing Speed by MitoNuclear Genotype at Ages 10 or 24 Days",
    subtitle = paste0("Females show more mtDNA effects and some genotypes get faster at older ages\n",
                      "Males get slower with age and show little mtDNA effects or Mito x Nuclear Interaction"),
    x = "Nuclear  Genotype",
    y = "Climbing Speed (slope)"
  )


## lines below error : data must be in data frame, not and S3 object.  Did this pass 'aes() to the 'data' argument?
ggplot(aes(climbm, slope)) +
  geom_line(size = 1.2, aes(group = geno, color = geno)) +
  geom_point(size = 2.6, aes(color = geno), shape = 15)

