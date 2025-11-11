# title: "Starvation Analysis"
# output: html_notebook

setwd('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mothers Curse Pervasive Absence/Data Analyses Mothers Curse/Starvation Analyses/')
#install.packages('googlesheet4')
install.packages('survival')
install.packages('survminer')
install.packages('tidyverse')
install.packages("magrittr")
install.packages("dplyr") 
install.packages("psych")

library(psych)
library(magrittr)
library(dplyr)
library(tidyverse)
library(gmodels)
library(ggplot2)
#library(googlesheets4)
library(survival)
library(survminer)

## starvation.data<- read_sheet('https://docs.google.com/spreadsheets/d/1khfioHB9rKGaau5ptkCNT8DZqWOhfmVwhgfgRI8Yv9k/edit#gid=1487910965', skip = 3)
## read downloaded googlsheet as a plain csv - from Shawns Google drive. Used CSV on hard drive in Shawn-David Drafts/Starvation folder
## starvation.data<- read_csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/MitoStarvationw1118.csv', skip = 3)

## combined the three MitoStarvation...csv files for 7744, 8469 and w1118 from Shawn-Kenny-Lindsay file to read alll three at once
starvation.data<- read_csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Mothers Curse Pervasive Absence/Data Analyses Mothers Curse/Starvation Analyses/MitoStarvation_7744_8469_w1118.csv', skip = 3)
starvation.data %>%
  mutate_at(.vars = c('mito','nuc', 'sex','rep'),.funs = factor) %>%
  gather(Time,Death, -c('mito','nuc', 'sex','rep')) %>%
  mutate(Time = sub(pattern = 'hr',replacement = '', x = Time))%>%
  mutate_at(.vars = 'Time',.funs = as.double) %>%
  filter(!is.na(Death)) %>%
  filter(Time != 0) -> stv
stv
head(stv)
dim(stv)

# write to a new csv for use in Lei's program that estimates survival data from HHMI 2014 project
write_csv(stv, '7744_8469_w1118.csv')

Dead = stv$Death # creates a variable with the number of deaths
survFormat.DF <- data.frame() # initialize empty data 
#frame to store the expanded tables 


for (iter in 1:dim(stv)[1]){
  i=0
  while (i < Dead[iter]){
    i = i+1 
    TMP<- stv[iter,]
    TMP$deaths=1
    survFormat.DF<- rbind(survFormat.DF,TMP) }
}
dim(survFormat.DF)
stv

## NOTE: This code was originally for w1118 only. Edited to include all 3 Dfs: 
## surv.resultsw1118 edited to surv.results7744_8469_w1118 and Survw1118.csv edited to Surv7744_8469_w1118.csv
survfit(Surv(Time, deaths == 1) ~ mito + nuc + sex, data= survFormat.DF) -> surv.results7744_8469_w1118
coxph(Surv(Time, deaths == 1) ~ mito + nuc + sex, data= survFormat.DF)
levels(survFormat.DF$mito)
levels(survFormat.DF$nuc)
levels(survFormat.DF$sex)
summary(surv.results7744_8469_w1118)$table %>% as_tibble(rownames = 'condition') %>% write_csv('Surv7744_8469_w1118.csv')
head('Surv7744_8469_w1118.csv')

## write_csv(surv.results7744_8469_w1118$
## write_excel_csv(surv.results7744_8469_w1118)$table
plot(surv.results7744_8469_w1118)
ggsurvplot(surv.results7744_8469_w1118)

survival_df <- function(x){
  x.sum <- summary(x) # summarize fit 
  Condition <- x.sum$strata %>%
    as.character() 
  
  # create data.frame with points to plot 
  df <- data.frame(Condition, time = x$time, survival = x$surv) 
  df <- mutate(df, Condition = factor(Condition))
  # add 0 time columns t df 
  for ( levs in levels(df$Condition)){
    df <- df %>%
      add_row(Condition = levs, time = 0, survival = 1 )
  }
  return(df)
}

surv.delevel<- function(x){
  sub('.*=','',x,perl=T)
}

test1<-survival_df(surv.results7744_8469_w1118) %>%
  separate(Condition, into = c('mito', 'nuc', 'sex'),sep = ", ") %>%
  mutate(mito = factor(surv.delevel(mito))) %>%
  mutate(nuc = factor(surv.delevel(nuc))) %>%
  mutate(sex = factor(surv.delevel(sex))) %>%
  
  ggplot(aes(x = time, y = survival, color = mito)) +
  geom_line(size = 0.5) +
  #xlab('Time (hrs)')+
  #ylab('Survival')+
  facet_wrap(~ sex + nuc, ncol = 3)+
  theme_bw(base_size = 11 )+
  theme(panel.grid = element_blank(), legend.position = 'right')

ggsave('mitoDfstarve.png', test1, dpi=300)

## Interaction plots of all three Deficiencies so far : 7744, 8469, w1118
## NOTE: this requires output CSVs for each nuclear genotype that generates the survival parameters
## 7744.csv, 8469.csv, w1118.csv. Code above generates these individually. Could be appended to do this once for : SurvStats3Dfs.csv.

starve3dfs<- read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/SurvStats3Dfs.csv')
starve3dfs$sex= as.factor(starve3dfs$sex)
starve3dfs$nuc= as.factor(starve3dfs$nuc)
starve3dfs$mito= as.factor(starve3dfs$mito)
#starve3dfs$replicate= as.factor(starve3dfs$replicate)

head(starve3dfs)
dim(starve3dfs)
starvem3dfs<-subset(starve3dfs, sex == 'Male')
starvef3dfs<-subset(starve3dfs, sex == 'Female')

interaction.plot(x.factor     = starvem3dfs$nuc,
                 trace.factor = starvem3dfs$mito,
                 response     = starvem3dfs$Mean,
                 fun = mean,
                 type="b",
                 col=c("black","red","green", "blue", "pink"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15, 16, 18),     ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 ylim = range(55:110, na.rm = TRUE),
                 leg.bty = "o",
                 las=1, xlab="Nuclear Genotype", main="Male Starvation Resistance by Mitonuclear Genotype", ylab="Survival (hours)")

## trying to facet wrap by variable 'sex'
starve3dfs %>% 
  group_by(sex, mito, nuc) -> starve3dfs2
#summarise(meanstarve = mean(Mean)) -> starve3dfs2

head(starve3dfs2)

starve3dfs2 %>% 
  ggplot() +
  aes(x = nuc, y = Mean, color = mito) +
  geom_line(aes(group = mito)) +
  geom_point() +
  facet_wrap(~sex) +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90)) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=12, angle=0),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=12, angle=0)) +
  labs(
    title = "Female and Male Starvation Resistance by MitoNuclear Genotype",
    x = "Nuclear Genotype",
    y = "Starvation Resistance (hours)"
  ) 
## Cut out this subtitle from ggplot frame;
##     subtitle = paste0("Females show considerable mtDNA and MitoNuclear interaction effects\n",
## "Males are more sensitive to starvation and show limited genetic effects"),

## Leah's Code for producing Mother's Cures biplots for Climbing and Flight raw csvs. Uses raw CSVs to calculate means and CVs
## Read in edited file for CVs 
cv_by_sex_starve <- read_csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/cv_by_sex_starve.csv')
ggplot(cv_by_sex_starve, aes(x = cv_female, y = cv_male, color = nuc))+
  geom_point(size=5) +
  geom_abline(intercept = 0, Y = 1) +
  xlim(0,.3) +
  ylim(0,.3) +
  theme_bw() +
  labs(color="Nuclear Genotype")

## Code from below: PLOT Female x Male CVs 
ggplot() + geom_point(size=5, aes(x=stvfmitomeans2.df$cv, y=stvmmitomeans2.df$cv)) +
  labs(
    title = "Female vs. Male mtDNA Coefficient of Varitation for Starvation",
    subtitle = paste0("3 mtDNAs in Nuclear Dfs 7744, 8469, w1118"),
    x = "Female Mitonuclear CV",
    y = "Male Mitonuclear CV)") +
  xlim(0,0.2) + ylim(0,0.2)


## CODE FOR CALCULATING Coefficient of Variation 
library(psych)
## MALE Code - NOTE SurvStats3Dfs_males.csv was sourcesd from SurvStats3Dfs.csv values providing means BUT NOT s.d. 
## flightmmeans <- describeBy(Y ~ mitonuc, data = flightm) ## Actually need means separated by mito to get mito-CV in each condition
## flightmmeans
## flightmmeans2 <- do.call("rbind", flightmmeans)
## flightmmeans2
## flightmmeans.df <- data.frame(flightmmeans2)
## flightmmeans.df$cv = flightmmeans.df$sd / flightmmeans.df$mean

## write to a file that allows taking means of mtDNAs across nuc and day
## write.csv(flightmmeans.df, '/Users/davidrand/Dropbox (Brown)/Manuscripts/Mito nuc Screen/Mito nuc Flight/batch_21_07_30/flightmmeans.csv')
## Edit this by hand in Excel = malemeans to have a NUC_DAY column added for CV calculation.
## READ this hand-edited file to allow mito means for mito-CV calculation 
stvmmeans<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/SurvStats3Dfs_males.csv')
stvmmitomeans <- describeBy(Mean ~ nuc, data = stvmmeans)
stvmmitomeans2 <- do.call("rbind", stvmmitomeans)
stvmmitomeans2.df <- data.frame(stvmmitomeans2)
stvmmitomeans2.df$cv = stvmmitomeans2.df$sd / stvmmitomeans2.df$mean
stvmmitomeans2.df

## FEMALE Code - NOTE SurvStats3Dfs_fmales.csv was sourcesd from SurvStats3Dfs.csv values providing means BUT NOT s.d. 
stvfmeans<-read.csv('/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/SurvStats3Dfs_females.csv')
stvfmitomeans <- describeBy(Mean ~ nuc, data = stvfmeans)
stvfmitomeans2 <- do.call("rbind", stvfmitomeans)
stvfmitomeans2.df <- data.frame(stvfmitomeans2)
stvfmitomeans2.df$cv = stvfmitomeans2.df$sd / stvfmitomeans2.df$mean
stvfmitomeans2.df

write.csv(stvmmitomeans2.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/stvmaleCVs.csv' )
write.csv(stvfmitomeans2.df, '/Users/davidrand/Brown Dropbox/David Rand/Manuscripts/Shawn Williams/Shawn-David Drafts/Starvation_KL_DRM/Starvation Analysis/stvfemaleCVs.csv' )


## PLOT Female x Male CVs 
ggplot() + geom_point(size=5, aes(x=stvfmitomeans2.df$cv, y=stvmmitomeans2.df$cv)) +
  labs(
    title = "Female vs. Male mtDNA Coefficient of Varitation for Starvation",
    subtitle = paste0("3 mtDNAs in Nuclear Dfs 7744, 8469, w1118"),
    x = "Female Mitonuclear CV",
    y = "Male Mitonuclear CV)") +
  xlim(0,0.2) + ylim(0,0.2)

