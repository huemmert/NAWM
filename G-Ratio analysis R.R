 #####                                 ####
###     MASTER THESIS PROJECT 2021    ###
####                                  ####


# G-RATIO analysis Con vs. MS

library(dplyr)
library(ggplot2)
library("gridExtra")
library(reshape2)
library("RColorBrewer")
library("ggsci")
library(car)
 library(pgirmess)
 library(ggpubr)

  ####                                      ####  
 ###             G-Ratio ANALYSIS         ###  
 ####                                      ####

 setwd("F:/TEM ANALYSIS R/data")

# load data 
xdata<-read.csv("G_ratio_MS_cleaned_complete_asc_age.csv")
xdata<-filter(xdata, !is.na(Axon_caliber))
xdata<-filter(xdata, !is.na(g_ratio))

# how many axons per patient were measured?
patient19_26<-filter(xdata, Patient=="19-26", !is.na(Axon_caliber))
patient19_106<-filter(xdata, Patient=="19-106", !is.na(Axon_caliber))
patient19_36<-filter(xdata, Patient=="19-36", !is.na(Axon_caliber))
patient19_47<-filter(xdata, Patient=="19-47", !is.na(Axon_caliber))
patient20_52<-filter(xdata, Patient=="20-52", !is.na(Axon_caliber))
patient20_54<-filter(xdata, Patient=="20-54", !is.na(Axon_caliber))
patient20_58<-filter(xdata, Patient=="20-58", !is.na(Axon_caliber))
patient20_60<-filter(xdata, Patient=="20-60", !is.na(Axon_caliber))
patient20_64<-filter(xdata, Patient=="20-64", !is.na(Axon_caliber))
patient20_65<-filter(xdata, Patient=="20-65", !is.na(Axon_caliber))
patient20_71<-filter(xdata, Patient=="20-71", !is.na(Axon_caliber))
patient20_77<-filter(xdata, Patient=="20-77", !is.na(Axon_caliber))
patient20_82<-filter(xdata, Patient=="20-82", !is.na(Axon_caliber))
patient20_92<-filter(xdata, Patient=="20-92", !is.na(Axon_caliber))
patient19_27<-filter(xdata, Patient=="19-27", !is.na(Axon_caliber))
patient19_31<-filter(xdata, Patient=="19-31", !is.na(Axon_caliber))
patient20_95<-filter(xdata, Patient=="20-95", !is.na(Axon_caliber))
patient19_77<-filter(xdata, Patient=="19-77", !is.na(Axon_caliber))
patient19_69<-filter(xdata, Patient=="19-69", !is.na(Axon_caliber))
patient19_98<-filter(xdata, Patient=="19-98", !is.na(Axon_caliber))
count(patient19_106, "axon_caliber") # n= 193
count(patient19_26, "axon_caliber") # n= 216
count(patient19_36, "axon_caliber") # n= 149
count(patient19_47, "axon_caliber") # n= 268
count(patient20_52, "axon_caliber") # n= 119
count(patient20_54, "axon_caliber") # n= 197
count(patient20_58, "axon_caliber") # n= 204
count(patient20_60, "axon_caliber") # n= 245
count(patient20_64, "axon_caliber") # n= 150
count(patient20_65, "axon_caliber") # n= 182
count(patient20_71, "axon_caliber") # n= 210
count(patient20_77, "axon_caliber") # n= 212
count(patient20_82, "axon_caliber") # n= 186
count(patient20_92, "axon_caliber") # n= 171
count(patient19_27, "axon_caliber") # n= 251
count(patient19_31, "axon_caliber") # n= 239
count(patient20_95, "axon_caliber") # n= 158
count(patient19_77, "axon_caliber") # n= 198
count(patient19_69, "axon_caliber") # n= 182
count(patient19_98, "axon_caliber") # n= 156


# GRAPH 1: POINT CLOUD
# do scatter plot of g-ratio (basic, just condition and g-ratio)
scatter_g_ratio <- ggplot(excluded_data, aes(Axon_caliber, g_ratio, col=Condition)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(aes(color=Condition),method = "lm", se = F)+
  scale_y_continuous()+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="g-ratio of human optic nerve tissue",
       x="Axon caliber [µm]", 
       y="g-ratio", col="condition")


# GRAPH 2: boxplot to show range of g-ratios per Patient
# reordering patients by hand to increasing age
desired_order<-c("19-31","19-26","20-60", "19-98", "19-27", "19-69", "19-47", "20-65", "20-64", "20-77","20-71", "19-106", "20-54", "19-36", "20-82","20-95", "20-92", "19-77", "20-58", "20-52")
xdata$Patient<-factor(as.character(xdata$Patient), levels=desired_order)

g_ratio_boxplot<-ggplot(xdata, aes(Patient, g_ratio, color=Condition))+
  stat_boxplot(aes(Patient, g_ratio),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Patient, g_ratio), outlier.shape=1)+
  #geom_violin(alpha = 0.5) +
  annotate(geom='text', 
           x=xdata$Patient, 
           y=0.07, 
           label=xdata$Age,
           size=2.5, 
           alpha=0.2, 
           fontface="plain")+
    geom_text(label="Patient age", # adds the patient age row below x-axis
              y=0.07, 
              x=0.3, 
              color="black", 
              size=2.5, 
              alpha=0.2)+ 
  coord_cartesian(ylim = c(0.1, 1), expand = FALSE, clip = "off") +
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.1)+
  scale_color_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + #complete cols "#003264",'#999999','#E69F00','#719c70'
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Individual g-ratio distributions", 
       x="Patient", 
       y="g-ratio",
       size=2.5, 
       fontface="bold")



# GRAPH 3: boxplot to show conditions vs. g-ratio 
compare_means(g_ratio~Condition, data=xdata, method="t.test") # used to implement significance into the graph below
my_comparison<-list(c("CON","MS"))#,c("CON","FTD"), c("CON", "PD"))

g_ratio_boxplot_per_condition<-ggplot(excluded_data, aes(Condition, g_ratio, color=Condition, palette="jco"))+
  stat_boxplot(aes(Condition, g_ratio),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, g_ratio), outlier.shape=1)+
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) + #complete colors'#999999','#719c70','#E69F00',"#003264"
  stat_compare_means(method = "t.test", label.y = 1.0)+
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="g-ratios per condition", 
       x="Condition", 
       y="g-ratio",
       size=2.5, 
       fontface="bold")

# GRAPH 4:do the same but look at axon calibers
compare_means(Axon_caliber~Condition, data=xdata, method="wilcox.test") # used to implement significance into the graph below
my_comparison1<-list(c("CON","MS"),
                     c("CON","FTD"), 
                     c("CON", "PD")) # which tests i want to visualize
# boxplot to show conditions vs. g-ratio 
Axon_caliber_boxplot_per_condition<-ggplot(excluded_data, aes(Condition, 
                                               Axon_caliber, 
                                               color=Condition, 
                                               palette="jco"))+
  stat_boxplot(aes(Condition, Axon_caliber),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, Axon_caliber), outlier.shape=1)+
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) +
  stat_compare_means(method="wilcox.test", label.y=6.9)+
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Axon calibers per condition", 
       x="Condition", 
       y="Axon caliber [µm]",
       size=2.5, 
       fontface="bold")

# plot g-ratio vs. age (ascending)
scatter_g_ratio_vsage <- ggplot(excluded_data, aes(Age, g_ratio, col=Condition)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(aes(color=Condition),method = "lm", se = F)+
  scale_y_continuous()+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="g-ratio correlation with patient age ",
       x="Age [years]", 
       y="g-ratio", col="condition")

# plot g-ratio vs. PMD time 
scatter_g_ratio_vs_PMD<- ggplot(excluded_data, aes(PMD_.min., g_ratio, col=Condition)) + 
  geom_point(alpha=0.4) + 
  geom_smooth(aes(color=Condition),method = "lm", se = F)+
  scale_y_continuous()+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="g-ratio correlation with PMD ",
       x="Post mortem delay [min]", 
       y="g-ratio", col="condition")
#-------------------------------------------------------------------------
###
### STATISTICS 
###

# subset into conditions and test for normality 
MS<-subset(xdata, Condition=="MS")
COntrol<-subset(xdata, Condition=="CON")
FTD<-subset(xdata, Condition=="FTD")
PD<-subset(xdata, Condition=="PD")

# check correlations for the cloud graph 

# exclude patient 19-31 because of remyelinated axons
MS_less<-subset(MS, Patient!="19-31")

qqnorm(MS_less$g_ratio)
qqnorm(COntrol$g_ratio) # look pretty normal 

shapiro.test(MS_less$g_ratio) # p<0.05 
shapiro.test(COntrol$g_ratio) # p<0.05
shapiro.test(FTD$g_ratio) # p<0.05
shapiro.test(PD$g_ratio) # p<0.05

qqnorm(MS$Axon_caliber) # look non-parametric
qqnorm(COntrol$Axon_caliber)

shapiro.test(MS$Axon_caliber) # p<0.05
shapiro.test(COntrol$Axon_caliber) # p<0.05
shapiro.test(FTD$Axon_caliber) # p<0.05
shapiro.test(PD$Axon_caliber) #p<0.05, everything is non-parametric

# do parametric test for g-ratios (look normal in qqplot)
t.test(MS_less$g_ratio, COntrol$g_ratio) # p<0.05, so significant 
t.test(COntrol$Axon_caliber, MS_less$Axon_caliber) # p>0.05, so NOT significant 
t.test()

#### TEST if the g-ratios are statistically significantly different between the condidtions FTD, PD, CON and MS

# do non-parametric Kruskal wallis 
kruskal.test(g_ratio~Condition, data=xdata) # p<0.05, g_ratios are significantly different 
kruskal.test(Axon_caliber~Condition, data=xdata) # p>0.05

# re-do kruskal test with now only CON (with FTD & PD) and MS 

# post-hoc analysis to determine which of the groups are different # use p adjust benjamini hochberg like jonas
pairwise.wilcox.test(xdata$g_ratio,
                     xdata$Patient,
                     p.adjust="BH")

g_ratio_wilcox_test <- as.data.frame(pairwise.wilcox.test(xdata$g_ratio,
                                                          xdata$Patient,
                                                          p.adjust="bonferroni")$p.value)

write.csv(g_ratio_wilcox_test, "Pairwise_wilcoxon_g_ratios_bonferroni.csv")

# Result: Significant difference between CON-MS, CON-PD, MS-FTD, MS-PD
pairwise.wilcox.test(xdata$Axon_caliber, 
                     xdata$Condition, 
                     p.adjust="bonferroni")
# Result: Significant difference only between FTD-CON, MS-FTD, FTD-PD

# do statistics like Jonas (Wilcoxon rank sum for each patient g-ratio vs. each other)
wilcox.test(xdata$g_ratio~xdata$Patient)


sophie_stinks <- pairwise.wilcox.test(xdata$g_ratio,
                                      xdata$Patient,
                                      p.adjust="bonferroni")

# statistics kolmogorov smirnov for each patient & Axon caliber distribuiton
ks.test(xdata$Patient, xdata$Axon_caliber)

# IT stats
shapiro.test(MS_less$InnerTongue_Area) # p<0.05
shapiro.test(COntrol$InnerTongue_Area)
qqnorm(MS$InnerTongue_Area) # look non-parametric
qqnorm(COntrol$InnerTongue_Area) # look non-parametric
kruskal.test(Condition~InnerTongue_Area, data=xdata) # p>0.05, not signif


#------------------------------------
####
#### RE-DO STATISTICS WITHOUT 19-31, 20-65, 19-26
####

excluded_data<-subset(xdata, Patient!="19-31" & Patient!="19-26" & Patient!="20-65")

MS<-subset(excluded_data, Condition=="MS")
CON<-subset(excluded_data, Condition=="CON")
qqnorm(MS$g_ratio) # looks normal
shapiro.test(MS$g_ratio) # p<0.05
qqnorm(CON$g_ratio) # looks normal
shapiro.test(CON$g_ratio) # p<0.05
shapiro.test(MS$Age) # p<0.05 (non normal)
shapiro.test(CON$Age) # p<0.05 
shapiro.test(MS$PMD_.min.) # p<0.05
shapiro.test(CON$PMD_.min.) # p<0.05

# check for equality of variance (F test)
var.test(MS$g_ratio, CON$g_ratio, alternative = "two.sided")

cor.test(MS$Axon_caliber, MS$g_ratio, method="spearman") # p<2.2e-16, rho = 0.55
cor.test(CON$Axon_caliber, CON$g_ratio, method="spearman") #p<2.2e-16, rho = 0.58

# check if ageing has a correlation to g.ratio
cor.test(MS$g_ratio, MS$Age, method = "spearman") # p=0.0493, rho=0.04
cor.test(CON$g_ratio, CON$Age, method="spearman") # p<0.05, rho =-0.101

# check if PMD has a correlation to g-ratio
cor.test(MS$g_ratio, MS$PMD_.min., method ="spearman") # p<0.05, rho = 0.234
cor.test(CON$g_ratio, CON$PMD_.min., method="spearman") # p>0.05, rho = 0.035


# do parametric t.test because g-ratios look normal 
t.test(MS$g_ratio, CON$g_ratio) # p<0.05 signif different 

# look at axon calibers
qqnorm(MS$Axon_caliber) # looks abnormal
qqnorm(CON$Axon_caliber) # also looks abnomal
shapiro.test(MS$Axon_caliber) #p<0.05 , non parametric
shapiro.test(CON$Axon_caliber) # p<0.05
wilcox.test(excluded_data$Axon_caliber~excluded_data$Condition) # wilcoxon rank sum for comparing two groups of non parametric values
# p=0.3788, not signif diff

# look at IT area
qqnorm(MS$InnerTongue_Area) # look non-parametric
qqnorm(CON$InnerTongue_Area) # look non-parametric
shapiro.test(MS$InnerTongue_Area) #p<0,05
shapiro.test(CON$InnerTongue_Area) # p<0,05
wilcox.test(excluded_data$InnerTongue_Area~excluded_data$Condition)

# look at PAS 
excluded_PAS<-subset(PAS, Patient!="19-31" & Patient!="19-26" & Patient!="20-65" )
MS_PAS<-subset(excluded_PAS, Condition=="MS")
CON_PAS<-subset(excluded_PAS, Condition=="CON ")
qqnorm(MS_PAS$Average_PS) # non-normal
qqnorm(CON_PAS$Average_PS) # non-normal
shapiro.test(MS_PAS$Average_PS) #non-normal
shapiro.test(CON_PAS$Average_PS) #non-normal
wilcox.test(excluded_PAS$Average_PS~excluded_PAS$Condition) # p<0.05, sign diff

# correlate to PMD
cor.test(MS_PAS$Average_PS, MS_PAS$PMD_.min., method ="spearman") # p = 0.1245, rho = -0.11
cor.test(CON_PAS$Average_PS, CON_PAS$PMD_.min., method ="spearman") # p = 0.0198, rho = 0.163

# correlate to Age
MS_PAS<-filter(MS_PAS, !is.na(Axons_per_Image) & !is.na(Age) & !is.na(Average_PS))
cor.test(MS_PAS$Age, MS_PAS$Average_PS, method ="spearman") # p < 0.05, rho = -0.407
cor.test(CON_PAS$Age, CON_PAS$Average_PS, method = "spearman") # p< 0.05, rho = 0.266


######
###   ANOVA task
#####

#### check influence of AGE
# test for equal variances
bartlett.test(g_ratio ~ Condition, data=excluded_data) # p <0.05, assumption of equal variance is violated
bartlett.test(g_ratio~Age, data = excluded_data) # p<0.05

# perform ANOVA

excluded_data$age_group<-findInterval(excluded_data$Age, c(60,70,80,90,100,110)) # groups data by age classes every 10 years
bartlett.test(g_ratio~as.factor(age_group), data = excluded_data) # p>0.05, assumption correct, variance is equal
g_ratio_age_group_aov<-aov(g_ratio~as.factor(age_group), data = excluded_data) 
summary(g_ratio_age_group_aov) # p<0.05, significantly different between age groups
plot(g_ratio_age_group_aov)

# do posthoc tukey test
TukeyHSD(g_ratio_age_group_aov, conf.level = .95) # p<0.05 for 3vs1, 4vs1, 5vs1, 3vs2, 5vs2, 5vs3 and 5vs4
plot(TukeyHSD(g_ratio_age_group_aov, conf.level = .95),las=2)

t.test(age_group~Condition, data=excluded_data) # p<0.05 significantly different.

### check influence of GENDER 
bartlett.test(g_ratio~Gender, data = excluded_data) # p>0.05 assumption met
g_ratio_gender_aov<-aov(g_ratio~Gender, data = excluded_data)
summary(g_ratio_gender_aov) # p=0.000559 significantly different
plot(g_ratio_age_group_aov)

dev.off()

### check influence of PMD
excluded_data$PMD_group<-findInterval(excluded_data$PMD_.min., c(250,300,350,400,450,500,550,600)) 
bartlett.test(g_ratio~PMD_group, data = excluded_data) # p<0.05, assumption violated
# do welch instead of anova
oneway.test(g_ratio~PMD_group, data = excluded_data, var.equal = F) # p<0.05, signif diff
t.test(PMD_group~Condition, data = excluded_data) # p<0.05, signif

####
###   LINEAR MODELs
# do LM to see average_g_ratio ~ age / PMD / Gender and look at R²
average_data<-read.csv("Average g-ratio_complete.csv")
excluded_average_data<-subset(average_data, Patient!="19-31" & Patient!="19-26" & Patient!="20-65")
lm_model_age<-lm(Average_g_ratio~Age, data=excluded_average_data)
summary(lm_model_age) # p>0.05, R² = 0.1198

lm_model_Gender<-lm(Average_g_ratio~Gender, data=excluded_average_data)
summary(lm_model_Gender) # Gender Male p>0.05, R² = 0.52

lm_model_condition<-lm(Average_g_ratio~Condition, data=excluded_average_data)
summary(lm_model_condition) # p>0.05, R²=0.1715

lm_model_experiment<-lm(Average_g_ratio~Condition+Age+PMD+Gender, data=excluded_average_data)
summary(lm_model_experiment) # p>0.05, R²=0.1627
lm_model_experiment1<-lm(g_ratio~Condition+Age+PMD_.min.+Gender, data=excluded_data)
summary(lm_model_experiment1) # all variables are significant?

lm_model_PMD<-lm(Average_g_ratio~PMD, data=excluded_average_data)
summary(lm_model_PMD) # p>0.05, R² = 0.04

# try again for PMD but with all g_ratio data points instead of averages
lm_model_whole<-lm(g_ratio~PMD_.min., data=excluded_data)
summary(lm_model_whole) #p<0.05, R² = 0.021 WHY IS RESULT DIFFERENT WHEN USING ALL DATA POINTS
#-------------------------------------------------------------------

##### 
###     Analyse INNER TONGUE #####
#####


# import the file from excel somehow

library(readr)
newdata<-read.csv("Average g-ratio_complete.csv")
newdata1<-newdata # reorders data into decreasing bars by copying original data and using reorder()
newdata1$Patient<-factor(newdata1$Patient, 
                         levels=newdata1$Patient[order(newdata1$Percentage_IT, 
                                                       decreasing=T)])

# GRAPH 5: and now plot the damn thing 
damnplot<-ggplot(newdata1, aes(Patient, Percentage_IT, fill=Condition))+
  geom_bar(stat="identity") + 
  #geom_smooth(aes(color=Condition),method = "lm", se = F)+
  #scale_y_continuous()+
  annotate(geom='text', 
           x=newdata1$Patient, 
           y=-3.8, 
           label=newdata1$Age,
           size=2.5, 
           alpha=0.8, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-3.8, 
            x=0.3, 
            color="black", 
            size=2.5, 
            alpha=0.2)+ 
  coord_cartesian(ylim = c(0, 80), expand = FALSE, clip = "off") +
  scale_fill_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + 
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Percentage of Inner tongue visible in myelinated Axons", 
       x="Patient",
       y="Visible Inner tongue [%]")
  
#-----------------------------------------------------------
####
####      plot the Area of Inner tongue 
####

Tongue_data<-filter(xdata, InnerTongue_Area>0) # exclude all dataPoints in which area = 0

# plot individual IT areas for each patient 
IT_area_boxplot<-ggplot(Tongue_data, aes(Patient, InnerTongue_Area, color=Condition))+
  stat_boxplot(aes(Patient, InnerTongue_Area),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Patient, InnerTongue_Area), outlier.shape=NA)+
  annotate(geom='text', 
           x=Tongue_data$Patient, 
           y=-0.5, 
           label=Tongue_data$Age,
           size=2.5, 
           alpha=0.2, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-0.5, 
            x=0.3, 
            color="black", 
            size=2.5, 
            alpha=0.2)+ 
  coord_cartesian(ylim = c(0, 11), expand = FALSE, clip = "off") +
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.25)+
  scale_color_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + #complete cols "#003264",'#999999','#E69F00','#719c70'
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Individual Inner tongue area distributions", 
       x="Patient", 
       y="Inner tongue area [µm]",
       size=2.5, 
       fontface="bold")

# now with just the conditions
Tongue_data_excluded<-subset(Tongue_data, Patient!="19-31"&Patient!="20-65"&Patient!="19-26")
IT_area_boxplot_CONMS<-ggplot(Tongue_data_excluded, aes(Condition, 
                                               InnerTongue_Area, 
                                               color=Condition, 
                                               palette="jco"))+
  stat_boxplot(aes(Condition, InnerTongue_Area),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, InnerTongue_Area), outlier.shape=NA)+
  stat_summary(fun=mean, geom="point", size=4)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) +
  stat_compare_means(method="wilcox.test", label.y=10, label.x=0.9)+
  coord_cartesian(ylim = c(-0.1, 11), expand = FALSE, clip = "off") +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Average Inner tongue area ", 
       x="Condition", 
       y="Inner tongue area [µm]",
       size=2.5, 
       fontface="bold")

# plot the Average IT size 
summaryMS<-summary(MS) # gives info of mean, median for each column for MS or Control (with PD & FTD patients)
summaryCON<-summary(COntrol)

# plot a distribution for the IT areas! 
histogram_IT<-ggplot(excluded_data, aes(x=InnerTongue_Area, fill=Condition))+
  geom_histogram(aes(y=..density..),alpha=0.8, position="dodge", binwidth=.8)+ # density means rel. frequency here (IMPORTANT)
  #geom_density(alpha=0.2)+
  #facet_grid(cols = vars(Condition))+
  scale_fill_manual(values = c('#999999','#E69F00')) +
  #geom_text(label="Kolmogorov-Smirnov, p=0.000193 ", 
            #y=-0.06, 
            #x=0.9, 
            #color="black", 
            #size=2.5, 
            #alpha=0.8)+ 
  coord_cartesian(xlim = c(-1, 11),ylim=c(0,1), expand = FALSE, clip = "off")+
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Frequency distribution of Inner tongue area",
       x="Inner tongue area [µm]", 
       y="Relative frequency")



# GRPAH 6: plot average g-ratio for each patient
# reorder data
newdata2<-newdata # reorders data into decreasing bars by copying original data and using reorder()
newdata2$Patient<-factor(newdata2$Patient, 
                         levels=newdata2$Patient[order(newdata1$Average_g_ratio,
                                                       decreasing=F)])
# plot this data 
average_g<-ggplot(newdata2, aes(Patient, Average_g_ratio, fill=Condition)) +
  geom_bar(stat="identity")+
  annotate(geom='text', 
           x=newdata2$Patient, 
           y=-0.05, 
           label=newdata2$Age,
           size=2.5, 
           alpha=0.8, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-0.05, 
            x=0.1, 
            color="black", 
            size=2.5, 
            alpha=0.8)+ 
  coord_cartesian(ylim = c(0, 0.7), expand = FALSE, clip = "off") +
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"),
        axis.title.x = element_text(vjust=-1.5))+
  labs(title="Average g-ratio per Patient (ordered increasing)", x="Patient", y="Average g-ratio")

# GRAPH 7: do the same but reorder in patient age (increasing)
newdata3<-newdata # reorders data into decreasing bars by copying original data and using reorder()
newdata3$Patient<-factor(newdata3$Patient, levels=newdata3$Patient[order(newdata3$Age, decreasing=F)])

average_g_age<-ggplot(newdata3, aes(Patient, Average_g_ratio, fill=Condition)) +
  geom_hline(yintercept=0, col="grey")+
  geom_bar(stat="identity")+
  annotate(geom='text', 
           x=newdata2$Patient,
           y=-0.06, 
           label=newdata2$Age, 
           size=2.5, 
           alpha=0.8,
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-0.06, 
            x=0.2, 
            color="black", 
            size=2.5, 
            alpha=0.8)+ 
  coord_cartesian(ylim = c(0, 0.7), expand = FALSE, clip = "off") +
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"),
        axis.title.x = element_text(vjust=-1.5))+
  labs(title="Average g-ratio per Patient ordered by increasing age", x="Patient", y="Average g-ratio")


# GRAPH 8: plot histogram for frequency distribution of axon caliber 
histogram<-ggplot(excluded_data, aes(x=Axon_caliber, fill=Condition))+
  geom_histogram(aes(y=..density..),alpha=0.8, position="dodge", binwidth=.5)+ # density means rel. frequency here (IMPORTANT)
  #geom_density(alpha=0.2)+
  #facet_grid(cols = vars(Condition))+
  scale_fill_manual(values = c('#999999','#E69F00')) +
  coord_cartesian(xlim = c(0, 8), ylim=c(0,1),expand = FALSE, clip = "off")+
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Frequency distribution of Axon calibers",
       x="Axon caliber [µm]", 
       y="Relative frequency")

# do kolmogorov-smirnov test here for the frequency distributions
MS<-subset(xdata, Condition=="MS")
CON<-subset(xdata, Condition=="CON")
ks.test(MS$Axon_caliber, CON$Axon_caliber) # p>0.05

#----------------------------------------------------------------------------


# PERIAXONAL SPACE ANALYSIS
PAS<-read.csv("Periaxonal_space_analysis.csv")

# plot PAS per patient 
# filter out NAs in Average_PS
PAS<-filter(PAS, !is.na(Average_PS))

desired_order<-c("19-31","19-26","20-60", "19-98", "19-27", "19-69", "19-47", "20-65", "20-64", "20-77","20-71", "19-106", "20-54", "19-36", "20-82","20-95", "20-92", "19-77", "20-58", "20-52")
PAS$Patient<-factor(as.character(PAS$Patient), levels=desired_order)

PAS_boxplot<-ggplot(PAS, aes(Patient, Average_PS, color=Condition))+
  stat_boxplot(aes(Patient, Average_PS),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Patient, Average_PS), outlier.shape=1)+
  #geom_violin(alpha = 0.5) +
  annotate(geom='text', 
           x=xdata$Patient, 
           y=-0.9, 
           label=xdata$Age,
           size=2.5, 
           alpha=0.2, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-0.9, 
            x=0.3, 
            color="black", 
            size=2.5, 
            alpha=0.2)+ 
  coord_cartesian(ylim = c(0.1, 22), expand = FALSE, clip = "off") +
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.1)+
  scale_color_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + #complete cols "#003264",'#999999','#E69F00','#719c70'
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Periaxonal space range per patient", 
       x="Patient", 
       y="Periaxonal space [nm]",
       size=2.5, 
       fontface="bold")

histogram_PAS_area<-ggplot(excluded_PAS, aes(x=Average_PS, fill=Condition))+
  geom_histogram(binwidth=1, position="dodge")+ # density means rel. frequency here (IMPORTANT)
  #geom_density(alpha=0.2)+
  #facet_grid(cols = vars(Condition))+
  scale_fill_manual(values = c('#999999','#E69F00')) +
  #coord_cartesian(xlim = c(0, 0.9), expand = FALSE, clip = "off")+
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Frequency distribution of Periaxonal space",
       x="Periaxonal space [nm]", 
       y="Relative frequency")

# con vs ms
compare_means(Average_PS~Condition, data=PAS, method="wilcox.test") # used to implement significance into the graph below
my_comparison<-list(c("CON","MS"))#,c("CON","FTD"), c("CON", "PD"))

PAS_area_boxplot_per_condition<-ggplot(excluded_PAS, aes(Condition, Average_PS, color=Condition, palette="jco"))+
  stat_boxplot(aes(Condition, Average_PS),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, Average_PS), outlier.shape=1)+
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) +
  stat_compare_means(method = "wilcox.test", label.y = 1.03)+
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Periaxonal space per condition", 
       x="Condition", 
       y="Peiaxonal space [nm]",
       size=2.5, 
       fontface="bold")

pairwise.wilcox.test(PAS$Average_PS, 
                     PAS$Condition, 
                     p.adjust="bonferroni") # significant

pairwise.wilcox.test(PAS$Average_PS,
                     PAS$Patient,
                     p.adjust="bonferroni")
PAS_wilcox_test <- as.data.frame(pairwise.wilcox.test(PAS$Average_PS,
                                                          PAS$Patient,
                                                          p.adjust="bonferroni")$p.value) # makes nice p-value table
# export p-value table 
write.csv(PAS_wilcox_test, "Pairwise_wilcoxon_PAS_new.csv")
