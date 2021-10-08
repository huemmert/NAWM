#####                           #####
###     Mitochondria Analysis    ###
#####                          #####  

library(dplyr)
library(ggplot2)
library(pgirmess)
library(ggpubr)

# working directory
setwd("F:/TEM ANALYSIS R/data")

# load data
View(Mito_Analysis)
  xdata<-read.csv("Mito_Analysis_norm_27072021.csv") # to analyse the counts
ydata<-read.csv("Mito_Analysis_norm_27072021.csv") # to analyse the Mito areas

# filter out NAs in Myelinated_Axon column & Outfolding column
xdata<-filter(xdata, !is.na(Myelinated_Axons))
ydata<-filter(ydata, !is.na(Mito_Percent_Axon))

# exclude MS-1,PD-1, FTD-1
ydata_excluded<-subset(ydata, Patient!="MS-1" & Patient!="FTD-1" & Patient!="PD-1")
xdata_excluded<-subset(xdata, Patient!="MS-1" & Patient!="FTD-1" & Patient!="PD-1")

#----------------------------------------------------
# Analyse the mito areas per patient and per condition
#

# make FREQUENCY distributions of mito area distributions
histogram_Mito_area<-ggplot(ydata_excluded, aes(x=Mito_Area, fill=Condition))+
  geom_histogram(aes(y=..density..),alpha=0.8, position="dodge", binwidth=0.2)+ # density means rel. frequency here (IMPORTANT)
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Frequency distribution of Mitochondria profile areas",
       x="Mitochondria area [µm²]", 
       y="Relative frequency")


newdata<-ydata_excluded # copy data so we can reorder the ages manually 
desired_order<-c("MS-1","PD-1","MS-2", "MS-3", "MS-4", "MS-5", "MS-6", "FTD-1", "MS-7", "MS-8","CON-1", "CON-2", "CON-3", "CON-4", "CON-5","MS-9", "CON-6", "CON-7", "CON-8", "CON-9")
newdata$Patient<-factor(as.character(newdata$Patient), levels=desired_order)

# plot mito areas boxplot per patient 
Mito_area_boxplot<-ggplot(ydata, aes(Patient, Mito_Percent_Axon, color=Condition))+
  stat_boxplot(aes(Patient, Mito_Percent_Axon),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Patient, Mito_Percent_Axon), outlier.shape=NA)+
  annotate(geom='text', 
           x=newdata$Patient, 
           y=-3.5, 
           label=newdata$Age,
           size=2.5, 
           alpha=0.2, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-3.5, 
            x=0.06, 
            color="black", 
            size=2.5, 
            alpha=0.2)+ 
  coord_cartesian(ylim = c(0, 100), expand = FALSE, clip = "off") +
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.25)+
  scale_color_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + #'#999999','#719c70','#E69F00',"#003264"complete cols "#003264",'#999999','#E69F00','#719c70'
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Individual Mitochondrial area distributions", 
       x="Patient", 
       y="Axon area covered by Mitochondria [%]",
       size=2.5, 
       fontface="bold")

# analyse conditions con vs. ms only 
compare_means(Mito_Percent_Axon~Condition, data=newdata, method="kruskal.test") # used to implement significance into the graph below
my_comparison<-list(c("CON","MS"))#,c("CON","FTD"), c("CON", "PD"))

mito_area_boxplot_per_condition<-ggplot(newdata, aes(Condition, Mito_Percent_Axon, color=Condition, palette="jco"))+
  stat_boxplot(aes(Condition, Mito_Percent_Axon),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, Mito_Percent_Axon), outlier.shape=1)+
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Mitochondria areas per condition", 
       x="Condition", 
       y="Axon area covered by Mitochondria [%]",
       size=2.5, 
       fontface="bold")

# which test


#-------------------------------------------------------------------------------
# Analyse the outfolding count per patient and condition
#

# reorder the data 
newdata2<-xdata
desired_order<-c("MS-1","PD-1","MS-2", "MS-3", "MS-4", "MS-5", "MS-6", "FTD-1", "MS-7", "MS-8","CON-1", "CON-2", "CON-3", "CON-4", "CON-5","MS-9", "CON-6", "CON-7", "CON-8", "CON-9")
newdata2$Patient<-factor(as.character(newdata2$Patient), levels=desired_order)

Outfoldings_boxplot<-ggplot(newdata2, aes(Patient, Outfoldings_norm, color=Condition))+
  stat_boxplot(aes(Patient, Outfoldings_norm),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Patient, Outfoldings_norm), outlier.shape=NA)+
  annotate(geom='text', 
           x=newdata2$Patient, 
           y=-0.0007, 
           label=newdata2$Age,
           size=2.5, 
           alpha=0.2, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-0.0007, 
            x=0.06, 
            color="black", 
            size=2.5, 
            alpha=0.2)+ 
  coord_cartesian(ylim = c(0, 0.02), expand = FALSE, clip = "off") +
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.25)+
  scale_color_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + #complete cols "#003264",'#999999','#E69F00','#719c70'
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Outfoldings per Patient", 
       x="Patient", 
       y="Outfoldings per 1 µm²",
       size=2.5, 
       fontface="bold")

#-- analyse outfoldings per condition
compare_means(Outfoldings~Condition, data=newdata2, method="kruskal.test") # used to implement significance into the graph below
my_comparison<-list(c("CON","MS"))#,c("CON","FTD"), c("CON", "PD"))

newdata2<-xdata_excluded

Outfoldings_boxplot_per_condition<-ggplot(newdata2, aes(Condition, Outfoldings_norm, color=Condition, palette="jco"))+
  stat_boxplot(aes(Condition, Outfoldings_norm),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, Outfoldings_norm), outlier.shape=1)+
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Outfoldings per condition", 
       x="Condition", 
       y="Outfoldings per 1 µm²",
       size=2.5,
       fontface="bold")

# frequency distribution of outfoldings
outfolding_distribution<-ggplot(xdata_excluded, aes(x=Outfoldings_norm, fill=Condition))+
  geom_histogram(aes(y = stat(count) / sum(count)),alpha=0.8, position="dodge", binwidth=0.002)+# density means rel. frequency here (IMPORTANT)
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Frequency distribution of Outfoldings",
       x="Outfoldings per 1 µm²", 
       y="Relative frequency")


#---------------------------------------------------------------
# do the same for myelinoid bodies
# boxplot for myelinoid bodies per patient 
newdata2<-xdata
desired_order<-c("MS-1","PD-1","MS-2", "MS-3", "MS-4", "MS-5", "MS-6", "FTD-1", "MS-7", "MS-8","CON-1", "CON-2", "CON-3", "CON-4", "CON-5","MS-9", "CON-6", "CON-7", "CON-8", "CON-9")
newdata2$Patient<-factor(as.character(newdata2$Patient), levels=desired_order)

myelinoid_bodies_boxplot<-ggplot(newdata2, aes(Patient, Myelinoid_bodies_norm, color=Condition))+
  stat_boxplot(aes(Patient, Myelinoid_bodies_norm),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Patient, Myelinoid_bodies_norm), outlier.shape=NA)+
  annotate(geom='text', 
           x=newdata2$Patient, 
           y=-0.001, 
           label=newdata2$Age,
           size=2.5, 
           alpha=0.2, 
           fontface="plain")+
  geom_text(label="Patient age", # adds the patient age row below x-axis
            y=-0.001, 
            x=0.001, 
            color="black", 
            size=2.5, 
            alpha=0.2)+ 
  coord_cartesian(ylim = c(0, 0.03), expand = FALSE, clip = "off") +
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.25)+
  scale_color_manual(values = c('#999999','#719c70','#E69F00',"#003264")) + #complete cols "#003264",'#999999','#E69F00','#719c70'
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Myelinoid bodies per Patient", 
       x="Patient", 
       y="Myelinoid bodies per 1 µm²",
       size=2.5, 
       fontface="bold")

compare_means(Myelinoid_bodies~Condition, data=newdata2, method="kruskal.test") # used to implement significance into the graph below
my_comparison<-list(c("CON","MS"))#,c("CON","FTD"), c("CON", "PD"))

# boxplot per condition
newdata2<-xdata_excluded
desired_order<-c("MS-1","PD-1","MS-2", "MS-3", "MS-4", "MS-5", "MS-6", "FTD-1", "MS-7", "MS-8","CON-1", "CON-2", "CON-3", "CON-4", "CON-5","MS-9", "CON-6", "CON-7", "CON-8", "CON-9")
newdata2$Patient<-factor(as.character(newdata2$Patient), levels=desired_order)
newdata2<-subset(newdata2, Patient!="MS-1" & Patient!="PD-1" & Patient!="FTD-1")

myelinoid_bodies_boxplot_per_condition<-ggplot(newdata2, aes(Condition, Myelinoid_bodies_norm, color=Condition, palette="jco"))+
  stat_boxplot(aes(Condition, Myelinoid_bodies_norm),
               geom="errorbar", 
               linetype=1, 
               width=0.5)+
  geom_boxplot(aes(Condition, Myelinoid_bodies_norm), outlier.shape=1)+
  stat_summary(fun=mean, geom="point", size=2)+
  geom_jitter(position=position_jitter(0.2), alpha=0.2)+
  scale_color_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5, vjust=-0.2), 
        axis.title.x = element_text(vjust=-1), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Myelinoid bodies per condition", 
       x="Condition", 
       y="Myelinoid bodies per 1 µm²",
       size=2.5,
       fontface="bold")

xdata_excluded<-filter(xdata_excluded, !is.na(Myelinoid_bodies_norm))
# frequency distribution of myelinoid bodies
myelinoid_bodies_distribution<-ggplot(xdata_excluded, aes(x=Myelinoid_bodies_norm, fill=Condition))+
  geom_histogram(aes(y = stat(count) / sum(count)),alpha=0.8, position="dodge", binwidth=0.003)+  # density means rel. frequency here (IMPORTANT)
  scale_fill_manual(values = c('#999999','#E69F00')) +
  theme_bw()+
  theme(panel.border=element_rect(colour="black"), 
        plot.title = element_text(hjust=0.5), 
        axis.ticks = element_line(colour="black"))+
  labs(title="Frequency distribution of Myelinoid bodies",
       x="Myelinoid bodies per 1 µm²", 
       y="Relative frequency")


######              ######
###     STATISTICS     ###
######              ######

MS_x<-subset(xdata, Condition=="MS")
MS_y<-subset(ydata_excluded, Condition=="MS")
CON_x<-subset(xdata, Condition=="CON")
CON_y<-subset(ydata_excluded, Condition=="CON")
shapiro.test(MS_y$Mito_Percent_Axon) # p<0.05, non-parametric
shapiro.test(CON_y$Mito_Percent_Axon) #p<0.05
qqnorm(MS_y$Mito_Percent_Axon) # non-normal
qqnorm(CON_y$Mito_Percent_Axon)

shapiro.test(MS_x$Myelinoid_bodies) # p>0.05
shapiro.test(CON_x$Myelinoid_bodies) # p<0.05

shapiro.test(MS_x$Outfoldings) # p<0.05
shapiro.test(CON_x$Outfoldings) # p<0.05

shapiro.test(MS_x$Myelinoid_bodies) # p>0.05
shapiro.test(CON_x$Myelinoid_bodies) # p<0.05

# test if mito areas are significantly different to each other 
kruskal.test(Mito_Percent_Axon~Condition, data=ydata) # p<0.05, significant

# post-hoc
pairwise.wilcox.test(ydata$Mito_Area,
                     ydata$Patient,
                     p.adjust="BH")
nito_area_wilcox_test <- as.data.frame(pairwise.wilcox.test(ydata$Mito_Area,
                                                          ydata$Patient,
                                                          p.adjust="BH")$p.value)
write.csv(nito_area_wilcox_test, "mito_area_wilcoxon.csv")

#---------------------
# RE-DO Statistics and exclude MS-1, PD-1, FTD-1

#-----------------------------------------------------------
# TEST IF EXCLUDE 19-31 AND FTD/PD
excluded_ydata<-subset(ydata, Patient!="MS-1"& Patient!="PD-1"&Patient!="FTD-1")
MS_y_less<-subset(excluded_ydata, Condition=="MS")
CON_y_less<-subset(excluded_ydata, Condition=="CON")
qqnorm(CON_y_less$Mito_Percent_Axon)
qqnorm(MS_y_less$Mito_Percent_Axon)

shapiro.test(MS_y_less$Mito_Percent_Axon) # non parametric
shapiro.test(CON_y_less$Mito_Percent_Axon) # non parametric

wilcox.test(Mito_Percent_Axon~Condition, data=excluded_ydata) # p=0.0753

# look at outfoldings
qqnorm(CON_y_less$Outfoldings) # not normal
qqnorm(MS_y_less$Outfoldings)
shapiro.test(CON_y_less$Outfoldings) # p<0.05
shapiro.test(MS_y_less$Outfoldings) # p<0.05
wilcox.test(Outfoldings~Condition, data=excluded_ydata) # p=0.1016

# look at myelinoid bodies 
qqnorm(CON_y_less$Myelinoid_bodies)
qqnorm(MS_y_less$Myelinoid_bodies)
shapiro.test(CON_y_less$Myelinoid_bodies) # p<0.05
shapiro.test(MS_y_less$Myelinoid_bodies) # p>0.05
wilcox.test(Myelinoid_bodies~Condition, data=excluded_ydata)
