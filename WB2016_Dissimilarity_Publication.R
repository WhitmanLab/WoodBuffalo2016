# Loading R packages
library(reshape)
library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(plotly)
library(vegan)

---------------------------------
#Read in bacteria data phyloseq object
ps = readRDS(file = "ps.WB2016")
theme_set(theme_bw())
ps

#Remove mitochondria and chloroplast from ps; ensure the only domain is Bacteria
ps <- ps %>%
  subset_taxa(
    Domain == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast"
  )
ps

#Remove samples 21 (bcause we don't really know its identity), sample 10, sample 8, and samples that don't have a pair
SamplesToKeep = c('WB16-01-L-M-2','WB16-01-L-O-2','WB16-01-S-M-2','WB16-01-S-O-1','WB16-02-L-M-2',
                  'WB16-02-S-M-2','WB16-03-L-M-1','WB16-03-L-O-2','WB16-03-S-M-2','WB16-03-S-O-2','WB16-04-L-M-2',
                  'WB16-04-L-O-1','WB16-04-S-M-2','WB16-04-S-O-1','WB16-05-L-M-2','WB16-05-L-O-2','WB16-05-S-M-2',
                  'WB16-05-S-O-1','WB16-06-L-M','WB16-06-L-O','WB16-06-S-M-1','WB16-06-S-O-1','WB16-07-L-M',
                  'WB16-07-L-O','WB16-07-S-M','WB16-07-S-O','WB16-09-L-M','WB16-09-L-O','WB16-09-S-M','WB16-09-S-O',
                  'WB16-10-L-M','WB16-10-S-M','WB16-11-L-M','WB16-11-L-O','WB16-11-S-M','WB16-11-S-O',
                  'WB16-12-L-M','WB16-12-L-O','WB16-12-S-M','WB16-12-S-O','WB16-13-L-M','WB16-13-L-O','WB16-13-S-M',
                  'WB16-13-S-O','WB16-14-L-M','WB16-14-S-M','WB16-15-L-M','WB16-15-S-M',
                  'WB16-16-L-M','WB16-16-L-O','WB16-16-S-M','WB16-16-S-O','WB16-17-L-M','WB16-17-L-O','WB16-17-S-M',
                  'WB16-17-S-O','WB16-18-L-O','WB16-18-S-O','WB16-19-L-M','WB16-19-L-O','WB16-19-S-M','WB16-19-S-O',
                  'WB16-20-L-O','WB16-20-S-O','WB16-22-L-M','WB16-22-L-O','WB16-22-S-M','WB16-22-S-O','WB16-23-L-O',
                  'WB16-23-S-O','WB16-24-L-O','WB16-24-S-O')

#Take out blanks, multicommunity, and samples from site 8 (allburn, burn, unburn)
ps.pruned =  prune_samples(SamplesToKeep, ps)
ps.pruned

#Transform to even sampling depth
ps.norm = transform_sample_counts(ps.pruned, function(x) x/sum(x))

#Run NMDS ordination
ps.ord = ordinate(ps.norm, method = "NMDS", distance = "bray", trymax = 1000, k=2)
ps.ord

otu_table(ps.norm) = t(otu_table(ps.norm))
row.names(otu_table(ps.norm))
sample_data(ps.norm)[1:40] %>% select(TSLF_At_Sample)

#Pull out the S and L disimilarity values for each Pair and put into a .csv (WB2016_analysis_dataframe.csv)
dist = as.matrix(vegdist(otu_table(ps.norm), method="bray", ncol=length(sample_names(ps.norm))))
#dist["SampleName1","SampleName2"]
dist['WB16-01-L-M-2','WB16-01-S-M-2']
dist['WB16-01-L-O-2','WB16-01-S-O-1']
dist['WB16-02-L-M-2','WB16-02-S-M-2']
dist['WB16-03-L-O-2','WB16-03-S-O-2']
dist['WB16-03-L-M-1','WB16-03-S-M-2']
dist['WB16-04-L-M-2','WB16-04-S-M-2']
dist['WB16-04-L-O-1','WB16-04-S-O-1']
dist['WB16-05-L-M-2','WB16-05-S-M-2']
dist['WB16-05-L-O-2','WB16-05-S-O-1']
dist['WB16-06-L-M','WB16-06-S-M-1']
dist['WB16-06-L-O','WB16-06-S-O-1']
dist['WB16-07-L-M','WB16-07-S-M']
dist['WB16-07-L-O','WB16-07-S-O']
dist['WB16-09-L-M','WB16-09-S-M']
dist['WB16-09-L-O','WB16-09-S-O']
dist['WB16-10-L-M','WB16-10-S-M']
dist['WB16-11-L-M','WB16-11-S-M']
dist['WB16-11-L-O','WB16-11-S-O']
dist['WB16-12-L-M','WB16-12-S-M']
dist['WB16-12-L-O','WB16-12-S-O']
dist['WB16-13-L-M','WB16-13-S-M']
dist['WB16-13-L-O','WB16-13-S-O']
dist['WB16-14-L-M','WB16-14-S-M']
dist['WB16-15-L-M','WB16-15-S-M']
dist['WB16-16-L-M','WB16-16-S-M']
dist['WB16-16-L-O','WB16-16-S-O']
dist['WB16-17-L-M','WB16-17-S-M']
dist['WB16-17-L-O','WB16-17-S-O']
dist['WB16-18-L-O','WB16-18-S-O']
dist['WB16-19-L-M','WB16-19-S-M']
dist['WB16-19-L-O','WB16-19-S-O']
dist['WB16-20-L-O','WB16-20-S-O']
dist['WB16-22-L-M','WB16-22-S-M']
dist['WB16-22-L-O','WB16-22-S-O']
dist['WB16-23-L-O','WB16-23-S-O']
dist['WB16-24-L-O','WB16-24-S-O']

#Make a new column combining the features for soil layer and drainage class
sample_data(ps.norm)$Layer_LandClass = paste(sample_data(ps.norm)$Soil_layer, sample_data(ps.norm)$Land_Class, sep = "-")
sample_data(ps.norm)$Layer_LandClass
sample_data(ps.norm)$Layer_LandClass = levels = c('M-Upland', 'M-Wetland','O-Upland', 'O-Wetland', ordered = TRUE)

# Read in data tables: 
df = read.csv('WB2016_analysis_dataframe.csv')
df

#Add column for Land Class, then make another column for land class and soil horizon to be able to map all four variables onto plot
df$LandClass = ifelse(df$Wetland == "no", "Upland", "Wetland")
df$Layer_LandClass = paste(df$Soil_Layer, df$LandClass, sep = "-")
df$Layer_LandClass = levels = c('M-Upland', 'M-Wetland','O-Upland', 'O-Wetland', ordered = TRUE)
df



# Figure 1: Non-metric multidimensional scaling plot of Bray-Curtis dissimilarities of the bacterial community composition 
# Black outline indicating low-intensity surface fires were added to figure manually.
p = plot_ordination(ps.norm, ps.ord, color = "Interval", shape = "Layer_LandClass", axes=c(1,2))
p = p + geom_point(size=5)
p = p + theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=18, face="bold"),
        axis.title = element_text(size=22, face="bold"),
        axis.text = element_text(size=18, face = "bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size = 22),
        strip.background = element_rect(colour="white", fill="white"))
p = p + theme(plot.title=element_text(size = 25))
p = p + theme(plot.subtitle=element_text(size = 18))
p = p + scale_colour_manual(name = "Interval", values = c("#009E73", "#D55E00","#999999","#E69F00"))
p = p + scale_shape_manual(name = "Soil Layer and \n Drainage Class", values=c(19, 21, 15, 22))
p

# Adding arrows joining SI and LI to the plot
x = data.frame(ps.ord$points)$MDS1
y = data.frame(ps.ord$points)$MDS2

sample_data(ps.norm)$x = x
sample_data(ps.norm)$y = y

# These are the co-ordinates of each point.
# We need to know how to join them.
colnames(sample_data(ps.norm))
sample_data(ps.norm)$Pair_ID
sample_data(ps.norm)$Soil_layer
sample_data(ps.norm)$Pairs = paste(sample_data(ps.norm)$Pair_ID,sample_data(ps.norm)$Soil_layer,sep="_")

# We will use geom_segment, which takes x, y, and xend, yend.
Segments = data.frame(sample_data(ps.norm)) %>%
  dplyr::select(Pairs,Pair_ID,Interval,Soil_layer,Layer_LandClass,Land_Class,x,y)

Segments = melt(Segments,id=c("Pairs","Pair_ID","Interval","Soil_layer","Layer_LandClass","Land_Class"))

# Using reshape2 package and dcast function
Segments = reshape2::dcast(Segments,Pairs+Pair_ID+Soil_layer+Layer_LandClass+Land_Class~Interval+variable,value.var = "value",fun=mean)

# Creating the same plot as above with the two years of data, linked by lines.
p = plot_ordination(ps.norm, ps.ord, color = "Interval", shape = "Layer_LandClass", axes=c(1,2))
p = p + geom_point(size=5)
p = p + theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=18, face="bold"),
        axis.title = element_text(size=22, face="bold"),
        axis.text = element_text(size=18, face = "bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size = 22),
        strip.background = element_rect(colour="white", fill="white"))
p = p + theme(plot.title=element_text(size = 25))
p = p + theme(plot.subtitle=element_text(size = 18))
p = p + scale_color_manual(name = "Interval", values = c("#009E73", "#D55E00","#999999","#E69F00"))
p = p + scale_shape_manual(name = "Soil Horizon and \n Drainage Class", values=c(19, 21, 15, 22))
p = p + geom_segment(data=Segments,aes(x=Long_x,y=Long_y,xend=Short_x,yend=Short_y),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="darkgrey",alpha=0.5)
p






## Table 2: ANOVA on stem count effect on bacterial dissimilarity

# Understory Stem density understory
anova.stem = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_understory_stems, data=df)
summary(anova.stem)

lm.stem = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_understory_stems, data=df)
anova(lm.stem)

summary(lm.stem)

#Stem density understory conifer only
anova.stem = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + diff_underCon_stems, data=df)
summary(anova.stem)

lm.stem = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + diff_underCon_stems, data=df)
anova(lm.stem)

summary(lm.stem)


#Stem density understory broadleaf only
anova.stem = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_underBL_stems, data=df)
summary(anova.stem)

lm.stem = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_underBL_stems, data=df)
anova(lm.stem)

summary(lm.stem)


## Figure 3: log of the difference (LI - SI) in understory stem count, plotted by all understory stems, understory conifer, and understory broadleaf

#3a: Difference in log understory stem count and dissimilarity
p = ggplot(df, aes(Diff_understory_stems_log, BrayCurtis_Dissimilarity, color = LandClass, shape = Soil_Layer))
p = p + geom_point(size=7)
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
#p = p + labs(title="Dissimilarity of Bacteria Community in Paired Sites by Difference in Number of Stems")
p = p + ylab("Bray-Curtis Dissimilarity")
p = p + xlab("Difference in Number of Understory Stems log(Long - Short)")
p = p + scale_colour_manual(name = "Land Class",values = c("#000000", "#56B4E9", "#009E73"))
p = p + scale_shape_manual(name = "Soil Layer", values=c(19, 15))
p = p + theme(plot.title=element_text(size = 25))
p


#3b: Difference in log Conifer understory stem count and dissimilarity
p = ggplot(df, aes(diff_underCon_stems_log, BrayCurtis_Dissimilarity, color = LandClass, shape = Soil_Layer))
p = p + geom_point(size = 7)
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
p = p + ylab("Bray-Curtis Dissimilarity")
p = p + xlab("Difference in Number of Understory Conifer Stems log(Long - Short)")
p = p + scale_colour_manual(name = "Land Class",values = c("#000000", "#56B4E9", "#009E73"))
p = p + scale_shape_manual(name = "Soil Layer", values=c(19, 15))
p = p + theme(plot.title=element_text(size = 25))
p


#3c: Difference in log broadleaf understory stem count and dissimilarity
p = ggplot(df, aes(Diff_underBL_stems_log, BrayCurtis_Dissimilarity, color = LandClass, shape = Soil_Layer))
p = p + geom_point(size = 7)
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
p = p + ylab("Bray-Curtis Dissimilarity")
p = p + xlab("Difference in Number of Understory Broadleaf Stems log(Long - Short)")
p = p + scale_colour_manual(name = "Land Class",values = c("#000000", "#56B4E9", "#009E73"))
p = p + scale_shape_manual(name = "Soil Layer", values=c(19, 15))
p




### Code for Supplemental Information

## Calulate the dissimlarity distances for the vegetation community
vegcom = read.csv("WBNPNWT_Vascular_Vegetation_Community_2016.csv", row.names = c("16-NT-11l", "16-NT-11s", "16-NT-12l", "16-NT-12s", "16-NT-13l", "16-NT-13s", "16-NT-16l", "16-NT-16s", "16-NT-17l", "16-NT-17s",
                                                                                  "16-NT-18l", "16-NT-18s", "16-NT-19l", "16-NT-19s", "16-NT-20l",  "16-NT-20s", "16-NT-22l",  "16-NT-22s",  "16-NT-23l", "16-NT-23s",
                                                                                  "16-NT-24l",  "16-NT-24s", "16-WB-01l", "16-WB-01s", "16-WB-02l", "16-WB-02s",  "16-WB-03l", "16-WB-03s", "16-WB-04l", "16-WB-04s",
                                                                                  "16-WB-05l",  "16-WB-05s",  "16-WB-06l", "16-WB-06s", "16-WB-07l", "16-WB-07s","16-WB-09l",  "16-WB-09s", "16-WB-10lA", "16-WB-10sA",
                                                                                  "16-WB-14l", "16-WB-14s", "16-WB-15l", "16-WB-15s", "16-WB-10lB", "16-WB-10sB", "16-NT-21l", "16-NT-21s", "16-WB-08l", "16-WB-08s"))
#want row.names to be sample ID, and col names to be plant ID
vegcom = vegcom[-c(1)]
ord.veg = metaMDS(vegcom)
ord.veg
#Pull out the S and L disimilarity values for each Pair
dist = as.matrix(vegdist(vegcom, method="bray", ncol=length(vegcom)))
dist["16-NT-11l", "16-NT-11s"]
dist["16-NT-12l", "16-NT-12s"]
dist["16-NT-13l", "16-NT-13s"]
dist["16-NT-16l", "16-NT-16s"]
dist["16-NT-17l", "16-NT-17s"]
dist["16-NT-18l", "16-NT-18s"]
dist["16-NT-19l", "16-NT-19s"]
dist["16-NT-20l",  "16-NT-20s"]
dist["16-NT-22l",  "16-NT-22s"]
dist["16-NT-23l", "16-NT-23s"]
dist["16-NT-24l",  "16-NT-24s"]
dist["16-WB-01l", "16-WB-01s"]
dist["16-WB-02l", "16-WB-02s"]
dist["16-WB-03l", "16-WB-03s"]
dist["16-WB-04l", "16-WB-04s"]
dist["16-WB-05l",  "16-WB-05s"]
dist["16-WB-06l", "16-WB-06s"]
dist["16-WB-07l", "16-WB-07s"]
dist["16-WB-09l",  "16-WB-09s"]
dist["16-WB-10lA", "16-WB-10sA"]
dist["16-WB-14l", "16-WB-14s"]
dist["16-WB-15l", "16-WB-15s"]

#These values were manually added into analysis csv


#Make datatable from ps.pruned to run PERMANOVA on
ps.df = as(sample_data(ps.pruned), "data.frame")
d = distance(ps.pruned, method = "bray")

##SI Table 1: PERMANOVA on Pair ID and Interval
d.adonis = adonis(d ~ Pair_ID + Interval, ps.df)
d.adonis

##SI Table 2: PERMANOVA on Drainage class, soil layer, and pH
d.adonis = adonis(d ~ Land_Class + Soil_layer + pH, ps.df)
d.adonis

##SI Table 3: ANOVA results for effect of TSLF, vegetation transition, land class, and difference FFI in individual models for their effects on paired site Bray-Curtis dissimilarities

#TSLF
anova.TSLF = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + TSLF_At_Sample, data=df)
summary(anova.TSLF)

lm.TSLF = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + TSLF_At_Sample, data=df)
anova(lm.TSLF)

summary(lm.TSLF)

#Vegetation Transition 
anova.Veg = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Veg_Transition, data=df)
summary(anova.Veg)

lm.veg = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Veg_Transition, data=df)
anova(lm.veg)

summary(lm.veg)

#Drainage Class 
anova.wet = aov(BrayCurtis_Dissimilarity ~ Soil_Layer + Wetland, data=df)
summary(anova.wet)

lm.wet = lm(BrayCurtis_Dissimilarity ~ Soil_Layer + Wetland, data=df)
anova(lm.wet)

summary(lm.wet)

head(df)

#Diff in FFI
anova.FFI = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_in_Recent_Interval, data=df)
summary(anova.FFI)

lm.FFI = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_in_Recent_Interval, data=df)
anova(lm.FFI)

summary(lm.FFI)

## SI Table 4: ANOVA results for effect of difference in stem count, difference in vegetation cover, and vegetation community dissimilarities on bacterial dissimilarities

#Vegetation Cover
anova.cov = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_veg_cover, data=df)
summary(anova.cov)

lm.cov = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_veg_cover, data=df)
anova(lm.cov)

summary(lm.cov)

#Stem count
df
anova.stem = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_total_stems, data=df)
summary(anova.stem)

lm.stem = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Diff_total_stems, data=df)
anova(lm.stem)

summary(lm.stem)

#Vegetation dissimilarity
anova.dis = aov(BrayCurtis_Dissimilarity ~ Layer_LandClass + Veg_Dissimilarity, data=df)
summary(anova.dis)

lm.dis = lm(BrayCurtis_Dissimilarity ~ Layer_LandClass + Veg_Dissimilarity, data=df)
anova(lm.dis)

summary(lm.dis)


#SI Figure 1: TSLF and bacterial dissimilarity
options(repr.plot.width=60, repr.plot.height=30)
p = ggplot(df, aes(TSLF_At_Sample, BrayCurtis_Dissimilarity, color = Veg_Transition, shape = Layer_LandClass))
p = p + geom_point(size = 7)
p = p + scale_shape_manual(name = "Soil Horizon and \nLand Class", values=c(19, 21, 15, 22))
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
#p = p + labs(title="Dissimilarity of Bacteria Community in Paired Sites by Time Since Last Fire")
p = p + ylab("Bray-Curtis Dissimilarity")
p = p + xlab("Time Since Last Fire (TSLF) at Date of Sampling")
p = p + scale_colour_manual(name = "Vegetation Transition",values = c("#56B4E9", "#000000"))
p = p + theme(plot.title=element_text(size = 25))
p

#SI Figure 2: Diff in recent interval and bacterial dissimilarity
options(repr.plot.width=60, repr.plot.height=30)
p = ggplot(df, aes(Diff_in_Recent_Interval, BrayCurtis_Dissimilarity, color = Veg_Transition, shape = Layer_LandClass))
p = p + geom_point(size = 7)
p = p + scale_shape_manual(name = "Soil Horizon and \nLand Class", values=c(19, 21, 15, 22))
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
#p = p + labs(title="Dissimilarity of Bacteria Community in Paired Sites by Differences in Fire Return Interval ")
p = p + ylab("Bray-Curtis Dissimilarity")
p = p + xlab("Difference in Fire Free Interval (FFI)")
p = p + scale_colour_manual(name = "Change in Leading Tree Species",values = c("#D55E00","#999999", "#009E73","#E69F00"))
p = p + theme(plot.title=element_text(size = 25))
p

#SI Figure 3: Vegetation cover and bacterial dissimilarity

p = ggplot(df, aes(Diff_veg_cover , BrayCurtis_Dissimilarity, color = Veg_Transition, shape = Layer_LandClass))
p = p + geom_point(size = 7)
p = p + scale_shape_manual(name = "Soil Horizon and \nLand Class", values=c(19, 21, 15, 22))
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
#p = p + labs(title="Dissimilarity of Bacteria Community in Paired Sites by Difference in Number of Stems")
p = p + ylab("Bray-Curtis Dissimilarity")
p = p + xlab("Difference in Vegetation Cover")
p = p + scale_colour_manual(name = "Vegetation Transition",values = c("#56B4E9", "#000000"))
p = p + theme(plot.title=element_text(size = 25))
p

# SI Figure 4: Vegetation dissimilarity and Bacterial dissimilarity

options(repr.plot.width=60, repr.plot.height=30)
p = ggplot(df, aes(Veg_Dissimilarity , BrayCurtis_Dissimilarity, color = Veg_Transition, shape = Layer_LandClass))
p = p + geom_point(size = 7)
p = p + scale_shape_manual(name = "Soil Horizon and \nLand Class", values=c(19, 21, 15, 22))
p = p + theme_bw() + 
  theme(axis.title.x = element_text(size=25, face="bold"),
        axis.title.y = element_text(size=25, face="bold"),
        axis.text = element_text(size=20, face = "bold", color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 22))
#p = p + labs(title="Dissimilarity of Bacteria Community in Paired Sites by Difference in Number of Stems")
p = p + ylab("Bray-Curtis Dissimilarity of (Bacteria)")
p = p + xlab("Bray-Curtis Dissimilarity of (Vegetation)")
p = p + scale_colour_manual(name = "Vegetation Transition",values = c("#56B4E9", "#000000"))
p = p + theme(plot.title=element_text(size = 25))
p

