# Loading R packages
library(reshape)
library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(plotly)
library(vegan)

#Read in phyloseq object (OTU table, metadata, and taxonomy table)
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

#Create a list of all the samples to keep; remove sample 21, 10B, 8, blanks, and samples that don't have a pair (i.e there was organic soil sample in LI, but not organic in SI)
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

# Prune samples that were not included in above list
ps.pruned =  prune_samples(SamplesToKeep, ps)
ps.pruned

#Transform to even sampling depth
ps.norm = transform_sample_counts(ps.pruned, function(x) x/sum(x))

#Make a new column in sample data that indicates the "class" of the leading tree species: Conifer, Broadleaf, or Mixed

Conifer = c('PICEMAR', 'PINUBAN','PICEGLA','PICEMAR/LARILAR','LARILAR' )
Broadleaf = c('POPUTRE','POPUBAL', 'POPUBAL/POPUTRE')
Mixed = c('PICEMAR/POPUTRE', 'PINUBAN/POPUTRE','POPUTRE/PINUBAN')


sample_data(ps.norm)$Leading_Type = ifelse(sample_data(ps.norm)$Leading_Species %in% Conifer, "Conifer",
                                           ifelse(sample_data(ps.norm)$Leading_Species %in% Broadleaf, "Broadleaf",
                                                  ifelse(sample_data(ps.norm)$Leading_Species %in% Mixed, "Mixed", NA)))

#Make a new column in sample data for layer and landclass so that each land class and soil layer can be plotted on same plot
sample_data(ps.norm)$Layer_LandClass = paste(sample_data(ps.norm)$Soil_layer, sample_data(ps.norm)$Land_Class, sep = "-")

# Reorder the landclass levels for graphing aesthetics
sample_data(ps.norm)$Layer_LandClass = ordered(sample_data(ps.norm)$Layer_LandClass, levels = c('M-Upland', 'M-Wetland','O-Upland', 'O-Wetland', ordered = TRUE))


#Run NMDS ordination
ps.ord = ordinate(ps.norm, method = "NMDS", distance = "bray", trymax = 1000, k=4)

## Look at ordination plots for the full community
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
p = p + scale_shape_manual(name = "Soil Layer and \n Land Class", values=c(19, 21, 15, 22))
p = p + geom_text(aes(label=Pair_ID),hjust=2, vjust=.5, size = 4, color = "black")
p

########### Addding arrows to the ordination plot ################
############ Arrows for factors related to samples

# Check sample info
colnames(sample_data(ps.norm))

# Which data are non-numeric?
charfactors = c("Moisture_Regime","Leading_Species","Land_Class","Interval","Indic_Class","Ecosite_Phase","Ecosite_Name","fwd_barcode","rev_barcode","Leading_Type","Layer_LandClass","Site_ID","Soil_layer")

# Prep sample data for envfit
samdat.envfit = data.frame(sample_data(ps.norm))

# Get only numerics and make sure they are formatted not as categorical
samdat.envfit.numeric = samdat.envfit[,-which(colnames(samdat.envfit) %in% charfactors)]
samdat.envfit.numeric = apply(samdat.envfit.numeric,2,paste)
samdat.envfit.numeric = apply(samdat.envfit.numeric,2,as.numeric)
samdat.envfit.numeric = data.frame(samdat.envfit.numeric)

# Get just the character columns
samdat.envfit.char = samdat.envfit[,which(colnames(samdat.envfit) %in% charfactors)]

# Join the numeric and characters together.
samdat.envfit[,-which(colnames(samdat.envfit) %in% charfactors)]=samdat.envfit.numeric
samdat.envfit[,which(colnames(samdat.envfit) %in% charfactors)]=samdat.envfit.char
head(samdat.envfit)

# If we want to cut the wordy variables
# "Moisture_Regime"
colnames(samdat.envfit)
Col.exclude = c("Site_ID","Indic_Class","rev_barcode","fwd_barcode","Leading_Species",
                "Ecosite_Phase","Ecosite_Name","Latitude","Longitude",
                "CMD_Ann_19912017","CMD_Y","CMD_Y.1","CMD_Y.1_Mean_Departure","CMD_Y1",                                 
                "CMD_Y1_Mean_Departure","CMD_Y2","CMD_Y2_Mean_Departure","CMD_Y_Mean_Departure",
                "Recent_Interval","Recent_Year","Old_Year","TSLF_At_Sample",
                "Pct_Conifer_Biomass","Understory_Conifer_Biomass_kg_Ha_Live",
                "Broadleaf_Biomass","Conifer_Biomass","Load_Rotten_nologs",
                "Load_Sound_nologs","Moisture_Regime","No_Trees","Old_Sites",
                "Overstory_Broadleaf_Biomass_kg_Ha_Live","Overstory_Conifer_Biomass_kg_Ha_Live",
                "Overstory_Total_Biomass_kg_Ha_Live","Pair_ID","Total_Biomass","Total_CWD_Load_Kg_m2",
                "Total_Load_Rotten","Understocked","Understory_Broadleaf_Biomass_kg_Ha_Live",
                "Understory_Total_Biomass_kg_Ha_Live","jitterx","jittery","Layer_LandClass")
samdat.envfit=samdat.envfit[,-which(colnames(samdat.envfit) %in% Col.exclude)]

# Run the envfit for axes 1 and 2
ps.ord.envfit.1.2 = envfit(ps.ord,samdat.envfit,choices=c(1,2),na.rm=TRUE)
# and for 3 and 4
ps.ord.envfit.3.4 = envfit(ps.ord,samdat.envfit,choices=c(3,4),na.rm=TRUE)
# and for 2 and 4
ps.ord.envfit.2.4 = envfit(ps.ord,samdat.envfit,choices=c(2,4),na.rm=TRUE)

# Working with the 2nd and 4th axes
ps.ord.envfit = ps.ord.envfit.2.4

# Grab the scores for each factor, where p<cutoff
cutoff = 0.05
ps.ord.envfit$vectors$pvals
ps.ord.envfit.scores.vectors = data.frame(scores(ps.ord.envfit,"vectors")[ps.ord.envfit$vectors$pvals<cutoff,])
ps.ord.envfit.scores.vectors$Factor = "Numeric"
ps.ord.envfit.scores.factors = data.frame(scores(ps.ord.envfit,"factors")[ps.ord.envfit$factors$pvals<cutoff,])
ps.ord.envfit.scores.factors$Factor = "Categorical"
ps.ord.envfit.scores = rbind(ps.ord.envfit.scores.factors,ps.ord.envfit.scores.vectors)

# Add rowname for label
ps.ord.envfit.scores$Label = row.names(ps.ord.envfit.scores)
colnames(ps.ord.envfit.scores)[1:2] = c("NMDSx","NMDSy")

sample_data(ps.norm)$Latitude = as.numeric(sample_data(ps.norm)$Latitude)

# Load the plot
p = plot_ordination(ps.norm, ps.ord, color = "Interval", axes=c(2,4))
p = p + geom_segment(data=ps.ord.envfit.scores,aes(x=0,y=0,xend=NMDSx,yend=NMDSy),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="darkgrey")
p = p + geom_point(size=5,aes(shape=sample_data(ps.norm)$Layer_LandClass))
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
p = p + scale_shape_manual(name = "Soil Layer and \n Land Class", values=c(19, 21, 15, 22))
p = p + geom_text(data=ps.ord.envfit.scores,aes(label=Label,x=NMDSx,y=NMDSy),hjust=0.5, vjust=.5, size = 2, color = "black")
p

######################## Arrow additions for long to short #############
# For NMDS
x = data.frame(ps.ord$points)$MDS2
y = data.frame(ps.ord$points)$MDS4

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
head(Segments)

Segments = melt(Segments,id=c("Pairs","Pair_ID","Interval","Soil_layer","Layer_LandClass","Land_Class"))
head(Segments)

# Using reshape2 package and dcast function
Segments = reshape2::dcast(Segments,Pairs+Pair_ID+Soil_layer+Layer_LandClass+Land_Class~Interval+variable,value.var = "value",fun=mean)
head(Segments)
colnames(sample_data(ps.norm))

# Creating the same plot as above with the two years of data, linked by lines.
p = plot_ordination(ps.norm, ps.ord, color = "Interval", shape = "Layer_LandClass", axes=c(2,4))
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
p = p + scale_shape_manual(name = "Soil Layer and \n Land Class", values=c(19, 21, 15, 22))
p = p + geom_segment(data=Segments,aes(x=Long_x,y=Long_y,xend=Short_x,yend=Short_y),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="darkgrey")
p

### Plot arrows only, using ggforce package

# Get just the arrows we want to focus on
head(ps.ord.envfit.scores)
ps.ord.envfit.scores
ps.ord.envfit.scores.plot = ps.ord.envfit.scores[c("PCt_Conifer_Seedlings","pH","IntervalLong"),]
ps.ord.envfit.scores.plot

p = ggplot(Segments)
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
p = p + geom_segment(data=Segments,aes(x=Long_x,y=Long_y,xend=Short_x,yend=Short_y),
                     arrow = arrow(type="open",length = unit(0.3,"cm")), 
                     color="darkgrey",size=1)
p = p + ggforce::geom_link(data=Segments,aes(x=Long_x,y=Long_y,xend=Short_x,yend=Short_y, colour = stat(index)), size = 1, show.legend = FALSE)
p = p + scale_colour_gradient(low = "#009E73", high="#D55E00")
p = p + xlab("Axis 2") + ylab("Axis 4")
# Add envfit arrows
p = p + geom_segment(data=ps.ord.envfit.scores.plot,aes(x=0,y=0,xend=NMDSx,yend=NMDSy),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="darkgrey")
# Code to see the envfit arrow identity
#p = p + geom_text(data=ps.ord.envfit.scores,aes(label=Label,x=NMDSx,y=NMDSy),hjust=0.5, vjust=.5, size = 2, color = "black")
p


##### Adding isolines with ordisurf #####
# Working from https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/

ordi <- ordisurf(ps.ord ~ samdat.envfit$PCt_Conifer_Seedlings, choices = c(2,4)) #created the ordisurf object
ordi.grid <- ordi$grid #extracts the ordisurf object
head(ordi.grid)
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.ps <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.ps$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.ps.na <- data.frame(na.omit(ordi.ps)) #gets rid of the nas
ordi.ps.na #looks ready for plotting!

points = data.frame(ps.ord$points)
points$Sample_ID = row.names(points)
samdat.envfit$Sample_ID = row.names(samdat.envfit)
ps.plot = join(points,samdat.envfit,by="Sample_ID")

p = ggplot()
p = p + theme_bw() +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),strip.text.x = element_text(size=18, face="bold"),
        axis.title = element_text(size=22, face="bold"),axis.text = element_text(size=18, face = "bold"),
        legend.text = element_text(size=18),legend.title = element_text(size = 22),
        strip.background = element_rect(colour="white", fill="white"))
p = p + theme(plot.title=element_text(size = 25))
p = p + theme(plot.subtitle=element_text(size = 18))
p = p + scale_fill_manual(name = "Interval", values = c("#000000", "#FFFFFF"))
p = p + stat_contour(data = ordi.ps.na, aes(x = x, y = y, z = z,colour = ..level..), binwidth =0.1,size=1)
p = p + geom_point(data=ps.plot, aes(x=MDS2,y=MDS4, fill=Interval),size=2,alpha=0.75,shape=21)
p = p + geom_segment(data=Segments,aes(x=Long_x,y=Long_y,xend=Short_x,yend=Short_y),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="black")
p = p + labs(colour = "Conifer\nSeedlings (%)") + xlab("NMDS2")+ylab("NMDS4")
p = p + scale_colour_gradient( low="#D55E00",high = "#009E73")
p


