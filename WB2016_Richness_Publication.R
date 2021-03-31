library(phyloseq)
library(tidyverse)
library(breakaway)
library(ggplot2)
library(dplyr)
library(PairedData)
library(ggpubr)

# Read in phyloseq object
ps = readRDS("ps.WB2016")
ps

# Get rid of taxa that are characterized as mitochondria or chloroplasts
ps <- ps %>%
  subset_taxa(
    Domain == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast"
  )
ps

sample_data(ps)

# Remove multi-community, sample 8, sample 21, and (samples that aren't paired) from ps object
Samples_To_Keep = c("WB16-01-L-M-2", "WB16-01-L-O-2", "WB16-02-L-M-2", "WB16-03-L-O-2", "WB16-03-L-M-1", "WB16-04-L-M-2", "WB16-04-L-O-1", "WB16-05-L-M-2", "WB16-05-L-O-2",
                    "WB16-06-L-M", "WB16-06-L-O", "WB16-07-L-M", "WB16-07-L-O", "WB16-09-L-M", "WB16-09-L-O", "WB16-10-L-M", "WB16-11-L-M", "WB16-11-L-O",  
                    "WB16-12-L-M", "WB16-12-L-O" ,"WB16-13-L-M", "WB16-13-L-O", "WB16-14-L-M", "WB16-15-L-M", "WB16-16-L-M", "WB16-16-L-O", "WB16-17-L-M",  
                    "WB16-17-L-O", "WB16-18-L-O", "WB16-19-L-M", "WB16-19-L-O", "WB16-20-L-O", "WB16-22-L-M", "WB16-22-L-O", "WB16-23-L-O", "WB16-24-L-O",
                    "WB16-01-S-M-2", "WB16-01-S-O-1", "WB16-02-S-M-2", "WB16-03-S-O-2", "WB16-03-S-M-2", "WB16-04-S-M-2", "WB16-04-S-O-1", "WB16-05-S-M-2", "WB16-05-S-O-1",
                    "WB16-06-S-M-1", "WB16-06-S-O-1", "WB16-07-S-M", "WB16-07-S-O", "WB16-09-S-M", "WB16-09-S-O", "WB16-10-S-M", "WB16-11-S-M", "WB16-11-S-O",
                    "WB16-12-S-M", "WB16-12-S-O", "WB16-13-S-M", "WB16-13-S-O", "WB16-14-S-M", "WB16-15-S-M", "WB16-16-S-M", "WB16-16-S-O", "WB16-17-S-M",
                    "WB16-17-S-O", "WB16-18-S-O", "WB16-19-S-M", "WB16-19-S-O",  "WB16-20-S-O", "WB16-22-S-M", "WB16-22-S-O",  "WB16-23-S-O", "WB16-24-S-O")

#Used 10A values for data
#Take out blanks, multicommunity, and samples from site 8 (allburn, burn, unburn)
ps.pruned =  prune_samples(Samples_To_Keep, ps)
ps.pruned

#Make a new column in sample data that indicates the "class" of the leading tree species 

Conifer = c('PICEMAR', 'PINUBAN','PICEGLA','PICEMAR/LARILAR','LARILAR' )
Broadleaf = c('POPUTRE','POPUBAL', 'POPUBAL/POPUTRE')
Mixed = c('PICEMAR/POPUTRE', 'PINUBAN/POPUTRE','POPUTRE/PINUBAN')


sample_data(ps.pruned)$Leading_Type = ifelse(sample_data(ps.pruned)$Leading_Species %in% Conifer, "Conifer",
                                           ifelse(sample_data(ps.pruned)$Leading_Species %in% Broadleaf, "Broadleaf",
                                                  ifelse(sample_data(ps.pruned)$Leading_Species %in% Mixed, "Mixed", NA)))

# Drop the OTUs that have zeros across all remaining samples
ps.nozero = subset_taxa(ps.pruned,taxa_sums(ps.pruned)>0)
ps.nozero

# Collect sample data
SamDat = data.frame(sample_data(ps.nozero))

### Making a plot of OTU abundances across full dataset
d = data.frame(taxa_sums(ps.nozero))
d$OTU = row.names(d)
colnames(d) = c("Total","OTU")

# Arrange by abundance (Total) of each OTU and add order variable (X) 
d = d %>%
  arrange(-Total)%>%
  mutate(X=c(1:dim(d)[1]))

# Plot abundances in order of abundances (need to log Total because so large)
p = ggplot(d,aes(x=X,y=log(Total)))
p = p + geom_bar(stat="identity")
p

(sample_data(ps.nozero))

# Create function to calculate rarefaction curves

# Code from https://github.com/joey711/phyloseq/issues/143
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# Run the function on our data...takes just a few minutes
rarefaction_curve_data = calculate_rarefaction_curves(ps.nozero, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))

# Summarize it
rarefaction_curve_data_summary = ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

# Fix the Sample naming issue
#rarefaction_curve_data_summary$Sample = substring(rarefaction_curve_data_summary$Sample,2)

#Fix . in sample names to - to match ps
string <- rarefaction_curve_data_summary$Sample
rarefaction_curve_data_summary$Sample = gsub('\\.', '-', string)
head(rarefaction_curve_data_summary$Sample)


# Join the sample data to it
rarefaction_curve_data_summary_verbose = merge(rarefaction_curve_data_summary, data.frame(sample_data(ps.nozero)), by.x = 'Sample', by.y = 'row.names')

# Just pull out the observed diversity for now (not looking at Shannon)
rarefaction_curve_data_summary_verbose.observed = rarefaction_curve_data_summary_verbose %>% filter(Measure=="Observed")

# Plot the data, wrapped by times mixed
ggplot(
  data = rarefaction_curve_data_summary_verbose.observed,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Soil_layer,
    shape = Land_Class,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(~Interval,scales="free_x") + theme_bw()

# Plot all the data on top of each other
ggplot(
  data = rarefaction_curve_data_summary_verbose.observed,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Interval,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
)  + theme_bw()


# Run Breakaway on WB2016 dataset to estimate bacterial richness

SD = data.frame(sample_data(ps.nozero))
row.names(SD)
SD$SampleID = row.names(SD)
SD

unique(SD$SampleID)
#ps.TM = prune_samples(sample_data(ps.nozero)$AgitControl != "Y",ps.nozero)
#ps.TM = merge_samples(ps.TM, "TimesMixed")
sample_data(ps.nozero)

# Set OTU table
otu_data = t(otu_table(ps.nozero))
# Set sample data
meta_data = sample_data(ps.nozero)
# Had to flip OTU table so rownames match sample data
head(row.names(otu_data) == row.names(meta_data))

# Run Breakaway's frequency table list function
frequencytablelist = build_frequency_count_tables(otu_data)

# Try Breakaway on a couple of samples
breakaway(frequencytablelist[[1]])
breakaway(frequencytablelist[[2]])

# Because no plot pops up, we know that we're dealing with the WLRM
# That's because dada2 won't allow singletons to pass through

# Run the richness estimator (breakaway) on all our samples (lists of frequency tables)
RichEsts = lapply(frequencytablelist,breakaway)

# Pull out the estimates, errors, and the model
Estimate = as.matrix(map(RichEsts, "estimate"))
Error = as.matrix(map(RichEsts, "error"))
Model = as.matrix(map(RichEsts, "model"))
df = data.frame(Estimate,Error,Model)

# Add sample ID column, estimate, and error
df$SampleID = row.names(df)
df$Estimate=as.numeric(df$Estimate)
df$Error=as.numeric(df$Error)

head(df)

SD$TSLF_At_Sample = factor(SD$TSLF_At_Sample, levels = c(1,2,3,8,12,21), ordered = TRUE)


# Merge the estimates with the sample data
head(SD)
head(df)
RichPlot = merge(SD,df,by="SampleID")
head(RichPlot)

unique(RichPlot$Pair_ID)

RichPlot$Pair_ID = factor(RichPlot$Pair_ID, levels = c("1","2","3","4","5","6","7","9","10A","11","12","13",
                                                             "14","15","16","17","18","19","20","22","23","24"), ordered = TRUE)
RichPlot

RichPlot$Site_Interval_Depth = paste(RichPlot$Pair_ID, RichPlot$Interval, RichPlot$Soil_layer, sep = "_")
head(RichPlot)
RichPlot$Model = as.character(RichPlot$Model)
write.csv(RichPlot,"WB2016-richness-df.csv", row.names = FALSE)

#############################################################################################

# Plot them
cbPalette <- c("#009E73", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#CC79A7")


df = read.csv("WB2016_analysis_dataframe.csv")
head(df)


#Plots for paper:
#More negative means that SI had higher richness

#add column in analysis dataframe for land class
df$LandClass = ifelse(df$Wetland == "no", "Upland", "Wetland")
df

df$LandClass = ifelse(df$Wetland == "no", "Upland", "Wetland")
df$Layer_LandClass = paste(df$Soil_Layer, df$LandClass, sep = "-")
df$Layer_LandClass = levels = c('M-Upland', 'M-Wetland','O-Upland', 'O-Wetland', ordered = TRUE)
df

## SI Table 5: ANOVA analysis on percent change in richness

#Vegetation Transition
anova.Veg = aov(Percent_change_richness ~ Layer_LandClass + Veg_Transition, data=df)
summary(anova.Veg)

lm.veg = lm(Percent_change_richness ~ Layer_LandClass + Veg_Transition, data=df)
anova(lm.veg)

summary(lm.veg)

#Wetland
anova.wet = aov(Percent_change_richness ~ Soil_Layer + Wetland, data=df)
summary(anova.wet)

lm.wet = lm(Percent_change_richness ~ Soil_Layer + Wetland, data=df)
anova(lm.wet)

summary(lm.wet)

#FRI
anova.FRI = aov(Percent_change_richness ~ Layer_LandClass + Diff_in_Recent_Interval, data=df)
summary(anova.FRI)

lm.FRI = lm(Percent_change_richness ~ Layer_LandClass + Diff_in_Recent_Interval, data=df)
anova(lm.FRI)

summary(lm.FRI)

#TSLF
anova.TSLF = aov(Percent_change_richness ~ Layer_LandClass + TSLF_At_Sample, data=df)
summary(anova.TSLF)

lm.TSLF = lm(Percent_change_richness ~ Layer_LandClass + TSLF_At_Sample, data=df)
anova(lm.TSLF)

summary(lm.TSLF)

#Number of stems
anova.stems = aov(Percent_change_richness ~ Layer_LandClass + Diff_total_stems, data=df)
summary(anova.stems)

lm.stems = lm(Percent_change_richness ~ Layer_LandClass + Diff_total_stems, data=df)
anova(lm.stems)

summary(lm.stems)



########################## Perform paired ttest on breakaway richness values ###############################################
RichPlot = read.csv("WB2016-richness-df.csv")
RichPlot$SampleID

#Make a column with the site ID and Depth to use later on in the ttest box plot
RichPlot$Site_Depth = paste(RichPlot$Pair_ID, RichPlot$Soil_layer, sep = "_")
head(RichPlot)

#Separate data by drainage class and soil layer
RichPlot.Wet = RichPlot %>%
  filter(Land_Class == "Wetland")
head(RichPlot.Wet)

RichPlot.Up = RichPlot %>%
  filter(Land_Class == "Upland")
head(RichPlot.Up)

RichPlot.Wet.M = RichPlot.Wet %>%
  filter(Soil_layer == "M")
RichPlot.Wet.M

RichPlot.Wet.O = RichPlot.Wet %>%
  filter(Soil_layer == "O")
RichPlot.Wet.O

RichPlot.Up.M = RichPlot.Up %>%
  filter(Soil_layer == "M")
RichPlot.Up.M

RichPlot.Up.O = RichPlot.Up %>%
  filter(Soil_layer == "O")
RichPlot.Up.O

RichPlot.M = RichPlot %>%
  filter(Soil_layer == "M")
RichPlot.M

RichPlot.O = RichPlot %>%
  filter(Soil_layer == "O")
RichPlot.O

colnames(RichPlot.O)

# Figure 4: Richness box plots
#Wetland Mineral
compare_means(Bacterial_Richness_Estimate ~ Interval, data = RichPlot.Wet.M, paired = TRUE, method = "t.test", ref.group = "Site_Depth")
options(repr.plot.width=10, repr.plot.height=15)
p = ggpaired(RichPlot.Wet.M, x = "Interval", y = "Bacterial_Richness_Estimate", id = "Site_Depth",
             color = "Interval", point.size = 3, line.color = "black", line.size = 0.4,
             font.label = list(size = 15, color = "black", style = "bold")) +
  stat_compare_means(paired = TRUE, method = "t.test")
p = p + ylab("Estimated Richness of Bacterial Community\n in Wetland Mineral Horizon")
p = p + xlab("Interval")
p = p + theme(axis.title.x = element_text(size=22, face="bold"),
              axis.title.y = element_text(size=22, face="bold"),
              axis.text.x = element_text(size=22, face = "bold"),
              axis.text.y = element_text(size=22, face = "bold"),
              legend.position = "none")
p = p + ylim(0,1500)
p = p + scale_colour_manual(values=cbPalette)
p

#Wetland Organic
p = ggpaired(RichPlot.Wet.O, x = "Interval", y = "Bacterial_Richness_Estimate", id = "Site_Depth",
             color = "Interval", point.size = 3, line.color = "black", line.size = 0.4,
             font.label = list(size = 15, color = "black", style = "bold")) +
  stat_compare_means(paired = TRUE, method = "t.test")
p = p + ylab("Estimated Richness of Bacterial Community\n in Wetland Organic Horizon")
p = p + xlab("Interval")
p = p + theme(axis.title.x = element_text(size=22, face="bold"),
              axis.title.y = element_text(size=22, face="bold"),
              axis.text.x = element_text(size=22, face = "bold"),
              axis.text.y = element_text(size=22, face = "bold"),
              legend.position = "none")
p = p + ylim(0,1500)
p = p + scale_colour_manual(values=cbPalette)
p

#Upland Mineral
p = ggpaired(RichPlot.Up.M, x = "Interval", y = "Bacterial_Richness_Estimate", id = "Site_Depth",
             color = "Interval", point.size = 3, line.color = "black", line.size = 0.4,
             font.label = list(size = 15, color = "black", style = "bold")) +
  stat_compare_means(paired = TRUE, method = "t.test")
p = p + ylab("Estimated Richness of Bacterial Community\n in Upland Mineral Horizon")
p = p + xlab("Interval")
p = p + theme(axis.title.x = element_text(size=22, face="bold"),
              axis.title.y = element_text(size=22, face="bold"),
              axis.text.x = element_text(size=22, face = "bold"),
              axis.text.y = element_text(size=22, face = "bold"),
              legend.position = "none")
p = p + ylim(0,1500)
p = p + scale_colour_manual(values=cbPalette)
p

#Upland Organic
p = ggpaired(RichPlot.Up.O, x = "Interval", y = "Bacterial_Richness_Estimate", id = "Site_Depth",
             color = "Interval", point.size = 3, line.color = "black", line.size = 0.4,
             font.label = list(size = 15, color = "black", style = "bold")) +
  stat_compare_means(paired = TRUE, method = "t.test")
p = p + ylab("Estimated Richness of Bacterial Community\n in Upland Organic Horizon")
p = p + xlab("Interval")
p = p + theme(axis.title.x = element_text(size=22, face="bold"),
              axis.title.y = element_text(size=22, face="bold"),
              axis.text.x = element_text(size=22, face = "bold"),
              axis.text.y = element_text(size=22, face = "bold"),
              legend.position = "none")
              #legend.position = "left",
              #legend.title = element_text(size = 20),
              #legend.text = element_text(size = 17))
p = p + ylim(0,1500)
p = p + scale_colour_manual(values=cbPalette)
p


#SI Table 6: Run model to predict seedling count using multiple variables
anova.seed.dens.bac = aov(Total_Understory_Stems_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                        Recent_Interval + Bacterial_Richness_Estimate,  data=RichPlot)
aov.seed.dens.bac = aov(Total_Understory_Stems_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                            Recent_Interval + Bacterial_Richness_Estimate,  data=RichPlot)

summary(aov.seed.dens.bac)
