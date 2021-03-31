# Using Corncob to detect differential abundance or differential variance
# Corncob: https://github.com/bryandmartin/corncob
# Paper: https://imstat.org/journals-and-publications/annals-of-applied-statistics/annals-of-applied-statistics-next-issues/

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install(version = "3.11")

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
#install.packages("ggplot2")
#install.packages("dplyr")


# To install:
# install.packages("devtools")
#devtools::install_github("bryandmartin/corncob")
library(corncob)
library(phyloseq)
library(dplyr)
library(ggplot2)

#?differentialTest

ps = readRDS("ps.WB2016")
ps

#Remove mitochondria and chloroplast from ps; ensure the only domain is Bacteria
ps <- ps %>%
  subset_taxa(
    Domain == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast"
  )
ps

# Remove multi-community, sample 8, sample 21, and (samples that aren't paired) from ps object
Samples_To_Keep = c("WB16-01-L-M-2", "WB16-01-L-O-2", "WB16-02-L-M-2", "WB16-03-L-O-2", "WB16-03-L-M-1", "WB16-04-L-M-2", "WB16-04-L-O-1", "WB16-05-L-M-2", "WB16-05-L-O-2",
                    "WB16-06-L-M", "WB16-06-L-O", "WB16-07-L-M", "WB16-07-L-O", "WB16-09-L-M", "WB16-09-L-O", "WB16-10-L-M", "WB16-11-L-M", "WB16-11-L-O",  
                    "WB16-12-L-M", "WB16-12-L-O" ,"WB16-13-L-M", "WB16-13-L-O", "WB16-14-L-M", "WB16-15-L-M", "WB16-16-L-M", "WB16-16-L-O", "WB16-17-L-M",  
                    "WB16-17-L-O", "WB16-18-L-O", "WB16-19-L-M", "WB16-19-L-O", "WB16-20-L-O", "WB16-22-L-M", "WB16-22-L-O", "WB16-23-L-O", "WB16-24-L-O",
                    "WB16-01-S-M-2", "WB16-01-S-O-1", "WB16-02-S-M-2", "WB16-03-S-O-2", "WB16-03-S-M-2", "WB16-04-S-M-2", "WB16-04-S-O-1", "WB16-05-S-M-2", "WB16-05-S-O-1",
                    "WB16-06-S-M-1", "WB16-06-S-O-1", "WB16-07-S-M", "WB16-07-S-O", "WB16-09-S-M", "WB16-09-S-O", "WB16-10-S-M", "WB16-11-S-M", "WB16-11-S-O",
                    "WB16-12-S-M", "WB16-12-S-O", "WB16-13-S-M", "WB16-13-S-O", "WB16-14-S-M", "WB16-15-S-M", "WB16-16-S-M", "WB16-16-S-O", "WB16-17-S-M",
                    "WB16-17-S-O", "WB16-18-S-O", "WB16-19-S-M", "WB16-19-S-O",  "WB16-20-S-O", "WB16-22-S-M", "WB16-22-S-O",  "WB16-23-S-O", "WB16-24-S-O")

#Take out blanks, multicommunity, and samples from site 8 (allburn, burn, unburn)
ps.pruned =  prune_samples(Samples_To_Keep, ps)

ps.wet = subset_samples(ps.pruned,  Land_Class == "Wetland")
ps.wet

ps.up = subset_samples(ps.pruned, Land_Class == "Upland")
ps.up

# Just make sure there are no OTUs that are left over from the other prairie
ps.up = prune_taxa(taxa_sums(ps.up)>0,ps.up)
ps.up

ps.wet = prune_taxa(taxa_sums(ps.wet)>0,ps.wet)
ps.wet

ps.relabund.wet = transform_sample_counts(ps.wet, function(x) x / sum(x))
AbundTaxa = taxa_names(filter_taxa(ps.relabund.wet, function(x) mean(x) > 0.00001, TRUE))
ps.wet = prune_taxa(AbundTaxa,ps.wet)
ps.wet

ps.relabund.up = transform_sample_counts(ps.up, function(x) x / sum(x))
AbundTaxa = taxa_names(filter_taxa(ps.relabund.up, function(x) mean(x) > 0.00001, TRUE))
ps.up = prune_taxa(AbundTaxa,ps.up)
ps.up

#Run the differential test on the upland samples controlling for Pair ID and Soil layer
dT.up.horizon = differentialTest(formula = ~ Pair_ID + Soil_layer + Interval,
                                 phi.formula = ~  Pair_ID + Soil_layer + Interval,
                                 formula_null = ~ Pair_ID + Soil_layer,
                                 phi.formula_null = ~  Pair_ID + Soil_layer + Interval,
                                 test = "Wald", boot = FALSE,
                                 data = ps.up,
                                 fdr_cutoff = 0.05)
dT.up.horizon

#Run the differential test on the wetland samples controlling for Pair ID and Soil layer
dT.wet.horizon = differentialTest(formula = ~ Pair_ID + Soil_layer + Interval,
                                  phi.formula = ~  Pair_ID + Soil_layer + Interval,
                                  formula_null = ~ Pair_ID + Soil_layer,
                                  phi.formula_null = ~  Pair_ID + Soil_layer + Interval,
                                  test = "Wald", boot = FALSE,
                                  data = ps.wet,
                                  fdr_cutoff = 0.05)
dT.wet.horizon

#Save the output of the differential tests so that you don't have to run the models again
saveRDS(dT.wet.horizon, "dT.wet.horizon")
saveRDS(dT.up.horizon, "dT.up.horizon")


#Read in the model 
dT.up.horizon = readRDS("dT.up.horizon")
#Check out the output for the significant models
TEST = dT.up.horizon$significant_models[2]
TEST
TEST[[1]]$coefficients
TEST[[1]]$coefficients[16,]

# Create an empty dataframe to store the results
df.up.horizon = data.frame()
# Get the dataset for the treatment
r = dT.up.horizon
# Run through all the significant models, pulling out their model coefficients
for (j in 1:length(r$significant_taxa)){
  # Get the significant models
  sig_models = r$significant_models[[j]]
  # Pull out the long-interval coefficients (intercept) for all plots
  mu.1.L.up = data.frame(t(as.matrix(sig_models$coefficients[1,])))
  mu.10.L.up = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  mu.11.L.up = data.frame(t(as.matrix(sig_models$coefficients[3,])))
  mu.14.L.up = data.frame(t(as.matrix(sig_models$coefficients[4,])))
  mu.15.L.up = data.frame(t(as.matrix(sig_models$coefficients[5,])))
  mu.16.L.up = data.frame(t(as.matrix(sig_models$coefficients[6,])))
  mu.17.L.up = data.frame(t(as.matrix(sig_models$coefficients[7,])))
  mu.19.L.up = data.frame(t(as.matrix(sig_models$coefficients[8,])))
  mu.2.L.up = data.frame(t(as.matrix(sig_models$coefficients[9,])))
  mu.3.L.up = data.frame(t(as.matrix(sig_models$coefficients[10,])))
  mu.4.L.up = data.frame(t(as.matrix(sig_models$coefficients[11,])))
  mu.5.L.up = data.frame(t(as.matrix(sig_models$coefficients[12,])))
  mu.6.L.up = data.frame(t(as.matrix(sig_models$coefficients[13,])))
  mu.7.L.up = data.frame(t(as.matrix(sig_models$coefficients[14,])))
  mu.9.L.up = data.frame(t(as.matrix(sig_models$coefficients[15,])))
  # Pull out the organic horizon coefficients (estimates)
  mu.O.up = data.frame(t(as.matrix(sig_models$coefficients[16,])))
  # Pull out the Short interval coefficients (estimates)
  mu.S.up = data.frame(t(as.matrix(sig_models$coefficients[17,])))
  # Pop all the coefficients together
  mu.up.hor = cbind(mu.1.L.up, mu.10.L.up, mu.11.L.up, mu.14.L.up, mu.15.L.up, mu.16.L.up, mu.17.L.up, mu.19.L.up,mu.2.L.up,mu.3.L.up,
                    mu.4.L.up, mu.5.L.up, mu.6.L.up, mu.7.L.up, mu.9.L.up, mu.O.up, mu.S.up)
  # Bring along the p_fdr for that model
  p_fdr = r$p_fdr[r$significant_taxa][j]
  mu.up.hor$p_fdr = p_fdr
  # Get the row names (OTU IDs)
  row.names(mu.up.hor)= paste(row.names(data.frame(p_fdr)))
  # Add an OTU column with that OTU ID
  mu.up.hor$OTU = row.names(data.frame(p_fdr))
  # Add this result to the dataframe
  df.up.horizon = rbind(df.up.horizon,mu.up.hor)
}

head(df.up.horizon)

#Make a new dataframe with the estimates for each of the burned plots and the unburned reference, along with the other values of interest

df.up.horizon.temp = df.up.horizon[,c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,62,63,64,65,66,67,68,69,70)]
head(df.up.horizon.temp)

# Fix column names
colnames(df.up.horizon.temp) = c("mu.1.L.up", "mu.10.L.up", "mu.11.L.up", "mu.14.L.up", "mu.15.L.up", "mu.16.L.up", "mu.17.L.up", "mu.19.L.up",
                                 "mu.2.L.up","mu.3.L.up", "mu.4.L.up", "mu.5.L.up", "mu.6.L.up", "mu.7.L.up", "mu.9.L.up", "mu.O.up","se.O.up",
                                 "t.value.O.up","p.O.up", "mu.S.up","se.S.up","t.value.S.up","p.S.up","p_fdr.S.up","OTU")
head(df.up.horizon.temp)

# Calculate the fold-change value (relative abundances) using corncob's inverse logit function and our estimates
#First make a new column with the relative abundance calculation for Burned
#use corncob invlogit on each of the plots; site 1 is the baseline, the other sample sites are added to the baseline
# Add all the plots together and take the average

#Average relative abundance in mineral horizons in long interval
df.up.horizon.temp$RelAbund.M.L = (corncob::invlogit(df.up.horizon.temp$mu.1.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.10.L.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.11.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.14.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.15.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.16.L.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.17.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.19.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.2.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.3.L.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.4.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.5.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.6.L.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.7.L.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.9.L.up))/15

df.up.horizon.temp$RelAbund.L

#Average relative abundance in organic horizons in short interval
df.up.horizon.temp$RelAbund.O.L = (corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.10.L.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.11.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.14.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.15.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.16.L.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.17.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.19.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.2.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.3.L.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.4.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.5.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.6.L.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.7.L.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.9.L.up+ df.up.horizon.temp$mu.O.up))/15


df.up.horizon.temp$RelAbund.O.L

#Do the same for the short intervals as done above, but add the SI baseline to each to get the SI relative abundance

df.up.horizon.temp$RelAbund.M.S = (corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.10.L.up + df.up.horizon.temp$mu.S.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.11.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.14.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.15.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.16.L.up + df.up.horizon.temp$mu.S.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.17.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.19.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.2.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.3.L.up + df.up.horizon.temp$mu.S.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.4.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.5.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.6.L.up + df.up.horizon.temp$mu.S.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.7.L.up + df.up.horizon.temp$mu.S.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.9.L.up + df.up.horizon.temp$mu.S.up))/15

df.up.horizon.temp$RelAbund.M.S

df.up.horizon.temp$RelAbund.O.S = (corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.S.up + df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.10.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.11.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.14.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.15.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.16.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.17.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.19.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.2.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.3.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.4.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.5.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.6.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) +
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.7.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up) + 
                                     corncob::invlogit(df.up.horizon.temp$mu.1.L.up + df.up.horizon.temp$mu.9.L.up + df.up.horizon.temp$mu.S.up+ df.up.horizon.temp$mu.O.up))/15

df.up.horizon.temp$RelAbund.O.S


#Divide the LI relative abundance by the SI relative abundance to calculate the change

df.up.horizon.temp$FC.O = df.up.horizon.temp$RelAbund.O.L/df.up.horizon.temp$RelAbund.O.S

df.up.horizon.temp$FC.M = df.up.horizon.temp$RelAbund.M.L/df.up.horizon.temp$RelAbund.M.S

df.up.horizon.temp$FC.O
df.up.horizon.temp$FC.M

#Take the average of the change
df.up.horizon.temp$FC.Ave = (df.up.horizon.temp$FC.M+ df.up.horizon.temp$FC.O)/2
df.up.horizon.temp$FC.Ave

# calculate the log2-fold change
df.up.horizon.temp$log2FC = log(df.up.horizon.temp$FC.Ave,base=2)

dim(df.up.horizon.temp)
head(df.up.horizon.temp)
# Final summary of results
# Make a new column with the average of the relative abundances
df.up.horizon.temp$RelAbundAve = (df.up.horizon.temp$RelAbund.M.L + df.up.horizon.temp$RelAbund.M.S + df.up.horizon.temp$RelAbund.O.L + df.up.horizon.temp$RelAbund.O.S)/4 

# Now lets add taxonomy information onto OTUs
SigOTUs = levels(as.factor(df.up.horizon.temp$OTU))
SigOTUs
pruned = prune_taxa(SigOTUs,ps.up)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
joined.df.up.horizon = merge(df.up.horizon.temp,taxtab,by=c("OTU"))
head(joined.df.up.horizon)

# This is some code for fixing weird long taxon names
# optional

ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium", "Ambiguous_taxa",
               "uncultured actinobacterium","uncultured planctomycete", "KD4-96", "IMCC26256",
               "Pir4_lineage", "Candidatus_Udaeobacter", "CL500-29_marine_group", "Clostridium_sensu_stricto_13", "CL500-29_marine_group",
               "SC-I-84", "Subgroup_17")

joined.df.up.horizon = joined.df.up.horizon %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Order %in% ignoreList | is.na(Order),
                                     ifelse(Class %in% ignoreList | is.na(Class),
                                            ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                                   paste(Phylum)),paste(Class)),paste(Order)),paste(Family)),paste(Genus)))%>%
  mutate(Name = ifelse(Name == "Burkholderia-Caballeronia-Paraburkholderia","Burkholderia-Paraburkholderia",Name))

joined.df.up.horizon

#Save the dataframe
write.csv(joined.df.up.horizon,"WB2016-L2FC-upland-horizon.csv", row.names = FALSE)

#Edit the dataframe so that only taxa that have a log2-fold change greater than |1| are included
#Read in new filtered dataframe
joined.df.up.horizon.filter = read.csv("WB2016-L2FC-upland-horizon-filtered.csv")

#Create variables for OTUs that responded positively, negatively, or neutrally in Whitman et al 2019
up.pos= c("27091ac6f4a1bc4de379aaff8b26c2e0", "28f945c073bb3d4053e9a94681f23f0f", "3a3a0d5306d35311a5f8df1484216895", "421ec41e37bd6025a62989ad0c82ee06",
          "4c31de261384e0253e628daad94ff8a1", "86f7e172d4a477ac548274d7cf30cf16", "c9f92118179f3d2a39da07e1ac99b054")
up.neg=c("06635b8dad601ec0288f8486f1eafd6b", "154a339ce8de9cb6ace77dbd1bbc6637", "4129a88675d97baeb6afd72e704e018b", "698d116b7adc4e587b5b367ed431019b",
         "7437ce99cc8d67f6624b53caceeaa10a", "7e53c50a7f23d2bfc2bbcb669f2456e6", "85d9a60195d46b8da185f8783826ffdb", "b2a901ef1fdeea2def2812fb9ad5e6ff",
         "e40f1dfae487b24206f0e4b0f07b0416")
up.neu = c("305924d9157cfa3dd3995fb10f25eca3","4c0cd2741e0e8766614831387ff3002b", "68c06b9a2a4f15666ba7a3eb5d0bf845", "72ece8855c482bc1386a25b176bd5d79",
           "75b058f5a171051f33aa9e79c89a1bfb", "75c5a7d1de33ad872da2331e8a605c45", "97111cec0076e6bd198dc03af22050f9", "aed7a81abdbcb1704d3ee64359942d92",
           "c454c3dda1f78304fe7da25b0e2a37a4","c6d87077f5ec9e5285efbbd3b4622427", "c818726be7b56994682bf42acc20ba13", "ce7f3027544b10996132606f20aa5164")

#Add a column for the response to fire
joined.df.up.horizon.filter = joined.df.up.horizon.filter %>%
  mutate(Response = ifelse(OTU %in% up.pos, "Positive",
                           ifelse(OTU %in% up.neg, "Negative",
                                  ifelse(OTU %in% up.neu, "Neutral", NA))))
head(joined.df.up.horizon.filter)
joined.df.up.horizon.filter$Response = factor(joined.df.up.horizon.filter$Response, levels = c("Positive", "Neutral", "Negative"))

# Figure 5a
p = ggplot(joined.df.up.horizon.filter,aes(y=log2FC,color=Response, shape =Response, fill = Response,  x= reorder(Name, log2FC)))
p = p + theme_bw()
p = p + geom_point(aes(size = RelAbundAve))
p = p + scale_shape_manual(values = c(24,21,25))
p = p + scale_fill_manual(values = c("#009E73", "#F0E442","#CC79A7"))
p = p + theme(axis.text.x = element_text(angle=90,face="italic",color = "black", vjust=0,hjust=1, size = 20),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 18))
p = p + ylab("Log2-Fold Change in Relative Abundance in Long-\n vs. Short-Interval Reburns for Upland Sites") + xlab("")
p = p + scale_x_discrete(breaks=joined.df.up.horizon.filter$Name,labels=joined.df.up.horizon.filter$Name)
p = p + scale_size(name = "Average Relative Abundance", range = c(4,12))
p = p + scale_colour_manual(name = "OTU Response to Fire \nin Whitman et al 2019a", values = c("#009E73", "#F0E442","#CC79A7"))
p = p + geom_hline(yintercept=0)
p = p + guides(color = guide_legend(reverse = TRUE))
p


## SI Figure 5
#Read in non-filtered datafile
joined.df.up.horizon = read.csv("WB2016-L2FC-upland-horizon.csv")
#Add a column for the response to fire using variable from above
joined.df.up.horizon = joined.df.up.horizon %>%
  mutate(Response = ifelse(OTU %in% up.pos, "Positive",
                           ifelse(OTU %in% up.neg, "Negative",
                                  ifelse(OTU %in% up.neu, "Neutral", NA))))
head(joined.df.up.horizon)
joined.df.up.horizon$Response = factor(joined.df.up.horizon$Response, levels = c("Positive", "Neutral", "Negative"))

p = ggplot(joined.df.up.horizon,aes(y=log2FC,color=Response, shape =Response, fill = Response,  x= reorder(Name, log2FC)))
p = p + theme_bw()
p = p + geom_point(aes(size = RelAbundAve))
p = p + scale_shape_manual(values = c(24,21,25))
p = p + scale_fill_manual(values = c("#009E73", "#F0E442","#CC79A7"))
p = p + theme(axis.text.x = element_text(angle=90,face="italic",color = "black", vjust=0,hjust=1, size = 20),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 18))
p = p + ylab("Log2-Fold Change in Relative Abundance in Long-\n vs. Short-Interval Reburns for Upland Sites") + xlab("")
p = p + scale_x_discrete(breaks=joined.df.up.horizon.filter$Name,labels=joined.df.up.horizon.filter$Name)
p = p + scale_size(name = "Average Relative Abundance", range = c(4,12))
p = p + scale_colour_manual(name = "OTU Response to Fire \nin Whitman et al 2019a", values = c("#009E73", "#F0E442","#CC79A7"))
p = p + geom_hline(yintercept=0)
p = p + guides(color = guide_legend(reverse = TRUE))
p





## Do the same as above for wetlands
################## Calculate L2FC for Wetland controlling for Soil Horizon
TEST = dT.wet.horizon$significant_models[2]
TEST
TEST[[1]]$coefficients
TEST[[1]]$coefficients[16,]

# Create an empty dataframe to store the results
df.wet.horizon = data.frame()
# Get the dataset for the treatment
r = dT.wet.horizon
# Run through all the significant models, pulling out their model coefficients
for (j in 1:length(r$significant_taxa)){
  # Get the significant models
  sig_models = r$significant_models[[j]]
  # Pull out the long-interval coefficients (intercept) for all plots
  mu.12.L.wet = data.frame(t(as.matrix(sig_models$coefficients[1,])))
  mu.13.L.wet = data.frame(t(as.matrix(sig_models$coefficients[2,])))
  mu.18.L.wet = data.frame(t(as.matrix(sig_models$coefficients[3,])))
  mu.20.L.wet = data.frame(t(as.matrix(sig_models$coefficients[4,])))
  mu.22.L.wet = data.frame(t(as.matrix(sig_models$coefficients[5,])))
  mu.23.L.wet = data.frame(t(as.matrix(sig_models$coefficients[6,])))
  mu.24.L.wet = data.frame(t(as.matrix(sig_models$coefficients[7,])))
  # Pull out the Organic horizon coefficients (estimates)
  mu.O.wet = data.frame(t(as.matrix(sig_models$coefficients[8,])))
  # Pull out the Short interval coefficients (estimates)
  mu.S.wet = data.frame(t(as.matrix(sig_models$coefficients[9,])))
  # Pop all the coefficients together
  mu.wet.hor = cbind(mu.12.L.wet, mu.13.L.wet, mu.18.L.wet, mu.20.L.wet, mu.22.L.wet, mu.23.L.wet, mu.24.L.wet, mu.O.wet, mu.S.wet)
  # Bring along the p_fdr for that model
  p_fdr = r$p_fdr[r$significant_taxa][j]
  mu.wet.hor$p_fdr = p_fdr
  # Get the row names (OTU IDs)
  row.names(mu.wet.hor)= paste(row.names(data.frame(p_fdr)))
  # Add an OTU column with that OTU ID
  mu.wet.hor$OTU = row.names(data.frame(p_fdr))
  # Add this result to the dataframe
  df.wet.horizon = rbind(df.wet.horizon,mu.wet.hor)
}

head(df.wet.horizon)

#Make a new dataframe with the estimates for each of the burned plots and the unburned reference, along with the other values of interest

df.wet.horizon.temp = df.wet.horizon[,c(1,5,9,13,17,21,25,29,30,31,32,33,34,35,36,37,38)]
head(df.wet.horizon.temp)

# Fix column names
colnames(df.wet.horizon.temp) = c('mu.12.L.wet', 'mu.13.L.wet', 'mu.18.L.wet', 'mu.20.L.wet', 'mu.22.L.wet', 'mu.23.L.wet', 'mu.24.L.wet',
                                  "mu.O.wet","se.O.wet","t.value.O.wet","p.O.wet","mu.S.wet","se.S.wet","t.value.S.wet","p.S.wet","p_fdr.S.wet","OTU")
head(df.wet.horizon.temp)

#Average relative abundance in mineral horizons in long interval
df.wet.horizon.temp$RelAbund.M.L = (corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.13.L.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.18.L.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.20.L.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.22.L.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.23.L.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.24.L.wet))/7

df.wet.horizon.temp$RelAbund.M.L
#Average relative abundance in organic horizons in long interval
df.wet.horizon.temp$RelAbund.O.L = (corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.13.L.wet+ df.wet.horizon.temp$mu.O.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.18.L.wet+ df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.20.L.wet+ df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.22.L.wet+ df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.23.L.wet+ df.wet.horizon.temp$mu.O.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.24.L.wet+ df.wet.horizon.temp$mu.O.wet))/7

df.wet.horizon.temp$RelAbund.O.L
#Average relative abundance in mineral horizons in short interval
df.wet.horizon.temp$RelAbund.M.S = (corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.S.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.13.L.wet + df.wet.horizon.temp$mu.S.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.18.L.wet + df.wet.horizon.temp$mu.S.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.20.L.wet + df.wet.horizon.temp$mu.S.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.22.L.wet + df.wet.horizon.temp$mu.S.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.23.L.wet + df.wet.horizon.temp$mu.S.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.24.L.wet + df.wet.horizon.temp$mu.S.wet))/7

df.wet.horizon.temp$RelAbund.M.S
#Average relative abundance in organic horizons in short interval
df.wet.horizon.temp$RelAbund.O.S = (corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.S.wet + df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.13.L.wet + df.wet.horizon.temp$mu.S.wet + df.wet.horizon.temp$mu.O.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.18.L.wet + df.wet.horizon.temp$mu.S.wet + df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.20.L.wet + df.wet.horizon.temp$mu.S.wet + df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.22.L.wet + df.wet.horizon.temp$mu.S.wet + df.wet.horizon.temp$mu.O.wet) +
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.23.L.wet + df.wet.horizon.temp$mu.S.wet + df.wet.horizon.temp$mu.O.wet) + 
                                      corncob::invlogit(df.wet.horizon.temp$mu.12.L.wet + df.wet.horizon.temp$mu.24.L.wet + df.wet.horizon.temp$mu.S.wet+ df.wet.horizon.temp$mu.O.wet))/7

df.wet.horizon.temp$RelAbund.O.S


#Divide the LI relative abundance by the SI relative abundance to calculate the change

df.wet.horizon.temp$FC.O = df.wet.horizon.temp$RelAbund.O.L/df.wet.horizon.temp$RelAbund.O.S

df.wet.horizon.temp$FC.M = df.wet.horizon.temp$RelAbund.M.L/df.wet.horizon.temp$RelAbund.M.S

df.wet.horizon.temp$FC.O
df.wet.horizon.temp$FC.M
#Get the average of the change
df.wet.horizon.temp$FC.Ave = (df.wet.horizon.temp$FC.M+ df.wet.horizon.temp$FC.O)/2
df.wet.horizon.temp$FC.Ave

# Calculate the log2-fold change
df.wet.horizon.temp$log2FC = log(df.wet.horizon.temp$FC.Ave,base=2)

dim(df.wet.horizon.temp)
head(df.wet.horizon.temp)
# Final summary of results
#Calculate the average relative abundances
df.wet.horizon.temp$RelAbundAve = (df.wet.horizon.temp$RelAbund.M.L + df.wet.horizon.temp$RelAbund.M.S + df.wet.horizon.temp$RelAbund.O.L + df.wet.horizon.temp$RelAbund.O.S)/4 


#Now let's add taxonomy data onto the dataframe
SigOTUs = levels(as.factor(df.wet.horizon.temp$OTU))
SigOTUs
pruned = prune_taxa(SigOTUs,ps.wet)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
joined.df.wet.horizon = merge(df.wet.horizon.temp,taxtab,by=c("OTU"))
head(joined.df.wet.horizon)

# This is some code for fixing weird long taxon names
ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium", "Ambiguous_taxa",
               "uncultured actinobacterium","uncultured planctomycete", "KD4-96", "Burkholderia-Caballeronia-Paraburkholderia", "IMCC26256",
               "Pir4_lineage", "Candidatus_Udaeobacter","SC-I-84", "Subgroup_17", "Rhizobiales_Incertae_Sedis", "Clostridium_sensu_stricto_13",
               "CL500-29_marine_group", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "TK10", "RB41")

joined.df.wet.horizon.filter = joined.df.wet.horizon.filter %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Order %in% ignoreList | is.na(Order),
                                     ifelse(Class %in% ignoreList | is.na(Class),
                                            ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                                   paste(Phylum)),paste(Class)),paste(Order)),paste(Family)),paste(Genus)))


#Save the dataframe
write.csv(joined.df.wet.horizon.filter,"WB2016-L2FC-wetland-horizon.csv", row.names = FALSE)

#Edit the dataframe so that only taxa that have a log2-fold change greater than |1| are included
#Read in new filtered dataframe
joined.df.wet.horizon.filter = read.csv("WB2016-L2FC-wetland-horizon-filtered.csv")

#Create variables for OTUs that responded positively, negatively, or neutrally in Whitman et al 2019
wet.pos= c("00c20c350c8c9aec561cc31befd8bb3f", "4854c719a46573ffb50f247c98af53a2", "4c31de261384e0253e628daad94ff8a1","674e9ef55b4e7d2965be9ed17b4e1a5a")
wet.neg=c("009b595c11c059ac75d090dafbb73e1c", "06635b8dad601ec0288f8486f1eafd6b", "06fe572905355a9b3fb61430a9279c1d", "18273717451643b4d481641dccee6348",
          "195bfb31f048e2f2dfbb0c7794fa147b", "1aa1e533395182d33feeebf0f1ca57ff", "2bc736f98d35194941c19df9a5490295", "4129a88675d97baeb6afd72e704e018b",
          "438986c28abec79e919563635f5d36c2", "4fe476b09a972d5737628a5182707144", "5748c96536a64033a6b47aa62a5baa9a", "7e53c50a7f23d2bfc2bbcb669f2456e6",
          "b5ac620241ecf15777fee8b68cea8861","d87c64b716251a82a270c8c6525c4a3d", "e0adeee47af9e320123cd353779f5bf7", "e40f1dfae487b24206f0e4b0f07b0416",
          "e43c145a42390c008bac59a5121e5992", "f8250040feff959fae4696e5f921cdaa")
wet.neu = c("04144d664c09655b570dac627bb18d0d", "08d77f4abd024cb55405c7e82020cf4c", "0b6df43ceaac857fa88f40e4867375e0", "0c00a12784a73646e3b6f7f494ab3b60",
            "305924d9157cfa3dd3995fb10f25eca3", "35f1139b7164644f9d4c97a5128b3c1a", "3aa82c0bfe1377d23d0ed48afefd8b19", "7558d6c7ef1f5cf9a1744b6e6ca6c743",
            "75c5a7d1de33ad872da2331e8a605c45", "767a51dc9c7bc79c71773c9f0c5f8c2c", "7923647be4c6141e042a28069bb18469", "7bb5668fbad7d45555c078be05f3b182",
            "804aea6b174f827ac3aeab08f2981e52", "97111cec0076e6bd198dc03af22050f9", "986580c4374e515d81f09226760dbbb7", "aed7a81abdbcb1704d3ee64359942d92",
            "c9ef89952e819877d3acc77face3dc14", "d21be3824021c7d43400132c1b9372e2", "eb65877471de7307a150ff363dcbc082")

#Add a column for the OTU response
joined.df.wet.horizon.filter = joined.df.wet.horizon.filter %>%
  mutate(Response = ifelse(OTU %in% wet.pos, "Positive",
                           ifelse(OTU %in% wet.neg, "Negative",
                                  ifelse(OTU %in% wet.neu, "Neutral", NA))))
unique(joined.df.wet.horizon.filter$OTU)

## Figure 5b
p = ggplot(joined.df.wet.horizon.filter,aes(y=log2FC,color=Response, shape =Response, fill = Response,  x= reorder(Name, log2FC)))
p = p + theme_bw()
p = p + geom_point(aes(size = RelAbundAve))
p = p + scale_shape_manual(values = c(25,21,24))
p = p + scale_fill_manual(values = c("#CC79A7", "#F0E442","#009E73"))
p = p + theme(axis.text.x = element_text(angle=90,face="italic",color = "black", vjust=0,hjust=1, size = 20),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 18))
p = p + ylab("Log2-Fold Change in Relative Abundance in Long-\n vs. Short-Interval Reburns for Wetland Sites") + xlab("")
p = p + scale_x_discrete(breaks=joined.df.wet.horizon.filter$Name,labels=joined.df.wet.horizon.filter$Name)
p = p + scale_size(name = "Average Relative Abundance", range = c(4,12))
p = p + scale_colour_manual(name = "OTU Response to Fire \nin Whitman et al 2019a", values = c("#CC79A7", "#F0E442","#009E73"))
p = p + geom_hline(yintercept=0)
p = p + guides(color = guide_legend(reverse = TRUE))
p

## SI Figure 6
#Read in non-filtered dataframe
joined.df.wet.horizon = read.csv("WB2016-L2FC-wetland-horizon.csv")
#Add a column for the OTU response using variable created above
joined.df.wet.horizon = joined.df.wet.horizon %>%
  mutate(Response = ifelse(OTU %in% wet.pos, "Positive",
                           ifelse(OTU %in% wet.neg, "Negative",
                                  ifelse(OTU %in% wet.neu, "Neutral", NA))))


p = ggplot(joined.df.wet.horizon,aes(y=log2FC,color=Response, shape =Response, fill = Response,  x= reorder(Name, log2FC)))
p = p + theme_bw()
p = p + geom_point(aes(size = RelAbundAve))
p = p + scale_shape_manual(values = c(25,21,24))
p = p + scale_fill_manual(values = c("#CC79A7", "#F0E442","#009E73"))
p = p + theme(axis.text.x = element_text(angle=90,face="italic",color = "black", vjust=0,hjust=1, size = 20),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 17),
              legend.title = element_text(size = 18))
p = p + ylab("Log2-Fold Change in Relative Abundance in Long-\n vs. Short-Interval Reburns for Wetland Sites") + xlab("")
p = p + scale_x_discrete(breaks=joined.df.wet.horizon.filter$Name,labels=joined.df.wet.horizon.filter$Name)
p = p + scale_size(name = "Average Relative Abundance", range = c(4,12))
p = p + scale_colour_manual(name = "OTU Response to Fire \nin Whitman et al 2019a", values = c("#CC79A7", "#F0E442","#009E73"))
p = p + geom_hline(yintercept=0)
p = p + guides(color = guide_legend(reverse = TRUE))
p








##### Post-hoc tests
#Pull out cale. taxa from wb2016-upland-melt to see if it has anything to do with conifer recruitment
library(vegan)

wb2016.melt= readRDS('wb2016-melt')
wb2016.melt
wb2016.melt.up = wb2016.melt %>%
  filter(Land_Class == "Upland")
wb2016.melt.wet = wb2016.melt %>%
  filter(Land_Class == "Wetland")

head(wb2016.melt.up)
wb2016.melt.up.cal = wb2016.melt.up %>%
  filter(OTU == "3a3a0d5306d35311a5f8df1484216895")
wb2016.melt.up.cal


anova.conseed.dens.cal = aov(Understory_Conifer_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                            Recent_Interval + Abundance, data=wb2016.melt.up.cal)
summary(anova.conseed.dens.cal)

anova.BLseed.dens.cal = aov(Understory_Broadleaf_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                               Recent_Interval + Abundance, data=wb2016.melt.up.cal)
summary(anova.BLseed.dens.cal)


sample_data(ps)

#Do with Blastococcus
#upland

wb2016.melt.blasto.up = wb2016.melt.up %>%
  filter(OTU == "4c31de261384e0253e628daad94ff8a1")
head(wb2016.melt.blasto.up)


anova.conseed.dens.blasto.up = aov(Understory_Conifer_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                               Recent_Interval + Abundance, data=wb2016.melt.blasto.up)
summary(anova.conseed.dens.blasto.up)

anova.BLseed.dens.blasto.up  = aov(Understory_Broadleaf_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                              Recent_Interval + Abundance, data=wb2016.melt.blasto.up)
summary(anova.BLseed.dens.blasto.up)

#wetland

wb2016.melt.blasto.wet = wb2016.melt.wet %>%
  filter(OTU == "4c31de261384e0253e628daad94ff8a1")
head(wb2016.melt.blasto.wet)

anova.allseed.dens.blasto.wet = aov(Understory_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                                      Recent_Interval + Abundance, data=wb2016.melt.blasto.wet)
summary(anova.allseed.dens.blasto.wet)

anova.conseed.dens.blasto.wet = aov(Understory_Conifer_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                                     Recent_Interval + Abundance, data=wb2016.melt.blasto.wet)
summary(anova.conseed.dens.blasto.wet)

anova.BLseed.dens.blasto.wet  = aov(Understory_Broadleaf_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                                     Recent_Interval + Abundance, data=wb2016.melt.blasto.wet)
summary(anova.BLseed.dens.blasto.wet)

#Do for Rhizobia (wetland only)

wb2016.melt.Rhizo.wet = wb2016.melt.wet %>%
  filter(OTU == "986580c4374e515d81f09226760dbbb7")
head(wb2016.melt.Rhizo.wet)


anova.conseed.dens.Rhizo.wet = aov(Understory_Conifer_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                                     Recent_Interval + Abundance, data=wb2016.melt.Rhizo.wet)
summary(anova.conseed.dens.Rhizo.wet)

anova.BLseed.dens.Rhizo.wet  = aov(Understory_Broadleaf_Stems_Ha_Live ~ Overstory_Stems_Ha_Live + Moisture_Numeric + TSLF_At_Sample +
                                     Recent_Interval + Abundance, data=wb2016.melt.Rhizo.wet)
summary(anova.BLseed.dens.Rhizo.wet)


sample_data(ps)
