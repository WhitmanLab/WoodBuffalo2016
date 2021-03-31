library(reshape)
library(ggplot2)
library(phyloseq)
library(dplyr)

ps = import_biom("WB2016_OTU_table/feature-table-metaD-tax_json.biom" , parseFunction = parse_taxonomy_default)
head(tax_table(ps))

##Fixing the tax_table object to "clean-up" the column names

x = data.frame(tax_table(ps))
# Making a dummy variable to store the taxonomy data

colnames(x) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Assigning the proper column names instead of SILVA ranks

x$Domain = gsub("d__", "", as.character(x$Domain))
x$Phylum = gsub("p__", "", as.character(x$Phylum))
x$Class = gsub("c__", "", as.character(x$Class))
x$Order = gsub("o__", "", as.character(x$Order))
x$Family = gsub("f__", "", as.character(x$Family))
x$Genus = gsub("g__", "", as.character(x$Genus))
x$Species = gsub("s__", "", as.character(x$Species))
# Substituting the characters we don't want with nothing in the taxonomy

x=tax_table(as.matrix(x,dimnames=list(row.names(x),colnames(x))))
# Turning it into a taxonomy table, while saving the rownames and column names
tax_table(ps)=x
# Reassigning the taxonomy table in ps_xxx to the new modified one

head(tax_table(ps))
# Check for success

# Save phyloseq object
#saveRDS(ps, "ps.WB2016")