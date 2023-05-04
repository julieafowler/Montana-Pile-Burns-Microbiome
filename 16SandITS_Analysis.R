#### Scripts for Montana Pile Burn Project - 16S & ITS Analyses ####

## 16S used for all examples, exclude copy number analysis for ITS

## Written and compiled by Julie A. Fowler and Amelia R. Nelson


#### Package Loading ####
library(phyloseq) ## if needed to download, will have to get from Bioconductor 
library(ggplot2)
library(grid)
library(plyr)
library(vegan)
library(Hmisc)
library(reshape2)
library(ggpubr)
library(skimr)
library(ggthemr) ## if needed to download, will have to get from GitHub, only needed if using ggthemr('dust')

ggthemr('dust')
swatch()


#### Set Working Directory ####
setwd("/Users/juliefowler/OneDrive - Colostate/PhD Work/Studies/Burn Piles/Montana Burn Pile Work/Montana-Burn-Piles-Project")


#### Import Files and Create Object ####

otus<-read.delim("16S-feature-table.tsv",header=T,row.names=1, check.names=FALSE) ## load in the feature table - this is your 'ASV' table with counts across samples
map_file<-read.delim("16S-map-file.txt",header = T,row.names=1,check.names=FALSE)


## Load in your map file 
## Make sure the chem variables in map file have outliers removed
all.equal(names(otus),row.names(map_file)) ## check to make sure row names and column names equal (aka the sample names), if this reads 'TRUE' you are good to move on!

## Ordering variable for ggplot
map_file$Layer <- factor(map_file$Layer, levels = c("Ash", "Organic", "Mineral", "Unburned Mineral"))


## Adding a column combining LWL and LWR under "Lonesome Wood" for a comparison of LWL&LWR vs. Rendes., if needed 
# library(stringr)
# library(dplyr)
# map_file$Site2<-map_file$Site
# map_file <- map_file %>% mutate(Site2 = ifelse(str_detect(map_file$Site, fixed("Lonesome Wood")), "Lonesome Wood", map_file$Site2))
# map_file <- map_file %>% mutate(Site2 = ifelse(str_detect(map_file$Site, fixed("Rendesvous")), "Rendesvous", map_file$Site2))

# Parsing data frame if only analyzing subset of data - for ANOSIM
# map_file<-subset(map_file, Layer == "Mineral" | Layer == "A")
# fin_sample_names<-row.names(map_file) #final list of subsetted sample names
# 
# totus<-t(otus) #Transpose the otu table
# totus.subset<-subset(totus,rownames(totus) %in% fin_sample_names)
# fin.otus<-t(totus.subset)
# all.equal(colnames(fin.otus),row.names(map_file))


## if parsing, re-do the stuff below to bray-curtis, replace otus with fin.otus
## otumat<-as.matrix(otus) if not parsing
## otumat<-as.matrix(fin.otus) if parsing 

otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)
taxa<-read.delim("16S-taxonomy.tsv",header=T,row.names=1) ## this is just a list of the ASV ids and the corresponding taxa strings, reordered so that they are listed in the same order as the ASV table
taxmat<-as.matrix(taxa)
all.equal(row.names(taxmat),row.names(otumat)) ## check to make sure they are in the same order
TAX = tax_table(taxmat)
physeq<-phyloseq(OTU,TAX) ## make your phyloseq object, which is basically just a R object that has all of your data for the analyses stored into it!
all.equal(row.names(map_file),sample_names(physeq)) ## LAST CHECK - is everything ordered all correctly?
sampledata<-sample_data(map_file)
mgd<-merge_phyloseq(physeq,sampledata) ## make the final phyloseq object


#### NMDS ####

## Transform to relative abundance, subset to match chem samples, subset by transects, make metadata df's: 
## Mention of relative abundance: https://rstudio-pubs-static.s3.amazonaws.com/496936_2cda5f07e6044d7e9d521893e484a558.html 
mgd_ge5K<-mgd
## mgd_ge5K<-subset_samples(mgd,sample_sums(mgd)>=5000) #filter if low counts
mgd_ge5K_relabund<-transform_sample_counts(mgd_ge5K,function(x)x/sum(x)) ## DOWNLOAD THIS, write.csv(mgd_ge5K_relabund@otu_table, "relative_abundance.csv") ; Calculate relative abundance from ASV count data - this is what you can save and use for other taxonomy analyses.
## write.csv(mgd_ge5K_relabund@otu_table, "relative_abundance.csv") 
mgd_ge5K_relabund.bray<-distance(mgd_ge5K_relabund,"bray") ## calculate distance object using Bray-Curtis dissimilarity
mgd_ge5K_relabund.bray.nmds<-ordinate(mgd_ge5K_relabund,"NMDS",mgd_ge5K_relabund.bray, trymax = 10000)
mgd_ge5K_relabund.bray

mgd_ge5K_relabund.bray.nmds$stress ## Use this line of code if you want to get out stress value

mgd_relabund_map=as(sample_data(mgd_ge5K_relabund),"data.frame") ## Create plain data frame 'map' of sample metadata
sample_tab<-mgd_relabund_map
head(sample_tab)
# sample_tab$NMDS1<-scores(mgd_ge5K_relabund.bray.nmds)[,1] ## add your NMDS values to the sample table! The 'scores' are the coordinates for your samples in NMDS space. Here you are adding the coordinate for NMDS1, so the x coordinate for your sample
# sample_tab$NMDS2<-scores(mgd_ge5K_relabund.bray.nmds)[,2] ## Here you are adding the coordinate for NMDS2, so the y coordinate for your sample

## Work around because above code wouldn't work for some reason - May 5th, 2022
sample_tab$NMDS1<-mgd_ge5K_relabund.bray.nmds$points[,1]
sample_tab$NMDS2<-mgd_ge5K_relabund.bray.nmds$points[,2]

plot(mgd_ge5K_relabund.bray.nmds) 

## Plot using ggplot
NMDS <- ggplot(sample_tab)+
  geom_point(aes(x=NMDS1, y=NMDS2, color=Layer), size=4)+
  theme(text=element_text(size = 20)) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group = Layer, color=Layer),linetype = 2) +
  ylim(-0.75, 0.50) + xlim(-0.75,1.00) +
  theme(legend.position = "none")

NMDS


#### NMDS with Envfit - Chemistry Variables #### 

ord <- mgd_ge5K_relabund.bray.nmds
fit <- envfit(ord, map_file[,1:25], perm = 10000, na.rm=TRUE)
fit


## P-value adjustment - from https://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:p.adjust.envfit 
p.adjust.envfit <- function (x, method = 'bonferroni', n)
{
  x.new <- x
  if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
  if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
  if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
  if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
  if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
  cat ('Adjustment of significance by', method, 'method')
  return (x.new)
}

fit_adjust <- p.adjust.envfit(fit)
fit_adjust

fit1 <- as.data.frame(scores(fit_adjust, display = "vectors"))
fit1 <- cbind(fit1, env.variables = rownames(fit1)) ## and then gives them their names
fit1 <- cbind(fit1, pval = fit_adjust$vectors$pvals) ## add pvalues to dataframe
env.scores.fit <- cbind(fit1, r = fit_adjust$vectors$r)
sig.env.scrs <- subset(env.scores.fit, pval<=0.05) ## add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.dune, r>=0.85) ## subset data to show variables significant at 0.05

library (ggrepel)

## Plot with ggplot
NMDS_Envfit_Chemistry <- NMDS +
  geom_segment(data=sig.env.scrs,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow=arrow(length=unit(0.5,'cm')),color='bisque4',inherit.aes=FALSE) + 
  coord_fixed() +
  geom_text_repel(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  theme(legend.position = "none") +
  ylim(-0.75, 0.6) + xlim(-0.75,1.0)

NMDS_Envfit_Chemistry


#### NMDS with Envfit - Taxa ####

## Made a spreadsheet with my OTU table where I kept the same layout and metadata as map_file 
## but with each column representing an OTU, but OTU code replaced with Phyla/Class of that OTU
## Now going to sum all columns of the same Phyla/Class

## Then need to normalize each row to one 
## Then move forward with envfit 

envfit_phyla <- read.delim("16S-feature-table-_NORMALIZED_forEnvfit_Phyla.txt",header = T)
envfit_class <- read.delim("16S-feature-table-_NORMALIZED_forEnvfit_Class.txt",header = T)

## Collapse all rows of the same phyla/class into one for each phyla/class

envfit_phyla_agg <- aggregate(. ~  Phylum, data = envfit_phyla, sum)
envfit_class_agg <- aggregate(. ~  Phylum, data = envfit_class, sum)

## Transpose the table so that each column is a phylum

envfit_phyla_pivot <- t(envfit_phyla_agg)
envfit_class_pivot <- t(envfit_class_agg)

## Export the dataframe to fix some things in excel
## Fixing ID labels, first column here is assumed to be row labels, making first row into headers,
## Add metadata columns
## Make each row sum to one 

## Export
write.table(envfit_phyla_pivot, file = "envfit_phyla_pivot.txt", sep = "\t")
write.table(envfit_class_pivot, file = "envfit_class_pivot.txt", sep = "\t")

## Import

envfit_phyla_edit <- read.delim("envfit_phyla_edited.txt",header=T)
envfit_class_edit <- read.delim("envfit_class_edited.txt",header=T)


## Do Envfit

options(ggrepel.max.overlaps = Inf)

## P-value adjustments
ord_phyla <- mgd_ge5K_relabund.bray.nmds
fit_phyla <- envfit(ord_phyla, envfit_phyla_edit, perm = 10000, na.rm=TRUE)
fit_phyla

p.adjust.envfit <- function (x, method = 'bonferroni', n)
{
  x.new <- x
  if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
  if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
  if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
  if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
  if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
  cat ('Adjustment of significance by', method, 'method')
  return (x.new)
}

fit_adjust <- p.adjust.envfit(fit_phyla)
fit_adjust

fit1 <- as.data.frame(scores(fit_adjust, display = "vectors"))
fit1 <- cbind(fit1, env.variables = rownames(fit1)) #and then gives them their names
fit1 <- cbind(fit1, pval = fit_adjust$vectors$pvals) # add pvalues to dataframe
env.scores.fit <- cbind(fit1, r = fit_adjust$vectors$r)
sig.env.scrs <- subset(env.scores.fit, pval<=0.05)# add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.dune, r>=0.85) #subset data to show variables significant at 0.05

library(ggrepel)

## Plot 
NMDS_Envfit_Phyla <- NMDS +
  geom_segment(data=sig.env.scrs,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow=arrow(length=unit(0.5,'cm')),color='bisque4',inherit.aes=FALSE) + 
  coord_fixed() +
  geom_text_repel(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  theme(legend.position = "none") +
  ylim(-0.75, 0.50) + xlim(-0.9,1.0)

NMDS_Envfit_Phyla


## Class 
ord_class <- mgd_ge5K_relabund.bray.nmds
fit_class <- envfit(ord_class, envfit_class_edit, perm = 10000, na.rm=TRUE)
fit_class

fit_adjust <- p.adjust.envfit(fit_class)
fit_adjust

fit1 <- as.data.frame(scores(fit_adjust, display = "vectors"))
fit1 <- cbind(fit1, env.variables = rownames(fit1)) ## and then gives them their names
fit1 <- cbind(fit1, pval = fit_adjust$vectors$pvals) ## add pvalues to dataframe
env.scores.fit <- cbind(fit1, r = fit_adjust$vectors$r)
sig.env.scrs <- subset(env.scores.fit, pval<=0.05) ## add pvalues to dataframe
#sig.env.scrs <- subset(env.scores.dune, r>=0.85) ## subset data to show variables significant at 0.05

library(ggrepel)

options(ggrepel.max.overlaps = Inf)

## Plot 
NMDS_Envfit_Class <- NMDS +
  geom_segment(data=sig.env.scrs,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow=arrow(length=unit(0.5,'cm')),color='bisque4',inherit.aes=FALSE) + 
  coord_fixed() +
  geom_text_repel(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  theme(legend.position = "none")

NMDS_Envfit_Class



#### Alpha Diversity Calculations and Histograms ####

## In the code below, you're making a new data frame 'alpha', calculating the four different diversity metrics 
## (Shannon, Simpson, Pielou, and Species Richness), and adding them to that data frame
totus=t(otus)
alpha = as.data.frame(matrix(data = NA, nrow = 43, ncol = 4))
row.names(alpha) = row.names(totus)

alpha[,1] = diversity(totus, index = "shannon")
alpha[,2] = diversity(totus, index = "simpson")
alpha[,3] = diversity(totus, index = "shannon")/log(specnumber(totus)) ## This is divided by the log because vegan defined Shannon's H with log, rather than ln
alpha[,4] = specnumber(totus)

colnames(alpha) = c("Shannon", "Simpson", "Pielou", "Species_Richness")

## Histogram of diversity metrics 
## write.csv(alpha,'sample_alpha_diversity.csv') #this will save the alpha data frame as a csv to your working directory 1
par(mfrow = c(2, 2)) ## make 2x2 environment
## Plot 
hist(alpha[,1], main="Shannon diversity", xlab="", breaks=10)
hist(alpha[,2], main="Simpson diversity", xlab="", breaks=10)
hist(alpha[,3], main="Pielou", xlab="", breaks=15)
hist(alpha[,4], main="Richness", xlab="", breaks=15)

## Adding columns to sample_tab 
map_file$Shannon <- alpha$Shannon
map_file$Species_Richness <- alpha$Species_Richness
## write.csv(map_file,'map_file_alphastats.csv')



#### Alpha Diversity Boxplots ####

## Boxplots w/ Tests of Significance 
## Note I use the Wilcoxon signed-rank test to test significance between groups, don't want to use pairwise t-test because assumes normally distributed data
## https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_nonparametric/bs704_nonparametric4.html
my_comparisons<-list(c('Ash','Organic'),c('Ash','Mineral'),c('Ash','Unburned Mineral'),c('Organic','Mineral'),c('Organic','Unburned Mineral'),c('Mineral','Unburned Mineral')) #list of all variables you're comparing
my_comparisons_nocontrols<-list(c('Ash','Organic'),c('Ash','Mineral'),c('Organic','Mineral'))
my_comparisons_sites<-list(c("Lonesome Wood Lake", "Lonesome Wood Ridge"), c("Lonesome Wood Lake", "Rendezvous"), c("Lonesome Wood Ridge", "Rendezvous"))


## Shannons

## By Site then Soil Layer
ggplot(map_file, aes(x = Layer, y = Shannon)) + geom_boxplot(stat = "boxplot", aes(fill=as.factor(Layer)), outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, 
                     size = 6, hide.ns = TRUE, tip.length = 0.035, bracket.size = 0.5) +
  facet_wrap(~Site) + theme(legend.position = "none") + geom_jitter(color = "black") +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Shannon Diversity Index") + theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"))

## By Soil Layer by Site
ggplot(map_file, aes(x = Site, y = Shannon)) + geom_boxplot(stat = "boxplot", aes(fill=as.factor(Site)), outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons_sites, method = "wilcox.test", label="p.signif", show.legend = FALSE, 
                     size = 6, hide.ns = TRUE, tip.length = 0.035, bracket.size = 0.5) +
  facet_wrap(~Layer) + theme(legend.position = "none") + geom_jitter(color = "black") +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Shannon Diversity Index") + theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"))

# Boxplot with Combined Sites by Soil Layer
ggplot(map_file, aes(x = Layer, y = Shannon)) + geom_boxplot(stat = "boxplot", aes(fill=as.factor(Layer)), outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, 
                     size = 6, hide.ns = TRUE, tip.length = 0.035, bracket.size = 0.5) +
  theme(legend.position = "none") +
  geom_jitter(color = "black") +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Shannon Diversity Index") 


## Species Richness

## By Site then Soil Layer 
ggplot(map_file, aes(x = Layer, y = Species_Richness)) + geom_boxplot(stat = "boxplot", aes(fill=as.factor(Layer)), outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, 
                     size = 6, hide.ns = TRUE, tip.length = 0.035, bracket.size = 0.5) +
  facet_wrap(~Site) + theme(legend.position = "none") + geom_jitter(color = "black") +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Species Richness") + theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"))


## By Soil Layer by Site
ggplot(map_file, aes(x = Site, y = Species_Richness)) + geom_boxplot(stat = "boxplot", aes(fill=as.factor(Site)), outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons_sites, method = "wilcox.test", label="p.signif", show.legend = FALSE, 
                     size = 6, hide.ns = TRUE, tip.length = 0.035, bracket.size = 0.5) +
  facet_wrap(~Layer) + theme(legend.position = "none") + geom_jitter(color = "black") +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Species Richness") + theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"))


# Boxplot with Combined Sites by Soil Layer
ggplot(map_file, aes(x = Layer, y = Species_Richness)) + geom_boxplot(stat = "boxplot", aes(fill=as.factor(Layer)), outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label="p.signif", show.legend = FALSE, 
                     size = 6, hide.ns = TRUE, tip.length = 0.035, bracket.size = 0.5) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  theme(legend.position = "none") +
  geom_jitter(color = "black") +
  ylab("Species Richness") 

 


#### ANOSIM and MRPP ####
## ANOSIM (analysis of similarities) statistically tests whether there are significant differences between two or more groups of samples
## ANOSIM will spit out an R and p-value: https://sites.google.com/site/mb3gustame/hypothesis-tests/anosim

## ANOSIM comparisons, what is compared depends on the way the data is parsed earlier in this document 
anosim(mgd_ge5K_relabund.bray,grouping=map_file$Burn_Severity) ## first input is the bray-curtis dissimilarity matrix made above, second input is the groups of samples you want to compare
## if want to subset, using parsing stuff above and re-do bray-curtis for this test


## MRPP (multi-response permutation procedure) is very similar to ANOSIM
## MRPP is usually used if >2 groups
## What is compared depends on the way the data is parsed earlier in this document 
mrpp(mgd_ge5K_relabund.bray,grouping=map_file$Burn_Severity) ## same inputs as anosim test




#### PERMANOVA ####

## Remember to subset data to answer questions you're looking at 

## Using adonis2 as adonis is deprecated and the use of by=margin allows for order not to matter
adonis2(mgd_ge5K_relabund.bray ~ Site, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Layer, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Code, data = sample_tab, permutations = 999, method = "bray")

adonis2(mgd_ge5K_relabund.bray ~ Site*Layer, data = sample_tab, permutations = 999, method = "bray", by = "margin")
adonis2(mgd_ge5K_relabund.bray ~ Site+Layer, data = sample_tab, permutations = 999, method = "bray", by = "margin")

adonis2(mgd_ge5K_relabund.bray ~ Burn_Severity*Layer, data = sample_tab, permutations = 999, method = "bray", by = "margin")
adonis2(mgd_ge5K_relabund.bray ~ Burn_Severity+Layer, data = sample_tab, permutations = 999, method = "bray", by = "margin")

adonis2(mgd_ge5K_relabund.bray ~ Burn_Severity, data = sample_tab, permutations = 999, method = "bray")

## Conduct p-value adjustments as needed
## P-value adjustments:
library(stats)
p = c(0.001, 0.001, 0.001)

p.adjusted = p.adjust(p, method = "BH")
p.adjusted


## Post hoc test to tell which pairs are sig different
## New method that uses adonis2 (not adonis)

# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Site, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Layer, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Code, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)

pairwise.adonis2(mgd_ge5K_relabund.bray ~ Layer*Site, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)

pairwise.adonis2(mgd_ge5K_relabund.bray ~ Burn_Severity, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)


## P-value adjustments:
library(stats)
p = c(, , , , , , , , , ,
      , )
p.adjusted = p.adjust(p, method = "BH")
p.adjusted



#### Relative Abundance Barplots - Stacked Barcharts #### 

## Phylum

setwd("/Users/juliefowler/OneDrive - Colostate/PhD Work/Studies/Burn Piles/Montana Burn Pile Work/16S/16S_Redo/Relative_Abundance_Taxa")

PivotTable_Phylum<-read.delim("16S_Phylum_PivotTable.txt",header = T,check.names=FALSE)

melted.data.phylum<-melt(PivotTable_Phylum,id.vars=c('Phylum'))

ggplot(melted.data.phylum,aes(x=variable,y=value,fill=Phylum)) +
  geom_bar(stat='identity',position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_light()


## Class

PivotTable_Class<-read.delim("16S_Class_PivotTable.txt",header = T,check.names=FALSE)

melted.data.Class<-melt(PivotTable_Class,id.vars=c('Class'))

ggplot(melted.data.Class,aes(x=variable,y=value,fill=Class)) +
  geom_bar(stat='identity',position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="bottom")



## Only including phyla whose relative abundances are >0.05% relative abundance averaged burned sites (all layers combined) or averaged unburned sites
## Also all normalized to sum to 1 within each site/layer combo

## Phylum - by Sites and Layers and Condition

setwd("/Users/juliefowler/OneDrive - Colostate/PhD Work/Studies/Burn Piles/Montana Burn Pile Work/16S/16S_Redo/Relative_Abundance_Taxa")


## Phyla

PivotTable_Phylum_CutOff<-read.delim("16S__Phylum_Cutoff_0.005_Sites.txt",header = T,check.names=FALSE)
view(PivotTable_Phylum_CutOff)

library(stringr)
library(dplyr)

melted.data.phylum.cutoff<-melt(PivotTable_Phylum_CutOff,id.vars=c('Phylum'))
view(melted.data.phylum.cutoff)
melted.data.phylum.cutoff$Site<-melted.data.phylum.cutoff$variable
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Site = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("LWL")), "Lonesome Wood Lake", melted.data.phylum.cutoff$Site))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Site = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("LWR")), "Lonesome Wood Ridge", melted.data.phylum.cutoff$Site))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Site = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("Rendesvous")), "Rendezvous", melted.data.phylum.cutoff$Site))

melted.data.phylum.cutoff$Layer<-melted.data.phylum.cutoff$variable
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("Ash")), "Ash", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("OrgC")), "Organic", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("Mineral_Burned")), "Mineral", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("mineral_Unburned")), "Unburned Mineral", melted.data.phylum.cutoff$Layer))

melted.data.phylum.cutoff$Layer <- factor(melted.data.phylum.cutoff$Layer, levels = c("Ash", "Organic", "Mineral", "Unburned Mineral"))

ggplot(melted.data.phylum.cutoff,aes(x=Layer,y=value,fill=Phylum)) +
  geom_bar(stat='identity', color = 'black', position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(cols = vars(Site)) +
  theme_classic() + theme(strip.text = element_text(size = 10)) +
  ylab("Value") +
  scale_fill_manual(values = coul)


# Load RColorBrewer
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(4, "Dark2") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(18)
pie(rep(1, length(coul)), col = coul , main="")



## Combined - Phylum
PivotTable_Phylum_CutOff<-read.delim("16S__Phylum_Cutoff_0.005_Combined.txt",header = T,check.names=FALSE)
melted.data.phylum.cutoff<-melt(PivotTable_Phylum_CutOff,id.vars=c('Phylum'))
melted.data.phylum.cutoff$Layer<-melted.data.phylum.cutoff$variable
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("Ash")), "Ash", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("Organic")), "Organic", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("Mineral")), "Mineral", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff <- melted.data.phylum.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.phylum.cutoff$variable, fixed("ControlMineral")), "Unburned Mineral", melted.data.phylum.cutoff$Layer))
melted.data.phylum.cutoff$Layer <- factor(melted.data.phylum.cutoff$Layer, levels = c("Ash", "Organic", "Mineral", "Unburned Mineral"))

ggplot(melted.data.phylum.cutoff,aes(x=Layer,y=value,fill=Phylum)) +
  geom_bar(stat='identity', color = 'black', position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 10)) +
  ylab("Value") +
  theme_bw() +
  scale_fill_manual(values = coul)



# Load RColorBrewer
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(4, "Dark2") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(37)
pie(rep(1, length(coul)), col = coul , main="")


## Combined - Class
PivotTable_Class_CutOff<-read.delim("16S__Class_Cutoff_0.005_Combined.txt",header = T,check.names=FALSE)
melted.data.Class.cutoff<-melt(PivotTable_Class_CutOff,id.vars=c('Class'))
melted.data.Class.cutoff$Layer<-melted.data.Class.cutoff$variable
melted.data.Class.cutoff <- melted.data.Class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.Class.cutoff$variable, fixed("Ash")), "Ash", melted.data.Class.cutoff$Layer))
melted.data.Class.cutoff <- melted.data.Class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.Class.cutoff$variable, fixed("Organic")), "Organic", melted.data.Class.cutoff$Layer))
melted.data.Class.cutoff <- melted.data.Class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.Class.cutoff$variable, fixed("Mineral")), "Mineral", melted.data.Class.cutoff$Layer))
melted.data.Class.cutoff <- melted.data.Class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.Class.cutoff$variable, fixed("ControlMineral")), "Unburned Mineral", melted.data.Class.cutoff$Layer))
melted.data.Class.cutoff$Layer <- factor(melted.data.Class.cutoff$Layer, levels = c("Ash", "Organic", "Mineral", "Unburned Mineral"))

ggplot(melted.data.Class.cutoff,aes(x=Layer,y=value,fill=Class)) +
  geom_bar(stat='identity', color = 'black', position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 10)) +
  ylab("Value") +
  theme_bw() +
  scale_fill_manual(values = coul)



## Class

PivotTable_Class_CutOff<-read.delim("16S__Class_Cutoff_0.005_Sites.txt",header = T,check.names=FALSE)

library(stringr)
library(dplyr)

melted.data.class.cutoff<-melt(PivotTable_Class_CutOff,id.vars=c('Class'))
melted.data.class.cutoff$Site<-melted.data.class.cutoff$variable
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Site = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("LWL")), "Lonesome Wood Lake", melted.data.class.cutoff$Site))
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Site = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("LWR")), "Lonesome Wood Ridge", melted.data.class.cutoff$Site))
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Site = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("Rendesvous")), "Rendezvous", melted.data.class.cutoff$Site))

melted.data.class.cutoff$Layer<-melted.data.class.cutoff$variable
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("Ash")), "Ash", melted.data.class.cutoff$Layer))
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("OrgC")), "Organic", melted.data.class.cutoff$Layer))
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("Mineral_Burned")), "Mineral", melted.data.class.cutoff$Layer))
melted.data.class.cutoff <- melted.data.class.cutoff %>% mutate(Layer = ifelse(str_detect(melted.data.class.cutoff$variable, fixed("mineral_Unburned")), "Unburned Mineral", melted.data.class.cutoff$Layer))

melted.data.class.cutoff$Layer <- factor(melted.data.class.cutoff$Layer, levels = c("Ash", "Organic", "Mineral", "Unburned Mineral"))

ggplot(melted.data.class.cutoff,aes(x=Layer,y=value,fill=Class)) +
  geom_bar(stat='identity', color = 'black', position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(cols = vars(Site)) + 
  theme_classic() + theme(strip.text = element_text(size = 10)) +
  ylab("Value")




#### Relative Abundance Barcharts - Only Pyrophilous Genera of Interest ####

# Load RColorBrewer
library(RColorBrewer)
library(stringr)
library(dplyr)

## Averaged values across all sites - "Average Relative Abundance"
Bacteria_Genera_Pyrophilous_Barchart<-read.delim("16S_Genera_Pyrophilous_Barchart.txt",header = T,check.names=FALSE)
Bacteria_Genera_Pyrophilous_Barchart<-Bacteria_Genera_Pyrophilous_Barchart[1:6,1:5]
Bacteria_Genera_Pyrophilous_Barchart<-melt(Bacteria_Genera_Pyrophilous_Barchart,id.vars=c('Genus'))
Bacteria_Genera_Pyrophilous_Barchart$Layer<-Bacteria_Genera_Pyrophilous_Barchart$variable
Bacteria_Genera_Pyrophilous_Barchart <- Bacteria_Genera_Pyrophilous_Barchart %>% mutate(Layer = ifelse(str_detect(Bacteria_Genera_Pyrophilous_Barchart$variable, fixed("Ash")), "Ash", Bacteria_Genera_Pyrophilous_Barchart$Layer))
Bacteria_Genera_Pyrophilous_Barchart <- Bacteria_Genera_Pyrophilous_Barchart %>% mutate(Layer = ifelse(str_detect(Bacteria_Genera_Pyrophilous_Barchart$variable, fixed("OrgC")), "Organic", Bacteria_Genera_Pyrophilous_Barchart$Layer))
Bacteria_Genera_Pyrophilous_Barchart <- Bacteria_Genera_Pyrophilous_Barchart %>% mutate(Layer = ifelse(str_detect(Bacteria_Genera_Pyrophilous_Barchart$variable, fixed("Mineral")), "Mineral", Bacteria_Genera_Pyrophilous_Barchart$Layer))
Bacteria_Genera_Pyrophilous_Barchart <- Bacteria_Genera_Pyrophilous_Barchart %>% mutate(Layer = ifelse(str_detect(Bacteria_Genera_Pyrophilous_Barchart$variable, fixed("UnburnedMineral")), "Unburned Mineral", Bacteria_Genera_Pyrophilous_Barchart$Layer))
Bacteria_Genera_Pyrophilous_Barchart$Layer <- factor(Bacteria_Genera_Pyrophilous_Barchart$Layer, levels = c("Ash", "Organic", "Mineral", "Unburned Mineral"))

ggplot(Bacteria_Genera_Pyrophilous_Barchart,aes(x=Layer,y=value,fill=Genus)) +
  geom_bar(stat='identity', color = 'black', position='stack') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 10)) +
  ylab("Value") +
  theme_bw() +
  scale_fill_brewer(palette="Dark2")




#### Core Microbiome ####

## Notes:

## Code template from Dr. Emily Bechtold

## Doing these analyses by both taxa (collapse to genus or species, rows combined so each taxa only represented once) and by unique ASV - one file for each of these analyses (genus, species, ASVs)
## This is based on presence/absence, rather than any abundance calculations

## Went forward using 50% ASV cutoff based on "Defining and quantifying the core microbiome: Challenges and prospects" by Neu et al. for the 50%, ASV as that works for our questions 


library(tidyverse) 
library(ggplot2)
library(writexl)
library(ggvenn)


merged <- read_csv("merged_all_genus_ASV.csv")
##metadata <- read_csv("metadata_genus_all_merged.csv")

merged[is.na(merged)] = 0

##metadata$Depth <- as.character(metadata$Depth)

d1 <- merged %>% 
  pivot_longer(2:44, names_to = "Sample", values_to = "Count") %>% 
  separate(Sample, into = c("Site", "Depth", "Burn_Status", "Sample_ID"), sep="_") %>% 
  mutate(Value=Count/Count)


d1[is.na(d1)] <- 0

d1$NewColumn <- 1


## Different Cutoffs to be Considered Core Microbiome

## See what other related papers do - using 100% and 50% as high/low demonstrations and as the two that seem common in the literature 


## 60% in a given combo of site/depth

d50 <- d1 %>% 
  group_by(Site, Depth, Burn_Status, ID) %>% 
  summarise(sum = sum(Value)/sum(NewColumn)) %>% 
  arrange(desc(sum)) %>% 
  ungroup %>% 
  filter(sum >= .5)

d60 <- d1 %>% 
  group_by(Site, Depth, Burn_Status, ID) %>% 
  summarise(sum = sum(Value)/sum(NewColumn)) %>% 
  arrange(desc(sum)) %>% 
  ungroup %>% 
  filter(sum >= .6)  

## 66% in a given combo of site/depth

d66 <- d1 %>% 
  group_by(Site, Depth, Burn_Status, ID) %>% 
  summarise(sum = sum(Value)/sum(NewColumn)) %>% 
  arrange(desc(sum)) %>% 
  ungroup %>% 
  filter(sum >= .66)  

## 80% in a given combo of site/depth

d80 <- d1 %>% 
  group_by(Site, Depth, Burn_Status, ID) %>% 
  summarise(sum = sum(Value)/sum(NewColumn)) %>% 
  arrange(desc(sum)) %>% 
  ungroup %>% 
  filter(sum >= .8)  

## 100% in a given combo of site/depth 

d100 <- d1 %>% 
  group_by(Site, Depth, Burn_Status, ID) %>% 
  summarise(sum = sum(Value)/sum(NewColumn)) %>% 
  arrange(desc(sum)) %>% 
  ungroup %>% 
  filter(sum >= 1)  

## 100% Cutoff - Looking for taxa that are present in at least 100% of a given site/depth combo in 100% of the sites

LWR_Ash_100 <- d100 %>% 
  filter(Site == "LWR" &  Depth == "Ash")

Rend_Ash_100 <- d100 %>% 
  filter(Site == "Rendezvous" &  Depth == "Ash")

Ash_all_100 <- LWR_Ash_100 %>% 
  select(ID) %>% 
  semi_join(Rend_Ash_100) %>% 
  left_join(LWR_Ash_100)

write_xlsx(Ash_all_100,"Ash_all_100_taxa_ASV.xlsx")


LWL_CharredOrgC_100 <- d100 %>% 
  filter(Site == "LWL" &  Depth == "CharredOrgC")


LWR_CharredOrgC_100 <- d100 %>% 
  filter(Site == "LWR" &  Depth == "CharredOrgC")


Rend_CharredOrgC_100 <- d100 %>% 
  filter(Site == "Rendezvous" &  Depth == "CharredOrgC")

CharredOrgC_all_100 <- LWL_CharredOrgC_100 %>% 
  select(ID) %>% 
  semi_join(LWR_CharredOrgC_100) %>% 
  semi_join(Rend_CharredOrgC_100) %>% 
  left_join(LWL_CharredOrgC_100)

write_xlsx(CharredOrgC_all_100,"CharredOrgC_all_100_taxa_ASV.xlsx")


LWL_Mineral_100 <- d100 %>% 
  filter(Site == "LWL" &  Depth == "Mineral")

LWR_Mineral_100 <- d100 %>% 
  filter(Site == "LWR" &  Depth == "Mineral")


Rend_Mineral_100 <- d100 %>% 
  filter(Site == "Rendezvous" &  Depth == "Mineral")


Mineral_all_100 <- LWL_Mineral_100 %>% 
  select(ID) %>% 
  semi_join(LWR_Mineral_100) %>% 
  semi_join(Rend_Mineral_100) %>% 
  left_join(LWL_Mineral_100)

write_xlsx(Mineral_all_100,"Mineral_all_100_taxa_ASV.xlsx")



LWL_ControlMineral_100 <- d100 %>% 
  filter(Site == "LWL" &  Depth == "ControlMineral")

LWR_ControlMineral_100 <- d100 %>% 
  filter(Site == "LWR" &  Depth == "ControlMineral")


Rend_ControlMineral_100 <- d100 %>% 
  filter(Site == "Rendezvous" &  Depth == "ControlMineral")


ControlMineral_all_100 <- LWL_ControlMineral_100 %>% 
  select(ID) %>% 
  semi_join(LWR_ControlMineral_100) %>% 
  semi_join(Rend_ControlMineral_100) %>% 
  left_join(LWL_ControlMineral_100)

write_xlsx(ControlMineral_all_100,"ControlMineral_all_100_taxa_ASV.xlsx")



AshCharMineral_all_100 <- Ash_all_100 %>% 
  select(ID) %>% 
  semi_join(CharredOrgC_all_100) %>% 
  semi_join(Mineral_all_100) %>%
  left_join(Ash_all_100)

write_xlsx(AshCharMineral_all_100,"AshCharMineral_all_100_taxa_ASV.xlsx")


AshChar_all_100 <- Ash_all_100 %>% 
  select(ID) %>% 
  semi_join(CharredOrgC_all_100) %>% 
  left_join(Ash_all_100)

write_xlsx(AshChar_all_100,"AshChar_all_100_taxa_ASV.xlsx")


AshMineral_all_100 <- Ash_all_100 %>% 
  select(ID) %>% 
  semi_join(Mineral_all_100) %>%
  left_join(Ash_all_100)

write_xlsx(AshMineral_all_100,"AshMineral_all_100_taxa_ASV.xlsx")



CharMineral_all_100 <- CharredOrgC_all_100 %>% 
  select(ID) %>% 
  semi_join(Mineral_all_100) %>% 
  left_join(CharredOrgC_all_100)

write_xlsx(CharMineral_all_100,"CharMineral_all_100_taxa_ASV.xlsx")



Mineral_BurnedandUnburned_all_100 <- Mineral_all_100 %>% 
  select(ID) %>% 
  semi_join(ControlMineral_all_100) %>% 
  left_join(Mineral_all_100)

write_xlsx(Mineral_BurnedandUnburned_all_100,"Mineral_BurnedandUnburned_all_100_taxa_ASV.xlsx")



## 50% Cutoff- Looking for taxa that are present in at least 50% of a given site/depth combo in 100% the sites 

LWR_Ash_50 <- d50 %>% 
  filter(Site == "LWR" &  Depth == "Ash")


Rend_Ash_50 <- d50 %>% 
  filter(Site == "Rendezvous" &  Depth == "Ash")


Ash_all_50 <- LWR_Ash_50 %>% 
  select(ID) %>% 
  semi_join(Rend_Ash_50) %>% 
  left_join(LWR_Ash_50)

write_xlsx(Ash_all_50,"Ash_all_50_taxa_ASV.xlsx")


LWL_CharredOrgC_50 <- d50 %>% 
  filter(Site == "LWL" &  Depth == "CharredOrgC")


LWR_CharredOrgC_50 <- d50 %>% 
  filter(Site == "LWR" &  Depth == "CharredOrgC")


Rend_CharredOrgC_50 <- d50 %>% 
  filter(Site == "Rendezvous" &  Depth == "CharredOrgC")


CharredOrgC_all_50 <- LWL_CharredOrgC_50 %>% 
  select(ID) %>% 
  semi_join(LWR_CharredOrgC_50) %>% 
  semi_join(Rend_CharredOrgC_50) %>% 
  left_join(LWL_CharredOrgC_50)

write_xlsx(CharredOrgC_all_50,"CharredOrgC_all_50_taxa_ASV.xlsx")


LWL_Mineral_50 <- d50 %>% 
  filter(Site == "LWL" &  Depth == "Mineral")


LWR_Mineral_50 <- d50 %>% 
  filter(Site == "LWR" &  Depth == "Mineral")


Rend_Mineral_50 <- d50 %>% 
  filter(Site == "Rendezvous" &  Depth == "Mineral")


Mineral_all_50 <- LWL_Mineral_50 %>% 
  select(ID) %>% 
  semi_join(LWR_Mineral_50) %>% 
  semi_join(Rend_Mineral_50) %>% 
  left_join(LWL_Mineral_50)

write_xlsx(Mineral_all_50,"Mineral_all_50_taxa_ASV.xlsx")



LWL_ControlMineral_50 <- d50 %>% 
  filter(Site == "LWL" &  Depth == "ControlMineral")


LWR_ControlMineral_50 <- d50 %>% 
  filter(Site == "LWR" &  Depth == "ControlMineral")


Rend_ControlMineral_50 <- d50 %>% 
  filter(Site == "Rendezvous" &  Depth == "ControlMineral")


ControlMineral_all_50 <- LWL_ControlMineral_50 %>% 
  select(ID) %>% 
  semi_join(LWR_ControlMineral_50) %>% 
  semi_join(Rend_ControlMineral_50) %>% 
  left_join(LWL_ControlMineral_50)

write_xlsx(ControlMineral_all_50,"ControlMineral_all_50_taxa_ASV.xlsx")



AshCharMineral_all_50 <- Ash_all_50 %>% 
  select(ID) %>% 
  semi_join(CharredOrgC_all_50) %>% 
  semi_join(Mineral_all_50) %>%
  left_join(Ash_all_50)

write_xlsx(AshCharMineral_all_50,"AshCharMineral_all_50_taxa_ASV.xlsx")



AshChar_all_50 <- Ash_all_50 %>% 
  select(ID) %>% 
  semi_join(CharredOrgC_all_50) %>% 
  left_join(Ash_all_50)

write_xlsx(AshChar_all_50,"AshChar_all_50_taxa_ASV.xlsx")



AshMineral_all_50 <- Ash_all_50 %>% 
  select(ID) %>% 
  semi_join(Mineral_all_50) %>%
  left_join(Ash_all_50)

write_xlsx(AshMineral_all_50,"AshMineral_all_50_taxa_ASV.xlsx")



CharMineral_all_50 <- CharredOrgC_all_50 %>% 
  select(ID) %>% 
  semi_join(Mineral_all_50) %>% 
  left_join(CharredOrgC_all_50)

write_xlsx(CharMineral_all_50,"CharMineral_all_50_taxa_ASV.xlsx")



Mineral_BurnedandUnburned_all_50 <- Mineral_all_50 %>% 
  select(ID) %>% 
  semi_join(ControlMineral_all_50) %>% 
  left_join(Mineral_all_50)

write_xlsx(Mineral_BurnedandUnburned_all_50,"Mineral_BurnedandUnburned_all_50_taxa_ASV.xlsx")




#### Venn Diagram of Core Microbiome Numbers ####

## install.packages("ggvenn")             
library("ggvenn")

Data_VennDiagram <- read.delim("Input_Data_GGVenn.txt",header = T,check.names=FALSE) 

ggvenn(Data_VennDiagram, c("Ash", "Organic","Mineral","Unburned Mineral"), 
       fill_color = c("#DB735C", "#EFA86E", "#B79667", "#7A6752"),
       text_size = 4,
       set_name_size = 5.5) +
  scale_x_continuous(expand = expansion(mult = 0.5))




#### Dot Plot ####

## Subset of phyla - same as the ones in the subsetted stacked barcharts
## Just deleted the rows with the phyla we don't want (so doesn't sum to 1 or 100) 
## So the numbers are accurate to the overall community 

boxplot_data_16S <- read.delim("boxplot-data-16S-subsetofphyla.txt",header = T,row.names=1,check.names=FALSE)

unique(boxplot_data_16S$Core)
boxplot_data_16S$Core <- factor(boxplot_data_16S$Core, levels = c("Ash", "Organic", "Mineral", "Ash & Organic", "Ash & Mineral",
                                                                  "Ash, Organic, & Mineral", "Ash, Organic, & Mineral; Unburned Mineral", 
                                                                  "Ash & Organic; Unburned Mineral",
                                                                  "Organic & Mineral", "Unburned Mineral", "None"))

unique(boxplot_data_16S$Core_Ash)
boxplot_data_16S$Core_Ash <- factor(boxplot_data_16S$Core_Ash, levels = c("Yes", "No"))
unique(boxplot_data_16S$Core_Organic)
boxplot_data_16S$Core_Organic <- factor(boxplot_data_16S$Core_Organic, levels = c("Yes", "No"))
unique(boxplot_data_16S$Core_Mineral)
boxplot_data_16S$Core_Mineral <- factor(boxplot_data_16S$Core_Mineral, levels = c("Yes", "No"))
unique(boxplot_data_16S$Core_Unburned_Mineral)
boxplot_data_16S$Core_Unburned_Mineral <- factor(boxplot_data_16S$Core_Unburned_Mineral, levels = c("Yes", "No"))

sort(unique(boxplot_data_16S$Phylum))
boxplot_data_16S$Phylum <- factor(boxplot_data_16S$Phylum, levels = c("Unknown",
                                                                      "WS2",
                                                                      "WPS-2",
                                                                      "Verrucomicrobiota",
                                                                      "Thermoplasmatota",
                                                                      "Sumerlaeota",
                                                                      "Spirochaetota",
                                                                      "SAR324_clade(Marine_group_B)",
                                                                      "RCP2-54",
                                                                      "Proteobacteria",
                                                                      "Planctomycetota",
                                                                      "Patescibacteria",
                                                                      "Opalinata",
                                                                      "Nitrospirota",
                                                                      "Nanoarchaeota",
                                                                      "Myxococcota",
                                                                      "Methylomirabilota",
                                                                      "Margulisbacteria",
                                                                      "Latescibacterota",
                                                                      "Hydrogenedentes",
                                                                      "Gemmatimonadota",
                                                                      "Fusobacteriota",
                                                                      "Firmicutes",
                                                                      "Fibrobacterota",
                                                                      "FCPU426",
                                                                      "Entotheonellaeota",
                                                                      "Elusimicrobiota",
                                                                      "Desulfobacterota",
                                                                      "Dependentiae",
                                                                      "Deinococcota",
                                                                      "Cyanobacteria",
                                                                      "Crenarchaeota",
                                                                      "Chloroflexi",
                                                                      "Bdellovibrionota",
                                                                      "Bacteroidota",
                                                                      "Armatimonadota",
                                                                      "Actinobacteriota",
                                                                      "Acidobacteriota",
                                                                      "Abditibacteriota"))


## Ash
boxplot_Ash <- ggplot(boxplot_data_16S, aes(x=Phylum, y=Ash_RA_100)) + 
  geom_jitter(aes(size = Ash_RA_100,color = Core_Ash), width = 0.30, alpha=0.7) +
  coord_flip() +
  scale_color_manual(values = c("#DB735C","#CDCBCB")) +
  scale_size(limits = c(0,20)) + 
  scale_alpha(limits = c(0,20)) +
  ylim(-.1,10) +
  theme(legend.position = "none")


boxplot_Ash


## Organic

boxplot_Organic <- ggplot(boxplot_data_16S, aes(x=Phylum, y=Organic_RA_100)) + 
  geom_jitter(aes(size = Organic_RA_100,color = Core_Organic), width = 0.30, alpha=0.7) +
  coord_flip() +
  scale_color_manual(values = c("#DB735C","#CDCBCB")) +
  scale_size(limits = c(0,20)) + 
  scale_alpha(limits = c(0,20)) +
  ylim(-.1,10) +
  theme(legend.position = "none")

boxplot_Organic


## Mineral

boxplot_Mineral <- ggplot(boxplot_data_16S, aes(x=Phylum, y=Mineral_RA_100)) + 
  geom_jitter(aes(size = Mineral_RA_100,color = Core_Mineral), width = 0.30, alpha=0.7) +
  coord_flip() +
  scale_color_manual(values = c("#DB735C","#CDCBCB")) +
  scale_size(limits = c(0,20)) + 
  scale_alpha(limits = c(0,20)) +
  ylim(-.1,10) +
  theme(legend.position = "none")

boxplot_Mineral


## Unburned Mineral

boxplot_Unburned_Mineral <- ggplot(boxplot_data_16S, aes(x=Phylum, y=Control_Mineral_RA_100)) + 
  geom_jitter(aes(size = Control_Mineral_RA_100,color = Core_Unburned_Mineral), width = 0.30, alpha=0.7) +
  coord_flip() +
  scale_color_manual(values = c("#DB735C","#CDCBCB")) +
  scale_size(limits = c(0,20)) + 
  scale_alpha(limits = c(0,20)) +
  ylim(-.1,10) +
  theme(legend.position = "none")

boxplot_Unburned_Mineral





#### Combined LEfSe/LDA Analysis ####

## Conducted online on the MicrobiomeAnalyst Server
## https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/ModuleView.xhtml 
## Overlaid core microbiome information using our above analysis and Adobe Illustrator



#### 16S Copy Numbers - Scatterplot ####


CopyNumber_Data_forGraph <- read.delim("CopyNumber_Data_forGraph.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Ash <- read.delim("CopyNumber_Data_forGraph_AshOnly.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Organic <- read.delim("CopyNumber_Data_forGraph_OrganicOnly.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Mineral <- read.delim("CopyNumber_Data_forGraph_MineralOnly.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_ControlMineral <- read.delim("CopyNumber_Data_forGraph_ControlMineralOnly.txt",header = T,row.names=1,check.names=FALSE)

## Ordering variable for ggplot
CopyNumber_Data_forGraph_Ash$Core <- factor(CopyNumber_Data_forGraph_Ash$Core, 
    levels = c("Ash", "Ash, Organic", "Ash, Organic, Mineral", "Ash, Organic, Unburned Mineral",
               "Organic", "Organic, Mineral", "Mineral", "Unburned Mineral", "None")) 
CopyNumber_Data_forGraph_Organic$Core <- factor(CopyNumber_Data_forGraph_Organic$Core, 
    levels = c("Ash", "Ash, Organic", "Ash, Organic, Mineral", "Ash, Organic, Unburned Mineral",
               "Organic", "Organic, Mineral", "Mineral", "Unburned Mineral", "None")) 
CopyNumber_Data_forGraph_Mineral$Core <- factor(CopyNumber_Data_forGraph_Mineral$Core, 
    levels = c("Ash", "Ash, Organic", "Ash, Organic, Mineral", "Ash, Organic, Unburned Mineral",
               "Organic", "Organic, Mineral", "Mineral", "Unburned Mineral", "None")) 
CopyNumber_Data_forGraph_ControlMineral$Core <- factor(CopyNumber_Data_forGraph_ControlMineral$Core, 
    levels = c("Ash, Organic", "Ash, Organic, Mineral", "Ash, Organic, Unburned Mineral",
               "Organic", "Organic, Mineral", "Unburned Mineral", "None")) 


library(ggpmisc)

## Ash

a <- mean(CopyNumber_Data_forGraph_Ash$Mean)

Ash_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Ash, aes(x = Ash_RA_100, y = Mean)) + 
  geom_point(aes(size = Core, shape = Core, color = Core)) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Ash Layer") +
  scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
                              '#000000','#000000','#000000','#939799')) +
  scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=a, linetype="longdash", colour = "black") +
  ylim(0,12.5)
  # stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  # stat_regline_equation(label.y = 9, aes(label = ..rr.label..))

Ash_CopyNumber <- Ash_CopyNumber + stat_cor(method = "spearman", 
                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                        label.x = 3)

Ash_CopyNumber


## Organic

b <- mean(CopyNumber_Data_forGraph_Organic$Mean)

Organic_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Organic, aes(x = Organic_RA_100, y = Mean)) + 
  geom_point(aes(size = Core, shape = Core, color = Core)) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Organic Layer") +
  scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
                              '#000000','#000000','#000000','#939799')) +
  scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=b, linetype="longdash", colour = "black") +
  ylim(0,12.5)
  # stat_regline_equation(label.y = 12, aes(label = ..eq.label..)) +
  # stat_regline_equation(label.y = 11, aes(label = ..rr.label..))

Organic_CopyNumber <- Organic_CopyNumber + stat_cor(method = "spearman", 
                                            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                            label.x = 5)

Organic_CopyNumber


## Mineral

c <- mean(CopyNumber_Data_forGraph_Mineral$Mean)

Mineral_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Mineral, aes(x = Mineral_RA_100, y = Mean)) + 
  geom_point(aes(size = Core, shape = Core, color = Core)) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Mineral Layer") +
  scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
                              '#000000','#000000','#000000','#939799')) +
  scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=c, linetype="longdash", colour = "black") +
  ylim(0,12.5)
  # stat_regline_equation(label.y = 12.5, aes(label = ..eq.label..)) +
  # stat_regline_equation(label.y = 11.5, aes(label = ..rr.label..))

Mineral_CopyNumber <- Mineral_CopyNumber + stat_cor(method = "spearman", 
                                                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                    label.x = 4)

Mineral_CopyNumber


## Control Mineral

d <- mean(CopyNumber_Data_forGraph_ControlMineral$Mean)

ControlMineral_CopyNumber <- ggplot(CopyNumber_Data_forGraph_ControlMineral, aes(x = ControlMineral_RA_100, y = Mean)) + 
  geom_point(aes(size = Core, shape = Core, color = Core)) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Control Mineral Layer") +
  scale_shape_manual(values=c(16,17,18,19,8,9,1)) +
  scale_color_manual(values=c('#000000','#000000','#000000','#000000',
                              '#000000','#000000','#939799')) +
  scale_size_manual(values=c(3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=d, linetype="longdash", colour = "black") +
  ylim(0,12.5)
  # stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  # stat_regline_equation(label.y = 9, aes(label = ..rr.label..))

ControlMineral_CopyNumber <- ControlMineral_CopyNumber + stat_cor(method = "spearman", 
                                                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                    label.x = 1.25)

ControlMineral_CopyNumber


## ggarrange

library(ggpubr)

figure <- ggarrange(Ash_CopyNumber + rremove("ylab"), Organic_CopyNumber + rremove("ylab"), 
                    Mineral_CopyNumber + rremove("ylab"), ControlMineral_CopyNumber + rremove("ylab"),
          labels = c("A. Ash", "B. Organic", "C. Mineral", "D. Control Mineral"),
          ncol = 2, nrow = 2, align = "hv",
          common.legend = TRUE, legend = "bottom")

annotate_figure(figure, left = text_grob("Mean 16S Copy Number at Lowest ID'd Taxanomic Level", 
                                        rot = 90, vjust = 1),
                        bottom = text_grob("% Relative Abundance of ASV in Given Layer", 
                                          rot = 0, vjust = 1))






## Version #2 - color code points by phyla, encircle points core to that layer ##


CopyNumber_Data_forGraph_Ash <- read.delim("CopyNumber_Data_forGraph_AshOnly.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Organic <- read.delim("CopyNumber_Data_forGraph_OrganicOnly.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Mineral <- read.delim("CopyNumber_Data_forGraph_MineralOnly.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_ControlMineral <- read.delim("CopyNumber_Data_forGraph_ControlMineralOnly.txt",header = T,row.names=1,check.names=FALSE)


## Ordering phylum variable for ggplot
CopyNumber_Data_forGraph_Ash$Phylum <- factor(CopyNumber_Data_forGraph_Ash$Phylum, 
    levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
               "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
               "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
               "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
               "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
               "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
               "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 
CopyNumber_Data_forGraph_Organic$Phylum <- factor(CopyNumber_Data_forGraph_Organic$Phylum, 
    levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
               "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
               "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
               "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
               "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
               "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
               "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 
CopyNumber_Data_forGraph_Mineral$Phylum <- factor(CopyNumber_Data_forGraph_Mineral$Phylum, 
    levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
               "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
               "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
               "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
               "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
               "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
               "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 
CopyNumber_Data_forGraph_ControlMineral$Phylum <- factor(CopyNumber_Data_forGraph_ControlMineral$Phylum, 
    levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
               "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
               "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
               "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
               "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
               "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
               "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 





library(ggpmisc)


# Load RColorBrewer
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(8, "Dark2") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(22)
pie(rep(1, length(coul)), col = coul , main="")


## Ash

a <- mean(CopyNumber_Data_forGraph_Ash$Mean)
CopyNumber_Data_forGraph_Ash$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_Ash$Core_ThisLayer)

Ash_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Ash, aes(x = Ash_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_Ash[CopyNumber_Data_forGraph_Ash$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=2, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Ash Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=a, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 9, aes(label = ..rr.label..))

Ash_CopyNumber <- Ash_CopyNumber + stat_cor(method = "spearman", 
                                            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                            label.x = 3)

Ash_CopyNumber


## Organic

b <- mean(CopyNumber_Data_forGraph_Organic$Mean)
CopyNumber_Data_forGraph_Organic$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_Organic$Core_ThisLayer)

Organic_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Organic, aes(x = Organic_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_Organic[CopyNumber_Data_forGraph_Organic$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=2, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Organic Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=b, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 12, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 11, aes(label = ..rr.label..))

Organic_CopyNumber <- Organic_CopyNumber + stat_cor(method = "spearman", 
                                                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                    label.x = 5)

Organic_CopyNumber


## Mineral

c <- mean(CopyNumber_Data_forGraph_Mineral$Mean)
CopyNumber_Data_forGraph_Mineral$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_Mineral$Core_ThisLayer)

Mineral_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Mineral, aes(x = Mineral_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_Mineral[CopyNumber_Data_forGraph_Mineral$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=2, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Mineral Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=c, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 12.5, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 11.5, aes(label = ..rr.label..))

Mineral_CopyNumber <- Mineral_CopyNumber + stat_cor(method = "spearman", 
                                                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                    label.x = 4)

Mineral_CopyNumber


## Control Mineral

d <- mean(CopyNumber_Data_forGraph_ControlMineral$Mean)
CopyNumber_Data_forGraph_ControlMineral$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_ControlMineral$Core_ThisLayer)

ControlMineral_CopyNumber <- ggplot(CopyNumber_Data_forGraph_ControlMineral, aes(x = ControlMineral_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_ControlMineral[CopyNumber_Data_forGraph_ControlMineral$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=2, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in ControlMineral Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=d, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 9, aes(label = ..rr.label..))

ControlMineral_CopyNumber <- ControlMineral_CopyNumber + stat_cor(method = "spearman", 
                                                                  aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                                  label.x = 1.25)

ControlMineral_CopyNumber


## ggarrange

library(ggpubr)

figure <- ggarrange(Ash_CopyNumber + rremove("ylab"), Organic_CopyNumber + rremove("ylab"), 
                    Mineral_CopyNumber + rremove("ylab"), ControlMineral_CopyNumber + rremove("ylab"),
                    labels = c("A. Ash", "B. Organic", "C. Mineral", "D. Control Mineral"),
                    ncol = 2, nrow = 2, align = "hv",
                    common.legend = TRUE, legend = "bottom")

annotate_figure(figure, left = text_grob("Mean 16S Copy Number at Lowest ID'd Taxanomic Level", 
                                         rot = 90, vjust = 1),
                bottom = text_grob("% Relative Abundance of ASV in Given Layer", 
                                   rot = 0, vjust = 1))




## Version #3 - same as version 2, but all ASVs less than 0.05% filtered out ##

CopyNumber_Data_forGraph_Ash <- read.delim("CopyNumber_Data_forGraph_AshOnly_Cutoff0.0005.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Organic <- read.delim("CopyNumber_Data_forGraph_OrganicOnly_Cutoff0.0005.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_Mineral <- read.delim("CopyNumber_Data_forGraph_MineralOnly_Cutoff0.0005.txt",header = T,row.names=1,check.names=FALSE)
CopyNumber_Data_forGraph_ControlMineral <- read.delim("CopyNumber_Data_forGraph_ControlMineralOnly_Cutoff0.0005.txt",header = T,row.names=1,check.names=FALSE)


## Ordering phylum variable for ggplot
CopyNumber_Data_forGraph_Ash$Phylum <- factor(CopyNumber_Data_forGraph_Ash$Phylum, 
                                              levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
                                                         "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                         "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
                                                         "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
                                                         "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
                                                         "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
                                                         "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 
CopyNumber_Data_forGraph_Organic$Phylum <- factor(CopyNumber_Data_forGraph_Organic$Phylum, 
                                                  levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
                                                             "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                             "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
                                                             "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
                                                             "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
                                                             "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
                                                             "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 
CopyNumber_Data_forGraph_Mineral$Phylum <- factor(CopyNumber_Data_forGraph_Mineral$Phylum, 
                                                  levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
                                                             "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                             "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
                                                             "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
                                                             "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
                                                             "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
                                                             "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 
CopyNumber_Data_forGraph_ControlMineral$Phylum <- factor(CopyNumber_Data_forGraph_ControlMineral$Phylum, 
                                                         levels = c("Abditibacteriota", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota",
                                                                    "Bdellovibrionota", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Deinococcota", "Dependentiae",
                                                                    "Desulfobacterota", "Elusimicrobiota", "Entotheonellaeota", "FCPU426", "Fibrobacterota", "Firmicutes",
                                                                    "Fusobacteriota", "Gemmatimonadota", "Hydrogenedentes", "Latescibacterota", "Margulisbacteria", "Methylomirabilota",
                                                                    "Myxococcota", "Nanoarchaeota", "Nitrospirota", "Opalinata", "Patescibacteria", "Planctomycetota",
                                                                    "Proteobacteria", "RCP2-54", "SAR324_clade(Marine_group_B)", "Spirochaetota", "Sumerlaeota",
                                                                    "Thermoplasmatota", "Verrucomicrobiota", "WPS-2","WS2")) 





library(ggpmisc)


# Load RColorBrewer
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(8, "Dark2") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(22)
pie(rep(1, length(coul)), col = coul , main="")


## Ash

a <- mean(CopyNumber_Data_forGraph_Ash$Mean)
CopyNumber_Data_forGraph_Ash$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_Ash$Core_ThisLayer)

Ash_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Ash, aes(x = Ash_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_Ash[CopyNumber_Data_forGraph_Ash$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=4, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Ash Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=a, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 9, aes(label = ..rr.label..))

Ash_CopyNumber <- Ash_CopyNumber + stat_cor(method = "spearman", 
                                            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                            label.x = 3)

Ash_CopyNumber


## Organic

b <- mean(CopyNumber_Data_forGraph_Organic$Mean)
CopyNumber_Data_forGraph_Organic$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_Organic$Core_ThisLayer)

Organic_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Organic, aes(x = Organic_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_Organic[CopyNumber_Data_forGraph_Organic$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=4, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Organic Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=b, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 12, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 11, aes(label = ..rr.label..))

Organic_CopyNumber <- Organic_CopyNumber + stat_cor(method = "spearman", 
                                                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                    label.x = 5)

Organic_CopyNumber


## Mineral

c <- mean(CopyNumber_Data_forGraph_Mineral$Mean)
CopyNumber_Data_forGraph_Mineral$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_Mineral$Core_ThisLayer)

Mineral_CopyNumber <- ggplot(CopyNumber_Data_forGraph_Mineral, aes(x = Mineral_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_Mineral[CopyNumber_Data_forGraph_Mineral$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=4, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in Mineral Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=c, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 12.5, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 11.5, aes(label = ..rr.label..))

Mineral_CopyNumber <- Mineral_CopyNumber + stat_cor(method = "spearman", 
                                                    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                    label.x = 4)

Mineral_CopyNumber


## Control Mineral

d <- mean(CopyNumber_Data_forGraph_ControlMineral$Mean)
CopyNumber_Data_forGraph_ControlMineral$Core_ThisLayer <- as.factor(CopyNumber_Data_forGraph_ControlMineral$Core_ThisLayer)

ControlMineral_CopyNumber <- ggplot(CopyNumber_Data_forGraph_ControlMineral, aes(x = ControlMineral_RA_100, y = Mean)) + 
  geom_point(aes(color = Phylum), size = 2) +
  geom_point(data=CopyNumber_Data_forGraph_ControlMineral[CopyNumber_Data_forGraph_ControlMineral$Core_ThisLayer %in% 'Yes', ],
             pch=21, fill=NA, size=4, colour="black", stroke=1) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme(axis.title.x = element_text(vjust = -0.3)) +
  theme(axis.title.y = element_text(vjust = 1.5)) +
  ylab("Mean 16S Copy Number at Lowest ID'd Taxanomic Level") +
  xlab("% Relative Abundance of ASV in ControlMineral Layer") +
  # scale_shape_manual(values=c(15,16,17,18,19,8,7,9,1)) +
  # scale_color_manual(values=c('#000000','#000000','#000000','#000000','#000000',
  #                             '#000000','#000000','#000000','#939799')) +
  # scale_size_manual(values=c(3,3,3,3,3,3,3,3,2)) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) +
  geom_hline(yintercept=d, linetype="longdash", colour = "black") +
  ylim(0,12.5) +
  scale_color_manual(values = coul)
# stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
# stat_regline_equation(label.y = 9, aes(label = ..rr.label..))

ControlMineral_CopyNumber <- ControlMineral_CopyNumber + stat_cor(method = "spearman", 
                                                                  aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                                                  label.x = 1.25)

ControlMineral_CopyNumber


## ggarrange

library(ggpubr)

figure <- ggarrange(Ash_CopyNumber + rremove("ylab"), Organic_CopyNumber + rremove("ylab"), 
                    Mineral_CopyNumber + rremove("ylab"), ControlMineral_CopyNumber + rremove("ylab"),
                    labels = c("A. Ash", "B. Organic", "C. Mineral", "D. Control Mineral"),
                    ncol = 2, nrow = 2, align = "hv",
                    common.legend = TRUE, legend = "bottom")

annotate_figure(figure, left = text_grob("Mean 16S Copy Number at Lowest ID'd Taxanomic Level", 
                                         rot = 90, vjust = 1),
                bottom = text_grob("% Relative Abundance of ASV in Given Layer", 
                                   rot = 0, vjust = 1))


