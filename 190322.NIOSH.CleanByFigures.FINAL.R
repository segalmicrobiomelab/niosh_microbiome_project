##Scrits for NIOSH New Badge (MSQ 44 TO 50)
#LNS/bck
#09/05/18

#Load required libraries
library(phyloseq)
library('ape')
library("plyr")
library("gplots")
library("RColorBrewer")
library("d3heatmap")
library("vegan")
library("Heatplus")
library("ade4")
library("RAM")
library("ggplot2")
library("mvabund")
library("reshape")
library("dplyr")
library("plotly")
library("devtools")
library("pairwiseAdonis")
library('DESeq2')
library("ggrepel")
library(data.table)
library(scales)
library(tidyr)
library(forcats)


#Code to load and save files 
load(file="190306.NIOSH.RData")
save.image(file="190306.NIOSH.RData")


####Set up code 
{
##Load the files needed
file = "/Users/segal-lab/Dropbox/LNS_Bianca/NIOSH.Project/Second.Part/NIOSH.New.Badge/Sequence.Data/Merged.Open/5_otus/otu_table_mc2_w_tax.filter.biom"
map = "/Users/segal-lab/Dropbox/LNS_Bianca/NIOSH.Project/Second.Part/NIOSH.New.Badge/Sequence.Data/Merged.Open/Map/Map.Combined.NIOSH.map.filter.sort.b.txt"

#leave A1 empty (the # make it not work) in the mapping file
treefile = "/Users/segal-lab/Dropbox/LNS_Bianca/NIOSH.Project/Second.Part/NIOSH.New.Badge/Sequence.Data/Merged.Open/5_otus/rep_set.tre"


# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))


###Merge
lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)


#Give a colnames to separate different taxonomic levels
colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Load the tree file
tree.obj = import_qiime(treefilename = treefile)

# Now merge the three separate phyloseq objects into a single object
otu.table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)

# Remove taxa with 0 abundance
otu.table = subset_taxa(otu.table, rowSums(otu_table(otu.table)) != 0)

otu_table(otu.table)

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)


rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))

#We will use this to get al of specific samples
sample_data(otu.relative.table)$Sample_Type_c_code

#We will use this as a code for each location 
sample_data(otu.relative.table)$Sample_Type_d_fluid_change_code


# Create phyllum and order tables (do it after normalization and out of the relative table) - we wil ony do Genus for now, but you can do others 
#Phylum.rel.table = tax_glom(otu.relative.table, taxrank = "Phylum")
#Class.rel.table = tax_glom(otu.relative.table, taxrank = "Class")
#Order.rel.table = tax_glom(otu.relative.table, taxrank = "Order")
#Family.rel.table = tax_glom(otu.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(otu.relative.table, taxrank = "Genus")
#OTU.rel.table = tax_glom(otu.relative.table, taxrank = "OTU")


}

####### MAIN FIGURES 

###Figure 1A
## Create PCOA of Skin, Oral, Nasal, Lung, and MWF 
{
#select just the right samples from the genus table
lung.mwf.human.otu.rel.table = subset_samples(Genus.rel.table, Sample_Type_d_fluid_change_code == 101 | Sample_Type_d_fluid_change_code == 102 | Sample_Type_d1_c == "MWF.NP" | Sample_Type_d1_c == "MWF.P" | Sample_Type_c_code == 20 | Sample_Type_c_code == 21 | Sample_Type_c_code == 22)
	
	
#this code creates a new variable based on how we'd like to name and identify the groups in out plotting later
#first get the closest variable 
var8 = as.character(get_variable(lung.mwf.human.otu.rel.table,"Sample_Type_d1_c"))
#now rename as you would like
var8[sample_data(lung.mwf.human.otu.rel.table)$Description == "Tissue_case"] <- "Tissue.case"
var8[sample_data(lung.mwf.human.otu.rel.table)$Description == "Tissue_control"] <- "Tissue.control"
var8[sample_data(lung.mwf.human.otu.rel.table)$Sample_Type_d_fluid_change_c == "Nasal.Machine.shop"] <- "Nasal.Machine.shop"
var8[sample_data(lung.mwf.human.otu.rel.table)$Sample_Type_d_fluid_change_c == "Nasal.Assembly"] <- "Nasal.Assembly"
var8[sample_data(lung.mwf.human.otu.rel.table)$Sample_Type_d_fluid_change_c == "Oral.Machine.Shop"] <- "Oral.Machine.Shop"
var8[sample_data(lung.mwf.human.otu.rel.table)$Sample_Type_d_fluid_change_c == "Oral.Assembly"] <- "Oral.Assembly"
var8[sample_data(lung.mwf.human.otu.rel.table)$Sample_Type_d_fluid_change_c == "Skin.Machine.Shop"] <- "Skin.Machine.Shop"
var8[sample_data(lung.mwf.human.otu.rel.table)$Sample_Type_d_fluid_change_c == "Skin.Assembly"] <- "Skin.Assembly"

#check names 
var8

#add the new variable back to the sample data of the phyloseq object, it will be called NewVar
sample_data(lung.mwf.human.otu.rel.table)$NewVar <- mapply(paste0, var8, collapse="_")
merge_samples(lung.mwf.human.otu.rel.table, "NewVar")

### Calculate Bray distance 
lung.mwf.human.bray.dist = distance(lung.mwf.human.otu.rel.table, method="bray")
lung.mwf.human.bray.pco = dudi.pco(cailliez(lung.mwf.human.bray.dist))

##Plot PCoA
pdf(file="Lung.MWF.Human.PCoA.Bray.pdf", width=10, height=10)
s.class(lung.mwf.human.bray.pco $li, interaction(sample_data(lung.mwf.human.otu.rel.table)$NewVar), col=c("grey", "grey", "red", "red", "red", "blue", "blue", "blue", "purple", "purple", "purple", "green", "green"))
dev.off()

# Calculate ordinate bray distance
ordinate.lung.mwf.human.bray.dist <- ordinate(lung.mwf.human.otu.rel.table, method="PCoA", distance="bray")

# Scatter plot , visualize to get axes #Axis1 [42.2%] #Axis2[9.1%]
plot_ordination(lung.mwf.human.otu.rel.table, ordinate.lung.mwf.human.bray.dist) 

#calculate pvalue from anova 
adonis(lung.mwf.human.bray.dist ~ NewVar, data=data.frame(sample_data(lung.mwf.human.otu.rel.table)))
#Result Call:
adonis(formula = lung.mwf.human.bray.dist ~ NewVar, data = data.frame(sample_data(lung.mwf.human.otu.rel.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
NewVar     13    139.04 10.6954   67.22 0.48203  0.001 ***
Residuals 939    149.41  0.1591         0.51797           
Total     952    288.44                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##Convert to ggplot fore more customizable plotting 
{
#create a new null variable to use for clustering - we want to plot by mucosa and color by work site so we need a clustering variable (NewVar2) and a dot color var (dots)
NewVar2 <- NULL
#Create a clustering variable from the one we already made earlier
NewVar2 = as.character(get_variable(lung.mwf.human.otu.rel.table,"NewVar"))
#check the variables, we need to change tissue case and control to just Tissue, etc. 
NewVar2 <- replace(NewVar2, NewVar2 == "Tissue.case", "Tissue")
NewVar2 <- replace(NewVar2, NewVar2 == "Tissue.control", "Tissue")

NewVar2 <- replace(NewVar2, NewVar2 == "Skin.Machine.Shop", "Skin")
NewVar2 <- replace(NewVar2, NewVar2 == "Skin.Assembly", "Skin")
NewVar2 <- replace(NewVar2, NewVar2 == "Skin.Administration", "Skin")

NewVar2 <- replace(NewVar2, NewVar2 == "Nasal.Machine.shop", "Nasal")
NewVar2 <- replace(NewVar2, NewVar2 == "Nasal.Assembly", "Nasal")
NewVar2 <- replace(NewVar2, NewVar2 == "Nasal.Administration", "Nasal")

NewVar2 <- replace(NewVar2, NewVar2 == "Oral.Machine.Shop", "Oral")
NewVar2 <- replace(NewVar2, NewVar2 == "Oral.Assembly", "Oral")
NewVar2 <- replace(NewVar2, NewVar2 == "Oral.Administration", "Oral")

NewVar2 <- replace(NewVar2, NewVar2 == "MWF.NP", "MWF")
NewVar2 <- replace(NewVar2, NewVar2 == "MWF.P", "MWF")


#NewVar has locations for dot plot
#NewVar2 fewer has groups for clustering
sample_data(lung.mwf.human.otu.rel.table)$NewVar2 <- NewVar2
#check the groups you would like to cluster are accurate
unique(sample_data(lung.mwf.human.otu.rel.table)$NewVar2)

#Now I am going to add a column of 'location' that we will use to color the dots. In this case we will use machine shop v assembly v admin
dots = as.character(get_variable(lung.mwf.human.otu.rel.table,"NewVar"))
#check th variables, we need to change
dots <- replace(dots, dots == "Tissue.case", "Case")
dots <- replace(dots, dots == "Tissue.control", "Control")

dots <- replace(dots, dots == "Skin.Machine.Shop", "Machine Shop")
dots <- replace(dots, dots == "Nasal.Machine.shop", "Machine Shop")
dots <- replace(dots, dots == "Oral.Machine.Shop", "Machine Shop")

dots <- replace(dots, dots == "Skin.Assembly", "Assembly")
dots <- replace(dots, dots == "Nasal.Assembly", "Assembly")
dots <- replace(dots, dots == "Oral.Assembly", "Assembly")

dots <- replace(dots, dots == "Skin.Administration", "Administration")
dots <- replace(dots, dots == "Nasal.Administration", "Administration")
dots <- replace(dots, dots == "Oral.Administration", "Administration")
#add dots to the sample data
sample_data(lung.mwf.human.otu.rel.table)$dots <- dots



# build ggplot dataframe with points (x,y), corresponding groups (cluster), and location for dots
gg.all <- data.frame(cluster=factor(sample_data(lung.mwf.human.otu.rel.table)$NewVar2), x2= lung.mwf.human.bray.pco$li$A1, y2= lung.mwf.human.bray.pco$li$A2, dots = factor(sample_data(lung.mwf.human.otu.rel.table)$dots))
# calculate group centroid locations
centroids.all <- aggregate(cbind(x2,y2)~cluster,data= gg.all,mean)
# merge centroid locations into ggplot dataframe
gg.all <- merge(gg.all, centroids.all,by="cluster",suffixes=c("",".centroid2"))


# plot with the lines to the points from the centroid
pdf(file="AllSamples_NoAir_Multicolor_ggplot.pdf", height = 18, width = 20) 
ggplot(gg.all, aes(x = gg.all $x2, y = gg.all $y2, color =cluster)) +
  #sets background
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #sets centroid points and draws lines to data points
  geom_point(data= centroids.all, aes(x=x2, y=y2), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
  #sets point size 
  geom_point(aes(x=x2,y=y2, color = dots), size = 2) +
  
  #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal (we did not use this)
  #stat_ellipse(type = "t") +

  #this sets the color scale all color inputs 
  scale_color_manual(values = c(
    	"Assembly" = "orange",
    	"Machine Shop" = "purple", 
    	"Administration"="palegreen", 
    	"Control"="blue", 
    	"Case"="red", 
    	"MWF.NP"="cyan", 
    	"MWF.P"="darkcyan",
    	"Oral" = "grey", 
    	"Skin"="grey",
    	"Nasal"="grey",
    	"MWF"="purple",
    	"Tissue"="blue")) +
  #labels centroids 
  geom_label(data = centroids.all, aes (x=x2, y=y2, label = cluster))
  
  dev.off()

}
}

### Figure 1 B-E - MWF vs. Human samples 
#This figure was developed in prism but the data files were processed in R as follows
{	

### MWF v Skin 
{
#create otu table of just the required samples
Skin.MWF.otu.rel.table <- subset_samples(Genus.rel.table, Sample_Type_d1_code =="25" | Sample_Type_d1_code == "26" |  Sample_Type_c_code == 22)
#calculate distance by bray
Skin.MWF.Bray.dist = distance(Skin.MWF.otu.rel.table, method="bray")

#melt the distance matrix into columns with each sample site and distance
b1 <- melt(as.matrix(Skin.MWF.Bray.dist))

#Then need to remove self distances and duplicated distances
p1    <- t(apply(b1[,c(1,2)],1,FUN=sort))
rmv11 <- which(p1[,1] == p1[,2])

p1    <- paste(p1[,1],p1[,2],sep="|")
rmv21 <- which(duplicated(p1))

#establish new data frame that removes those values 
b1.df   <- b1[-c(rmv11,rmv21),] 
head(b1.df)

##Now we need to replace variable rows with topo 
#set up new data frame
new1.df <- NULL
new1.df <- b1.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new1.df[] <- lapply(b1.df, function(x) sample_data(Skin.MWF.otu.rel.table)$Sample_Type_d1[match(x, rownames(sample_data(Skin.MWF.otu.rel.table)))])

head(new1.df)

#create two lists of the group variable 
topo.var11 <- new1.df[,1]
topo.var21 <-new1.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b1.var.df <- cbind(b1.df, topo.var11, topo.var21)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b1.var.df <- b1.var.df
#set row names to re-zero
rownames(btw.b1.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b1.var.df)) {
	if (btw.b1.var.df$topo.var11[i] == btw.b1.var.df$topo.var21[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b1.var.df <- btw.b1.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b1.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw1 = paste(btw.b1.var.df$topo.var11, "to", btw.b1.var.df$topo.var21)
new.cat.btw1.df <- data.frame(btw.b1.var.df$topo.var11, btw.b1.var.df$topo.var21)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw1.df, 1, sort))
unique.new.cat.btw1 <- unique(new.cat.btw1.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw1) <- NULL
rownames(unique.new.cat.btw1) <- NULL
unique.new.cat.btw1 <- paste(unique.new.cat.btw1[,1], "to", unique.new.cat.btw1[,2])


#create new data frame 
clean.btw.b1.var.df <- btw.b1.var.df

#reset row names
rownames(clean.btw.b1.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b1.var.df)){
	if (paste(clean.btw.b1.var.df$topo.var21[i], "to", clean.btw.b1.var.df$topo.var11[i]) %in% unique.new.cat.btw1) {
		clean.btw.b1.var.df$topo.var11[i] <- btw.b1.var.df$topo.var21[i]
		clean.btw.b1.var.df$topo.var21[i] <- btw.b1.var.df$topo.var11[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw1.clean = paste(clean.btw.b1.var.df$topo.var11, "to", clean.btw.b1.var.df$topo.var21)

#confirm permutations 
unique(new.cat.btw1.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
inter.skin.mwf.data.bray <- data.frame(Sample1 = clean.btw.b1.var.df$topo.var11, Sample2 = clean.btw.b1.var.df$topo.var21, Category = new.cat.btw1.clean, Distance = as.numeric(clean.btw.b1.var.df$value))

#because the above also plots between skin and between mwf categories want to subset out some 
MWF.Skin.Subset= inter.skin.mwf.data.bray[c(which(inter.skin.mwf.data.bray$Category == "MWF.P to Skin.Machine.Shop"), which(inter.skin.mwf.data.bray$Category == "MWF.NP to Skin.Machine.Shop"), which(inter.skin.mwf.data.bray$Category == "MWF.NP to Skin.Assembly"), which(inter.skin.mwf.data.bray$Category == "MWF.P to Skin.Assembly"), which(inter.skin.mwf.data.bray$Category == "MWF.NP to Skin.Administration"),  which(inter.skin.mwf.data.bray $Category == "MWF.P to Skin.Administration")),]

head(MWF.Skin.Subset)
#save as txt file to take into prism
write.csv(MWF.Skin.Subset, file = "Inter_Bray_Skin_MWF.txt")
}

### MWF v Oral 
{
#create otu table ot just samples we'd like
Oral.MWF.otu.rel.table <- subset_samples(Genus.rel.table, Sample_Type_d1_code == "25" | Sample_Type_d1_code =="26"| Sample_Type_c_code == 21)
#calculate bray distance
Oral.MWF.Bray.dist = distance(Oral.MWF.otu.rel.table, method="bray")

#melt the distance matrix into columns with each sample site and distance
b1 <- melt(as.matrix(Oral.MWF.Bray.dist))
head(b1)

#Then need to remove self distances and duplicated distances
p1    <- t(apply(b1[,c(1,2)],1,FUN=sort))
rmv11 <- which(p1[,1] == p1[,2])

p1    <- paste(p1[,1],p1[,2],sep="|")
rmv21 <- which(duplicated(p1))

#establish new data frame that removes those values 
b1.df   <- b1[-c(rmv11,rmv21),] 
head(b1.df)

##Now we need to replace variable rows with topo 
#set up new data frame
new1.df <- b1.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new1.df[] <- lapply(b1.df, function(x) sample_data(Oral.MWF.otu.rel.table)$Sample_Type_d1[match(x, rownames(sample_data(Oral.MWF.otu.rel.table)))])

head(new.df)

#create two lists of the group variable 
topo.var11 <- new1.df[,1]
topo.var21 <-new1.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b1.var.df <- cbind(b1.df, topo.var11, topo.var21)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b1.var.df <- b1.var.df
#set row names to re-zero
rownames(btw.b1.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b1.var.df)) {
	if (btw.b1.var.df$topo.var11[i] == btw.b1.var.df$topo.var21[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b1.var.df <- btw.b1.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b1.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw1 = paste(btw.b1.var.df$topo.var11, "to", btw.b1.var.df$topo.var21)
new.cat.btw1.df <- data.frame(btw.b1.var.df$topo.var11, btw.b1.var.df$topo.var21)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw1.df, 1, sort))
unique.new.cat.btw1 <- unique(new.cat.btw1.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw1) <- NULL
rownames(unique.new.cat.btw1) <- NULL
unique.new.cat.btw1 <- paste(unique.new.cat.btw1[,1], "to", unique.new.cat.btw1[,2])


#create new data frame 
clean.btw.b1.var.df <- btw.b1.var.df

#reset row names
rownames(clean.btw.b1.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b1.var.df)){
	if (paste(clean.btw.b1.var.df$topo.var21[i], "to", clean.btw.b1.var.df$topo.var11[i]) %in% unique.new.cat.btw1) {
		clean.btw.b1.var.df$topo.var11[i] <- btw.b1.var.df$topo.var21[i]
		clean.btw.b1.var.df$topo.var21[i] <- btw.b1.var.df$topo.var11[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw1.clean = paste(clean.btw.b1.var.df$topo.var11, "to", clean.btw.b1.var.df$topo.var21)

#confirm permutations 
unique(new.cat.btw1.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
inter.oral.mwf.data.bray <- data.frame(Sample1 = clean.btw.b1.var.df$topo.var11, Sample2 = clean.btw.b1.var.df$topo.var21, Category = new.cat.btw1.clean, Distance = as.numeric(clean.btw.b1.var.df$value))

#because the above also plots between oral and between mwf categories want to subset out some 
MWF.Oral.Subset= inter.oral.mwf.data.bray[c(which(inter.oral.mwf.data.bray$Category == "MWF.P to Oral.Machine.Shop"), which(inter.oral.mwf.data.bray $Category == "MWF.NP to Oral.Machine.Shop"), which(inter.oral.mwf.data.bray $Category == "MWF.NP to Oral.Assembly"), which(inter.oral.mwf.data.bray$Category == "MWF.P to Oral.Assembly"), which(inter.oral.mwf.data.bray $Category == "Oral.Administration to MWF.NP"),  which(inter.oral.mwf.data.bray $Category == "Oral.Administration to MWF.P")),]
unique(MWF.Oral.Subset $Category)

MWF.Oral.Subset.c <-MWF.Oral.Subset

#oral admin label is backwards so going to correct that now 
MWF.Oral.Subset.c $Category <- as.character(MWF.Oral.Subset.c$Category)
MWF.Oral.Subset.c $Category[MWF.Oral.Subset.c $Category == "Oral.Administration to MWF.P"] <- "MWF.P to Oral.Administration"
MWF.Oral.Subset.c $Category[MWF.Oral.Subset.c $Category == "Oral.Administration to MWF.NP"] <- "MWF.NP to Oral.Administration"

#export for prism
write.csv(MWF.Oral.Subset.c, file = "Inter_Bray_Oral_MWF.txt")
}

### MWF v Nasal
{
#create otu table of just the samples we want 
Nasal.MWF.otu.rel.table <- subset_samples(Genus.rel.table, Sample_Type_d1_code == "25" |Sample_Type_d1_code == "26" | Sample_Type_c_code ==20)
#calculate bray distance 
Nasal.MWF.Bray.dist = distance(Nasal.MWF.otu.rel.table, method="bray")

#melt the distance matrix into columns with each sample site and distance
b1 <- melt(as.matrix(Nasal.MWF.Bray.dist))
head(b1)

#Then need to remove self distances and duplicated distances
p1    <- t(apply(b1[,c(1,2)],1,FUN=sort))
rmv11 <- which(p1[,1] == p1[,2])

p1    <- paste(p1[,1],p1[,2],sep="|")
rmv21 <- which(duplicated(p1))

#establish new data frame that removes those values 
b1.df   <- b1[-c(rmv11,rmv21),] 
head(b1.df)

##Now we need to replace variable rows with topo 
#set up new data frame
new1.df <- b1.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new1.df[] <- lapply(b1.df, function(x) sample_data(Nasal.MWF.otu.rel.table)$Sample_Type_d1[match(x, rownames(sample_data(Nasal.MWF.otu.rel.table)))])

head(new1.df)

#create two lists of the group variable 
topo.var11 <- new1.df[,1]
topo.var21 <-new1.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b1.var.df <- cbind(b1.df, topo.var11, topo.var21)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b1.var.df <- b1.var.df
#set row names to re-zero
rownames(btw.b1.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b1.var.df)) {
	if (btw.b1.var.df$topo.var11[i] == btw.b1.var.df$topo.var21[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b1.var.df <- btw.b1.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b1.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw1 = paste(btw.b1.var.df$topo.var11, "to", btw.b1.var.df$topo.var21)
new.cat.btw1.df <- data.frame(btw.b1.var.df$topo.var11, btw.b1.var.df$topo.var21)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw1.df, 1, sort))
unique.new.cat.btw1 <- unique(new.cat.btw1.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw1) <- NULL
rownames(unique.new.cat.btw1) <- NULL
unique.new.cat.btw1 <- paste(unique.new.cat.btw1[,1], "to", unique.new.cat.btw1[,2])


#create new data frame 
clean.btw.b1.var.df <- btw.b1.var.df

#reset row names
rownames(clean.btw.b1.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b1.var.df)){
	if (paste(clean.btw.b1.var.df$topo.var21[i], "to", clean.btw.b1.var.df$topo.var11[i]) %in% unique.new.cat.btw1) {
		clean.btw.b1.var.df$topo.var11[i] <- btw.b1.var.df$topo.var21[i]
		clean.btw.b1.var.df$topo.var21[i] <- btw.b1.var.df$topo.var11[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw1.clean = paste(clean.btw.b1.var.df$topo.var11, "to", clean.btw.b1.var.df$topo.var21)

#confirm permutations 
unique(new.cat.btw1.clean)

# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
inter.nasal.mwf.data.bray <- data.frame(Sample1 = clean.btw.b1.var.df$topo.var11, Sample2 = clean.btw.b1.var.df$topo.var21, Category = new.cat.btw1.clean, Distance = as.numeric(clean.btw.b1.var.df$value))

#because the above also plots between skin and between mwf categories want to subset out some 
MWF.Nasal.Subset= inter.nasal.mwf.data.bray[c(which(inter.nasal.mwf.data.bray $Category == "MWF.P to Nasal.Machine.Shop"), which(inter.nasal.mwf.data.bray $Category == "MWF.NP to Nasal.Machine.Shop"), which(inter.nasal.mwf.data.bray $Category == "MWF.NP to Nasal.Assembly"), which(inter.nasal.mwf.data.bray $Category == "MWF.P to Nasal.Assembly"), which(inter.nasal.mwf.data.bray $Category == "MWF.P to Nasal.Administration"),  which(inter.nasal.mwf.data.bray $Category == "MWF.NP to Nasal.Administration")),]

write.csv(MWF.Nasal.Subset, file = "Inter_Bray_Nasal_MWF.txt")
}

### MWF v Lung 
{
#select just the samples you want for an otu table	
lung.mwf.otu.rel.table = subset_samples(Genus.rel.table, Sample_Type_d_fluid_change_code == "101" | Sample_Type_d_fluid_change_code == "102" | Sample_Type_d1_code == "25" | Sample_Type_d1_code == "26")
#since there is no variable that encompasses both lung and mwf variable labels, we will create one called newvar and use that 
sample_data(lung.mwf.otu.rel.table)$Description
var1 = as.character(get_variable(lung.mwf.otu.rel.table,"Description"))
var2 = as.character(get_variable(lung.mwf.otu.rel.table,"Sample_Type_d1_c"))
sample_data(lung.mwf.otu.rel.table)$NewVar <- mapply(paste0, var1, var2, collapse="_")
merge_samples(lung.mwf.otu.rel.table, "NewVar")


## # Calculate Bray distance
lung.mwf.bray.dist = distance(lung.mwf.otu.rel.table, method="bray")
lung.mwf.bray.pco = dudi.pco(cailliez(lung.mwf.bray.dist))

#melt the distance matrix into columns with each sample site and distance
b2 <- melt(as.matrix(lung.mwf.bray.dist))
head(b2)

#Then need to remove self distances and duplicated distances
p2    <- t(apply(b2[,c(1,2)],1,FUN=sort))
rmv2 <- which(p2[,1] == p2[,2])

p2    <- paste(p2[,1],p2[,2],sep="|")
rmv3 <- which(duplicated(p2))

#establish new data frame that removes those values 
b2.df   <- b2[-c(rmv2,rmv3),] 
head(b2.df)

##Now we need to replace variable rows with loca 
#set up new data frame
new2.df <- b2.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new2.df[] <- lapply(b2.df, function(x) sample_data(lung.mwf.otu.rel.table)$NewVar[match(x, rownames(sample_data(lung.mwf.otu.rel.table)))])

head(new2.df)

#create two lists of the group variable 
topo.var11 <- new2.df[,1]
topo.var21 <-new2.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b2.var.df <- cbind(b2.df, topo.var11, topo.var21)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b2.var.df <- b2.var.df
#set row names to re-zero
rownames(btw.b2.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b2.var.df)) {
	if (btw.b2.var.df$topo.var11[i] == btw.b2.var.df$topo.var21[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b2.var.df <- btw.b2.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b2.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw2 = paste(btw.b2.var.df$topo.var11, "to", btw.b2.var.df$topo.var21)
new.cat.btw2.df <- data.frame(btw.b2.var.df$topo.var11, btw.b2.var.df$topo.var21)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw2.df, 1, sort))
unique.new.cat.btw2 <- unique(new.cat.btw2.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw2) <- NULL
rownames(unique.new.cat.btw2) <- NULL
unique.new.cat.btw2 <- paste(unique.new.cat.btw2[,1], "to", unique.new.cat.btw2[,2])


#create new data frame 
clean.btw.b2.var.df <- btw.b2.var.df

#reset row names
rownames(clean.btw.b2.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b2.var.df)){
	if (paste(clean.btw.b2.var.df$topo.var21[i], "to", clean.btw.b2.var.df$topo.var11[i]) %in% unique.new.cat.btw2) {
		clean.btw.b2.var.df$topo.var11[i] <- btw.b2.var.df$topo.var21[i]
		clean.btw.b2.var.df$topo.var21[i] <- btw.b2.var.df$topo.var11[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw2.clean = paste(clean.btw.b2.var.df$topo.var11, "to", clean.btw.b2.var.df$topo.var21)

#confirm permutations 
unique(new.cat.btw2.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
inter.lung.mwf.data.bray <- data.frame(Sample1 = clean.btw.b2.var.df$X1, Sample2 = clean.btw.b2.var.df$X2, Category = new.cat.btw2.clean, Distance = as.numeric(clean.btw.b2.var.df$value))

#because the above also plots between tissue and between mwf categories want to subset out some, select for the ones you want to keep 
MWF.Lung.Subset= inter.lung.mwf.data.bray[c(which(inter.lung.mwf.data.bray$Category == "MWFMWF.P to Tissue_control"), which(inter.lung.mwf.data.bray$Category == "MWFMWF.P to Tissue_case"), which(inter.lung.mwf.data.bray$Category == "MWFMWF.NP to Tissue_control"), which(inter.lung.mwf.data.bray$Category == "MWFMWF.NP to Tissue_case")),]

#save as txt file 
write.csv(MWF.Lung.Subset, file = "Inter_Bray_Lung_MWF.txt")	
#Then take this file to prism for creating the image desired 	
}	
}	


###Figure 2A Volcano plot and DeSeq analysis as included in the supplement 
{
#Choose Alpha/FDR (here you can go more or less conservative, choose one)
alpha = 0.5

##Lung 
{
## Initial DeSeq processing and generation of data tables 
{
# Subset OTU table, nonrelative 
Lung.otu.table = subset_samples(otu.table, Sample_Type_d_fluid_change_code == 101 | Sample_Type_d_fluid_change_code == 102)

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(Lung.otu.table, ~Sample_Type_d_fluid_change)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)


#Drop rows with no data in your comparison variable
diagdds$Sample_Type_d_fluid_change <- droplevels(diagdds$Sample_Type_d_fluid_change)

#Choose which is the 'control' in your comparison variable, in this case, within the category Sample_Type_d_fluid_change 
diagdds$Sample_Type_d_fluid_change <- relevel(diagdds$Sample_Type_d_fluid_change, ref ="Control")

# Run DESEQ
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Get Results of DESeq into a table
res = results(diagdds, cooksCutoff = FALSE)


# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]
res.lung.otu = res

##Generate tables 

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(Lung.otu.table)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table as Taxa Names
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res$OTU <- res$Gene.symbol
res$Gene.symbol <- res$row2

otu.to.save <-as.character(res$OTU)

#get abundnace data - RAW
abundance.lung <- as.data.frame(taxa_sums(Lung.otu.table))
abundance.lung$OTU <-rownames(abundance.lung)
abundance.lung <-abundance.lung[abundance.lung$OTU %in% otu.to.save,]
abundance.lung <- as.data.frame(abundance.lung)


#get abundance data - mean relative - use otu.relative.table
#subset lung
Lung.relative.otu.table = subset_samples(otu.relative.table, Sample_Type_d_fluid_change_code == 101 | Sample_Type_d_fluid_change_code == 102)

#from relative table we should get the mean across the row of the otu table
lung.otu.rel.df <- data.frame(otu_table(Lung.relative.otu.table))
lung.meanRA <- rowMeans(lung.otu.rel.df)

#need to subset AND reorder just the otus that we have 
lung.meanRA.save <- lung.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- lung.meanRA.save

#Keep only the variables you need for pathway analysis
res.PA <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.IPA <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
#For Pathway Analysis
write.table(res.PA,file=  "Lung_OTU_DeSeq.txt", sep="\t", col.names = NA, row.names = TRUE)
#For IPA
write.table(res.IPA,file="Lung_OTU_DeSeq_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
}


##DeSeq plot in supplement 
{
alpha = 0.2

#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])

#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Lung.Genus.table))), ntaxa(Lung.Genus.table)), Lung.Genus.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))


# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comparing Status to a character vector from a factor
dat$Location <- as.character(dat$Sample_Type_d_fluid_change)

# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$Genus <- dat$o

#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "Tissue.Case", "Tissue.Control"))

#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]


#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c('#FF0000','#0000FF')) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c('#FF0000','#0000FF')) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("OTU") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())


pdf("DESeq2_Lung_NIOSH.pdf", height = 2, width = 10)
print(p)
dev.off()

}


##Violin Plot 
{
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")


#Lung
res.lung = res
res.lung$sample = "Lung"
res.lung$otu = row.names(res.lung)
res.lung = cbind(as(res.lung, "data.frame"), as(tax_table(Lung.Genus.table)[rownames(res.lung), ], "matrix"))
res.lung$row2 <- paste(res.lung$Domain,res.lung$Phylum,res.lung$Class,res.lung$Order,res.lung$Family,res.lung$Genus)
#Replace Spaces with .
res.lung$row2 <- gsub('\\s+', '.', res.lung$row2)

res <- res.lung

#Fix Labeling of Taxa Trail
res$Genus[res$Genus =="g__"] <- NA
res$o <-ifelse(is.na(res$Genus), paste(as.character(res$Family),"(u.g.)",sep=" "), as.character(res$Genus))
res$Family[res$Family =="f__"] <- NA
res$o <-ifelse(is.na(res$Family), paste(as.character(res$Order),"(u.g.)",sep=" "), as.character(res$o))
res$Order[res$Order =="o__"] <- NA
res$o <-ifelse(is.na(res$Order), paste(as.character(res$Class),"(u.g.)",sep=" "), as.character(res$o))
res$Class[res$Class =="c__"] <- NA
res$o <-ifelse(is.na(res$Class), paste(as.character(res$Phylum),"(u.g.)",sep=" "), as.character(res$o))
res$Phylum[res$Phylum =="p__"] <- NA
res$o <-ifelse(is.na(res$Phylum), paste(as.character(res$Domain),"(u.g.)",sep=" "), as.character(res$o))


#label all top genuses, can see what else is significant
    ggplot(res, aes(x= sample, y=sig, label=otu, fill=sample)) + 
    geom_violin()+
    geom_jitter(shape=1, position=position_jitter(0.2))+
 geom_text_repel(aes(label=ifelse(res$sig>1.3 , as.character(res$o),'')),size=3,force=25) +
     #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme

## Final figure - label just pseudomonas 813945 IF it is signifncant and add dot size by scaled mean relative abundance 
pdf("DESeq2_VPlot_LungOTU.Abudance.pdf", height = 8, width = 5)
    ggplot(res, aes(x= sample, y=sig, label=otu,)) + 
    geom_violin(fill = "gold3")+
    geom_point(position=position_jitter(0.2), size = 50* res$abundance, color = "black")+
      geom_text_repel(aes(label=ifelse(res$o=="f__Pseudomonadaceae (u.g.)" & res$sig>1.3, as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
 	#add horizontal line at p=0.05, dashed and red in color 
 	geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    theme
dev.off()

}
}

##Skin 
{
## Initial DeSeq processing and generation of data tables 
{
#OTU 
Skin.otu.table = subset_samples(otu.table,   Sample_Type_c_code ==22)

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(Skin.otu.table, ~ Sample_Type_d_fluid_change_c)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)

# Run DESEQ
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Get Results of DESeq into a table
res = results(diagdds, cooksCutoff = FALSE)
res.skin = res
# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(Skin.otu.table)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table as Taxa Names
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res$OTU <- res$Gene.symbol
res$Gene.symbol <- res$row2

otu.to.save <-as.character(res$OTU)



#get abundance data - mean relative - use otu.relative.table
#subset lung
Skin.relative.otu.table = subset_samples(otu.relative.table, Sample_Type_c_code ==22)

#from relative table we should get the mean across the row of the otu table
skin.otu.rel.df <- data.frame(otu_table(Skin.relative.otu.table))
skin.meanRA <- rowMeans(skin.otu.rel.df)

#need to subset AND reorder just the otus that we have 
skin.meanRA.save <- skin.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- skin.meanRA.save

#Keep only the variables you need for pathway analysis
res.PA <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.IPA <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
#For Pathway Analysis
write.table(res.PA,file=  "Skin_OTU_DeSeq.txt", sep="\t", col.names = NA, row.names = TRUE)
#For IPA
write.table(res.IPA,file="Skin_OTU_DeSeq_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
	
}

#DeSeq plots in supp
{
#Set FDR Threshold
alpha = 0.2

#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])
str9
#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Skin.Genus.table))), ntaxa(Skin.Genus.table)), Skin.Genus.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comparing Status to a character vector from a factor
dat$Location <- as.character(dat$Sample_Type_d_fluid_change_c)

# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$Genus <- dat$o

#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "Skin.Machine.Shop", "Skin.Assembly", "Skin.Administration"))

#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]


#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c('#98FB98','#FFA500', "#800080")) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c('#98FB98','#FFA500', "#800080")) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("OTU") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())


pdf("DESeq2_Skin_NIOSH.pdf", height = 10, width = 10)
print(p)
dev.off()


}


##Violin Plot
{
	theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")



#Skin
res.skin <- res
res.skin$sample = "Skin"
res.skin$otu = row.names(res.skin)
res.skin = cbind(as(res.skin, "data.frame"), as(tax_table(Skin.Genus.table)[rownames(res.lung), ], "matrix"))
res.skin$row2 <- paste(res.skin$Domain,res.skin$Phylum,res.skin$Class,res.skin$Order,res.skin$Family,res.skin$Genus)
#Replace Spaces with .
res.skin$row2 <- gsub('\\s+', '.', res.skin$row2)
res<-res.skin


#Fix Labeling of Taxa Trail
res$Genus[res$Genus =="g__"] <- NA
res$o <-ifelse(is.na(res$Genus), paste(as.character(res$Family),"(u.g.)",sep=" "), as.character(res$Genus))
res$Family[res$Family =="f__"] <- NA
res$o <-ifelse(is.na(res$Family), paste(as.character(res$Order),"(u.g.)",sep=" "), as.character(res$o))
res$Order[res$Order =="o__"] <- NA
res$o <-ifelse(is.na(res$Order), paste(as.character(res$Class),"(u.g.)",sep=" "), as.character(res$o))
res$Class[res$Class =="c__"] <- NA
res$o <-ifelse(is.na(res$Class), paste(as.character(res$Phylum),"(u.g.)",sep=" "), as.character(res$o))
res$Phylum[res$Phylum =="p__"] <- NA
res$o <-ifelse(is.na(res$Phylum), paste(as.character(res$Domain),"(u.g.)",sep=" "), as.character(res$o))


#label all top
    ggplot(res, aes(x= sample, y=sig, label=otu, fill=sample)) + 
    geom_violin()+
    geom_jitter(shape=1, position=position_jitter(0.2))+
 geom_text_repel(aes(label=ifelse(res$sig>1.3 , as.character(res$o),'')),size=3,force=25) +
     #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme

## size by abundance 
pdf("DESeq2_VPlot_SkinOTU.Abudance.pdf", height = 8, width = 5)
    ggplot(res, aes(x= sample, y=sig, label=otu)) + 
    geom_violin(fill="deeppink1")+
    geom_point(position=position_jitter(0.2), size = 50* res$abundance, color = "black")+
      geom_text_repel(aes(label=ifelse(res$o=="f__Pseudomonadaceae (u.g.)" & res$sig>1.3, as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
 	#add horizontal line at p=0.05, dashed and red in color 
 	geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    theme
dev.off()

} 


}

##Nasal 
{
## Initial DeSeq processing and generation of data tables 
{
#OTU 
Nasal.otu.table = subset_samples(otu.table,  Sample_Type_c_code ==20)

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(Nasal.otu.table, ~ Sample_Type_d_fluid_change_c)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)

# Run DESEQ
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Get Results of DESeq into a table
res = results(diagdds, cooksCutoff = FALSE)
res.nasal = res
# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(Nasal.otu.table)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table as Taxa Names
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res$OTU <- res$Gene.symbol
res$Gene.symbol <- res$row2

otu.to.save <-as.character(res$OTU)



#get abundance data - mean relative - use otu.relative.table
#subset lung
Nasal.relative.otu.table = subset_samples(otu.relative.table, Sample_Type_c_code ==20)

#from relative table we should get the mean across the row of the otu table
nasal.otu.rel.df <- data.frame(otu_table(Nasal.relative.otu.table))
nasal.meanRA <- rowMeans(nasal.otu.rel.df)

#need to subset AND reorder just the otus that we have 
nasal.meanRA.save <- nasal.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- nasal.meanRA.save

#Keep only the variables you need for pathway analysis
res.PA <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.IPA <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
#For Pathway Analysis
write.table(res.PA,file=  "Nasal_OTU_DeSeq.txt", sep="\t", col.names = NA, row.names = TRUE)
#For IPA
write.table(res.IPA,file="Nasal_OTU_DeSeq_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
}

##DeSeq barplot - in supplmenet 
{

#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])

#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Nasal.Genus.table))), ntaxa(Nasal.Genus.table)), Nasal.Genus.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comparing Status to a character vector from a factor
dat$Location <- as.character(dat$Sample_Type_d_fluid_change_c)

# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$Genus <- dat$o

#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "Nasal.Machine.shop", "Nasal.Assembly", "Nasal.Administration"))


#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]


#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c('#98FB98','#FFA500', "#800080")) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c('#98FB98','#FFA500', "#800080")) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("OTU") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())


pdf("DESeq2_Nasal_NIOSH.pdf", height = 10, width = 10)
print(p)
dev.off()	
	
}


##Violin Plot
{
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")



#Nasal
res.nasal <- res
res.nasal$sample = "Nasal"
res.nasal$otu = row.names(res.nasal)
res.nasal = cbind(as(res.nasal, "data.frame"), as(tax_table(Nasal.Genus.table)[rownames(res.nasal), ], "matrix"))
res.nasal$row2 <- paste(res.nasal$Domain,res.nasal$Phylum,res.nasal$Class,res.nasal$Order,res.nasal$Family,res.nasal$Genus)
#Replace Spaces with .
res.nasal$row2 <- gsub('\\s+', '.', res.nasal$row2)
res<-res.nasal


#Fix Labeling of Taxa Trail
res$Genus[res$Genus =="g__"] <- NA
res$o <-ifelse(is.na(res$Genus), paste(as.character(res$Family),"(u.g.)",sep=" "), as.character(res$Genus))
res$Family[res$Family =="f__"] <- NA
res$o <-ifelse(is.na(res$Family), paste(as.character(res$Order),"(u.g.)",sep=" "), as.character(res$o))
res$Order[res$Order =="o__"] <- NA
res$o <-ifelse(is.na(res$Order), paste(as.character(res$Class),"(u.g.)",sep=" "), as.character(res$o))
res$Class[res$Class =="c__"] <- NA
res$o <-ifelse(is.na(res$Class), paste(as.character(res$Phylum),"(u.g.)",sep=" "), as.character(res$o))
res$Phylum[res$Phylum =="p__"] <- NA
res$o <-ifelse(is.na(res$Phylum), paste(as.character(res$Domain),"(u.g.)",sep=" "), as.character(res$o))


#label all top
    ggplot(res, aes(x= sample, y=sig, label=otu, fill=sample)) + 
    geom_violin()+
    geom_jitter(shape=1, position=position_jitter(0.2))+
 geom_text_repel(aes(label=ifelse(res$sig>1.3 , as.character(res$o),'')),size=3,force=25) +
     #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme

## size by abundance  
pdf("DESeq2_VPlot_NasalOTU.Abundance.pdf", height = 8, width = 5)
    ggplot(res, aes(x= sample, y=sig, label=otu,)) + 
    geom_violin(fill="cyan2")+
    geom_point(position=position_jitter(0.2), size = 50* res$abundance, color = "black")+
      geom_text_repel(aes(label=ifelse(res$o=="f__Pseudomonadaceae (u.g.)" & res$sig>1.3, as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
 	#add horizontal line at p=0.05, dashed and red in color 
 	geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    theme
dev.off()
} 


}

##Oral 
{
## Initial DeSeq processing and generation of data tables 
{
#OTU 
Oral.otu.table = subset_samples(otu.table,  Sample_Type_c_code ==21)

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(Oral.otu.table, ~ Sample_Type_d_fluid_change_c)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)

# Run DESEQ
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Get Results of DESeq into a table
res = results(diagdds, cooksCutoff = FALSE)
res.oral = res
# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(Oral.otu.table)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table as Taxa Names
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res$OTU <- res$Gene.symbol
res$Gene.symbol <- res$row2

otu.to.save <-as.character(res$OTU)



#get abundance data - mean relative - use otu.relative.table
#subset lung
Oral.relative.otu.table = subset_samples(otu.relative.table, Sample_Type_c_code ==21)

#from relative table we should get the mean across the row of the otu table
oral.otu.rel.df <- data.frame(otu_table(Oral.relative.otu.table))
oral.meanRA <- rowMeans(oral.otu.rel.df)

#need to subset AND reorder just the otus that we have 
oral.meanRA.save <- oral.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- oral.meanRA.save

#Keep only the variables you need for pathway analysis
res.PA <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.IPA <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
#For Pathway Analysis
write.table(res.PA,file=  "Oral_OTU_DeSeq.txt", sep="\t", col.names = NA, row.names = TRUE)
#For IPA
write.table(res.IPA,file="Oral_OTU_DeSeq_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
}


##DeSeq bar plot in Supplmenet 
{
#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])

#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Oral.Genus.table))), ntaxa(Oral.Genus.table)), Oral.Genus.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# agglomerate taxa
#glom <- tax_glom(phy, taxrank = 'Genus')

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comparing Status to a character vector from a factor
dat$Location <- as.character(dat$Sample_Type_d_fluid_change_c)

# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$Genus <- dat$o

#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "Oral.Machine.Shop", "Oral.Assembly", "Oral.Administration"))

#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]

#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c('#98FB98','#FFA500', "#800080")) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c('#98FB98','#FFA500', "#800080")) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("OTU") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())


pdf("DESEq2_Oral_NIOSH.pdf", height = 10, width = 10)
print(p)
dev.off()	
}

##Violin Plot 
{
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")



#Oral
res.oral <- res
res.oral$sample = "Oral"
res.oral$otu = row.names(res.oral)
res.oral = cbind(as(res.oral, "data.frame"), as(tax_table(Oral.Genus.table)[rownames(res.oral), ], "matrix"))
res.oral$row2 <- paste(res.oral$Domain,res.oral$Phylum,res.oral$Class,res.oral$Order,res.oral$Family,res.oral$Genus)
#Replace Spaces with .
res.oral$row2 <- gsub('\\s+', '.', res.oral$row2)
res<-res.oral


#Fix Labeling of Taxa Trail
res$Genus[res$Genus =="g__"] <- NA
res$o <-ifelse(is.na(res$Genus), paste(as.character(res$Family),"(u.g.)",sep=" "), as.character(res$Genus))
res$Family[res$Family =="f__"] <- NA
res$o <-ifelse(is.na(res$Family), paste(as.character(res$Order),"(u.g.)",sep=" "), as.character(res$o))
res$Order[res$Order =="o__"] <- NA
res$o <-ifelse(is.na(res$Order), paste(as.character(res$Class),"(u.g.)",sep=" "), as.character(res$o))
res$Class[res$Class =="c__"] <- NA
res$o <-ifelse(is.na(res$Class), paste(as.character(res$Phylum),"(u.g.)",sep=" "), as.character(res$o))
res$Phylum[res$Phylum =="p__"] <- NA
res$o <-ifelse(is.na(res$Phylum), paste(as.character(res$Domain),"(u.g.)",sep=" "), as.character(res$o))


#label all top
    ggplot(res, aes(x= sample, y=sig, label=otu, fill=sample)) + 
    geom_violin()+
    geom_jitter(shape=1, position=position_jitter(0.2))+
 geom_text_repel(aes(label=ifelse(res$sig>1.3 , as.character(res$o),'')),size=3,force=25) +
     #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme

## size by abundance 
pdf("DESeq2_VPlot_OralOTU.Abundance.pdf", height = 8, width = 5)
    ggplot(res, aes(x= sample, y=sig, label=otu,)) + 
    geom_violin(fill="dodgerblue2")+
    geom_point(position=position_jitter(0.2), size = 50* res$abundance, color = "black")+
      geom_text_repel(aes(label=ifelse(res$o=="f__Pseudomonadaceae (u.g.)" & res$sig>1.3, as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
 	#add horizontal line at p=0.05, dashed and red in color 
 	geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    theme
dev.off()	
}

}

##MWF 
{
## Initial DeSeq processing and generation of data tables 
{
# subset otu table 
MWF.otu.table = subset_samples(otu.table,  Sample_Type_d1_code %in% c("25","26"))

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(MWF.Genus.table, ~Sample_Type_d1)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)

# Run DESEQ
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Get Results of DESeq into a table
res = results(diagdds, cooksCutoff = FALSE)
res.mwf = res
# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]


#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(MWF.otu.table)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table as Taxa Names
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res$OTU <- res$Gene.symbol
res$Gene.symbol <- res$row2

#save the list of OTUs in the res dataframe for selection later
otu.to.save <-as.character(res$OTU)

#get abundance data - mean relative - use otu.relative.table
#subset lung
MWF.relative.otu.table = subset_samples(otu.relative.table, Sample_Type_d1_code %in% c("25","26"))

#from relative table we should get the mean across the row of the otu table
mwf.otu.rel.df <- data.frame(otu_table(MWF.relative.otu.table))
mwf.meanRA <- rowMeans(mwf.otu.rel.df)

#need to subset AND reorder just the otus that we have 
mwf.meanRA.save <- mwf.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- mwf.meanRA.save

#Keep only the variables you need for pathway analysis
res.PA <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.IPA <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
#For Pathway Analysis
write.table(res.PA,file=  "MWF_OTU_DeSeq.txt", sep="\t", col.names = NA, row.names = TRUE)
#For IPA
write.table(res.IPA,file="MWF_OTU_DeSeq_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
}

##DeSeq bar plot - in supplement 
{
#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])

#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(MWF.Genus.table))), ntaxa(MWF.Genus.table)), MWF.Genus.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comparing Status to a character vector from a factor
dat$Location <- as.character(dat$Sample_Type_d1)

# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$Genus <- dat$o

#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "MWF.NP", "MWF.P"))

#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]


#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c('#F08080','#008B8B')) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c('#F08080','#008B8B')) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("OTU") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())


pdf("DESeq2_MWF_NIOSH.pdf", height = 10, width = 10)
print(p)
dev.off()	
}

##Violin Plot 
{
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")



#MWF
res.mwf <- res
res.mwf$sample = "MWF"
res.mwf$otu = row.names(res.mwf)
res.mwf = cbind(as(res.mwf, "data.frame"), as(tax_table(MWF.Genus.table)[rownames(res.mwf), ], "matrix"))
res.mwf$row2 <- paste(res.mwf$Domain,res.mwf$Phylum,res.mwf$Class,res.mwf$Order,res.mwf$Family,res.mwf$Genus)
#Replace Spaces with .
res.mwf$row2 <- gsub('\\s+', '.', res.mwf$row2)
res<-res.mwf


#Fix Labeling of Taxa Trail
res$Genus[res$Genus =="g__"] <- NA
res$o <-ifelse(is.na(res$Genus), paste(as.character(res$Family),"(u.g.)",sep=" "), as.character(res$Genus))
res$Family[res$Family =="f__"] <- NA
res$o <-ifelse(is.na(res$Family), paste(as.character(res$Order),"(u.g.)",sep=" "), as.character(res$o))
res$Order[res$Order =="o__"] <- NA
res$o <-ifelse(is.na(res$Order), paste(as.character(res$Class),"(u.g.)",sep=" "), as.character(res$o))
res$Class[res$Class =="c__"] <- NA
res$o <-ifelse(is.na(res$Class), paste(as.character(res$Phylum),"(u.g.)",sep=" "), as.character(res$o))
res$Phylum[res$Phylum =="p__"] <- NA
res$o <-ifelse(is.na(res$Phylum), paste(as.character(res$Domain),"(u.g.)",sep=" "), as.character(res$o))


#label all top genuses for visualisation 
    ggplot(res, aes(x= sample, y=sig, label=otu, fill=sample)) + 
    geom_violin()+
    geom_jitter(shape=1, position=position_jitter(0.2))+
 	geom_text_repel(aes(label=ifelse(res$sig>1.3 , as.character(res$o),'')),size=3,force=25) +
    ylab("-log10(adjusted p-value)") + 
    theme

##size by abundance 
pdf("DESeq2_VPlot_MWFOTU.Abudance.pdf", height = 8, width = 5)
    ggplot(res, aes(x= sample, y=sig, label=otu,)) + 
    geom_violin(fill="green4")+
    geom_point(position=position_jitter(0.2), size = 25* res$abundance, color = "black")+
    geom_text_repel(aes(label=ifelse(res$o=="f__Pseudomonadaceae (u.g.)" & res$sig>1.3, as.character(res$o),'')),size=3,force=25) +
    ylab("-log10(adjusted p-value)") + 
 	geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    theme
dev.off()

}

}

##Air 
{
## Initial DeSeq Processing and generation of data datables 
{
#subset from nonrelative OTU 
Air.otu.table = subset_samples(otu.table,  Sample_Type_d1_code == 41 | Sample_Type_d1_code == 42 | Sample_Type_d1_code == 43)

#Convert Phyloseq Object to DESEq, correncting for any potential confounders
diagdds <- phyloseq_to_deseq2(Air.otu.table, ~ Sample_Type_d1)

#Calculate geometric means prior to estimate size factor
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)

# Estimate Size, Disperssion and Variance
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)

# Run DESEQ
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Get Results of DESeq into a table
res = results(diagdds, cooksCutoff = FALSE)

# Reorder Results based on FDR
res = res[order(res$padj, na.last = NA), ]
res.air <-res

#Get Taxa Names from Phyloseq Object
res = cbind(as(res, "data.frame"), as(tax_table(Air.otu.table)[rownames(res), ], "matrix"))

#Replace OTU with Taxa
res$row2 <- paste(res$Domain,res$Phylum,res$Class,res$Order,res$Family,res$Genus)
#Replace Spaces with .
res$row2 <- gsub('\\s+', '.', res$row2)

#Convert Resuts table into a data.frame
res <- as.data.frame(res)

#Set Names of Results Table as Taxa Names
res <- setNames(cbind(rownames(res), res, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","Domain","Phylum","Class","Order","Family","Genus","OTU","row2"))
res$OTU <- res$Gene.symbol
res$Gene.symbol <- res$row2

otu.to.save <-as.character(res$OTU)



#get abundance data - mean relative - use otu.relative.table
#subset lung
Air.relative.otu.table = subset_samples(otu.relative.table, Sample_Type_d1_code == 41 | Sample_Type_d1_code == 42 | Sample_Type_d1_code == 43)

#from relative table we should get the mean across the row of the otu table
air.otu.rel.df <- data.frame(otu_table(Air.relative.otu.table))
air.meanRA <- rowMeans(air.otu.rel.df)

#need to subset AND reorder just the otus that we have 
air.meanRA.save <- air.meanRA[otu.to.save]

#add the abundnace data for the res dataframe
res$abundance <- air.meanRA.save

#Keep only the variables you need for pathway analysis
res.PA <- res[,c("Gene.symbol","logFC","adj.P.Val")]
res.IPA <- res[,c("Gene.symbol","logFC","pvalue","adj.P.Val")]

#Write Tables to TXT file
#For Pathway Analysis
write.table(res.PA,file=  "Air_OTU_DeSeq.txt", sep="\t", col.names = NA, row.names = TRUE)
#For IPA
write.table(res.IPA,file="Air_OTU_DeSeq_IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
}

## DeSeq plot in the supplement 
{
#Set FDR Threshold
alpha = 0.5

#Select top 50 Significant OTUS
keepOTUs = rownames(res[res$padj < alpha, ])

#Create New Labels for Taxa
    x10 = prune_taxa(tail(names(sort(taxa_sums(Air.Genus.table))), ntaxa(Air.Genus.table)), Air.Genus.table)
    tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
    myranks = c( "Family", "Genus", "Strain")
    mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
    # Add concatenated labels as a new rank after strain
    tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)

#Keep OTUs that are significant in DESeq
phy =  prune_taxa(keepOTUs, x10)

# get abundance in % (relative abundance)
phy <- transform_sample_counts(phy, function(x) x/sum(x))

# create dataframe from phyloseq object
dat <- data.table(psmelt(phy))

# convert Genus to a character vector from a factor
dat$Genus <- as.character(dat$catglab)

# convert Comparing Status to a character vector from a factor
dat$Location <- as.character(dat$Sample_Type_d1)

# Fix taxa name for genus or OTU with (u.g.) if applicable
dat <- separate(data=dat, col=Genus, into= c("f", "g", "o"), sep="__")
dat <- separate(data=dat, col=o, into= c("o", "otu"), sep="_")
dat <- separate(data=dat, col=g, into= c("f2", "g"), sep="_")
dat[dat==""] <- NA
dat$o[is.na(dat$o)] <- paste(as.character(dat$f2[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$o <- ifelse(is.na(dat$f2), dat$f2, dat$o)
dat$o[is.na(dat$o)] <- paste("OTU", as.character(dat$otu[is.na(dat$o)]),"(u.g.)", sep=" ")
dat$Genus <- dat$o

#Calculate Median and IQR
data <- setDT(dat)[,list(Abundance=as.numeric(median(Abundance, na.rm=TRUE)), iqr=as.numeric(quantile(Abundance, probs=.75, na.rm=TRUE))), by=c("Location", "Genus")]

data.table <-as.data.table(data)
new <- data.table[data.table[ , .I[Abundance == max(Abundance)], by = Genus]$V1]
new2 <- new[new[ , .I[iqr == max(iqr)], by = Genus]$V1]
new2$facet <- new2$Location
new2 <- mutate(new2, facet=fct_relevel(facet, "Machine.Shop", "Assembly", "Administration"))


#for each genus need to assign facet in original data frame
data$facet <- new2$facet[match(data$Genus, new2$Genus)]
dat$facet <- new2$facet[match(dat$Genus, new2$Genus)]


#Plot Median + IQR
p <- ggplot(data, aes(x=reorder(Genus,+Abundance), y=Abundance, fill=Location)) +
    geom_blank() +
    facet_grid(facet ~., scales = "free_y", space = "free") +
    geom_jitter(data=dat, aes(x=reorder(Genus, +Abundance), y=Abundance, fill=Location, color=Location), shape=16, position=position_dodge(.9)) +
    geom_rect(aes(xmin=as.numeric(reorder(Genus, +Abundance))-0.4,xmax=as.numeric(reorder(Genus, +Abundance))+.4, ymax=Abundance, ymin=0, fill=Location),color="black",position=position_dodge(0.9)) +
    scale_fill_manual(values=c('#98FB98','#FFA500', "#800080")) +
    geom_errorbar(aes(ymin=Abundance, ymax=iqr), width=.2, position=position_dodge(.9)) +
    scale_color_manual(values=c('#98FB98','#FFA500', "#800080")) +
    scale_y_continuous(name="Relative Abundance",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
    xlab("OTU") +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.y = element_blank(), strip.background = element_blank())

pdf("DESeq2_Air_NIOSH.pdf", height = 10, width = 10)
print(p)
dev.off()
}

##Violin Plot 
{
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")



#Air
res.air <- res
res.air$sample = "Air"
res.air$otu = row.names(res.air)
res.air = cbind(as(res.air, "data.frame"), as(tax_table(Air.otu.table)[rownames(res.air), ], "matrix"))
res.air$row2 <- paste(res.air$Domain,res.air$Phylum,res.air$Class,res.air$Order,res.air$Family,res.air$Genus)
#Replace Spaces with .
res.air$row2 <- gsub('\\s+', '.', res.air$row2)
res<-res.air


#Fix Labeling of Taxa Trail
res$Genus[res$Genus =="g__"] <- NA
res$o <-ifelse(is.na(res$Genus), paste(as.character(res$Family),"(u.g.)",sep=" "), as.character(res$Genus))
res$Family[res$Family =="f__"] <- NA
res$o <-ifelse(is.na(res$Family), paste(as.character(res$Order),"(u.g.)",sep=" "), as.character(res$o))
res$Order[res$Order =="o__"] <- NA
res$o <-ifelse(is.na(res$Order), paste(as.character(res$Class),"(u.g.)",sep=" "), as.character(res$o))
res$Class[res$Class =="c__"] <- NA
res$o <-ifelse(is.na(res$Class), paste(as.character(res$Phylum),"(u.g.)",sep=" "), as.character(res$o))
res$Phylum[res$Phylum =="p__"] <- NA
res$o <-ifelse(is.na(res$Phylum), paste(as.character(res$Domain),"(u.g.)",sep=" "), as.character(res$o))


#label all top
    ggplot(res, aes(x= sample, y=sig, label=otu, fill=sample)) + 
    geom_violin()+
    geom_jitter(shape=1, position=position_jitter(0.2))+
 geom_text_repel(aes(label=ifelse(res$sig>1.3 , as.character(res$o),'')),size=3,force=25) +
     #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
    #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
    #geom_point(color=cols) +
    #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
    theme

## size by abundance 
pdf("DESeq2_VPlot_AirOTU.Abudance.pdf", height = 8, width = 5)
    ggplot(res, aes(x= sample, y=sig, label=otu,)) + 
    geom_violin(fill="tomato1")+
    geom_point(position=position_jitter(0.2), size = 25* res$abundance, color = "black")+
      geom_text_repel(aes(label=ifelse(res$o=="f__Pseudomonadaceae (u.g.)" & res$sig>1.3, as.character(res$o),'')),size=3,force=25) +
    #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
    #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
    ylab("-log10(adjusted p-value)") + 
 	#add horizontal line at p=0.05, dashed and red in color 
 	geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
    theme
dev.off()

}

}

}



###Figure 2B Relative Abundance
# This figure was developed in prism but the data files were process in R as follows
{		
#create genus table of just 813945 abundance in the lung to add to abundance data in prism already compiled from biomfile 
pseudomonas.lung.genus <- get_sample(lung.otu.rel.table, "813945")
#obtain what each sample is for labeling
code <- sample_data(lung.otu.rel.table)$Sample_Type_d_fluid_change_code
#bind the abundance with the label
lung.pseudomonas <- rbind(pseudomonas.lung.genus, code)
#write to a txt for export to prism
write(lung.pseudomonas, file="lung.pseudomonas.rel abundance.txt")
}

###Figure 3 - not processed in R 

####### SUPPPLEMENTARY FIGURES 

###Supp Fig 1
{
#1A was processed in prism but shannon div was export from R 
#1B - PCOA of MWF 

##(1A)Alpha diversity of  MWF 
{
#First choose environmental samples 
Enviromental.Genus.rel.table = subset_samples(Genus.rel.table, Sample_Type_a_code ==2)

# exclude NYU experiment
NIOSH.Enviromental.Genus.rel.table = subset_samples(Enviromental.Genus.rel.table, No_NYU_Exp ==1)


#Now choose just fluid samples 
Fluid.NIOSH.Enviromental.Genus.rel.table = subset_samples(NIOSH.Enviromental.Genus.rel.table, Sample_Type_E_code ==2)

#calculate shannon
Fluid.NIOSH.Enviromental.Shannon_diversity = diversity(otu_table(Fluid.NIOSH.Enviromental.Genus.rel.table), index = "shannon", MARGIN = 2, base = exp(1))

#write it to a txt file and take to prism for plotting for figure 
write.table(Fluid.NIOSH.Enviromental.Shannon_diversity, file ="Fluid.NIOSH.Enviromental.Shannon_diversity.A.txt" )

pdf(file="MWF.NIOSH.Shannon1.pdf", width=14, height=6)
boxplot(Fluid.NIOSH.Enviromental.Shannon_diversity ~ sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$Sample_Type_d1, col=c("purple", "plum2", "plum4", "black", "red", "pink1", "pink3", "grey", "blue"), outline=FALSE)
stripchart(Fluid.NIOSH.Enviromental.Shannon_diversity ~ sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$Sample_Type_d1, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
dev.off()
}

##(1B)PCOA of MWF 
{
#all MWF 
#Calculate Bray curtis for FLUID samples

#using sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$Sample_Type_d1
#this code creates a new variable based on how we'd like to name and identify the groups 
var8 = as.character(get_variable(Fluid.NIOSH.Enviromental.Genus.rel.table,"Sample_Type_d1"))

var8[sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$Sample_Type_d1 == "MWF.P"] <- "In.Use.MWF"
var8[sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$Sample_Type_d1 == "MWF.NP"] <- "In.Use.MWF"

#check names 
var8

#add the new variable back to the sample data of the phyloseq object 
sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$NewVar <- mapply(paste0, var8, collapse="_")
merge_samples(Fluid.NIOSH.Enviromental.Genus.rel.table, "NewVar")

#calculate bray 
Fluid.NIOSH.Enviromental.Bray.dist = distance(Fluid.NIOSH.Enviromental.Genus.rel.table, method="bray")
Fluid.NIOSH.Enviromental.Bray.pco = dudi.pco(cailliez(Fluid.NIOSH.Enviromental.Bray.dist))

#plot s.class
s.class(Fluid.NIOSH.Enviromental.Bray.pco$li, interaction(sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$NewVar), col=c("lightcoral", "plum2", "plum4", "black", "darkcyan", "pink1", "pink3", "grey", "blue"))

##GET AXIS
Fluid.NIOSH.Enviromental.ordinate.bray <- ordinate(Fluid.NIOSH.Enviromental.Genus.rel.table, method="PCoA", distance="bray")

# Scatter plot - Axis1 [37.1%] Axis2[7.6%]
plot_ordination(Fluid.NIOSH.Enviromental.Genus.rel.table, Fluid.NIOSH.Enviromental.ordinate.bray, color = "NewVar", shape = "NewVar")

#get pval 
adonis(Fluid.NIOSH.Enviromental.Bray.dist ~ NewVar, data=data.frame(sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)))
#Call:
adonis(formula = Fluid.NIOSH.Enviromental.Bray.dist ~ NewVar,      data = data.frame(sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sample_Type_d1  8    9.4697 1.18372  6.7595 0.45413  0.001 ***
Residuals      65   11.3827 0.17512         0.54587           
Total          73   20.8524                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#### Convert to ggplot
# build ggplot dataframe with points (x,y), corresponding groups (cluster), and location for dots
gg.all.mwf <- data.frame(cluster=factor(sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$NewVar), x2= Fluid.NIOSH.Enviromental.Bray.pco $li$A1, y2= Fluid.NIOSH.Enviromental.Bray.pco $li$A2, dots = factor(sample_data(Fluid.NIOSH.Enviromental.Genus.rel.table)$Sample_Type_d1))
# calculate group centroid locations
centroids.all.mwf <- aggregate(cbind(x2,y2)~cluster,data= gg.all.mwf,mean)
# merge centroid locations into ggplot dataframe
gg.all.mwf <- merge(gg.all.mwf, centroids.all.mwf,by="cluster",suffixes=c("",".centroid2"))


pdf(file="MWF.All.PCA_ggplot.pdf", height = 4, width = 5) 
ggplot(gg.all.mwf, aes(x = gg.all.mwf $x, y = gg.all.mwf $y2, label = cluster, color = cluster)) +
  #sets backgrounds
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2,color= cluster)) +

  #sets centroid points and draws lines to data points
  geom_point(data=centroids.all.mwf, aes(x=x2, y=y2, color=cluster), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
  #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  #stat_ellipse(level = 0.5, type = "t", show.legend = FALSE) + 
  scale_color_manual(values =c("maroon", "plum2", "plum4", "black", "pink1", "pink3", "grey", "blue")) + 
  #labels centroids 
  geom_label(data = centroids.all.mwf, aes (x=x2, y=y2))+ 
  #remove legend 
  theme(legend.position="none")
  
 dev.off()


#Now do the same for JUST the in use MWF 
#subset samples
MWF.NP.P.otu.relative.table = subset_samples(Genus.rel.table, Sample_Type_d1_code %in% c("25","26"))
#calculate bray 
MWF.NP.P.Bray.dist = distance(MWF.NP.P.otu.relative.table, method="bray")
MWF.NP.P.Bray.pco = dudi.pco(cailliez(MWF.NP.P.Bray.dist))

#Check the s.class plot 
pdf(file="MWF.NP.P.NIOSH.PCoA.Bray.pdf", width=4, height=4)
s.class(MWF.NP.P.Bray.pco$li, interaction(sample_data(MWF.NP.P.otu.relative.table)$Sample_Type_d1), col=c("purple", "red"))
dev.off()

#Calculate the p value from the anova 
adonis(MWF.NP.P.Bray.dist ~ Sample_Type_d1, data=data.frame(sample_data(MWF.NP.P.otu.relative.table)))
#Call:
adonis(formula = MWF.NP.P.Bray.dist ~ Sample_Type_d1, data = data.frame(sample_data(MWF.NP.P.otu.relative.table))) 
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Sample_Type_d1  1   0.48737 0.48737  6.4538 0.19292  0.003 **
Residuals      27   2.03897 0.07552         0.80708          
Total          28   2.52635                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##GET AXIS
MWF.NP.P.ordinate.bray <- ordinate(MWF.NP.P.otu.relative.table, method="PCoA", distance="bray")

#Check scatter plot --> AXIS1 [40.1%], AXIS 2 [19.1%]
pdf(file="MWF.NP.P.Ordinate.PCoA.Bray.pdf") 
plot_ordination(MWF.NP.P.otu.relative.table, MWF.NP.P.ordinate.bray, color = "Sample_Type_d1", shape = "Sample_Type_d1")
dev.off()

##Plotting
# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
gg.mwf <- data.frame(cluster=factor(sample_data(MWF.NP.P.otu.relative.table)$Sample_Type_d1), x=MWF.NP.P.Bray.pco$li$A1, y=MWF.NP.P.Bray.pco$li$A2)
# calculate group centroid locations
centroids.mwf <- aggregate(cbind(x,y)~cluster,data= gg.mwf,mean)
# merge centroid locations into ggplot dataframe
gg.mwf <- merge(gg.mwf,centroids.mwf,by="cluster",suffixes=c("",".centroid"))


pdf(file="MWF.NP.P.PCA_ggplot.pdf", height = 4, width = 5) 
ggplot(gg.mwf, aes(x = gg.mwf $x, y = gg.mwf $y, label = cluster, color = cluster)) +
  #sets backgrounds
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x,y=y,color= cluster)) +

  #sets centroid points and draws lines to data points
  geom_point(data=centroids.mwf, aes(x=x, y=y, color=cluster), size=0) +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, color=cluster))+ 
  #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  #stat_ellipse(level = 0.5, type = "t", show.legend = FALSE) + 
  scale_color_manual(values = c("MWF.NP" = "lightcoral", "MWF.P" = "darkcyan")) + 
  #labels centroids 
  geom_label(data = centroids.mwf, aes (x=x, y=y))+ 
  #remove legend 
  theme(legend.position="none")
  
 dev.off()
 




}
}	


###SuppFig 2
{
#2A was processed in prism but shannon div was exported from R 
#2B -PCOA of Air 

##(2A) Alpha div of Air 
{
##Alpha diveristy of all seq data 
length(All_Shannon_diversity_wOld)
all.otu.matrix <-otu_table(otu.relative.table)
All_Shannon_diversity_wOld = diversity(all.otu.matrix, index = "shannon", MARGIN = 2, base = exp(1))

All_Shannon_diversity_2_wOld = estimate_richness(otu.relative.table, measures="Shannon")

write.table(All_Shannon_diversity_wOld,file ="All_Shannon_diversity.A1.txt" )
write.table(All_Shannon_diversity_2_wOld,file ="All_Shannon_diversity.B1.txt")
}

##(2B) PCOA of Air 
{
#subset samples
Air.Genus.rel.table = subset_samples(Genus.rel.table, Sample_Type_c_code == 15)	

#this code creates a new variable based on how we'd like to name and identify the groups 
var8 = as.character(get_variable(Air.Genus.rel.table,"Sample_Type_d1"))

var8[sample_data(Air.Genus.rel.table)$Sample_Type_d1 == "Administration"] <- "In.Use.Air"
var8[sample_data(Air.Genus.rel.table)$Sample_Type_d1 == "Machine.Shop"] <- "In.Use.Air"
var8[sample_data(Air.Genus.rel.table)$Sample_Type_d1 == "Assembly"] <- "In.Use.Air"

#check names 
var8

#add the new variable back to the sample data of the phyloseq object 
sample_data(Air.Genus.rel.table)$NewVar <- mapply(paste0, var8, collapse="_")
merge_samples(Air.Genus.rel.table, "NewVar")


#calculate bray 
Air.Bray.dist = distance(Air.Genus.rel.table, method="bray")
Air.Bray.pco = dudi.pco(cailliez(Air.Bray.dist))

#plot s.class and adonis
s.class(Air.Bray.pco$li, interaction(sample_data(Air.Genus.rel.table)$NewVar), col=c("darkgreen", "purple", "grey", "red"))

adonis(Air.Bray.dist ~ Sample_Type_d1, data=data.frame(sample_data(Air.Genus.rel.table)))

                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Sample_Type_d1   5     1.695 0.33898  1.6873 0.05359  0.001 ***
Residuals      149    29.933 0.20090         0.94641          
Total          154    31.628                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1	
	
Air.ordinate.bray <- ordinate(Air.Genus.rel.table, method="PCoA", distance="bray")

# Scatter plot #Axis1 [18.2%], Axis2[9.1%]
plot_ordination(Air.Genus.rel.table, Air.ordinate.bray, color = "Sample_Type_d1", shape = "Sample_Type_d1")


#### Convert to ggplot
# build ggplot dataframe with points (x,y), corresponding groups (cluster), and location for dots
gg.all.air <- data.frame(cluster=factor(sample_data(Air.Genus.rel.table)$NewVar), x2= Air.Bray.pco $li$A1, y2= Air.Bray.pco $li$A2, dots = factor(sample_data(Air.Genus.rel.table)$NewVar))
# calculate group centroid locations
centroids.all.air <- aggregate(cbind(x2,y2)~cluster,data= gg.all.air,mean)
# merge centroid locations into ggplot dataframe
gg.all.air <- merge(gg.all.air, centroids.all.air,by="cluster",suffixes=c("",".centroid2"))


pdf(file="Air.All.PCA_ggplot.pdf", height = 4, width = 5) 
ggplot(gg.all.air, aes(x = gg.all.air $x, y = gg.all.air $y2, label = cluster, color = cluster)) +
  #sets backgrounds
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2,color= cluster)) +

  #sets centroid points and draws lines to data points
  geom_point(data= centroids.all.air, aes(x=x2, y=y2, color=cluster), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
  #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  #stat_ellipse(level = 0.5, type = "t", show.legend = FALSE) + 
  scale_color_manual(values =c("grey", "maroon", "deepskyblue","blue")) + 
  #labels centroids 
  geom_label(data = centroids.all.air, aes (x=x2, y=y2))+ 
  #remove legend 
  theme(legend.position="none")
  
 dev.off()
 

###Just work site air for three sites 
#subset samples
Ad.Assem.Mach.Air.genus.relative.table = subset_samples(Air.Genus.rel.table, Sample_Type_d1_code %in% c("41", "42", "43"))
Ad.Assem.Mach.Air.Bray.dist = distance(Ad.Assem.Mach.Air.genus.relative.table, method="bray")
Ad.Assem.Mach.Air.Bray.pco = dudi.pco(cailliez(Ad.Assem.Mach.Air.Bray.dist))

##Plot s.class and adonis
pdf(file="Ad.Assem.Mach.Air.PCoA.Bray.pdf", width=8, height=8)
s.class(Ad.Assem.Mach.Air.Bray.pco$li, interaction(sample_data(Ad.Assem.Mach.Air.genus.relative.table)$Sample_Type_d1), col=c("palegreen", "orange", "purple"))
dev.off()

adonis(Ad.Assem.Mach.Air.Bray.dist ~ Sample_Type_d1_code, data=data.frame(sample_data(Ad.Assem.Mach.Air.genus.relative.table)))
#Call:
adonis(formula = Ad.Assem.Mach.Air.Bray.dist ~ Sample_Type_d1_code,      data = data.frame(sample_data(Ad.Assem.Mach.Air.genus.relative.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sample_Type_d1_code   1    0.2418 0.24182  1.1473 0.00917   0.27
Residuals           124   26.1359 0.21077         0.99083       
Total               125   26.3777                 1.00000 

#ordinate
Used.Air.ordinate.bray <- ordinate(Ad.Assem.Mach.Air.genus.relative.table, method="PCoA", distance="bray")
#Scatter plot Axis1 [17.2%] v Axis 2 [9.7%]
plot_ordination(Ad.Assem.Mach.Air.genus.relative.table, Used.Air.ordinate.bray, color = "Sample_Type_d1", shape = "Sample_Type_d1")


#### Convert to ggplot
# build ggplot dataframe with points (x,y), corresponding groups (cluster), and location for dots
gg.used.air <- data.frame(cluster=factor(sample_data(Ad.Assem.Mach.Air.genus.relative.table)$Sample_Type_d1), x2= Ad.Assem.Mach.Air.Bray.pco $li$A1, y2= Ad.Assem.Mach.Air.Bray.pco $li$A2, dots = factor(sample_data(Ad.Assem.Mach.Air.genus.relative.table)$Sample_Type_d1))
# calculate group centroid locations
centroids.used.air <- aggregate(cbind(x2,y2)~cluster,data= gg.used.air,mean)
# merge centroid locations into ggplot dataframe
gg.used.air <- merge(gg.used.air, centroids.used.air,by="cluster",suffixes=c("",".centroid2"))


pdf(file="Air.Used.PCA_ggplot.pdf", height = 4, width = 5) 
ggplot(gg.used.air, aes(x = gg.used.air $x, y = gg.used.air $y2, label = cluster, color = cluster)) +
  #sets backgrounds
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2,color= cluster)) +

  #sets centroid points and draws lines to data points
  geom_point(data= centroids.used.air, aes(x=x2, y=y2, color=cluster), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
  #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  #stat_ellipse(level = 0.5, type = "t", show.legend = FALSE) + 
  scale_color_manual(values =c("palegreen", "orange", "purple")) + 
  #labels centroids 
  geom_label(data = centroids.used.air, aes (x=x2, y=y2))+ 
  #remove legend 
  theme(legend.position="none")
  
 dev.off()



}
}

##Supplemental Figure 3 - PCOA plots processed in R as follows 
{
### Lung Tissue
{
#subset lung samples
lung.otu.rel.table = subset_samples(Genus.rel.table, Sample_Type_d_fluid_change_code == 101 | Sample_Type_d_fluid_change_code == 102)

 ## # Calculate Bray distance
lung.bray.dist = distance(lung.otu.rel.table, method="bray")
lung.bray.pco = dudi.pco(cailliez(lung.bray.dist))

##Plot PCoA Lung
pdf(file="Lung.PCoA.Bray.pdf", width=10, height=10)
s.class(lung.bray.pco $li, interaction(sample_data(lung.otu.rel.table)$Sample_Type_d_fluid_change_code), col=c("red", "blue"))
dev.off()

# Calculate ordinate bray distance
ordinate.lung.bray.dist <- ordinate(lung.otu.rel.table, method="PCoA", distance="bray")

# Scatter plot
plot_ordination(lung.otu.rel.table, ordinate.lung.bray.dist) #Axis1 [29.2%], #Axis2 [17.7%]

#check adonis
adonis(lung.bray.dist ~ Sample_Type_d_fluid_change_code, data=data.frame(sample_data(lung.otu.rel.table)))
Call:
adonis(formula = lung.bray.dist ~ Sample_Type_d_fluid_change_code,      data = data.frame(sample_data(lung.otu.rel.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sample_Type_d_fluid_change_code  1    0.2932 0.29324  1.8383 0.09267  0.048 *
Residuals                       18    2.8712 0.15951         0.90733         
Total                           19    3.1645                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#### Convert to ggplot for easier plotting
# build ggplot dataframe with points (x,y), corresponding groups (cluster), and location for dots
gg.lung <- data.frame(cluster=factor(sample_data(lung.otu.rel.table)$Sample_Type_d_fluid_change_code), x2= lung.bray.pco $li$A1, y2= lung.bray.pco $li$A2, dots = factor(sample_data(lung.otu.rel.table)$Sample_Type_d_fluid_change_code))
# calculate group centroid locations
centroids.lung <- aggregate(cbind(x2,y2)~cluster,data= gg.lung,mean)
# merge centroid locations into ggplot dataframe
gg.lung <- merge(gg.lung, centroids.lung,by="cluster",suffixes=c("",".centroid2"))

# plot with the lines to the points from the centroid
pdf(file="LungPCA_ggplot.b.pdf", height = 4, width = 5) 
ggplot(gg.lung, aes(x = gg.lung $x2, y = gg.lung $y2, color =cluster)) +
  #sets background
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2, color = dots), size = 2) +
 

  #sets centroid points and draws lines to data points
  geom_point(data= centroids.lung, aes(x=x2, y=y2), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
   #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  stat_ellipse(type = "norm") +

  #this sets the color scale all color inputs 
  scale_color_manual(values = c(
    	"102"="blue", 
    	"101"="red" 
    	)) +
  #labels centroids 
  geom_label(data = centroids.lung, aes (x=x2, y=y2, label = c("101"="Case", "102" ="Control")))
  
  dev.off()
}

### Skin 
{
#subset samples	
Skin.Human.otu.relative.table = subset_samples(Genus.rel.table, Sample_Type_c_code ==22)

## Calculate Bray distance
Skin.Human.Bray.dist = distance(Skin.Human.otu.relative.table, method="bray")
Skin.Human.Bray.pco = dudi.pco(cailliez(Skin.Human.Bray.dist))

##Plot PCoA
pdf(file="Skin.Human.PCoA.Bray.b.pdf", width=10, height=10)
s.class(Skin.Human.Bray.pco$li, interaction(sample_data(Skin.Human.otu.relative.table)$Sample_Type_d_fluid_change_c), col=c("palegreen", "orange", "purple"))
dev.off()

#Check pval from anova 
adonis(Skin.Human.Bray.dist ~ Sample_Type_d_fluid_change_c, data=data.frame(sample_data(Skin.Human.otu.relative.table)))

#Cal
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sample_Type_d_fluid_change_c   2     1.101 0.55048  2.2663 0.01508  0.001 ***
Residuals                    296    71.897 0.24290         0.98492           
Total                        298    72.998                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Get Axis via Ordination
Skin.Human.ordinate.bray <- ordinate(Skin.Human.otu.relative.table, method="PCoA", distance="bray")

# Scatter plot #Axis1 [15.9%] #Axis2 [12.8%]
pdf(file="Skin.Ordinate.PCoA.Bray.pdf", width=10, height=10) 
plot_ordination(Skin.Human.otu.relative.table, Skin.Human.ordinate.bray, color = "Sample_Type_d_fluid_change_c", shape = "Sample_Type_d_fluid_change_c")
dev.off()


###Plotting
# build ggplot dataframe with points (x,y), corresponding groups (cluster)
gg.skin <- data.frame(cluster=factor(sample_data(Skin.Human.otu.relative.table)$Sample_Type_d_fluid_change_c), x2= Skin.Human.Bray.pco $li$A1, y2= Skin.Human.Bray.pco $li$A2, dots = factor(sample_data(Skin.Human.otu.relative.table)$Sample_Type_d_fluid_change_c))
# calculate group centroid locations
centroids.skin <- aggregate(cbind(x2,y2)~cluster,data= gg.skin,mean)
# merge centroid locations into ggplot dataframe
gg.skin <- merge(gg.skin, centroids.skin,by="cluster",suffixes=c("",".centroid"))


# plot with the lines to the points from the centroid
pdf(file="SkinPCA_ggplot.b.pdf", height = 4, width = 4) 
ggplot(gg.skin, aes(x = gg.skin$x2, y = gg.skin$y2, color =cluster)) +
  #sets background
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(legend.position="none") + 
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2, color = dots), size = 2) +
  #sets centroid points and draws lines to data points
  geom_point(data= centroids.skin, aes(x=x2, y=y2), size=0) +
  geom_segment(aes(x=x2.centroid, y=y2.centroid, xend=x2, yend=y2, color=cluster))+ 
   #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  stat_ellipse(type = "norm") +
  #this sets the color scale all color inputs 
  scale_color_manual(values = c("orange",  "purple","palegreen")) +
  #labels centroids 
  geom_label(data = centroids.skin, aes (x=x2, y=y2, label = c("Skin.Assembly" = "Assembly", "Skin.Machine.Shop" = "Machine Shop", "Skin.Administration"="Administration")))
dev.off()



	
}	
  
### Oral
{ 
#select oral samples
Oral.Human.otu.relative.table = subset_samples(Genus.rel.table, Sample_Type_c_code ==21)

## # Calculate Bray distance
Oral.Human.Bray.dist = distance(Oral.Human.otu.relative.table, method="bray")
Oral.Human.Bray.pco = dudi.pco(cailliez(Oral.Human.Bray.dist))

##Plot PCoA
pdf(file="Oral.Human.PCoA.Bray.b.pdf", width=10, height=10)
s.class(Oral.Human.Bray.pco$li, interaction(sample_data(Oral.Human.otu.relative.table)$Sample_Type_d_fluid_change_c), col=c("palegreen", "orange", "purple"))
dev.off()

#calc adonis
adonis(Oral.Human.Bray.dist ~ Sample_Type_d_fluid_change_c, data=data.frame(sample_data(Oral.Human.otu.relative.table)))

Terms added sequentially (first to last)

                              Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
Sample_Type_d_fluid_change_c   2    0.1547 0.077340  1.2122 0.0081  0.296
Residuals                    297   18.9492 0.063802         0.9919       
Total                        299   19.1039                  1.0000      

#Get Axis via Ordination
Oral.Human.ordinate.bray <- ordinate(Oral.Human.otu.relative.table, method="PCoA", distance="bray")

# Scatter plot #Axis1 [30.6%] #Axis2 [21.3%]
pdf(file="Oral.Ordinate.PCoA.Bray.pdf", width=10, height=10) 
plot_ordination(Oral.Human.otu.relative.table, Oral.Human.ordinate.bray, color = "Sample_Type_d_fluid_change_c", shape = "Sample_Type_d_fluid_change_c")
dev.off()


##Convert to ggplot 
# build ggplot dataframe with points (x,y), corresponding groups (cluster)
gg.oral <- data.frame(cluster=factor(sample_data(Oral.Human.otu.relative.table)$Sample_Type_d_fluid_change_c), x2= Oral.Human.Bray.pco $li$A1, y2= Oral.Human.Bray.pco $li$A2, dots = factor(sample_data(Oral.Human.otu.relative.table)$Sample_Type_d_fluid_change_c))
# calculate group centroid locations
centroids.oral <- aggregate(cbind(x2,y2)~cluster,data= gg.oral,mean)
# merge centroid locations into ggplot dataframe
gg.oral <- merge(gg.oral,centroids.oral,by="cluster",suffixes=c("",".centroid2"))

# plot with the lines to the points from the centroid
pdf(file="OralPCA_ggplot.b.pdf", height = 4, width = 4) 
ggplot(gg.oral, aes(x = gg.oral $x2, y = gg.oral $y2, color =cluster)) +
  #sets background
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(legend.position="none") +

  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2, color = dots), size = 2) +
 

  #sets centroid points and draws lines to data points
  geom_point(data=centroids.oral, aes(x=x2, y=y2), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
   #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  stat_ellipse(type = "norm") +

  #this sets the color scale all color inputs 
  scale_color_manual(values = c("Oral.Administration"="palegreen", "Oral.Assembly" ="orange","Oral.Machine.Shop" = "purple")) +
  #labels centroids 
  geom_label(data=centroids.oral, aes (x = x2, y = y2, label = c(
    	"Oral.Administration"="Administration",
    	"Oral.Assembly" = "Assembly",
    	"Oral.Machine.Shop" = "Machine Shop" 
    	)))
dev.off()

}
 
 ### Nasal 
 {
 #subset nasal samples
 Nasal.Human.otu.relative.table = subset_samples(Genus.rel.table, Sample_Type_c_code ==20)

##Calculate Bray distance
Nasal.Human.Bray.dist = distance(Nasal.Human.otu.relative.table, method="bray")
Nasal.Human.Bray.pco = dudi.pco(cailliez(Nasal.Human.Bray.dist))

##Plot PCoA
pdf(file="Nasal.Human.PCoA.Bray.b.pdf", width=10, height=10)
s.class(Nasal.Human.Bray.pco$li, interaction(sample_data(Nasal.Human.otu.relative.table)$Sample_Type_d_fluid_change_c), col=c("palegreen", "orange", "purple"))
dev.off()

adonis(Nasal.Human.Bray.dist ~ Sample_Type_d_fluid_change_c, data=data.frame(sample_data(Nasal.Human.otu.relative.table)))
#adonis(formula = Nasal.Human.Bray.dist ~ Sample_Type_d_fluid_change_c,      data = data.frame(sample_data(Nasal.Human.otu.relative.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sample_Type_d_fluid_change_c   2     0.718 0.35921   2.028 0.01338  0.027 *
Residuals                    299    52.962 0.17713         0.98662         
Total                        301    53.680                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Calculate ordinate bray distance
ordinate.nasal.bray.dist <- ordinate(Nasal.Human.otu.relative.table, method="PCoA", distance="bray")

# Scatter plot
plot_ordination(Nasal.Human.otu.relative.table, ordinate.nasal.bray.dist) #Axis1 [27.7%], #Axis2 [20.7%]


##Convert to ggplot 
# build ggplot dataframe with points (x,y), corresponding groups (cluster)
gg.nasal <- data.frame(cluster=factor(sample_data(Nasal.Human.otu.relative.table)$Sample_Type_d_fluid_change_c), x2= Nasal.Human.Bray.pco $li$A1, y2= Nasal.Human.Bray.pco $li$A2, dots = factor(sample_data(Nasal.Human.otu.relative.table)$Sample_Type_d_fluid_change_c))
# calculate group centroid locations
centroids.nasa <- aggregate(cbind(x2,y2)~cluster,data= gg.nasal,mean)
# merge centroid locations into ggplot dataframe
gg.nasal <- merge(gg.nasal,centroids.nasa,by="cluster",suffixes=c("",".centroid2"))


# plot with the lines to the points from the centroid
pdf(file="NasalPCA_ggplot.b.pdf", height = 4, width = 4) 
ggplot(gg.nasal, aes(x = gg.nasal $x2, y = gg.nasal $y2, color =cluster)) +
  #sets background
  xlab(NULL) + ylab(NULL) +
  theme_linedraw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(legend.position="none") +
  #sets point size to pseudomonas genus and scales
  geom_point(aes(x=x2,y=y2, color = dots), size = 2) +
 

  #sets centroid points and draws lines to data points
  geom_point(data=centroids.nasa, aes(x=x2, y=y2), size=0) +
  geom_segment(aes(x=x2.centroid2, y=y2.centroid2, xend=x2, yend=y2, color=cluster))+ 
   #draws ellipse and sets colors, type t is a multivatiate t-distribution, "norm" is multivariate normal
  stat_ellipse(level = 0.9, type = "t") +

  #this sets the color scale all color inputs 
  scale_color_manual(values = c(
    	"Nasal.Administration"="palegreen", "Nasal.Assembly" = "orange", "Nasal.Machine.shop" = "purple")) +
  #labels centroids 
  geom_label(data = centroids.nasa, aes (x=x2, y=y2, label = c(
    	"Nasal.Administration"="Administration", "Nasal.Assembly" = "Assembly","Nasal.Machine.shop" = "Machine Shop")))
dev.off()
}
}

##Supp Fig 4
{
#processed in R and then taken to prism for figure development
#Air v. Lung - this can be repeated for each sample type to produce the figures 
{
air.lung.otu.rel.table = subset_samples(Genus.rel.table, Sample_Type_d_fluid_change_code == "101" | Sample_Type_d_fluid_change_code == "102" | Sample_Type_d_fluid_change_code_c == "41" | Sample_Type_d_fluid_change_code_c == "42"| Sample_Type_d_fluid_change_code_c == "43")

sample_data(air.lung.otu.rel.table)$Sample_Type_d_fluid_change_code
## # Calculate Bray distance

air.lung.bray.dist = distance(air.lung.otu.rel.table, method="bray")
air.lung.bray.pco = dudi.pco(cailliez(air.lung.bray.dist))

#melt the distance matrix into columns with each sample site and distance
b2 <- melt(as.matrix(air.lung.bray.dist))
head(b2)

#Then need to remove self distances and duplicated distances
p2    <- t(apply(b2[,c(1,2)],1,FUN=sort))
rmv2 <- which(p2[,1] == p2[,2])

p2    <- paste(p2[,1],p2[,2],sep="|")
rmv3 <- which(duplicated(p2))

#establish new data frame that removes those values 
b2.df   <- b2[-c(rmv2,rmv3),] 
head(b2.df)

##Now we need to replace variable rows with loca 
#set up new data frame
new2.df <- b2.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new2.df[] <- lapply(b2.df, function(x) sample_data(air.lung.otu.rel.table)$Sample_Type_d_fluid_change_code[match(x, rownames(sample_data(air.lung.otu.rel.table)))])

head(new2.df)

#create two lists of the group variable 
topo.var11 <- new2.df[,1]
topo.var21 <-new2.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b2.var.df <- cbind(b2.df, topo.var11, topo.var21)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b2.var.df <- b2.var.df
#set row names to re-zero
rownames(btw.b2.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b2.var.df)) {
	if (btw.b2.var.df$topo.var11[i] == btw.b2.var.df$topo.var21[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b2.var.df <- btw.b2.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b2.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw2 = paste(btw.b2.var.df$topo.var11, "to", btw.b2.var.df$topo.var21)
new.cat.btw2.df <- data.frame(btw.b2.var.df$topo.var11, btw.b2.var.df$topo.var21)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw2.df, 1, sort))
unique.new.cat.btw2 <- unique(new.cat.btw2.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw2) <- NULL
rownames(unique.new.cat.btw2) <- NULL
unique.new.cat.btw2 <- paste(unique.new.cat.btw2[,1], "to", unique.new.cat.btw2[,2])


#create new data frame 
clean.btw.b2.var.df <- btw.b2.var.df

#reset row names
rownames(clean.btw.b2.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b2.var.df)){
	if (paste(clean.btw.b2.var.df$topo.var21[i], "to", clean.btw.b2.var.df$topo.var11[i]) %in% unique.new.cat.btw2) {
		clean.btw.b2.var.df$topo.var11[i] <- btw.b2.var.df$topo.var21[i]
		clean.btw.b2.var.df$topo.var21[i] <- btw.b2.var.df$topo.var11[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw2.clean = paste(clean.btw.b2.var.df$topo.var11, "to", clean.btw.b2.var.df$topo.var21)

#confirm permutations 
unique(new.cat.btw2.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
inter.air.lung.data.bray <- data.frame(Sample1 = clean.btw.b2.var.df$X1, Sample2 = clean.btw.b2.var.df$X2, Category = new.cat.btw2.clean, Distance = as.numeric(clean.btw.b2.var.df$value))

#because the above also plots between tissue and between mwf categories want to subset out some, select for the ones you want to keep 
Air.Lung.Subset= inter.air.lung.data.bray[c(which(inter.air.lung.data.bray $Category == "101 to 42"), which(inter.air.lung.data.bray $Category == "101 to 41"), which(inter.air.lung.data.bray $Category == "41 to 102"), which(inter.air.lung.data.bray $Category == "101 to 43"), which(inter.air.lung.data.bray $Category == "43 to 102"), which(inter.air.lung.data.bray $Category == "42 to 102")),]

#save as txt file 
write.csv(Air.Lung.Subset, file = "Inter_Bray_Air_Lung.txt")	
#Then take this file to prism for creating the image desired 	
}

}


