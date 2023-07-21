# DATA ANALYSIS AND VISUALIZATION OF METAGENOMIC ANALYSIS - SAMPLES (EXAMPLE)
## Author: Victor Guillermo Cornejo Villanueva

# Set your working directory
setwd("~/Escritorio/R/Curso/input_files")
list.files() #Check files inside the working directory

# Create a phyloseq object with the .qza files
## Load necessary libraries
library(xfun)
library(qiime2R)
library(phyloseq)

## Check R Documentation in reading Artifacts (.qza) with qiime2R
?read_qza

## 1. Read a table of sequence variants
SVs <- read_qza("table.qza")
names(SVs)

### To access the actual data stored within the object
SVs$data[1:5, 1:5] #show first 5 samples and first 5 taxa

#### To look the unique identifier of the object
SVs$uuid

#### To see the type of artifact
SVs$type

#### To obtain a complete list of the files within the artifact and their sizes
SVs$contents

#### To print the providence
print_provenance(SVs)

## 2. Reading Metadata (read_q2metadata)
metadata <- read_q2metadata("metadata.tsv") #2nd line have to start with #q2:types
head(metadata) #show top lines of metadata

## 3. Reading Taxonomy
taxonomy <- read_qza("taxonomy.qza")
head(taxonomy$data)

## 4. Reading Tree
tree <- read_qza("rooted_tree.qza")
tree_file <- tree$data

## Creating a Phyloseq Object
physeq_O <- qza_to_phyloseq(
  features = "table.qza",
  tree = "rooted_tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

physeq_O

# Removing Archaea and Eukaryota
physeq_OF <- subset_taxa(physeq_O, Kingdom != "d__Archaea")
physeq_OF
physeq_OF <- subset_taxa(physeq_OF, Kingdom != "d__Eukaryota")
physeq_OF
tax_table(physeq_OF)
taxa <- tax_table(physeq_OF)
setwd("~/Escritorio/R/Curso")
write.csv(taxa, "taxa_phyloseq.csv")


# 1_Rarefaction Curves
library(vegan)

## Set the file for the results
setwd("~/Escritorio/R/Curso/1_Rarefaction_curves")

tab <- otu_table(physeq_OF)@.Data
class(tab) #it has to be a matrix
rarecurve(t(tab), step = 50, cex = 0.5)
### Save as "Rarefaction_curve_simple.png"

raremax <- min(rowSums(t(tab)))
raremax
col <- c("black", "darkred", "forestgreen", "orange3", "blue", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
rarecurve(t(tab), step = 20, col = col, lty = lty, label = TRUE)
### Save as "Rarefaction_curve_samples_colors.png"

write.csv(tab, "sample_size_per_species.csv", row.names = FALSE)
## Import and Excel with the sample size per groups
tab_groups <- read_excel('sample_size_per_groups.xlsx') 
raremax2 <- min(rowSums(t(tab_groups)))
raremax2
rarecurve(t(tab_groups), step = 20, col = col, lty = lty, label = TRUE)
abline(v = raremax2)
### Save as "Rarefaction_curve_groups.png"


# 2_Alpha Diversity
library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gplots)

setwd("~/Escritorio/R/Curso/input_files")

shannon <- read_qza("shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("ID") #to move the sample names to a new column

## Set the file for the results
setwd("~/Escritorio/R/Curso/2_Alpha_diversity")

## Save the Shannon output
capture.output(shannon, file = "Shannon_entropy.txt")

## Check if there's an extra sample not assigned to a Shannon diversity value
gplots::venn(list(metadata = metadata$SampleID, shannon = shannon$ID)) 
### Save as "venn_samplesxmetadata.png"


## Estimate Richness
richness <- estimate_richness(physeq_OF)
head(richness)
write.table(richness, file = "richness.tsv", sep = ",", quote = FALSE, row.names = T)

## Plot richness
plot_richness(physeq_OF)
### Save as "plot_richness_general.png"

a_my_comparisons_treatment <- list( c("ANDINA-AL", "ANDINA-AY"))

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_OF, x = "Group", color = "Group", measures = c("Observed", "Chao1","Shannon"), title = "Alpha diversity indices of samples") + 
  geom_boxplot(alpha = 0.6) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons_treatment, label = "p.signif", symnum.args = symnum.args)
ggsave("alpha_diversity_treatment.png", height = 7, width = 8, dpi = 400, device = "png")

plot_richness(physeq_OF, x = "Group", color = "Group", measures = c("Observe", "Chao1", "Shannon"), title = "Alpha diversity indices of samples") + 
  geom_boxplot(alpha = 0.6) +
  theme(legend.position = "right", axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 12))
ggsave("alpha_diversity_treatment_without_NS.png", height = 7, width = 8, dpi = 400, device = "png")

## Violin plot for alpha diversity
library(MicrobiotaProcess)

alphaobj <- get_alphaindex(physeq_OF)
head(as.data.frame(alphaobj))

p_alpha_treatment <- ggbox(alphaobj, geom = "violin", factorNames = "Group", p_textsize = 3) +
  scale_fill_manual(values = c("#00AED7", "#FD9347")) +
  theme(strip.background = element_rect(colour = NA, fill = "grey"))
p_alpha_treatment
ggsave("plot_alpha_violin-treatment.png", height = 10, width = 22, device = "png")


## Alpha Plot richness
p1 = plot_richness(physeq_OF, x = "Group", color = "Group", measures = c("Chao1", "Shannon"))
p1 + geom_point(size=5, alpha=0.7)
ggsave("alpha_plot_richness_treatment.png", height = 7, width = 6, dpi = 400, device = "png")


## Create an histogram to check if the Shannon values estimated for the metagenomes are normally distributed
hist(richness$Shannon, main = "Shannon index", xlab = "")
### Save as "histogram_Shannon-index.png"

## Analysis of variance (ANOVA)
anova.sh.treatment = aov(richness$Shannon ~ Group, metadata)
summary_ANOVA_shannon_treatment = summary(anova.sh.treatment)
summary_ANOVA_shannon_treatment
capture.output(summary_ANOVA_shannon_treatment, file = "ANOVA_Shannon_Treatment.txt")

anova.ch1.treatment = aov(richness$Chao1 ~ Group, metadata)
summary_ANOVA_chao1_treatment = summary(anova.ch1.treatment)
summary_ANOVA_chao1_treatment
capture.output(summary_ANOVA_chao1_treatment, file = "ANOVA_Chao1_Treatment.txt")

### Tukey Honest Significant Differences (Tukey HSD)
TukeyHSD_shannon_treatment = TukeyHSD(anova.sh.treatment)
TukeyHSD_shannon_treatment
capture.output(TukeyHSD_shannon_treatment, file = "TukeyHSD_Shannon_Treatment.txt")

TukeyHSD_chao1_treatment = TukeyHSD(anova.ch1.treatment)
TukeyHSD_chao1_treatment
capture.output(TukeyHSD_chao1_treatment, file = "TukeyHSD_Chao1_Treatment.txt")

## Kruskall-Wallis Rank Sum Test (for non-normally distributed data)
Kruskall_treatment = kruskal.test(richness$Shannon ~ Group, metadata)
Kruskall_treatment
capture.output(Kruskall_treatment, file = "Kruskall_Treatment.txt")

## A list with the p-values resulted of the Wilcoxon Tests considering each pair of groups
Wilcoxon_treatment = pairwise.wilcox.test(richness$Shannon, metadata$Group, p.adj = "bonf")
Wilcoxon_treatment
capture.output(Wilcoxon_treatment, file = "Wilcoxon_Treatment.txt")



# 3_Beta diversity
library(vegan)

## Set the file for the results
setwd("~/Escritorio/R/Curso/3_Beta_diversity")

## OTU-based metrics
OTUs = physeq_OF@otu_table@.Data #To use the “Group” column as the row names so that it will match our metadata
OTU_final = as.data.frame(t(OTUs)) #transpose
OTU_final
### (In case you have to delete a column:)
# OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))] # to remove the “label”, “numOTUs”, and “Group” columns as they are not OTU counts like the rest of the table

### Ordination scatterplots (Bray-Curtis metric)
BC.nmds = metaMDS(OTU_final, distance = "bray", k = 2, trymax = 20) #calculate the nMDS values
# check the stress is lower than 0.2

#### Plot the nMDS
par(mfrow = c(1, 1))
##### Create a blank plot for the NMDS
plot(BC.nmds, type = "n", main = "Bray-Curtis")
##### Add the points colored by age
points(BC.nmds, display = "sites", pch = 20, col = c("blue", "green")[metadata$Group])
##### Add a legend
legend(-0.75, 0.38, legend = c("ANDINA-AL","ANDINA-AY"), col = c("blue","green"), pch = 20)
###### Save as "OTU_based_NMDS_BC_treatment.png"

### Ordination scatterplots (Jaccard metric)
J.nmds = metaMDS(OTU_final, distance = "jaccard", k = 2, trymax = 1000)

plot(J.nmds, type = "n", main = "Jaccard")
points(J.nmds, display = "sites", pch = 20, col = c("blue", "green")[metadata$Group])
legend(-0.65, 0.32, legend = c("ANDINA-AL","ANDINA-AY"), col = c("blue","green"), pch = 20)
#### Save as "OTU_based_NMDS_J_treatment.png"


### To create two separate 2D plots to show the 3D data
library(plotly)

#### Calculate the Bray-Curtis nMDS for 3-axis
BC.nmds.3D = metaMDS(OTU_final, distance = "bray", k = 3, trymax = 1000)
##### Extract x-y-z values for this nMDS
BCxyz = scores(BC.nmds.3D, display="sites")
######This is a table that looks like 
BCxyz

par(mfrow = c(1,2))
#### Axis 1 and 2 (x and y)
plot(BCxyz[,1], BCxyz[,2], main = "Bray-Curtis 1:2", pch = 20, col = c("blue", "green")[metadata$Group])
#### Axis 1 and 3 (x and z)
plot(BCxyz[,1], BCxyz[,3], main = "Bray-Curtis 1:3", pch = 20, col = c("blue", "green")[metadata$Group])
##### Save as "3Din2D_OTU_based_NMDS_BC_treatment.png"


J.nmds.3D = metaMDS(OTU_final, distance = "jaccard", k = 3, trymax = 1000)
Jxyz = scores(J.nmds.3D, display = "sites")
Jxyz

plot(Jxyz[,1], Jxyz[,2], main = "Jaccard 1:2", pch = 20, col = c("blue", "green")[metadata$Group])
plot(Jxyz[,1], Jxyz[,3], main = "Jaccard 1:3", pch = 20, col = c("blue", "green")[metadata$Group])
##### Save as "3Din2D_OTU_based_NMDS_J_treatment.png"


## Phylogenetic-based metrics

### Ordination scatterplots
#### Weighted unifrac
wUF.ordu = ordinate(physeq_OF, method = "NMDS", distance = "unifrac", weighted = TRUE)

plot_ordination(physeq_OF, wUF.ordu, type = "sites", color = "Group") + 
  scale_colour_manual(values = c("ANDINA-AL" = "green", "ANDINA-AY" = "red")) + 
  theme_bw() + 
  stat_ellipse() +
  ggtitle("Weighted UniFrac")
ggsave("Phylo_based_wUF_NMDS_treatment.png", height = 7, width = 8, dpi = 400, device = "png")

### Unweighted unifrac
uwUF.ordu = ordinate(physeq_OF, method = "NMDS", distance = "unifrac", weighted = FALSE)

plot_ordination(physeq_OF, uwUF.ordu, type = "sites", color = "Group") + 
  scale_colour_manual(values = c("ANDINA-AL" = "green", "ANDINA-AY" = "red")) + 
  stat_ellipse() +
  theme_bw() + 
  ggtitle("Unweighted UniFrac")
ggsave("Phylo_based_uwUF_NMDS_treatment.png", height = 7, width = 8, dpi = 400, device = "png")

## Statistically test beta-diversity

### PERMANOVA (permutational analysis of variance)
#### Calculate distance and save as a matrix
BC.dist = vegdist(OTU_final, distance = "bray")
#### Run PERMANOVA on distances
PERMANOVA_BC = adonis2(BC.dist ~ Group, data = metadata, permutations = 1000)
PERMANOVA_BC
capture.output(PERMANOVA_BC, file = "PERMANOVA_BC.txt")

J.dist = vegdist(OTU_final, distance = "jaccard")
PERMANOVA_J = adonis2(J.dist ~ Group, data = metadata, permutations = 1000)
PERMANOVA_J
capture.output(PERMANOVA_J, file = "PERMANOVA_J.txt")


wUF.dist = UniFrac(physeq_OF, weighted = TRUE, normalized = TRUE)
PERMANOVA_wUF = adonis2(wUF.dist ~ Group, data = metadata, permutations = 1000)
PERMANOVA_wUF
capture.output(PERMANOVA_wUF, file = "PERMANOVA_wUF.txt")

uwUF.dist = UniFrac(physeq_OF, weighted = FALSE, normalized = TRUE)
PERMANOVA_uwUF = adonis2(uwUF.dist ~ Group, data = metadata, permutations = 1000)
PERMANOVA_uwUF
capture.output(PERMANOVA_uwUF, file = "PERMANOVA_uwUF.txt")

### ANOSIM (Analysis of similarities)
ANOSIM_BC_Treatment = anosim(BC.dist, metadata$Group, permutations = 1000)
ANOSIM_BC_Treatment
capture.output(ANOSIM_BC_Treatment, file = "ANOSIM_BC_Treatment.txt")

ANOSIM_J_Treatment = anosim(J.dist, metadata$Group, permutations = 1000)
ANOSIM_J_Treatment
capture.output(ANOSIM_J_Treatment, file = "ANOSIM_J_Treatment.txt")

ANOSIM_wUF_Treatment = anosim(wUF.dist, metadata$Group, permutations = 1000)
ANOSIM_wUF_Treatment
capture.output(ANOSIM_wUF_Treatment, file = "ANOSIM_wUF_Treatment.txt")

ANOSIM_uwUF_Treatment = anosim(uwUF.dist, metadata$Group, permutations = 1000)
ANOSIM_uwUF_Treatment
capture.output(ANOSIM_uwUF_Treatment, file = "ANOSIM_uwUF_Treatment.txt")


### Beta dispersion
disp.BC.treatment = betadisper(BC.dist, metadata$Group)
test.disp.BC.treatment = permutest(disp.BC.treatment, pairwise = TRUE, permutations = 1000)
capture.output(test.disp.BC.treatment, file = "Beta_dispersion_BC_Treatment.txt")


disp.J.treatment = betadisper(J.dist, metadata$Group)
test.disp.J.treatment = permutest(disp.J.treatment, pairwise = TRUE, permutations = 1000)
capture.output(test.disp.J.treatment, file = "Beta_dispersion_J_Treatment.txt")


disp.wUF.treatment = betadisper(wUF.dist, metadata$Group)
test.disp.wUF.treatment = permutest(disp.wUF.treatment, pairwise = TRUE, permutations = 1000)
capture.output(test.disp.wUF.treatment, file = "Beta_dispersion_wUF_Treatment.txt")


disp.uwUF.treatment = betadisper(uwUF.dist, metadata$Group)
test.disp.uwUF.treatment = permutest(disp.uwUF.treatment, pairwise = TRUE, permutations = 1000)
capture.output(test.disp.uwUF.treatment, file = "Beta_dispersion_uwUF_Treatment.txt")

### To realize the calculation with Euclidean method and Bonferroni correction
library(pairwiseAdonis)
pairwise.adonis_BC_Euc_Bonf_treatment <- pairwise.adonis(BC.dist, metadata$Group, sim.method = "euclidean",
                                                                p.adjust.m = "bonferroni")
pairwise.adonis_BC_Euc_Bonf_treatment
capture.output(pairwise.adonis_BC_Euc_Bonf_treatment, file = "(n)pairwise_Treatment_BC_euclidean_bonferroni.txt")


pairwise.adonis_J_Euc_Bonf_treatment <- pairwise.adonis(J.dist, metadata$Group, sim.method = "euclidean",
                                                               p.adjust.m = "bonferroni")
pairwise.adonis_J_Euc_Bonf_treatment
capture.output(pairwise.adonis_J_Euc_Bonf_treatment, file = "(n)pairwise_Treatment_J_euclidean_bonferroni.txt")


pairwise.adonis_wUF_Euc_Bonf_treatment <- pairwise.adonis(wUF.dist, metadata$Group, sim.method = "euclidean",
                                                                 p.adjust.m = "bonferroni")
pairwise.adonis_wUF_Euc_Bonf_treatment
capture.output(pairwise.adonis_wUF_Euc_Bonf_treatment, file = "(n)pairwise_Treatment_wUF_euclidean_bonferroni.txt")


pairwise.adonis_uwUF_Euc_Bonf_treatment <- pairwise.adonis(uwUF.dist, metadata$Group, sim.method = "euclidean",
                                                                  p.adjust.m = "bonferroni")
pairwise.adonis_uwUF_Euc_Bonf_treatment
capture.output(pairwise.adonis_uwUF_Euc_Bonf_treatment, file = "(n)pairwise_Treatment_uwUF_euclidean_bonferroni.txt")


#### With the best method, to realize the calculation without method or standard corrections
pairwise.adonis_uwUF_treatment <- pairwise.adonis(uwUF.dist, metadata$Group)
pairwise.adonis_uwUF_treatment
capture.output(pairwise.adonis_uwUF_treatment, file = "pairwise_uwUF_Treatment_default.txt")

#### To realize the calculation without method and with Bonferroni correction
pairwise.adonis_uwUF_Bonf_treatment <- pairwise.adonis(uwUF.dist, metadata$Group,
                                                              p.adjust.m = "bonferroni")
pairwise.adonis_uwUF_Bonf_treatment
capture.output(pairwise.adonis_uwUF_Bonf_treatment, file = "pairwise_uwUF_Treatment_bonferroni.txt")

# To realize the calculation with BC method and Bonferroni correction
pairwise.adonis_uwUF_BC_Bonf_treatment <- pairwise.adonis(uwUF.dist, metadata$Group, sim.method = "bray",
                                                                 p.adjust.m = "bonferroni")
pairwise.adonis_uwUF_BC_Bonf_treatment
capture.output(pairwise.adonis_uwUF_BC_Bonf_treatment, file = "pairwise_uwUF_Treatment_bray_bonferroni.txt")

#Constrained Ordinations -1
# Remove data points with missing metadata
library(dplyr) #To use the function %>%
erie_not_na <- physeq_OF %>%
  subset_samples(
    !is.na(Group)
  )

bray_not_na <- phyloseq::distance(physeq = erie_not_na, method = "unifrac")

# CAP ordinate
#An alternative (more PERMANOVA-like) would be to use a canonical analysis of principal coordinates (CAP) to find axes that best 
#discriminate among groups of interest (like a distance-based discriminant analysis) and then overlay Spearman rank correlation 
#vectors of species abundances with the 2 canonical axes. This is also a nice way of graphically looking at patterns.

cap_ord <- ordinate(
  physeq = erie_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ Group
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = erie_not_na, 
  ordination = cap_ord, 
  color = "Group", 
  axes = c(1,2)
) + 
  aes(shape = Group) + 
  geom_point(aes(colour = Group), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("blue3","red1")) + theme_bw()
cap_plot

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = MDS1, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.4 * CAP1, 
                 y = 1.4 * MDS1, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap2_bray <- cap_plot +
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )
cap2_bray
ggsave(file="beta_diversity_cap_uwUF.png", plot = cap2_bray, width = 20, height = 9, dpi = 600, limitsize = FALSE)


## Additional 2D plot beta diversity
library(microbial)
plotbeta(physeq_OF, group = "Group", distance = "bray", method = "PCoA", ellipse = TRUE)
ggsave("plot_beta_Bray_PCoA-treatment.png", height = 9, width = 10, device = "png")

ordination = ordinate(physeq_OF, "DPCoA", "bray")
plot_ordination(physeq_OF, ordination, "Group", color ="Group")
ggsave("plot_beta_ordination_DPCoA_bray_treatment.png", height = 5, width = 5, device = "png")


# Taxonomic barplots with ABSOLUTE ABUNDANCE (PART 1)
setwd("~/Escritorio/R/Curso/input_files")
SVs <- physeq_OF@otu_table@.Data
taxonomy <- read_qza("taxonomy.qza")$data %>% parse_taxonomy()
taxasums <- summarize_taxa(SVs, taxonomy)$Genus
taxasums2 <- summarize_taxa(SVs, taxonomy)$Phylum

setwd("~/Escritorio/R/Curso/4_Taxonomic_barplot/Absolute_abundance")

taxa_barplot(taxasums, metadata, "Group", normalize = "none")
ggsave("abundance_barplot_treatment.png", height = 8, width = 11, device = "png") # save a PNG 4 inches by 8 inches

setwd("~/Escritorio/R/Curso/4_Taxonomic_barplot")
write.csv(taxasums, "Taxasums_genus.csv")
write.csv(taxasums2, "Taxasums_phylum.csv")

library(readxl)
taxasums3 <- read_excel('Taxasums_treatment_genus.xlsx')
taxasums_genus <- taxasums3 %>%
  tibble::column_to_rownames("Taxa")

taxasums4 <- read_excel('Taxasums_treatment_phylum.xlsx')
taxasums_phylum <- taxasums4 %>%
  tibble::column_to_rownames("Taxa")

MyPalette3 <- c("dodgerblue3","green4","red1","mediumpurple1","yellow2",
                "blueviolet", "chartreuse3", "cadetblue2","chocolate1", 
                "burlywood2", "red3" ,"blue3", "brown3", "cadetblue4", 
                "bisque3","coral", "aquamarine3", "brown","deeppink3", 
                "gold3","darkturquoise","darkslateblue", 
                "darksalmon", "forestgreen", "pink2", "aliceblue",
                "wheat3", "plum", "olivedrab1", "grey81",
                "orangered", "turquoise1", "magenta2", "red",
                "coral3","springgreen1", "orchid1","yellowgreen",
                "darkgoldenrod1","chocolate4", "royalblue1","tan4",
                "darkorange","darkkhaki","pink","firebrick1","seagreen1",
                "peru","lightgoldenrod1")

setwd("~/Escritorio/R/Curso/4_Taxonomic_barplot/Absolute_abundance")

absolute_genus_35 = taxa_barplot(taxasums_genus, normalize = "none", ntoplot = 35)
absolute_genus_35_figure = absolute_genus_35 + scale_fill_manual(values = MyPalette3) 
absolute_genus_35_figure
ggsave("absolute_35genus_treatment.png", height = 7, width = 14, device = "png")

absolute_phylum_10 = taxa_barplot(taxasums_phylum, normalize = "none", ntoplot = 10)
absolute_genus_10_figure = absolute_phylum_10 + scale_fill_manual(values = MyPalette3) 
absolute_genus_10_figure
ggsave("absolute_10phylum_treatment.png", height = 5, width = 4, device = "png")

absolute_genus_45 = taxa_barplot(taxasums, normalize = "none", ntoplot = 45)
absolute_genus_45_figure = absolute_genus_45 + scale_fill_manual(values = MyPalette3) 
absolute_genus_45_figure
ggsave("absolute_45genus_sample.png", height = 7, width = 22, device = "png")


# Taxonomic barplots with ABSOLUTE ABUNDANCE (PART 2)
prune.dat <- prune_taxa(taxa_sums(physeq_OF) > 2, physeq_OF) #remove less than 2%

dat.aglo = tax_glom(prune.dat, taxrank = "Genus")
dat.dataframe = psmelt(dat.aglo)

ggplot(dat.dataframe, aes(x = Group, y = Abundance, fill = Genus)) + geom_bar(stat = "identity", position = "stack")

## Obtain only the TOP 10 results
top13 <- names(sort(taxa_sums(prune.dat), decreasing=TRUE)[1:13])
top13 #shows 10 results

prune.dat.two = prune_taxa(top13, dat.aglo)
dat.dataframe13 = psmelt(prune.dat.two)

ggplot(dat.dataframe13, aes(x = Group, y = Abundance, fill = Genus)) + geom_bar(stat = "identity", position = "stack")
ggsave("absolute_barplot_genus_treatment.png", height = 6, width = 7, dpi = 400, device = "png")


top260 <- names(sort(taxa_sums(prune.dat), decreasing = TRUE)[1:260])
dat.aglo2 = tax_glom(prune.dat, taxrank = "Phylum")
prune.dat.three = prune_taxa(top260, dat.aglo2)
dat.dataframe2 = psmelt(prune.dat.three)
ggplot(dat.dataframe2, aes(x = Group, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity", position = "stack")
ggsave("absolute_barplot_phylum_treatment.png", height = 6, width = 7, dpi = 400, device = "png")


# Taxonomic barplots with RELATIVE ABUNDANCE (with 'Microbial' library)
setwd("~/Escritorio/R/Curso/4_Taxonomic_barplot/Relative_abundance")

library(microbial)

## default normalize method is relative
plotbar(physeq_OF,level = "Phylum", top = 30, return = FALSE, fontsize.x = 10, fontsize.y = 10)
ggsave("plot_bar-phylum-all.png", height = 9, width = 10, device = "png")

plotbar(physeq_OF,level = "Phylum", top = 10, group = "Group", return = FALSE, fontsize.x = 10, fontsize.y = 10)
ggsave("plot_bar-phylum-treatment_10.png", height = 8, width = 7, device = "png")

plotbar(physeq_OF,level = "Genus", top = 10, fontsize.x = 10, fontsize.y = 10)
ggsave("plot_bar-genus-all.png", height = 9, width = 10, device = "png")

plotbar(physeq_OF,level = "Genus", top = 45, fontsize.x = 10, fontsize.y = 10)
ggsave("plot_bar-genus-all_45.png", height = 9, width = 14, device = "png")

plotbar(physeq_OF,level = "Genus", group = "Group", top = 35, fontsize.x = 10, fontsize.y = 10)
ggsave("plot_bar-genus-treatment_35.png", height = 9, width = 14, device = "png")


# Taxonomic barplots with RELATIVE ABUNDANCE (PART 2)
prune.dat <- prune_taxa(taxa_sums(physeq_OF) > 2, physeq_OF) #remove less than 2%

dat.aglo = tax_glom(prune.dat, taxrank = "Genus")
dat.dataframe = psmelt(dat.aglo)

ggplot(dat.dataframe, aes(x = Group, y = Abundance, fill = Genus)) + geom_bar(stat = "identity", position = "fill")

## Obtain only the TOP 10 results
top13 <- names(sort(taxa_sums(prune.dat), decreasing = TRUE)[1:13])
top13 #shows 10 results

prune.dat.two = prune_taxa(top13, dat.aglo)
dat.dataframe13 = psmelt(prune.dat.two)

ggplot(dat.dataframe13, aes(x = Group, y = Abundance, fill = Genus)) + geom_bar(stat = "identity", position = "fill")
ggsave("relative_barplot_genus_treatment.png", height = 6, width = 7, dpi = 400, device = "png")

top260 <- names(sort(taxa_sums(prune.dat), decreasing = TRUE)[1:260])
dat.aglo2 = tax_glom(prune.dat, taxrank = "Phylum")
prune.dat.three = prune_taxa(top260, dat.aglo2)
dat.dataframe2 = psmelt(prune.dat.three)
ggplot(dat.dataframe2, aes(x = Group, y = Abundance, fill = Phylum)) + geom_bar(stat = "identity", position = "fill")
ggsave("relative_barplot_phylum_treatment.png", height = 6, width = 7, dpi = 400, device = "png")


# Taxonomic barplots with RELATIVE ABUNDANCE (PART 3)
taxa_barplot(taxasums,metadata,"Group")
ggsave("relative_abundance_barplot_treatment.png", height = 4, width = 10, device = "png") # save a PNG 4 inches by 8 inches


relative_genus_35 = taxa_barplot(taxasums_genus, ntoplot = 35)
relative_genus_35_figure = relative_genus_35 + scale_fill_manual(values = MyPalette3) 
relative_genus_35_figure
ggsave("relative_35genus_treatment.png", height = 7, width = 14, device = "png")

relative_phylum_10 = taxa_barplot(taxasums_phylum, ntoplot = 10)
relative_genus_10_figure = relative_phylum_10 + scale_fill_manual(values = MyPalette3) 
relative_genus_10_figure
ggsave("relative_10phylum_treatment.png", height = 5, width = 4, device = "png")

relative_genus_45 = taxa_barplot(taxasums, ntoplot = 45)
relative_genus_45_figure = relative_genus_45 + scale_fill_manual(values = MyPalette3) 
relative_genus_45_figure
ggsave("relative_45genus_treatment.png", height = 7, width = 22, device = "png")

## That's all! Thank you so much for your attention! 