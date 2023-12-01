# Weeks 13 and 14 Project

library(ggplot2) # install packages
otu <- read.table("HMP_rarefied_otu_with_taxa.xlsx") # contains the data from the rarefied OTU table
metadata <- read.table('HMP_metadata.tsv') # contains the data from the metadata file 
colnames(metadata)    
samples1 <- colnames(metadata)[1:ncol(metadata)-1] # defines all the samples in the OTU table minus the last column
IDs_Keep <- intersect(samples1, rownames(metadata)) # intersecting with metadata
metadata <- metadata[IDs_Keep,] # filters metadata to only keep the intersecting samples 
otu2 <- metadata[IDs_Keep] # filters the OTU table to keep only the intersecting samples
alpha <- read.table('HMP_alpha_diversity.tsv') 
beta <- as.matrix(read.table("HMP_unweighted_unifrac_distance-matrix.tsv"))
otu2$taxonomy <- metadata$taxonomy 
alpha <- alpha [IDs_Keep, ] # filters alpha diversity to keep the samples
as.character(rownames(metadata)) == colnames(otu2)[1:(ncol(otu2)-1)]


combined_alphadata <- metadata # make copy of metadata to work with
alpha <- as.data.frame(alpha) # creating new column of Shannon indecies 
combined_alphadata$shannon <- alpha$shannon

ggplot() + geom_boxplot(data=combined_alphadata, aes(x= , y= shannon))
colnames(alpha)

ggplot() + geom_boxplot(data=combined_alphadata)
males.ix <- metadata$V3 == 'male' # find samples that are male in metadata
males <- alpha[males.ix,] # subset the alpha data 
females.ix <- metadata$V3 == 'female' # find samples that are female in metadata 
females <- alpha[females.ix,] # subset the alpha data 
hist(females$shannon, xlab = "Alpha Diversity", main = "Females") # creates histogram 
hist(males$shannon, xlab = "Alpha Diversity", main = "Males")

# Shapiro-Wilk Normality Test, p values less than 0.1 indicate the data are significantly different from normal distribution 
# Run this test for each group in the covariate of interest
shapiro.test(females$shannon)
shapiro.test(males$shannon)

# Mann-Whitney U Test i.e., Wilcoxon Rank-Sum test 
# Tests for differences, does not require normal distribution 
wilcox.test(females$shannon, males$shannon, na.rm=TRUE) #na.rm = TRUE removes missing data

# Plot Shannon Index by sex 
library(ggplot2)
colnames(alpha2) 
alpha2 <- merge(alpha, metadata)
ggplot(data = alpha2, aes(x=V3, y = shannon)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, aes(color = V3))
names(metadata) <- c("SampleID", "BodySite", "Sex", "Country", "Description", "Age", "Smoker", "BMI", "RUN_PREFIX", "BarcodeSequence", "RUN_CENTER")

# Set all groups available, use bodysite instead of sex
groups <- unique(metadata$BodySite)
# Create empty vectors to store the pair-wise p-values and the groups tested
pw.pvalues <- NULL
pw.names <- NULL

# Set 2 counters, 'i' starts at 1 and goes until one less than the number of groups; 'j' will start at 2, and go up tot he full number of groups
for(i in 1:(length(groups)-1)){
  for(j in (i+1):length(groups)){
    # Use this to pick the groups assigned to 'i' 
    ix.metric.i <- metadata$BodySite == groups[i]
    # use this to pick the groups assigned to 'j'
    ix.metric.j <- metadata$BodySite == groups[j]
    # stores the p-value from the test 
    pvalue <- wilcox.test(alpha[ix.metric.i, "shannon"],
                          alpha[ix.metric.j, "shannon"])$p.value
    # appends the new p-value to the list 
    pw.pvalues <- c(pw.pvalues, pvalue)
    # sets the names of the groups tested
    test.name <- paste(groups[i], "_vs_", groups[j], sep = '')
    # appends the names of the groups tested to the list 
    pw.names <- c(pw.names, test.name)
  }
}
names(pw.pvalues) <- pw.names

# Correct for type I errors, which is rejecting the true null hypothesis (false positive)
# Corrects using 'fdr' (false discovery rate)
fdr.pvalues <- p.adjust(pw.pvalues, 'fdr')

sink("alpha_stats.txt")
is.na(alpha$shannon)
is.na(metadata$BodySite)

cat("\nNumber of samples in each group:\n")
print(table(metadata$BodySite)) # prints a table of the number of samples at each body site 

cat("\nMean Alpha Diversity:\n")
print(tapply(alpha$shannon, metadata$BodySite, mean)) # gets the mean of the alpha diversity at each body siteby using tapply() to apply the mean function across the alpha table subsetted into body site groups 

cat("nMedian Alpha Diveristy:\n")
print(tapply(alpha$shannon, metadata$BodySite, median))
# gets median alpha diversity at each body site 

cat("\nStandard Deviation:\n")
print(tapply(alpha$shannon, metadata$BodySite,sd))

cat("\nPairwise Mann-Whitney-Wilcoxon Tests were performed.\n")
cat("Pairwise p-values are :/n")
print(pw.pvalues)

cat("\nFDR-corrected pairwise p-values are:\n")
print(p.adjust(pw.pvalues, 'fdr'))

sink()

ggplot(data=alpha2, aes(x=BodySite, y= shannon)) +
  geom_boxplot() +
  geom_jitter(width= 0.1, aes(color=BodySite)) + 
  theme_bw()

plot_output <- ggplot(data = alpha2, aes(x = BodySite, y= shannon)) + 
  geom_boxplot()+
  geom_jitter(width = 0.1, aes(color=BodySite))+
  theme_bw()+
  scale_x_discrete(labels=c("ear fold", "vagina", "saliva", "stool", "plaque")) +
  guides(color= "none") # because they are labelled at the x-axis 
pdf("Alpha_Diversity.pdf", height = 4, width = 6)
plot(plot_output)
dev.off()




