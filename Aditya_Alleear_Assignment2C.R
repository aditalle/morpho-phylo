#STUDY OF THE MORPHO GENUS----
#The Morpho genus is a group of butterflies found in Central America, Mexico, and South America. The data for the Morpho genus will be obtained from the BOLD. A markercode will be picked according to counts and all the entries of that markercode will then be clustered. After clustering, OTUs will be formed and then represented in a phylogenetic tree. The phylogenetic tree will show the branch distance through colour amongst the OTUs/species. To check that the branch lengths are properly calculated, a heatmap will be used as well as a phylogenetic tree with highlighted alignment section of the markercode beside each species name in the tree. These visualization should allow any reader to quickly determine the relatedness between the species. 

# Paula's comment: I am just wondering what is the research question specifically? 

#PACKAGES ----
library(tidyverse)
library(Biostrings)
library(ape)
library(DECIPHER)
library(plyr)
library(ggtree)
library(gplots)
#sessionInfo()  #checking session info
#PICKING AND GETTING DATA FOR STUDY----
##Downloading the data from the BOLD
Morpho.BOLD <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Morpho&format=tsv")
#Looking at what we're dealing with 
summary(Morpho.BOLD)

#Next step is to find what marker to use for analysis down the road

#Filtering: 

# Paula's Edit: When I run the original code, I found out the output is only 1 obs. of 1 variable which is different from normal expectation. I think it is because of the installation of the package plyr.It is possible that there are some overlapping function names between the package plyr and other packages. So at this step, I uninstall the package plyr. 

count.by.marker <- Morpho.BOLD %>%
  group_by(markercode) %>%
  summarize(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

#We find that COI-5P to have the largest counts and so we subset the data using the filter function based on that 
Morpho.COI <- Morpho.BOLD %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[ACGT]"))

#ANALYSIS CODE ----
#Converting the COI sequence from str to biostrings variable in order to work with alignment package muscle 
Morpho.COI$nucleotides <- DNAStringSet(Morpho.COI$nucleotides)

#Using the muscle package, alignment is run. The settings used were maxiter at 20 which is the number of iterations to run in the algorithm. And diags was turned on to allow optimization. 
Morpho.COI.alignment <- DNAStringSet(muscle::muscle(Morpho.COI$nucleotides, maxiters = 20, diags1 = TRUE))

#Paula's Edit: I set the Diag1 to TRUE to allow MUSCLE adopting k-mer extension to firstly find the "diagonals" - those short regions of high similarity between two sequences. The number 1 placed right after the option Diag indicates the diagonal optimization takes place in the first iteration. 

#Naming the alignments in order to keep track of aligned sequences and to merge data frames later on. 
names(Morpho.COI.alignment) <- Morpho.COI$processid

#Using the alignment, we can make a distance matrix using the dist.dna function from the ape package. The model used was K80 because it considers the substitution frequence between certain nucleotides. 
DNA.bin.Morpho <- as.DNAbin(Morpho.COI.alignment)
distanceMatrix.Morpho <- dist.dna(DNA.bin.Morpho, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

# Paula's edit: I changed the model of evolution from K81 to TN93 because TN93 not only assume unequal transition/transversion rate but also consider unequal base frequencies in one sequence. 

#Using the distanceMatrix, we can perform clustering using the IdCluster function from the DECIPHER package. The setting for the model was set to neibhbour joining. The cutoff was set at 0.04 because the it is around the optimal cutoff value for the COI gene. 
Morpho.COI.cluster <- IdClusters(distanceMatrix.Morpho, method = "NJ", cutoff= 0.03, showPlot = TRUE, type = "both", verbose = TRUE)

# Paula's edit: For the OTU assignment for both 16S rRNA and COI, I think it is still better to cluster sequences with a similarity threshold of 97%. Lower threshold raises the possibility to clustering sequences from different species to one cluster. 

########################################################

# Paula's Edit: By observing the dendrogram, I found that the dendrogram height is at 0.2 and the height indicates there are outliers in the dataset. So I will remove the outliers and re-do the clustering process. 

# Check the Singleton Cluster

Cluster.Info <- data.frame(Morpho.COI.cluster[[1]])

Cluster.Info <- Cluster.Info %>% 
  mutate (Processid = Morpho.COI$processid) 

Cluster.size <- Cluster.Info %>% 
  group_by(cluster) %>% 
  summarize(n=length(Processid)) %>% 
  arrange(desc(n)) 

Cluster.singleton <- Cluster.size %>% 
  filter(n==1)  

# From here, I know that there are 14 singleton clusters.

# Find the sequence processid in each singleton cluster
Outlier.Seq <- merge(Cluster.singleton, Cluster.Info, by ="cluster",  x.all=TRUE)

# Outlier Removal 
grep("NCM262-10",Morpho.COI$processid)
grep("NCM145-10",Morpho.COI$processid)
grep("GBGL10987-12",Morpho.COI$processid)
grep("GBLN3376-13",Morpho.COI$processid)
grep("GBGL10978-12",Morpho.COI$processid)
grep("GBLN3348-13",Morpho.COI$processid)
grep("GBLN3376-13",Morpho.COI$processid)
grep("GBLN3368-13",Morpho.COI$processid)
grep("NCM148-10",Morpho.COI$processid)
grep("GBLN3352-13",Morpho.COI$processid)
grep("GBLN3325-13",Morpho.COI$processid)
grep("GBLN5120-14",Morpho.COI$processid)
grep("GBLN3426-13",Morpho.COI$processid)
grep("GBGL10979-12",Morpho.COI$processid)
grep("GBLN3305-13",Morpho.COI$processid)

Morpho.COI.no.outlier <- Morpho.COI[-c(62,411,155,370,74,179,375,74,311,348,16,13,378,187,418,62),]

# Check if the outlier is removed. 
length(Morpho.COI.no.outlier$processid)

# New Alignment and Cluster (without outlier)

Morpho.COI.no.outlier.alignment <- DNAStringSet(muscle::muscle(Morpho.COI.no.outlier$nucleotides, maxiters = 20, diags1 = TRUE))

names(Morpho.COI.no.outlier.alignment) <- Morpho.COI.no.outlier$processid

#Using the alignment, we can make a distance matrix using the dist.dna function from the ape package. The model used was K80 because it considers the substitution frequence between certain nucleotides. 
DNA.bin.Morpho.2 <- as.DNAbin(Morpho.COI.no.outlier.alignment)
distanceMatrix.Morpho.2 <- dist.dna(DNA.bin.Morpho.2, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

Morpho.COI.cluster.2 <- IdClusters(distanceMatrix.Morpho.2, method = "NJ", cutoff= 0.03, showPlot = TRUE, type = "both", verbose = TRUE)

# After removing the outlier, we have found that the dendrogram height is effectively reduced to 0.08. 

########################################################
#Paula's Edit: Here, I install the package plyr back for the further work. 
library(plyr)

#Adding processid column
Morpho.COI.cluster.2[[1]]$processid <- Morpho.COI.no.outlier$processid

#Picking one sequence per OTU 
#First, we randomize the rows
OTU.Morpho.all <- Morpho.COI.cluster.2[[1]][sample(nrow(Morpho.COI.cluster.2[[1]])),] 

#Here, I make a list that will store the processids 
OTU.Morpho.list <- c()

#I run a for loop to go through the unique cluster numbers
for (i in unique(OTU.Morpho.all[[1]])) {
  #i make a temporary data frame to store the current i cluster
  Morpho.temp <- subset(OTU.Morpho.all, cluster == i) 
  #i pick the first row of the subsetted data and put it in the list that i made earlier
  OTU.Morpho.list <- c(OTU.Morpho.list, Morpho.temp[1,]$processid)
}

#I take the list and make a data frame using ldplyr 
OTU.Morpho.seq <- ldply(OTU.Morpho.list, data.frame)

#I make sure that the column has no name
names(OTU.Morpho.seq) <- NULL 

#I name the column
colnames(OTU.Morpho.seq) <- "processid"

#Using the column from the new data frame, I left merge the data frame with the original data frame that contains the metadata
OTU.Morpho.tree <- merge(OTU.Morpho.seq, Morpho.COI, by ="processid",  x.all=TRUE)

#removing duplicates of species name in data frame 
OTU.Morpho.tree = OTU.Morpho.tree[!duplicated((OTU.Morpho.tree$species_name)),]

#I run an alignment for the OTUs 
OTU.Morpho.tree.alignment <- DNAStringSet(muscle::muscle(OTU.Morpho.tree$nucleotides, maxiters = 30, diags = TRUE))

#I write the MSA to a fasta which will be used for another function downstream
writeXStringSet(OTU.Morpho.tree.alignment, file="OTU_Morpho_alignment.fasta", format = "fasta", width = 585)

#I give the alignment species name
names(OTU.Morpho.tree.alignment) <- OTU.Morpho.tree$species_name

#I transformed the aligment object so that it can be used with dist.dna
dnabin.OTU.Morpho <- as.DNAbin(OTU.Morpho.tree.alignment)

#I make a distance matrix using the same model as before
distanceMatrix.OTU.Morpho <- dist.dna(dnabin.OTU.Morpho, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#Paula's Edit: I also changed the model to "TN93" here. 

#I used the single method for clustering because it uses single linkage which would be appropriate for phylogenetic tree
hclust.OTU <- hclust(dist(distanceMatrix.OTU.Morpho), method = "single")

#using as.phylo, I transform the hclust object to a tree object
Morpho.treeobj <- as.phylo(hclust.OTU)

#this is the first phylogenetic tree which uses blue and red to display branch length. xlim was used to make sure nothing gets cut off. theme_tree2() puts the x-scale which is evolutionary distance. geom_tiplab is what displays the species name
Morpho.tree.branch <- ggtree(Morpho.treeobj, aes(color=branch.length)) + geom_tiplab(size=3) + geom_tree() + theme_tree2() + scale_color_continuous(low="blue", high="red") + theme(legend.position="bottom") + labs(title = "Morpho Phylogenetic Tree showing variability in branch length") + xlim(NA, 0.15)
Morpho.tree.branch

# Paula's Edit: Here, I change the xlim from 0.2 to 0.15 - to fit the graph better. 

#i made a heatmap that uses a RGB colour scale. Red means 100% similarity and green is 0% similarity. The heatmap.2 function also by default shows a histogram and colour scale which allows a better understanding of the heatmap. The margins setting was used for resizing. 
scaleRYG <- colorRampPalette(c("red","yellow","green"), space = "rgb")(50)
heatmap.2(distanceMatrix.OTU.Morpho, col=scaleRYG, margins = c(9, 9), main="Morpho COI-5P") 
#the following code is for making a phylogenetic tree with aligned sequences shown next to species name. The fasta file is used in the msaplot function. The offset setting was used to move the MSA a bit to the right. The window setting was used to show a select area of the sequence. 
tree.msa <- ggtree(Morpho.treeobj) + geom_tiplab(size=3) + xlim(NA, 0.2) + labs(title = "Morpho Phylogenetic Tree using COI-5P with MSA")

#Paula's Edit: 
# Here I change the xlim value from 0.4 to 0.2 - to fit the graph better. 

Morpho.tree.highlight <- msaplot(tree.msa, fasta = "OTU_Morpho_alignment.fasta", offset = 0.05, window = c(250, 400)) +  theme(legend.position = "right") 
Morpho.tree.highlight

#Paula's Edit: 
# Here I change the offset value and side window value - to fit the graph better.  

#CONCLUSION ---- 
# From the first visual representation, we can see that most of the species in the phylogenetic tree have very similar branch lengths with the exception of two species (Morpho marcus, Morpho achilles) which have the longest branch lengths. The second visual representation which was the heatmap shows that Morpho archilles has the largest dissimilarity when comparing to the other species. The heatmap also helps to confirm that the clades at the top are closely related as the heatmap is more orange than green. The last visual representation just shows the phylogenetic tree with the alignment. From all these figures, we can see that the genus Morpho contains species that are closely related. After looking more into the genus, it is not surprising given that the species are found in the same climate and environment. The genus Morpho is typically found in the warmer climate of South America, Central America, and Mexico. To get a better understanding, the next step would be to look at the subfamily Satyrinae where the Morpho genus belongs. It would be interesting to look at how drastic the diversity would be given that the subfamily is found all over the globe. There is definitely more work that can be done to the study, The first would be to use either whole genome to make the phylogenetic trees but that would be too lengthy of a process. I would like to have used several genetic markers and built several genetic trees. The phylogenetic trees can be combined into one through consensus. 
