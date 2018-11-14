#STUDY OF THE MORPHO GENUS----
#The Morpho genus is a group of butterflies found in Central America, Mexico, and South America. The data for the Morpho genus will be obtained from the BOLD. A markercode will be picked according to counts and all the entries of that markercode will then be clustered. After clustering, OTUs will be formed and then represented in a phylogenetic tree. The phylogenetic tree will show the branch distance through colour amongst the OTUs/species. To check that the branch lengths are properly calculated, a heatmap will be used as well as a phylogenetic tree with highlighted alignment section of the markercode beside each species name in the tree. These visualization should allow any reader to quickly determine the relatedness between the species. 

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
#count.by.marker <- Morpho.BOLD %>%
#group_by(markercode) %>%
 # summarize(n = length(processid)) %>%
  #arrange(desc(n)) %>%
  #print()

count.by.marker <- Morpho.BOLD %>%  #Comment: I do not understand the function of codes you wrote for this part
  group_by(markercode) %>%          #I assume you want to find the frequency of eahc marker, and this is the easiest way to do this. 
  tally(sort = T)                 

#We find that COI-5P to have the largest counts and so we subset the data using the filter function based on that 
Morpho.COI <- Morpho.BOLD %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[ACGT]"))
#ANALYSIS CODE ----
#Converting the COI sequence from str to biostrings variable in order to work with alignment package muscle 
Morpho.COI$nucleotides <- DNAStringSet(Morpho.COI$nucleotides)

#Using the muscle package, alignment is run. The settings used were maxiter at 20 which is the number of iterations to run in the algorithm. And diags was turned on to allow optimization. 
Morpho.COI.alignment <- DNAStringSet(muscle::muscle(Morpho.COI$nucleotides, maxiters = 20, diags = TRUE))

#Naming the alignments in order to keep track of aligned sequences and to merge data frames later on. 
names(Morpho.COI.alignment) <- Morpho.COI$processid

#Using the alignment, we can make a distance matrix using the dist.dna function from the ape package. The model used was K80 because it considers the substitution frequence between certain nucleotides. 
distanceMatrix.Morpho <- dist.dna(as.DNAbin(Morpho.COI.alignment), model = "K81", as.matrix = TRUE, pairwise.deletion = TRUE)

#Using the distanceMatrix, we can perform clustering using the IdCluster function from the DECIPHER package. The setting for the model was set to neibhbour joining. The cutoff was set at 0.04 because the it is around the optimal cutoff value for the COI gene. 
Morpho.COI.cluster <- IdClusters(distanceMatrix.Morpho, method = "NJ", cutoff= 0.04, showPlot = TRUE, type = "both", verbose = TRUE)

#Adding processid column
Morpho.COI.cluster[[1]]$processid <- Morpho.COI$processid

#Picking one sequence per OTU 
#First, we randomize the rows
OTU.Morpho.all <- Morpho.COI.cluster[[1]][sample(nrow(Morpho.COI.cluster[[1]])),] 

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
writeXStringSet(OTU.Morpho.tree.alignment, file="OTU_Morpho_alignment.fas", format = "fasta", width = 585)

#I give the alignment species name
names(OTU.Morpho.tree.alignment) <- OTU.Morpho.tree$species_name

#I transformed the aligment object so that it can be used with dist.dna
dnabin.OTU.Morpho <- as.DNAbin(OTU.Morpho.tree.alignment)

#I make a distance matrix using the same model as before
distanceMatrix.OTU.Morpho <- dist.dna(dnabin.OTU.Morpho, model = "K81", as.matrix = TRUE, pairwise.deletion = TRUE)

#I used the single method for clustering because it uses single linkage which would be appropriate for phylogenetic tree
hclust.OTU <- hclust(dist(distanceMatrix.OTU.Morpho), method = "single")

plot(hclust.OTU)
#Comment: I can not plot the Morpho.treeobj after runnung heatmap2. It only works after I plot hclst.OTU. 
#using as.phylo, I transform the hclust object to a tree object
Morpho.treeobj <- as.phylo(hclust.OTU)


#this is the first phylogenetic tree which uses blue and red to display branch length. xlim was used to make sure nothing gets cut off. theme_tree2() puts the x-scale which is evolutionary distance. geom_tiplab is what displays the species name
Morpho.tree.branch <- ggtree(Morpho.treeobj, aes(color=branch.length)) + 
  geom_tiplab(size=3) + 
  xlim(NA, 0.4) +            #Comment: I adjusted the limit of x axis to make the graph visible for me
  geom_tree() + 
  theme_tree2() + 
  scale_color_continuous(low="blue", high="red") + 
  theme(legend.position="bottom") + 
  labs(title = "Morpho Phylogenetic Tree showing variability in branch length")

Morpho.tree.branch


#i made a heatmap that uses a RGB colour scale. Red means 100% similarity and green is 0% similarity. The heatmap.2 function also by default shows a histogram and colour scale which allows a better understanding of the heatmap. The margins setting was used for resizing. 
scaleRYG <- colorRampPalette(c("red","yellow","green"), space = "rgb")(50)
heatmap.2(distanceMatrix.OTU.Morpho, col=scaleRYG, margins = c(9, 9), main="Morpho COI-5P") 

#Comment: ggtree also has a function to plot heatmap along with phylogeny tree and easy to use.
new.tree <- ggtree(Morpho.treeobj) + 
  geom_tiplab(size=3)+
  xlim(0,0.8)
  
gheatmap(new.tree, distanceMatrix.OTU.Morpho, offset = 0.1, width=1, colnames_angle = "90", colnames_offset_y = -3, font.size = 3, low = "red", high = "green", hjust = 0.5)+
  ylim(-5,30) +
  labs(title = "Morpho COI-5P")

plot(hclust.OTU)
#the following code is for making a phylogenetic tree with aligned sequences shown next to species name. The fasta file is used in the msaplot function. The offset setting was used to move the MSA a bit to the right. The window setting was used to show a select area of the sequence. 
tree.msa <- ggtree(Morpho.treeobj) + 
  geom_tiplab(size=3) + 
  xlim(NA, 1) +    #Comment: I adjusted the limit of x axis to make the graph visible for me
  labs(title = "Morpho Phylogenetic Tree using COI-5P with MSA")

# I also can not see anything for this plot, so I made few changes
Morpho.tree.highlight <- msaplot(tree.msa, fasta = "OTU_Morpho_alignment.fasta", offset = 0.02, window = c(250, 400)) +  theme(legend.position = "right") 
Morpho.tree.highlight

#install.packages("seqinr")
library(seqinr)
write.fasta(OTU.Morpho.tree.alignment,names(OTU.Morpho.tree.alignment),file.out = "OTU_Morpho_alignment.fasta")
# Comment : I get error when I run this code: arguments imply differing number of rows: 1, 151, 0. 
# Comment : I only get figure after I overwrite your fasta file. I use write.fasta function from seqinr package to write the file.
#CONCLUSION ---- 
# From the first visual representation, we can see that most of the species in the phylogenetic tree have very similar branch lengths with the exception of two species (Morpho marcus, Morpho achilles) which have the longest branch lengths. The second visual representation which was the heatmap shows that Morpho archilles has the largest dissimilarity when comparing to the other species. The heatmap also helps to confirm that the clades at the top are closely related as the heatmap is more orange than green. The last visual representation just shows the phylogenetic tree with the alignment. From all these figures, we can see that the genus Morpho contains species that are closely related. After looking more into the genus, it is not surprising given that the species are found in the same climate and environment. The genus Morpho is typically found in the warmer climate of South America, Central America, and Mexico. To get a better understanding, the next step would be to look at the subfamily Satyrinae where the Morpho genus belongs. It would be interesting to look at how drastic the diversity would be given that the subfamily is found all over the globe. There is definitely more work that can be done to the study, The first would be to use either whole genome to make the phylogenetic trees but that would be too lengthy of a process. I would like to have used several genetic markers and built several genetic trees. The phylogenetic trees can be combined into one through consensus.
