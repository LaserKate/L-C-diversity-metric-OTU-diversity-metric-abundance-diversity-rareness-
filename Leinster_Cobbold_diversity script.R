
#This code calculates the diversity of a community based on OTU similarity and abundance.
#It was originally written for species and has been adapted for OTUs, including DNA sequence similarity for OTUs using the package ape
#It also includes rarity, as the parameter q, where the importance of the rare species or OTU is q=0 (rare) to q->infinity (very abundant)

#Please direct questions to: K.Quigley (katemarie.quigley@my.jcu.edu.au)

#Citations----
#(1)
#Quigley, K. M., Warner, P. A., Bay, L. K., & Willis, B. L. (2018). 
#Unexpected mixed-mode transmission and moderate genetic regulation of Symbiodinium communities in a brooding coral. 
#Heredity, 121(6), 524-536.

#(2)
#Leinster, Tom, and Christina A. Cobbold. 
#"Measuring diversity: the importance of species similarity." 
#Ecology 93.3 (2012): 477-489.


###################################################
###1) OTUs aka species matrix----
###################################################

#observations = species abundances in samples, comma delimited text file (.csv)
#species names as column headersplot
#example dataset
obs <-read.csv("DiversityNormalizedAbundance.csv", header=T)
#fyi, there are no row names in this initial abundance file. The script is ASSUMING that the order of your OTUs is the same as what you provide later on, so make sure it is. 

samples = length(obs[1,]) #names of samples are in the first row. Check that you have all the samples you think you should have. In this example, it is 105
taxa = length(obs[,1]) #Same here, the taxa are in each column? 161 OTUs (excluding first row, which is your sample names)


# Initialise the abundance matrix
p=matrix(0,taxa,samples)

# Convert the data table to a matrix
for (k in 1:samples){
  p[,k]<-obs[,k]/sum(obs[,k])
}


###################################################
###2) dissimilarity matrix----
###################################################

#You need to provide ape with a .fasta file of sequences that have all been formatted to the same length (use D to fill in the sizes to the longest sequence length)
library(ape)
x <-read.dna("otusnhystrixreduced_samelength.fasta", "fasta")
#from ape manual: boot.phylo Tree Bipartition and Bootstrapping Phylogenies
#lots of diff models of evolution you can chose from, using the "raw" model here. Read ape documenation for more info

matrix<-dist.dna(x, model="raw", variance=TRUE, gamma=FALSE, pairwise.deletion=FALSE, base.freq=NULL, as.matrix=TRUE)
Z<-matrix #raw evolution model, diagnol is zero
#Z = species similarity information, comma delimited text file with entries between 0 and 1.


#################################################
#3)Similarity matrix-
#################################################
#ape pairwise spits out dissimilarity, needed to change it to similarity. I did this in excel but can be done in R

####Import new similarity matrix back into R
Z <- read.csv("Z similarity matrix.csv",header=T) #note: you will now replace the above Z object 

#, eqivalent ot exp(Shannon) if q=1 100 percent dissimilarity matrix, this should be equalivanet to not having a Z matrix because biologist used to treat species as if they had no similarities
#from leinster paper: naive measure of diversity is "where commanalities between species ignored. similarities between distinct species is ZERO." therefore dissimilarity is 1., all 100% dissimilar
#all zero. all species have zero similarity
#so along diagnols needs to be all ones. this is a similarity matrix. so each species is 100% similar to itself. however to have a naive diversity measure, "the similarites between distinct species is zero"
#below is naive similarity matrix


# Specify the x-axis on the diversity profile.  lenq should be at least 5 (probably larger) so the you can see the effect of the different weightings of abundant and rare species.
lenq = 50
qq <- seq(length=lenq, from=0, by=.11)



# Initialise the Zp matrix to zero
Zp=matrix(0,taxa,samples)

# Compute Zp
for (k in 1:samples) {
  for (i in 1:taxa){
    for (j in 1:taxa){
      Zp[i,k]<-Zp[i,k]+Z[i,j]*p[j,k]
    }
  }}



# Initialise the Diversity matrix to zero
Dqz = matrix(0, lenq ,samples)

#  Loop to calculate the diversity Dqz for each value of q (iq) and each sample (k)


for (k in 1:samples) {
  
  for (iq in 1:lenq)  {
    q<-qq[iq];
    
    
    for (zpi in 1:length(Zp[,k])){
      if (Zp[zpi,k]>0)(
        Dqz[iq,k]<-Dqz[iq,k]+ p[zpi,k]*(Zp[zpi,k])^(q-1))
    }
    
    Dqz[iq,k] <- Dqz[iq,k]^(1/(1-q));
  }
}

#THIS can take a while (up to 20 min) depending on OS and number of OUTs/samples

#Plot the diversity profiles for all the samples (e.g. understorey and canopy) on the same graph.
#Plot the columns of one matrix against the columns of another.
#in Dqz matrix, each column is a sample (I guess in the correct order, but check this)
#the rows are for each value of q used. So the first row is each sample at value q=0
matplot(qq,Dqz, type="l",
        xlab="Sensitivity parameter q",
        ylab="Diversity DqZ(p)",
        main="Diversity profile using a similarity matrix")
#or do a ggplot graph so you can add legend
head(Dqz)
head(qq)
#write.xlsx(Dqz, "Dqz.xlsx")
#write.xlsx(qq, "qq.xlsx")
dqz_qs<-read.csv("Dqz_qs.csv", header=T)
head(dqz_qs)
class(dqz_qs)
dqz_qs$q<-factor(dqz_qs$q)


#matplot is great because you can directly plot matrix columns, but it provides no default legend
library(reshape2)
#the diversity metric with all otus, normalized abundnace and APE raw similarity matrix
dataL<-melt(dqz_qs, id="q")                        
ggplot(dataL, aes_string(x="q", y="value", colour="variable", group="variable"))+geom_line()+labs(x ="Sensitivity parameter q", y="Diversity DqZ(p)", title="Diversity profile using a similarity matrix")+theme_bw()


#of those, I will only choose q=0 (where rare species are of most importance, so carry on to heritability index)
qzero<-(Dqz[1:2,])
class(qzero)

library(xlsx)
write.csv(qzero, "qzero.csv")

