##### Adapted from Dr Andrew Nguyen's workshop I found here https://adnguyen.github.io/2017_Ecological_Genomics/Tutorial/2017-03-22_pop4_followalong.html #####

#Packages used
library('dartR')
library('adegenet')
library('PopGenReport')
library('poppr')

#Load your genomic data (in my case a genlight object)
load("glBenmoref5.rdata")

#Perform a PCA - I used dartR to do this
pcBenmore <- gl.pcoa(glBenmoref5)

## Which SNPs load most strongly on the 1st PC axis?

loadingplot(abs(pcBenmore$loadings[,1]),
            threshold=quantile(abs(pcBenmore$loadings), 0.99))

## Get their locus names

threshold<-quantile(abs(pcBenmore$loadings[,1]),0.99)

loadings <- glBenmoref5$loc.names[which(abs(pcBenmore$loadings[,1])>threshold)] 

## To check that these were the large effect loci, I ran a PCA with only those loci, as well as a PCA without those loci to see the difference removing them made
#I used the gl.keep.loc() and gl.drop.loc() functions in dartR to do this

## Check if these loci have private alleles in your populations, which may explain why they are contributin to so much differentiation in your PCA

#Convert gl to genind and use poppr package
gi.pc1loadings <- gl2gi(gl.loadings)

load1pa <- private_alleles(gi.pc1loadings)

## Using popgenreport I looked at the allele frequencies of these loci of interest

pdf(file = "allele frequencies Benmore popgenreport PC1 loci.pdf")
allele.dist(gi.pc1loadings)
dev.off()

## Repeat the above steps but for the 2nd PC axis

loadingplot(abs(pcBenmore$loadings[,2]),
            threshold=quantile(abs(pcBenmore$loadings), 0.99))

threshold2<-quantile(abs(pcBenmore$loadings[,2]),0.99)

loadings2 <- glBenmoref5$loc.names[which(abs(pcBenmore$loadings[,2])>threshold)] #21 loci, 1 overlapping with PC1

gl.loadings2 <- gl.keep.loc(glBenmoref5, loc.list = loadings2) #gl of only loci with PC2 99th percentile loadings
