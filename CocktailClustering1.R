# R Code as published in Appendix S2
# Bruelheide, H. (2016): Cocktail clustering - a new hierarchical agglomerative 
# cluster algorithm for extracting species groups in vegetation databases. 
# Journal of Vegetation Science (submitted 22.10.2015, revised version submitted 
# 22.11.2015, decided for a minor revision 17.02.2016, JVS-RA-03393).

# Names of variables follow the definitions in 
# Bruelheide, H. 2000. A new measure of fidelity and its application to defining 
# species groups. Journal of Vegetation Science 11: 167-178.
# N = number of all plots (=relev?s) in the data set.
# n = number of all species in the data set.
# k = total number of species in a species group
# m = required minimum number a plot must have to belong to a species group

# The program returns several objects

# Cluster.species: 
# a binary cluster by species matrix that states which species occurs 
# in which cluster (i.e. in which species group)

# Cluster.info:
# a cluster by 2 matrix holding the information of k = cluster size k and m = the
# required minimum number of species in a cluster

# Plot.cluster:
# a binary cluster by plot matrix that states which plot is assigned to which cluster 
# according to the required criteria of that cluster (i.e. species group)

# Cluster.merged:
# a cluster by 2 matrix holding the information of which species or clusters are merged
# Row i of Cluster.merged describes the merging of clusters at step i of the clustering.
# A negative sign  is given if the element to be merged is a species.
# A positive value is given if the element to be merged is a cluster, formed at an 
# earlier stage of the clusterin galgorithm. 
# This notation correspondes to the "merge" value in the hclust procedure in R stats.
# Thus, negative entries indicate agglomerations of singletons, and positive entries 
# indicate agglomerations of non-singletons, i.e. clusters.

# Cluster.height: 
# A cluster x 1 matrix holding the phi value at which the cluster was formed.
# This notation correspondes to the "height" value in the hclust procedure in R stats.
 
### 1. Test data ###
library(vegan)
# Using Ter Braak's dune data set as the first test example
data(dune)
str(dune) # 20 plots x 30 species
vegmatrix <- dune
# vegmatrix holds the plots by species data

### 2. Calculate expected frequencies ###
vegmatrix[vegmatrix[,]>0] <- 1 # transform the plot x species matrix into a p/a matrix
N <- dim(vegmatrix)[[1]] # number of plots
n <- dim(vegmatrix)[[2]] # number of species !!!! later  be renamed to k
col.sums <- apply(vegmatrix, 2, sum) # frequency of every species in plots
row.sums <- apply(vegmatrix, 1, sum) # species richness of every plot
p.freq <- col.sums/N   # proportional frequency of occurrence of every species
q.freq <- 1-p.freq    # proportional frequency of non-occurrence of every species

### 3. Function Expected.plot.freq ###
# According to Bruelheide (2000):
# Define a function to calculated the expected frequencies of occurrences in plots, 
# comprising the classes of 0, 1, ... n species a plot might contain.
# Given the observed frequencies of species in the whole data set, 
# the numbers of relev?s are calculated that are expected to have
# >=0, >=1, ... >=k species (cumulative frequency distribution).
# The principle of the calculation of expected frequencies of plots for a species group 
# containing n species, resulting in 0, 1, ... k (i.e. k+1 classes) frequency classes,
# is to start with the first species, where the expected frequencies in class 0 and 1 
# are q.freq and p.freq.
# With the next species the expected frequency in class 0 is the previous value of 
# class 0 (Exp.no.of.plots.inter[1]) multiplied by q.freq of that species.
# Similarly, the expected frequency of the highest class is obtained by using the 
# previous value of that class (Exp.no.of.plots.inter[j]) multiplied by p.freq of 
# that species.
# The classes in between with 1, ... k-1 species are obtained by adding the product of
# the previous value of the class of 1 species lower (Exp.no.of.plots.inter[k]) 
# and p.freq and the product of the previous value of that class 
# (Exp.no.of.plots.inter[k+1]) and q.freq.

Expected.plot.freq <- function(species.in.cluster){
  No.of.spec.in.cluster <- length(species.in.cluster)
  Exp.no.of.plots <- array(0,No.of.spec.in.cluster+1, dimnames=list(seq(0,No.of.spec.in.cluster[1],1)))
  Exp.no.of.plots.inter <- array(0,No.of.spec.in.cluster+1,dimnames=list(seq(0,No.of.spec.in.cluster[1],1)))
  Exp.no.of.plots.inter[1] <- 1 # set to 1 in order to enable the loop run correctly
  for (j in 1:No.of.spec.in.cluster){
    Exp.no.of.plots[1] <- Exp.no.of.plots.inter[1]*q.freq[species.in.cluster[j]]
    # j = 1 refers to class >=0, i.e. the class with none or at least one species of
    # the species group, which is the total number of species n in the data set.
    if (j>1){
      for (k in 1:(j-1)){
        Exp.no.of.plots[k+1] <- Exp.no.of.plots.inter[k] * p.freq[species.in.cluster[j]] +
          Exp.no.of.plots.inter[k+1] * q.freq[species.in.cluster[j]] 
      }
    }    
    Exp.no.of.plots[j+1] <- Exp.no.of.plots.inter[j]*p.freq[species.in.cluster[j]]
    for (k in 1:(j+1)){
      Exp.no.of.plots.inter[k] <- Exp.no.of.plots[k]
    }
  }
  Exp.no.of.plots
}

### 4. Function Compare.obs.exp.freq ###
# This function compares expected and observed counts of how many plots occur 
# in the species frequency classes of 0, 1, .. k species of a species group and
# returns m, the minimum number of species a plot must have
# to belong to a species group.
# Addition of frequencies is done from frequency classes with highest to lowest
# species number.
Compare.obs.exp.freq <- function(Obs.freq, Exp.freq){
  No.of.spec.in.cluster <- length(Obs.freq) - 1
  Cum.obs.no.of.plots <- array(0,No.of.spec.in.cluster+1, 
                        dimnames=list(seq(0,No.of.spec.in.cluster[1],1)))
  Cum.exp.no.of.plots <- array(0,No.of.spec.in.cluster+1, 
                        dimnames=list(seq(0,No.of.spec.in.cluster[1],1)))
  # The frequency class of highest species number (i.e. all k species in the species group)
  # is set at the beginning
  Cum.obs.no.of.plots[No.of.spec.in.cluster+1] <- Obs.freq[No.of.spec.in.cluster+1]
  Cum.exp.no.of.plots[No.of.spec.in.cluster+1] <- Exp.freq[No.of.spec.in.cluster+1]
  m = 1
  # m is the minimum number of species any cluster should contain.
  # Setting this variable is done to produce the cluster at the highest hierarchical level,
  # which do not follow the species group criteria (as here species do not co-occur any
  # longer with higher probability than expected by chance). When that level is reached,
  # m is set to 0, and a species group is then defined by 1 out of k species.
  m.found <- -1
  # m.found is a variable initially set to -1.
  # After the first observed occurrences are found (i.e. Cum.obs.no.of.plots[No.of.spec.in.cluster]>0),
  # m.found is set to 0. This shows that m is now settable.
  # After m has been finally determined, m.found is set to 1.
  if (Cum.obs.no.of.plots[No.of.spec.in.cluster+1] > 
        Cum.exp.no.of.plots[No.of.spec.in.cluster+1]) {
    m <- No.of.spec.in.cluster
    m.found <- 0
    # This is usually the case for small clusters with two or three species 
    # which have occurrences in the class of k species out of k in the species group.
  }
  for (j in seq(No.of.spec.in.cluster,1,-1)){
    Cum.obs.no.of.plots[j] <- Cum.obs.no.of.plots[j+1]+Obs.freq[j]
    if (m.found==-1 & Cum.obs.no.of.plots[j] > 0){
       m.found <- 0
    }
    Cum.exp.no.of.plots[j] <- Cum.exp.no.of.plots[j+1]+Exp.freq[j]
    if (j>1 & m.found==0 & Cum.exp.no.of.plots[j] > Cum.obs.no.of.plots[j]){
        m = j
        m.found <- 1
    }
  }
  m
  # The minimum number of species a species group should contain is returned
}

### 5. Prepare arrays ###
Cluster.species <-  array(0,c(n-1,n)) 
# Assignment matrix, dim: cluster x species
# a 0/1 matrix showing which species belongs to which cluster (i.e. species group)
dimnames(Cluster.species)[[2]] <- names(vegmatrix)
# The dimnames in the first dimension are cluster numbers. They range from 1 to n-1
# and thus, are not required to be set here.

Cluster.info <-  array(0,c(n-1,2), dimnames=list(c(1:(n-1)),c("k","m"))) 
# k = cluster size, i.e. total number of species in that species group
# m = required minimum number of species a plot must have to belong to this cluster
# Matrix holding cluster information, dim: cluster x 2
Plot.cluster <-  array(0,c(N,n-1), dimnames=list(c(1:N),c(1:(n-1))))
# a binary matrix holding the information which plot belongs to which cluster, 
# dim: number of plots N x number of clusters n-1
# Plots are renumbered from 1 to N. No original plot names are kept.
# The dimnames in the second dimension are cluster numbers. They range from 1 to n-1

vegmatrix2 <- vegmatrix 
# a copy of the plot by species matrix is made that stepwise will be reduced  
# and transformed into a cluster by plot matrix
# The current dimnames are species and plot names. Species names will be replaced by
# cluster numbers. These are formed by pasting "c_" with the respective cluster numbers.
n2 <- n   
# a copy of number of species n in the data set is made that will be reduced stepwise
Cluster.merged <- array(0,c(n-1,2))
# a cluster by 2 matrix holding the information of which species or clusters are merged
Cluster.height <- array(0,n-1)
# A cluster x 1 matrix holding the phi value at which the cluster was formed.
i <- 0           # index for clusters to be defined, is increased in the while loop
multiple.max <- 1 
# multiple.max holds the information of how many maxima there are
# if there are multiple maxima, they are fused in subsequent steps without 
# recalculating the phi similarity

### 6. cluster loop ###
pb <- txtProgressBar(title = "progress bar", min = 0, max =n-1, width = 300)
while (i <= (n-2) ){
  # Calculate the whole phi distance matrix
  phi.index <- designdist(t(vegmatrix2), method = "(a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))",
                     terms = c("binary"), abcd = T, "phi")
  # The distance defined by phi.index gives NA if one of the factors
  # (a+c),(b+d),(a+b),(c+d) = 0,
  # which occurs if all plots have been assigned to one group
  phi.index[is.na(phi.index)]<- 0
  # Then NAs are replaced by 0
  allmax <- which(phi.index==max(phi.index))
  # All maxima are extracted. This is necessary to allow for simultaneous
  # clustering of all maxima
  elements.per.col <- array(0,n2-1)
  # used to extract phi values for certain pairs from the distance matrix
  for (x in 1:n2-1){
    elements.per.col[x] <- (x-1)*n2-((x-1)*x/2)
  }
  multiple.max <- length(allmax)
  e1 <- array(0,multiple.max)
  e2 <- array(0,multiple.max)
  for (j in 1:multiple.max){
    # find the elements to be combined in the new cluster
    e1[j] <- max(which(elements.per.col<allmax[j]))
    e2[j] <- allmax[j] - elements.per.col[e1[j] ] + e1[j]
  }
  
  # Test for sequential clustering of the same nodes
  compare1 <- as.vector(t(as.matrix(cbind(e1,e2))))
  if (anyDuplicated(compare1)>2) {
    # then more than two nodes are fused at the same time, which
    # causes circularity, unnecessary fusions have to be eleminated
    tobekept <- rep(TRUE,multiple.max)
    for (j in 2:multiple.max){
      compare2 <- compare1[1:(2*(j-1))]
      k1 <- length(match(compare2,e1[j])[!is.na(match(compare2,e1[j]))])
      k2 <- length(match(compare2,e2[j])[!is.na(match(compare2,e2[j]))])
      if (k1>0 & k2 >0){
        # then this node has to be excluded
        tobekept[j] <- FALSE
      }
    }
    e1 <- e1[tobekept]
    e2 <- e2[tobekept]
    multiple.max <- length(e1)
  }
  
  i1 <- i+1 # variable remembering the value of i at the start of this loop
  for (j in 1:multiple.max){
    # go through all elements to be combined in the new cluster
    i <- i + 1 # increment of loop variable
    Cluster.height[i] <- phi.index[allmax][j]
    if(substr(names(vegmatrix2)[e1[j]],1,2)=="c_"){
      # then the element merged into a cluster is a cluster itself
      cluster.no.1 <- as.numeric(strsplit(names(vegmatrix2)[e1[j]],split="_")[[1]][2])
      Cluster.merged[i,1] <- cluster.no.1
      # a positive sign is given if the element to be merged is a cluster
      Cluster.species[i,which(Cluster.species[cluster.no.1,]==1)] <- 1 
      # The species of the cluster are filled into the Cluster.species matrix 
      # which holds the information which species belongs to which cluster
    } else {
      f1 <- which(names(vegmatrix)==names(vegmatrix2)[e1[j]])
      Cluster.merged[i,1] <- - f1
      # a negative sign is given if the element to be merged is a species
      Cluster.species[i,f1] <- 1
      # The species of the cluster are filled into the Cluster.species matrix 
      # which holds the information which species belongs to which cluster
    }      
    if(substr(names(vegmatrix2)[e2[j]],1,2)=="c_"){
      # then the element merged into a cluster is a cluster itself
      cluster.no.2 <- as.numeric(strsplit(names(vegmatrix2)[e2[j]],split="_")[[1]][2])
      Cluster.merged[i,2] <- cluster.no.2
      # a positive sign is given if the element to be merged is a cluster
      Cluster.species[i,which(Cluster.species[cluster.no.2,]==1)] <- 1 
      # The species of the cluster are filled into the Cluster.species matrix 
      # which holds the information which species belongs to which cluster
    } else {
      f2 <- which(names(vegmatrix)==names(vegmatrix2)[e2[j]])
      Cluster.merged[i,2] <- - f2
      # a negative sign is given if the element to be merged is a species
      Cluster.species[i,f2] <- 1
      # The species of the cluster are filled into the Cluster.species matrix 
      # which holds the information which species belongs to which cluster
    }
    # now the information is kept that the species now belongs to a cluster
    names(vegmatrix2)[c(e1,e2)][names(vegmatrix2)[c(e1,e2)]==names(vegmatrix2)[e1[j]] ] <- 
        paste("c_",i, sep="")  
    names(vegmatrix2)[c(e1,e2)][names(vegmatrix2)[c(e1,e2)]==names(vegmatrix2)[e2[j]] ] <- 
        paste("c_",i, sep="") 
    # Species names in the dimnames of vegmatrix2 are replaced by cluster numbers. 
    # These are formed by pasting "c_" with the respective cluster numbers.
    
    Cluster.info[i,1] <- sum(Cluster.species[i,])
    # The size of the spcies group is assigned to Cluster.info
    
    # Calculate "Observed.plot.freq"
    species.in.plot <- rowSums(vegmatrix[,which(Cluster.species[i,]>0)])  
    species.in.plot.table <- as.matrix(table(species.in.plot))
    Obs.plot.freq <- matrix(0, Cluster.info[i,1]+1, dimnames=list(0:Cluster.info[i,1]))
    index <- match(dimnames(Obs.plot.freq)[[1]],dimnames(table(species.in.plot))$species.in.plot)
    Obs.plot.freq[,1] <- as.numeric(table(species.in.plot))[index]
    Obs.plot.freq <- ifelse(is.na(Obs.plot.freq),0,Obs.plot.freq)
    # the results are absolute number of plots per frequency class
    # Call the function "Expected.plot.freq", defined above
    Exp.plot.freq <- Expected.plot.freq(which(Cluster.species[i,]>0)) * N
    # proportional frequencies are multiplied by total plot number N to 
    # obtain  expected number of plots per frequency class
    
    # call function "Compare.obs.exp.freq", defined above
    # The result is the minimum number of species required for a plot to be assigned
    # to the species group
    Cluster.info[i,2] <- Compare.obs.exp.freq(Obs.plot.freq, Exp.plot.freq)
    # The minimum number m is assigned to the second column of Cluster.info

    # Calculate which plots belong to that cluster, i.e. which transgress the 
    # required minimum number m 
    Plot.cluster[species.in.plot >= Cluster.info[i,2],i]  <- 1
    
  }  
  if (n2==3) {
    # required to remember the column names when only three cluster remain
    name.last.cluster <- names(vegmatrix2)[-c(e1[1],e2[1])]
  }
  # only those clusters have to be added that have not already been fused
  # at the same level of phi.index
  index.e1 <- !duplicated(names(vegmatrix2)[e1],fromLast=T)
  index.e2 <- !duplicated(names(vegmatrix2)[e2],fromLast=T)
  #index.e1 <- !duplicated(c(e1,e2),fromLast=T)[1:multiple.max]
  #index.e2 <- !duplicated(c(e1,e2),fromLast=T)[(multiple.max+1):(2*multiple.max)]
  index.e <- c(i1:i)[index.e1 & index.e2]
  # those clusters with index.e==F will be omitted
  vegmatrix2 <- cbind(vegmatrix2[,-unique(c(e1,e2))],Plot.cluster[,index.e])
  n2 <- dim(vegmatrix2)[[2]]
  for (j in 1:length(index.e)){
    # those clusters with index.e==F will be omitted
    names(vegmatrix2)[n2-length(index.e)+j] <- paste("c_",index.e[j], sep="") 
  }
  
  if (n2==2) {
    # required to remember the column names when only two cluster remain
    names(vegmatrix2)[1] <- name.last.cluster
  }  
      
  if (sum(Plot.cluster[,i])==N){   
    # Stop criterion, reached when all plots are merged
    # Then all clusters are combined at phi = 0
    for(j in (i+1):(n-1)){
      g1 <- dim(vegmatrix2)[[2]] # the latest cluster is the rightmost one in vegmatrix2
      cluster.no.1 <- as.numeric(strsplit(names(vegmatrix2)[g1],split="_")[[1]][2])
      Cluster.merged[j,1] <- cluster.no.1
      g2 <- 1
      # the leftmost cluster is combined first
      cluster.no.2 <- as.numeric(strsplit(names(vegmatrix2)[g2],split="_")[[1]][2])
      Cluster.merged[j,2] <- cluster.no.2
      Cluster.height[j]=0
      Plot.cluster[,j] <- 1
      Cluster.info[j,1] <- sum(Cluster.species[j,])
      Cluster.info[j,2] <- 1
      Cluster.species[j,which(Cluster.species[cluster.no.1,]==1)] <- 1 
      Cluster.species[j,which(Cluster.species[cluster.no.2,]==1)] <- 1 
      vegmatrix2 <- cbind(vegmatrix2,Plot.cluster[,j])
      names(vegmatrix2)[dim(vegmatrix2)[[2]]] <- paste("c_",j, sep="")
      vegmatrix2 <- vegmatrix2[,-c(g1,g2)]
      n2 <- dim(vegmatrix2)[[2]]
    }
    i <- j
  }
  setTxtProgressBar(pb, i, title=paste("Cluster ",i,"of ", n-1, " is processed"))
}
close(pb)

# write output files for subsequent use
write.csv(Cluster.species,"Cluster_species_dune.csv",row.names=F)
write.csv(Plot.cluster,"Plot_cluster_dune.csv",row.names=F)

### 7. Plotting the tree of the Cocktail clustering results ###
Species.sort <- array("",n)
# holds the information in which sequence the species are arranged along the tips
# of the Tree. Species.sort holds a string of binary numbers for all n-1 clusters
# that show the assignment to a cluster in ascending order, i.e. from clusters formed last 
# to those formed first
for (i in 1:n){
  for (j in 1:n-1){
    Species.sort[i] <- paste(Species.sort[i],Cluster.species[order(Cluster.height)[j],i],sep="")
  }
}
plot(c(1:(n)),c(seq(0.1,-1,(-0.1-1)/(n-1))),type="n", xaxt = "n", yaxt = "n", xlab="", ylab=expression(paste(phi," coefficient")))
axis(1, las=2,at=seq(1:n), labels=names(vegmatrix)[order(Species.sort)], cex.axis=1)
axis(2, las=2,at=seq(0,-1,-0.2), labels=seq(0,1,0.2))  
Cluster.position <- array(NA, c(n-1,6),dimnames=list(c(1:(n-1)), c("left.leg.x", "left.leg.y0", "left.leg.y1","right.leg.x", "right.leg.y0","right.leg.y1")))
# Array with x and y coordinates for all clusters. Coordinates have a negative
# sign to arrange the tree in descending values of the phi.index
# y0: lower part of the leg
# y1: upper part of the leg
for (i in 1:(n-1)){
   # loop for all clusters
   if (Cluster.merged[i,1] < 0){
     # left leg
     Cluster.position[i,1] <- which(order(Species.sort)==-Cluster.merged[i,1])
     Cluster.position[i,2] <- -1
     Cluster.position[i,3] <- -Cluster.height[i]
   } else {
     Cluster.position[i,1] <- 
       mean(match(which(Cluster.species[Cluster.merged[i,1],]==1),order(Species.sort)))
     Cluster.position[i,2] <- -Cluster.height[Cluster.merged[i,1]]
     Cluster.position[i,3] <- -Cluster.height[i]
   }
   if (Cluster.merged[i,2] < 0){
     # right leg
     Cluster.position[i,4] <- which(order(Species.sort)==-Cluster.merged[i,2])
     Cluster.position[i,5] <- -1
     Cluster.position[i,6] <- -Cluster.height[i]
   } else {
     Cluster.position[i,4] <- 
       mean(match(which(Cluster.species[Cluster.merged[i,2],]==1),order(Species.sort)))
     Cluster.position[i,5] <- -Cluster.height[Cluster.merged[i,2]]
     Cluster.position[i,6] <- -Cluster.height[i]
   }
   # labels
   text(mean(c(Cluster.position[i,1],Cluster.position[i,4])),
        Cluster.position[i,6]-0.03,i, cex=0.6)   
   # left leg
   lines(c(Cluster.position[i,1],Cluster.position[i,1]),
         c(Cluster.position[i,2],Cluster.position[i,3]))
   # right leg
   lines(c(Cluster.position[i,4],Cluster.position[i,4]),
         c(Cluster.position[i,5],Cluster.position[i,6]))
   # horizontal line
   lines(c(Cluster.position[i,1],Cluster.position[i,4]),
         c(Cluster.position[i,6],Cluster.position[i,6]))
}

### 8. Plotting UPGMA clusters using the phi coefficient ###
phi.index <- designdist(t(dune), method = "(a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))",
                        terms = c("binary"), abcd = T, "phi")
phi.complement <- as.dist(1-phi.index) 
# the complement to phi is used, as a distance matrix instead of a similarity
# matrix is needed
dune.clust.average <- hclust(phi.complement, "average")
plot(dune.clust.average)

# Cluster diagram drawn manually for a direct comparison to Cocktail clustering
A2 <-  array(0,c(n-1,n)) 
# Assignment matrix, dim: cluster x species
for (i in 1:(n-1)){
  if (dune.clust.average$merge[i,1] < 0){
    A2[i,-dune.clust.average$merge[i,1]] <- 1 
  } else {
    A2[i,] <- A2[i,]+A2[dune.clust.average$merge[i,1],]
  }
  if (dune.clust.average$merge[i,2] < 0){
    A2[i,-dune.clust.average$merge[i,2]] <- 1 
  } else {
    A2[i,] <- A2[i,]+A2[dune.clust.average$merge[i,2],]
  }
}
Species.sort <- array("",n)
for (i in 1:n){
  Species.sort[i]  <- paste(A2[order(1-dune.clust.average$height),i],collapse="")
}
plot(c(1:(n)),c(seq(0.2,-1,(-0.2-1)/(n-1))),type="n", xaxt = "n", yaxt = "n", xlab="", ylab=expression(paste(phi," coefficient")))
axis(1, las=2,at=seq(1:n), labels=names(vegmatrix)[order(Species.sort)], cex.axis=1)
axis(2, las=2,at=seq(0,-1,-0.2), labels=seq(0,1,0.2))  
Cluster.position <- array(NA, c(n-1,6),dimnames=list(c(1:(n-1)), c("left.leg.x", "left.leg.y0", "left.leg.y1","right.leg.x", "right.leg.y0","right.leg.y1")))
for (i in 1:(n-1)){
  if (dune.clust.average$merge[i,1] < 0){
    # left leg
    Cluster.position[i,1] <- which(order(Species.sort)==-dune.clust.average$merge[i,1])
    Cluster.position[i,2] <- -1
    Cluster.position[i,3] <- -(1-dune.clust.average$height[i])
  } else {
    Cluster.position[i,1] <- 
      mean(match(which(A2[dune.clust.average$merge[i,1],]==1),order(Species.sort)))
    Cluster.position[i,2] <- -(1-dune.clust.average$height)[dune.clust.average$merge[i,1]]
    Cluster.position[i,3] <- -(1-dune.clust.average$height[i])
  }
  if (dune.clust.average$merge[i,2] < 0){
    # right leg
    Cluster.position[i,4] <- which(order(Species.sort)==-dune.clust.average$merge[i,2])
    Cluster.position[i,5] <- -1
    Cluster.position[i,6] <- -(1-dune.clust.average$height[i])
  } else {
    Cluster.position[i,4] <- 
      mean(match(which(A2[dune.clust.average$merge[i,2],]==1),order(Species.sort)))
    Cluster.position[i,5] <- -(1-dune.clust.average$height)[dune.clust.average$merge[i,2]]
    Cluster.position[i,6] <- -(1-dune.clust.average$height[i])
  }
  # labels
  text(mean(c(Cluster.position[i,1],Cluster.position[i,4])),
       Cluster.position[i,6]-0.03,i, cex=0.6)   
  # left leg
  lines(c(Cluster.position[i,1],Cluster.position[i,1]),
        c(Cluster.position[i,2],Cluster.position[i,3]))
  # right leg
  lines(c(Cluster.position[i,4],Cluster.position[i,4]),
        c(Cluster.position[i,5],Cluster.position[i,6]))
  # horizontal line
  lines(c(Cluster.position[i,1],Cluster.position[i,4]),
        c(Cluster.position[i,6],Cluster.position[i,6]))
}

