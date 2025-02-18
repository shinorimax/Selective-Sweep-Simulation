### Compute pairwise distances 
library("devtools")
install_github("JuliaPalacios/phylodyn",force=TRUE)
library("phylodyn")
library("ape")

trees1<-read.tree("~/Documents/Introgression/data/Ioutput(3).txt")
seq_length<-read.delim("~/Documents/Introgression/data/Itree_lengths(3).csv",sep=",")
#true_reg<-read.delim("~/Documents/Introgression/data/True_introgressed_regions(3).csv",header=TRUE,sep=",",nrows=3)

start<-as.numeric(gsub(",","",true_reg[,2]))
end<-as.numeric(gsub(",","",true_reg[,3]))
seq_length[,3]<-cumsum(seq_length[,2])


start<-c(167674,422107,858128,945109)
end<-c(198267,644020,865878,1000000)
#pairwise differences
n_tr<-length(trees1)
d2.dmat <- matrix(0, nrow=n_tr, ncol=n_tr)
for (i in 1:n_tr) {
  for (j in i:n_tr) {
    d2.dmat[i,j] <- dist_pairwise(trees1[[i]],trees1[[j]],dist.method="l2",weighted=TRUE,tol=3)
  }
}

d2.dmat_unw <- matrix(0, nrow=n_tr, ncol=n_tr)
for (i in 1:n_tr) {
  for (j in i:n_tr) {
    d2.dmat_unw[i,j] <- dist_pairwise(trees1[[i]],trees1[[j]],dist.method="l2",weighted=FALSE,tol=3)
  }
}

saveRDS(list(unw=d2.dmat_unw,ww=d2.dmat),file="distances_intro.rds")
d2.dmatsimUNW<-d2.dmat_unw+t(d2.dmat_unw)
d2.dmatsim<-d2.dmat+t(d2.dmat)
total_original<-apply(d2.dmatsim,1,sum)
total_original_unw<-apply(d2.dmatsimUNW,1,sum)

par(mfrow=c(1,2))
plot(seq_length[,3],total_original,main="Introgressed",ylab="distance to all",xlab="Genomic position")
plot(seq_length[,3],total_original_unw)


for (j in 1:length(start)){
  abline(v=start[j],col="red")
  abline(v=end[j],col="red")
}


for (j in 1:nrow(true_reg)){
  abline(v=start[j]/1000,col="red")
  abline(v=end[j]/1000,col="blue")
}

##compute total lenght of trees
treelength<-rep(0,length(trees1))
firstcoal<-rep(0,length(trees1))
for (j in 1:length(trees1)){
  phy <- summarize_phylo(trees1[[j]])
  summ<-phylodyn:::gen_INLA_args(samp_times = phy$samp_times, n_sampled = phy$n_sampled,
                coal_times = phy$coal_times)
  treelength[j]<-sum(summ$lineages*diff(summ$s))
  firstcoal[j]<-phy$coal_times[1]
}

plot(seq_length[,3],treelength)
cor(treelength,total_original)
#points(seq_length[,3],firstcoal,col="red")

plot(seq_length[,3],firstcoal,col="red")

plot(total_original,treelength)
abline(c(0,1))

##This code is to detect coalescence of 2 and 3 with coal time < 15000
##Do not use this, tip labels are weird
# n<-trees1[[1]]$Nnode+1
# list_max<-0
# coal_times<-0
# 
# for (j in 1:n_tr){
# off<-trees1[[j]]$edge[,1]
# ##want to find true 2 and 3
# l1<-seq(1,n)[trees1[[j]]$tip.label==2]
# l2<-seq(1,n)[trees1[[j]]$tip.label==3]
# 
# if (sum(sort(trees1[[j]]$edge[off==5,2])==c(1,2))==2) {
#   print("it is a pair")
#   print(j)
#   #it is a coalescence of 1 and 2
#   if (coalescent.intervals(trees1[[j]])$interval.length[1]<15000){
#     print("potential")
#     print(j)
#     list_max<-c(list_max,j)
#     coal_times<-c(coal_times,coalescent.intervals(trees1[[j]])$interval.length[1])
#   }
# }
# }
# 
# list_max<-list_max[-1]
# coal_times<-coal_times[-1]
# #what is the lenght?

#plot(seq_length[list_max,3],coal_times)


            

#Now we will discretize the genome to make the histogram comparisons "fair"
##we will pick the length of the grid at which we discretize the genome
ngrid<-1000
lengthgenome<-sum(seq_length[,2])

lgrid<-1000 #length of the grid = lengthgenome/ngrid
aux_bins<-diff(seq_length[,3]%/%lgrid)
ind_bins<- which(aux_bins>0) ## selects the trees
trees_bin <- c(1,rep(ind_bins ,aux_bins[which(aux_bins>0)])) #adds the first tree and 
#repeats the number of times a single tree occupies different trees
newd2 <- d2.dmat[trees_bin,trees_bin]
newd2unw <- d2.dmat_unw[trees_bin,trees_bin]


#Let's look at the distribution of distances
#list.d.all<-as.vector(newd2[upper.tri(d2.dmat, diag = FALSE)]) 

list.d<- as.vector(newd2[upper.tri(newd2, diag = FALSE)])
list.dunw<- as.vector(newd2unw[upper.tri(newd2unw, diag = FALSE)])
#list.d2<- as.vector(newd2[upper.tri(newd2.mat, diag = FALSE)]) 

apply(list.d,1,sum)
p1<-hist(list.d,main="With introgression",freq=FALSE)
p1unw<-hist(list.dunw,main="With introgression",freq=FALSE)

par(mfrow=c(3,1)list.dpar(mfrow=c(3,1), mar=c(0,0,3,3))
hist(list.d.all,main="All trees")
hist(list.d,main="Bins")
hist(list.d2,main="Bins")

#MDS plot of binned weighted distance introgressed
newd2sim<-newd2+t(newd2)
mds <- cmdscale(newd2sim)
plot(mds[,1],mds[,2],main="MDS Weighted")
for (j in 1:length(start)){
points(mds[(start[j]%/%lgrid):(end[j]%/%lgrid),1],mds[(start[j]%/%lgrid):(end[j]%/%lgrid),2],col="red",pch=16)
}
totals<-apply(newd2sim,1,sum)


hist(totals)
plot(totals, main="with introgression")
indic<-rep(mean(totals),nrow(newd2))
indic[list_max]<-max(totals)
points(indic,col="red")
hist(totals[list_max])
plot(sort(totals))

##Without introgression


treesB<-read.tree("~/Documents/Introgression/data/output(3).txt")
seq_lengthB<-read.delim("~/Documents/Introgression/data/tree_lengths(3).csv",sep=",")
#true_reg<-read.delim("~/Documents/Introgression/data/True_introgressed_regions(3).csv",header=TRUE,sep=",",nrows=3)

seq_lengthB[,3]<-cumsum(seq_lengthB[,2])


#pairwise differences
n_tr<-length(treesB)
d2.dmatB <- matrix(0, nrow=n_tr, ncol=n_tr)
for (i in 1:n_tr) {
  for (j in i:n_tr) {
    d2.dmatB[i,j] <- dist_pairwise(treesB[[i]],treesB[[j]],dist.method="l2",weighted=TRUE,tol=3)
  }
}

d2.dmat_unwB <- matrix(0, nrow=n_tr, ncol=n_tr)
for (i in 1:n_tr) {
  for (j in i:n_tr) {
    d2.dmat_unwB[i,j] <- dist_pairwise(treesB[[i]],treesB[[j]],dist.method="l2",weighted=FALSE,tol=3)
  }
}

saveRDS(list(unw=d2.dmat_unwB,ww=d2.dmatB),file="distances_introB.rds")
d2.dmatsimUNWB<-d2.dmat_unwB+t(d2.dmat_unwB)
d2.dmatsimB<-d2.dmatB+t(d2.dmatB)
total_originalB<-apply(d2.dmatsimB,1,sum)
total_original_unwB<-apply(d2.dmatsimUNWB,1,sum)

par(mfrow=c(1,2))
plot(seq_lengthB[,3],total_originalB,main="Standard",ylab="distance to all",xlab="Genomic position")
plot(seq_lengthB[,3],total_original_unwB)



##compute total lenght of trees
treelengthB<-rep(0,length(treesB))
firstcoalB<-rep(0,length(treesB))
for (j in 1:length(treesB)){
  phy <- summarize_phylo(treesB[[j]])
  summ<-phylodyn:::gen_INLA_args(samp_times = phy$samp_times, n_sampled = phy$n_sampled,
                                 coal_times = phy$coal_times)
  treelengthB[j]<-sum(summ$lineages*diff(summ$s))
  firstcoalB[j]<-phy$coal_times[1]
}

plot(seq_lengthB[,3],treelengthB)
cor(treelengthB,total_originalB)
#points(seq_length[,3],firstcoal,col="red")

#plot(seq_lengthB[,3],firstcoalB,col="red")

plot(total_originalB,treelengthB)
abline(c(0,1))


ngrid<-1000
lengthgenomeB<-sum(seq_lengthB[,2])

lgrid<-1000 #length of the grid = lengthgenome/ngrid
aux_bins<-diff(seq_lengthB[,3]%/%lgrid)
ind_bins<- which(aux_bins>0) ## selects the trees
trees_binB <- c(1,rep(ind_bins ,aux_bins[which(aux_bins>0)])) #adds the first tree and 
#repeats the number of times a single tree occupies different trees
newd2B <- d2.dmatB[trees_binB,trees_binB]
newd2Bunw <- d2.dmat_unwB[trees_binB,trees_binB]


#Let's look at the distribution of distances
#list.d.all<-as.vector(newd2[upper.tri(d2.dmat, diag = FALSE)]) 

list.dB<- as.vector(newd2B[upper.tri(newd2B, diag = FALSE)])
#list.d2<- as.vector(newd2[upper.tri(newd2.mat, diag = FALSE)]) 
list.dBunw<- as.vector(newd2Bunw[upper.tri(newd2Bunw, diag = FALSE)])

p2<-hist(list.dB,main="Without introgression",freq=FALSE)
p2unw<-hist(list.dBunw,main="Without introgression",freq=FALSE)

plot( p2, col=rgb(0,0,1,1/4),freq=FALSE,ylim=c(0,.00001),main="Pairwise tree distance",xlab="Distance")  # first histogram
plot( p1, col=rgb(1,0,0,1/4), add=T,freq=FALSE)
legend('topright',c('standard','introgressed'),
       fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), bty = 'n',
       border = NA)

plot( p2unw, col=rgb(0,0,1,1/4),freq=FALSE,ylim=c(0,.5),main="Pairwise tree distance",xlab="Distance")  # first histogram
plot( p1unw, col=rgb(1,0,0,1/4), add=T,freq=FALSE)
legend('topright',c('standard','introgressed'),
       fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), bty = 'n',
       border = NA)
#How about unweighted?

