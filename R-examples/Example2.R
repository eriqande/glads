#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#
#                      Examples for submitted manuscript:
# The genomic landscape of divergence across the speciation continuum in an island-colonising bird
#
# Questions to:
# Claudio S. Quilodrán. Department of Zoology, University of Oxford.
#         Email: claudio.quilodran@zoo.ox.ac.uk; claudio.quilodran@unige.ch
# Eric C. Anderson. Department of Ecology and Evolutionary Biology, University of California.
#         Email: eric.anderson@noaa.gov
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
rm(list=ls())
graphics.off()
library(gids)
library(progress)


##################
### Functions  ###
##################
# mating function
resample <- function(x, ...) x[sample.int(length(x), ...)]
mating <- function(fitn, struct){
  females      <- which(fitn[,"sex"] %in% 1)
  males        <- which(fitn[,"sex"] %in% 2)
  try(if(length(females)<1 || length(males)<1 ) stop("Population extinction. There is not enough females or males for mating.", call.=F))
  female.fitn  <- fitn[females,]
  female.stru  <- struct[females,,, drop = F]
  male.fitn    <- fitn[males,]
  mum.stru     <- female.stru[c(rep(1:length(females),times=female.fitn[,"fit"])),,, drop = F]
  dads         <- resample(males,sum(female.fitn[,"fit"]),TRUE,prob=male.fitn[,"fit"])
  dad.stru     <- struct[c(dads),,, drop = F]
  return(list(mum.stru,dad.stru))
}


######
doer <- function(x){
  struct <- x[[1]]
  sex.ratio <- x[[2]]
  mean.fitness <- x[[3]]
  loci.pos <-  x[[4]]
  chromo_mb <- x[[5]]
  ve <- x[[6]]
  crossover <- x[[7]]

  z <- as.data.frame(cbind(sex = rbinom(nrow(struct), 1, sex.ratio)+1)) #zs
  rownames(z)=rownames(struct)

  n = nrow(struct)
  z$fit <- round( pmax(rpois(n, lambda=mean.fitness) - ve*n, 0) )

  pairs <- mating(z,struct)

  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))


  struct.rt[,,1] <- rcpp_recombo_segregate_expo(mums, dim(mums), loci.pos, chromo_mb, crossover)
  struct.rt[,,2] <- rcpp_recombo_segregate_expo(dads, dim(dads), loci.pos, chromo_mb, crossover)

  return(struct.rt)
}


###
##########################################
#' disperal function
#' @param pop1 pop struct 1
#' @param pop2 pop struct 2
#' @param rate1to2  rate at which indivs in pop1 migrate to pop2
#' @param rate2tp1 rate at which indivds in pop2 migrate to pop1
#' @return Returns a list of pop structs
#' @export
fast_dispersal <- function(pop1,pop2,rate1to2,rate2to1){
  n1 <- dim(pop1)[1]
  n2 <- dim(pop2)[1]
  L <- dim(pop1)[2]
  g <- 2 #diploidia?
  p1 <- ifelse(runif(n1)<rate1to2,2,1)
  p2 <- ifelse(runif(n2)<rate2to1,1,2)

  ret <- rcpp_dispersal_placement(pop1, pop2, dim(pop1), dim(pop2), p1, p2);
  dim(ret[[1]]) <- c(sum(c(p1, p2)==1), L, g)
  dim(ret[[2]]) <- c(sum(c(p1, p2)==2), L, g)

  ret
}


######
main.function<- function(start.1, start.2, sex.ratio, mean.fitness, time, loci.pos, disp, chromo_mb, ve, crossover){


  n.gens = time

  struct.1 <- start.1
  struct.2 <- start.2

  sex.ratio <- sex.ratio
  mean.fitness <- mean.fitness

  loci.pos <-  loci.pos
  chromo_mb <- chromo_mb

  disp <- disp
  ve <- ve

  crossover <- crossover

  ############################################
  pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)

  res <- list()
  for (i in 1:n.gens){

    x1 <- list(struct.1, sex.ratio, mean.fitness, loci.pos, chromo_mb, ve, crossover)
    x2 <- list(struct.2, sex.ratio, mean.fitness, loci.pos, chromo_mb, ve, crossover)


    x <- list(x1,x2)

    out <- lapply(x,doer)
    struct.1 <- out[[1]]
    struct.2 <- out[[2]]

    outd <- fast_dispersal(struct.1,struct.2,disp,disp)
    struct.1 <- outd[[1]]
    struct.2 <- outd[[2]]
    if (i %% n.gens==0) res<- list(struct.1,struct.2)
    pb$tick()
  }

  return(res)
}


#########################
### Orginizing datset ###
#########################
Geno<-read.table("./R-examples/Data/165-129K_ZF_ified_17.geno", header=F, comment.char = "")

Geno.t<-t(Geno[,-1])
Geno.t[1,1] <- "ID"
colnames(Geno.t) <- Geno.t[1,]
Geno.t<-Geno.t[-1,]
rownames(Geno.t) <- Geno.t[,1]
Geno.t<-Geno.t[,-1]


n1<-apply(Geno.t,1:2,function(i) strsplit(i,split='')[[1]][1])
n2<-apply(Geno.t,1:2,function(i) strsplit(i,split='')[[1]][2])

n1[n1 %in% "A"] <- 1
n1[n1 %in% "C"] <- 2
n1[n1 %in% "G"] <- 3
n1[n1 %in% "T"] <- 4
n1[n1 %in% "N"] <- NA

n1<-apply(n1,2, as.numeric)
rownames(n1)<-rownames(n2)


n2[n2 %in% "A"] <- 1
n2[n2 %in% "C"] <- 2
n2[n2 %in% "G"] <- 3
n2[n2 %in% "T"] <- 4
n2[n2 %in% "N"] <- NA

n2<-apply(n2,2, as.numeric)
rownames(n2)<-rownames(n1)

struct.t0<-list(n1,n2)

L<-list(struct.t0[[1]], struct.t0[[2]])
struct.t00<-array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)), dimnames=list(rownames(L[[1]]),colnames(L[[1]]), 1:2 ))

#########################
###  Parameter values ###
#########################

loci.pos<-as.numeric(colnames(struct.t00) )
loci.pos<-loci.pos-loci.pos[1] + 1
chromo_mb<-max(loci.pos)

crossover <- 3/100000000.0 #recombination rate (cM/MB)
sex.ratio = 0.5 			#ratio of males and females
Fledged = 1.9				#number of fledglings by clutch
Survival.summer = 0.96 		#Survival rate at summer
Survival.autoumn = 0.76 		#Survival rate at autoumn
Survival.winter = 0.63 		#Survival rate at winter

Generation.time = 3 			#generation time (years)
mean.fitness = Fledged*Survival.summer*Survival.autoumn*Survival.winter*Generation.time

n.gens = 100	#number of simulated generations
ve= 0.0013		#density-dependent demographic effect
ini.ind = 400 	#initial number of individuals at each population
disp=0.001 		#migration rate

############################################
##### Generating initial individuals #######
############################################

dat0<-as.data.frame(rbind(struct.t0[[1]],struct.t0[[2]]) )
rownames(dat0)<-1:nrow(dat0)
snp<-lapply(1:ncol(dat0), function(i){ table(dat0[[i]], exclude=NA) })

set.seed(1)
Pop<-sapply(1:ini.ind, function(i){
	print(i)
		sapply(1:ncol(struct.t0[[1]]), function(j){
			n1<-sample(names(snp[[j]]),1,T, prob=snp[[j]])
			n2<-sample(names(snp[[j]]),1,T, prob=snp[[j]])
	return(list(n1, n2))
			})
	 })


Pop<-t(Pop)
even<-seq_len(ncol(Pop)) %% 2 ==0

n1<-Pop[, even]
n2<-Pop[, !even]

n1<-apply(n1,2, as.numeric)
n2<-apply(n2,2, as.numeric)

L<-list(n1,n2)

start<-array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)), dimnames=list(1:ini.ind,colnames(dat0), 1:2 ))

start.1=start #Initial population 1
start.2=start #Initial population 2

#########################
###    Simulations    ###
#########################

out<-main.function(start.1, start.2, sex.ratio, mean.fitness, n.gens, loci.pos, disp, chromo_mb, ve, crossover)


###################################################
###               Analyzing Fst         ###########
###               Python Script         ###########
###        modified from evomics.org        #######
# https://github.com/simonhmartin/genomics_general#
###################################################

###Orginizing dataset
#First population
P1<-out[[1]]
P1<-sapply(1:dim(P1)[2], function(i){
	paste(P1[,i,1], P1[,i,2], sep="")
})
colnames(P1)<-colnames(start)
P1 <-gsub("NA", "N", P1)
P1 <-chartr("1234", "ACGT ", P1)
P1 <-t(P1)
P1<-as.data.frame( cbind(rep("Chr_17", nrow(P1)), colnames(start), P1) )
colnames(P1)<-c("#CHROM", "POS", paste("P1_", 1:(ncol(P1)-2), sep="") )
rownames(P1)<-1:dim(P1)[1]

#Second population
P2<-out[[2]]
P2<-sapply(1:dim(P2)[2], function(i){
	paste(P2[,i,1], P2[,i,2], sep="")
})
colnames(P2)<-colnames(start)
P2 <-gsub("NA", "N", P2)
P2 <-chartr("1234", "ACGT ", P2)
P2 <-t(P2)
colnames(P2)<-paste("P2_", 1:ncol(P2), sep="")

#both populations in a single file
Data_P1_P2<-as.data.frame(cbind(P1,  P2 ) )
filename=paste("Output.geno", sep="")
write.table(Data_P1_P2, filename, sep="\t", row.name=F, quote=F)

#Calling the Python script
pop1<-dim(out[[1]])[1]
pop2<-dim(out[[2]])[1]
filename2<-paste(filename, ".gz", sep="")

	system(paste("gzip", filename))
	system( paste("python -W ignore ./R-examples/Data/calculate_FST.py", filename2, pop1+1, pop2 +1) )

	filename3=paste("FST_",  filename2, ".csv", sep="")
	Fst<-read.csv(filename3 , na.strings = "nan")[,6]

	system( paste("rm ", filename2) )
	system( paste("rm ", filename3) )
	system("rm ./R-examples/Data/genomics.pyc")


#########################
#####    Figure    ######
#########################

h <- hist(Fst, breaks=seq(0, 1, length=15), freq=T)
h$counts=h$counts/sum(h$counts)
plot(h, col="grey", xlim=c(0,1),  ylim=c(0, 1), xlab="Fst", ylab="Proportion", main=NULL, bty="l", xaxs="i",yaxs="i", las=1, border=F)
