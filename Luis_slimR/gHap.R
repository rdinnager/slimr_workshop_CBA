#Function: ghap.loadphase
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load phased genotypes
  #  samples.file <-  samples_haplo
  # markers.file <- markers_haplo
  # phase.file <-phase_haplo
  # verbose = TRUE

ghap.loadphase<-function(
  samples.file,
  markers.file,
  phase.file,
  verbose = TRUE
){
  
  #Check files
  # if(file.exists(phase.file) == FALSE){
  #   stop("Could not find phased genotypes file")
  # }
  # if(file.exists(markers.file) == FALSE){
  #   stop("Could not find marker map file")
  # }
  # if(file.exists(samples.file) == FALSE){
  #   stop("Could not find sample file")
  # }
  # 
  #Load marker map file
  if(verbose == TRUE){
    cat("\nReading in marker map information... ")
  }
  marker<- markers.file
    # read.table(markers.file,header=FALSE,colClasses = c("character","character","numeric","character","character"))
  
  #Check if the map file contains correct dimension
  if(ncol(marker) != 5){
    stop("Marker map contains wrong number of columns (expected 3)")
  }
  
  #Check if the file contains information on a single chromosome
  chr<-unique(marker[,1])
  if(length(chr) != 1){
    stop("Your marker map file contains information on more than one chromosome")
  }
  
  #Check for duplicated marker ids
  if(length(unique(marker[,2])) < nrow(marker)){
    stop("Your marker map file contains duplicated ids!")
  }
  
  #Check if markers are sorted by bp
  if(identical(marker[,3],sort(marker[,3])) == FALSE){
    stop("Markers are not sorted by base pair position")
  }
  
  #Check for duplicated bp
  if(length(unique(marker[,3])) != nrow(marker)){
    warning("Your marker map file contains duplicated ids! Be careful in your analysis!")
  }
  
  #Map passed checks
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(marker), " markers were found for chromosome ",chr,".\n",sep=""))
  }
  
  #Load sample file
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample<- samples.file
    # read.table(samples.file,header=FALSE,colClasses = "character")
  
  #Check if the sample file contains correct dimension
  if(ncol(sample) != 2){
    stop("Sample file contains wrong number of columns (expected 2)")
  }
  
  #Check for duplicated ids
  if(length(unique(sample[,2])) < nrow(sample)){
    stop("Sample file contains duplicated ids!")
  }
  
  pop <- rep(sample[,1],each=2)
  ids <- rep(sample[,2],each=2)
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ", length(unique(pop)), " populations.\n",sep=""))
  }
  
  #Create GHap.phase object
  phase<-NULL
  phase$chr<-chr
  phase$nsamples<-nrow(sample)
  phase$nmarkers<-nrow(marker)
  phase$nsamples.in<-nrow(sample)
  phase$nmarkers.in<-nrow(marker)
  phase$pop<-pop
  phase$id<-ids
  phase$id.in<-rep(TRUE,times=length(phase$id))
  phase$marker<-marker[,2]
  phase$marker.in<-rep(TRUE,times=length(phase$marker))
  phase$bp<-marker[,3]
  phase$A0<-marker[,4]
  phase$A1<-marker[,5]
  if(verbose == TRUE){
    cat("Reading in phased genotypes... (may take a few minutes for large datasets)\n")
  }
  phase$phase<- phase.file
    # read.big.matrix(filename = phase.file,sep = " ",header = FALSE,type = "char")
  
  #Check phase file dimensions
  
  if(nrow(phase$phase) != phase$nmarkers & ncol(phase$phase) != 2*phase$nsamples){
    stop("Your phased genotypes file contains wrong dimensions")
  }
  
  #Return ghap object
  class(phase) <- "GHap.phase"
  if(verbose == TRUE){
    cat("Your GHap.phase object was successfully loaded without apparent errors.\n\n")
  }
  return(phase)
  
}

#Function: ghap.haplotyping
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Output haplotype genotype matrix for user-defined haplotype blocks

ghap.haplotyping<-function(
  phase,
  blocks,
  outfile,
  freq=c(0,1),
  drop.minor=FALSE,
  batchsize=500,
  ncores=1,
  verbose=TRUE
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  #Insert suffix to outfile
  hapgenotypes <- paste(outfile,"hapgenotypes",sep=".")
  hapalleles <- paste(outfile,"hapalleles",sep=".")
  hapsamples <- paste(outfile,"hapsamples",sep=".")
  
  #Check if output will overwrite existing files before opening connection
  # if(file.exists(hapsamples) == TRUE){
  #   stop(paste("File", hapsamples, "already exists!"))
  # }
  # if(file.exists(hapgenotypes) == TRUE){
  #   stop(paste("File", hapgenotypes, "already exists!"))
  # }
  # if(file.exists(hapalleles) == TRUE){
  #   stop(paste("File", hapalleles, "already exists!"))
  # }
  
  #Identify activated samples
  ids.in <- which(phase$id.in)
  id <- phase$id[ids.in]
  id <- id[1:length(id) %% 2 == 0]
  pop <- phase$pop[ids.in]
  pop <- pop[1:length(pop) %% 2 == 0]
  ids.n <- length(id)
  
  #Output hapsamples file
  write.table(x = cbind(pop,id),file = hapsamples,quote = FALSE,row.names = FALSE,col.names=FALSE)
  
  #Generate batch index
  if(batchsize > nrow(blocks)){
    batchsize <- nrow(blocks)
  }
  id1<-seq(1,nrow(blocks),by=batchsize)
  id2<-(id1+batchsize)-1
  id1<-id1[id2<=nrow(blocks)]
  id2<-id2[id2<=nrow(blocks)]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,nrow(blocks))
  if(id1[length(id1)] > nrow(blocks)){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", nrow(blocks), " blocks in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
  }
  
  
  block.iter.FUN<-function(i){
    
    outline<-NULL
    
    #Get block info
    block.info <- blocks[i,c("BLOCK","CHR","BP1","BP2")]
    
    #SNPs in the block
    snps <- which(phase$bp >= .subset2(block.info,3) & phase$bp <= .subset2(block.info,4) & phase$marker.in == TRUE)
    phase.A0 <- phase$A0[snps]
    phase.A1 <- phase$A1[snps]
    
    if(length(snps) >= 1){
      
      #Subset block
      if(length(snps) == 1){
        haplotypes <- as.character(phase$phase[snps,ids.in])
      }else{
        block.subset <- phase$phase[snps,ids.in]
        haplotypes <- apply(block.subset,MARGIN = 2, paste, collapse="")
      }
      
      #Haplotype library
      lib <- table(haplotypes)/(2*ids.n)
      lib <- sort(lib)
      if(drop.minor == TRUE & length(lib) > 1){
        lib <- lib[-1]
      }
      lib <- lib[lib >= freq[1] & lib <= freq[2]]
      
      #Output genotype matrix
      if(length(lib) >= 1){
        for(j in names(lib)){
          phase.geno <- as.integer(haplotypes[1:length(haplotypes) %% 2 == 0] == j)
          phase.geno <- phase.geno + as.integer(haplotypes[1:length(haplotypes) %% 2 == 1] == j)
          j <- unlist(strsplit(j,""))
          phase.allele <- rep(NA,times=length(j))
          phase.allele[j == "0"] <- phase.A0[j == "0"]
          phase.allele[j == "1"] <- phase.A1[j == "1"]
          phase.allele <- paste(phase.allele,collapse="")
          phase.info<-paste(paste(block.info,collapse=" "),phase.allele,sep=" ")
          phase.geno<-paste(phase.geno,collapse=" ")
          outline <- c(outline,phase.info,phase.geno)
        }
      }
      
    }
    return(outline)
  }
  
  
  #Iterate blocks
  nblocks.done <- 0
  for(i in 1:length(id1)){
    
    hapgenotypes.con  <- file(hapgenotypes, open = "a")  
    hapalleles.con  <- file(hapalleles, open = "a")  
    
    #Compute blocks
    ncores <- ncores
      # min(c(detectCores(),ncores))
    if(Sys.info()["sysname"] == "Windows"){
      cat("\nParallelization not supported yet under Windows (using a single core).\n")
      mylines <- unlist(lapply(FUN = block.iter.FUN, X = id1[i]:id2[i]))
    }else{
      mylines <- unlist(mclapply(FUN = block.iter.FUN, X = id1[i]:id2[i], mc.cores = ncores))
    }
    
    #Write batch to files
    if(is.null(mylines) == F){
      writeLines(text = mylines[1:length(mylines) %% 2 == 1],con=hapalleles.con)
      writeLines(text = mylines[1:length(mylines) %% 2 == 0],con=hapgenotypes.con)
    }
    #Log message
    if(verbose == TRUE){
      nblocks.done <- nblocks.done + (id2[i]-id1[i]) + 1
      cat(nblocks.done, "blocks written to file\r")
    }
    
    #Close connections
    close(con = hapalleles.con)
    close(con = hapgenotypes.con)
    
  }
  
  
}

#Function: ghap.blockgen
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Generate blocks based on sliding windows

ghap.blockgen<-function(
  phase,
  windowsize=10,
  slide=5,
  unit="marker",
  nsnp=2
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  if(unit %in% c("marker","kbp") == FALSE){
    stop("Unit must be specified as 'marker' or 'kbp'")
  }
  
  #Initialize vectors
  BP1 <- rep(NA,times=phase$nmarkers.in)
  BP2 <- rep(NA,times=phase$nmarkers.in)
  SIZE <- rep(NA,times=phase$nmarkers.in)
  NSNP <- rep(NA,times=phase$nmarkers.in)
  
  if(unit == "kbp"){
    
    #Transform windowsize to bp
    windowsize <- windowsize*1e+3
    slide <- slide*1e+3
    
    #Get bp from active markers
    markersin <- which(phase$marker.in)
    bp <- phase$bp[markersin]
    
    #Generate indices
    id1 <- seq(1,max(bp),by=slide)
    id2 <- id1 + windowsize
    id1 <- id1[id2 <= max(bp)]
    id2 <- id2[id2 <= max(bp)]
    for(i in 1:length(id1)){
      slice <- which(bp >= id1[i] & bp <= id2[i])
      BP1[i] <- id1[i]
      BP2[i] <- id2[i]
      NSNP[i] <- length(slice)
    }
    
  }else if(unit == "marker"){
    
    id1<-seq(1,phase$nmarkers.in,by=slide);
    id2<-id1+(windowsize-1);
    id1<-id1[id2<=phase$nmarkers.in]
    id2<-id2[id2<=phase$nmarkers.in]
    markersin <- which(phase$marker.in)
    for(i in 1:length(id1)){
      slice <- markersin[id1[i]:id2[i]]
      BP1[i] <- phase$bp[slice[1]]
      BP2[i] <- phase$bp[slice[length(slice)]]
      NSNP[i] <- length(slice)
    }
    
  }
  
  results <- data.frame(BP1,BP2,NSNP,stringsAsFactors = FALSE)
  results <- unique(results)
  results$SIZE <- results$BP2 - results$BP1
  results$SIZE[results$NSNP == 1] <- 1
  results <- results[order(results$BP1,results$BP2),]
  results <- results[results$NSNP >= nsnp,]
  if(nrow(results) == 0){
    stop("No blocks could be generated with the specified options. Try setting the nsnp argument to a smaller value.")
  }
  results <- na.exclude(results)
  results$BLOCK <- paste("CHR",phase$chr,"_B",1:nrow(results),sep="")
  results$CHR <- phase$chr
  results <- results[,c("BLOCK","CHR","BP1","BP2","SIZE","NSNP")]
  return(results);
  
}

#Function: ghap.loadhaplo
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Load haplotype genotypes

ghap.loadhaplo<-function(
  hapsamples.file,
  hapalleles.file,
  hapgenotypes.file,
  verbose = TRUE
){
  
  #Check files
  if(file.exists(hapgenotypes.file) == FALSE){
    stop("Could not find haplotype genotypes file")
  }
  if(file.exists(hapalleles.file) == FALSE){
    stop("Could not find haplotype alleles file")
  }
  if(file.exists(hapsamples.file) == FALSE){
    stop("Could not find sample file")
  }
  
  #Load haplotype alleles file
  if(verbose == TRUE){
    cat("\nReading in haplotype allele information... ")
  }
  hapalleles<-read.table(hapalleles.file,header=FALSE,colClasses = c("character","character","numeric","numeric","character"))
  
  #Check if the haplotype alleles file contains correct dimension
  if(ncol(hapalleles) != 5){
    stop("Haplotype alleles file contains wrong number of columns (expected 5)")
  }else{
    if(verbose == TRUE){
      cat("Done.\n")
      cat(paste("A total of ", nrow(hapalleles), " haplotype alleles were found.\n",sep=""))
    }
  }
  
  #Load sample file
  if(verbose == TRUE){
    cat("Reading in sample information... ")
  }
  sample<-read.table(hapsamples.file,header=FALSE,colClasses = "character")
  
  #Check if the sample file contains correct dimension
  if(ncol(sample) != 2){
    stop("Sample file contains wrong number of columns (expected 2)")
  }
  
  #Check for duplicated ids
  if(length(unique(sample[,2])) < nrow(sample)){
    stop("Sample file contains duplicated ids!")
  }
  
  pop<-sample[,1]
  ids<-sample[,2]
  if(verbose == TRUE){
    cat("Done.\n")
    cat(paste("A total of ", nrow(sample), " individuals were found in ", length(unique(pop)), " populations.\n",sep=""))
  }
  
  #Create GHap.hap object
  hap<-NULL
  hap$nsamples<-nrow(sample)
  hap$nalleles<-nrow(hapalleles)
  hap$nsamples.in<-nrow(sample)
  hap$nalleles.in<-nrow(hapalleles)
  hap$pop<-pop
  hap$id<-ids
  hap$id.in<-rep(TRUE,times=length(hap$id))
  hap$chr<-hapalleles[,2]
  hap$block<-hapalleles[,1]
  hap$bp1<-hapalleles[,3]
  hap$bp2<-hapalleles[,4]
  hap$allele<-hapalleles[,5]
  hap$allele.in<-rep(TRUE,times=length(hap$allele))
  if(verbose == TRUE){
    cat("Reading in haplotype genotypes... (may take a few minutes for large datasets)\n")
  }
  hap$genotypes<- read.table(file = hapgenotypes.file,sep = " ",header = FALSE,colClasses = "character")
    # read.big.matrix(filename = hapgenotypes.file,sep = " ",header = FALSE,type = "char")
  
  #Check haplotype genotypes file dimensions
  if(nrow(hap$genotypes) != hap$nalleles & ncol(hap$genotypes) != 2*hap$nsamples){
    stop("Your haplotype genotypes file contains wrong dimensions")
  }
  
  #Return GHap.haplo object
  class(hap) <- "GHap.haplo"
  if(verbose == TRUE){
    cat("Your GHap.haplo object was successfully loaded without apparent errors.\n\n")
  }
  return(hap)
  
}

#Function: ghap.hapstats
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Summary statistics for haplotype alleles

ghap.hapstats<-function(
  haplo,
  alpha=c(1,1),
  only.active.samples=F,
  only.active.alleles=F,
  ncores=1
){
  
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if inactive alleles and samples should be reactived
  if(only.active.alleles == FALSE){
    haplo$allele.in <- rep(TRUE,times=haplo$nalleles)
    haplo$nalleles.in<-length(which(haplo$allele.in))
  }
  if(only.active.samples == FALSE){
    haplo$id.in <- rep(TRUE,times=haplo$nsamples)
    haplo$nsamples.in<-length(which(haplo$id.in))
  }
  
  #Hapstats iterate function
  hapstats.FUN <- function(j){
    hap.geno <- as.numeric(haplo$genotypes[j,haplo$id.in])
    N<-sum(hap.geno)
    FREQ<-N/(2*haplo$nsamples.in)
    O.HOM <- length(which(hap.geno == 2))
    O.HET <- length(which(hap.geno == 1))
    return(c(N,FREQ,O.HOM,O.HET))
  }
  
  #Compute haplotype allele statistics
  ncores <- min(c(detectCores(),ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    a <- lapply(X = which(haplo$allele.in), FUN = hapstats.FUN)
  }else{
    a <-  mclapply(FUN=hapstats.FUN,X=which(haplo$allele.in),mc.cores = ncores)
  }
  a <- data.frame(matrix(unlist(a), nrow=haplo$nalleles.in, byrow=TRUE))
  hapstats <- NULL
  hapstats$BLOCK <- haplo$block[haplo$allele.in]
  hapstats$CHR <- haplo$chr[haplo$allele.in]
  hapstats$BP1 <- haplo$bp1[haplo$allele.in]
  hapstats$BP2 <- haplo$bp2[haplo$allele.in]
  hapstats$ALLELE<- haplo$allele[haplo$allele.in]
  hapstats$N <- a[,1]
  hapstats$FREQ <- a[,2]
  hapstats$O.HOM <- a[,3]
  hapstats$O.HET <- a[,4]
  hapstats$E.HOM <- (hapstats$FREQ^2)*haplo$nsamples.in
  hapstats$RATIO <- (hapstats$E.HOM+alpha[1])/(hapstats$O.HOM+alpha[2])
  hapstats$BIN.logP <- -1*pbinom(q = hapstats$O.HOM,size = haplo$nsamples.in,prob = hapstats$FREQ^2,lower.tail=TRUE,log.p = TRUE)/log(10)
  hapstats$POI.logP <- -1*ppois(q = hapstats$O.HOM,lambda = hapstats$E.HOM,lower.tail=TRUE,log.p = TRUE)/log(10)
  hapstats <- data.frame(hapstats,stringsAsFactors = FALSE)
  hapstats$TYPE <- NA
  for(i in unique(hapstats$BLOCK)){
    slice <- which(hapstats$BLOCK == i)
    freq <- hapstats$FREQ[slice]
    sumfreq <- sum(freq)
    nalleles <- length(slice)
    type <- rep("REGULAR",times=nalleles)
    type[freq == 0] <- "ABSENT"
    minfreq <- min(freq[freq != 0])
    maxfreq <- max(freq)
    if(sumfreq == 1 & nalleles > 2){
      type[which(freq == minfreq)[1]] <- "MINOR"
      type[which(freq == maxfreq)[1]] <- "MAJOR"
    }else if(sumfreq == 1 & nalleles == 2){
      if(freq[1] == freq[2]){
        type[1] <- "MINOR"
        type[2] <- "MAJOR"
      }else if(freq[1] != freq[2] & pmin(freq[1],freq[2]) != 0){
        type[which(freq == minfreq)] <- "MINOR"
        type[which(freq == maxfreq)] <- "MAJOR"
      }else if(maxfreq == 1){
        type[which(freq == 1)] <- "SINGLETON"
      }
    }else if(sumfreq == 1 & nalleles == 1){
      type <- "SINGLETON"
    # }else if(sumfreq != 1 & nalleles == 1){
    #   type <- "REGULAR"
    # }else{
    #   type[which(freq == maxfreq)[1]] <- "MAJOR"
    # }
    }else{
      type <- "REGULAR"
    }
    hapstats$TYPE[slice] <- type
  }

  #Return object
  return(hapstats)
  
}

#Function: ghap.blockstats
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Calculate block summary statistics

ghap.blockstats <- function(hapstats, ncores = 1){
  
  #Get unique blocks
  blocks <- unique(hapstats[,c("BLOCK","CHR","BP1","BP2")])
  
  #Set of internal functions
  my.fun <- function(block){
    freq <- hapstats$FREQ[hapstats$BLOCK==block & hapstats$TYPE != "ABSENT"]
    if(sum(freq) == 1){
      exp.het <- 1-(sum((freq)^2))
    }else{
      exp.het <- NA
    }
    return(c(exp.het,length(freq)))
  }
    
  #Calculation of expected heterozygosity
  ncores <- min(c(detectCores(),ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    temp <- lapply(FUN = my.fun, X = blocks$BLOCK)
  }else{
    temp <- mclapply(FUN = my.fun, X = blocks$BLOCK, mc.cores = ncores)
  }
  temp <- unlist(temp)
  blocks$EXP.H <- temp[1:length(temp) %% 2 == 1]
  blocks$N.ALLELES <- temp[1:length(temp) %% 2 == 0]
  
  #Return output
  return(blocks)
  
}




