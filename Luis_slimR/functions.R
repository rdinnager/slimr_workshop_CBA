slimr_output_nucleotides_2 <- function(name = "seqs"){
  slimr_output(paste(sim.subpopulations.individuals.genome1.nucleotides(),sim.subpopulations.individuals.genome2.nucleotides()), name, type = "slim_nucleotides", expression = "slimr_output_nucleotides()")
}

ped <- function(df_ped){
    chromosome1 <- as.character(df_ped[3])
    chromosome2 <- as.character(df_ped[4])
    # chromosome1 <- gsub("[-^1-9]", "0", data[3])
    # chromosome2 <- gsub("[-^1-9]", "0", data[4])
  split_seqs <- strsplit(c(chromosome1, chromosome2), split = "")
  genotypes <- as.data.frame(matrix(nrow =nchar(chromosome1),ncol =2 ))
  genotypes$V1 <-split_seqs[[1]]
  genotypes$V2 <-split_seqs[[2]]
  genotypes_final <- paste0(paste(genotypes[,1],genotypes[,2]),collapse = " ")
  # return(paste(df_ped[2],paste0("o",df_ped["id"],collapse = ""),df_ped[5],df_ped[6],df_ped[1],1,genotypes_final))
  
  return(paste(df_ped[2],paste0("O",df_ped[2],collapse = ""),as.character(1),as.character(1),df_ped[1],1,genotypes_final))
}

# this code returns the same LD output as the software haploview
# this is the input file for the LD analyses
LD_fun <- function(df_LD_fun,pop_size ,show_warnings=T){
  df_LD_fun$V1[df_LD_fun$V1=="Male"]   <- 1
  df_LD_fun$V1[df_LD_fun$V1=="Female"] <- 2
  df_LD_fun$id <- df_LD_fun[,2]
  
  # df_LD_fun$id <- paste0(df_LD_fun[,2],"_",1:(pop_size))
  plink_ped <- apply(df_LD_fun,1,ped)
  haploview <- gsub("G", "1", plink_ped) # converting allele names to numbers
  haploview <- gsub("C", "2", haploview)
  haploview <- gsub("A", "3", haploview) # converting allele names to numbers
  haploview <- gsub("T", "4", haploview)
  write.table(haploview,file = "haploview.ped",quote = F,row.names = F,col.names = F)
  snp_stats  <- read.pedfile_b("haploview.ped",sep = " ",snps = as.character(1:nchar(df_LD_fun[1,3])))
  genotype_pops <- snp_stats$genotypes
  genotype_pop1 <- genotype_pops[1:Ne,]
  # genotype_pop2 <- genotype_pops[(ind_pop1+1):nrow(genotype_pops),]
  snpsum.col_pop1 <- col.summary(genotype_pop1)
  colnames(snpsum.col_pop1) <-  paste0(colnames(snpsum.col_pop1),"_pop1")
  # snpFst <- Fst(genotype_pops,group=substr(rownames(genotype_pops),2,5))[1]
  plink_map_2 <- plink_map
  colnames(plink_map_2) <- c("chr_name","row_name","loc_cM","loc_bp")
  #This is mutual information foe each SNP
  # MI_snps <- genotype_pops@.Data
  # MI_snps <- as.data.frame(matrix(as.double(MI_snps),nrow = nrow(MI_snps),ncol = ncol(MI_snps) ))
  # MI_snps[MI_snps==1] <- 11
  # MI_snps[MI_snps==2] <- 12
  # MI_snps[MI_snps==3] <- 22
  pop_names <- as.data.frame(rep(1, Ne))
  # MI_snps_2 <- as.data.frame(cbind(pop_names,MI_snps))
  # allele_matrix_2 <- allele.count(MI_snps_2)
  # MI <- as.numeric(paste(lapply(allele_matrix_2, mutual_information))) 
  # # AFD_snps <- as.numeric(paste(lapply(allele_matrix_2, AFD_fun))) 
  #This is Shannon index for each SNP
  # sha_snps_pops <- lapply(allele_matrix_2,shannon)
  # shannon_pop1 <- unname(unlist(lapply(sha_snps_pops,"[",1)))
  # shannon_pop2 <- unname(unlist(lapply(sha_snps_pops,"[",2)))
  snp_final <- cbind(plink_map_2,snpsum.col_pop1[,4:9])
  # minor allele frequency (MAF) threshold
  minor <- MAF 
  # Filter on MAF
  use <- with(snp_final, (!is.na(z.HWE_pop1) & MAF_pop1 >= minor)) 
  if(sum(use)<=9){return(list(NA,NA))}
  if(sum(use)>=10){
    # Remove NA's
    use[is.na(use)] <- FALSE 
    snp_final <- snp_final[use,]
    # Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
    genotype <- genotype_pop1[,use]
    # genotype_sample <- sample(1: ncol(genotype),size = ncol(genotype)/2)
    # genotype <- genotype[,genotype_sample]
    snp_loc  <- plink_map[use,] # this is the location of each snp in cM and in bp
    # snp_loc <- snp_loc[genotype_sample,]
    snp_loc$V2 <- as.numeric(snp_loc$V2)
    #this is the mean distance between each snp which is used to determine the depth at which 
    # LD analyses are performed 
    mean_dis <- mean(diff(snp_loc$V4))
    # ld_depth_b <- ceiling((ld_max_pairwise/mean_dis)*2)
    ld_depth_b <- ld_depth
    ld_snps <- ld(genotype,depth = ld_depth_b,stats = "R.squared") #function to calculate LD
    ld_matrix <- as.matrix(ld_snps)
    return(list(ld_matrix,snp_final))
  }
}

read.pedfile_b <- function (file, n, snps, which, split = "\t| +", sep = ".", na.strings = "0", lex.order = FALSE,show_warnings=T){
  r0 <- as.raw(0)
  r1 <- as.raw(1)
  r2 <- as.raw(2)
  r3 <- as.raw(3)
  con <- gzfile(file)
  open(con)
  if (missing(n)) {
    n <- 0
    repeat {
      line <- readLines(con, n = 1)
      if (length(line) == 0) 
        break
      n <- n + 1
    }
    if (n == 0) 
      stop("Nothing read")
    seek(con, 0)
  }
  gen <- missing(snps)
  map <- NULL
  if (!gen) {
    m <- length(snps)
    if (m == 1) {
      map <- read.table(snps, comment.char = "")
      m <- nrow(map)
      if (missing(which)) {
        which <- 1
        repeat {
          snps <- map[, which]
          if (!any(duplicated(snps))) 
            break
          if (which == ncol(map)) 
            stop("No unambiguous snp names found on file")
          which <- which + 1
        }
      }
      else {
        snps <- map[, which]
      }
    }
  }
  else {
    line <- readLines(con, n = 1)
    fields <- strsplit(line, split)[[1]]
    nf <- length(fields)
    if (nf%%2 != 0) 
      stop("Odd number of fields")
    m <- (nf - 6)/2
    seek(con, 0)
  }
  nf <- 6 + 2 * m
  result <- matrix(raw(n * m), nrow = n)
  ped <- character(n)
  mem <- character(n)
  pa <- character(n)
  ma <- character(n)
  sex <- numeric(n)
  aff <- numeric(n)
  rownms <- character(n)
  a1 <- a2 <- rep(NA, m)
  a1m <- a2m <- rep(TRUE, m)
  mallelic <- rep(FALSE, m)
  for (i in 1:n) {
    line <- readLines(con, n = 1)
    fields <- strsplit(line, "\t| +")[[1]]
    to.na <- fields %in% na.strings
    fields[to.na] <- NA
    ped[i] <- fields[1]
    mem[i] <- fields[2]
    pa[i] <- fields[3]
    ma[i] <- fields[4]
    sex[i] <- as.numeric(fields[5])
    aff[i] <- as.numeric(fields[6])
    alleles <- matrix(fields[7:nf], byrow = TRUE, ncol = 2)
    one <- two <- rep(FALSE, m)
    for (k in 1:2) {
      ak <- alleles[, k]
      akm <- is.na(ak)
      br1 <- !akm & a1m
      a1[br1] <- ak[br1]
      a1m[br1] <- FALSE
      br2 <- !akm & (a1 == ak)
      one[br2] <- TRUE
      br3 <- !akm & !a1m & (a1 != ak)
      br4 <- br3 & a2m
      a2[br4] <- ak[br4]
      a2m[br4] <- FALSE
      br5 <- br3 & (a2 == ak)
      two[br5] <- TRUE
      mallelic <- mallelic | !(akm | one | two)
    }
    gt <- rep(r0, m)
    gt[one & !two] <- r1
    gt[one & two] <- r2
    gt[two & !one] <- r3
    result[i, ] <- gt
  }
  close(con)
  if (any(a1m & show_warnings==T)){ 
    warning("no data for ", sum(a1m), " loci")
  }
  mono <- (a2m & !a1m)
  if (any(mono & show_warnings==T)){ 
    warning(sum(mono), " loci were monomorphic")
  }
  if (any(mallelic & show_warnings==T)) {
    result[, mallelic] <- r0
    warning(sum(mallelic), " loci were multi-allelic --- set to NA")
  }
  if (gen) 
    snps <- paste("locus", 1:m, sep = sep)
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      rnames <- paste(ped, mem, sep = sep)
      if (any(duplicated(rnames))) 
        stop("could not create unique subject identifiers")
    }
    else rnames <- mem
  }
  else rnames <- ped
  dimnames(result) <- list(rnames, snps)
  result <- new("SnpMatrix", result)
  if (lex.order) {
    swa <- (!(is.na(a1) | is.na(a2)) & (a1 > a2))
    switch.alleles(result, swa)
    a1n <- a1
    a1n[swa] <- a2[swa]
    a2[swa] <- a1[swa]
    a1 <- a1n
  }
  fam <- data.frame(row.names = rnames, pedigree = ped, member = mem, 
                    father = pa, mother = ma, sex = sex, affected = aff)
  if (is.null(map)) 
    map <- data.frame(row.names = snps, snp.name = snps, 
                      allele.1 = a1, allele.2 = a2)
  else {
    map$allele.1 <- a1
    map$allele.2 <- a2
    names(map)[which] <- "snp.names"
  }
  list(genotypes = result, fam = fam, map = map)
}

# 
# do_LD_analysis <- function(path_folder,list_ld,list_snps,list_files){
#   for( i in 1:nrow(list_files)){
#     ld_columns <- read.table(paste0(getwd(),"/","ld_pop1.csv"),sep = ",",row.names = 1)
#     # ld_columns <- read.table(paste0(path_folder,"/",list_ld[list_files[i,3]]),sep = ",",row.names = 1)
#     colnames(ld_columns) <- rownames(ld_columns)
#     snpsum.col <- read.table(paste0(getwd(),"/","snps_pop1.csv"),header=T,sep = ",")
#     # if(is.na(snpsum.col)){next()}
#     ld_columns <- as.data.frame(as.table(as.matrix(ld_columns)))
#     ld_columns <- ld_columns[-ld_columns$Freq < 0,] #remove cases where LD was not calculated
#     ld_columns$Var1 <- as.numeric(as.character(ld_columns$Var1))
#     ld_columns$Var2 <- as.numeric(as.character(ld_columns$Var2))
#     #determine the distance at which LD was calculated
#     ld_columns$dis <- ld_columns$Var2 - ld_columns$Var1 
#     #remove pairwise LD results that were calculated at larger distances than the required in the 
#     # settings and then filtering and rearraging dataframes to match each other and then merge them
#     test_linkageb <-ld_columns[which(ld_columns$dis<=ld_max_pairwise),] 
#     ldtb <- data.table( test_linkageb , key = "Var1" )
#     ldtc <- data.table( test_linkageb , key = "Var2" )
#     use <- with(snpsum.col, (!is.na(MAF_pop1)))
#     # Remove NA's
#     use[is.na(use)] <- FALSE 
#     # Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
#     snp_loc  <- snpsum.col[use,c("chr_name","row_name","loc_cM","loc_bp")] # this is the location of each snp in cM and in bp
#     dtb <- data.table(snp_loc, key = "loc_bp" )
#     t_locationb <- ldtb[dtb, "row_name", nomatch=0]
#     t_locationc <- ldtc[dtb, "row_name", nomatch=0]
#     test_linkageb <- test_linkageb[order(test_linkageb$Var1),]
#     test_linkageb$locus_a <- t_locationb
#     test_linkageb <- test_linkageb[order(test_linkageb$Var2),]
#     test_linkageb$locus_b <- t_locationc
#     test_linkageb <- test_linkageb[order(test_linkageb$Var1),]
#     test_linkage <- cbind(test_linkageb[,5],test_linkageb[,3], test_linkageb[,4], test_linkageb[,1])
#     colnames(test_linkage) <- c("L1","r^2","Dist","location")
#     test_split <- split(test_linkage,test_linkage$L1) #separating snps into individual dataframes
#     test_split <- lapply(test_split,"[",,-1)
#     # function to calculate the mean r^2 by locus in windows of the length specified in the settings
#     func_rsqr <- function(dfx){stats.bin(dfx$Dist, dfx$`r^2`,breaks = seq(0,ld_max_pairwise,ld_resolution))}
#     # this function is just to obtain the same output as the above funtion so the two outputs match each 
#     # other and therefeore possible to merge them
#     func_loc  <- function(dfx){stats.bin(dfx$Dist, dfx$location,breaks = seq(0,ld_max_pairwise,ld_resolution))}
#     #applying functions and merging dataframes
#     t_rsqr <- lapply(test_split,func_rsqr)
#     t_rsqr_b <-lapply(t_rsqr,"[[",3)
#     t_rsqr_c <- lapply(t_rsqr_b,"[",2,)
#     t_rsqr_d <- lapply(t_rsqr_c, as.data.frame)
#     t_loc <- lapply(test_split,func_loc)
#     t_loc_b <- lapply(t_loc,"[[",3)
#     t_loc_c <- lapply(t_loc_b,"[",2,)
#     t_loc_d <- lapply(t_loc_c, as.data.frame)
#     linkage <- Map(cbind,t_loc_d, t_rsqr_d)
#     linkage_b <- rbindlist(linkage)
#     linkage_b <- as.data.frame(linkage_b)
#     # LD has been now calculated at the same maximum distance and at the same resolution for each snp. 
#     # A vector (based on the maximum distance of LD and its resolution) is used to separate LD 
#     # resolutions into bins so LD can be calculated in each genomic region
#     linkage_b$factor <- 1:(ld_max_pairwise/ld_resolution)
#     linkage_b <- linkage_b[complete.cases(linkage_b),]
#     colnames(linkage_b) <- c("loc","rsqr","factor")
#     linkage_c<-linkage_b[order(linkage_b$factor),]
#     linkage_c$factor<- factor(linkage_c$factor)
#     linkage_c<- as.data.frame(linkage_c)
#     # separating the groups of LD resolutions into individual dataframes
#     linkage_split<-split(linkage_c,linkage_c$factor) 
#     linkage_split<-lapply(linkage_split,as.data.frame)
#     # function to calculate the mean r^2 of each genomic region in windows of length of the region_size
#     fun_rsqr<- function(df){stats.bin(df$loc, df$rsqr,  breaks = seq(0,chromosome_length,region_size))}
#     fun_loc<- function(df){stats.bin(df$loc, df$factor,  breaks = seq(0,chromosome_length,region_size))}
#     link_rsqr<-lapply(linkage_split,fun_rsqr)
#     link_rsqr_b<-lapply(link_rsqr,"[[",3)
#     link_rsqr_c<- lapply(link_rsqr_b,"[",2,)
#     link_rsqr_d<-lapply(link_rsqr_c, as.data.frame)
#     link <- list.cbind(link_rsqr_d)
#     # this is the number of non-synonymous mutations at the end of the simulations
#     use_b <- with(snpsum.col,(RAF_pop1>0)) 
#     # Remove NA's
#     use_b[is.na(use_b)] <- FALSE 
#     # Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
#     NS_sim <- snpsum.col[use_b,"loc_bp"] 
#     # number_nonsyn[,i] <- as.data.frame(table(cut(NS_sim, breaks=seq(0,chromosome_length,region_size))))[,2]
#     # this is to  calculate the means of the SNP's statistics in bins of size of the region_size 
#     snpsum.col_2 <- snpsum.col[use_b,] 
#     snpsum.col_2 <- split(snpsum.col_2,cut(snpsum.col_2$loc_bp,breaks=seq(0,chromosome_length,region_size)))
#     snpsum.col_2 <- as.data.frame(snpsum.col_2[,5:ncol(snpsum.col_2)])
#     # snpsum.col_2 <- lapply(snpsum.col_2,"[",,5:24)
#     snpsum.col_2 <- colMeans(snpsum.col_2,,na.rm=T)
#     snpsum.col_2 <- list.rbind(snpsum.col_2)
#     snpsum.col_mean[[i]] <- snpsum.col_2
#     for (ld_distance in 1: ncol(link) ) {
#       link_iterations[[ld_distance]][,i] <- link[,ld_distance]
#     }
#   }
#   return(list(link_iterations,snpsum.col_mean))
#   # return(list(link_iterations,number_nonsyn,snpsum.col_mean))
# }
# 
# 
# 
