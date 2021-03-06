
library(slimr)
library(dartR)
library(vcfR)
library(parallel)
source("./bernd_Ne/gl2genepop.r")
library(tictoc)
library(future)
tic()
## number of generations to simulate
simngen<-50
## number of snps to simulate
simnsnps <- 5000


#sample GED chromosom
gedchromo <- read.csv("./bernd_Ne/GED_chromosomes.csv") #averaged length of sex chromosome.


total.chromo <- sum(gedchromo$csize)* 1e6  #length in bases

chromo <- sort(sample(0:(total.chromo-1),simnsnps+1, replace=FALSE))
recrates <- (diff(chromo)/1e6)*0.01 #convert to recombination rate (1 cM=0.01 per 1e6 basis)

chromends <- cumsum(gedchromo$csize)*1e6

for (i in 1:length(chromends))
{
  index <- which(chromo[1:(simnsnps-1)]<chromends[i]  & chromo[2:simnsnps]>chromends[i])
  
  recrates[index]<- 0.5
}

#plot(chromo, rep(1,length(chromo)), pch=15)
#abline(v=chromends)


snpends <- 0:4999





slim_script(
  slim_block(initialize(),
             {
               defineConstant("alpha1", 0.1);
               defineConstant("alpha2", 1e-8);
               ##defineConstant("Ne",10);
               ## set the overall mutation rate
               #initializeMutationRate(0.05); 
               
               ##defineConstant("test_vector", slimr_inline(recrates[[slimr_template("index", 0)]], delay = TRUE));
               defineConstant("rates", slimr_inline(recrates));
               defineConstant("ends", slimr_inline(snpends));
               
               
               ##[[slimr_template("index", 0)]], delay = TRUE));
               
               
               ## m1 mutation type: neutral
               defineConstant("L", !!simnsnps);
               defineConstant("simlen", !!simngen);
               
               #how many outputs (simlen=only once at the end)
               defineConstant("every", simlen);
               initializeSLiMOptions(nucleotideBased=T);
               initializeAncestralNucleotides(randomNucleotides(L));
               initializeMutationTypeNuc("m1", 0.5, "f",0.0 );
               m1.mutationStackPolicy = "l";
               m1.convertToSubstitution=T;
               ## g1 genomic element type: uses m1 for all mutations
               defineConstant("mm", matrix(c(0.0,alpha1,0.0,0.0, alpha1,0.0,0.0,0.0, 0.0,0.0,0.0,alpha1, 0.0,0.0,alpha1,0.0 ),nrow=4, ncol=4));
               defineConstant("mm2", matrix(c(0.0,alpha2,0.0,0.0, alpha2,0.0,0.0,0.0, 0.0,0.0,0.0,alpha2, 0.0,0.0,alpha2,0.0 ),nrow=4, ncol=4));
               initializeGenomicElementType("g1", m1, 1.0, mm);
               #initializeGenomicElementType("g1", m1, 1.0);
               ## uniform chromosome of length 100 kb
               initializeGenomicElement(g1, 0, L-1);
               ## uniform recombination along the chromosome
               #initializeRecombinationRate(slimr_template("RR"));
               #calculated recombination rates
               
               crates = rates
               if (slimr_template("RR")!=0.99)  crates = rep(slimr_template("RR"),5000); #change to constant or use GED version
                 
               initializeRecombinationRate(rates = crates, ends = ends)     
               
               ##initializeRecombinationRate(slimr_template("RR")));
               
             }),
  slim_block(1,
             {
               sim.addSubpop("p1",slimr_template("Ne"));
             }),
  

  
  slim_block(1,simngen,late(), {
         
               if (sim.generation / every== round(sim.generation/every)) {
               pre = ""
               if (sim.generation<100) pre="0" 
               if (sim.generation<10) pre="00" 
               nn = paste0( c("./bernd_Ne/Neest/out","_Ne_",slimr_template("Ne"),"_rep_",slimr_template("rep"),"_RR_",slimr_template("RR"),"_",pre,sim.generation,".vcf")) ;
               
               p1.individuals.genomes.outputVCF(nn, simplifyNucleotides=T);
              
               
               ##c = sim.chromosome;
               ##catn("Ancestral: " + c.ancestralNucleotides());
               
               
              ## g = p1.genomes[0];
              ## catn("Derived:   " + g.nucleotides());
               
              ## catn("SNPs: " + paste(g.mutations.nucleotide));
              } ##end every
             }),
  
  slim_block(20, early(), {
    g1.setMutationMatrix(mm2);
  }),
  
  
  
  slim_block(simngen, late(), {
    ##sim.outputFull("model_output.txt");
               sim.simulationFinished();
              ## sim.outputFull();
              ## allIndividuals = p1.individuals;
               
             
             })
) -> script_temp

##slimr_template_info(script_temp)
##script_dummy <- slimr_script_render(script_temp, template=(data.frame(rep=1, Ne=50)))
##slimr_write(script_dummy[[1]], "script_dummy.txt")
##system("./slim script_dummy.txt")


##slimr::slim_run(script_1)

#df <- expand.grid(rep=1:100,  Ne = c(10,20,30,40,50,100), RR=c(0.0001,0.001,0.005, 0.01,0.1,0.50,0.99))


#df <- expand.grid(rep=1:10,  Ne = c(10,50), RR=c(0.001,0.005, 0.01,0.5,0.99))
#df <- expand.grid(rep=1:15,  Ne = c(10,20,30,40,50,100), RR=round(mean(recrates),4)) #just to have a number
#df <- expand.grid(rep=1:3,  Ne = c(10,20), RR=c(0.0001,0.001))
df <- expand.grid(rep=1:3,  Ne = c(10,20), RR=c(0.5,0.99))
#RR = 0.99 = GED recrates!!!
#df <- expand.grid(rep=1:5,  Ne = c(10,20,30,50), RR=c(0.005,0.99))


script_temp_df <- slimr_script_render(script_temp, template =df)
#slim_run(script_temp_df[[1]])

plan(multisession(workers = 2))
sr <- slim_run(script_temp_df,parallel = TRUE) #based on future plan

#####
#for (i in 1:nrow(df))
#{
#slimr_write(script_temp_df[[i]], "script_dummy.txt")
#system("./slim script_dummy.txt")
#}
  

files <- list.files(path="./bernd_Ne/Neest/",pattern="out" )
#gls <- list()
#for (i in 1:length(files)){
#  dumvcf<- read.vcfR(file.path("~/slimsim/out",files[i]),verbose=FALSE)
#  gls[[i]]<- vcfR2genlight(dumvcf)
#  names(gls)[i] <- files[i]
#}  

#unlist(lapply(gls,nLoc))
#plot(1:length(gls),unlist(lapply(gls, function(x) mean(gl.Ho(x)))))


#snplevels <- c(100,500,1000,2000)
snplevels <- c(100,1000)
#sslevels <- c(10,20,30,40,50,100)
sslevels <- c(5,10)
sim <- function(i)

  
  
#for (i in 1:length(gls))
{
  
  dumvcf<- read.vcfR(file.path("./bernd_Ne/Neest/",files[i]),verbose=FALSE)
  gls<- vcfR2genlight(dumvcf)
  namesgls <- files[i]  
  
  
res <- data.frame(rep=NA, Nesim=NA, nGen=NA, nSnp=NA, sample=NA,RR=NA, Nee0=NA, Nee0_low=NA, Nee0_up=NA )
cc <-1
  

Nesim=as.numeric(gsub("out_Ne_([0-9]+).*$", "\\1", namesgls))

nrep <-as.numeric(gsub(".*?rep_([0-9]+).*$", "\\1", namesgls))
nGen <-  as.numeric(gsub(".*_([0-9]+).vcf$", "\\1", namesgls))
nRR <- as.numeric(gsub(".*?RR_(0.[0-9]+).*$", "\\1", namesgls))

glout <- gls
glout <- gl.filter.monomorphs(glout, verbose = 0)
pop(glout)<- rep("A", nInd(glout))

for (ii in 1:length(snplevels))
{
 
  nSnp=snplevels[ii] 
  if(nSnp<=nLoc(glout)) {
for (iii in 1:length(sslevels))
{
  
  ss = sslevels[iii]
  if (ss<=nInd(glout))
  {

gldummy <- glout[sample(1:nInd(glout), ss,replace = FALSE),sample(1:nLoc(glout), nSnp, replace = FALSE)]


gl2genepop(gldummy,paste0("./bernd_Ne//Neest/slimne",i,".gen"))

inf <- readLines("./bernd_Ne/Neest/info")
inf[3]<-paste0("slimne",i,".gen")

inf[6]<- paste0("slimLD",i,".txt")
con <- file(paste0("./bernd_Ne/Neest/info",i),"w")
writeLines(inf, con)
close(con)

system(paste0("./bernd_Ne/Neest/Ne2-1L i:./bernd_Ne/Neest/info",i))

ll <- readLines(paste0("./bernd_Ne/Neest/slimLD",i,".txt"))
index <- grep("Estimated Ne", ll)
nes <- unlist(strsplit(ll[index],split = "[ ]+"))
index <- grep("JackKnife on Samples", ll)
neslow <- unlist(strsplit(ll[index],split = "[ ]+"))
nesup  <- unlist(strsplit(ll[index+1],split = "[ ]+"))
Nee0 <- as.numeric(nes[length(nes)])
Nee0_low <- as.numeric(neslow[length(neslow)])
Nee0_up <- as.numeric(nesup[length(nesup)])

res[cc,]<- c(nrep, Nesim, nGen, nSnp, ss, nRR, Nee0, Nee0_low, Nee0_up)
cc <- cc+1
#res<- c(nrep, Nesim, nGen, nSnp, ss, Nee0, Nee0_low, Nee0_up)


#cat(paste("run",cc,"\n"))
#flush.console()
} #endif ss<= nInd
} #end iii
  } #if snplevels < nLoc
} #end ii (nSnps)
return(res)
} #gls files

system.time(res <- mclapply(1:length(files), function(x) sim(x), mc.cores = 1))
res<- do.call("rbind",res)

#write.csv(res,"res_all100.csv")

#res <- res[complete.cases(res[,1:7]),]


#res <- res[complete.cases(res),]
ggplot(res, aes(x=factor(Nesim), y=Nee0/Nesim,fill=factor(sample)))+geom_boxplot()+facet_grid(RR~nSnp)+ylim(c(0,1))


#ggplot(resall, aes(x=factor(Nesim), y=Nee0,fill=factor(run)))+geom_boxplot()+facet_wrap(.~nSnp)+ylim(c(0,60))

toc()



