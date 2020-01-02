#!/applications/R/R-3.3.2/bin/Rscript

# Author: Xiaohui Zhao (xz289@cam.ac.uk)
# Description: For reads that have valid alignments to multiple
#  positions in the genome, extract the best alignment
#  args[1] is the summary stats file for multi_fq10.sam lines
#  args[2] is the input sam file multi_fq10.sam
#  args[3] is halfway saved the .RData file for multi-unique
#  args[4] is the final output sam file without header.
#  This script is invoked by PE_ChIPseq_MNaseSeq_mapping_pipeline_bowtie2.sh

args <- commandArgs(TRUE)
dim1 <- ceiling(read.table(args[1], header=F)[1,1]/20000)
infile <- args[2]
dat_out <- NULL
for(s in 1:dim1){
  print(s)
  dat_in <-  matrix(scan(infile, what=" ", quote="\"", skip=9+(s-1)*20000, nlines=20000), ncol=21, byrow=T)
  colnames(dat_in) <- c("qname", "flag", "chr", "locs", "mapq", paste("tag",c(1:16), sep=""))
  qnames.len <- length(unique(dat_in[,1]))
  test <-  by(dat_in, list(qname=dat_in[,1]), function(x){
    y <- subset(x, select= -qname)
    y})
  testfn <- function(ldat, ind){
    newd <- ldat[[ind]][1:2,]
    find <- as.matrix(newd,ncol=20)
  }
  #ntest <- sapply(c(1:qnames.len), function(x) testfn(test,x))
  finaltest <- NULL
  for(i in 1:qnames.len){
    ntest <- testfn(test, i)
    finaltest <- rbind(finaltest, ntest)
    finaltest
  }
  newnames <- rep(names(test), each=2)
  totaldat <- cbind(newnames, finaltest)
  dat_out <- rbind(dat_out, totaldat)
  dat_out
}
colnames(dat_out) <- c("qname", "flag", "chr", "locs", "mapq", paste("tag",c(1:16), sep=""))
rpio <- dat_out
flags <- rpio[,2]
findex <- which(is.na(flags))
indexfb <- c(findex,findex-1)
rpi <- rpio[-indexfb,]
save(rpi, file=args[3])
check.len0 <- dim(rpi)[1]/2
check.len1 <- length(unique(rpi[,1]))
s0 <- table(rpi[,1])

if(check.len0==check.len1){
  final.rpi <- rpi
}else{
  multi <- names(which(s0!=2))
  dindex <- rpi[,1]%in%multi
  subdat <- rpi[dindex==T,]
  olddat <- rpi[dindex==F,]
  qnames.len <- length(unique(subdat[,1]))
  newsubdat <- subdat[order(subdat[,1]),]
  test <-  by(subdat, list(qname=subdat[,1]), function(x){
    y <- subset(x, select= -qname)
    y})
  subtest <- NULL
  for(i in 1:qnames.len){
    stest <- testfn(test,i)
    subtest <- rbind(subtest,stest)
    subtest
  }
  newnames <- rep(names(test), each=2)
  newsubdat <- cbind(newnames, subtest)
  colnames(newsubdat) <- c("qname", "flag", "chr", "locs", "mapq", paste("tag",c(1:16), sep=""))
  final.rpi <- rbind(olddat,newsubdat)
 }
write.table(final.rpi, file = args[4], row.names = F, col.names=F, quote=F, sep="\t")



