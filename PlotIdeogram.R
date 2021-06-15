args = commandArgs(trailingOnly=TRUE)
faiFile <- args[1]
bedFile <- args[2]
outFile <- args[3]
#faiFile <- "assembly.orig.fasta.fai"
#bedFile <- "sedef_out/final.sorted.bed.final.filt"
fai <- read.table(faiFile, comment.char="")
bed <- read.table(bedFile, comment.char="")
minSize <- 2000000
contigs <- which(fai$V2 > minSize)
fai <- fai[contigs,]
nContigs <- length(contigs)
maxLen <- max(fai$V2)
fai$V1 <- as.character(fai$V1)
c1 <- which(bed$V1==fai$V1[1])
i <- 1
bin <- 500000

getText <- function(v) {
  if (v == 0) {
    return(paste(1))
  }
  else {
    return(sprintf("%dM",floor(v/1e6)));
  }
}
    
plotChrom <- function(bed, fai, i, bin) {
  print(i)
  maxLen <- max(fai$V2)
  c1 <- which(bed$V1==fai$V1[i])
  c1his <- hist(bed$V2[c1], breaks=seq(0,ceiling(fai$V2[i]/bin))*bin, plot=F)
  plot(c(), ylim=c(0,max(c1his$counts)+2), xlim=c(0,maxLen), ylab="Number of duplications per 0.5Mbp", xlab="Position", xaxt='n')
  rect(c1his$mids-bin/2, 0, c1his$mids+bin/2,c1his$counts)
  mtext(fai$V1[i],side=3,line=-0.5,cex=0.5)
  segments(-1, 1, -1, fai$V2[i])

  if (fai$V2[i] < 50E6) {
    sapply(c(1,pretty(fai$V2[i],n=1)[1]), function(v) mtext(getText(v), side=1, at=v, cex=0.5))
  } else {
  sapply(pretty(c(1,min(fai$V2[i],150E6)),n=3), function(v) mtext(getText(v), side=1, at=v, cex=0.5))
  }   
}
print(nContigs)
pdf(outFile, width=7,height=8)
maxChromPerFile <- 12
nPages <- ceiling(nContigs/maxChromPerFile)
par(mfrow=c(maxChromPerFile/3,3))
par(mar=c(1,1,1,1))
par(bty='n')
for (p in 1:nPages) {
  nToDraw <- min(maxChromPerFile, nContigs-maxChromPerFile*p)
  print(seq((p-1)*maxChromPerFile+1, min(p*maxChromPerFile,nContigs)))
  sapply(seq((p-1)*maxChromPerFile+1, min(p*maxChromPerFile,nContigs)), function(i) plotChrom(bed, fai, i, bin))
}
dev.off()

