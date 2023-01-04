library(recluster)
library (iodatabase)
#Open the datacoi
datacoi<-read.table("datacoi.txt", h=T, sep="\t")
nrow(datacoi)
sequences1<-read.FASTA("sequencesGLQ.fas")
length(sequences1)

