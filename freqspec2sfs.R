#!/usr/bin/env Rscript

# [output of freqspec] [number of individuals in graph] [output pdf]


args=(commandArgs(TRUE))



data<-read.table(args[1]);

N<-as.numeric(args[2]);
outputpdf<-args[3];


d<-round(N*(data$V4/(data$V3+data$V4)))

pdf(outputpdf);
hist(d,main="Site frequency spectrum",xlab="Allele frequency")
dev.off();
