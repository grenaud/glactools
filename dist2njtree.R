#!/usr/bin/env Rscript

#DO NOT use, not mature yet

library(stringr)
library("phangorn")

args=(commandArgs(TRUE))

tmf<-tempfile();

cmd1<-paste("grep -n  \"\\-ALL\\-\"  ",args[1],sep="");
print(cmd1);
o<-system(cmd1,intern=TRUE);

n<-strsplit(o,":")[[1]][1]
system(paste("tail -n ",n," ",args[1]," > ",tmf,sep=""))

d<-read.table(tmf);


c<-cbind(str_split_fixed(d$V1, "-", 2),d$V20)

an<-unique(c[,1])

M <- array(0, c(length(an), length(an)), list(an, an))

i <- match(c[,1], an)
j <- match(c[,2], an)


M[cbind(i,j)] <- M[cbind(j,i)] <- as.numeric(c[,3])
tree<-NJ(as.dist(M))
pdf("nj.pdf")
plot(tree, , "unrooted", main = "Neighbor Joining")
dev.off();


write.tree(tree,file="/dev/stdout");
