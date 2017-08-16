#!/usr/bin/env Rscript

library("ggplot2");
library("gplots");
library("reshape2");

args=(commandArgs(TRUE))
#paircoacompute2heatmap.R [output pairwise coa] [samples to include, comma delim] [pdf out prefix] [pdf size]




cmd1<-paste("grep -n  \"^$\" ",args[1]," |tail -1 |sed \"s/://g\" ",sep="");
linenumber <- system(cmd1, intern = TRUE)
tmpf<-tempfile();
cmd2<-paste("awk '{if(NR>",linenumber,"){print $0}}' ",args[1]," > ",tmpf);
system(cmd2);
data <- read.table(tmpf,header=FALSE);

#data <- read.table(args[1],header=FALSE);

t1<-gsub("-([0-9][0-9]?)$","_\\1",data$V1)
t2<-sub("-([0-9][0-9]?)-","_\\1-",t1)

#print(t1);
#print(t2);
#quit();
data$V1<-t2

vectorinString<-(strsplit( args[2]  , "," )[[1]]);
#vectorinString<-(strsplit( "15302,15304,S_BantuHerero-2,S_BantuKenya-2,S_Biaka-2,S_Bougainville-2"  , "," )[[1]]);
t1<-gsub("-([0-9][0-9]?)$","_\\1",vectorinString);
vectorinString<-paste(t1,collapse=",");

#vectoroutString<-(strsplit( args[3]  , "," )[[1]]);
#t1<-gsub("-([0-9][0-9]?)$","_\\1",vectoroutString);
#vectoroutString<-paste(t1,collapse=",");



ind1<-unlist(strsplit( as.character(data$V1) , "-" ))[c(TRUE, FALSE)]
ind2<-unlist(strsplit( as.character(data$V1) , "-" ))[c(FALSE,TRUE) ]

#print(ind1);
#quit();
m<-as.matrix(data)

dat<-as.data.frame(cbind(as.matrix(ind1),as.matrix(ind2), m[,seq(2,(length(data)-1))]) ,stringsAsFactors=FALSE)


cols = seq(3,length(dat)-2);
dat[,cols] = apply(dat[,cols], 2, function(x) as.numeric(x));


#lengthindinclude<-length(strsplit( args[2]  , "," )[[1]])

vectorin<-rep(FALSE,length(dat$V1));
#for(indinc in  strsplit( "15302,15304,Andaman,Dinka_A,Mbuti_A,French_A,Papuan_A,Han_A,Yoruba_A,San_A,Mandenka_A,Stuttgart,Loschbour,French_B,Han_B,Mandenka_B,Mbuti_B,Papuan_B,San_B,Yoruba_B,Australian1_B,Dinka_B,Ust"     , "," )[[1]]){
for(indinc in  strsplit( vectorinString     , "," )[[1]]){
    #print(indinc)
    vectorin <- vectorin | dat$V1==indinc;
    #print(vectorin);
}
dat<-dat[vectorin  , ]

#print(dat);
#quit();

#vectorin<-rep(FALSE,length(dat$V1));
vectorin<-rep(FALSE,length(dat$V1));

for(indinc in  strsplit( vectorinString    , "," )[[1]]){
    #print(indinc)
    vectorin <- vectorin | dat$V2==indinc;
}
dat<-dat[vectorin  , ]

#print(dat);
#print(dat$V1)
#print(dat$V2)



#quit();

val<-acast(   dat, V1~V2, value.var="V65")
valmin<-acast(dat, V1~V2, value.var="V66")
valmax<-acast(dat, V1~V2, value.var="V67")
class(val)   <-"numeric"
class(valmin)<-"numeric"
class(valmax)<-"numeric"


val<-    round(val*100,2)
valmin<- round(valmin*100,2)
valmax<- round(valmax*100,2)


                                        #mirror
val[lower.tri(val,diag = FALSE)]<-NA
valmin[lower.tri(valmin,diag = FALSE)]<-NA
valmax[lower.tri(valmax,diag = FALSE)]<-NA



for(j in 1:length(val[1,])){
   for(i in 1:length(val[1,])){
       val[i,j]    = val[j,i] ;
       valmin[i,j] = valmin[j,i] ;
       valmax[i,j] = valmax[j,i] ;       
   }
}

#    val[lower.tri(val,diag = FALSE)]   <- val[upper.tri(val,diag = FALSE)]
#valmin[lower.tri(valmin,diag = FALSE)] <- valmin[upper.tri(valmin,diag = FALSE)]
#valmax[lower.tri(valmax,diag = FALSE)] <- valmax[upper.tri(valmax,diag = FALSE)]

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)



pdf(paste(args[3],"heat.pdf",sep=""),width = as.integer(args[4]), height = as.integer(args[4]));
#png("filename.png",width=2000,height=800)
#pdf("filename.pdf",width=24,height=8)

p<-heatmap.2(val,
  cellnote = val,  # same data set for cell labels
  main = "Pairwise coalescence", # heat map title
  xlab="Individual#1",
  ylab="Individual#2",
  
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  key=FALSE, 
  Colv="NA",
          Rowv="NA",
             lwid=c(0.1,4), lhei=c(0.1,0.4)
          )            # turn off column clustering
dev.off();


pdf(paste(args[3],"reorderheat.pdf",sep=""),width = as.integer(args[4]), height = as.integer(args[4]));
#png("filename.png",width=2000,height=800)
#pdf("filename.pdf",width=24,height=8)

p<-heatmap.2(val,
  cellnote = val,  # same data set for cell labels
  main = "Pairwise coalescence", # heat map title
  xlab="Individual#1",
  ylab="Individual#2",

  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  key=FALSE, 
  lwid=c(0.1,4),
  lhei=c(0.1,0.4)
          )            # turn off column clustering
dev.off();
