#!/usr/bin/env Rscript

library(ggplot2)

args=(commandArgs(TRUE))

#./dstats2pdf [dstat file] [output prefix pdf]



cmd1<-paste("grep -n  \"^$\" ",args[1]," |tail -1 |sed \"s/://g\" ",sep="");
linenumber <- system(cmd1, intern = TRUE)
tmpf<-tempfile();
cmd2<-paste("awk '{if(NR>",linenumber,"){print $0}}' ",args[1]," > ",tmpf);
system(cmd2);
data <- read.table(tmpf,header=FALSE);

#data <- read.table(args[1]);
outprefix<-args[2];



sourcepops<-unique(unlist(strsplit(as.character(data$V1),"@"))[ c(FALSE,TRUE) ]);

for(sourcepop in sourcepops){



    dat<-data[grepl(paste("@",sourcepop,"$",sep=""), data$V1),]
    uniqueids<-unique(unlist(strsplit(as.character(dat$V1),"-"))[ c(TRUE,FALSE) ]);

    for(pop in uniqueids){
        filename<-paste(outprefix,"@",sourcepop,"_",pop,".pdf",sep="");
        print(filename);

        patternname<-paste("^",pop,"-",sep="");
        da<-dat[grepl(patternname, dat$V1),]
                                        #V34,V35,V36
        minG<-min(da$V34,da$V35,da$V36);
        maxG<-max(da$V34,da$V35,da$V36);
                                        #sub(patternname,"",da$V1)
        da$V1<-sub( paste("@",sourcepop,"$",sep="") ,"" ,   sub(patternname,"",da$V1))
        da2<-da;
        da3<-da;
        
        da3$V1 <-factor(da2$V1, levels=da2[order(da$V34), "V1"]);


        p<-ggplot(da3, aes(y=V1, x=V34))+geom_point(stat="identity",color="darkblue") + geom_errorbarh(aes(xmax = da$V35, xmin = da$V36),height = .2) + labs(x = paste("D(X,",pop,";",sourcepop,",Chimp") )+ labs(y = "X" )+theme(text = element_text(size=15));
        ggsave(filename, plot = p)

        
        
    }
}

