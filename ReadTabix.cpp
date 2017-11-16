/*
 * ReadTabix
 * Date: Aug-13-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "ReadTabix.h"

// using namespace std;

ReadTabix::ReadTabix(string file,string indexForFile,string chrName,int start,int end){
    if(!isFile(file)){
	cerr<<"File "<<file<<" does not exist"<<endl;
	exit(1);
    }

    if(!isFile(indexForFile)){
	cerr<<"File "<<indexForFile<<" does not exist"<<endl;
	exit(1);
    }

    tb = new Tabix(file);
    string region = chrName+":"+stringify(start)+"-"+stringify(end);
    if( !tb->setRegion( region) ){
	cerr<<"Cannot relocated to  "<<chrName<<":"<<start<<"-"<<end<<endl;
	exit(1);
    }

    //fpTab=ti_open(file.c_str(),indexForFile.c_str());	       

    // -1 is substracted from the start because, for some reason that is unknown to me (probably for BED coordinates)
    // Heng Li does that in index.c in int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end)
    // to the get the coordinates right
    // I do this for consistency with the command line program   
    //if(start>0)
    //start--;

    //iteratorTab=ti_query(fpTab,chrName.c_str(),start,end); 
}

ReadTabix::ReadTabix(string file,string indexForFile,string chrName){
    if(!isFile(file)){
	cerr<<"File "<<file<<" does not exist"<<endl;
	exit(1);
    }

    if(!isFile(indexForFile)){
	cerr<<"File index "<<indexForFile<<" does not exist"<<endl;
	exit(1);
    }


    // if ((fpTab=ti_open(file.c_str(),indexForFile.c_str() ))== 0) {
    // 	cerr<<"Cannot open file "<<indexForFile<<""<<endl;
    // 	exit(1);	
    // }

    // if (ti_lazy_index_load(fpTab) < 0){
    // 	cerr<<"Cannot load index for file "<<indexForFile<<""<<endl;
    // 	exit(1);	
    // }

    // int start;
    // int end;
    // int tid;
    
    // if (ti_parse_region(fpTab->idx, chrName.c_str(), &tid, &start, &end) == -1) {
    // 	cerr<<"Cannot locate chromosome "<<chrName<<" "<<endl;
    // 	exit(1);
    // }
    if(!tb->setRegion(chrName)){
	cerr<<"Cannot locate chromosome "<<chrName<<" "<<endl;
	exit(1);
    }

    //iteratorTab=ti_queryi(fpTab,tid,start,end);
}

ReadTabix::~ReadTabix(){
    // ti_iter_destroy(iteratorTab);
    // ti_close(fpTab);
    delete tb;
}

const kstring_t * ReadTabix::getKstringPtr(){
    return tb->getKstringPtr();
}

string ReadTabix::getHeader(){
    //adapted from main.c from tabix
    //cout<<"get header"<<endl;
    // const ti_conf_t *idxconf = ti_get_conf(fpTab->idx);
    // ti_iter_t iter = ti_query(fpTab, 0, 0, 0);
    // const char *s;
    // int len;
    // string toreturn="";
    // while ((s = ti_read(fpTab, iter, &len)) != 0) {
    // 	if ((int)(*s) != idxconf->meta_char) 
    // 	    break;
    // 	//toreturn+=string(s);
    // 	toreturn.append(s);
    // 	toreturn.append("\n");
    // 	//fputs(s, stdout); fputc('\n', stdout);
    // }
    // ti_iter_destroy(iter);
    string toreturn;
    tb->getHeader(toreturn);
    return toreturn;
}


void ReadTabix::repositionIterator(string chrName,int start,int end){
    // ti_iter_destroy(iteratorTab);
    // // -1 is substracted from the start because, for some reason that is unknown to me
    // // Heng Li does that in index.c in int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end)
    // // to the get the coordinates right
    // // I do this for consistency with the command line program   
    // if(start>0)
    // 	start--;
    // iteratorTab=ti_query(fpTab,chrName.c_str(),start,end); 
    string region = chrName+":"+stringify(start)+"-"+stringify(end);
    if( !tb->setRegion(region) ){
	cerr<<"Cannot relocated to  "<<chrName<<":"<<start<<"-"<<end<<endl;
	exit(1);
    }
    
}


void ReadTabix::repositionIterator(string chrName){
    // ti_iter_destroy(iteratorTab);
    // int start;
    // int end;
    // int tid;
    
    // if (ti_parse_region(fpTab->idx, chrName.c_str(), &tid, &start, &end) == -1) {
    // 	cerr<<"Cannot locate chromosome "<<chrName<<" "<<endl;
    // 	exit(1);
    // }

    // iteratorTab=ti_queryi(fpTab,tid,start,end); 
    if( !tb->setRegion(chrName)){
	cerr<<"Cannot locate chromosome "<<chrName<<" "<<endl;
	exit(1);
    }

}

bool ReadTabix::readLine(string & line){
    // int length; //useless ?
    // const char *buffer;

    // bool toReturn=(buffer=ti_read(fpTab,iteratorTab,&length));
    // if(toReturn){
    // 	line=string(buffer);
    // }

    
    bool toReturn=tb->getNextLine(line);
    if(!toReturn){
      	line="";
    }


    return toReturn;
}



//The user is responsible for getting the addr of the kstring using getKstringPtr 
bool ReadTabix::readLineKS(){
    return tb->getNextLineKS();
}

