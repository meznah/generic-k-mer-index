#include "db.h"

void ReadNode::init(){
  _pageNum=0;
  _nameLength=_seqLength=0;
  for (int i=0; i<maxNameLength; i++) {
	_name[i]='*';
	}
  for (int i=0; i<maxSeqLength; i++) {
	_seq[i]='*';
	}
}

ReadNode::ReadNode()  {init();}

ReadNode::~ReadNode(){;};

SeqIdType ReadNode::getPageNum(){return _pageNum;}

UI ReadNode::getNameLength() {return _nameLength;}

UI ReadNode::getSeqLength() {return _seqLength;}

std::string ReadNode::getName() {
  std::string name="";
  for(UI i=0; i<_nameLength; i++) { name+=_name[i]; }
  return name;
}

std::string ReadNode::getSeq() {
  std::string seq="";
  for(UI i=0; i<_seqLength; i++) { seq+=_seq[i]; }
  return seq;
}

void ReadNode::setPageNum(SeqIdType newPageNum)      {
  _pageNum=newPageNum;
}

void ReadNode::setNameLength(UI newNameLength) {
  if (newNameLength>maxNameLength) {
	std::cout<<"Warning newNameLength > maxNameLength, set to max\n";
	_nameLength=newNameLength;
	}
  else {_nameLength=newNameLength;}
}

void ReadNode::setSeqLength(UI newSeqLength)        {
  if (newSeqLength > maxSeqLength) {
	std::cout<<"Warning newSeqLength > maxSeqLength,"
			"set newSeqLength=maxSeqLength\n";
	_seqLength = maxSeqLength;
	}
	else _seqLength=newSeqLength;
}

void ReadNode::setName(UI newNameLength,std::string  newName)       {
  for (UI i=0; (i<newNameLength) && (i<maxNameLength); i++) {
	_name[i]=newName[i];
        }
  if (newNameLength > maxNameLength) {
	std::cout<<"Warning newNameLength ("<<newNameLength<< 
			") > maxNameLength, only first part is saved\n";
	newNameLength=maxNameLength;
        }
  this->setNameLength(newNameLength);
}

void ReadNode::setSeq(UI newSeqLength,std::string  newSeq) {
  for (UI i=0; (i<newSeqLength && i<maxSeqLength ); i++) {
	 _seq[i]=newSeq[i];
        }
  if (newSeqLength > maxSeqLength) {
	std::cout<<"Warning newSeqLength ("<<newSeqLength << 
		") > maxSeqLength, only the first part is saved\n";
                newSeqLength = maxSeqLength;
        }
  this->setSeqLength(newSeqLength);
}

void ReadNode::nullfy(){this->init();}

void ReadNode::printPageNum() {
  std::cout<<_pageNum;
}

void ReadNode::printName()    {
  for (UI i=0; i<_nameLength; i++)      {
	std::cout<<_name[i];
        }
}

void ReadNode::printSeq()             {
  for (UI i=0; i<_seqLength ; i++)      {
	std::cout<<_seq[i];
	}
}

void ReadNode::printPageNum(std::fstream& outFile) {
  outFile<<_pageNum;
}

void ReadNode::printName(std::fstream& outFile)       {
  for (UI i=0; i<_nameLength; i++)      {
	outFile<<_name[i];
	}
}

void ReadNode::printSeq(std::fstream& outFile)        {
   for (UI i=0; i<_seqLength ; i++)      {
	outFile<<_seq[i];
        }
}

//--- create db ----------
void createDB(std::string readsfilename,std::string dbfilename) {
  //save db in binary format 
  dbfilename=dbfilename+".db";
  std::ofstream dbfile;
  dbfile.open(dbfilename.c_str(),std::ios::out|std::ios::binary);
  assert(dbfile.is_open());
  dbfile.clear();
  //save the length disrbution of DB entries
  std::string statfilename=""; statfilename=dbfilename+".lengths";
  std::ofstream statfile;
  statfile.open(statfilename.c_str(),std::ios::out);
  assert(statfile.is_open());
  //start the building the DB and Query set
  ULL numReads = 0, numChars =0, sum =0;
  UI min, max;
  double avg;


  IParser * parser  = IParser::get_parser(readsfilename.c_str());
  Read read;
  ReadNode readnode;
  std::string seq = "";
  
  max= 0; min = 10^6;

  while(!parser->is_complete())  {
	read = parser->get_next_read();
 	numReads++; //read address in db
	readnode.setPageNum(numReads);
        readnode.setName(read.name.size(),read.name);
	readnode.setSeq(read.seq.size(),read.seq);
	writeToDiskRead(dbfile,&readnode,numReads);

	statfile<<read.seq.length()<<std::endl;
	numChars+=read.seq.size();
	sum+=read.seq.size();
	if (read.seq.size() > max) max=read.seq.size();
	if (read.seq.size() < min) min=read.seq.size();
	readnode.nullfy();
    }
  avg = sum / numReads; 
  std::cout<<"\t db is saved in file :"<<dbfilename<<std::endl;

  //write header of the bin file
  writeToDiskHeader(dbfile,&numReads, &numChars, &min, &max, &avg);

  dbfile.close();
  statfile.close();

  return;
}

//------- db IO ----------
void writeToDiskHeader(std::ofstream& dbfile, ULL* numReads,
			ULL* numChars, UI* min, UI* max, double* avg){
  dbfile.seekp(0);
  dbfile.write((char*) numReads, sizeof(ULL));
  dbfile.write((char*) numChars, sizeof(ULL));
  dbfile.write((char*) min, sizeof(UI));
  dbfile.write((char*) max, sizeof(UI));
  dbfile.write((char*) avg, sizeof(double));
}

void readFromDiskHeader(std::ifstream& dbfile,ULL* numReads,
			ULL* numChars, UI* min, UI* max, double* avg){
  dbfile.seekg(0);
  dbfile.read((char*)numReads,sizeof(ULL));
  dbfile.read((char*)numChars,sizeof(ULL));
  dbfile.read((char*) min, sizeof(UI));
  dbfile.read((char*) max, sizeof(UI));
  dbfile.read((char*) avg, sizeof(double));  
}

void writeToDiskRead(std::ofstream& dbfile,ReadNode* read,
				SeqIdType pageNum){
  dbfile.seekp(pageNum*blockSizeRead);
  dbfile.write((char*) read, sizeof(ReadNode));
}

void readFromDiskRead(std::ifstream& dbfile,ReadNode* read,
				SeqIdType pageNum){
  dbfile.seekg(pageNum*blockSizeRead);
  dbfile.read((char*)read,sizeof(ReadNode));
}

//------ db query ---------
void retrieveSeqById(std::ifstream& dbfile, std::unordered_set<SeqIdType>& reads_ids, 
		std::unordered_map<SeqIdType,std::string>& mymap) {
    ReadNode read;
    std::string  seq;
    for (std::unordered_set<SeqIdType>::iterator pi = reads_ids.begin(); 
			pi != reads_ids.end(); pi++) {
        readFromDiskRead(dbfile, &read, *pi);
        seq=read.getSeq();
	mymap.insert ( std::pair<SeqIdType,std::string>(*pi,seq) );
	read.nullfy();
	seq = "";
        }
}

void retrieveSeqNameById(std::ifstream& dbfile, std::unordered_set<SeqIdType>& reads_ids,
                std::unordered_map<SeqIdType,std::string>& mymap) { 
    ReadNode read;
    std::string  Name;
    for (std::unordered_set<SeqIdType>::iterator pi = reads_ids.begin(); 
                        pi != reads_ids.end(); pi++) {
        readFromDiskRead(dbfile, &read, *pi);
        Name=read.getName();
        mymap.insert ( std::pair<SeqIdType,std::string>(*pi,Name) );
        read.nullfy();
        Name = "";
        }
}

std::string retrieveSeqNameById(std::ifstream& dbfile, SeqIdType seqId){

  ReadNode read;
  std::string  name;
  readFromDiskRead(dbfile, &read, seqId);
  name=read.getName();
  return name;

}

std::string retrieveSeqById(std::ifstream& dbfile, SeqIdType seqId){
  ReadNode read;
  std::string  seq;
  readFromDiskRead(dbfile, &read, seqId);
  seq=read.getSeq();
  return seq;

}

void printdbinfo(std::string dbfilename){
  std::cout<<"dbfilename: "<<dbfilename<<std::endl;
  std::ifstream dbfile;
  dbfile.open(dbfilename.c_str(),std::ios::in|std::ios::binary);
  assert(dbfile.is_open());

  ULL* numReads = new ULL;
  ULL* numChars = new ULL;
  UI*	min	= new UI;
  UI*	max	= new UI;
  double* avg	= new double;

  readFromDiskHeader(dbfile,numReads,numChars,min,max,avg);
  std::cout<<"num sequences:\t "<<*numReads<<std::endl;
  std::cout<<"num char:\t"<<*numChars<<std::endl;   
  std::cout<<"min seq:\t"<<*min<<std::endl;
  std::cout<<"max seq:\t"<<*max<<std::endl;
  std::cout<<"avg seq:\t"<<*avg<<std::endl;
  dbfile.close();
 
}
