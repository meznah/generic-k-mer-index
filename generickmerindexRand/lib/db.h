#ifndef DB_H
#define DB_H

#include <string>
#include <assert.h>
#include <unordered_map>
#include <unordered_set> 
#include"parsers.h"
#include"basic.h"


typedef unsigned long long SeqIdType;
#define maxNameLength 100 //to be user defined later 
#define maxSeqLength 1000 //to be user defined later


class ReadNode
{

public:
  ReadNode();
  ~ReadNode();
  SeqIdType  getPageNum();
  UI        getNameLength();
  UI        getSeqLength();
  std::string  getName();
  std::string  getSeq();
  void setPageNum(SeqIdType);
  void setNameLength(UI);
  void setSeqLength(UI);
  void setName(UI,std::string);  
  void setSeq(UI,std::string); 
  void nullfy();
  void printPageNum();
  void printName();
  void printSeq();
  void printPageNum(std::fstream&);
  void printName(std::fstream&);
  void printSeq(std::fstream&);

private:
  SeqIdType _pageNum;
  char _name[maxNameLength];
  char _seq[maxSeqLength];
  UI _nameLength;
  UI _seqLength;
  void init();
};

#define blockSizeRead sizeof(ReadNode)

void createDB(std::string, std::string);
//------ db IO -----------
void writeToDiskHeader(std::ofstream&, ULL*, ULL*, UI*, UI*, double*);
void readFromDiskHeader(std::ifstream&, ULL*, ULL*, UI*, UI*, double*);

void writeToDiskRead(std::ofstream&, ReadNode*, SeqIdType);
void readFromDiskRead(std::ifstream&, ReadNode*, SeqIdType);
//----- db query -------
void retrieveSeqById(std::ifstream&, std::unordered_set<SeqIdType>&, 
	std::unordered_map<SeqIdType,std::string>&);
void retrieveSeqNameById(std::ifstream&, std::unordered_set<SeqIdType>&,
        std::unordered_map<SeqIdType,std::string>&);

std::string retrieveSeqNameById(std::ifstream&, SeqIdType);
std::string retrieveSeqById(std::ifstream&, SeqIdType);

void printdbinfo(std::string dbfilename);

#endif
