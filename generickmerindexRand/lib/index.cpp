
#include "index.h"

void makeIndexNoncompact(std::string dbfilename, std::string idxfilename,
                UI ksize, UI samplingStep,
                INDEX_MODE mode){
  std::ifstream dbfile;
  dbfile.open(dbfilename.c_str(),std::ios::in|std::ios::binary);
  assert(dbfile.is_open());
  ULL numReads = 0, numChars = 0;
  UI min, max;
  double avg;
  readFromDiskHeader(dbfile,&numReads, &numChars, &min, &max, &avg);
  ReadNode readnode;
  UI seqLen=0, cntError=0, cnt =0;
  std::string  seq="";
 
  std::unordered_map<KmerType, SeqPosSet > mymap;
  std::unordered_map<KmerType, SeqPosSet >::iterator it;
  KmerPosSet newseeds;
  std::unordered_set<KmerType> allseeds;
  KmerType kmer;
  SeqPosPair p;

  for (ULL i=1; i <= numReads; i++) {
        if (i % 250000 == 0 ) std::cout<<"process read # "<<i<<std::endl;
	//if (i % 50000 == 0 ) std::cout<<"process read # "<<i<<std::endl;
        readFromDiskRead(dbfile,&readnode,i);
        if (readnode.getSeqLength() < ksize ) {
                cntError++;continue;}
        seqLen=readnode.getSeqLength();
        seq=readnode.getSeq();
	//std::cout<<seq<<std::endl;
        newseeds.clear();
        if (mode == FIX) {
           newseeds = extractSeqSeedFixedWithLocations(seq, ksize, samplingStep);
        	} 
	else if (mode == SMALLEST){
           newseeds = extractSeqSeedSmallesWithLocations(seq, ksize, samplingStep);
        	}
    	else if (mode == SMALLESTOPT){
	   newseeds = extractSeqSeedSmallesWithLocationsOPT(seq, ksize, samplingStep);
		}
	//std::cout<<"#newseeds " << newseeds.size() << std::endl;
	cnt +=newseeds.size();
        for (KmerPosSet::iterator pi=newseeds.begin(); pi!=newseeds.end(); ++pi){
                kmer = (*pi).first;
                p.first = i; p.second = (*pi).second;
                if (mymap.count(kmer) > 0 ) {
                        it = mymap.find(kmer);
                        (it->second).insert(p);
                        }
                 if (mymap.count(kmer) == 0 ) {
                        SeqPosSet p2;
                        p2.insert(p);
                        mymap.insert(std::pair<KmerType,
                                SeqPosSet > (kmer,p2));
                        }
                }

        readnode.nullfy();      seq="";
        }
  if (cntError>0) {
        std::cout << "# reads with length < k: "<<cntError<<std::endl;
        }
  std::cout << "#kmers: "<< mymap.size() << '\n';
  std::cout << "#positions: "<< cnt << '\n';
  saveIndexNoncompact(idxfilename, ksize, samplingStep, mode, mymap);
  dbfile.close();
}


void saveIndexNoncompact(std::string idxfilename, UI ksize, UI samplingStep,
        INDEX_MODE mode, std::unordered_map<KmerType, SeqPosSet > &mymap) {

   std::ofstream indexfile(idxfilename.c_str(), std::ios::binary);
   std::string statfilename= idxfilename+".stat";
   std::ofstream statfile(statfilename.c_str());
   indexfile.seekp(0); statfile.seekp(0);

   UI _ksize=ksize;
   indexfile.write((const char *) &_ksize, sizeof(_ksize));

   UI _samplingStep=samplingStep;
   indexfile.write((const char *) &_samplingStep,sizeof(_samplingStep));

   INDEX_MODE _mode = mode;
   indexfile.write((const char *) &_mode, sizeof(_mode));

   INDEX_TYPE _type = NONCOMPACT;
   indexfile.write((const char *) &_type, sizeof(_type));

   ULL  _dicsize = mymap.size();
   indexfile.write((const char *) &_dicsize, sizeof(_dicsize));

   KmerType * kmerbuf = new KmerType[_dicsize];
   ULL * acclistsizebuf = new ULL [_dicsize];
   ULL listssize=0;
   unsigned int i = 0;

   std::unordered_map<KmerType, SeqPosSet > ::iterator it;
   for (it = mymap.begin() ; it != mymap.end(); ++it,++i) {
        statfile<<it->first<<"\t"<<it->second.size()<<"\n";
        kmerbuf[i]= it->first;
                if ((ULLONG_MAX-listssize) < it->second.size() ) {
                 std::cout<<"ERROR the sum exceed max value\n";
                 if( remove( idxfilename.c_str() ) != 0 )
                        perror( "Error deleting index file" );
                 else   puts( "Index file successfully deleted" );
                 if( remove( statfilename.c_str() ) != 0 )
                          perror( "Error deleting index stat file" );
                 else     puts( "Index stat file successfully deleted" );
                 return;
                 }
        listssize += it->second.size();
        acclistsizebuf[i]=listssize;
        }


   indexfile.write((const char *) &listssize, sizeof(listssize));
   indexfile.write((const char *) kmerbuf, sizeof(KmerType) * _dicsize);
   delete kmerbuf;
   indexfile.write((const char *) acclistsizebuf, sizeof(ULL) * _dicsize);
   delete acclistsizebuf;
  SeqPosSet seq_pos_pairs;
  SeqPosSet::iterator pi;
  SeqPosPair seq_pos_pair;
  for (it = mymap.begin(); it != mymap.end(); ++it) {
        seq_pos_pairs.clear();
        seq_pos_pairs = (*it).second;
        for (pi = seq_pos_pairs.begin(); pi != seq_pos_pairs.end(); pi++) {
        seq_pos_pair = *pi;
        indexfile.write((const char *) &seq_pos_pair, sizeof(seq_pos_pair));
                }
        }

  indexfile.close(); // *.idx     
  statfile.close(); // *.idx.stat  
}

void loadIndexHeader(std::string idxfilename, IndexHeader& idxheader) {

  std::ifstream idxfile;
  idxfile.open(idxfilename.c_str(),std::ios::out|std::ios::binary);
  assert(idxfile.is_open());
  idxfile.seekg(0);

  idxheader.clear();
  idxfile.seekg(0);
  UI ksize;
  idxfile.read((char *) &ksize, sizeof(ksize));
  idxheader.set_ksize(ksize);
  UI samplingStep;
  idxfile.read((char *) &samplingStep, sizeof(samplingStep));
  idxheader.set_samplingStep(samplingStep);
  INDEX_MODE mode;
  idxfile.read((char *) &mode, sizeof(mode));
  idxheader.set_mode(mode);
  INDEX_TYPE type;
  idxfile.read((char *) &type, sizeof(type));
  idxheader.set_type(type);
  ULL  dicsize;
  idxfile.read((char *) &dicsize, sizeof(dicsize));
  idxheader.set_dicsize(dicsize);
  ULL listssize;
  idxfile.read((char *) &listssize, sizeof(listssize));
  idxheader.set_listssize(listssize);

  idxfile.close();
}

void loadIndexNoncompact(std::ifstream& idxfile, IndexHeader& idxheader,
        std::unordered_map<KmerType, SeqPosSet > & index) {

  idxheader.clear();
  index.clear();

  //read header 

  idxfile.seekg(0);
  UI ksize;
  idxfile.read((char *) &ksize, sizeof(ksize));
  idxheader.set_ksize(ksize);
  UI samplingStep;
  idxfile.read((char *) &samplingStep, sizeof(samplingStep));
  idxheader.set_samplingStep(samplingStep);
  INDEX_MODE mode;
  idxfile.read((char *) &mode, sizeof(mode));
  idxheader.set_mode(mode);
  INDEX_TYPE type;
  idxfile.read((char *) &type, sizeof(type));
  idxheader.set_type(type);
  ULL  dicsize;
  idxfile.read((char *) &dicsize, sizeof(dicsize));
  idxheader.set_dicsize(dicsize);
  ULL listssize;
  idxfile.read((char *) &listssize, sizeof(listssize));
  idxheader.set_listssize(listssize);

  //read dictionary & accumulated lists sizes

  KmerType * kmerbuf = new KmerType[dicsize];
  idxfile.read((char *) kmerbuf, sizeof(KmerType) * dicsize);
  ULL * acclistsizebuf = new ULL [dicsize];
  idxfile.read((char *) acclistsizebuf, sizeof(ULL) * dicsize);

  //read copmact lists & pair kmers with it list and insert into index
   SeqPosSet seq_pos_pairs;
  SeqPosPair seq_pos_pair, _seq_pos_pair;
  ULL size;

 for (ULL i = 0; i < dicsize ; i++) {
    seq_pos_pairs.clear();
    if (i == 0 ) { size = acclistsizebuf[i]; }
    else  {size = acclistsizebuf[i] - acclistsizebuf [i-1];}

    for (ULL j = 0; j < size ; j++) {
        idxfile.read((char *) &seq_pos_pair, sizeof(SeqPosPair));
	seq_pos_pairs.insert(seq_pos_pair);
        }
    index.insert(std::pair<KmerType, SeqPosSet > (kmerbuf[i],seq_pos_pairs));
    }

  delete kmerbuf;
  delete acclistsizebuf;
}

void printIndexInfo(std::string idxfilename) {

  std::ifstream idxfile(idxfilename.c_str(), std::ios::binary);

  UI _ksize;
  idxfile.read((char *) &_ksize, sizeof(_ksize));
  UI _samplingStep;
  idxfile.read((char *) &_samplingStep, sizeof(_samplingStep));
  INDEX_MODE _mode;
  idxfile.read((char *) &_mode, sizeof(_mode));
  INDEX_TYPE _type;
  idxfile.read((char *) &_type, sizeof(_type));
  ULL  _dicsize;
  idxfile.read((char *) &_dicsize, sizeof(_dicsize));
  ULL _listssize;
  idxfile.read((char *) &_listssize, sizeof(_listssize));

  std::cout << "index file: " << idxfilename <<std::endl;
  std::cout << "kmer size: " << _ksize << std::endl;
  std::cout << "sampling step: " << _samplingStep << std::endl;
  std::cout << "mode: " ;
  if (_mode == FIX ) std::cout << "FIX" << std::endl;
  if (_mode == SMALLEST) std::cout << "SMALLEST" << std::endl;
  if (_mode == SMALLESTOPT) std::cout << "SMALLESTOPT" << std::endl;
  std::cout << "type: " ;
  //if (_type == COMPACT)  std::cout << "COMPACT" << std::endl;
  if (_type == NONCOMPACT)  std::cout << "NONCOMPACT" << std::endl;
  std::cout << "dictionary size: " << _dicsize << std::endl;
  std::cout << "total list size: " << _listssize << std::endl;
  idxfile.close();
}


