# include "query.h"

void query(std::string qfilename, std::string dbfilename, std::string idxfilename,
                 UI ksize, bool allow_sampling, std::string hitsfilename ) {

  if (idxfilename == "" ) {
    std::cout<<"ERROR must send the index name\n";
    return ;
    }

  printIndexInfo(idxfilename);
  IndexHeader idxHeader;
  loadIndexHeader(idxfilename, idxHeader);
  if (idxHeader.get_ksize() > ksize) {
    std::cout << "Error index kmer > query kmer \n";
    return ;
    }

  queryNoncompactWithMEM(qfilename, dbfilename, idxfilename,
        ksize, allow_sampling, hitsfilename);
}

//------------ MEM search ------
void queryNoncompactWithMEM(std::string qfilename, std::string dbfilename, std::string idxfilename,
        UI ksize, bool allow_sampling, std::string hitsfilename ) {

  std::cout<<"MEM search \n";
  std::ofstream hitsfile;
  hitsfile.open(hitsfilename.c_str(),std::ios::out);
  std::string statfilename = hitsfilename + ".stat";
  std::ofstream statfile;
  statfile.open(statfilename.c_str(),std::ios::out);

  std::ifstream dbfile;
  dbfile.open(dbfilename.c_str(),std::ios::in|std::ios::binary);
  assert(dbfile.is_open());


  ReadNode dbreadnode;
  UI dbcntError = 0;

  std::ifstream idxfile;
  idxfile.open(idxfilename.c_str(),std::ios::in|std::ios::binary);
  assert(idxfile.is_open());

  IndexHeader idxheader;
  std::unordered_map<KmerType, SeqPosSet > index;
  clock_t t; t = clock();
  loadIndexNoncompact(idxfile, idxheader, index);
  t = clock() - t;
  printf ("index time %Lf seconds \n",((long double)t)/CLOCKS_PER_SEC);
  std::unordered_map<KmerType, SeqPosSet >::iterator got;
  UI kstar = ksize;
  UI kprime = idxheader.get_ksize();

  IParser * parser  = IParser::get_parser(qfilename.c_str());
  Read qread;
  UI qcntError = 0;
  SeqPosSet seq_pos_set;
  SeqPosPair set_pos_pair;

  KmerType qkmer;
  KmerPosSet qkmerslist;

  UI samplingStep;
  UI qpos , dbpos;
  SeqIdType dbseq;
  PosPairSet posset; //to aggregate kmers for a hit in map
  std::map<SeqIdType,PosPairSet> mymap; //collect all kmers from one hit
  std::map<SeqIdType,PosPairSet>::iterator map_it;
  PosPair pos_pair;
  if (!allow_sampling) { samplingStep = 1;}
  else if (allow_sampling == true) { samplingStep = idxheader.get_samplingStep(); }
  //std::cout<<"query samplingStep = "<<samplingStep<<std::endl;
  UI tp = 0, fp = 0;
  t = clock();
  while(!parser->is_complete())  {
   qread = parser->get_next_read();
   if (qread.seq.length() < kstar) {qcntError++; continue;}

   //step1 : extract all query kmers
      qkmerslist.clear();
   if (idxheader.get_mode() == FIX) {
        qkmerslist = extractSeqSeedFixedWithLocations(qread.seq, kprime, samplingStep);
   }
   else if ((idxheader.get_mode() == SMALLEST) || (idxheader.get_mode()== SMALLESTOPT)) {
        qkmerslist = extractSeqSeedSmallesWithLocations (qread.seq, kprime, samplingStep);
   }

  //step2: find all shared kmers and aggregate them based on seq-id
     mymap.clear();
   for ( auto it = qkmerslist.begin(); it != qkmerslist.end(); ++it ){
     qkmer = (*it).first;       qpos  = (*it).second;
     got = index.find(qkmer);
     if ( got != index.end()) {
        seq_pos_set.clear();
        seq_pos_set = got -> second;
        for ( auto pi = seq_pos_set.begin(); pi != seq_pos_set.end(); ++pi) {
           dbseq = (*pi).first; dbpos=(*pi).second;
           map_it = mymap.find(dbseq);
           if (map_it != mymap.end()) {map_it->second.insert(std::make_pair(qpos,dbpos)); }
           else {
                  posset.clear(); posset.insert(std::make_pair(qpos,dbpos));
                  mymap.insert ( std::pair<SeqIdType,PosPairSet>(dbseq,posset) );
                }
           } //collect kmers per hit into mymap 
         }// current shared qkmer
        }   //current qkmer

  //step3: for all hits, check MEMs
  //hitsfile << "> " << qread.name <<std::endl;
    for ( auto it = mymap.begin(); it != mymap.end(); ++it) {
     dbseq = it->first; posset = it->second;
     readFromDiskRead(dbfile, &dbreadnode, dbseq);
     if (dbreadnode.getSeqLength() < ksize ) {dbcntError++; continue;}
     filterbyMEM(qread.seq, dbreadnode.getSeq(),posset,kstar,kprime,
                qread.name,dbreadnode.getName(),hitsfile,tp, fp, statfile);
        } // end currnt hit 
   } //current query seq 
  t = clock() - t;
  printf ("query time %Lf seconds \n",((long double)t)/CLOCKS_PER_SEC);
  std::cout<<"#tp " << tp << " #fp " << fp << std::endl;
  if (qcntError>0) { std::cout<<"# qseq with length < k:"<<qcntError<<std::endl;}

  hitsfile.close();     statfile.close();
  idxfile.close();      dbfile.close();

}

void filterbyMEM(std::string seq, std::string cseq, PosPairSet psopair,
                UI kstar, UI kprime, std::string seqname, std::string cseqname,
                std::ofstream& hitsfile, UI& tp_cnt, UI& fp_cnt, std::ofstream& statfile) {

  bool is_discovered_kmer;
  std::vector<MEM> allMEM;
  MEM currMEM;
  //UI tp_cnt=0,fp_cnt=0;
  allMEM.clear();
  for ( auto it = psopair.begin(); it != psopair.end(); ++it) {
        is_discovered_kmer = false;
        for (std::vector<MEM>::iterator pi = allMEM.begin() ; pi != allMEM.end(); ++pi) {
                currMEM = *pi;
                if ( (it->first - currMEM[0]) == (it->second - currMEM[1]) )
                   {is_discovered_kmer = true; tp_cnt++;break;} //?check this if works all times
                }
        if (!is_discovered_kmer) {
                currMEM = kmerExtendToMEM (seq, cseq, kstar, kprime, it->first, it->second);
                if (currMEM[2] >= kstar) {allMEM.push_back(currMEM); tp_cnt++;}
                else {fp_cnt++; }
                }
        }
    //statfile << seqname << "\t"<< cseqname << "\t"  << tp_cnt << "\t" << fp_cnt << std::endl;
    for (std::vector<MEM>::iterator pi = allMEM.begin() ; pi != allMEM.end(); ++pi) {
       currMEM = *pi;
       hitsfile << seqname << "\t"<< cseqname  << "\t"<< currMEM[0] << "\t" << currMEM[1] << "\t" << currMEM[2] << "\n" ;
	//hitsfile << cseqname <<  "\t" << currMEM[1]+1<< "\t"<< currMEM[0]+1 << "\t" <<currMEM[2] << std::endl ;
       		}
}


MEM kmerExtendToMEM(std::string seq, std::string cseq, UI kstar, UI kprime,
                UI seqpos, UI cseqpos){
  MEM newMEM;
  newMEM[0]=seqpos; newMEM[1]=cseqpos; newMEM[2]=kprime;
  std::string substr1, substr2;
  substr1 = seq.substr (seqpos, kprime);
  substr2 = cseq.substr(cseqpos,kprime);

  bool is_right_extendable = true,  is_left_extendable = true;
  UI qindexl, qindexr, dbindexl, dbindexr;

  qindexl = seqpos;  qindexr  = seqpos + kprime -1;
  dbindexl= cseqpos; dbindexr = cseqpos+ kprime -1;

  if (seqpos == 0 || cseqpos == 0) is_left_extendable = false ;
  if ((seqpos+kprime == seq.length()) || (cseqpos+kprime == cseq.length()))
        is_right_extendable = false ;

  while ( is_right_extendable || is_left_extendable ) {
  if (is_right_extendable) {
        if ((qindexr == seq.length()-1 )|| dbindexr == cseq.length()-1) {
                is_right_extendable = false;
                }
        else {
                qindexr++; dbindexr++;
                if  (seq[qindexr] != cseq[dbindexr]) {
                        is_right_extendable = false;
                        qindexr--; dbindexr--; }
                }
        }
   if (is_left_extendable) {
        if (qindexl == 0 || dbindexl == 0) {
                is_left_extendable = false;
                }
        else {
                qindexl --; dbindexl--;
                if (seq[qindexl] != cseq[dbindexl]) {
                        is_left_extendable = false;
                        qindexl++; dbindexl++;}
             }
        }
   }
  newMEM[0]=qindexl ; newMEM[1]=dbindexl;
  newMEM[2]=(qindexr-qindexl)+1;

  return newMEM;
}

