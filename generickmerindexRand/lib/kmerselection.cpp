
#include "kmerselection.h"

//--- extract kmers end at every wth position

KmerPosSet extractSeqSeedFixedWithLocations (std::string seq,
                UI ksize, UI samplingStep) {

  KmerPosSet newseeds;
  newseeds.clear();
  if (seq.length() < ksize) {
        std::cout<<"ERROR cannot select 1st seed\n";
        return newseeds;
        }

  SeqKMerIterator kmers(seq.c_str(), ksize);
  UI pos = 0 ;
  KmerPosPair p;
  while(!kmers.done()) {
        p.first  = kmers.next();
        p.second = pos;
	//std::cout<<"process kmer "<< _revhash(p.first,ksize)<<" @pos "<< p.second <<std::endl;
        if ((pos+ksize-1) % (samplingStep) == 0) {
                newseeds.insert(p);
		//std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
                }
        pos ++;
        }
  return newseeds;
}

//-------------------------------------
//------  save the right most ties -----------------
KmerPosSet extractSeqSeedSmallesWithLocationsOPT (std::string seq,
        UI ksize, UI samplingStep) {

  if (samplingStep == 1) {
        return extractSeqSeedFixedWithLocations(seq,ksize,samplingStep);
        }
  KmerPosSet newseeds;
  if (seq.length() < ksize) {
        std::cout<<"ERROR cannot select 1st seed\n";
        return newseeds;
        }
  newseeds.clear();

  SeqKMerIterator kmers(seq.c_str(), ksize);
  KmerPosPair p;
  if (seq.length() < samplingStep+ksize-1 ) {
        p.first = kmers.next(); p.second = 0;
        newseeds.insert(p);
        return newseeds;
        }
  //std::cout<<"str="<<seq<<std::endl;
  //std::cout<<"k = " << ksize <<" w = "<< samplingStep << std::endl;
  UI  cur_pos, min_pos;
  PosSet min_pos_set;
  KmerType kmer, min;
  kmer = kmers.next(); cur_pos = 0;
  min = kmer; min_pos = cur_pos;
  //std::cout<<"process kmer "<< _revhash(kmer,ksize)<<" pos "<<cur_pos<<std::endl;
  
  std::map<KmerType,PosSet> mymap;
  std::map<KmerType,PosSet>:: iterator it;
  PosSet posset; 
  mymap.clear();
  
  min = kmer; min_pos_set.clear(); min_pos_set.insert(cur_pos);

  while (cur_pos < samplingStep-1) {
        kmer = kmers.next();    cur_pos++;
        //std::cout<<"process kmer "<< _revhash(kmer,ksize)<<" pos "<<cur_pos<<std::endl;
	        if (kmer > min) {
                it=mymap.find(kmer);
                if (it!=mymap.end()) { 
                   (it -> second).insert(cur_pos);
                        }
                else {
                   posset.clear(); posset.insert(cur_pos);
                   mymap.insert ( std::pair<KmerType,PosSet> (kmer,posset) );
                        }
                }
        if (kmer <= min ) {
                if (kmer == min) { min_pos_set.clear();	min_pos_set.insert(cur_pos); }
                if (kmer < min)  { min = kmer; min_pos_set.clear(); min_pos_set.insert(cur_pos);}
                mymap.clear();
                        }
        }// end 1st window
  for (auto it = min_pos_set.begin(); it!=min_pos_set.end(); ++it ) {
        p.first = min; p.second = *it;
        newseeds.insert(p); min_pos = p.second;
        //std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
		}
    it=mymap.begin();
  for (it=mymap.begin(); it!=mymap.end(); ++it) {
    posset = it -> second;
    for (auto pi = posset.begin(); pi!=posset.end(); ++pi) {
        if (*pi < p.second) { posset.erase(pi);}
                }
   if (posset.size() == 0)  { mymap.erase(it);}
        }
  //std::cout<<"handle the remaining windows\n";
  while(!kmers.done()) {
    while ((!kmers.done()) && (cur_pos - min_pos < samplingStep)) {
        kmer = kmers.next();    cur_pos++; 
        //std::cout<<"process kmer "<< _revhash(kmer,ksize)<<" pos "<<cur_pos<<std::endl;
	    if (kmer > min) {
               it=mymap.find(kmer);
               if (it!=mymap.end()) {(it -> second).insert(cur_pos);}
                else { posset.clear(); posset.insert(cur_pos);
                   mymap.insert ( std::pair<KmerType,PosSet> (kmer,posset) );
                        }
                } //end track potintial mins
            if (kmer == min){
                //std::cout<<"a min" <<  _revhash(p.first,ksize)  << "@pos" << cur_pos <<std::endl;
		//std::cout<<"\tit is added as min but saved for future considertions\n"; 
		it=mymap.find(kmer);
		if (it!=mymap.end()) {(it -> second).insert(cur_pos);}
		else { posset.clear(); posset.insert(cur_pos);
			mymap.insert ( std::pair<KmerType,PosSet> (kmer,posset) );
                        }
                   }
            if (kmer < min) {
                min = kmer;  min_pos = cur_pos;
                p.first = min; p.second = min_pos;
                newseeds.insert(p);
                min_pos = p.second;
                //std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
		mymap.clear();
                        } // update the min
        } // reaches the end of window

   if (cur_pos - min_pos == samplingStep) {
        it=mymap.begin();
        min = it->first; min_pos_set = it->second;
	auto pi = min_pos_set.rbegin(); //get the last element
	p.first = min; p.second = *pi;
	newseeds.insert(p);
  	min_pos = p.second;
	//std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
        mymap.erase (it);
        for (it=mymap.begin(); it!=mymap.end(); ++it) {
          posset = it -> second;
          for (auto pi = posset.begin(); pi!=posset.end(); ++pi) {
            if (*pi < p.second) {posset.erase(pi);}
            if (posset.size() == 0) {mymap.erase(it);}
                } //clean pos older than cur min
              } //chelc pos for every kmer
            } // end adding new min 
        } // end process the read

  return newseeds;
}


//------  save both ties ---------------------------
KmerPosSet extractSeqSeedSmallesWithLocations (std::string seq,
        UI ksize, UI samplingStep) {

  if (samplingStep == 1) {
        return extractSeqSeedFixedWithLocations(seq,ksize,samplingStep);
        }
  KmerPosSet newseeds;
  if (seq.length() < ksize) {
        std::cout<<"ERROR cannot select 1st seed\n";
        return newseeds;
        }
  newseeds.clear();

  SeqKMerIterator kmers(seq.c_str(), ksize);
  KmerPosPair p;
  if (seq.length() < samplingStep+ksize-1 ) {
        p.first = kmers.next(); p.second = 0;
        newseeds.insert(p);
        return newseeds;
        }

  //std::cout<<"str="<<seq<<std::endl;
  //std::cout<<"k = " << ksize <<" w = "<< samplingStep << std::endl;
  UI  cur_pos, min_pos;
  PosSet min_pos_set;
  KmerType kmer, min;
  kmer = kmers.next(); cur_pos = 0;
  min = kmer; min_pos = cur_pos;
  //std::cout<<"process kmer "<< _revhash(kmer,ksize)<<" pos "<<cur_pos<<std::endl;
  std::map<KmerType,PosSet> mymap;
  std::map<KmerType,PosSet>:: iterator it;
  PosSet posset; 
  mymap.clear();
  
  min = kmer; min_pos_set.clear(); min_pos_set.insert(cur_pos);

  while (cur_pos < samplingStep-1) {
        kmer = kmers.next();    cur_pos++;
        //std::cout<<"process kmer "<< _revhash(kmer,ksize)<<" pos "<<cur_pos<<std::endl;
	if (kmer > min) {
                it=mymap.find(kmer);
                if (it!=mymap.end()) { 
		   (it -> second).insert(cur_pos);
			}
		else {
		   posset.clear(); posset.insert(cur_pos);
                   mymap.insert ( std::pair<KmerType,PosSet> (kmer,posset) );
			}
                }
        if (kmer <= min ) {
                if (kmer == min) { min_pos_set.insert(cur_pos); }
                if (kmer < min)  { min = kmer; min_pos_set.clear(); min_pos_set.insert(cur_pos);}
                mymap.clear();
                        }
        }// end 1st window
  for (auto it = min_pos_set.begin(); it!=min_pos_set.end(); ++it ) {
  	p.first = min; p.second = *it;
  	newseeds.insert(p); min_pos = p.second;
	//std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
  		}
  it=mymap.begin();
  for (it=mymap.begin(); it!=mymap.end(); ++it) {
    posset = it -> second;
    for (auto pi = posset.begin(); pi!=posset.end(); ++pi) {
	if (*pi < p.second) { posset.erase(pi);}
		}
   if (posset.size() == 0)  { mymap.erase(it);}
	}
  //std::cout<<"handle the remaining windows\n";
  while(!kmers.done()) {
    while ((!kmers.done()) && (cur_pos - min_pos < samplingStep)) {
	kmer = kmers.next();    cur_pos++; 
	//std::cout<<"process kmer "<< _revhash(kmer,ksize)<<" pos "<<cur_pos<<std::endl;
	if (kmer > min) {
            it=mymap.find(kmer);
            if (it!=mymap.end()) {(it -> second).insert(cur_pos);}
	    else { posset.clear(); posset.insert(cur_pos);
                   mymap.insert ( std::pair<KmerType,PosSet> (kmer,posset) );
			}
            	} 
        if (kmer <= min ) {
            if (kmer == min){
		p.first = min; p.second = cur_pos; 
		newseeds.insert(p); 
		min_pos = p.second;
		//std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
			}
            if (kmer < min) {
		min = kmer;  min_pos = cur_pos;
		p.first = min; p.second = min_pos;
                newseeds.insert(p);
		min_pos = p.second;
		//std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
                	}
            mymap.clear();
		} // update the min
	} // reaches the end of window
   if (cur_pos - min_pos == samplingStep) {
        it=mymap.begin();
        min = it->first; min_pos_set = it->second;
	for (auto pi = min_pos_set.begin(); pi!=min_pos_set.end(); ++pi ) {
	    p.first = min; p.second = *pi;
            newseeds.insert(p);
	    min_pos = p.second;
	    //std::cout<<"a min " <<  _revhash(p.first,ksize)  << "@pos" << p.second <<std::endl;
		}
        mymap.erase (it);
	for (it=mymap.begin(); it!=mymap.end(); ++it) {
    	  posset = it -> second;
    	  for (auto pi = posset.begin(); pi!=posset.end(); ++pi) {
            if (*pi < p.second) {posset.erase(pi);}
   	  if (posset.size() == 0) {mymap.erase(it);}
        	}
	      }
	    } 
        } // end process the read

  return newseeds;
}

