#ifndef INDEX_H
#define INDEX_H

#include <iostream>
#include <string>
#include <climits>              //to compute the max value
#include <time.h>
#include <assert.h>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include "parsers.h"
#include "basic.h"
#include "db.h"
#include "kmerselection.h"

enum INDEX_MODE {FIX, SMALLEST, SMALLESTOPT};
enum INDEX_TYPE {NONCOMPACT};

typedef std::pair<KmerType, UI> KmerPosPair;	// also def in kmerselction.h
typedef std::pair<SeqIdType, UI> SeqPosPair;
typedef std::set<KmerPosPair> KmerPosSet;	// also def in kmerselction.h   
typedef std::set<SeqPosPair > SeqPosSet;
typedef std::pair<KmerType, SeqPosSet> KmerSeqPosSetPair;

class IndexHeader {
  protected:
    UI _ksize;
    UI _samplingStep;
    INDEX_MODE _mode;
    INDEX_TYPE _type;
    ULL  _dicsize;
    ULL _listssize;

  public:
    IndexHeader() { clear();}
    UI          get_ksize () {return _ksize;}
    UI          get_samplingStep () {return _samplingStep; }
    INDEX_MODE  get_mode () {return _mode; }
    INDEX_TYPE  get_type () {return _type;}
    ULL         get_dicsize () {return _dicsize; }
    ULL         get_listssize () {return _listssize ;}
    void        set_ksize(UI ksize) { _ksize = ksize; }
    void        set_samplingStep (UI samplingStep) { _samplingStep = samplingStep; }
    void        set_mode (INDEX_MODE mode) {_mode = mode;}
    void        set_type (INDEX_TYPE type) {_type = type; }
    void        set_dicsize (ULL dicsize) { dicsize =_dicsize ;}
    void        set_listssize (ULL listssize) { _listssize =listssize; }
    void        clear() {
        _ksize = _samplingStep = _dicsize = _listssize =0;
        _mode = FIX; _type = NONCOMPACT;}
};

//---- one pass algorithms-------
void makeIndexNoncompact(std::string, std::string, UI , UI, INDEX_MODE);
//----- save index -----------
void saveIndexNoncompact(std::string,UI, UI, INDEX_MODE,
        std::unordered_map<KmerType, SeqPosSet > &);
//----- load index ---------
void loadIndexHeader(std::string, IndexHeader&);
void loadIndexNoncompact(std::ifstream&, IndexHeader&,
        std::unordered_map<KmerType, SeqPosSet > & );
//---- index info -------
void printIndexInfo(std::string);

#endif
