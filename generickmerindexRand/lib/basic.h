//the file is modeifed from GED lab @ MSU


#ifndef BASIC_H
#define BASIC_H

#include <iostream>
#include <string.h>
#include <algorithm>   
#include <cmath>	//for round function in hashing large k
#include <assert.h>  
#include <math.h>       // for pow function
#include <set>		// for test examples
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


//this varialbe is a global varialbe  to randmoize hahsing


// test validity
#define is_valid_dna(ch) ((toupper(ch)) == 'A' || (toupper(ch)) == 'C' || \
                            (toupper(ch)) == 'G' || (toupper(ch)) == 'T')
// bit representation of A/T/C/G.
/* KHMER encoding
#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
                            (toupper(ch)) == 'T' ? 1LL : \
                            (toupper(ch)) == 'C' ? 2LL : 3LL)

#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'T' : \
                           (n) == 2 ? 'C' : 'G')
*/
// BLAST encoding
#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
                            (toupper(ch)) == 'C' ? 1LL : \
                            (toupper(ch)) == 'G' ? 2LL : 3LL)

#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'C' : \
                           (n) == 2 ? 'G' : 'T')

//typedef unsigned long long	KmerType; //hashed kmer value
typedef uint64_t		KmerType; //hashed kmer value
typedef unsigned long long 	ULL; //tot #read, tot #char
typedef unsigned int		UI;  //ksize,sampeling_step,locations,min read

//GLOBAL_VARS
extern unsigned int  GLOBAL_SEED;
extern KmerType  GLOBAL_KMER;

KmerType  _hash(const char * kmer, const unsigned int k);
std::string _revhash(KmerType  hash, unsigned int k);


class SeqKMerIterator {
  protected:
    const char * _seq;
    const unsigned char _ksize;
    KmerType _kmer;
    unsigned int  index, length;
    bool initialized;
  public:
   //----define constructor-----
   SeqKMerIterator (const char * seq, unsigned char k) : _seq(seq), _ksize(k) {
	index = 0;
	length = strlen(seq);
        _kmer = 0 ;
 	initialized = false;
	}
   //--define first ----
   KmerType first(KmerType& f) {
	_kmer = _hash(_seq, _ksize);
	index ++;
	return _kmer;
	}
   //--define next -----
   KmerType next(KmerType& f) {
	if (done()) {
        throw std::exception();
      }

      if (!initialized) {
        initialized = true;
        return first(f);
      }
     std::string subseq(_seq, index, index + _ksize); 
     const char * c = subseq.c_str();
     _kmer =  _hash(c, _ksize);
     index++;
     return _kmer;
     }

    KmerType first()  { return first(_kmer);}
    KmerType next()   { return next(_kmer); }
    bool done()                 { return  index + _ksize > length; }
};

//---- hamming distance- -----
//This is better when most bits in x are 0
//It uses 3 arithmetic operations and one comparison/branch per "1" bit in x.
unsigned int popcount(unsigned long long x);
bool popcount_bounded(unsigned long long x,unsigned int error);
bool is_similar(unsigned long long x,unsigned long long y,unsigned int error);


#endif
