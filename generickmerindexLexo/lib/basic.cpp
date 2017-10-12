//the file is modeifed from GED lab @ MSU

#include "basic.h"

//----- two bit represntation -------------
KmerType  _hash(const char * kmer, const unsigned int k)
//KmerType  _hash(char * kmer, const unsigned int k)
{
   unsigned int limit = strlen(kmer);
   KmerType h = 0;
   // if k<= 32 then use bit representation
   if (k <= 32) {
   h |= twobit_repr(kmer[0]);
   for (unsigned int i = 1; (i < k) && (i < limit); i++) {
        h = h << 2;
        h |= twobit_repr(kmer[i]);
   		}
   return h;
	}
   // if k > 32, then do hasing
   // you may need to add differnt function with differnt possible hasing 
   // 1- XOR 2-pick the minimzer  3- universal hasing 4-division or multipcation hashing!
   else  {
   unsigned int num_h = 0 , index = 0;
   num_h = std::round( k / 32.0 );
   KmerType h_ptr[num_h], h;
   for (unsigned int n = 0; n < num_h; n++) {
	h = 0 ;
	h |= twobit_repr(kmer[index]);
	for (unsigned int i = 1; (i < 32) && (index+i < limit); i++) {
        	h = h << 2;
       		h |= twobit_repr(kmer[index+i]);
   		 }
	index +=32;
	h_ptr[n] = h ;
	}
   //final hash is XOR for all parts
   h = h_ptr[0] ;
   for (unsigned int n = 1 ; n < num_h; n++) {
	h ^= h_ptr[n];
   		}
   return h;
	}
}

std::string _revhash(KmerType  hash, unsigned int k)
{
	//this works only for k < 32, add excpetion
    std::string s = "";
    unsigned int val = hash & 3;
    s += revtwobit_repr(val);
    for (unsigned int i = 1; i < k; i++) {
        hash = hash >> 2;
        val = hash & 3;
        s += revtwobit_repr(val);
    }
    std::reverse(s.begin(), s.end());
    return s;
}

//------ hamming distance ------------
unsigned int popcount(unsigned long long x)
{
    unsigned int count;
    for (count=0; x; count++)
        x &= x-1;
    return count;
}

bool popcount_bounded(unsigned long long x,unsigned int error)
{
    unsigned int count;
    count = popcount(x);
    return (count < error+1);
}

bool is_similar(unsigned long long x,unsigned long long y,unsigned int error)
{
  unsigned long long z;
  z = x ^ y; 
  std::cout<<"Z:"<<z<<std::endl;
   std::cout<<"Z:"<< _revhash(z,4)<<std::endl;
  return popcount_bounded(z,error);
 
}
