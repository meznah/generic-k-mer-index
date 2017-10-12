
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>      //for atoi 
#include "query.h"


int main(int argc,char *argv[])
{
  std::string qfilename = "", dbfilename ="", idxfilename="", 
		hitsfilename ="" ;
  UI ksize = 1; // samplingStep=1;
  bool allow_sampling = false; 
  int optind = 1;  
  while ((optind < argc) && (argv[optind][0] == '-')) {
        std::string sw = argv[optind];
        if (sw == "-db") {
            optind++;
            dbfilename = argv[optind];
        } else if (sw == "-idx") {
            optind++;
            idxfilename = argv[optind];
 	} else if (sw == "-qf") {
            optind++;
            qfilename = argv[optind];
        } else if (sw == "-hf") {
            optind++;
            hitsfilename = argv[optind];
	} else if (sw == "-k") {
            optind++;
            ksize = atoi(argv[optind]);
	} else if (sw == "-sample") {
           allow_sampling = true; 
	} else if (sw == "-help") {
        std::cout<<"this program finds all shared MEM between a query sequence(s) and a database using k-mer index \n";
	std::cout<<"\t-qf    query sequence file name, mandatory input \n";
        std::cout<<"\t-db    database file name, mandatory input \n";
        std::cout<<"\t-idx   index file name, mandatory input \n";
	std::cout<<"\t-rf    MEM results file name \n";
        std::cout<<"\t-l     MEM minimum length \n";
	std::cout<<"\t-help  print this help message \n";
       return 0;

        } else { std::cout << "Unknown switch: " << argv[optind] << std::endl; }
        optind++;
    }

  if (hitsfilename == "" ) { hitsfilename = "h." + qfilename; }  

  if (qfilename == "" || dbfilename == "" || idxfilename == "") {
    std::cout << "Error query and db file names must be give \n";
    return -1 ;
    }
  
  clock_t t; t = clock();
  
  query(qfilename, dbfilename+".db", idxfilename+".idx", ksize,
	allow_sampling, hitsfilename); 

  t = clock() - t;
  printf ("total time %Lf seconds \n",((long double)t)/CLOCKS_PER_SEC);
  
  std::cout<<"done!\n";

  return 0;

}
