
#include <iostream>
#include <string>
#include <cstring>      //for atoi 
#include "index.h"

int main(int argc,char *argv[])
{
  INDEX_MODE mode=FIX;
  std::string dbfilename = "", idxfilename="";
  int optind = 1 , method = 1;
  UI ksize=12, samplingStep=1;
  bool printIdxInfo = false;

  while ((optind < argc) && (argv[optind][0] == '-')) {
	std::string sw = argv[optind];
        if (sw == "-db") {
	    optind++;
            dbfilename = argv[optind];
        } else if (sw == "-idx") {
            optind++;
            idxfilename = argv[optind];
	} else if (sw == "-k") {
            optind++;
            ksize = atoi(argv[optind]);
        } else if (sw == "-w") {
            optind++;
            samplingStep = atoi(argv[optind]);
	} else if (sw == "-m") {
            optind++;
            method = atoi(argv[optind]);
	} else if (sw == "-info") {
	    printIdxInfo = true;
	} else if (sw == "-help") {
        std::cout<<"this program creates a RAM based kmer index \n";
        std::cout<<"\t-db    database File Name, mandatory input\n";
        std::cout<<"\t-idx   index name (the default is the database name) \n";
        std::cout<<"\t-k     kmer size (must be in [12, 32], the default is 12)\n";
	std::cout<<"\t-w     sampling step (must be a positive integer, the default is 1)\n";
	std::cout<<"\t-m     kmer selection methods are 1-Fixed, 2-Minimizer(Lex,many), and 3-Minimizer(lex,one) (the default is 1)\n";
	std::cout<<"\t-info  print an index file information, (must give an index file name) \n"; 
	std::cout<<"\t-help  print this help message\n";
       return 0;

        } else { std::cout << "Unknown switch: " << argv[optind] << std::endl; }
        optind++;
    }
   
  if (printIdxInfo == true) {
	if (idxfilename == "") { 
		std::cout<<"Error, must give an existing index file name\n";
		return -1; 
		}
	printIndexInfo(idxfilename+".idx"); 
 	std::cout<<"done!\n";
	return 0; 
	}

  if ( dbfilename == "" ) {
	std::cout<<"database filename is mandatory \n";
	return -1;
	}

  if ( ksize <= 0 || samplingStep <= 0 ) { 
	std::cout<<"ERROR k and/or w <=0 \n"; 
	return -2; 
	}
  if ( ksize < 12 || ksize >  32 ) {
        std::cout<<"ERROR k must be in [12,32] \n";
        return -2;
	}

  if (method == 0 ) {method = 1; std::cout << "set selection method to 1 \n";}
  switch (method) {
	case 1: mode=FIX;   break;
	case 2: mode=SMALLEST; break;
	case 3: mode=SMALLESTOPT; break;
	default:
		std::cout << "Invalid selection method (options are 1(Fix), 2(Minimizer(Lex,many)), or 3 (Minimizer(Lex,one))" 
		<< std::endl;
	}
  
  if (idxfilename == "") { idxfilename=dbfilename;}
 
  std::cout<<"create index with (k,w,mode) = ("<<ksize<<","<<samplingStep
		<<","<<method<<")"<<std::endl;

  clock_t t; t = clock();
  makeIndexNoncompact(dbfilename+".db", idxfilename+".idx", ksize, samplingStep, mode);
  t = clock() - t;
  printf ("time  %Lf seconds \n",((long double)t)/CLOCKS_PER_SEC);

  std::cout<<"done!\n";
  return 0;
}
