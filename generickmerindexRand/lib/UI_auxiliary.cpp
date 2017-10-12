
#include <iostream>
#include <string>
#include <cstring>	//for atoi 
#include "auxiliary.h"

void printMesage();

int main(int argc,char *argv[]) 
{
  std::string infilename = "", outfilename = "", 
		outfilename2 = "";
  int optind = 1;
  bool parse_g=true;
  UI len=100, ovlp=0; 
  while ((optind < argc) && (argv[optind][0] == '-')) {
        std::string sw = argv[optind];
  	if (sw == "-i") {
            optind++;
            infilename = argv[optind];
	} else if (sw == "-o") {
	    optind++;
            outfilename = argv[optind];
	} else if (sw == "-len") {
             optind++;
             len = atoi(argv[optind]);
        } else if (sw == "-ovlp") {
             optind++;
             ovlp = atoi(argv[optind]);
        } else if (sw == "-help") {
	std::cout<<"this program splits long sequences in fasta format to smaller sequences \n";
	std::cout<<"\t-i     input file name\n";
	std::cout<<"\t-o     output file name \n"; 
	std::cout<<"\t-len   upper bound on sequence length (the defult 100) \n";
	std::cout<<"\t-ovlp  size of overlapping (the default 0, no overlapping)\n";
	std::cout<<"\t-help  print this message \n";
	return 0;
        } else { std::cout << "Unknown switch: " << argv[optind] << std::endl; }
        optind++;
    }
      
   if (infilename == "" || outfilename == "") {
	std::cout<<"Error, must pass infile and outfile names\n";
	return -1;
	} 

   if (parse_g == true) {
        parseGenome(infilename,outfilename,len,ovlp);
        return 0;
        }

  std::cout<<"done!\n";
  return 0;
}

