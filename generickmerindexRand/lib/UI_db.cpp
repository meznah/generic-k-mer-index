
#include <iostream>
#include <string>
#include <cstring>	//for atoi 
#include "db.h"

int main(int argc,char *argv[]) 
{
  std::string readsfilename = "", dbfilename = "";
  int optind=1;
  bool info = false;

  while ((optind < argc) && (argv[optind][0] == '-')) {
        std::string sw = argv[optind];
  	if (sw == "-i") {
            optind++;
            readsfilename = argv[optind];
	} else if (sw == "-o") {
	    optind++;
            dbfilename = argv[optind];
  	} else if (sw == "-info") {
             info = true ;
        } else if (sw == "-help") {
	std::cout<<"this program creates a disk-based database\n";
	std::cout<<"\t-i     sequence fasta-file name, mandatory input\n";
	std::cout<<"\t-o     database name (the default is sequence file name)\n"; 
	std::cout<<"\t-info  print info about database file \n";
	std::cout<<"\t-help  print this message\n";
	return 0;
        } else { std::cout << "Unknown switch: " << argv[optind] << std::endl; }
        optind++;
    }
  
  if (info) {
 	if (dbfilename == "") {std::cout << "ERROR must pass dbfile name\n"; return -1;}
   	printdbinfo(dbfilename+".db");
	std::cout<<"done!\n";
	return 0;
	}
  if ( readsfilename == "" ) { 
	std::cout<<"readsfilename is mandatory \n"; 
	return -1; 
	}
  if (dbfilename == "") { dbfilename=readsfilename;}

  createDB(readsfilename,dbfilename);
 
  std::cout<<"done!\n";
  return 0;
}
