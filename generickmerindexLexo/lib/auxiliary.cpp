#include "auxiliary.h"

//--- preprocess fasta files-----
// split very long physical seq to logical seq of length len
// where the the logical seq are overlapping with ovlp char
// any ambiguous characters are removed 

void parseGenome(std::string infilename, std::string outfilename,
                UI len, UI ovlp) {

  std::fstream infile;
  infile.open(infilename.c_str(),std::ios::in);
  assert(infile.is_open());

  std::fstream outfile;
  outfile.open(outfilename.c_str(),std::ios::out);

 //code modified from GED lab at MSU
  std::string line, seq, next_name = "";
  Read current_read;
  bool done, is_last=false;
  UI seq_cnt = 0, err_cnt= 0, cnt;
  std::string substr;
  ULL start;
  while (!infile.eof() && !is_last) {
    done = false;
    while (!done && !is_last)  {
      line = ""; seq = "";
      if (next_name == "")  {
        getline(infile, current_read.name);
        assert(current_read.name[0] == '>');
        current_read.name = current_read.name.substr(1);
        }
      getline(infile, line);
      while(line[0] != '>' && !infile.eof()) {
         if ((int)line.find('N') == -1) {
                seq += line;
                }
           else {err_cnt++;}
         getline(infile, line);
        }
     if (line[0] == '>') {
        next_name = line.substr(1);
        done = true;
        }
    if (infile.eof()) {is_last=true;}
    } // finish reading the current seq
  
   current_read.seq=seq;        seq_cnt++;
   done = false; start = 0 , cnt = 0;
   while ( !done) {
        if ( (start + len) < current_read.seq.size()){
           substr = current_read.seq.substr (start, len);
           outfile << ">" << current_read.name << "p" << cnt <<"len" << len << std::endl;
           outfile << substr << std::endl;
           cnt ++; start = (start+len)-ovlp;
                }
        if ( (start + len) >= current_read.seq.size()) {
           substr = current_read.seq.substr (start, current_read.seq.size()-1);
           outfile << ">" << current_read.name << "p" << cnt <<"len"
                << current_read.seq.size() - start << std::endl;
           outfile << substr << std::endl;
           done = true;
                }
        } // finish split the current seq
   current_read.name = next_name;
  } // finish parase the file

  //handel last seq
  done =false; start = 0 , cnt = 0;
  while ( !done) {
        if ( (start + len) < current_read.seq.size()){
           substr = current_read.seq.substr (start, len);
           outfile << ">" << current_read.name << "p" << cnt <<"len" << len << std::endl;
           outfile << substr << std::endl;
           cnt ++; start = (start+len)-ovlp;
                }
        if ( (start + len) >= current_read.seq.size()) {
           substr = current_read.seq.substr (start, current_read.seq.size()-1);
           outfile << ">" << current_read.name << "p" << cnt <<"len"
                 << current_read.seq.size() - start<< std::endl;
           outfile << substr << std::endl;
           done = true;
                }
        } // finish split the current seq

  outfile.close();
  infile.close();

}



