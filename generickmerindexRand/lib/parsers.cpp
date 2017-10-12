//the file is modeifed from GED lab @ MSU

#include "parsers.h"

IParser* IParser::get_parser(const std::string &inputfile)
{
   std::string filename(inputfile);
   int found = filename.find_last_of(".");

   std::string type = filename.substr(found+1);
   //currentely we are only suporting fasta format 
   //if (type == "fa" || type == "fasta") 
	return new FastaParser(inputfile);
}

FastaParser::FastaParser(const std::string &inputfile) : 
                         infile(inputfile.c_str())
{
   std::string line, seq;
   next_name = "";
   
   assert(infile.is_open());

   bool valid_read = 0;

   while (!valid_read)  {
      line = "";
      seq = "";
      if (next_name == "")  {
        getline(infile, current_read.name);
        assert(current_read.name[0] == '>');
        current_read.name = current_read.name.substr(1);
      }
      else  {
         current_read.name = next_name;
         next_name = "";
      }
      
      while(line[0] != '>' && !infile.eof()) {
         getline(infile, line);
         if (line[0] != '>') {
            seq += line;
         }
      }

      if ((int)seq.find('N') == -1)  {
         valid_read = 1;
      }

      if (line[0] == '>') {
         next_name = line.substr(1);
      } else {
         seq += line;
      }

   }

   one_read_left = false;

   if (infile.eof())
   {
      one_read_left = true;
   }

   current_read.seq = seq;
}

Read FastaParser::get_next_read()
{
   std::string line = "", seq = "";
   Read next_read = current_read;

   if (one_read_left)  {
      one_read_left = false;
      return next_read;
   }

   bool valid_read = 0;

   while (!valid_read)  {
      current_read.name = next_name;
      next_name = "";
      current_read.seq = "";

      getline(infile, line);

      while(line[0] != '>' && !infile.eof())
      {

         if (line[0] != '>')  {
            seq += line;
         }
         getline(infile, line);
      }

      if (line[0] == '>')  {
         next_name = line.substr(1);
      }

      if ((int)seq.find('N') == -1)  {
         valid_read = 1;
      }  else if (infile.eof())  {
         one_read_left = false;
         break;
      }

      current_read.seq = seq;
      seq = "";

      if (infile.eof()) {
         one_read_left = true;
      }

   }

   return next_read;
}

