#ifndef QUERY_H
#define QUERY_H

#include <iostream>
#include <string>
#include <assert.h>
#include <array>
#include <vector>
#include <set>
#include "parsers.h"
#include "basic.h"
#include "db.h"
#include "kmerselection.h"
#include "index.h"

typedef std::pair<UI, UI> PosPair;
typedef std::set<PosPair> PosPairSet;
typedef std::array<ULL,3> MEM; //(qpos, dbpos, size)
// ------------------

void query (std::string, std::string, std::string, UI,
                bool, std::string);
//---- MEM search ------
void queryNoncompactWithMEM(std::string, std::string, std::string,
                UI, bool, std::string);
void filterbyMEM(std::string, std::string, PosPairSet, UI, UI,
	 std::string, std::string, std::ofstream&,
	 UI&, UI&, std::ofstream&);
MEM kmerExtendToMEM(std::string, std::string, UI, UI, UI, UI);

#endif
