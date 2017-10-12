#ifndef KMERSELECTION
#define KMERSELECTION

#include <iostream>
#include <string>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "parsers.h"
#include "basic.h"
#include "db.h"

typedef std::pair<KmerType, UI> KmerPosPair;
typedef std::set<KmerPosPair> KmerPosSet;
typedef std::set<UI> PosSet;

KmerPosSet extractSeqSeedFixedWithLocations(std::string , UI, UI);
KmerPosSet extractSeqSeedSmallesWithLocations(std::string , UI, UI);
KmerPosSet extractSeqSeedSmallesWithLocationsOPT(std::string , UI, UI);
#endif
