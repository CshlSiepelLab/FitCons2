#pragma once

#include "butils/butils.h"
// define Coord as a unsigned32 or signed 64 quantity
// for performance unsigned32 is better, but requires some care 
//	e.g. don't subtract 1 from valeus that might be 0...
//	soe for Coord Index, use index < (a.size() +1) rather than (index-1<a.size())
//typedef long long Coord;

using butils::ulong;
using butils::llong;
using butils::ullong;
//typedef llong Coord;
typedef ulong Coord;
