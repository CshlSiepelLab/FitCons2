#pragma once

#include <cstdio>
#include <vector>

#include "butils/butils.h"		// fastfile namespace
#include "insdb/insdb.h"		// Base class Locus

// For vectors of objects derived from the locus class, get num loci and number of positions covered by elements in vector.
template<typename T> 
void WLociInfo( const std::vector<T> &V, butils::ulong &NumLoci, butils::ulong &NumPos, butils::ulong &Weight) {
	ulong tmp;
	NumLoci = NumPos= Weight = (ulong) 0;
	for (std::size_t i = 0; i < V.size(); i++) {
		NumLoci++; tmp = (ulong)(V[i].en - V[i].st); NumPos += tmp; Weight += (tmp*V[i].counts);	}
}


class WLocus : public insdb::Locus {
public:
	typedef butils::ulong ulong;
	uint32_t counts;
public:
	WLocus() { clear(); };
	~WLocus(void) { clear(); };		// prevent vtable allocation
	inline void clear() { counts = 0; Locus::clear(); };

	FFStatus read( FILE *f , bool ignoreWeight = false ) { 
		using namespace butils::fastfile;
		// read a single rcord from a stream
		FFStatus s=Locus::read(f);
		if ( (s == SEOF) || (s == SERR) ) return(s);
		if (s == SEOL) counts = 1; 
		if (s == SEOR && !ignoreWeight ) {
			ulong t=0;
			s = getInt(f, &t); counts = t;
		}
		return(s);
	};

	FFStatus readln( FILE *f , bool ignoreWeight = false ) {
		using namespace butils::fastfile;
		// read a single rcord and rest of line... from a stream
		FFStatus s=read(f, ignoreWeight );
		if (ignoreWeight) counts = 1;
		if (s==SEOF) return( s );
		if (s==SEOR) s=clearLine(f);
		return( s );
	}

	bool	writeln (FILE *f ) { // slow method used for debuffing
		if (fprintf(f,"%s\t%lu\t%lu\t%lu\n",chrom.c(),(ulong)st,(ulong)en,(ulong)counts)<1) return(false); 
		fflush(f); return(true);
	}

};

typedef std::vector<WLocus> vWLocus;
