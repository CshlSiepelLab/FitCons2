#pragma once

#include "util.h"
#include "butils/butils.h"

#include <cstring>
#include <vector>
#include <string>

class Locus {
public:
	typedef butils::stringT stringT;				// sigh, why wont using work here?
	typedef butils::fastfile::FFStatus FFStatus;	// sigh, why wont using work here?
	stringT chrom;
	Coord	st;
	Coord	en;			// end is a reserved word....
	inline void clear() {chrom[0]='\0'; st=en=(Coord)0; };
public:
	Locus()	{ clear(); };
	virtual ~Locus(void) { clear(); };
	inline bool isEmpty( void ) const { return( st >= en ); }
	inline bool isBefore( const char *Ch, Coord pos ) const	{ 
		return( compare( Ch,pos ) < 0 ); }
	inline bool OverlapsPos( const char *Ch, Coord pos ) const	{ 
		return( compare(Ch,pos)==0); }
	inline bool Overlaps(const Locus &L ) const {
		return(compare(L) == 0);
	}

	// return -1, 0, 1 if this is entirely before, overlapping or entirely after position
	inline int compare( const char *Ch, Coord Pos ) const {
		int s=strcmp(chrom.c(), Ch);
		if (s!=0) return(s);			
//		if (s<0 ) return( -1 );			// pos is on an earlier chromosome
//		if (s>0 ) return(  1 );			// pos is in a later chromosome
		if (en <= Pos ) return( -1 );	// pos is after the end of the locus (position en is NOT prt of locus)
		if (st >  Pos ) return(  1 );	// pos is before start of locus (position st IS part of locus)
		return 0;						// pos is contained in locus
	}

	// return -1, 0, 1 if this is entirely before, overlapping or entirely after L.
	inline int compare( const Locus &L ) const {
		int s = strcmp( chrom.c(), L.chrom.c());
		if (s != 0 ) return( s );
		if (en<=L.st) return( -1 );
		if (st>=L.en) return( 1 );
		return 0;
	}

	// advance the start of this th the end of L, if possible
	inline bool advanceBy( const Locus &L ) {
		int ci = strcmp( chrom.c(), L.chrom.c());
		// wrong chromosome, removing L from this just this...
		if (ci==0) advanceBy_unsafe( L );
		if (ci<0) st = en;
		return (st<en);
	}

	// Unsafe operations assume that *this and Locus &L are on the same chromosome (that higher level code has already checked)
	//	and that there IS an sctual overlap (non null result).
	inline void advanceBy_unsafe( const Locus &L ) {
		if (L.en>st) { st = (L.en<en?L.en:en); };
		return;
	}

	inline void intersect_unsafe( const Locus &L, Locus &Overlap ) const {
		char *a = Overlap.chrom.c();	const char *b = chrom.c();
		if (strcmp(a,b)) // only returns 0 (false) if they are ==
			strcpy(a, b);
		Overlap.st=(st < L.st ? L.st : st   );	// max
		Overlap.en=(en < L.en ? en   : L.en );	// min
	}

	inline bool intersect( const Locus &L, Locus &Overlap ) const {
		if (compare(L)!=0 ) return( false );
		intersect_unsafe(L,Overlap);
		return( true );
	}


	FFStatus read( FILE *f )	{ 
		// read a single rcord from a stream
		using namespace butils;
		clear(); FFStatus s= fastfile::getString(f,chrom.c());
		if (s==fastfile::SEOF) return( s );
		if (s						!=fastfile::SEOR ) return fastfile::SERR;
		st = en = 0;
		if (fastfile::getInt(f, &st) != fastfile::SEOR) return fastfile::SERR;
		return( fastfile::getInt(f,&en) );
	};

	FFStatus readln( FILE *f ) {
		// read a single rcord and rest of line... from a stream
		using namespace butils;
		FFStatus s=read(f);
		if (s==fastfile::SEOF) return( s );
		if (s==fastfile::SEOR) s=fastfile::clearLine(f);
		return( s );
	}

	bool	writeln (FILE *f ) const { // slow method used for debuffing
		if (fprintf(f,"%s\t%lu\t%lu\n",chrom.c(),(ulong)st,(ulong)en)<1) return(false); 
		fflush(f); return(true);
	}

	// slow method used for debugging
	std::string toString(void) const {
		using std::string;
		return(string(chrom.s) + " " + std::to_string(st) + " " + std::to_string(en));
	}

	// For vectors of objects derived from the locus class, get num loci and number of positions covered by elements in vector.
	template<typename T>
	static void LociInfo(const std::vector<T> &V, ulong &NumLoci, ulong &NumPos) {
		NumLoci = NumPos = (ulong)0;
		for (std::size_t i = 0; i<V.size(); i++) { NumLoci++; NumPos += (ulong)(V[i].en - V[i].st); }
	}

	// For vectors of objects derived from the locus class, advance to the first index of a locus that is not entirely before Loc. 
	template<typename T>
	static bool SkipPrevious(const std::vector<T> &V, const char *Chrom, Coord Pos, ulong &Index) {
		ulong imax = (ulong)V.size();
		// if vector is empty return false
		if (imax == 0) { Index = 0;  return false; }
		imax--;
		// if Index is already at the end of the vector, just check the last item.
		if (Index >= imax) { Index = imax; return (V[Index].compare(Chrom,Pos) >= 0); };
		int c = 0;
		// relies on short circuiting && for safety. Binary search might be faster for random access, but we are 
		//	often incrementing by small amounts as we scan through genomic positions in order.
		while ((Index < imax) && ((c = V[Index].compare(Chrom,Pos)) < 0)) Index++;
		// If weseek a position in the last block, the penultimate block increments Index to imax, but then short curcuts before testing the block, so we miss it! 0.16d...
		if (Index == imax) c = V[Index].compare(Chrom, Pos);	// if the last, force a test
		return (c >= 0); // only false when Index >= Imax and compare is <0; This means No more data.
	}

	// For vectors of objects derived from the locus class, advance to the first index of a locus that is not entirely before Loc. 
	template<typename T>
	static bool SkipPrevious( const std::vector<T> &V, const Locus &L, ulong &Index)  {
		return(SkipPrevious(V, L.chrom.c(), L.st, Index));
	}
};

typedef ::std::vector<Locus> vLocus;
