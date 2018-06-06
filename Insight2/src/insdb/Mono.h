#pragma once

#include <cstdio>	// FILE
#include <climits>
#include <cstdlib>		// malloc / free
#include <cstring>		// strcpy etcv...

#include <iostream>		// low performance IO
#include <string>		// low performance string manips
#include <vector>
#include <map>

#include "util.h"
#include "Locus.h"
#include "butils/butils.h"

// allows for 255 tags, plus a missing data value (0). If we need mroe tags, we can change dbElem to ushorts, which
//	doubles memory footprint to 6gb, but allows for 65536 tags... (or 1/2 precision prob class, phalf=16bits).
class Mono  {
public:
	typedef	unsigned char dbElem;						// databe element type... unsigned is important as chars are signed by default.
	typedef butils::stringS	stringS;
	typedef butils::RTable<stringS, dbElem>	RTableS;	// high performance map from index to value, oen for string one for double...
	typedef butils::RTable<double, dbElem>	RTableD;
private:
	typedef butils::ulong	ulong;
	typedef	dbElem *dbPtr;								// pointer into databsae
	typedef insdb::Coord	Coord;
	typedef insdb::Locus	Locus;
	//
	dbElem	*db;										// database of values, genome wide, indexed to 0. 1 entry per genomic position.
	Coord	dbSize;										// number of elements in DB
	RTableS	mapString;									// map an entry in DB (index) to an output value. Must be fast, used often...
	RTableD mapDouble;

	// used for rare / low speed data access
	std::vector<Locus>			loci;				// list of chromosome extents... maintained for debugging.
	std::map<std::string,Coord>	chromDbOffset;		// Lookup, chromosome to position, can be slow as it is rarely used

	// used for high speed data access, unsafe!
	dbPtr		chromBase;							// pointer to position 0 of "current" chromosome
	dbPtr		curPos;								// "current" index into database....
	stringS		curChrom;							// string copy odf "current" might be replaced by pointer to table, if needed...

	void clear( bool first = false );
	long long read64(FILE *Handle, void *buf, long long Bytes, int verbose = 0 );
	// long long read64(int Handle, void *buf, long long Bytes, int verbose = 0 );

public:
	typedef butils::fastfile::FFStatus FFStatus;
	static const std::size_t dbElemMax = std::numeric_limits<dbElem>::max();
	Mono() { clear( true ); };
	virtual	~Mono(void) {clear( false ); };

	// retrives pointers to tables, fast, but potentially unstable
	void getTables( void **tString, void **tDouble) {
		if (tString) *tString = (void *) &mapString; 
		if (tDouble) *tDouble= (void *) &mapDouble;	
	}

	// copy actual tables, slower, but stable and tables are small...
	void copyPriors(RTableS *tString, RTableD *tDouble) {
		if (tString) *tString = mapString;
		if (tDouble) *tDouble = mapDouble;
	}

	FFStatus		fRead( std::string &DBBase, int Verbosity = 0 );

	// find the db position of the start of chrom, add Start as an offset and return pointer to appropriate tag, or NULL.
	// increment current pos nad reutrn the next tag... maintain these for backwards compatbility
	inline const char		*getFirst(const char *chrom, Coord Start)	{ return((const char *)elemToString(getFirstElem(chrom, Start))); };
	inline const char		*getNext( void )							{ return((const char *)elemToString(getNextElem())); };

	// return smaller, and mor informative, dbElem values, these are indices into table, use maps to covert indices to values, very fast.
	dbElem					getFirstElem(const char *chrom, Coord Start);
	inline dbElem			getNextElem(void) { curPos++; return(*curPos ); }

	const char	*			elemToString(dbElem E)  const { return(E==0?NULL:&( (*(mapString[E]))[0]) ); };
	inline const double	*	elemToDouble(dbElem E)  const { return(E==0?NULL:mapDouble[E]); };

	inline double			majDouble(dbElem E)		const { return(E==0?0.0:*(mapDouble[E])); };
	inline const char *		majString(dbElem E)		const { return(E==0?"0":elemToString(E)); };


	static void		LociInfo( Mono &m, ulong &NumLoci, ulong &NumPos ) { Locus::LociInfo<Locus>(m.loci,NumLoci,NumPos); return; }
	// could be made faster but requires FastFile support for writing lu.
	inline static bool writeInp( FILE *f, const char *Chr, Coord Pos, const char *Tag ){
		if (!f) return(false); return( fprintf(f,"site\t%s:%lu\tM\t%s\n",Chr,(ulong)Pos,Tag)>0); }
	// used for debugging...
	bool writeLn( FILE *f ) { 
		if (!f) return(false); fprintf(f,"%s\tptr=%llx\tTagIDX=%lu\tTag=%s\n",curChrom.c(),(long long int)curPos,(ulong)*curPos,(*curPos==0?"null":elemToString(*curPos))); return( true); }
};



void Mono::clear(bool first) {
	db = (first ? NULL : db);
	if (db != NULL) { free(db); db = NULL; };
	dbSize = (Coord)0;
	mapString.clear(); mapDouble.clear();
	loci.clear();
	chromDbOffset.clear();
	chromBase = curPos = db;
	curChrom[0] = (dbElem)0;
}

// deal with fact that in MS land x64 ints AND longs are 32 bits....
// WHY? int ->32, long -> 64.... + room on upside for longlong -> 64 or 128... sigh...
long long
Mono::read64(FILE *Handle, void *buf, long long BytesToRead, int verbose) {
	//long long	max_read= (long long) INT_MAX;
	using std::cerr;
	using std::endl;
	long long	bytes_read = 0;
	//int			read_sz, s;
	//char		*p = (char *) buf;
	if (verbose > 0) { cerr << "Read64, handle " << Handle << " requested " << BytesToRead << " bytes." << endl; }
	if (BytesToRead < 0) return(-1);
	if (verbose > 1) { cerr << "\t requesting " << BytesToRead << " bytes... "; cerr.flush(); }
	bytes_read = fread(buf, 1, (size_t) BytesToRead, Handle);
	if (verbose > 1) { cerr << "\t got " << bytes_read << " bytes. " << endl; cerr.flush(); }
	return(bytes_read);
}


Mono::FFStatus
Mono::fRead(std::string &DBbase, int v) {
	using std::cerr;
	using std::endl;
	std::string fname;
	using namespace butils;
	clear();

	//read chromosome extent list
	fname = DBbase + ".chroms";
	{	
		if (v > 0)  cerr << "Reading chromosome file " << fname << endl;
		FFStatus s =fastfile::FileToVec(fname.c_str(), loci); if (s != fastfile::SEOF) return fastfile::SERR;
		if (v > 0)  cerr << "\tFile read, identifying offset " << fname << endl;
		for (unsigned int i = 0; i< loci.size(); i++) {
			// use a std::map to map chromosome names to starting offset in the database. This is slow, 
			//	but we always traverse chromosomes in order, so lookups are rare, and only occur when chromosomes
			//	change during a whole genome scan change...
			if (v > 1) cerr << "\t" << loci[i].chrom.c() << "\t" << dbSize << endl;
			chromDbOffset[std::string(loci[i].chrom.c())] = dbSize; dbSize += loci[i].en;
		}
		if (v > 0) cerr << "Total positions : " << dbSize << endl;
	}

	// Allocate the database and read the main data, one (byte) for each position in the gneome
	fname = DBbase + ".db";
	// magic number 16 is just to add a little at the end of the long buffer. We should never access it, but it provides a defensive margin.
	if (v > 0)  cerr << "Allocating  " << dbSize *sizeof(dbElem) << "bytes for mono array." << endl;
	db = (dbPtr)calloc(dbSize + 16, sizeof(dbElem)); if (db == (dbPtr)NULL) return (fastfile::SERR);
	{
		if (v > 0)  cerr << "\t reading from filesd based at " << fname << endl;
		FILE *fi = fastfile::ffopen(fname.c_str()); if (fi == NULL) return(fastfile::SERR);
		long long nread = read64(fi, db, (long long)dbSize*sizeof(dbElem), v); fastfile::ffclose(&fi); 
		if ((long long unsigned int)nread != dbSize*sizeof(dbElem)) return(fastfile::SERR);
		if (v > 0)  cerr << "Done. Got " << nread << " bytes " << endl;
	}

	// Read a list of values, up to one for each possible data base value, generally range is [1-255] with 0 beign a null (missing value)
	// this is slow, but these tables are small.
	fname = DBbase + ".tags";
	if (!mapString.readSimple(fname.c_str(),v)) return (fastfile::SERR);
	if (!mapDouble.readSimple(fname.c_str(),v)) return (fastfile::SERR);

	return(fastfile::SEOF);
}

Mono::dbElem
Mono::getFirstElem(const char *chrom, Coord Offset) {
	//	strings are short, so strcmp is fast and usually returns a match (0)
	if (strcmp(chrom, curChrom.c()) != 0) {
		// only lookup starting position for a chromosome when chromosome changes, this is rare
		strcpy(curChrom.c(), chrom);
		Coord index = chromDbOffset[chrom];
		chromBase = &(db[index]);					// get new base for chromosome
	}
	curPos = &(chromBase[Offset]);					// get position for offset from desired chromosome in the db
	return(*curPos );								// the db returns an index into the tags list, return pointer to corrsponding tags string... 0-len string returns NULL pointer (no value)
}
