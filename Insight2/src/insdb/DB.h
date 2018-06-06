#pragma once

#include <cstdio>
#include <cstring>
#include <cstdlib>		// atoi
#include <cassert>		// debugging

#include <vector>
#include <iostream>		// low performance IO
#include <string>		// low performance string manips
#include <map>			// low performance associateive array
#include <limits>

#include "butils/butils.h"
#include "insdb/insdb.h"

class DB {
public:
	typedef insdb::Locus	Locus;
	typedef insdb::Mono		Mono;
	typedef insdb::Poly		Poly;
	typedef insdb::Block	Block;
	typedef	insdb::Coord	Coord;
	typedef butils::ulong	ulong;
	typedef butils::ullong	ullong;
private:
	Mono				m;		// monomorphic sites (large footprint / 3GB)
	std::vector<Block>	vb;		// block list with theta & lambdas
	std::vector<Poly>	vp;		// polymorphisms
	std::vector<Poly>	vpn;	// neutral polymorphisms
public:
	enum  streamtype	{ AllData, NeutralPolysOnly };
	enum  cumestatus	{ NewLocus, NewBlock, BlocksDone };
	struct polycount	{ insdb::Poly *poly;  uint32_t count; };
	//struct monocount	{ Mono::dbElem iMono; uint8_t  count;};
	struct bigsite {
		const insdb::Poly  *polyptr; insdb::Mono::dbElem monoval; Coord monopos; bool ismono;
		void clear() { polyptr = NULL; monoval = 0; monopos = 0; ismono = false; };	// used to generate inp files...
		bigsite() { clear(); }
	};
	struct bigblock {
		insdb::Block			*head;
		insdb::DB	 const		*db;
		std::vector<bigsite>	sites;	// 21 bytes per site, 10-20 x larger than fastblock....
		inline bigblock()		{ clear(); };
		inline bigblock(insdb::DB const *DB) { clear(false ); db = DB; };
		inline bool numPos( ulong &Mono, ulong &PolyL, ulong &PolyH );
		inline void clear(bool SitesOnly = false) { if (!SitesOnly) { head = NULL; }; sites.clear(); }; // NEVER clear the DB....
		bool dump(FILE *fout, bool PolyOnly=false );	// Mono lookup is needed to resolve indiced into Mono table
	};
	struct fastblock {
		// these poitners are to persistant objects. Its ok to just clear the pointers.
		insdb::Block				*head;
		uint32_t					items;
		std::vector<uint32_t>		mono;	// 1 <= 4 bytes / site
		// std::vector<uint32_t>		mono(Mono::dbElemMax + 1);	// 1 <= 4 bytes / site
		std::vector<polycount>		poly;	// 12 bytes per site 
		//inline ullong dataCount(void)	{ return(mono.size() + poly.size()); };
		inline void  clear(void) { head = NULL; items = 0;  poly.clear(); for (ulong i = 0; i < mono.size(); i++) mono[i] = 0; }
		fastblock() : mono(Mono::dbElemMax + 1) { clear(); }
	};
	struct sitecounter {			// used to identify the number of informative sites encountered durign a stream...
		ulong		M, L, H;	// weighted and unweighted counts.
		uint64_t	wM,wL,wH;	// weighted and unweighted counts.
		void clear() { M=L=H=0;wM=wH=wL=0; return; }
		uint64_t	numSites()	const { return(  M +  L +  H ); };
		uint64_t	numSitesW()	const { return( wM + wL + wH ); };
		sitecounter() { clear(); }
	};
private:
	// values used during streaming
	ulong				iblock, ipoly;			// position indices
	streamtype          stype;
	Locus				t_roi;		// temporary variable used in Stream accumulatuion. Class scope to prevent iterated construction.
	sitecounter			infSites;
public:
	int				ReadDB(std::string DBdir, int Verbose = 0);
	inline void		locStreamInit(streamtype t) { stype = t;  iblock = ipoly = 0; infSites.clear(); return; };				// prepare to iterate over sorted loci
	cumestatus		locStreamFast(Locus &L, fastblock &Cumulant, uint32_t Count = 1);	// just return sufficient stats for monos. fast.
	cumestatus		locStreamBig( Locus &L, bigblock  &Cumulant);	// return stats and positions... test every position... slower...
	const sitecounter &infSiteCount( ) { return( infSites ); };		// read only - number of informative sites found in last stream, cleared by StreamInit
	// need access to tables to lookup mono values as mono lookups return index into table, poly lookups return pointer to complete poly object....
	inline double	monoMajDouble(Mono::dbElem e) { return(m.majDouble(e)); };
	inline const char *monoMajString(Mono::dbElem e) { return(m.majString(e)); };
	// hack......
	std::vector<Poly> &NeutralPolys() { return vpn; };
	void			copyAncPri(Mono::RTableS *TableS, Mono::RTableD *TableD) { m.copyPriors(TableS, TableD); };
	inline const double *AncPriVal(bool IsMono, Mono::dbElem Index) {
		if ( Index == 0 ) return NULL;	// special value!
		// OK, assume that Mono and Poly share index tables... true for now...
		return m.elemToDouble(Index);
	}
};

bool DB::bigblock::numPos(ulong &Mono, ulong &PolyL, ulong &PolyH) {
	Mono = PolyL = PolyH = 0;
	for (ulong i = 0; i < sites.size(); i++) {
		bigsite &s = sites[i];
		if (s.ismono) {	Mono++;
		} else {
			if (s.polyptr->freq == 'L') PolyL++;
			if (s.polyptr->freq == 'H') PolyH++;
		}
	}
	return true;
}

bool DB::bigblock::dump( FILE *fout, bool PolyOnly ) {
	bool ok = true; bigsite *s;
	if (sites.size() < 1) return(true);	// block w/ no sites, this is not an error, just dump nothing
	if ((!head) || (head->writeInp(fout) != butils::fastfile::SEOL)) return(false);	// write the header
	for (ulong i = 0; (ok && (i < sites.size())); i++) {							// write each site, in order...
		s = &(sites[i]); ok = true;
		if (!s->ismono) ok = s->polyptr->writeInp(fout);
		if (s->ismono && !PolyOnly) ok = Mono::writeInp(fout, head->chrom.c(), s->monopos, db->m.elemToString(s->monoval));
	}
	return(ok);
}

// Use: Call locStreamInit - Clear Cumulant, 
//	read a locus, call locStream()
//		== BlocksDone, dump cumulant and terminate
//		== NewLocus, read a new locus and continue
//		== NewBlock, dump and clear cumulant then call again
//  untill out of loci.... then dump cumulant if nonempty and exit.
DB::cumestatus DB::locStreamBig(Locus &L, bigblock &Cumulant) {

	// find the first block that does not end before Locus starts
	if (!Locus::SkipPrevious<Block>(vb, L, iblock)) return DB::BlocksDone;	// no more blocks? then we are all done..

	//Locus t_roi;	// put in DB class to avoid iterated construction....
	if (!L.intersect(vb[iblock], t_roi)) return(DB::NewLocus);				// no overlap between locus and next block? get the next locus...

	if (Cumulant.head == NULL) Cumulant.head = &(vb[iblock]);				// Freshly cleared cumulant! 
	if (Cumulant.head != &(vb[iblock])) return DB::NewBlock;				// new roi is not in cumulant's block, allow caller to clear cumulant

	assert(!t_roi.isEmpty());
	bigsite bs; const char *chrom=Cumulant.head->chrom.c();
	std::vector<Poly> &t_vp = (stype == NeutralPolysOnly ? vpn : vp);
	Mono::dbElem v = m.getFirstElem(t_roi.chrom.c(), t_roi.st);
	for (Coord pos = t_roi.st; pos < t_roi.en; pos++) {
		// process monos, if requested
		if ((stype == AllData) && (v != 0)) {
			bs.clear();
			bs.ismono = true; bs.monopos = pos; bs.monoval = v; Cumulant.sites.push_back(bs); infSites.M++; infSites.wM += 1;
		};
		// process polys , 
		Locus::SkipPrevious<Poly>(t_vp, chrom, pos, ipoly);
		if (t_vp[ipoly].OverlapsPos(chrom, pos)) {
			bs.clear();
			bs.ismono = false; bs.polyptr = &t_vp[ipoly]; Cumulant.sites.push_back(bs);
			if (bs.polyptr->freq == 'L') { infSites.L++; infSites.wL += 1; } else { infSites.H++; infSites.wH += 1; };
		};
		v = m.getNextElem();
	};

	// Remove the roi from the current locus. If remainder is empty, we need a new locus
	//	Otherwise the locus has spanned blocks, and we need a new block....
	L.advanceBy(t_roi);
	return(L.isEmpty() ? DB::NewLocus : DB::NewBlock);
}

DB::cumestatus DB::locStreamFast(Locus &L, fastblock &Cumulant, uint32_t Count) {

	// find the first block that does not end before Locus starts
	if (!Locus::SkipPrevious<Block>(vb, L, iblock)) return DB::BlocksDone;	// no more blocks? then we are all done..
	if (!L.intersect(vb[iblock], t_roi)) return(DB::NewLocus);				// no overlap between locus and next block? get the next locus...

	// Make sure current block head matches cumulant head...
	if (Cumulant.head == NULL) Cumulant.head = &(vb[iblock]);				// Freshly cleared cumulant! 
	if (Cumulant.head != &(vb[iblock])) return DB::NewBlock;				// new roi is not in cumulant's block, allow caller to clear cumulant

	// get the right set of polys, all, or neutral only...
	assert(!t_roi.isEmpty());
	std::vector<Poly> &t_vp = (stype == NeutralPolysOnly ? vpn : vp);

	polycount tmp_poly;
	// Get polys in roi. No more polys? That is ok.
	if (Locus::SkipPrevious<Poly>(t_vp, t_roi, ipoly)) {
		while (ipoly < t_vp.size() && t_roi.Overlaps(t_vp[ipoly])) {
			tmp_poly.poly = &(t_vp[ipoly]); tmp_poly.count = Count;
			if (tmp_poly.poly->freq == 'L') { infSites.L++; infSites.wL += Count; } else { infSites.H++; infSites.wH += Count; };
			Cumulant.poly.push_back(tmp_poly); ipoly++;
			Cumulant.items += Count;
		};
	};

	// Only get monos if we are asking for all data, skip if we are asking for neutral polys & blocks only.
	Mono::dbElem v = 0;
	if (stype == AllData) {
		// Get monos in roi. Monos span entire genome, should should cover any ROI.
		// only pushback non-0 values, 0 represents missing data.
		if ((v = m.getFirstElem(t_roi.chrom.c(), t_roi.st)) != 0) { 
			Cumulant.mono[v] += Count; Cumulant.items += Count; infSites.M++; infSites.wM += Count; };
		for (Coord i = t_roi.st + 1; i < t_roi.en; i++)
			if ((v = m.getNextElem()) != 0) { 
				Cumulant.mono[v] += Count; Cumulant.items += Count; infSites.M++; infSites.wM += Count;};
	};

	// Remove the roi from the current locus. If remainder is empty, we need a new locus
	//	Otherwise the locus has spanned blocks, and we need a new block....
	L.advanceBy(t_roi);
	return(L.isEmpty() ? DB::NewLocus : DB::NewBlock);
}


int  DB::ReadDB(std::string DBdir, int Verbose) {
	using std::endl;
	using std::cerr;
	using butils::VecArray;
	int v = Verbose; bool ok; std::string fname;

	// Accelerated fetch of the Block file (small)
	fname = DBdir + "/block.bedg";
	if (v > 0) { cerr << "Processing Block file " << fname << "..."; cerr.flush(); }
	ok = VecArray<Block>::FileToVecCached(fname, vb);
	if (!ok) {
		auto i = vb.size();  cerr << "Error processing file, exiting. Num Lines read " << i << endl;
		if (i > 1) { cerr << "\tN-1  line : " << vb[i - 2].toString() << endl; }
		if (i > 0) { cerr << "\tLast line : " << vb[i - 1].toString() << endl; }
		cerr << endl << endl;
		return(-1);
	}
	if (v > 0) { ulong l = 0, p = 0; Locus::LociInfo<Block>(vb, l, p); cerr << "found " << l << " loci covering " << p << " genomic positions." << endl; };

	// Accelerated fetch of the monomorphic site files (large 3GB)
	fname = DBdir + "/monoDB";
	if (v > 0) { cerr << "Processing Monomorphism database " << fname << "..."; cerr.flush(); }
	ok = (m.fRead(fname, v) == butils::fastfile::SEOF);
	if (v > 0) { ulong l, p; Mono::LociInfo(m, l, p); cerr << "found " << l << " loci covering " << p << " genomic positions." << endl; };
	if (!ok) {
		cerr << "Error processing file, exiting.\n\n"; return(-4);
	}
	// HACK: share the tables from mono to poly before reading POLY files. Tables are coppied, but small.
	{	void *a, *b; m.getTables(&a, &b); Poly::setTables(a, b); }

	// Accelerated fetch of Polymorphism file (medium)
	fname = DBdir + "/poly.bedg";
	if (v > 0) { cerr << "Processing poly file " << fname << "..."; cerr.flush(); }
	ok = VecArray<Poly>::FileToVecCached(fname, vp);
	if (!ok) {
		auto i = vp.size();  cerr << "Error processing file, exiting. Num Lines read " << i << endl;
		if (i > 1) { cerr << "\tN-1  line : " << vp[i - 2].toString() << endl; }
		if (i > 0) { cerr << "\tLast line : " << vp[i - 1].toString() << endl; }
		cerr << endl << endl;
		return(-2);
	}
	if (v > 0) { ulong l = 0, p = 0; Locus::LociInfo<Poly>(vp, l, p); cerr << "found " << l << " loci covering " << p << " genomic positions." << endl; };

	// Accelerated fetch of Neutral Polymorphism file (medium)
	fname = DBdir + "/polyn.bedg";
	if (v > 0) { cerr << "Processing poly file " << fname << "..."; cerr.flush(); }
	ok = VecArray<Poly>::FileToVecCached(fname, vpn);
	if (!ok) {
		auto i = vpn.size();  cerr << "Error processing file, exiting. Num Lines read " << i << endl;
		if (i > 1) { cerr << "\tN-1  line : " << vpn[i - 2].toString() << endl; }
		if (i > 0) { cerr << "\tLast line : " << vpn[i - 1].toString() << endl; }
		cerr << endl << endl;
		return(-3);
	}
	if (v > 0) { ulong l = 0, p = 0; Locus::LociInfo<Poly>(vpn, l, p); cerr << "found " << l << " loci covering " << p << " genomic positions." << endl; };

	return(0);
}
