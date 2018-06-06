#pragma once
#include <vector>
#include <algorithm>	// std::remove()
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "stringer.h"

// Creates a map from a set of names chomosomes and positions to a linear index.
//	Generally used to simplify access to sorted bed files the span the entire range
//	of a reference genome, in sorted order.
//	
// typedef Tpos		pos;	// enough positions to specify every location hg18 is about 3billion. Larger chromosomes might require promotion to uint64_t
// typedef Tid		chID;	// must contain number of chromosomes 255 should be enough... for most uses, increase if more chroms / uniquely named loci....

template< typename Tpos = uint32_t, typename Tid = uint8_t >
class posMapper
{
public:
	typedef stringer	string;
private:
	static constexpr const char *hg19auto =	\
		"chr1    0       249250621\nchr10   0       135534747\nchr11   0       135006516\n"
		"chr12   0       133851895\nchr13   0       115169878\nchr14   0       107349540\n"
		"chr15   0       102531392\nchr16   0       90354753\n chr17   0       81195210\n"
		"chr18   0       78077248\n chr19   0       59128983\n chr2    0       243199373\n"
		"chr20   0       63025520\n	chr21   0       48129895\n chr22   0       51304566\n"
		"chr3    0       198022430\nchr4    0       191154276\nchr5    0       180915260\n"
		"chr6    0       171115067\nchr7    0       159138663\nchr8    0       146364022\n"
		"chr9    0       141213431";
	std::vector<string>	chrName;	// locus (chrom) name
	std::vector<Tpos>	chrStart;	// linear position of start of chromosome
	std::vector<Tpos>	chrLen;		// Length of locus chrEnd-chrOffset+1, also number of positiosn occupuied in linear map.
	std::vector<Tpos>	chrOst;		// Offset of start of locus into chrom, usually 0 (bed is 0 indexed, start closed)
	std::vector<Tpos>	chrOen;		// Offset of end of locus into chrom, bed end position -1 (bed is end open)
	// Tpos				lastPos;	// last linear position
	inline bool idOK( Tid ID ) const { return( ID<numIDs()); };
	inline bool idBad(Tid ID) const { return(ID>=numIDs()); };
	inline bool posBad(Tpos Pos) const { return(Pos> lastPos()); };
public:
	void clear() { chrName.clear(); chrOst.clear(); chrOen.clear(); chrStart.clear(); chrLen.clear(); };
	posMapper() { clear(); };
	~posMapper() { };
	// Number of loci / chromosomes in map
	inline Tid  numIDs( ) const { return (Tid) chrName.size(); }
	// Get the linear positions based on ID, or sumamry for entire chromosome set.
	inline Tpos  firstPos(Tid ID) const {  if (idBad(ID)) return 0; return(chrStart[ID]); }
	inline Tpos  lastPos( Tid ID ) const { if (idBad(ID)) return 0; return(chrStart[ID] + chrLen[ID] - 1); }
	inline Tpos  lastPos() const { Tid n=numIDs(); return (n<1?0:lastPos(n-1)); };
	// Defaults to autosomal HG19. Low performance IO.. but OK, this is typically < 255 short lines....
	bool read( const stringer &fname="") {
		std::istream *sin;
		sin = (fname == "" ? (std::istream *) new std::istringstream(posMapper::hg19auto) : (std::istream *) new std::ifstream(fname, std::ios_base::in));
		string aline; std::vector<string> toks;
		while (sin->good()) {
			std::getline( *sin, aline ); aline.trim(); 
			if (aline.length() < 1 ) continue;	// skip blanks
			if (aline[0]=='#') continue;		// skip comments
			aline.tokenize(toks); toks.erase( std::remove(toks.begin(), toks.end(), ""), toks.end());	// tokenize and remove blanks
			Tid tok_num=(Tid)chrName.size(); Tpos st = (Tpos) std::stoull(toks[1]), en = (Tpos) std::stoull(toks[2]);
			if (numIDs()==0) { chrStart.push_back(0); } else { chrStart.push_back(chrStart[tok_num-1]+chrLen[tok_num-1]); };
			chrName.push_back(toks[0]); chrOst.push_back(st); chrOen.push_back(en-1); chrLen.push_back(en-st);
		};
		if (sin->eof()) return true;
		return false;
	}

	// These are fast, OK to use as accessor methods...
	inline bool idInfo(Tid ID, Tpos &FirstPos, string &Chr, Tpos &Start, Tpos &Length ) const {
		if (idBad(ID)) return false; FirstPos=chrStart[ID]; Chr=chrName[ID]; Start=chrOst[ID]; Length=chrLen[ID]; return true; };
	inline bool posInCh(Tid ID, Tpos Pos ) const {
		if (idBad(ID) ) return false;
		if ( Pos < firstPos(ID) ) return false;
		if ( Pos > lastPos(ID)  ) return false; 
		return true; };

	// these are slow, avoid using when possible
	bool posToID( Tpos pos, Tid &ID ) const { if (posBad(pos)) return false; ID=0; while( !posInCh( ID,pos) ) ID++; return true; };
	bool chrToID( const string &chr, Tid &ID ) const {ID=0; while (idOK(ID)) { if (chr==chrName[ID]) return true; ID++; }; return false; };
};

